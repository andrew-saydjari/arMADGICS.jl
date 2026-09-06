## E5: bright sky-line DETECTOR (AKS design, 2026-09-06) and its calibration
## against DR17's actual bright mask.
#
# AKS: "figure out what flux in our new units reproduces the bright_msk from DR17 best
# (in terms of specific overlap in the lines masked). because we have not throughput
# divided, it might be that this mask needs to be computed on a median smoothed spectrum
# or something like that. since we are just trying to mask bright skylines that are an
# order of magnitude too bright." ... "something closest to C, but ... identifying the
# brightest lines via some moving median, then outliers relative to quantiles, then
# probably grow the mask via dilation to catch the wings of the line."
#
# Three steps, all operating on the FULL wavetarg grid (NaN outside valid pixels) so
# that wavelength adjacency is physical — the builder's compressed submsk space splices
# chip gaps together, which would corrupt any moving window:
#   (1) LOCAL CONTINUUM: NaN-aware running median of median_sky over `window` pixels.
#       This absorbs the throughput/blaze shape we have NOT divided out.
#   (2) OUTLIERS vs QUANTILES: residual = median_sky - continuum, scaled by a running
#       robust spread (MAD, or IQR) of that residual; flag residual > k * spread.
#       Mode :ratio instead flags median_sky > f * continuum ("order of magnitude too
#       bright") — kept so the calibration can choose between them on evidence.
#   (3) DILATION: grow the flagged set by `dilation` pixels to cover line wings.
#
# The output is a MASK, not a scalar threshold: bright lines cannot be separated from a
# throughput-shaped continuum by any single flux value (that is exactly why the
# inherited constant failed, finding #35).
# Author - Andrew Saydjari (E5 pass 1)

"""
    running_median_nan(x, window) -> Vector{Float64}

NaN-aware running median over a centred window (odd length enforced). Positions whose
window holds no finite value return NaN.
"""
function running_median_nan(x::AbstractVector, window::Int)
    n = length(x)
    r = max(window ÷ 2, 1)
    out = fill(NaN, n)
    buf = Float64[]
    for i in 1:n
        empty!(buf)
        @inbounds for j in max(1, i - r):min(n, i + r)
            v = x[j]
            isfinite(v) && push!(buf, v)
        end
        isempty(buf) || (out[i] = median(buf))
    end
    return out
end

"""
    running_spread_nan(x, window; kind=:mad) -> Vector{Float64}

NaN-aware running robust spread of `x`: MAD scaled to sigma (1.4826) or IQR/1.349.
"""
function running_spread_nan(x::AbstractVector, window::Int; kind::Symbol=:mad)
    n = length(x)
    r = max(window ÷ 2, 1)
    out = fill(NaN, n)
    buf = Float64[]
    for i in 1:n
        empty!(buf)
        @inbounds for j in max(1, i - r):min(n, i + r)
            v = x[j]
            isfinite(v) && push!(buf, v)
        end
        length(buf) < 5 && continue
        if kind == :mad
            m = median(buf)
            out[i] = 1.4826 * median(abs.(buf .- m))
        else
            out[i] = (quantile(buf, 0.75) - quantile(buf, 0.25)) / 1.349
        end
    end
    return out
end


"""
    running_spread_fast(x, window; kind=:mad, stride=50) -> Vector{Float64}

Running robust spread evaluated on a coarse grid (every `stride` pixels) and linearly
interpolated. The quantity it estimates -- the local noise/throughput envelope -- varies
smoothly over hundreds of pixels, so the coarse grid is faithful while being ~`stride`x
cheaper than the per-pixel version (which is prohibitive at window ~1000).
"""
function running_spread_fast(x::AbstractVector, window::Int; kind::Symbol=:mad, stride::Int=50)
    n = length(x)
    r = max(window ÷ 2, 1)
    nodes = unique(vcat(collect(1:stride:n), n))
    vals = fill(NaN, length(nodes))
    buf = Float64[]
    for (t, i) in enumerate(nodes)
        empty!(buf)
        @inbounds for j in max(1, i - r):min(n, i + r)
            v = x[j]
            isfinite(v) && push!(buf, v)
        end
        length(buf) < 5 && continue
        if kind == :mad
            m = median(buf)
            vals[t] = 1.4826 * median(abs.(buf .- m))
        else
            vals[t] = (quantile(buf, 0.75) - quantile(buf, 0.25)) / 1.349
        end
    end
    # linear interpolation over the finite nodes
    out = fill(NaN, n)
    good = findall(isfinite, vals)
    isempty(good) && return out
    for idx in 1:(length(good)-1)
        a, b = good[idx], good[idx+1]
        i0, i1 = nodes[a], nodes[b]
        v0, v1 = vals[a], vals[b]
        @inbounds for i in i0:i1
            out[i] = v0 + (v1 - v0) * (i - i0) / max(i1 - i0, 1)
        end
    end
    @inbounds for i in 1:nodes[good[1]]
        out[i] = vals[good[1]]
    end
    @inbounds for i in nodes[good[end]]:n
        out[i] = vals[good[end]]
    end
    return out
end

"binary dilation of a Bool vector by `rad` pixels (grow TRUE regions)."
function dilate_msk(msk::AbstractVector{Bool}, rad::Int)
    rad <= 0 && return copy(msk)
    n = length(msk)
    out = falses(n)
    @inbounds for i in 1:n
        if any(view(msk, max(1, i - rad):min(n, i + rad)))
            out[i] = true
        end
    end
    return out
end

"""
    detect_bright_lines(x_full; window=101, k=8.0, dilation=4, mode=:mad,
                        spread_kind=:mad, ratio=10.0, spread_window=nothing)

Bright sky-line mask on the full grid. `x_full` is `median_sky` mapped onto wavetarg
with NaN at invalid pixels. Returns `(mask, continuum, residual, spread)`.
"""
function detect_bright_lines(x_full::AbstractVector; window::Int=101, k::Float64=8.0,
        dilation::Int=4, mode::Symbol=:mad, spread_kind::Symbol=:mad,
        ratio::Float64=10.0, spread_window=nothing)
    cont = running_median_nan(x_full, window)
    resid = x_full .- cont
    sw = isnothing(spread_window) ? window : spread_window
    spread = running_spread_nan(resid, sw; kind=spread_kind)
    base = falses(length(x_full))
    if mode == :ratio
        @inbounds for i in eachindex(x_full)
            if isfinite(x_full[i]) && isfinite(cont[i]) && cont[i] > 0 && x_full[i] > ratio * cont[i]
                base[i] = true
            end
        end
    else
        @inbounds for i in eachindex(x_full)
            if isfinite(resid[i]) && isfinite(spread[i]) && spread[i] > 0 && resid[i] > k * spread[i]
                base[i] = true
            end
        end
    end
    msk = dilate_msk(base, dilation)
    # never flag invalid pixels
    @inbounds for i in eachindex(msk)
        isfinite(x_full[i]) || (msk[i] = false)
    end
    return msk, cont, resid, spread
end

## ---------------------------------------------------------------- line-level metrics

"contiguous TRUE runs of a Bool vector, as UnitRanges (one per 'line')."
function mask_runs(msk::AbstractVector{Bool})
    runs = UnitRange{Int}[]
    i = 1
    n = length(msk)
    while i <= n
        if msk[i]
            j = i
            while j < n && msk[j+1]
                j += 1
            end
            push!(runs, i:j)
            i = j + 1
        else
            i += 1
        end
    end
    return runs
end

"""
    line_overlap_stats(ref_msk, new_msk, valid) -> NamedTuple

Line-level agreement between a reference bright mask (DR17) and a candidate mask,
restricted to pixels valid in both. A reference line counts as RECOVERED when the
candidate overlaps it at all; its IoU uses the union of candidate runs touching it.
Candidate lines touching no reference line are SPURIOUS. Also returns pixel IoU.
"""
function line_overlap_stats(ref_msk::AbstractVector{Bool}, new_msk::AbstractVector{Bool},
        valid::AbstractVector{Bool})
    R = ref_msk .& valid
    N = new_msk .& valid
    rruns = mask_runs(R)
    nruns = mask_runs(N)
    recovered = 0
    ious = Float64[]
    for rr in rruns
        touching = filter(nr -> !isempty(intersect(nr, rr)), nruns)
        if isempty(touching)
            push!(ious, 0.0)
            continue
        end
        recovered += 1
        # range arithmetic only: candidate runs are disjoint, so
        # |rr ∪ (∪touching)| = |rr| + Σ|nr| - |rr ∩ (∪touching)|
        inter = sum(nr -> length(intersect(nr, rr)), touching)
        uni = length(rr) + sum(length, touching) - inter
        push!(ious, uni == 0 ? 0.0 : inter / uni)
    end
    spurious = count(nr -> all(i -> !R[i], nr), nruns)
    pix_inter = count(R .& N)
    pix_union = count(R .| N)
    recall = isempty(rruns) ? NaN : recovered / length(rruns)
    precision = isempty(nruns) ? NaN : (length(nruns) - spurious) / length(nruns)
    f1 = (isnan(recall) || isnan(precision) || (recall + precision) == 0) ? NaN :
         2recall * precision / (recall + precision)
    return (n_ref_lines=length(rruns), n_new_lines=length(nruns), recovered=recovered,
        spurious=spurious, recall=recall, precision=precision, f1=f1,
        pixel_iou=pix_union == 0 ? NaN : pix_inter / pix_union,
        mean_line_iou=isempty(ious) ? NaN : mean(ious),
        ref_pix=count(R), new_pix=count(N), line_ious=ious)
end
