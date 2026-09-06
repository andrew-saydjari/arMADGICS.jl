## E5: calibrate the bright-line detector against DR17's ACTUAL bright mask.
# Target (AKS): maximum agreement in the specific LINES masked.
#
# PASS-2 CALIBRATION (pass 1 was unsound — its results are discarded, not reported):
#  - objective is now PIXEL IoU, with line recall/precision as diagnostics. Optimising
#    line-level F1 alone is gameable: dilation MERGES adjacent runs, raising line
#    "precision" while flagging 3x too many pixels (pass 1 chose bright=27% vs DR17 8.3%).
#  - the SCALE is global (per fiber), not a running MAD. A running MAD collapses in
#    line-free regions — where the residual is pure noise — so k*MAD there sits at the
#    noise floor and flags everything. That was the over-flagging mechanism in pass 1.
#  - k range extended far beyond pass 1, which pinned at its boundary (k=20, dil=8).
#  - :ratio mode DROPPED as structurally degenerate: it tests x > f*continuum, but the
#    samples are already continuum-subtracted upstream (sample_sky_defs.jl:69,
#    `fnew = fvec .- x_comp_lst[1]`), so the local continuum is ~0 (MEASURED: median
#    -0.036, global MAD 0.497 on f010; -0.025 / 0.166 on f600) and the ratio test
#    divides by ~zero. AKS's premise "we have not throughput divided" is handled
#    upstream: the throughput/blaze shape is removed by the per-spectrum continuum fit
#    against the starCont basis, not left in this quantity. The moving median of step
#    (1) is retained but is nearly a no-op here — window=0 (no continuum step) is in the
#    grid so the data decides whether it earns its place.
# Usage: julia --project=<plots env> e5_bright_calibrate.jl
using HDF5, Statistics, StatsBase, Printf, Dates

const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
include(ARM * "scripts/prior_build/e5_bright_line_detect.jl")

const OUT = get(ENV, "E5_OUT", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1")
const DR17 = get(ENV, "E5_OLD_PRIORS", "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/sky_priors")
const SCREENS = joinpath(OUT, "screens")

wavetarg = 10 .^ range((4.179 - 125 * 6.0e-6), step=6.0e-6, length=8575 + 125)

function load_fiber(adjfib)
    n = lpad(adjfib, 3, "0")
    p = joinpath(SCREENS, "median_sky_$n.h5")
    isfile(p) || return nothing
    ms = h5read(p, "median_sky")
    sub = Bool.(h5read(p, "submsk"))
    x = fill(NaN, length(wavetarg))
    x[sub] .= ms
    pb = joinpath(DR17, "APOGEE_skyline_bright_svd_120_f$n.h5")
    pf = joinpath(DR17, "APOGEE_skyline_faint_svd_120_f$n.h5")
    (isfile(pb) && isfile(pf)) || return nothing
    db = Bool.(h5read(pb, "submsk"))
    df = Bool.(h5read(pf, "submsk"))
    valid = sub .& (db .| df)
    return (adjfib=adjfib, x=x, sub=sub, dr17_bright=db, valid=valid)
end

fibers = [parse(Int, x) for x in split(get(ENV, "E5_CAL_FIBERS",
    "10,40,76,100,150,200,250,295,340,400,450,500,550,570,598,600"), ",")]
data = filter(!isnothing, [load_fiber(f) for f in fibers])
apo = filter(d -> d.adjfib <= 300, data)
lco = filter(d -> d.adjfib > 300, data)
@printf("calibration set: %d APO + %d LCO fibers\n", length(apo), length(lco))

windows = [0, 101, 301]          # 0 = no moving-median continuum step
ks = [3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0, 200.0]
qs = [0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08]
dils = [0, 1, 2, 3, 4, 5, 6, 8]

# standardized residual z = (x - continuum) / global_MAD, per (fiber, window)
zs = Dict{Tuple{Int,Int},Vector{Float64}}()
for d in data, w in windows
    cont = w == 0 ? zeros(length(d.x)) : running_median_nan(d.x, w)
    r = d.x .- cont
    fin = filter(isfinite, r)
    mad = 1.4826 * median(abs.(fin .- median(fin)))
    zs[(d.adjfib, w)] = r ./ max(mad, 1e-12)
end
println("precomputed standardized residuals"); flush(stdout)

function config_mask(d, w, dil, mode, k, q)
    z = zs[(d.adjfib, w)]
    cut = mode == :kmad ? k : quantile(filter(isfinite, z), 1 - q)
    base = falses(length(z))
    @inbounds for i in eachindex(z)
        (isfinite(z[i]) && z[i] > cut) && (base[i] = true)
    end
    msk = dilate_msk(base, dil)
    @inbounds for i in eachindex(msk)
        isfinite(d.x[i]) || (msk[i] = false)
    end
    return msk
end

function eval_config(set, w, dil; mode, k=0.0, q=0.0)
    ious = Float64[]; recs = Float64[]; precs = Float64[]; lious = Float64[]; bf = Float64[]
    for d in set
        msk = config_mask(d, w, dil, mode, k, q)
        st = line_overlap_stats(d.dr17_bright, msk, d.valid)
        push!(ious, st.pixel_iou); push!(recs, st.recall)
        push!(precs, st.precision); push!(lious, st.mean_line_iou)
        push!(bf, count(msk .& d.sub) / count(d.sub))
    end
    ff(v) = mean(filter(isfinite, v))
    return (pixiou=ff(ious), recall=ff(recs), precision=ff(precs), lineiou=ff(lious), bfrac=ff(bf))
end

results = Dict{String,Vector}()
for (name, set) in (("APO", apo), ("LCO", lco))
    isempty(set) && continue
    rows = []
    for w in windows, dil in dils
        for k in ks
            push!(rows, (mode=:kmad, w=w, k=k, q=NaN, dil=dil,
                eval_config(set, w, dil; mode=:kmad, k=k)...))
        end
        for q in qs
            push!(rows, (mode=:quant, w=w, k=NaN, q=q, dil=dil,
                eval_config(set, w, dil; mode=:quant, q=q)...))
        end
    end
    sort!(rows, by=x -> -(isfinite(x.pixiou) ? x.pixiou : -1))
    results[name] = rows
    println("\n===== $name : top 10 by PIXEL IoU vs DR17 bright mask =====")
    @printf("%-6s %4s %6s %6s %4s | %7s %7s %7s %7s %8s\n",
        "mode", "win", "k", "q", "dil", "pixIoU", "recall", "prec", "lineIoU", "bright%")
    for row in first(rows, 10)
        @printf("%-6s %4d %6s %6s %4d | %7.3f %7.3f %7.3f %7.3f %7.2f%%\n",
            row.mode, row.w, isnan(row.k) ? "-" : string(row.k),
            isnan(row.q) ? "-" : string(row.q), row.dil,
            row.pixiou, row.recall, row.precision, row.lineiou, 100row.bfrac)
    end
end

open(joinpath(SCREENS, "e5_bright_calibration.txt"), "w") do io
    println(io, "E5 bright-line detector calibration vs DR17 bright mask (pass 2) — $(Dates.now())")
    println(io, "objective: pixel IoU; scale: per-fiber global MAD; see script header for rationale")
    for (name, set) in (("APO", apo), ("LCO", lco))
        haskey(results, name) || continue
        for (rank, row) in enumerate(first(results[name], 5))
            @printf(io, "%s #%d: mode=%s window=%d k=%s q=%s dil=%d -> pixIoU=%.3f recall=%.3f prec=%.3f lineIoU=%.3f bright=%.2f%%\n",
                name, rank, row.mode, row.w, string(row.k), string(row.q), row.dil,
                row.pixiou, row.recall, row.precision, row.lineiou, 100row.bfrac)
        end
        best = results[name][1]
        println(io, "\n$name per-fiber detail at best config:")
        for d in set
            msk = config_mask(d, best.w, best.dil, best.mode, best.k, best.q)
            st = line_overlap_stats(d.dr17_bright, msk, d.valid)
            @printf(io, "  f%03d: DR17 lines=%d recovered=%d(%.0f%%) new=%d spurious=%d pixIoU=%.3f bright=%.2f%% (DR17 %.2f%%)\n",
                d.adjfib, st.n_ref_lines, st.recovered, 100st.recall, st.n_new_lines,
                st.spurious, st.pixel_iou, 100count(msk .& d.sub) / count(d.sub),
                100st.ref_pix / count(d.valid))
            rr = mask_runs(d.dr17_bright .& d.valid)
            low = [i for i in eachindex(st.line_ious) if st.line_ious[i] < 0.2]
            if !isempty(low)
                @printf(io, "      DR17 lines poorly matched (lambda, npix, IoU): %s\n",
                    string([(round(wavetarg[first(rr[i])], digits=1), length(rr[i]),
                        round(st.line_ious[i], digits=2)) for i in first(low, 10)]))
            end
        end
    end
end
println("\nwrote $(joinpath(SCREENS, "e5_bright_calibration.txt"))")
