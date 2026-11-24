# 02_solve_offsets.jl
#
# CHAIN POSITION: fiber offsets step 2 of 5 (see ../README.md).
# INPUTS:  ARGS[1] WORKDIR (reads build.h5)
# OUTPUTS: WORKDIR/offsets_raw.csv, WORKDIR/solve.h5
# SOURCE:  work/2026_09_16/fiber_offsets/s2_solve.jl (0.4.0). Logic unchanged; hard-coded paths became
#          ARGS (WORKDIR = the fiber-offset working directory shared by steps 1-5).
#
# s2_solve.jl -- the per-(telescope, fiber, era) RV zero-point solve.
#
# ------------------------------------------------------------------------------------------
# THE ESTIMATOR, and why it is not the diluted one
#
# Model (AJ 169, 167, arXiv:2408.07126, Sec. 4): every exposure m of star i taken on
# block j carries
#                        y_m = mu_i + f_j + noise
# with mu_i the star's true velocity and f_j the zero point of block j.  Blocks are indexed
#   b = adjfiberindx + 600*(era-1)  in 1..1200,
# so ONE vector f absorbs the per-fiber term, the APO/LCO telescope term and the plate->FPS
# era step simultaneously.  That is what was asked for.
#
# The naive estimator -- take each star's median, difference, average per fiber -- is DILUTED.
# 65 % of stars are seen on ONE block only; their median already contains that block's own
# offset, so they contribute a residual of exactly zero and drag the per-fiber average toward
# zero (2026_09_13: naive correlates r = 0.98 with the correct answer but at slope 0.34).
#
# This script instead solves the FULL joint least-squares problem for mu and f together, and
# then profiles mu out analytically.  With w_m the weight, W_ij = sum of star i's weights on
# block j, T_i = sum_j W_ij, Y_ij = sum of w_m y_m on block j and Z_i = sum_j Y_ij,
#
#     M_jk =  delta_jk * sum_i W_ij  -  sum_i W_ij W_ik / T_i
#     b_j  =  sum_i ( Y_ij - W_ij Z_i / T_i )          ,        M f = b
#
# M is the weighted Laplacian of the star-sharing graph over blocks.  A star confined to one
# block gives W_ij = T_i and its contribution cancels EXACTLY (W_j - W_j^2/T = 0): it neither
# informs nor dilutes.  This is the undiluted estimator, and it is algebraically identical to
# the paper's ubercalibration once its alternating iteration is run to convergence.
#
# GAUGE.  M's null space is the constant vector on each connected component (the absolute
# zero point is unconstrained -- the paper's "one eigenvalue 6-10 orders of magnitude below
# the next").  We solve with a pseudo-inverse, which returns the minimum-norm solution, i.e.
# sum_j f_j = 0 over each component.  On the giant component that is the paper's own gauge
# ("enforcing the mean fiber offset to be zero") and it means the correction introduces NO
# net velocity shift: v_rad_corr is on the same absolute scale as v_rad.
#
# WEIGHTS.  Pass 1 is unweighted (C^-1 = I), exactly as the paper does, to break the
# circularity between the offsets and the error model.  The residuals of pass 1 then define
# an EMPIRICAL single-exposure precision curve sigma_emp(S/N) per (telescope, era) -- robust
# MAD scatter in 24 log-S/N bins, monotone-free, floored at 5 m/s -- and passes 2..n use
# w = 1/sigma_emp^2.  e_v_rad is NOT used: it has no systematic floor at high S/N, so it
# would over-weight exactly the bright stars whose achieved scatter turns back up.
#
# SAMPLE.  Release default mask, category == science, snr >= 15, and an RV-STABILITY cut on
# the star (robust scatter <= 1.0 km/s and max deviation from its own median <= 2.0 km/s) to
# keep binaries and pulsators out of the fit.  1.0 km/s is >= 2.5x the largest possible
# offset DIFFERENCE (+-200 m/s offsets -> <= 400 m/s spread), so the cut cannot clip the
# signal it is measuring; s6 reports the sensitivity to moving it.  Residuals are then
# sigma-clipped at 4 sigma_emp each iteration.
#
# TIME DEPENDENCE IS DELIBERATELY LEFT OUT.  2026_09_13 measured a median 21.6 m/s epoch
# drift.  Folding a per-epoch term in would make the table time-resolved and much harder to
# validate, and the epoch term is ~2x smaller than the per-fiber rms it would sit on top of.
# This first table is time-independent within an era.  Stated as a limitation, not a claim.
#
# ERRORS.  Two, both reported: the analytic inverse curvature diag(pinv(M))^1/2 (the paper's
# Fig. 4 bottom-right), which assumes independent exposures and therefore UNDERSTATES, and a
# BOOTSTRAP OVER STARS (64 replicates, stars drawn with replacement), which is the shipped
# `err` because it carries the within-star and within-visit correlation correctly.
# ------------------------------------------------------------------------------------------
using HDF5, Statistics, Printf, LinearAlgebra, Random, Base.Threads, DataFrames, CSV

const DIR = ARGS[1]
const NB  = 1200
const SNRMIN   = 15.0
const SSTAB    = 1.0      # km/s, star robust-scatter cut
const SMAXDEV  = 2.0      # km/s, star max |y - median| cut
const CLIP     = 4.0      # sigma_emp
const NITER    = 4
const NBOOT    = 64
const SEED     = 20260916

say(s) = (println(s); flush(stdout))
madsig(v) = length(v) < 2 ? NaN : (length(v) == 2 ? abs(v[1]-v[2])/sqrt(2) : 1.4826*median(abs.(v .- median(v))))

# ---------------------------------------------------------------- normal equations
"""Accumulate M, b from per-star block aggregates.  `mult` is the star multiplicity
(1 everywhere for the point estimate; Poisson/bootstrap counts for a replicate)."""
function accumulate!(M, b, soff, sblk, sW, sY, mult)
    fill!(M, 0.0); fill!(b, 0.0)
    ns = length(soff) - 1
    for i in 1:ns
        c = mult[i]; c == 0 && continue
        lo = soff[i]; hi = soff[i+1] - 1
        nbk = hi - lo + 1
        nbk >= 2 || continue
        T = 0.0; Z = 0.0
        @inbounds for k in lo:hi; T += sW[k]; Z += sY[k]; end
        T > 0 || continue
        @inbounds for k in lo:hi
            j = sblk[k]
            b[j] += c * (sY[k] - sW[k]*Z/T)
            M[j,j] += c * sW[k]
            for l in lo:hi
                M[j, sblk[l]] -= c * sW[k]*sW[l]/T
            end
        end
    end
    return nothing
end

"""Minimum-norm solve of M f = b, per connected component, dropping null directions."""
function solve_mn(M, b, live; want_curv::Bool = true)
    f = fill(NaN, NB); ferr = fill(NaN, NB)
    idx = findall(live)
    isempty(idx) && return f, ferr, Int[]
    # connected components of the off-diagonal sparsity
    lab = zeros(Int, NB); nc = 0
    for s in idx
        lab[s] == 0 || continue
        nc += 1; stack = [s]; lab[s] = nc
        while !isempty(stack)
            u = pop!(stack)
            for v in idx
                (lab[v] == 0 && v != u && abs(M[u,v]) > 0) || continue
                lab[v] = nc; push!(stack, v)
            end
        end
    end
    comps = [findall(==(c), lab) for c in 1:nc]
    for cc in comps
        length(cc) >= 2 || continue
        Ms = Symmetric(M[cc, cc]); bs = b[cc]
        F = eigen(Ms)
        λ = F.values; V = F.vectors
        tol = maximum(λ) * 1e-9
        keep = λ .> tol
        xs = V[:, keep] * ((V[:, keep]' * bs) ./ λ[keep])
        xs .-= mean(xs)                     # min-norm gauge within the component
        f[cc] = xs
        if want_curv   # analytic inverse curvature (paper Fig. 4, bottom right)
            Vk = V[:, keep]; iv = 1.0 ./ λ[keep]
            ferr[cc] = sqrt.(max.([sum(iv .* (@view Vk[r, :]).^2) for r in 1:length(cc)], 0.0))
        end
    end
    return f, ferr, lab
end

function main()
    t0 = time()
    say("loading working table ...")
    h = h5open(joinpath(DIR, "build.h5"), "r")
    sid = read(h["exp/sdss_id"]); blk = Int.(read(h["exp/block"]))
    y   = read(h["exp/v_rad"]);   snr = read(h["exp/snr"])
    cat = read(h["exp/category"]); tele = read(h["exp/tele"]); era = read(h["exp/era"])
    close(h)
    n0 = length(sid)
    say("  usable exposures = $n0")

    sel = (cat .== 0x01) .& (snr .>= SNRMIN) .& (sid .> 0)
    say(@sprintf("  science & snr >= %.0f & sdss_id > 0 : %d (%.2f %%)", SNRMIN, count(sel), 100count(sel)/n0))
    E = findall(sel)
    sidE = sid[E]; blkE = blk[E]; yE = y[E]; snrE = snr[E]
    teleE = tele[E]; eraE = era[E]
    ne = length(E)

    # ------------------------------------------------------------ star grouping
    say("grouping by star ...")
    ord = sortperm(sidE)
    sidE = sidE[ord]; blkE = blkE[ord]; yE = yE[ord]; snrE = snrE[ord]
    teleE = teleE[ord]; eraE = eraE[ord]
    gb = Int32[1]
    for k in 2:ne; sidE[k] == sidE[k-1] || push!(gb, Int32(k)); end
    push!(gb, Int32(ne+1))
    NS = length(gb) - 1
    say("  distinct stars = $NS")

    # ------------------------------------------------------------ star stability cut
    say("RV-stability cut on stars ...")
    keepobs = falses(ne)
    nstab = 0; nmulti = 0
    starscat = fill(NaN, NS); starnblk = zeros(Int16, NS)
    buf = Float64[]
    for i in 1:NS
        lo = gb[i]; hi = gb[i+1]-1
        n = hi - lo + 1
        resize!(buf, n); @inbounds for k in 1:n; buf[k] = yE[lo+k-1]; end
        nbk = length(unique(@view blkE[lo:hi]))
        starnblk[i] = Int16(nbk)
        if n == 1
            starscat[i] = 0.0
        else
            m = median(buf); s = madsig(buf)
            starscat[i] = s
            (isfinite(s) && s <= SSTAB && maximum(abs.(buf .- m)) <= SMAXDEV) || continue
        end
        nstab += 1
        nbk >= 2 && (nmulti += 1)
        @inbounds for k in lo:hi; keepobs[k] = true; end
    end
    say(@sprintf("  stars passing stability: %d / %d (%.2f %%);  of those %d see >= 2 blocks (%.2f %%)",
        nstab, NS, 100nstab/NS, nmulti, 100nmulti/nstab))
    say(@sprintf("  exposures kept: %d / %d", count(keepobs), ne))

    # block census, and the adjacency matrix (paper Fig. 4 bottom-left)
    adj = zeros(Int32, NB, NB)
    nexp_blk = zeros(Int64, NB); nstar_blk = zeros(Int32, NB); nstar_info = zeros(Int32, NB)
    for i in 1:NS
        lo = gb[i]; hi = gb[i+1]-1
        keepobs[lo] || continue
        bs = unique(@view blkE[lo:hi])
        for j in bs; nstar_blk[j] += 1; end
        if length(bs) >= 2
            for j in bs
                nstar_info[j] += 1
                for k in bs; j == k || (adj[j,k] += 1); end
            end
        end
        for k in lo:hi; nexp_blk[blkE[k]] += 1; end
    end
    say(@sprintf("  blocks with any data: %d / %d ; with >= 1 informing star: %d",
        count(>(0), nexp_blk), NB, count(>(0), nstar_info)))

    # ------------------------------------------------------------ iterate
    w    = ones(Float64, ne)
    good = copy(keepobs)
    f    = zeros(Float64, NB)
    ferr_an = fill(NaN, NB)
    sigmod = Dict{Tuple{Int8,Int8}, Tuple{Vector{Float64},Vector{Float64}}}()  # (tele,era) -> (logsnr edges, sigma)
    sblk = Int32[]; sW = Float64[]; sY = Float64[]; soff = Int32[1]
    M = zeros(NB, NB); bvec = zeros(NB)
    live = falses(NB)
    resid = fill(NaN, ne)

    for it in 1:NITER
        # --- per-star block aggregates over the currently-good obs
        empty!(sblk); empty!(sW); empty!(sY); resize!(soff, 1); soff[1] = 1
        for i in 1:NS
            lo = gb[i]; hi = gb[i+1]-1
            base = length(sblk)
            for k in lo:hi
                good[k] || continue
                j = blkE[k]
                p = 0
                @inbounds for q in (base+1):length(sblk); sblk[q] == j && (p = q; break); end
                if p == 0
                    push!(sblk, Int32(j)); push!(sW, w[k]); push!(sY, w[k]*yE[k])
                else
                    sW[p] += w[k]; sY[p] += w[k]*yE[k]
                end
            end
            length(sblk) - base >= 2 || (resize!(sblk, base); resize!(sW, base); resize!(sY, base))
            push!(soff, Int32(length(sblk)+1))
        end
        nsused = count(i -> soff[i+1] > soff[i], 1:NS)
        accumulate!(M, bvec, soff, sblk, sW, sY, ones(Int32, NS))
        live .= [M[j,j] > 0 for j in 1:NB]
        f, ferr_an, lab = solve_mn(M, bvec, live)

        # --- residuals with the star mean profiled out
        fill!(resid, NaN)
        for i in 1:NS
            lo = gb[i]; hi = gb[i+1]-1
            sw = 0.0; sz = 0.0
            for k in lo:hi
                (good[k] && isfinite(f[blkE[k]])) || continue
                sw += w[k]; sz += w[k]*(yE[k] - f[blkE[k]])
            end
            sw > 0 || continue
            mu = sz/sw
            for k in lo:hi
                (keepobs[k] && isfinite(f[blkE[k]])) || continue
                resid[k] = yE[k] - f[blkE[k]] - mu
            end
        end

        # --- empirical error model from the residuals, per (tele, era)
        empty!(sigmod)
        edges = collect(range(log10(SNRMIN), 3.0, length = 25))
        NQ = length(edges) - 1
        bins = [Float64[] for _ in 1:2, _ in 1:2, _ in 1:NQ]   # single pass, capped
        for k in 1:ne
            isfinite(resid[k]) || continue
            q = clamp(searchsortedlast(edges, log10(snrE[k])), 1, NQ)
            v = bins[teleE[k], eraE[k], q]
            length(v) < 200_000 && push!(v, resid[k])
        end
        for tl in Int8.(1:2), er in Int8.(1:2)
            sig = fill(NaN, NQ)
            for q in 1:NQ
                vv = bins[tl, er, q]
                length(vv) >= 200 && (sig[q] = 1.4826*median(abs.(vv .- median(vv))))
            end
            fin = findall(isfinite, sig)
            for q in 1:NQ
                isfinite(sig[q]) && continue
                sig[q] = isempty(fin) ? 0.1 : sig[fin[argmin(abs.(fin .- q))]]
            end
            sigmod[(tl,er)] = (edges, max.(sig, 0.005))
        end
        mm = [sigmod[(tl,er)][2] for tl in Int8.(1:2), er in Int8.(1:2)]
        say(@sprintf("  it %d: stars informing = %d ; good obs = %d ; sigma_emp median (m/s) apo-plate %.1f apo-fps %.1f lco-plate %.1f lco-fps %.1f",
            it, nsused, count(good), 1000median(mm[1,1]), 1000median(mm[1,2]), 1000median(mm[2,1]), 1000median(mm[2,2])))

        # --- weights and clipping for the next iteration
        nclip = 0
        for k in 1:ne
            keepobs[k] || (good[k] = false; continue)
            edg, sg = sigmod[(teleE[k], eraE[k])]
            q = clamp(searchsortedlast(edg, log10(snrE[k])), 1, length(sg))
            s = sg[q]
            w[k] = 1.0 / s^2
            gk = isfinite(resid[k]) ? abs(resid[k]) <= CLIP*s : true
            gk || (nclip += 1)
            good[k] = gk
        end
        say(@sprintf("       clipped at %.0f sigma: %d obs (%.3f %%)", CLIP, nclip, 100nclip/ne))
    end

    # ------------------------------------------------------------ bootstrap over stars
    say("bootstrap over stars ($NBOOT replicates) ...")
    Random.seed!(SEED)
    fb = fill(NaN, NB, NBOOT)
    Mb = zeros(NB, NB); bb = zeros(NB)
    for r in 1:NBOOT
        mult = zeros(Int32, NS)
        for _ in 1:NS; mult[rand(1:NS)] += Int32(1); end
        accumulate!(Mb, bb, soff, sblk, sW, sY, mult)
        lv = [Mb[j,j] > 0 for j in 1:NB]
        fr, _, _ = solve_mn(Mb, bb, lv; want_curv = false)
        fb[:, r] = fr
        r % 16 == 0 && say("   ... $r / $NBOOT  ($(round(time()-t0,digits=1)) s)")
    end
    ferr = fill(NaN, NB)
    for j in 1:NB
        v = filter(isfinite, @view fb[j, :])
        length(v) >= div(NBOOT, 2) && (ferr[j] = std(v))
    end

    # ------------------------------------------------------------ write
    say("writing the offset table ...")
    rows = DataFrame(block=Int[], tele=String[], fiber=Int[], era=Int[],
                     offset_kms=Float64[], err_kms=Float64[], err_curv_kms=Float64[],
                     n_stars_used=Int[], n_stars_any=Int[], n_exp=Int[], calibratable=Bool[])
    for j in 1:NB
        er = div(j - 1, 600) + 1
        af = j - 600*(er-1)
        tl = af <= 300 ? "apo" : "lco"
        fb_ = af <= 300 ? af : af - 300
        push!(rows, (j, tl, fb_, er, f[j], ferr[j], ferr_an[j],
                     Int(nstar_info[j]), Int(nstar_blk[j]), Int(nexp_blk[j]), false))
    end
    CSV.write(joinpath(DIR, "offsets_raw.csv"), rows)

    h5open(joinpath(DIR, "solve.h5"), "w") do g
        g["offset"] = f; g["err_boot"] = ferr; g["err_curv"] = ferr_an
        g["boot"] = fb; g["adj"] = adj
        g["n_stars_info"] = nstar_info; g["n_stars_any"] = nstar_blk; g["n_exp"] = nexp_blk
        g["star_scatter"] = starscat; g["star_nblk"] = starnblk
        g["resid"] = resid; g["resid_snr"] = snrE; g["resid_tele"] = teleE
        g["resid_era"] = eraE; g["resid_blk"] = Int16.(blkE)
    end
    chmod(joinpath(DIR, "solve.h5"), 0o644); chmod(joinpath(DIR, "offsets_raw.csv"), 0o644)

    for tl in ("apo","lco"), er in 1:2
        m = (rows.tele .== tl) .& (rows.era .== er) .& isfinite.(rows.offset_kms)
        v = 1000 .* rows.offset_kms[m]
        isempty(v) && continue
        say(@sprintf("  %s era %d : n=%3d  mean %+7.1f  rms-about-mean %6.1f  range %+8.1f .. %+8.1f  median err %5.1f m/s",
            tl, er, count(m), mean(v), sqrt(mean((v .- mean(v)).^2)), minimum(v), maximum(v),
            1000median(filter(isfinite, rows.err_kms[m]))))
    end
    say("DONE  $(round(time()-t0,digits=1)) s")
end

main()
