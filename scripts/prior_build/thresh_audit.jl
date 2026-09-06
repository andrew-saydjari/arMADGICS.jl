## E5 bright/faint threshold audit (integration-agent finding #35).
# Reproduces build_skyLines' median_sky EXACTLY for given fibers, reports the
# flux-unit distribution, what the DR17-inherited thresholds (apo 2000 /
# lco 645) flag, and whether the exported submsk equals the no-split mask.
# Usage: julia --project=<arM-E5b> thresh_audit.jl <adjfib> [adjfib ...]
using HDF5, Statistics, Printf

# verbatim copies of the two builder helpers (arM src/utils.jl:11-37,
# scripts/prior_build/prior_utils.jl:178-190) — copied rather than included so
# the audit does not drag in the Distributed/Suppressor deps of prior_utils.jl
isnanorzeroorinf(x) = isnan(x) || iszero(x) || isinf(x)
nanzeromedian(x) = all(isnanorzeroorinf, x) ? NaN : median(filter(!isnanorzeroorinf, x))
nanzeromedian(x, y) = mapslices(nanzeromedian, x, dims=y)
function expand_msk(msk; rad=1)
    lmsk = length(msk)
    msknew = ones(Bool, lmsk)
    for i in 1:lmsk
        lindx = maximum((1, i - rad))
        rindx = minimum((i + rad, lmsk))
        if any(.!view(msk, lindx:rindx))
            msknew[i] = false
        end
    end
    return msknew
end

const OUT = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1"
const CHIPGAP = "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5"
default_thresh(adjfib) = (adjfib <= 300) ? 2000 : 645

read_sample(prefix, adjfib) = h5read(joinpath(OUT, "samples",
    prefix * "_" * lpad(adjfib, 3, "0") * ".h5"), prefix)

for adjfib in parse.(Int, ARGS)
    n = lpad(adjfib, 3, "0")
    tele = adjfib <= 300 ? "apo" : "lco"
    skyline = read_sample("skyline", adjfib)
    skymsk = convert.(Float64, read_sample("skymsk", adjfib))
    chipgapmsk = h5open(CHIPGAP, "r") do f
        adjfib > 300 ? read(f["lco"]) : read(f["apo"])
    end

    # ---- exactly the builder's pre-split block ----
    specsum = dropdims(sum(skyline, dims=1), dims=1)
    obscnt = dropdims(sum(skymsk, dims=2), dims=2)
    submsk = (obscnt .>= 10) .& chipgapmsk
    Vred = skyline[submsk, specsum .> 0]
    skymsked = skymsk[submsk, specsum .> 0]
    Vred .*= skymsked
    median_sky = dropdims(nanzeromedian(Vred, 2), dims=2)

    thr = default_thresh(adjfib)
    mskflux = .!expand_msk(median_sky .< thr, rad=4)   # true = BRIGHT
    nbright = count(mskflux)
    npix = length(median_sky)

    qs = quantile(filter(isfinite, median_sky), [0.5, 0.9, 0.99, 0.999, 1.0])
    @printf("f%s (%s) nsamp=%d submsk_pix=%d  median_sky: med=%.1f p90=%.1f p99=%.1f p99.9=%.1f MAX=%.1f\n",
        n, tele, size(Vred, 2), npix, qs...)
    @printf("    DR17-inherited thresh=%d -> BRIGHT pixels %d (%.4f%%)%s\n",
        thr, nbright, 100nbright / npix, nbright == 0 ? "   <-- SPLIT IS A NO-OP" : "")

    # candidate thresholds: fraction flagged bright (post expand rad=4)
    print("    candidates: ")
    for c in [20, 25, 30, 35, 40, 45, 50, 75, 100, 200, 500]
        nb = count(.!expand_msk(median_sky .< c, rad=4))
        @printf("%d->%.2f%%  ", c, 100nb / npix)
    end
    println()
    # threshold achieving target bright fractions (DR17 deployed ~8.35% apo / 8.15% lco)
    for target in [0.0835, 0.05, 0.02]
        lo, hi = 1.0, 3000.0
        for _ in 1:40
            mid = sqrt(lo * hi)
            fr = count(.!expand_msk(median_sky .< mid, rad=4)) / npix
            fr > target ? (lo = mid) : (hi = mid)
        end
        @printf("    target bright %.2f%% -> threshold %.1f\n", 100target, sqrt(lo * hi))
    end
    h5open(joinpath(OUT, "screens", "median_sky_$n.h5"), "w") do fh
        fh["median_sky"] = median_sky
        fh["submsk"] = Vector{Bool}(submsk)
    end

    # exported submsk in the built faint prior vs the no-split mask
    fp = joinpath(OUT, "built", "APOGEE_skyline_faint_svd_120_f$n.h5")
    if isfile(fp)
        sm = Bool.(h5read(fp, "submsk"))
        nosplit = submsk
        @printf("    built submsk: %d pix; identical to (obscnt>=10 & chipgap) no-split mask: %s\n",
            count(sm), count(sm .!= nosplit) == 0)
    else
        println("    built faint prior not present yet")
    end
    flush(stdout)
end
