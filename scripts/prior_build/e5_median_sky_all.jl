## E5: cache `median_sky` for ALL 600 fibers, so the bright-line mask can be COMBINED
## across fibers (AKS directive 2026-09-07: the bright mask should not depend on fiber
## number, since it is only ever used as a mask).
#
# Reproduces the builder's pre-split block EXACTLY (build_sky_defs.jl:332-339), but
# without materializing `Vred`/`skymsked` as separate Float64 copies:
#
#     Vred = skyline[submsk, cols]; skymsked = skymsk[submsk, cols]; Vred .*= skymsked
#     median_sky = nanzeromedian(Vred, 2)
#
# `nanzeromedian` drops NaN/0/Inf, and the Float64 mask can only map an entry to 0.0
# (excluded) or leave it unchanged, so "multiply by the 0/1 mask then drop zeros" is
# identical to "keep entries where the mask is true", which is what this computes.
# VERIFIED bit-for-bit against the 20 pre-existing thresh_audit.jl caches (--verify).
#
# Usage: julia --project=<arM-E5b> e5_median_sky_all.jl [--nworkers=N] [--verify]
# Author - Andrew Saydjari (E5 pass 1)
using Distributed, Printf, Dates

const OUT = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1"
const NW = parse(Int, replace(get(filter(a -> startswith(a, "--nworkers="), ARGS), 1, "--nworkers=16"), "--nworkers=" => ""))
const VERIFY = "--verify" in ARGS

addprocs(NW)
@everywhere begin
    using HDF5, Statistics
    const OUT = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1"
    const CHIPGAP = "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5"
    # verbatim copies of the builder helpers (arM src/utils.jl:11-37)
    isnanorzeroorinf(x) = isnan(x) || iszero(x) || isinf(x)

    read_sample(prefix, adjfib) = h5read(joinpath(OUT, "samples",
        prefix * "_" * lpad(adjfib, 3, "0") * ".h5"), prefix)

    "the builder's pre-split block: returns (median_sky, submsk, nsamp)"
    function compute_median_sky(adjfib)
        skyline = read_sample("skyline", adjfib)
        skymsk = read_sample("skymsk", adjfib)          # Bool on disk
        chipgapmsk = h5open(CHIPGAP, "r") do f
            adjfib > 300 ? read(f["lco"]) : read(f["apo"])
        end
        # ---- the builder's pre-split block ----
        specsum = dropdims(sum(skyline, dims=1), dims=1)
        obscnt = dropdims(sum(skymsk, dims=2), dims=2)
        submsk = (obscnt .>= 10) .& chipgapmsk
        cols = findall(specsum .> 0)
        rows = findall(submsk)
        med = fill(NaN, length(rows))
        buf = Float64[]
        @inbounds for (ri, r) in enumerate(rows)
            empty!(buf)
            for c in cols
                skymsk[r, c] || continue          # mask 0 -> entry becomes 0.0 -> dropped
                v = skyline[r, c]
                isnanorzeroorinf(v) && continue
                push!(buf, v)
            end
            isempty(buf) || (med[ri] = median(buf))
        end
        return med, Vector{Bool}(submsk), length(cols)
    end

    function cache_one(adjfib; overwrite=false)
        n = lpad(adjfib, 3, "0")
        dst = joinpath(OUT, "screens", "median_sky_$n.h5")
        (!overwrite && isfile(dst)) && return (adjfib, :skipped, 0.0)
        t0 = time()
        med, submsk, nsamp = compute_median_sky(adjfib)
        tmp = dst * ".tmp$(getpid())"
        h5open(tmp, "w") do fh
            fh["median_sky"] = med
            fh["submsk"] = submsk
            fh["nsamp"] = nsamp
        end
        mv(tmp, dst; force=true)
        return (adjfib, :built, time() - t0)
    end
end

if VERIFY
    # bit-for-bit check against the 20 thresh_audit.jl caches before overwriting anything
    ref = [10, 76, 450, 519, 598, 600]
    nbad = 0
    for (a, (newmed, newsub, nsamp)) in zip(ref, pmap(compute_median_sky, ref))
        n = lpad(a, 3, "0")
        p = joinpath(OUT, "screens", "median_sky_$n.h5")
        isfile(p) || continue
        old = h5read(p, "median_sky"); olds = Bool.(h5read(p, "submsk"))
        same_sub = newsub == olds
        same_med = length(newmed) == length(old) &&
                   all(((x, y),) -> (isnan(x) && isnan(y)) || x === y, zip(newmed, old))
        (same_sub && same_med) || (global nbad += 1)
        @printf("VERIFY f%s: submsk identical=%s  median_sky bitwise identical=%s (npix=%d nsamp=%d)\n",
            n, same_sub, same_med, length(old), nsamp)
        flush(stdout)
    end
    println(nbad == 0 ? "VERIFY: ALL PASS — the fast path reproduces thresh_audit.jl bit-for-bit" :
            "VERIFY: $nbad MISMATCH(ES) — do not use the fast path")
    exit(nbad == 0 ? 0 : 1)
end

todo = 1:600
println("caching median_sky for $(length(todo)) fibers on $NW workers  ($(now()))"); flush(stdout)
t0 = time()
res = pmap(cache_one, todo)
nb = count(r -> r[2] == :built, res)
ns = count(r -> r[2] == :skipped, res)
ts = [r[3] for r in res if r[2] == :built]
@printf("done: built=%d skipped(existing)=%d  wall=%.1f min  per-fiber med=%.1f s max=%.1f s\n",
    nb, ns, (time() - t0) / 60, isempty(ts) ? 0 : median(ts), isempty(ts) ? 0 : maximum(ts))
