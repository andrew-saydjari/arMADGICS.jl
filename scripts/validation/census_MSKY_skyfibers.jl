# M-SKY impact census: how many real 2026_05_01 object exposures carry sky fibers
# that (pre-guard) poison or (post-guard) get excluded from the exposure-level sky
# prior (getSkyRough/combine_sky_fibers)?
#
# Samples object-exposure ar1Dunical files across mjds/telescopes, applies the EXACT
# sky-fiber selection the sky path uses (get_fibTargDict + fibtype[1:3]=="sky"), and
# classifies each exposure:
#   - pre-guard poisoning: a sky fiber with non-finite flux that PASSES the pre-guard
#     scale z-cut -> its NaNs/Infs enter VLocSkyLines (counted; also counted restricted
#     to chebmsk pixels, which are the ones that can reach a target's simplemsk)
#   - all-NaN (A4-shape) sky fibers: silently dropped pre-guard by the z-cut accident
#     (NaN scale), excluded + flagged post-guard
#   - post-guard exclusions / low-sky-count fallback (skyBit)
#
# Usage:
#   julia --project=. scripts/validation/census_MSKY_skyfibers.jl <apredDir> <almanacFile> <chebmskFile> <out_h5> [nexp_target] [seed]
#
# Read-only on all inputs; writes only <out_h5>.

using LinearAlgebra, Statistics, StatsBase, HDF5, JLD2, FITSIO, Random
using EllipsisNotation, ShiftedArrays, SparseArrays, Interpolations, LowRankOps

src_dir = joinpath(dirname(@__DIR__), "..") # repo root
include(joinpath(src_dir, "src/utils.jl"))
include(joinpath(src_dir, "src/fileNameHandling.jl"))
include(joinpath(src_dir, "src/ingest.jl"))
logUniWaveAPOGEE = 10 .^ range((start = 4.179 - 125 * 6.0e-6), step=6.0e-6, length=8575 + 125)

length(ARGS) >= 4 || error("usage: census_MSKY_skyfibers.jl <apredDir> <almanacFile> <chebmskFile> <out_h5> [nexp_target] [seed]")
apredDir = abspath(ARGS[1])
almanacFile = abspath(ARGS[2])
chebmskFile = abspath(ARGS[3])
out_h5 = abspath(ARGS[4])
nexp_target = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 350
seed = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 20260902
max_per_mjd_tele = 6

chebmsk = convert.(Bool, h5read(chebmskFile, "chebmsk_exp"))

rng = MersenneTwister(seed)
mjds = sort(readdir(apredDir))
mjds_shuffled = shuffle(rng, mjds)

# collect the sample of (tele, mjd, expnum, path)
sample = NamedTuple{(:tele, :mjd, :expnum, :path),Tuple{String,String,Int,String}}[]
nscanned_mjds = 0
for mjd in mjds_shuffled
    global nscanned_mjds
    length(sample) >= nexp_target && break
    dir = joinpath(apredDir, mjd)
    isdir(dir) || continue
    nscanned_mjds += 1
    files = filter(x -> startswith(x, "ar1Dunical_") && endswith(x, "_object.h5"), readdir(dir))
    isempty(files) && continue
    for tele in ("apo", "lco")
        tf = filter(x -> occursin("_$(tele)_", x), files)
        isempty(tf) && continue
        pick = shuffle(rng, tf)[1:min(max_per_mjd_tele, length(tf))]
        for fname in pick
            expnum = parse(Int, split(fname, "_")[4])
            push!(sample, (tele=tele, mjd=mjd, expnum=expnum, path=joinpath(dir, fname)))
        end
    end
end
println("sampled $(length(sample)) object exposures from $nscanned_mjds mjd dirs (target $nexp_target, seed $seed)")

n = length(sample)
tele_v = [s.tele for s in sample]
mjd_v = [parse(Int, s.mjd) for s in sample]
expnum_v = [s.expnum for s in sample]
n_skycand = fill(-1, n)          # candidate sky fibers in the config
n_allnan = zeros(Int, n)         # A4-shape: all-NaN/zero flux
n_partial_nonfinite = zeros(Int, n) # non-finite somewhere but not all-NaN
n_excluded_postguard = zeros(Int, n)
n_poisoning_preguard = zeros(Int, n) # non-finite fibers that PASS the pre-guard z-cut
npix_poisoned = zeros(Int, n)        # V pixels poisoned pre-guard (any col non-finite, incremental over healthy)
npix_poisoned_cheb = zeros(Int, n)   # ... restricted to chebmsk
nSkyFibers_post = fill(-1, n)
skyBit_v = fill(-1, n)
status_v = fill("", n)

falm = h5open(almanacFile)
t_start = time()
for (i, s) in enumerate(sample)
    try
        fibtargDict, _ = get_fibTargDict(falm, s.tele, s.mjd, s.expnum)
        fibtypelist = map(x -> fibtargDict[x], 1:300)
        skyfibIndxs = findall(map(x -> x[1:3] == "sky", fibtypelist))
        n_skycand[i] = length(skyfibIndxs)
        if isempty(skyfibIndxs)
            skyBit_v[i] = SKY_NO_FIBERS_BIT | SKY_TOO_FEW_FIBERS_BIT
            nSkyFibers_post[i] = 0
            status_v[i] = "nosky"
            continue
        end
        # jldopen, exactly like getSkyRough (HDF5.jl lacks fancy column indexing)
        fh = jldopen(s.path)
        skyspec = fh["flux_1d"][:, skyfibIndxs]
        skyivar = fh["ivar_1d"][:, skyfibIndxs]
        skymsk = fh["mask_1d"][:, skyfibIndxs]
        close(fh)
        ncand = length(skyfibIndxs)
        bits = [validate_sky_fiber(view(skyspec, :, j), view(skyivar, :, j), view(skymsk, :, j)) for j in 1:ncand]
        allnan = [(b & SKYFIB_ALLNANZERO_BIT) != 0 for b in bits]
        nonfin = [(b & SKYFIB_NONFINITE_BIT) != 0 for b in bits]
        n_allnan[i] = count(allnan)
        n_partial_nonfinite[i] = count(nonfin .& .!allnan)
        n_excluded_postguard[i] = count(bits .!= 0)

        # pre-guard inclusion: the original z-cut over ALL candidates
        skyScale = dropdims(nanzeromedian(skyspec, 1), dims=1)
        skyMed = nanzeromedian(skyScale)
        skyIQR = nanzeroiqr(skyScale)
        skyZ = (skyScale .- skyMed) ./ skyIQR
        mskSky_pre = (abs.(skyZ) .< 10)
        poisoners = findall(nonfin .& mskSky_pre)
        n_poisoning_preguard[i] = length(poisoners)
        if !isempty(poisoners)
            badpix = vec(any(.!isfinite.(skyspec[:, poisoners]), dims=2))
            npix_poisoned[i] = count(badpix)
            npix_poisoned_cheb[i] = count(badpix .& chebmsk)
        end

        # post-guard summary
        nSky, _, _, _, skyBit, _, _ = combine_sky_fibers(skyspec, skyivar, skymsk)
        nSkyFibers_post[i] = nSky
        skyBit_v[i] = skyBit
        status_v[i] = "ok"
    catch e
        status_v[i] = "ERR: " * sprint(showerror, e)[1:min(end, 200)]
    end
    if i % 50 == 0
        println("  $i / $n  ($(round(time()-t_start, digits=1))s)")
        flush(stdout)
    end
end
close(falm)

rm(out_h5; force=true)
mkpath(dirname(out_h5))
h5open(out_h5, "w") do fh
    fh["tele"] = tele_v
    fh["mjd"] = mjd_v
    fh["expnum"] = expnum_v
    fh["n_skycand"] = n_skycand
    fh["n_allnan"] = n_allnan
    fh["n_partial_nonfinite"] = n_partial_nonfinite
    fh["n_excluded_postguard"] = n_excluded_postguard
    fh["n_poisoning_preguard"] = n_poisoning_preguard
    fh["npix_poisoned"] = npix_poisoned
    fh["npix_poisoned_cheb"] = npix_poisoned_cheb
    fh["nSkyFibers_post"] = nSkyFibers_post
    fh["skyBit"] = skyBit_v
    fh["status"] = status_v
    attrs(fh)["apredDir"] = apredDir
    attrs(fh)["seed"] = seed
end

# summary with binomial (Wald) errors
binerr(k, m) = m == 0 ? NaN : sqrt(max(k / m * (1 - k / m), 1 / m) / m)
function report(sel, label)
    m = count(sel)
    m == 0 && return println("$label: no exposures")
    ok = sel .& (status_v .== "ok") .| (sel .& (status_v .== "nosky"))
    mok = count(ok)
    poisoned = count(ok .& (npix_poisoned_cheb .> 0))
    poisonedany = count(ok .& (npix_poisoned .> 0))
    flagged = count(ok .& (skyBit_v .> 0))
    excl = count(ok .& (n_excluded_postguard .> 0))
    a4 = count(ok .& (n_allnan .> 0))
    lowsky = count(ok .& ((skyBit_v .& SKY_TOO_FEW_FIBERS_BIT) .> 0))
    nosky = count(ok .& (n_skycand .== 0))
    println("$label: n=$mok (of $m sampled)")
    println("  pre-guard POISONED (nonfinite in-chebmsk V pixels): $poisoned  frac=$(round(poisoned/mok, sigdigits=3)) +/- $(round(binerr(poisoned, mok), sigdigits=2))")
    println("  pre-guard poisoned incl. out-of-chebmsk-only:        $poisonedany  frac=$(round(poisonedany/mok, sigdigits=3))")
    println("  >=1 sky fiber excluded post-guard:                   $excl  frac=$(round(excl/mok, sigdigits=3)) +/- $(round(binerr(excl, mok), sigdigits=2))")
    println("  >=1 all-NaN (A4) sky fiber:                          $a4  frac=$(round(a4/mok, sigdigits=3)) +/- $(round(binerr(a4, mok), sigdigits=2))")
    println("  low-sky fallback (<$(SKY_MIN_FIBERS) survive):       $lowsky")
    println("  zero sky candidates in config:                       $nosky")
end
println("\n=== census summary ===")
report(trues(n), "ALL")
report(tele_v .== "apo", "APO")
report(tele_v .== "lco", "LCO")
nerr = count(startswith.(status_v, "ERR"))
println("errors: $nerr")
for i in findall(startswith.(status_v, "ERR"))[1:min(end, 10)]
    println("  $(tele_v[i]) $(mjd_v[i]) $(expnum_v[i]): $(status_v[i])")
end
println("results -> $out_h5")
