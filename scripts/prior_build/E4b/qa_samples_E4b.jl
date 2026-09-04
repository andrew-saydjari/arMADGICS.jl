# E4b step 3: QA the regenerated LCO starCont samples (adjfib 301-600, C2_LCO=3000
# tfunlist) against the retained 2026_09_03 E4 generation. Follows E4's
# qa_starCont.jl conventions.
#  (a) sweep: 300/300 present, shape (8700,10000), all finite, per-fiber stats.
#  (d) full per-fiber per-column comparison vs the 20260903 LCO samples, checked
#      against the EXACT RNG-stream prediction (E4b_run/rng_predict_final.jdat):
#      - fibers with unchanged tfunlist (113): file must be BIT-IDENTICAL;
#      - affected fibers (187): the set of differing columns must EQUAL the
#        predicted set (changed Tfunindx draws), all other columns bit-identical.
#      This simultaneously proves the leaked columns are GONE (the prediction is
#      built from the leak-free list and asserts no leak entry is drawn).
#  (b) routing audit: deterministic bit-exact regeneration of the first nchk stored
#      columns from the telescope-routed inputs incl. the NEW lcoC2 tfunlist, for
#      spot fibers (one unchanged; 450/595 headline; 519 = collateral-trimmed).
# Run (after sampling completes):
#   cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run && nice -n 10 \
#   julia +1.11.6 --project=/mnt/home/asaydjari/gitcode/worktrees/arM-E4b -p 8 \
#   /mnt/home/asaydjari/gitcode/worktrees/arM-E4b/scripts/prior_build/E4b/qa_samples_E4b.jl
using Distributed
@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using Serialization, Statistics, StatsBase, Printf
end
using HDF5, Dates, Random

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run"
@everywhere const NEW = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_20260904_lcoC2/tell_prior_disk"
@everywhere const OLD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/prior_outputs/starCont_20260903/tell_prior_disk"

pred = deserialize(joinpath(RUND, "rng_predict_final.jdat"))
identset = Set(pred.identfib)                       # fiberindx (1:300)
diffcols = pred.diffcols                            # fiberindx => predicted changed cols

io = open(joinpath(RUND, "QA_SAMPLES_REPORT.txt"), "w")
logmsg(args...) = (println(args...); println(io, args...); flush(stdout); flush(io))
logmsg("E4b LCO sample QA — ", now())

@everywhere function check_fiber(adjfib, predcols_or_nothing)
    fn = joinpath(NEW, "starCont_" * lpad(adjfib, 3, "0") * ".jdat")
    isfile(fn) || return (adjfib=adjfib, ok=false, msg="MISSING")
    X = deserialize(fn)
    size(X) == (8700, 10000) || return (adjfib=adjfib, ok=false, msg="BAD SHAPE $(size(X))")
    all(isfinite, X) || return (adjfib=adjfib, ok=false, msg="NONFINITE")
    smed = median(X, dims=1)
    fibmed = median(smed); fibiqr = iqr(vec(smed))
    Xo = deserialize(joinpath(OLD, "starCont_" * lpad(adjfib, 3, "0") * ".jdat"))
    changed = findall(j -> (@view X[:, j]) != (@view Xo[:, j]), 1:10000)
    expected = predcols_or_nothing === nothing ? Int[] : predcols_or_nothing
    okpred = (changed == expected)
    msg = okpred ? "" :
        @sprintf("PREDICTION MISMATCH: %d changed (expected %d); symdiff %d",
            length(changed), length(expected), length(symdiff(Set(changed), Set(expected))))
    (adjfib=adjfib, ok=okpred, msg=msg, nchanged=length(changed), nexpected=length(expected),
     fibmed=fibmed, fibiqr=fibiqr)
end

args = [(adjfib, (adjfib - 300) in identset ? nothing : diffcols[adjfib - 300]) for adjfib in 301:600]
res = pmap(a -> check_fiber(a...), args)

nfail = count(r -> !r.ok, res)
for r in res
    r.ok || logmsg("FAIL adjfib $(r.adjfib): $(r.msg)")
end
nident = count(r -> r.ok && r.nexpected == 0, res)
nchanged_tot = sum(r -> get(r, :nchanged, 0), res)
meds = [r.fibmed for r in res if r.ok]
logmsg(@sprintf("(a) sweep+compare: %d/300 pass; per-fiber median-of-sample-medians range [%.4g, %.4g]",
    300 - nfail, minimum(meds), maximum(meds)))
logmsg(@sprintf("(d) bit-identical fibers (unchanged tfunlist): %d (predicted %d)", nident, length(pred.identfib)))
logmsg(@sprintf("(d) affected fibers: %d; total changed columns %d — all EXACTLY matching the RNG-stream prediction: %s",
    300 - nident, nchanged_tot, nfail == 0 ? "YES" : "NO"))
logmsg("(d) leak-gone: prediction was generated from the leak-free lcoC2 list with an assertion that rows {2578,3830} are never drawn; prediction match above certifies the stored samples contain no leaked-domeflat draws.")

# ── (b) routing audit: bit-exact regeneration from routed inputs ─────────────────────────
using FITSIO, LowRankOps, EllipsisNotation, ShiftedArrays, JLD2, FileIO
using Interpolations, SparseArrays, AstroTime, Suppressor, DataFrames, ProgressMeter
using SortFilters, BasisFunctions, DustExtinction, DelimitedFiles
using ApogeeReduction: adjFiberIndx2FiberIndx, get_lsf_matrix
proj_path = "/mnt/home/asaydjari/gitcode/worktrees/arM-E4b/"
for f in ["src/utils.jl","src/gridSearch.jl","src/componentAndPosteriors.jl","src/fileNameHandling.jl",
          "src/ingest.jl","src/lowRankPrescription.jl","src/marginalizeEW.jl","src/spectraInterpolation.jl",
          "src/chi2Wrappers.jl","scripts/prior_build/prior_utils.jl"]
    include(proj_path*f)
end
prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
LSF_LCO = prior_dir*"2026_08_31/prior_inputs/lsf_20260427/fpiLSFparams_REGULARIZED_lco_60861.h5"
TFUN_LCO = prior_dir*"2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5"
LST_LCO = prior_dir*"2026_09_04/tfunlists_lcoC2/20260904_lco_tfunlist.jdat"
FRAC_LCO = prior_dir*"2026_04_26/outsamptell_lco.jdat"
nsamp = 10_000
Teff_rng = 4_000:1:10_000; Av_rng = 0:1e-4:5; Rv_rng = 2.6:1e-4:3.6
x_model = 15000:0.01:17000

function regen_check(adjfib; nchk=8, seed=203)
    fiberindx = adjFiberIndx2FiberIndx(adjfib)
    logmsg("(b) adjfib=$adjfib (LCO fiberindx $fiberindx) consumes:")
    for p in (LSF_LCO, TFUN_LCO, LST_LCO, FRAC_LCO); logmsg("      ", p); end
    Ksp = get_lsf_matrix(adjfib, LSF_LCO)
    Atell = h5open(TFUN_LCO, "r") do f
        permutedims(read(f["design_matrix"]), [2,1])
    end
    lst = deserialize(LST_LCO)
    frac = deserialize(FRAC_LCO)
    rng = MersenneTwister(seed)
    Teff_lst = rand(rng,Teff_rng,nsamp); Av_lst = rand(rng,Av_rng,nsamp); Rv_lst = rand(rng,Rv_rng,nsamp)
    Tfunindx_lst = rand(rng,lst[fiberindx],nsamp); Tfracindx_lst = rand(rng,1:size(frac,2),nsamp)
    X = deserialize(joinpath(NEW, "starCont_"*lpad(adjfib,3,"0")*".jdat"))
    theta = h5open(TFUN_LCO,"r") do f; f["theta"][:,fiberindx,:]; end
    nexact = 0; maxrel = 0.0
    for s in 1:nchk
        bbs = blackbody.(Ref(Teff_lst[s]), x_model*1e-8)
        rvec = redden_mult(x_model,Av_lst[s],Rv_lst[s])
        Tfunsample = exp.(Atell*theta[:,Tfunindx_lst[s]])
        v = frac[:,Tfracindx_lst[s]].*Tfunsample./nanzeromedian(Tfunsample).*(Ksp*(rvec.*bbs))
        rel = maximum(abs.(v .- X[:,s]) ./ max.(abs.(X[:,s]), 1e-30))
        maxrel = max(maxrel, rel)
        nexact += (v == X[:,s])
    end
    logmsg(@sprintf("(b) adjfib=%d: %d/%d columns bit-exact; max rel dev %.3g -> %s",
        adjfib, nexact, nchk, maxrel, (nexact==nchk || maxrel < 1e-12) ? "ROUTING VERIFIED" : "MISMATCH"))
    (nexact==nchk || maxrel < 1e-12)
end
first_unchanged = 300 + minimum(pred.identfib)
verdict_b = all(regen_check.([first_unchanged, 450, 519, 595]))
logmsg("(b) VERDICT: ", verdict_b ? "PASS" : "FAIL")

logmsg(nfail == 0 && verdict_b ? "E4B-SAMPLE-QA-PASS" : "E4B-SAMPLE-QA-FAIL", "  ", now())
close(io)
