# E4b QA (b) FIX: the routing audit inside qa_samples_E4b.jl ran under default
# OpenBLAS and showed the known ~5e-6 Float32-gemm deviations (same signature as
# E4's first-pass audit; see E4 QA_ROUTING_FIX.txt). Rerun (b) under run-identical
# numerics: BLISBLAS + 1 BLAS thread. Appends to QA_SAMPLES_REPORT.txt.
# Run: cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run && nice -n 10 \
#   julia +1.11.6 --project=/mnt/home/asaydjari/gitcode/worktrees/arM-E4b \
#   /mnt/home/asaydjari/gitcode/worktrees/arM-E4b/scripts/prior_build/E4b/qa_routing_fix_E4b.jl
using BLISBLAS
using LinearAlgebra
BLAS.set_num_threads(1)
using HDF5, Serialization, Statistics, StatsBase, Printf, Dates, Random, SparseArrays
using FITSIO, LowRankOps, EllipsisNotation, ShiftedArrays, JLD2, FileIO
using Interpolations, AstroTime, Suppressor, DataFrames, ProgressMeter
using SortFilters, BasisFunctions, DustExtinction, DelimitedFiles
using ApogeeReduction: adjFiberIndx2FiberIndx, get_lsf_matrix

proj_path = "/mnt/home/asaydjari/gitcode/worktrees/arM-E4b/"
for f in ["src/utils.jl","src/gridSearch.jl","src/componentAndPosteriors.jl","src/fileNameHandling.jl",
          "src/ingest.jl","src/lowRankPrescription.jl","src/marginalizeEW.jl","src/spectraInterpolation.jl",
          "src/chi2Wrappers.jl","scripts/prior_build/prior_utils.jl"]
    include(proj_path*f)
end

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run"
const NEW = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_20260904_lcoC2/tell_prior_disk"
prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
LSF_LCO = prior_dir*"2026_08_31/prior_inputs/lsf_20260427/fpiLSFparams_REGULARIZED_lco_60861.h5"
TFUN_LCO = prior_dir*"2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5"
LST_LCO = prior_dir*"2026_09_04/tfunlists_lcoC2/20260904_lco_tfunlist.jdat"
FRAC_LCO = prior_dir*"2026_04_26/outsamptell_lco.jdat"
nsamp = 10_000
Teff_rng = 4_000:1:10_000; Av_rng = 0:1e-4:5; Rv_rng = 2.6:1e-4:3.6
x_model = 15000:0.01:17000

io = open(joinpath(RUND, "QA_SAMPLES_REPORT.txt"), "a")
logmsg(args...) = (println(args...); println(io, args...); flush(stdout); flush(io))
logmsg("\n--- (b) ROUTING AUDIT RERUN under BLISBLAS + 1 BLAS thread — ", now(), " ---")
logmsg("BLAS config: ", BLAS.get_config())

function regen_check(adjfib; nchk=8, seed=203)
    fiberindx = adjFiberIndx2FiberIndx(adjfib)
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
    logmsg(@sprintf("(b) adjfib=%d (lcoC2 list): %d/%d columns bit-exact; max rel dev %.3g -> %s",
        adjfib, nexact, nchk, maxrel, (nexact==nchk || maxrel < 1e-12) ? "ROUTING VERIFIED" : "MISMATCH"))
    (nexact==nchk || maxrel < 1e-12)
end
verdict_b = all(regen_check.([334, 450, 519, 595]))
logmsg("(b) RERUN VERDICT: ", verdict_b ? "PASS — bit-exact under run-identical numerics; earlier MISMATCH lines were the known OpenBLAS-vs-BLIS Float32 gemm rounding (~6e-6), not a routing failure" : "FAIL")
logmsg(verdict_b ? "E4B-SAMPLE-QA-PASS (with (b) rerun)" : "E4B-SAMPLE-QA-FAIL", "  ", now())
close(io)
