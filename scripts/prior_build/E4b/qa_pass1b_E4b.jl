# pass-1b QA (intelligent policy delta): the 4 resampled fibers (adjfib
# 388/448/459/519) and their rebuilt priors.
#  (s1) samples: shape/finiteness; per-column comparison vs the ORIGINAL
#       2026_09_03 files must equal the RNG-stream prediction exactly
#       (leak-gone certification: prediction asserts rows 2578/3830 undrawable).
#  (s2) routing: bit-exact regeneration of first 8 columns under BLISBLAS.
#  (b1) built priors: structural invariants; λ comparison vs the superseded
#       C2=3000 builds (variance restored on the FPI fibers is EXPECTED);
#       FPI kept-entry counts reported prominently.
# Run: cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/E4c_run && nice -n 10 \
#   julia +1.11.6 --project=/mnt/home/asaydjari/gitcode/worktrees/arM-E4b qa (this file)
using BLISBLAS, LinearAlgebra
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

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/E4c_run"
const NEW = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_20260905_intelligent/tell_prior_disk"
const ORIG = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/prior_outputs/starCont_20260903/tell_prior_disk"
const NEWB = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_20260905_intelligent/built"
const OLDB = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_pass1/built_lco"
const CHANGED = [388, 448, 459, 519]
const LST = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_intelligent/20260905_lco_tfunlist.jdat"

io = open(joinpath(RUND, "QA_PASS1B_REPORT.txt"), "w")
logmsg(args...) = (println(args...); println(io, args...); flush(stdout); flush(io))
logmsg("pass-1b delta QA — ", now())

pred = deserialize("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/c2_analysis/rng_predict_intelligent.jdat")
lst = deserialize(LST)
logmsg("FPI/changed-fiber kept tfun entries (intelligent policy): ",
    join(["adjfib $(a): $(length(lst[a-300]))" for a in CHANGED], "  "))

# (s1)
allok = true
for a in CHANGED
    Xn = deserialize(joinpath(NEW, "starCont_" * lpad(a, 3, "0") * ".jdat"))
    ok = size(Xn) == (8700, 10000) && all(isfinite, Xn)
    Xo = deserialize(joinpath(ORIG, "starCont_" * lpad(a, 3, "0") * ".jdat"))
    changed = findall(j -> (@view Xn[:, j]) != (@view Xo[:, j]), 1:10000)
    match = changed == pred[a-300]
    global allok &= ok && match
    logmsg(@sprintf("(s1) adjfib %d: shape/finite %s; changed cols vs ORIGINAL %d, == RNG prediction: %s",
        a, ok, length(changed), match))
end
logmsg("(s1) leak-gone: prediction built from the leak-free intelligent list with an assertion rows {2578,3830} are never drawn; matches above certify the stored samples. VERDICT: ", allok ? "PASS" : "FAIL")

# (s2) routing bit-exact
prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
LSF = prior_dir*"2026_08_31/prior_inputs/lsf_20260427/fpiLSFparams_REGULARIZED_lco_60861.h5"
TF = prior_dir*"2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5"
FRAC = prior_dir*"2026_04_26/outsamptell_lco.jdat"
x_model = 15000:0.01:17000
frac = deserialize(FRAC)
Atell = h5open(TF, "r") do f
    permutedims(read(f["design_matrix"]), [2,1])
end
sok = true
for a in CHANGED
    fb = adjFiberIndx2FiberIndx(a)
    Ksp = get_lsf_matrix(a, LSF)
    rng = MersenneTwister(203)
    Te = rand(rng, 4_000:1:10_000, 10_000); Av = rand(rng, 0:1e-4:5, 10_000); Rv = rand(rng, 2.6:1e-4:3.6, 10_000)
    Ti = rand(rng, lst[fb], 10_000); Fi = rand(rng, 1:size(frac, 2), 10_000)
    X = deserialize(joinpath(NEW, "starCont_" * lpad(a, 3, "0") * ".jdat"))
    theta = h5open(TF, "r") do f; f["theta"][:, fb, :]; end
    nex = 0
    for s in 1:8
        bbs = blackbody.(Ref(Te[s]), x_model * 1e-8)
        rvec = redden_mult(x_model, Av[s], Rv[s])
        Ts = exp.(Atell * theta[:, Ti[s]])
        v = frac[:, Fi[s]] .* Ts ./ nanzeromedian(Ts) .* (Ksp * (rvec .* bbs))
        nex += (v == X[:, s])
    end
    global sok &= (nex == 8)
    logmsg("(s2) adjfib $a: $nex/8 bit-exact")
end
logmsg("(s2) VERDICT: ", sok ? "PASS" : "FAIL")

# (b1) built priors
msk_lco = h5open(prior_dir*"2026_04_25/StarContChipGapMsk.h5", "r") do f
    Bool.(read(f["lco"]))
end
lp(p) = h5open(p, "r") do f; (read(f["Vmat"]), read(f["λv"]), Bool.(read(f["chipgapmsk"]))); end
bok = true
for a in CHANGED
    Vn, lam, m = lp(joinpath(NEWB, "APOGEE_starcont_svd_60_f" * lpad(a, 3, "0") * ".h5"))
    G = transpose(Vn) * Vn
    od = maximum(abs(G[i, j]) / sqrt(lam[i] * lam[j]) for i in 1:60, j in 1:60 if i != j)
    ok = size(Vn) == (8700, 60) && all(isfinite, Vn) && m == msk_lco &&
         all(Vn[.!msk_lco, :] .== 0) && issorted(lam, rev=true) && all(lam .> 0) && od < 1e-8
    global bok &= ok
    Vo, lo, _ = lp(joinpath(OLDB, "APOGEE_starcont_svd_60_f" * lpad(a, 3, "0") * ".h5"))
    logmsg(@sprintf("(b1) adjfib %d: structural %s; λ1 new/C2build = %.4f  λ2 = %.3f  λ3 = %.3f (variance restored on FPI is expected)",
        a, ok ? "PASS" : "FAIL", lam[1] / lo[1], lam[2] / lo[2], lam[3] / lo[3]))
end
logmsg("(b1) VERDICT: ", bok ? "PASS" : "FAIL")
logmsg(allok && sok && bok ? "PASS1B-QA-PASS" : "PASS1B-QA-FAIL", "  ", now())
close(io)
