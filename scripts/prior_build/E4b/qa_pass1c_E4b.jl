# pass-1c QA (FINAL AKS-approved policy, 2026-09-05): all 300 resampled APO
# fibers and their rebuilt priors.
#  (s1) samples: shape/finiteness; per-column comparison vs the ORIGINAL
#       2026_09_03 APO files must equal the RNG-stream prediction exactly for
#       EVERY fiber (removed-entry-undrawable certification lives in the
#       prediction script's assertion).
#  (s2) routing: bit-exact regeneration of the first 8 columns under BLISBLAS
#       for spot fibers (28, 76, 150, 171, 226, 300).
#  (s3) LCO side (pass-1b reused wholesale): bit-exact 2-column regeneration
#       for an unchanged fiber (adjfib 350, lcoC2 corpus) and a pass-1b
#       resampled fiber (adjfib 388) against the FINAL lco list (byte-identical
#       to the intelligent list, so streams must reproduce the stored files).
#  (b1) built priors: structural invariants on all 300; λ comparison vs the
#       superseded pass-1 ceiling builds — variance restored on FPI f76/f226
#       is EXPECTED and is the headline.
# Run: cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/pass1c_run && nice -n 10 \
#   JULIA_PROJECT=/mnt/home/asaydjari/gitcode/worktrees/arM-E4b julia +1.11.6 -p 3 qa (this file)
# ARM_QA_PHASE=samples runs (s1)(s2)(s3) only; =built runs (b1) only; default
# all. Each phase appends to QA_PASS1C_REPORT.txt (samples phase overlaps the
# prior rebuild; (b1) runs once all 300 built h5s exist).
using Distributed
const QA_PHASE = get(ENV, "ARM_QA_PHASE", "all")
const DO_S = QA_PHASE in ("all", "samples")
const DO_B = QA_PHASE in ("all", "built")
using BLISBLAS, LinearAlgebra
BLAS.set_num_threads(1)
using HDF5, Serialization, Statistics, StatsBase, Printf, Dates, Random, SparseArrays
using FITSIO, LowRankOps, EllipsisNotation, ShiftedArrays, JLD2, FileIO
using Interpolations, AstroTime, Suppressor, DataFrames, ProgressMeter
using SortFilters, BasisFunctions, DustExtinction, DelimitedFiles
using ApogeeReduction: adjFiberIndx2FiberIndx, get_lsf_matrix

proj_path = "/mnt/home/asaydjari/gitcode/worktrees/arM-E4b/"
for f in ["src/utils.jl", "src/gridSearch.jl", "src/componentAndPosteriors.jl", "src/fileNameHandling.jl",
    "src/ingest.jl", "src/lowRankPrescription.jl", "src/marginalizeEW.jl", "src/spectraInterpolation.jl",
    "src/chi2Wrappers.jl", "scripts/prior_build/prior_utils.jl"]
    include(proj_path * f)
end

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/pass1c_run"
const NEW = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_20260905_final/tell_prior_disk"
const ORIG = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/prior_outputs/starCont_20260903/tell_prior_disk"
const NEWB = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starCont_20260905_final/built"
const OLDB = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_pass1/built_apo"
const LSTA = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_final/20260905_apo_tfunlist.jdat"
const LSTL = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_final/20260905_lco_tfunlist.jdat"
const PRED = joinpath(RUND, "rng_predict_final_apo.jdat")
const prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"

io = open(joinpath(RUND, "QA_PASS1C_REPORT.txt"), QA_PHASE == "built" ? "a" : "w")
logmsg(args...) = (println(args...); println(io, args...); flush(stdout); flush(io))
logmsg("pass-1c QA (final policy, phase=$(QA_PHASE)) — ", now())

lst = deserialize(LSTA)
logmsg("FPI kept tfun entries (final policy): f76: $(length(lst[76]))  f226: $(length(lst[226]))")

s1ok = sok = s3ok = bok = true

if DO_S
# ── (s1) all 300 fibers: per-column delta == RNG prediction ─────────────────
@everywhere begin
    using Serialization
    const NEW_W = $NEW
    const ORIG_W = $ORIG
    function s1_check(fb, predcols)
        f3 = lpad(fb, 3, "0")
        Xn = deserialize(joinpath(NEW_W, "starCont_$(f3).jdat"))
        ok = size(Xn) == (8700, 10000) && all(isfinite, Xn)
        Xo = deserialize(joinpath(ORIG_W, "starCont_$(f3).jdat"))
        changed = findall(j -> (@view Xn[:, j]) != (@view Xo[:, j]), 1:10000)
        (ok, changed == predcols, length(changed))
    end
end
pred = deserialize(PRED)
predv = [pred[fb] for fb in 1:300]   # positional args, so pmap ships only the per-task column list
res = @showprogress pmap(s1_check, 1:300, predv)
s1ok = all(r -> r[1] && r[2], res)
nch = [r[3] for r in res]
logmsg(@sprintf("(s1) 300/300 fibers: shape/finite %d ok; prediction-match %d ok; changed cols min %d med %.0f max %d",
    count(r -> r[1], res), count(r -> r[2], res), minimum(nch), median(nch), maximum(nch)))
for fb in findall(r -> !(r[1] && r[2]), res)
    logmsg("(s1) FAIL fiber $fb: ok=$(res[fb][1]) match=$(res[fb][2])")
end
logmsg("(s1) VERDICT: ", s1ok ? "PASS" : "FAIL")

# ── (s2) bit-exact regeneration, first 8 columns, APO spot fibers ───────────
x_model = 15000:0.01:17000
function regen_check(adjfib, tele, lstT, ncols, srcdir)
    fb = adjFiberIndx2FiberIndx(adjfib)
    LSF = prior_dir * "2026_08_31/prior_inputs/lsf_20260427/fpiLSFparams_REGULARIZED_$(tele)_60861.h5"
    TF = prior_dir * "2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_$(tele).h5"
    FRAC = tele == "apo" ? prior_dir * "2026_04_25/outsamptell_apo.jdat" :
           prior_dir * "2026_04_26/outsamptell_lco.jdat"
    frac = deserialize(FRAC)
    Atell = h5open(TF, "r") do f
        permutedims(read(f["design_matrix"]), [2, 1])
    end
    Ksp = get_lsf_matrix(adjfib, LSF)
    rng = MersenneTwister(203)
    Te = rand(rng, 4_000:1:10_000, 10_000); Av = rand(rng, 0:1e-4:5, 10_000); Rv = rand(rng, 2.6:1e-4:3.6, 10_000)
    Ti = rand(rng, lstT[fb], 10_000); Fi = rand(rng, 1:size(frac, 2), 10_000)
    X = deserialize(joinpath(srcdir, "starCont_" * lpad(adjfib, 3, "0") * ".jdat"))
    theta = h5open(TF, "r") do f
        f["theta"][:, fb, :]
    end
    nex = 0
    for s in 1:ncols
        bbs = blackbody.(Ref(Te[s]), x_model * 1e-8)
        rvec = redden_mult(x_model, Av[s], Rv[s])
        Ts = exp.(Atell * theta[:, Ti[s]])
        v = frac[:, Fi[s]] .* Ts ./ nanzeromedian(Ts) .* (Ksp * (rvec .* bbs))
        nex += (v == X[:, s])
    end
    nex
end
sok = true
for a in (28, 76, 150, 171, 226, 300)
    nex = regen_check(a, "apo", lst, 8, NEW)
    global sok &= (nex == 8)
    logmsg("(s2) adjfib $a: $nex/8 bit-exact")
end
logmsg("(s2) VERDICT: ", sok ? "PASS" : "FAIL")

# ── (s3) LCO reuse certification (pass-1b corpus under the final list) ──────
lstl = deserialize(LSTL)
s3ok = true
for (a, src) in ((350, prior_dir * "2026_09_04/prior_outputs/starCont_20260904_lcoC2/tell_prior_disk"),
    (388, prior_dir * "2026_09_05/prior_outputs/starCont_20260905_intelligent/tell_prior_disk"))
    nex = regen_check(a, "lco", lstl, 2, src)
    global s3ok &= (nex == 2)
    logmsg("(s3) adjfib $a (reused corpus): $nex/2 bit-exact under FINAL lco list")
end
logmsg("(s3) VERDICT: ", s3ok ? "PASS" : "FAIL")
end  # DO_S

if DO_B
# ── (b1) built priors: structural invariants on 300; λ restoration ──────────
msk_apo = h5open(prior_dir * "2026_04_25/StarContChipGapMsk.h5", "r") do f
    Bool.(read(f["apo"]))
end
lp(p) = h5open(p, "r") do f
    (read(f["Vmat"]), read(f["λv"]), Bool.(read(f["chipgapmsk"])))
end
bok = true
lr2 = zeros(300)
for a in 1:300
    Vn, lam, m = lp(joinpath(NEWB, "APOGEE_starcont_svd_60_f" * lpad(a, 3, "0") * ".h5"))
    G = transpose(Vn) * Vn
    od = maximum(abs(G[i, j]) / sqrt(lam[i] * lam[j]) for i in 1:60, j in 1:60 if i != j)
    ok = size(Vn) == (8700, 60) && all(isfinite, Vn) && m == msk_apo &&
         all(Vn[.!msk_apo, :] .== 0) && issorted(lam, rev=true) && all(lam .> 0) && od < 1e-8
    global bok &= ok
    ok || logmsg("(b1) FAIL structural adjfib $a")
    _, lo, _ = lp(joinpath(OLDB, "APOGEE_starcont_svd_60_f" * lpad(a, 3, "0") * ".h5"))
    lr2[a] = lam[2] / lo[2]
    if a in (28, 76, 171, 226)
        logmsg(@sprintf("(b1) adjfib %d: structural %s; λ1 new/pass1 = %.4f  λ2 = %.3f  λ3 = %.3f%s",
            a, ok ? "PASS" : "FAIL", lam[1] / lo[1], lam[2] / lo[2], lam[3] / lo[3],
            a in (76, 226) ? "  (FPI — variance restoration EXPECTED)" : ""))
    end
end
logmsg(@sprintf("(b1) structural: %s on 300/300; λ2 ratio vs pass-1 fleet: med %.4f  q90 %.3f  f76 %.3f  f226 %.3f",
    bok ? "PASS" : "FAIL", median(lr2), quantile(lr2, 0.9), lr2[76], lr2[226]))
logmsg("(b1) VERDICT: ", bok ? "PASS" : "FAIL")
end  # DO_B

phase_ok = (!DO_S || (s1ok && sok && s3ok)) && (!DO_B || bok)
logmsg(phase_ok ? "PASS1C-QA-$(uppercase(QA_PHASE))-PASS" : "PASS1C-QA-$(uppercase(QA_PHASE))-FAIL", "  ", now())
close(io)
