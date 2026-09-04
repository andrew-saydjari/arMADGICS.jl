# E4b step 6: QA the 300 per-fiber LCO starCont SVD priors built from the
# regenerated (C2_LCO=3000) samples. E6 structural invariants on all 300, plus
# headline leak-deflation checks against the E6 reference builds:
#   - f450: λ2 must deflate ~4.2x vs the leaked 20260903-sample build
#     (E6 built_new/f450) and match the E6 drop-variant (f450_dropleak) closely.
#   - f595: mid-mode (k~6-10) deflation consistent with f595_dropleak.
#   - f350 (fiberindx 50 tfunlist unchanged -> samples bit-identical): build must
#     reproduce E6 built_new/f350 to numerical precision (cross-Manifest control).
# Run: cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run && nice -n 10 \
#   julia +1.11.0 --project=/mnt/home/asaydjari/gitcode/worktrees/arM-E6 \
#   /mnt/home/asaydjari/gitcode/worktrees/arM-E4b/scripts/prior_build/E4b/qa_built_E4b.jl
using HDF5, LinearAlgebra, Statistics, Printf, Serialization, Dates

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run"
const BUILT = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_pass1/built_lco"
const E6D = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E6_run"
const NSUB = 60

io = open(joinpath(RUND, "QA_BUILT_REPORT.txt"), "w")
logmsg(args...) = (println(args...); println(io, args...); flush(stdout); flush(io))
logmsg("E4b built-prior QA — ", now())

msk_lco = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5", "r") do f
    Bool.(read(f["lco"]))
end
logmsg("lco chipgapmsk count = ", count(msk_lco), " (expect 7833)")

loadprior(path) = h5open(path, "r") do f
    read(f["Vmat"]), read(f["λv"]), Bool.(read(f["chipgapmsk"]))
end
signaligned_cor(a, b) = abs(dot(a, b) / (norm(a) * norm(b)))

# ── structural invariants on all 300 ─────────────────────────────────────────────────────
fails = String[]
chk(cond, msg) = (cond || push!(fails, msg); cond)
lam1 = zeros(300); lam2 = zeros(300)
for adjfib in 301:600
    fn = joinpath(BUILT, "APOGEE_starcont_svd_60_f" * lpad(adjfib, 3, "0") * ".h5")
    if !isfile(fn)
        push!(fails, "f$adjfib: MISSING"); continue
    end
    V, lam, m = loadprior(fn)
    chk(size(V) == (8700, NSUB) && eltype(V) == Float64, "f$adjfib: Vmat dims/dtype")
    chk(length(lam) == NSUB && all(lam .> 0) && issorted(lam, rev=true), "f$adjfib: λv")
    chk(m == msk_lco, "f$adjfib: chipgapmsk != 2026_04_25 lco mask")
    chk(all(V[.!msk_lco, :] .== 0), "f$adjfib: nonzero outside mask")
    chk(all(isfinite, V) && all(isfinite, lam), "f$adjfib: NaN/Inf")
    G = transpose(V) * V
    od = maximum(abs(G[i, j]) / sqrt(lam[i] * lam[j]) for i in 1:NSUB, j in 1:NSUB if i != j)
    dg = maximum(abs.(diag(G) .- lam) ./ lam)
    chk(od < 1e-8 && dg < 1e-8, "f$adjfib: orthogonality od=$od dg=$dg")
    lam1[adjfib-300] = lam[1]; lam2[adjfib-300] = lam[2]
end
logmsg(isempty(fails) ? "STRUCTURAL: PASS on all 300 (S1,S3,S4,S5 + mask)" :
    "STRUCTURAL: FAIL — $(length(fails)) violations")
foreach(m -> logmsg("  ", m), fails[1:min(end, 20)])
logmsg(@sprintf("λ1 across 300 lco fibers: min %.4g med %.4g max %.4g", minimum(lam1), median(lam1), maximum(lam1)))
logmsg(@sprintf("λ2 across 300 lco fibers: min %.4g med %.4g max %.4g", minimum(lam2), median(lam2), maximum(lam2)))

# ── headline checks vs E6 references ─────────────────────────────────────────────────────
function compare_pair(tag, fnew, fref; klead=10)
    Vn, ln_, _ = loadprior(fnew)
    Vr, lr, _ = loadprior(fref)
    rat = ln_ ./ lr
    cors = [signaligned_cor(Vn[:, k], Vr[:, k]) for k in 1:klead]
    logmsg(@sprintf("%s: λ ratio k=1..6: %s", tag, join([@sprintf("%.4f", r) for r in rat[1:6]], " ")))
    logmsg(@sprintf("%s: |cor| k=1..%d: %s", tag, klead, join([@sprintf("%.4f", c) for c in cors], " ")))
    rat, cors
end

f450new = joinpath(BUILT, "APOGEE_starcont_svd_60_f450.h5")
rat_leak, _ = compare_pair("f450 REGEN vs LEAKED build (E6 built_new)", f450new,
    joinpath(E6D, "built_new", "APOGEE_starcont_svd_60_f450.h5"))
logmsg(@sprintf("HEADLINE f450: λ2 deflation vs leaked build = %.2fx (expect ~4.2x)", 1 / rat_leak[2]))
rat_drop, cors_drop = compare_pair("f450 REGEN vs E6 drop-variant", f450new,
    joinpath(E6D, "built_new", "APOGEE_starcont_svd_60_f450_dropleak.h5"))
hl450 = 3.0 < 1 / rat_leak[2] < 6.0 && 0.8 < rat_drop[2] < 1.25
logmsg("HEADLINE f450 VERDICT: ", hl450 ? "PASS" : "FAIL")

f595new = joinpath(BUILT, "APOGEE_starcont_svd_60_f595.h5")
rat595_leak, _ = compare_pair("f595 REGEN vs LEAKED build (E6 built_new)", f595new,
    joinpath(E6D, "built_new", "APOGEE_starcont_svd_60_f595.h5"))
rat595_drop, _ = compare_pair("f595 REGEN vs E6 drop-variant", f595new,
    joinpath(E6D, "built_new", "APOGEE_starcont_svd_60_f595_dropleak.h5"))
# drop-variant's own deflation prediction: λ(dropleak)/λ(leaked full)
_, ldrop595, _ = loadprior(joinpath(E6D, "built_new", "APOGEE_starcont_svd_60_f595_dropleak.h5"))
_, lfull595, _ = loadprior(joinpath(E6D, "built_new", "APOGEE_starcont_svd_60_f595.h5"))
logmsg(@sprintf("f595: mid-mode λ6-10 REGEN/LEAKED: %s (drop-variant predicted %s)",
    join([@sprintf("%.3f", r) for r in rat595_leak[6:10]], " "),
    join([@sprintf("%.3f", r) for r in (ldrop595 ./ lfull595)[6:10]], " ")))
hl595 = rat595_leak[1] > 0.99 && rat595_leak[1] < 1.01
logmsg("f595 λ1 stable: ", hl595 ? "PASS" : "FAIL")

# control: f350 (fiberindx 50) — if its tfunlist was unchanged its samples are
# bit-identical, so the build must reproduce E6 built_new/f350 to numerical precision
pred = deserialize(joinpath(RUND, "rng_predict_final.jdat"))
if 50 in pred.identfib
    Vn, ln_, _ = loadprior(joinpath(BUILT, "APOGEE_starcont_svd_60_f350.h5"))
    Vr, lr, _ = loadprior(joinpath(E6D, "built_new", "APOGEE_starcont_svd_60_f350.h5"))
    dlam = maximum(abs.(ln_ .- lr) ./ lr)
    dv = maximum([1 - signaligned_cor(Vn[:, k], Vr[:, k]) for k in 1:10])
    logmsg(@sprintf("CONTROL f350 (unchanged samples): max|Δλ|/λ = %.3g, max(1-|cor|) k=1..10 = %.3g -> %s",
        dlam, dv, (dlam < 1e-6 && dv < 1e-6) ? "REPRODUCED" : "DEVIATES (check numerics)"))
else
    logmsg("CONTROL f350: fiberindx 50 tfunlist CHANGED — control not applicable")
end

allpass = isempty(fails) && hl450 && hl595
logmsg(allpass ? "E4B-BUILT-QA-PASS" : "E4B-BUILT-QA-FAIL", "  ", now())
serialize(joinpath(RUND, "qa_built_results.jdat"),
    (lam1=lam1, lam2=lam2, fails=fails))
close(io)
