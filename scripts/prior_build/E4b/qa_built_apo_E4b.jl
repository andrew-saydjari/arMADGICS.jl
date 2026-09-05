# E4b scope extension QA: the 300 per-fiber APO starCont SVD priors built from the
# canonical 2026_09_03 E4 samples. E6 structural invariants on all 300, plus:
#   - EXACT-REPRODUCTION control on f101/f245/f295/f335: E6 built_new used the SAME
#     samples + builder + mask + env (arM-E6, julia 1.11.0), so these four must
#     reproduce E6's built_new files to numerical precision (bit-comparable).
#   - f101: λ2 vs E6 built_old (2026_04_26 samples) should show the ~0.23x cut
#     effect E6 measured (new APO bright/chi2 tfunlist cuts removed variance).
# Run: cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run && nice -n 10 \
#   julia +1.11.0 --project=/mnt/home/asaydjari/gitcode/worktrees/arM-E6 \
#   /mnt/home/asaydjari/gitcode/worktrees/arM-E4b/scripts/prior_build/E4b/qa_built_apo_E4b.jl
using HDF5, LinearAlgebra, Statistics, Printf, Serialization, Dates

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run"
const BUILT = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_pass1/built_apo"
const E6D = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E6_run"
const NSUB = 60

io = open(joinpath(RUND, "QA_BUILT_APO_REPORT.txt"), "w")
logmsg(args...) = (println(args...); println(io, args...); flush(stdout); flush(io))
logmsg("E4b APO built-prior QA — ", now())

msk_apo = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5", "r") do f
    Bool.(read(f["apo"]))
end
logmsg("apo chipgapmsk count = ", count(msk_apo), " (expect 7742)")

loadprior(path) = h5open(path, "r") do f
    read(f["Vmat"]), read(f["λv"]), Bool.(read(f["chipgapmsk"]))
end
signaligned_cor(a, b) = abs(dot(a, b) / (norm(a) * norm(b)))

fails = String[]
chk(cond, msg) = (cond || push!(fails, msg); cond)
lam1 = zeros(300); lam2 = zeros(300)
for adjfib in 1:300
    fn = joinpath(BUILT, "APOGEE_starcont_svd_60_f" * lpad(adjfib, 3, "0") * ".h5")
    if !isfile(fn)
        push!(fails, "f$adjfib: MISSING"); continue
    end
    V, lam, m = loadprior(fn)
    chk(size(V) == (8700, NSUB) && eltype(V) == Float64, "f$adjfib: Vmat dims/dtype")
    chk(length(lam) == NSUB && all(lam .> 0) && issorted(lam, rev=true), "f$adjfib: λv")
    chk(m == msk_apo, "f$adjfib: chipgapmsk != 2026_04_25 apo mask")
    chk(all(V[.!msk_apo, :] .== 0), "f$adjfib: nonzero outside mask")
    chk(all(isfinite, V) && all(isfinite, lam), "f$adjfib: NaN/Inf")
    G = transpose(V) * V
    od = maximum(abs(G[i, j]) / sqrt(lam[i] * lam[j]) for i in 1:NSUB, j in 1:NSUB if i != j)
    dg = maximum(abs.(diag(G) .- lam) ./ lam)
    chk(od < 1e-8 && dg < 1e-8, "f$adjfib: orthogonality od=$od dg=$dg")
    lam1[adjfib] = lam[1]; lam2[adjfib] = lam[2]
end
logmsg(isempty(fails) ? "STRUCTURAL: PASS on all 300 (S1,S3,S4,S5 + mask)" :
    "STRUCTURAL: FAIL — $(length(fails)) violations")
foreach(m -> logmsg("  ", m), fails[1:min(end, 20)])
logmsg(@sprintf("λ1 across 300 apo fibers: min %.4g med %.4g max %.4g", minimum(lam1), median(lam1), maximum(lam1)))
logmsg(@sprintf("λ2 across 300 apo fibers: min %.4g med %.4g max %.4g", minimum(lam2), median(lam2), maximum(lam2)))
for (nm, v) in (("λ1", lam1), ("λ2", lam2))
    top = sortperm(v, rev=true)[1:5]
    logmsg("  top-5 $nm fibers: ", join(["f$(lpad(i,3,"0"))=$(round(v[i],sigdigits=4))" for i in top], "  "))
end

# exact-reproduction control vs E6 built_new (same samples/builder/mask/env)
repro_ok = true
for fib in (101, 245, 295)   # NOTE: E6's fourth fiber 335 is LCO (adjfib>300), not APO
    global repro_ok
    Vn, ln_, _ = loadprior(joinpath(BUILT, "APOGEE_starcont_svd_60_f" * lpad(fib, 3, "0") * ".h5"))
    Vr, lr, _ = loadprior(joinpath(E6D, "built_new", "APOGEE_starcont_svd_60_f" * lpad(fib, 3, "0") * ".h5"))
    bit = (Vn == Vr) && (ln_ == lr)
    dlam = maximum(abs.(ln_ .- lr) ./ lr)
    dv = maximum([1 - signaligned_cor(Vn[:, k], Vr[:, k]) for k in 1:10])
    ok = bit || (dlam < 1e-10 && dv < 1e-10)
    repro_ok &= ok
    logmsg(@sprintf("CONTROL f%03d vs E6 built_new: bit-identical=%s  max|Δλ|/λ=%.3g  max(1-|cor|)=%.3g -> %s",
        fib, bit, dlam, dv, ok ? "REPRODUCED" : "DEVIATES"))
end

# f101 cut-effect check vs E6 built_old (2026_04_26 samples): λ2 ratio ~0.23
V101, l101, _ = loadprior(joinpath(BUILT, "APOGEE_starcont_svd_60_f101.h5"))
_, l101old, _ = loadprior(joinpath(E6D, "built_old", "APOGEE_starcont_svd_60_f101.h5"))
r2 = l101[2] / l101old[2]
hl101 = 0.15 < r2 < 0.35
logmsg(@sprintf("HEADLINE f101: λ2 new/old(2026_04_26) = %.4f (E6 measured ~0.23) -> %s", r2, hl101 ? "PASS" : "FAIL"))

allpass = isempty(fails) && repro_ok && hl101
logmsg(allpass ? "E4B-BUILT-APO-QA-PASS" : "E4B-BUILT-APO-QA-FAIL", "  ", now())
serialize(joinpath(RUND, "qa_built_apo_results.jdat"), (lam1=lam1, lam2=lam2, fails=fails))
close(io)
