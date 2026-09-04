# E6 step 3: structural regression of NEW vs OLD built starCont priors (+ rough control)
# Usage: julia --project=<arM-E6> compare_structural.jl
using HDF5, LinearAlgebra, Statistics, Printf, Serialization

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E6_run"
const FIBERS = [101, 245, 295, 335, 350, 450, 595]
const NSUB = 60
const KLEAD = 10

cgf = "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5"
msk_apo, msk_lco = h5open(cgf, "r") do f
    Bool.(read(f["apo"])), Bool.(read(f["lco"]))
end

roughf = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/APOGEE_starcont_svd_60_rough.h5"
V_rough, lam_rough = h5open(roughf, "r") do f
    read(f["Vmat"]), read(f["λv"])
end

loadprior(path) = h5open(path, "r") do f
    read(f["Vmat"]), read(f["λv"]), Bool.(read(f["chipgapmsk"]))
end

# principal angles' cosines between the K leading left-singular subspaces
function subspace_cosines(Va, Vb, K)
    Ua = Matrix(qr(Va[:, 1:K]).Q)
    Ub = Matrix(qr(Vb[:, 1:K]).Q)
    svdvals(transpose(Ua) * Ub)
end

signaligned_cor(a, b) = abs(dot(a, b) / (norm(a) * norm(b)))

results = Dict{Int,Dict{String,Any}}()
fails = String[]
chk(cond, msg) = (cond || push!(fails, msg); cond)

println("=== E6 structural regression: NEW (2026_09_03 samples) vs OLD (2026_04_26 samples) builds ===\n")
for fib in FIBERS
    tele = fib > 300 ? "lco" : "apo"
    mskref = fib > 300 ? msk_lco : msk_apo
    fo = joinpath(RUND, "built_old", "APOGEE_starcont_svd_60_f" * lpad(fib, 3, "0") * ".h5")
    fn = joinpath(RUND, "built_new", "APOGEE_starcont_svd_60_f" * lpad(fib, 3, "0") * ".h5")
    Vo, lo, mo = loadprior(fo)
    Vn, ln_, mn = loadprior(fn)

    # S1 schema
    chk(size(Vo) == (8700, NSUB) && size(Vn) == (8700, NSUB), "f$fib: Vmat dims")
    chk(eltype(Vn) == Float64, "f$fib: dtype")
    chk(length(ln_) == NSUB && all(ln_ .> 0) && issorted(ln_, rev=true), "f$fib: λv")
    # S2 mask
    chk(mo == mskref && mn == mskref, "f$fib: chipgapmsk != 2026_04_25 $tele mask")
    # S3 zero-padding
    chk(all(Vo[.!mskref, :] .== 0) && all(Vn[.!mskref, :] .== 0), "f$fib: nonzero outside mask")
    # S4 finite
    chk(all(isfinite, Vo) && all(isfinite, Vn) && all(isfinite, lo) && all(isfinite, ln_), "f$fib: NaN/Inf")
    # S5 orthogonality-with-scale
    for (V, lam, tag) in ((Vo, lo, "old"), (Vn, ln_, "new"))
        G = transpose(V) * V
        od = maximum(abs(G[i, j]) / sqrt(lam[i] * lam[j]) for i in 1:NSUB, j in 1:NSUB if i != j)
        dg = maximum(abs.(diag(G) .- lam) ./ lam)
        chk(od < 1e-8 && dg < 1e-8, "f$fib $tag: orthogonality od=$od dg=$dg")
    end
    # S6 scale
    lamratio = ln_ ./ lo
    chk(0.5 < lamratio[1] < 2.0, "f$fib: λ1 ratio $(lamratio[1]) outside [0.5,2]")

    cors = [signaligned_cor(Vo[:, k], Vn[:, k]) for k in 1:KLEAD]
    cosang = subspace_cosines(Vo, Vn, KLEAD)
    cos_rough_o = subspace_cosines(V_rough, Vo, KLEAD)
    cos_rough_n = subspace_cosines(V_rough, Vn, KLEAD)

    results[fib] = Dict("lam_old" => lo, "lam_new" => ln_, "cors" => cors,
        "cos_oldnew" => cosang, "cos_rough_old" => cos_rough_o, "cos_rough_new" => cos_rough_n,
        "V_old_lead" => Vo[:, 1:KLEAD], "V_new_lead" => Vn[:, 1:KLEAD])

    @printf("fiber %3d (%s): λ1 new/old = %.4f   λ2 = %.4f   λ10 = %.4f\n",
        fib, tele, lamratio[1], lamratio[2], lamratio[10])
    @printf("   |cor(k)| k=1..5: %s\n", join([@sprintf("%.5f", c) for c in cors[1:5]], " "))
    @printf("   subspace cos(theta) old-new (min of 10): %.5f ; full: %s\n",
        minimum(cosang), join([@sprintf("%.4f", c) for c in cosang], " "))
    @printf("   control cos(theta) rough-old min: %.4f | rough-new min: %.4f\n",
        minimum(cos_rough_o), minimum(cos_rough_n))
end

serialize(joinpath(RUND, "structural_results.jdat"), results)

println()
if isempty(fails)
    println("STRUCTURAL-REGRESSION-PASS: all invariants S1-S6 hold for $(length(FIBERS)) fibers x {old,new}")
else
    println("STRUCTURAL-REGRESSION-FAIL:")
    foreach(println, fails)
end
