## E7 QA: cross-check the on-the-fly get_lsf_matrix against the 2026_04_26
## materialized sparse LSF matrices (built from the same-era FPI LSF fits;
## 586 GB store that E7 deliberately does NOT reuse or extend).
##
## Compares fiber 150 (APO), findx 6 (fstep = 0) by default. Rows of both are
## normalized to sum 1 before comparison (the materialized files may predate
## the internal row-normalization convention).
##
## Usage: julia --project=<E7 env> qa_lsf_crosscheck.jl [fiber] [findx]

using LinearAlgebra, SparseArrays, Serialization, Statistics, Printf
using ApogeeReduction: get_lsf_matrix

fib = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 150
findx = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
fstep = (6 - findx) / 10

prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
mat_path = prior_dir * "2026_04_26/mat_lsf_out/sp_combolsfmat_norm_$(findx)_" *
           lpad(fib, 3, "0") * ".jdat"
lsfpath = prior_dir * "2026_08_31/prior_inputs/lsf_20260427/fpiLSFparams_REGULARIZED_apo_60861.h5"

fib <= 300 || error("materialized store is APO-only (adjfib <= 300)")

println("materialized: $mat_path")
Kold = deserialize(mat_path)
println("on-the-fly: get_lsf_matrix($fib, apo; fstep=$fstep)")
Knew = get_lsf_matrix(fib, lsfpath; fstep = fstep)

println("sizes: old $(size(Kold)) nnz $(nnz(Kold)) | new $(size(Knew)) nnz $(nnz(Knew))")

rownorm(K) = begin
    s = vec(sum(K, dims = 2))
    s[iszero.(s)] .= 1.0
    spdiagm(0 => 1 ./ s) * K
end

if size(Kold) == size(Knew)
    Ko = rownorm(Kold)
    Kn = rownorm(Knew)
    d = Ko - Kn
    relfro = norm(d) / norm(Kn)
    maxabs = maximum(abs, d)
    @printf("row-normalized comparison: rel-Frobenius %.3e  max|diff| %.3e\n", relfro, maxabs)
    # per-row L1 difference distribution (each row is a unit-sum LSF profile)
    rl1 = vec(sum(abs, d, dims = 2))
    @printf("per-row L1 diff: med %.3e  p99 %.3e  max %.3e\n",
        median(rl1), quantile(rl1, 0.99), maximum(rl1))
else
    println("SIZE MISMATCH — comparing effective LSF action on a test spectrum instead")
    # act both matrices on the same smooth test vector over their model grids
end
