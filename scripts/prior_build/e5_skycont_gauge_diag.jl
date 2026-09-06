## Diagnose WHY the skyCont rebuild is not bitwise identical: real difference,
## SVD sign/rotation gauge freedom, or BLAS non-determinism?
using HDF5, LinearAlgebra, Printf, Statistics
const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
const OUT = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1"
const CHIPGAP = "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5"
include(ARM * "scripts/prior_build/build_sky_defs.jl")

scratch = mktempdir()
for adjfib in parse.(Int, ARGS)
    n = lpad(adjfib, 3, "0")
    ship = joinpath(OUT, "built", "APOGEE_skycont_svd_30_f$n.h5")
    fresh = build_skyCont(adjfib; sample_dir=joinpath(OUT, "samples"),
        chipgap_msk_path=CHIPGAP, out_dir=scratch, nsub=30)
    V1 = h5read(ship, "Vmat"); V2 = h5read(fresh, "Vmat")
    l1 = h5read(ship, "λv"); l2 = h5read(fresh, "λv")

    @printf("\nf%s  nsub=%d\n", n, size(V1, 2))
    @printf("  eigenvalues: max REL diff = %.3e   (l1[1]=%.6e l2[1]=%.6e)\n",
        maximum(abs.(l1 .- l2) ./ max.(abs.(l1), 1e-300)), l1[1], l2[1])

    # per-column, sign-corrected relative difference (SVD column signs are arbitrary)
    relsign = Float64[]
    for j in 1:size(V1, 2)
        a = view(V1, :, j); b = view(V2, :, j)
        nrm = max(norm(a), 1e-300)
        push!(relsign, min(norm(a .- b), norm(a .+ b)) / nrm)
    end
    @printf("  per-column sign-corrected rel diff: med=%.3e max=%.3e (worst col %d)\n",
        median(relsign), maximum(relsign), argmax(relsign))
    nflip = count(j -> norm(view(V1, :, j) .+ view(V2, :, j)) < norm(view(V1, :, j) .- view(V2, :, j)), 1:size(V1, 2))
    @printf("  columns whose sign flipped between builds: %d / %d\n", nflip, size(V1, 2))

    # subspace agreement (gauge-invariant): principal angles
    for k in (5, 10, 30)
        k > size(V1, 2) && continue
        Q1 = Matrix(qr(V1[:, 1:k]).Q); Q2 = Matrix(qr(V2[:, 1:k]).Q)
        s = svdvals(transpose(Q1) * Q2)
        @printf("  span k=%2d: worst principal angle = %.3e deg\n", k, acosd(clamp(minimum(s), 0, 1)))
    end
end
rm(scratch, recursive=true, force=true)
