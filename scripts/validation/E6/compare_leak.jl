# E6 step 5c: quantify the leaked-exposure imprint — full NEW build vs drop-variant
using HDF5, LinearAlgebra, Statistics, Printf

const BND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E6_run/built_new"

msk_lco = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5", "r") do f
    collect(Bool, read(f["lco"]) .!= 0)
end

for fib in (450, 595)
    fm = joinpath(BND, "APOGEE_starcont_svd_60_f$(fib).h5")
    fd = joinpath(BND, "APOGEE_starcont_svd_60_f$(fib)_dropleak.h5")
    # append mask if the failed final write left it off
    h5open(fd, "r+") do f
        haskey(f, "chipgapmsk") || (f["chipgapmsk"] = msk_lco)
    end
    Vm, lm = h5open(fm, "r") do f; read(f["Vmat"]), read(f["λv"]); end
    Vd, ld = h5open(fd, "r") do f; read(f["Vmat"]), read(f["λv"]); end
    println("=== f$fib: NEW full (10000 cols) vs dropleak ===")
    @printf("  λ ratio drop/full - 1: k=1..6: %s\n",
        join([@sprintf("%+.3e", ld[k] / lm[k] - 1) for k in 1:6], " "))
    @printf("  max_k |λ ratio - 1| over 60: %.3e at k=%d\n",
        maximum(abs.(ld ./ lm .- 1)), argmax(abs.(ld ./ lm .- 1)))
    cors = [abs(dot(Vm[:, k], Vd[:, k]) / (norm(Vm[:, k]) * norm(Vd[:, k]))) for k in 1:10]
    @printf("  |cor| matched k=1..10: %s\n", join([@sprintf("%.4f", c) for c in cors], " "))
    # where does the change live? RMS of component delta per k (sign aligned)
    for k in 1:4
        s = sign(dot(Vm[:, k], Vd[:, k]))
        @printf("  k=%d: rms(ΔV)=%.3e  max|ΔV|=%.3e  (col norm %.3e)\n",
            k, sqrt(mean((s .* Vd[:, k] .- Vm[:, k]) .^ 2)),
            maximum(abs.(s .* Vd[:, k] .- Vm[:, k])), norm(Vm[:, k]))
    end
end
