# E6 step 0a: inspect the production rough prior, chip-gap mask, tellurics lco metadata
using HDF5, Statistics, LinearAlgebra

rough = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/APOGEE_starcont_svd_60_rough.h5"
h5open(rough, "r") do f
    println("rough keys: ", keys(f))
    V = read(f["Vmat"]); lam = read(f["λv"]); msk = read(f["cont_msk"])
    println("Vmat ", size(V), " ", eltype(V), " finite=", all(isfinite, V),
        " sqrt(lam1)=", sqrt(lam[1]), " normcol1=", norm(V[:, 1]))
    println("lam[1:5] = ", lam[1:5])
    println("lam[56:60] = ", lam[56:60])
    println("nmask = ", count(Bool.(msk)), " / ", length(msk))
    println("rows outside mask all zero: ", all(V[.!Bool.(msk), :] .== 0))
    G = transpose(V) * V
    od = maximum(abs(G[i, j]) / sqrt(lam[i] * lam[j]) for i in 1:60, j in 1:60 if i != j)
    dg = maximum(abs.(diag(G) .- lam) ./ lam)
    println("max offdiag(VtV)/sqrt(li lj) = ", od, "   max |diag-lam|/lam = ", dg)
end

cg = "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5"
h5open(cg, "r") do f
    for k in keys(f)
        println("chipgap key ", k, "  ", size(f[k]), " ntrue=", count(Bool.(read(f[k]))))
    end
end

tf = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5"
h5open(tf, "r") do f
    for k in keys(f)
        println("tell lco key: ", k, "  ", size(f[k]), "  ", eltype(f[k]))
    end
end
