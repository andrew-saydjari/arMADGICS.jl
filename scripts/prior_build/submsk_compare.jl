using HDF5, Printf, Statistics
NEW="/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1/built"
OLD="/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/sky_priors"
CG="/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5"
cg = h5open(CG,"r") do f; (apo=count(read(f["apo"])), lco=count(read(f["lco"]))); end
@printf("chipgap mask pixels: apo=%d lco=%d\n", cg.apo, cg.lco)
newn=Int[]; oldn=Int[]; tel=String[]
for adjfib in vcat(1:20:300, 301:20:600)
    n=lpad(adjfib,3,"0")
    pn=joinpath(NEW,"APOGEE_skyline_faint_svd_120_f$n.h5")
    po=joinpath(OLD,"APOGEE_skyline_faint_svd_120_f$n.h5")
    (isfile(pn) && isfile(po)) || continue
    push!(newn, count(Bool.(h5read(pn,"submsk"))))
    push!(oldn, count(Bool.(h5read(po,"submsk"))))
    push!(tel, adjfib<=300 ? "apo" : "lco")
end
for t in ("apo","lco")
    m = tel .== t
    any(m) || continue
    @printf("%s (n=%d): E5 faint submsk med=%.0f (range %d-%d) | DR17-era dump med=%.0f (range %d-%d)\n",
        t, count(m), median(newn[m]), minimum(newn[m]), maximum(newn[m]),
        median(oldn[m]), minimum(oldn[m]), maximum(oldn[m]))
end
