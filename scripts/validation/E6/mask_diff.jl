# E6: quantify support difference between the production rough prior mask and the
# current 2026_04_25 apo chip-gap mask (both apply to the apo fixture).
using HDF5, Statistics, LinearAlgebra
r = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/APOGEE_starcont_svd_60_rough.h5", "r") do f
    (V=read(f["Vmat"]), m=read(f["cont_msk"]) .!= 0)
end
apo = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5", "r") do f
    read(f["apo"]) .!= 0
end
println("rough mask true: ", count(r.m), " ; apo2026 mask true: ", count(apo))
println("in rough but not apo2026 (pixels LOSING prior support): ", count(r.m .& .!apo))
println("in apo2026 but not rough (pixels GAINING support): ", count(apo .& .!r.m))
# where are the losing pixels?
loss = findall(r.m .& .!apo)
function run_ranges(idx)
    runs = Tuple{Int,Int}[]
    isempty(idx) && return runs
    s = idx[1]; p = idx[1]
    for i in idx[2:end]
        if i == p + 1; p = i; else; push!(runs, (s, p)); s = i; p = i; end
    end
    push!(runs, (s, p))
    runs
end
println("lost-support pixel runs (index ranges): ", run_ranges(loss))
# how much variance does rough carry on those pixels?
if !isempty(loss)
    v_loss = sum(abs2, r.V[loss, :]); v_tot = sum(abs2, r.V)
    println("fraction of rough prior power on lost pixels: ", v_loss / v_tot)
end
