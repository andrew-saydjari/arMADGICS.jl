# 01_seqmap.jl -- map every telescope-night to its almanac OBJECT SEQUENCES.
#
# CHAIN POSITION: allVisitRV step 1 of 4 (also read by fiber_offsets/01_build_working_tables.jl).
# INPUTS:  ARGS[1]  the almanac file the workup used (raw/<obs>/<mjd>/sequences/objects)
# OUTPUTS: ARGS[2]  seqmap.jls -- Serialization of Dict{(tele::Int8, mjd::Int32) =>
#                   Vector{(first_expnum, last_expnum)}}
# SOURCE:  work/2026_09_16/allvisitrv/s0_seqmap.jl (0.4.0). Logic unchanged; the almanac path
#          and output path became ARGS[1], ARGS[2].
# Usage:   julia --project=scripts/release_tables 01_seqmap.jl ALMANAC.h5 seqmap.jls
using HDF5, Printf
A=ARGS[1]
t0=time()
f=h5open(A,"r"); r=f["raw"]
const RES = Dict("nseq"=>0, "nnight"=>0)
seqmap = Dict{Tuple{Int8,Int32},Vector{Tuple{Int32,Int32}}}()
for (ti,obs) in ((Int8(1),"apo"),(Int8(2),"lco"))
    g=r[obs]
    for m in keys(g)
        nd=g[m]
        haskey(nd,"sequences") || continue
        s=nd["sequences"]
        haskey(s,"objects") || continue
        arr=read(s["objects"])
        size(arr,2)==0 && continue
        v=[(Int32(arr[1,j]), Int32(arr[2,j])) for j in 1:size(arr,2)]
        seqmap[(ti, parse(Int32,m))] = v
        RES["nseq"] += length(v); RES["nnight"] += 1
    end
end
close(f)
println("nights with object sequences = $(RES["nnight"]), total object sequences = $(RES["nseq"])  ($(round(time()-t0,digits=1)) s)")
using Serialization
serialize(ARGS[2], seqmap)
