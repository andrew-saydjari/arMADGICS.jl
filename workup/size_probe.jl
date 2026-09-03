## size_probe.jl — data-side inputs for run_workup.sh's WORKUP_RANKS=auto
## sizing: per-batch payload bytes (all keys of a sample batch, the same
## discovery RowContract's startup performs) and the number of batch files
## (work units) in the requested fiber window.
#
#   julia --project=<workup dir> size_probe.jl <rawdir> <f1> <f2>
#
# stdout (machine-parsed by run_workup.sh):
#   BATCH_BYTES=<sum over non-hdr datasets of prod(shape)*sizeof(eltype)>
#   NBATCH=<batch-file count in fiber dirs f1:f2>
#   SAMPLE=<basename of the probed batch>
#
# Deliberately does NOT load full_list_info.h5 (the ~1-minute identity load
# belongs to the writer, not to sizing); directory listings are bounded and
# non-recursive (one readdir per 3-digit fiber dir).

using HDF5

length(ARGS) == 3 || (println("SIZE_PROBE_ERROR=usage"); exit(1))
rawdir = ARGS[1]
f1 = parse(Int, ARGS[2])
f2 = parse(Int, ARGS[3])

const BATCH_RE = r"^arMADGICS_fiber_\d{3}_batch_\d{7}\.h5$"

nbatch = 0
sample = ""
for entry in sort(readdir(rawdir))
    occursin(r"^\d{3}$", entry) || continue
    f = parse(Int, entry)
    (f1 <= f <= f2) || continue
    sub = joinpath(rawdir, entry)
    isdir(sub) || continue
    for fn in readdir(sub)
        occursin(BATCH_RE, fn) || continue
        global nbatch += 1
        isempty(sample) && (global sample = joinpath(sub, fn))
    end
end

if nbatch == 0
    println("SIZE_PROBE_ERROR=no_batches_in_fiber_window")
    exit(1)
end

bytes = 0
h5open(sample, "r") do io
    for k in keys(io)
        k == "hdr" && continue
        obj = io[k]
        obj isa HDF5.Dataset || continue
        global bytes += prod(size(obj)) * sizeof(eltype(obj))
    end
end

println("BATCH_BYTES=$bytes")
println("NBATCH=$nbatch")
println("SAMPLE=$(basename(sample))")
