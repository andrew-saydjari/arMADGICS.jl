## workup_mpi.jl — W2 MPI tier: rank-strided batch readers, all ranks writing
## the SAME preallocated output files through the parallel-HDF5 MPIO driver.
#
# Shares ALL contract/planning/layout logic with the serial tier:
# RowContract.jl (W1 identity/expected-set/integrity) and WorkupSerial.jl
# (plan_workup, output preallocation, sentinel fills, slab placement). The
# only difference is the I/O topology:
#   * rank 0 plans (identity load + reconcile) and broadcasts the plan;
#   * rank 0 creates every arMADGICS_out_<key>.h5 serially (datasets
#     preallocated, hdr provenance attrs, missing_row_mask + sentinel fill);
#   * ALL ranks then collectively open every output file with
#     Drivers.MPIO(comm) and each rank streams its strided share of batches
#     (present[rank+1 : nranks : end], the gist's assignment) — reading each
#     batch file ONCE (all keys) and writing 37 independent hyperslabs at the
#     contract-derived row ranges. Row placement never depends on rank count.
#   * any batch failing the integrity check aborts the whole job (MPI.Abort).
#
# Environment: MUST run under workup/mpi_env (module-provided MPI + parallel
# HDF5; see setup_mpi_env.sh — cephtweaks + openmpi/5.0.6 + hdf5/mpi-1.14.5,
# HDF5_USE_FILE_LOCKING=FALSE on ceph).
#
#   source mpi_env/setup_mpi_env.sh
#   mpiexec -n 4 julia --project=mpi_env workup_mpi.jl \
#       --rawdir RAW --outdir OUT [--fibers F1:F2] [--allow-missing] \
#       [--batchsize 100] [--progress-every 50]
#
# Raw batch files are NEVER deleted (W5).

using MPI
using HDF5
using Dates
using Printf

include(joinpath(@__DIR__, "RowContract.jl"))
include(joinpath(@__DIR__, "WorkupSerial.jl"))
using .RowContract
using .WorkupSerial

function getopt(flag; default = nothing, required = false)
    i = findfirst(==(flag), ARGS)
    if isnothing(i)
        required && error("$flag is required")
        return default
    end
    return ARGS[i + 1]
end
getflag(flag) = flag in ARGS

function parse_fibers(s::AbstractString)
    m = match(r"^(\d+):(\d+)$", strip(s))
    isnothing(m) && error("--fibers must be F1:F2 (got '$s')")
    return parse(Int, m.captures[1]):parse(Int, m.captures[2])
end

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    rawdir = getopt("--rawdir"; required = true)
    outdir = getopt("--outdir"; required = true)
    fibers = parse_fibers(getopt("--fibers"; default = "1:600"))
    allow_missing = getflag("--allow-missing")
    batchsize = parse(Int, getopt("--batchsize"; default = "100"))
    progress_every = parse(Int, getopt("--progress-every"; default = "50"))

    log0(s) = rank == 0 ? (println("[$(now())] $s"); flush(stdout)) : nothing

    HDF5.has_parallel() ||
        (log0("FATAL: HDF5.jl not built against parallel HDF5"); MPI.Abort(comm, 1))

    # ---- rank 0 plans + creates outputs; everyone else waits -----------------
    t0 = MPI.Wtime()
    plan = nothing
    if rank == 0
        plan = plan_workup(rawdir, outdir; fibers = fibers,
            allow_missing = allow_missing, batchsize = batchsize, log = log0)
        files = WorkupSerial.open_outputs(plan)
        try
            WorkupSerial.write_missing_fills!(files, plan)
        finally
            foreach(close, values(files))
        end
        log0("Outputs preallocated ($(length(plan.keyinfo)) files); broadcasting plan ...")
    end
    plan = MPI.bcast(plan, comm; root = 0)

    # ---- collective MPIO open of every output file ---------------------------
    okeys = sort(collect(keys(plan.keyinfo)))
    files = Dict{String, HDF5.File}()
    for k in okeys   # identical order on every rank: opens are collective
        files[k] = h5open(out_file_path(plan.outdir, k), "r+";
            driver = HDF5.Drivers.MPIO(comm))
    end

    # ---- rank-strided streaming ---------------------------------------------
    myids = plan.present[(rank + 1):nranks:end]
    log0("Streaming $(length(plan.present)) batches over $nranks ranks (~$(length(myids))/rank) ...")
    ndone = 0
    for id in myids
        nrow = length(RowContract.batch_within_range(plan.fidx, id))
        data = try
            WorkupSerial.read_batch(plan.paths[id], id, nrow, plan.keyinfo)
        catch e
            println("rank $rank: ABORT: $(sprint(showerror, e))")
            flush(stdout)
            MPI.Abort(comm, 3)
            error("unreachable")
        end
        rng = out_row_range(plan, id)
        for k in okeys
            WorkupSerial._write_slab!(files[k][k], rng, data[k])
        end
        ndone += 1
        if rank == 0 && (ndone % progress_every == 0 || ndone == length(myids))
            el = MPI.Wtime() - t0
            log0(@sprintf("  rank0 wrote %d / %d of its batches (%.2f batches/s/rank, elapsed %.0fs)",
                ndone, length(myids), ndone / el, el))
        end
    end
    foreach(close, values(files))   # collective closes
    MPI.Barrier(comm)

    if rank == 0
        el = MPI.Wtime() - t0
        rowbytes = sum(prod(v.shape[1:(end - 1)]) * sizeof(v.eltype)
                       for v in values(plan.keyinfo))
        covered = sum(length(out_row_range(plan, id)) for id in plan.present)
        log0(@sprintf("DONE: %d keys × %d rows in %.1fs (%.2f GB written, %.2f GB/s aggregate)",
            length(plan.keyinfo), plan.nsamp_out, el,
            covered * rowbytes / 1e9, covered * rowbytes / 1e9 / el))
    end
    MPI.Finalize()
end

main()
