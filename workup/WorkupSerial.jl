## WorkupSerial.jl — W2 serial tier: Distributed readers + a single streaming
## writer, built on the W1 row contract (RowContract.jl).
#
# Design:
#   * Startup is EXACTLY RowContract's expected-set construction: identity from
#     full_list_info.h5, ordering contract asserted (verify_fiber_blocks),
#     expected batches derived from identity, discovery reconciled against it.
#     Missing batches are a hard error unless `allow_missing`, which writes a
#     `missing_row_mask` dataset into every output file and fills the missing
#     rows with type-appropriate sentinels (RowContract.sentinel_for).
#   * One output file per key, `<outdir>/arMADGICS_out_<key>.h5`, preallocated
#     to the full identity-derived row count (nsamp of the fiber window) and
#     filled by slab writes at each batch's contract-derived row range. Row
#     placement NEVER depends on discovery order or arrival order.
#   * ALL keys are discovered from a reference batch (no hardcoded key tuple);
#     the `hdr` DATASET attrs (git provenance) are copied onto the `hdr`
#     dataset of every output file, matching the historical workup.jl layout.
#   * Every batch passes an integrity check (openable, expected key set,
#     eltypes, leading axes, last axis == identity-derived batch length) on the
#     reader BEFORE its rows are written; any problem aborts the whole run.
#   * Optional contiguous fiber window (`fibers=f1:f2`) for subset workups
#     (e.g. regression tests): output rows are the identity rows of those
#     fibers, i.e. global rows (offset+1):(offset+nsamp_out); the offset is
#     recorded as a file-level attribute `workup_row_offset`.
#   * Raw batch files are NEVER deleted (W5: deletion is a separate, later,
#     validation-gated step).
#
# Distribution model: worker processes read + integrity-check batches and ship
# them to the master over a RemoteChannel; the master is the single writer.
# With no workers (nworkers()==1) everything runs inline in the master — the
# code path unit tests exercise.

if !isdefined(@__MODULE__, :RowContract)
    include(joinpath(@__DIR__, "RowContract.jl"))
end

module WorkupSerial

using ..RowContract
using Distributed
using Dates
using HDF5
using Printf

export WorkupPlan, plan_workup, run_workup, out_file_path, out_row_range

const HDR_PLACEHOLDER = "This is only a header"   # historical workup.jl layout
const CKPT_NAME = "workup_serial.ckpt"

# ------------------------------------------------------------------- planning

struct WorkupPlan
    rawdir::String
    outdir::String
    fibers::UnitRange{Int}
    batchsize::Int
    allow_missing::Bool
    nsamp_out::Int
    row_offset::Int              # global row offset of the fiber window
    fidx::RowContract.FiberIndex
    present::Vector{BatchId}
    missing::Vector{BatchId}
    paths::Dict{BatchId, String}
    keyinfo::Dict{String, <:NamedTuple}
    hdr::Dict{String, Any}
    refbatch::BatchId
end

"Row range of batch `id` in the OUTPUT datasets (identity-derived, window-shifted)."
out_row_range(plan::WorkupPlan, id::BatchId) =
    batch_row_range(plan.fidx, id) .- plan.row_offset

out_file_path(outdir::AbstractString, key::AbstractString) =
    joinpath(outdir, "arMADGICS_out_" * key * ".h5")

"""
    plan_workup(rawdir; fibers=1:600, allow_missing=false, batchsize=100, log=...)

W1 startup: load identity, assert the ordering contract, build the expected
batch set (restricted to the contiguous fiber window), discover + reconcile,
and pull key schema + `hdr` provenance from the first present batch.
"""
function plan_workup(rawdir::AbstractString, outdir::AbstractString;
        fibers::AbstractUnitRange = 1:RowContract.NFIBER_MAX,
        allow_missing::Bool = false, batchsize::Int = RowContract.DEFAULT_BATCHSIZE,
        log = default_log)
    rawdir = abspath(expanduser(rawdir))
    outdir = abspath(expanduser(outdir))
    f1, f2 = first(fibers), last(fibers)
    (1 <= f1 <= f2 <= RowContract.NFIBER_MAX) ||
        error("fiber window $fibers outside 1:$(RowContract.NFIBER_MAX)")

    log("Loading identity source full_list_info.h5 ...")
    fli = load_full_list_info(joinpath(rawdir, "full_list_info.h5"))
    fidx = verify_fiber_blocks(fli; batchsize = batchsize)
    log("  nsamp(all fibers) = $(fli.nsamp); ordering contract VERIFIED")

    row_offset = fidx.offsets[f1]
    nsamp_out = (fidx.offsets[f2] + fidx.counts[f2]) - row_offset
    nsamp_out > 0 || error("fiber window $fibers contains no rows")

    expected = [id for id in expected_batches(fidx) if f1 <= id.adjfiberindx <= f2]
    discovered = discover_batches(rawdir; fibers = f1:f2)
    rec = reconcile_batches(discovered, expected; allow_missing = allow_missing)
    log("  fiber window $f1:$f2 → $nsamp_out rows, " *
        "$(length(rec.present)) present / $(length(expected)) expected batches" *
        (isempty(rec.missing) ? "" : " ($(length(rec.missing)) MISSING, allow_missing)"))
    isempty(rec.present) && error("no batch files present in fiber window $fibers")

    refbatch = rec.present[1]
    keyinfo = discover_keys(rec.paths[refbatch])
    hdr = read_hdr_attrs(rec.paths[refbatch])
    log("  reference batch $(batch_filename(refbatch)): $(length(keyinfo)) keys; " *
        "hdr attrs: $(join(sort(collect(keys(hdr))), ", "))")

    return WorkupPlan(rawdir, outdir, f1:f2, batchsize, allow_missing,
        nsamp_out, row_offset, fidx, rec.present, rec.missing, rec.paths,
        keyinfo, hdr, refbatch)
end

default_log(s) = (println("[$(now())] $s"); flush(stdout))

# ------------------------------------------------------------------- outputs

"""
Create (or, on resume, reopen and re-verify) the per-key output files with
preallocated datasets sized (leading..., nsamp_out). Layout matches historical
workup.jl: dataset named after the key, plus a `hdr` string dataset carrying
the producer git-provenance attrs. Workup provenance goes in file-level attrs.
"""
function open_outputs(plan::WorkupPlan; resume::Bool = false)
    mkpath(plan.outdir)
    files = Dict{String, HDF5.File}()
    for (k, v) in plan.keyinfo
        p = out_file_path(plan.outdir, k)
        dshape = (v.shape[1:(end - 1)]..., plan.nsamp_out)
        if resume
            isfile(p) || error("--resume: expected output file missing: $p")
            f = h5open(p, "r+")
            haskey(f, k) || (close(f); error("--resume: $p lacks dataset '$k'"))
            size(f[k]) == dshape || (close(f);
                error("--resume: $p dataset '$k' has size $(size(f[k])) != $dshape"))
            files[k] = f
        else
            f = h5open(p, "w")
            write(f, "hdr", HDR_PLACEHOLDER)
            for (a, val) in plan.hdr
                write_attribute(f["hdr"], a, val)
            end
            attrs(f)["workup_rawdir"] = plan.rawdir
            attrs(f)["workup_fibers"] = [first(plan.fibers), last(plan.fibers)]
            attrs(f)["workup_row_offset"] = plan.row_offset
            attrs(f)["workup_nsamp"] = plan.nsamp_out
            attrs(f)["workup_generated"] = string(now())
            attrs(f)["workup_writer"] = "workup_serial.jl"
            create_dataset(f, k, datatype(v.eltype), dataspace(dshape))
            if plan.allow_missing
                msk = zeros(UInt8, plan.nsamp_out)
                for id in plan.missing
                    msk[out_row_range(plan, id)] .= 0x1
                end
                write(f, "missing_row_mask", msk)
            end
            files[k] = f
        end
    end
    return files
end

"Sentinel-fill the rows of every missing batch (allow_missing mode)."
function write_missing_fills!(files, plan::WorkupPlan)
    isempty(plan.missing) && return
    fs = fill_spec(plan.keyinfo)
    for id in plan.missing
        rng = out_row_range(plan, id)
        for (k, v) in plan.keyinfo
            lead = v.shape[1:(end - 1)]
            block = fill(v.eltype(fs[k]), lead..., length(rng))
            _write_slab!(files[k][k], rng, block)
        end
    end
end

@inline function _write_slab!(dset::HDF5.Dataset, rng::UnitRange{Int}, arr)
    nd = ndims(dset)
    idx = ntuple(_ -> Colon(), nd - 1)
    dset[idx..., rng] = arr
    return nothing
end

"Slab-write one batch's data at its contract-derived output row range."
function write_batch!(files, plan::WorkupPlan, id::BatchId, data::Dict{String, Any})
    rng = out_row_range(plan, id)
    for (k, v) in plan.keyinfo
        haskey(data, k) || error("batch $(batch_filename(id)): reader returned no data for key $k")
        _write_slab!(files[k][k], rng, data[k])
    end
    return length(rng)
end

# ------------------------------------------------------------------- readers

"Integrity-check then read ALL keys of one batch file."
function read_batch(path::AbstractString, id::BatchId, nrow_expect::Int, ref_keyinfo)
    problems = check_batch_integrity(path, nrow_expect; ref_keyinfo = ref_keyinfo)
    isempty(problems) ||
        error("batch $(batch_filename(id)) FAILED integrity check:\n    " *
              join(problems, "\n    "))
    data = Dict{String, Any}()
    h5open(path, "r") do f
        for k in keys(ref_keyinfo)
            data[k] = read(f[k])
        end
    end
    return data
end

"""
Worker loop: pull (id, path, nrow_expect) jobs, integrity-check + read, push
(:ok, id, data) or (:err, id, msg) results. Runs until the jobs channel is
closed and drained.
"""
function reader_loop(jobs::RemoteChannel, results::RemoteChannel, ref_keyinfo)
    while true
        job = try
            take!(jobs)
        catch
            break   # channel closed + drained → done
        end
        id, path, nrow = job
        try
            data = read_batch(path, id, nrow, ref_keyinfo)
            put!(results, (:ok, id, data))
        catch e
            put!(results, (:err, id, sprint(showerror, e)))
        end
    end
    return nothing
end

"take! that also watches reader tasks, so a dead reader can't deadlock the writer."
function _take_monitored!(results::RemoteChannel, rtasks)
    while true
        isready(results) && return take!(results)
        if all(istaskdone, rtasks)
            isready(results) && return take!(results)
            for t in rtasks
                istaskfailed(t) && fetch(t)   # rethrows the reader failure
            end
            error("all reader tasks exited but results are incomplete")
        end
        sleep(0.02)
    end
end

# ------------------------------------------------------------------ main run

"""
    run_workup(rawdir, outdir; fibers=1:600, allow_missing=false, batchsize=100,
               resume=false, ckpt_every=200, progress_every=100, log=...)

The W2 serial-tier workup. Uses all currently-added Distributed workers as
readers (inline single-process mode when there are none). Returns the plan.

Checkpointing: after every `ckpt_every` written batches all output files are
flushed and the written batch ids appended to `<outdir>/workup_serial.ckpt`;
`resume=true` skips those batches (rewrites are idempotent — slabs land at
identity-derived ranges — so a crash between data flush and checkpoint append
costs only rework, never corruption).
"""
function run_workup(rawdir::AbstractString, outdir::AbstractString;
        fibers::AbstractUnitRange = 1:RowContract.NFIBER_MAX,
        allow_missing::Bool = false, batchsize::Int = RowContract.DEFAULT_BATCHSIZE,
        resume::Bool = false, ckpt_every::Int = 200, progress_every::Int = 100,
        log = default_log)
    t0 = time()
    plan = plan_workup(rawdir, outdir; fibers = fibers,
        allow_missing = allow_missing, batchsize = batchsize, log = log)

    ckptp = joinpath(plan.outdir, CKPT_NAME)
    done = Set{BatchId}()
    do_resume = resume && isfile(ckptp)
    if do_resume
        for ln in eachline(ckptp)
            id = parse_batch_filename(strip(ln))
            isnothing(id) || push!(done, id)
        end
        log("Resuming: $(length(done)) batches already written per checkpoint")
    end
    todo = [id for id in plan.present if !(id in done)]

    files = open_outputs(plan; resume = do_resume)
    nbytes_written = 0
    try
        write_missing_fills!(files, plan)

        ckpt_io = open(ckptp, do_resume ? "a" : "w")
        pending = String[]                     # written since last flush
        function checkpoint!(force::Bool = false)
            (force || length(pending) >= ckpt_every) || return
            isempty(pending) && return
            for f in values(files)
                flush(f)
            end
            foreach(l -> println(ckpt_io, l), pending)
            flush(ckpt_io)
            empty!(pending)
        end

        nwritten = 0
        rowbytes = sum(prod(v.shape[1:(end - 1)]) * sizeof(v.eltype)
                       for v in values(plan.keyinfo))
        function after_write!(id, nrows)
            nwritten += 1
            nbytes_written += nrows * rowbytes
            push!(pending, batch_filename(id))
            checkpoint!()
            if nwritten % progress_every == 0 || nwritten == length(todo)
                el = time() - t0
                log(@sprintf("  written %d / %d batches (%.2f batches/s, %.2f MB/s out, ETA %ds)",
                    nwritten, length(todo), nwritten / el, nbytes_written / el / 1e6,
                    round(Int, (length(todo) - nwritten) * el / max(nwritten, 1))))
            end
        end

        if nworkers() > 1
            log("Streaming with $(nworkers()) distributed readers, single writer ...")
            jobs = RemoteChannel(() -> Channel{Tuple{BatchId, String, Int}}(max(length(todo), 1)))
            results = RemoteChannel(() -> Channel{Any}(2 * nworkers()))
            for id in todo
                put!(jobs, (id, plan.paths[id],
                    length(batch_within_range(plan.fidx, id))))
            end
            close(jobs)
            rtasks = [@async remotecall_fetch(reader_loop, w, jobs, results, plan.keyinfo)
                      for w in workers()]
            for _ in 1:length(todo)
                status, id, payload = _take_monitored!(results, rtasks)
                if status != :ok
                    close(results)   # unblocks any reader mid-put!; they exit
                    error("ABORT (nothing further written): $payload")
                end
                nrows = write_batch!(files, plan, id, payload)
                after_write!(id, nrows)
            end
            foreach(wait, rtasks)
        else
            log("Streaming inline (no distributed workers) ...")
            for id in todo
                data = read_batch(plan.paths[id], id,
                    length(batch_within_range(plan.fidx, id)), plan.keyinfo)
                nrows = write_batch!(files, plan, id, data)
                after_write!(id, nrows)
            end
        end
        checkpoint!(true)
        close(ckpt_io)
    finally
        foreach(close, values(files))
    end
    rm(ckptp; force = true)   # complete → checkpoint no longer meaningful

    covered = sum(length(out_row_range(plan, id)) for id in plan.present)
    miss = plan.nsamp_out - covered
    log(@sprintf("DONE: %d keys × %d rows (%d sentinel-filled missing rows) in %.1fs; %.2f GB written to %s",
        length(plan.keyinfo), plan.nsamp_out, miss, time() - t0,
        nbytes_written / 1e9, plan.outdir))
    return plan
end

end # module WorkupSerial
