## RowContract.jl — W1: the row-index contract for the arMADGICS workup
#
# Output row ranges are derived from IDENTITY, never from file-list position.
#
# Producer contract (pipeline.jl, verified against the 2026_05_01 run):
#   * run_lst is filtered per adjfiberindx in 1:600, in increasing fiber order;
#     each fiber's spectra are enumerated with a within-fiber `linear_index`
#     (1-based) and partitioned into contiguous batches of `batchsize` (=100).
#   * full_list_info.h5 stores (sdss_id, tele, mjd, expnum, adjfiberindx) for
#     every spectrum in that exact global ordering; it is the identity source.
#   * batch_info.txt lists the same rows in the same order as
#     "linear_index, tele, mjd, expnum, adjfiberindx" (linear_index is
#     WITHIN-fiber, restarting at 1 for each fiber).
#   * A batch file arMADGICS_fiber_FFF_batch_SSSSSSS.h5 holds the spectra with
#     within-fiber linear_index in startind:(startind+batchsize-1), clipped to
#     the fiber's count. Its global output row range is therefore
#         fiber_offset(F) .+ (S : min(S + batchsize - 1, fiber_count(F)))
#     where fiber_offset(F) = number of rows of all fibers ordered before F.
#
# Nothing here assumes the discovered file list is complete or sorted: the
# expected batch set is built from full_list_info.h5 identity, discovery is
# reconciled against it, and any discrepancy is a hard error unless the caller
# explicitly opts into `allow_missing`, which yields a missing-rows mask and a
# sentinel fill specification instead of silent row shifts.

module RowContract

using HDF5
using Printf

export BatchId, parse_batch_filename, batch_filename, batch_relpath,
    FullListInfo, load_full_list_info, tele_code, tele_name,
    FiberIndex, verify_fiber_blocks, fiber_row_range,
    expected_batches, batch_row_range, batch_within_range, batch_of_linear_index,
    crosscheck_batch_info, BatchInfoCrossCheck,
    discover_batches, reconcile_batches, BatchReconciliation,
    missing_rows_mask, sentinel_for, fill_spec,
    discover_keys, read_hdr_attrs, check_batch_integrity

const DEFAULT_BATCHSIZE = 100
const NFIBER_MAX = 600  # adjfiberindx 1:300 apo, 301:600 lco

# ------------------------------------------------------------------ BatchId --

"""
Identity of one batch file: the fiber it belongs to and the within-fiber
1-based `linear_index` of its first row (the producer's `startind`).
"""
struct BatchId
    adjfiberindx::Int
    startind::Int
end

Base.isless(a::BatchId, b::BatchId) =
    (a.adjfiberindx, a.startind) < (b.adjfiberindx, b.startind)

const BATCH_FNAME_RE = r"^arMADGICS_fiber_(\d{3})_batch_(\d{7})\.h5$"

"""
    parse_batch_filename(name) -> Union{BatchId, Nothing}

Parse `(adjfiberindx, startind)` from a batch file basename. Returns `nothing`
for non-batch files. The producer writes
`arMADGICS_fiber_FFF_batch_SSSSSSS.h5` with FFF = zero-padded adjfiberindx and
SSSSSSS = zero-padded within-fiber linear_index of the batch's first row.
"""
function parse_batch_filename(name::AbstractString)
    m = match(BATCH_FNAME_RE, basename(String(name)))
    isnothing(m) && return nothing
    return BatchId(parse(Int, m.captures[1]), parse(Int, m.captures[2]))
end

batch_filename(id::BatchId) =
    "arMADGICS_fiber_" * lpad(id.adjfiberindx, 3, "0") *
    "_batch_" * lpad(id.startind, 7, "0") * ".h5"

"Relative path of a batch file below the raw dir (fiber subdirectory + name)."
batch_relpath(id::BatchId) =
    joinpath(lpad(id.adjfiberindx, 3, "0"), batch_filename(id))

# ------------------------------------------------------------- FullListInfo --

const TELE_CODES = Dict{String, Int8}("apo" => Int8(1), "lco" => Int8(2))
const TELE_NAMES = Dict{Int8, String}(Int8(1) => "apo", Int8(2) => "lco")

tele_code(t::AbstractString) = TELE_CODES[String(t)]
tele_name(c::Integer) = TELE_NAMES[Int8(c)]

"""
Compact in-memory copy of full_list_info.h5 (the identity source).
`tele` is coded (1=apo, 2=lco) and `mjd` parsed to Int32 to keep ~26.5M rows
small; use `tele_name` to decode.
"""
struct FullListInfo
    nsamp::Int
    tele::Vector{Int8}
    mjd::Vector{Int32}
    expnum::Vector{Int32}
    adjfiberindx::Vector{Int16}
    sdss_id::Vector{Int64}
end

"""
    load_full_list_info(path; chunk=2_000_000) -> FullListInfo

Read full_list_info.h5 (written by pipeline.jl via jldsave; plain HDF5
datasets) chunk-wise, converting the variable-length string columns
(tele, mjd) to compact codes on the fly.
"""
function load_full_list_info(path::AbstractString; chunk::Int = 2_000_000)
    h5open(path, "r") do f
        for k in ("tele", "mjd", "expnum", "adjfiberindx", "sdss_id")
            haskey(f, k) || error("full_list_info file $path missing dataset '$k'")
        end
        n = length(f["tele"])
        for k in ("mjd", "expnum", "adjfiberindx", "sdss_id")
            length(f[k]) == n || error("full_list_info dataset '$k' length $(length(f[k])) != $n")
        end
        tele = Vector{Int8}(undef, n)
        mjd = Vector{Int32}(undef, n)
        expnum = Vector{Int32}(undef, n)
        adjfiberindx = Vector{Int16}(undef, n)
        sdss_id = Vector{Int64}(undef, n)
        dt = f["tele"]; dm = f["mjd"]; de = f["expnum"]; da = f["adjfiberindx"]; ds = f["sdss_id"]
        i = 1
        while i <= n
            j = min(i + chunk - 1, n)
            tchunk = dt[i:j]
            mchunk = dm[i:j]
            @inbounds for (kk, idx) in enumerate(i:j)
                tele[idx] = tele_code(tchunk[kk])
                mjd[idx] = parse(Int32, mchunk[kk])
            end
            expnum[i:j] .= Int32.(de[i:j])
            adjfiberindx[i:j] .= Int16.(da[i:j])
            sdss_id[i:j] .= Int64.(ds[i:j])
            i = j + 1
        end
        return FullListInfo(n, tele, mjd, expnum, adjfiberindx, sdss_id)
    end
end

# --------------------------------------------------------------- FiberIndex --

"""
Per-fiber row accounting derived from (and verified against) the identity
ordering: `counts[f]` rows for adjfiberindx f, starting at global row
`offsets[f] + 1`. Constructed ONLY through `verify_fiber_blocks`, which
asserts the contract rather than assuming it.
"""
struct FiberIndex
    counts::Vector{Int}    # length NFIBER_MAX
    offsets::Vector{Int}   # offsets[f] = rows before fiber f
    nsamp::Int
    batchsize::Int
end

"""
    verify_fiber_blocks(fli::FullListInfo; batchsize=100) -> FiberIndex

The one-time verification that the global ordering defined by
full_list_info.h5 is exactly "fibers in increasing adjfiberindx order, each
fiber's rows contiguous" — the precondition for batches being contiguous
`batchsize`-row slices of that ordering. Hard-errors (with the offending row)
if any fiber's rows are non-contiguous or fibers are out of order. This is an
assertion of the producer contract, not an assumption.
"""
function verify_fiber_blocks(fli::FullListInfo; batchsize::Int = DEFAULT_BATCHSIZE)
    counts = zeros(Int, NFIBER_MAX)
    prev = 0
    @inbounds for i in 1:fli.nsamp
        fib = Int(fli.adjfiberindx[i])
        (1 <= fib <= NFIBER_MAX) ||
            error("row $i: adjfiberindx $fib outside 1:$NFIBER_MAX")
        if fib != prev
            fib > prev || error(
                "identity-ordering contract violated at global row $i: " *
                "adjfiberindx $fib follows $prev (fiber blocks must be " *
                "contiguous and in increasing fiber order)")
            prev = fib
        end
        counts[fib] += 1
    end
    offsets = zeros(Int, NFIBER_MAX)
    acc = 0
    for f in 1:NFIBER_MAX
        offsets[f] = acc
        acc += counts[f]
    end
    acc == fli.nsamp || error("internal accounting error: $acc != $(fli.nsamp)")
    return FiberIndex(counts, offsets, fli.nsamp, batchsize)
end

"Global row range of one fiber (empty range when the fiber has no rows)."
fiber_row_range(fidx::FiberIndex, fib::Int) =
    (fidx.offsets[fib] + 1):(fidx.offsets[fib] + fidx.counts[fib])

# -------------------------------------------------------- expected batches --

"""
    expected_batches(fidx) -> Vector{BatchId}

The EXPECTED batch set implied by identity: for each fiber with rows, batches
start at within-fiber linear_index 1, 1+batchsize, 1+2*batchsize, ...
"""
function expected_batches(fidx::FiberIndex)
    out = BatchId[]
    for f in 1:NFIBER_MAX
        c = fidx.counts[f]
        c == 0 && continue
        for s in 1:fidx.batchsize:c
            push!(out, BatchId(f, s))
        end
    end
    return out
end

"""
    batch_within_range(fidx, id) -> UnitRange

Within-fiber linear_index range covered by batch `id` (clipped to the fiber's
row count). Hard-errors if `id` is not a legal batch of this index (fiber
empty, startind out of range, or startind not on the batch grid).
"""
function batch_within_range(fidx::FiberIndex, id::BatchId)
    (1 <= id.adjfiberindx <= NFIBER_MAX) ||
        error("batch $(batch_filename(id)): adjfiberindx out of range")
    c = fidx.counts[id.adjfiberindx]
    c > 0 || error("batch $(batch_filename(id)): fiber $(id.adjfiberindx) has no rows in the identity list")
    (1 <= id.startind <= c) ||
        error("batch $(batch_filename(id)): startind $(id.startind) outside 1:$c for fiber $(id.adjfiberindx)")
    mod(id.startind - 1, fidx.batchsize) == 0 ||
        error("batch $(batch_filename(id)): startind $(id.startind) not on the batchsize=$(fidx.batchsize) grid")
    return id.startind:min(id.startind + fidx.batchsize - 1, c)
end

"""
    batch_row_range(fidx, id) -> UnitRange

GLOBAL output row range of batch `id`, derived purely from identity:
`fiber_offset .+ batch_within_range`. This function is the whole point of W1.
"""
function batch_row_range(fidx::FiberIndex, id::BatchId)
    wr = batch_within_range(fidx, id)
    off = fidx.offsets[id.adjfiberindx]
    return (off + first(wr)):(off + last(wr))
end

"Batch that contains within-fiber `linear_index` li of fiber `fib`."
batch_of_linear_index(fidx::FiberIndex, fib::Int, li::Int) =
    BatchId(fib, fld(li - 1, fidx.batchsize) * fidx.batchsize + 1)

# --------------------------------------------------- batch_info crosscheck --

struct BatchInfoCrossCheck
    nrows::Int
    n_mismatch::Int
    first_mismatches::Vector{String}   # up to `maxreport` human-readable lines
    header_total_batches::Union{Int, Nothing}
end

"""
    crosscheck_batch_info(path, fli, fidx; maxreport=20) -> BatchInfoCrossCheck

Row-for-row comparison of batch_info.txt against full_list_info.h5:
same row count, and for every row identical (tele, mjd, expnum, adjfiberindx)
plus `linear_index == within-fiber position` implied by the fiber blocks.
Streams the text file; never loads it whole.
"""
function crosscheck_batch_info(path::AbstractString, fli::FullListInfo,
        fidx::FiberIndex; maxreport::Int = 20)
    n_mismatch = 0
    lines = String[]
    header_total = nothing
    row = 0
    open(path, "r") do io
        for ln in eachline(io)
            s = strip(ln)
            isempty(s) && continue
            if startswith(s, "#")
                m = match(r"Total batches:\s*(\d+)", s)
                isnothing(m) || (header_total = parse(Int, m.captures[1]))
                continue
            end
            row += 1
            if row > fli.nsamp
                n_mismatch += 1
                length(lines) < maxreport &&
                    push!(lines, "row $row: batch_info.txt has more rows than full_list_info.h5 ($(fli.nsamp))")
                continue
            end
            parts = split(s, ',')
            if length(parts) != 5
                n_mismatch += 1
                length(lines) < maxreport && push!(lines, "row $row: unparseable line: $s")
                continue
            end
            li = parse(Int, strip(parts[1]))
            tel = strip(parts[2])
            mj = parse(Int32, strip(parts[3]))
            ex = parse(Int32, strip(parts[4]))
            fib = parse(Int, strip(parts[5]))
            exp_li = row - fidx.offsets[fib > 0 && fib <= NFIBER_MAX ? fib : 1]
            ok = (fib == fli.adjfiberindx[row]) &&
                 (tele_code(tel) == fli.tele[row]) &&
                 (mj == fli.mjd[row]) &&
                 (ex == fli.expnum[row]) &&
                 (li == exp_li)
            if !ok
                n_mismatch += 1
                length(lines) < maxreport && push!(lines,
                    "row $row: batch_info=(li=$li,$tel,$mj,$ex,fib=$fib) vs " *
                    "full_list=(li=$exp_li,$(tele_name(fli.tele[row]))," *
                    "$(fli.mjd[row]),$(fli.expnum[row]),fib=$(fli.adjfiberindx[row]))")
            end
        end
    end
    if row != fli.nsamp
        n_mismatch += 1
        push!(lines, "row-count mismatch: batch_info.txt has $row data rows, full_list_info.h5 has $(fli.nsamp)")
    end
    return BatchInfoCrossCheck(row, n_mismatch, lines, header_total)
end

# ------------------------------------------------------ discovery/reconcile --

"""
    discover_batches(rawdir; fibers=1:NFIBER_MAX) -> Dict{BatchId, String}

Enumerate batch files on disk (one readdir per 3-digit fiber subdirectory —
bounded, non-recursive). Files that do not match the batch-name pattern are
ignored. Discovery order is irrelevant by construction: identity comes from
the parsed filename, never from list position. `fibers` restricts discovery
to fiber subdirectories in that range (for fiber-subset workups; the expected
set must be restricted identically or reconciliation will flag extras).
"""
function discover_batches(rawdir::AbstractString; fibers::AbstractUnitRange = 1:NFIBER_MAX)
    out = Dict{BatchId, String}()
    for entry in sort(readdir(rawdir))
        occursin(r"^\d{3}$", entry) || continue
        parse(Int, entry) in fibers || continue
        sub = joinpath(rawdir, entry)
        isdir(sub) || continue
        for fn in readdir(sub)
            id = parse_batch_filename(fn)
            isnothing(id) && continue
            id.adjfiberindx == parse(Int, entry) || error(
                "batch file $fn found in wrong fiber directory $entry")
            haskey(out, id) && error("duplicate batch id for $fn")
            out[id] = joinpath(sub, fn)
        end
    end
    return out
end

struct BatchReconciliation
    present::Vector{BatchId}     # sorted, = expected ∩ discovered
    missing::Vector{BatchId}     # expected but not on disk
    extra::Vector{BatchId}       # on disk but not expected
    paths::Dict{BatchId, String}
end

function _fmt_batch_list(ids::Vector{BatchId}; maxn::Int = 15)
    strs = [batch_filename(id) for id in first(ids, maxn)]
    length(ids) > maxn && push!(strs, "... ($(length(ids) - maxn) more)")
    return join(strs, "\n    ")
end

"""
    reconcile_batches(discovered, expected; allow_missing=false) -> BatchReconciliation

Compare the discovered batch set against the identity-derived expected set.
Any difference is a HARD ERROR with an explicit missing/extra listing, unless
`allow_missing=true`, in which case missing batches are tolerated (recorded in
the result so callers can build a missing-rows mask via `missing_rows_mask`
and fill with sentinels via `fill_spec`). Extra (unexpected) batch files are
always a hard error: they indicate an identity/contract violation, not a
partial run.
"""
function reconcile_batches(discovered::Dict{BatchId, String},
        expected::Vector{BatchId}; allow_missing::Bool = false)
    expset = Set(expected)
    missing_ids = sort([id for id in expected if !haskey(discovered, id)])
    extra_ids = sort([id for id in keys(discovered) if !(id in expset)])
    present = sort([id for id in expected if haskey(discovered, id)])
    if !isempty(extra_ids)
        error("Discovered batch files NOT in the expected (identity-derived) set — " *
              "refusing to proceed ($(length(extra_ids)) extra):\n    " *
              _fmt_batch_list(extra_ids))
    end
    if !isempty(missing_ids) && !allow_missing
        error("Expected batch files MISSING from disk ($(length(missing_ids)) of " *
              "$(length(expected))) and allow_missing=false — refusing to " *
              "proceed (silently skipping them would shift every subsequent " *
              "output row):\n    " * _fmt_batch_list(missing_ids))
    end
    return BatchReconciliation(present, missing_ids, extra_ids, discovered)
end

"""
    missing_rows_mask(fidx, missing) -> BitVector

Length-nsamp mask, `true` where the global output row belongs to a missing
batch (rows the workup must fill with sentinels and flag).
"""
function missing_rows_mask(fidx::FiberIndex, missing::Vector{BatchId})
    msk = falses(fidx.nsamp)
    for id in missing
        msk[batch_row_range(fidx, id)] .= true
    end
    return msk
end

"Sentinel fill value for missing rows of a dataset of element type T."
sentinel_for(::Type{T}) where {T <: AbstractFloat} = T(NaN)
sentinel_for(::Type{T}) where {T <: Signed} = typemin(T)
sentinel_for(::Type{T}) where {T <: Unsigned} = typemax(T)
sentinel_for(::Type{Bool}) = false

"""
    fill_spec(keyinfo) -> Dict{String, Any}

Sentinel fill specification (key => fill value) for every discovered key, for
`allow_missing` workups: missing rows are written with these values AND
flagged in the `missing_row_mask` dataset the workup writer must emit.
"""
fill_spec(keyinfo::AbstractDict) =
    Dict{String, Any}(k => sentinel_for(v.eltype) for (k, v) in keyinfo)

# -------------------------------------------------------- keys / hdr attrs --

"""
    discover_keys(batchfile) -> Dict{String, NamedTuple}

Discover ALL datasets (key => (shape, eltype)) from one batch file — no
hardcoded key tuple; scalars, RV-level matrices and {n,8700} arrays alike.
The `hdr` provenance dataset is excluded (it is metadata, handled by
`read_hdr_attrs`). Shapes are as seen by Julia/HDF5.jl (column-major): the
batch-row axis is the LAST axis.
"""
function discover_keys(batchfile::AbstractString)
    out = Dict{String, NamedTuple{(:shape, :eltype), Tuple{Tuple{Vararg{Int}}, DataType}}}()
    h5open(batchfile, "r") do f
        for k in keys(f)
            k == "hdr" && continue
            obj = f[k]
            obj isa HDF5.Dataset || continue
            out[k] = (shape = size(obj), eltype = eltype(obj))
        end
    end
    return out
end

"""
    read_hdr_attrs(batchfile) -> Dict{String, Any}

Read the git-provenance attributes stored ON the `hdr` DATASET (not file-level
attributes — the mistake in workup_streaming.py). These must be copied onto
every workup output file.
"""
function read_hdr_attrs(batchfile::AbstractString)
    out = Dict{String, Any}()
    h5open(batchfile, "r") do f
        haskey(f, "hdr") || error("$batchfile has no 'hdr' dataset")
        for a in keys(attributes(f["hdr"]))
            out[a] = read_attribute(f["hdr"], a)
        end
    end
    return out
end

"""
    check_batch_integrity(path, id, fidx; ref_keyinfo=nothing) -> Vector{String}

Per-batch integrity check before any workup write: file openable, has `hdr`,
every dataset's LAST axis equals the identity-derived batch length, and (when
`ref_keyinfo` from a reference batch is given) the key set and leading axes
match. Returns a list of problem strings (empty = clean).
"""
check_batch_integrity(path::AbstractString, id::BatchId,
    fidx::FiberIndex; ref_keyinfo = nothing) =
    check_batch_integrity(path, length(batch_within_range(fidx, id));
        ref_keyinfo = ref_keyinfo)

"""
    check_batch_integrity(path, nrow_expect::Int; ref_keyinfo=nothing) -> Vector{String}

Method taking the identity-derived batch length directly (for callers — e.g.
distributed readers — that carry the expected length rather than the whole
`FiberIndex`).
"""
function check_batch_integrity(path::AbstractString, nrow_expect::Int;
        ref_keyinfo = nothing)
    problems = String[]
    local ki
    try
        ki = discover_keys(path)
        h5open(path, "r") do f
            haskey(f, "hdr") || push!(problems, "missing 'hdr' provenance dataset")
        end
    catch e
        push!(problems, "unreadable: $(sprint(showerror, e))")
        return problems
    end
    for (k, v) in ki
        nlast = isempty(v.shape) ? -1 : v.shape[end]
        nlast == nrow_expect || push!(problems,
            "key $k: last-axis length $nlast != expected batch length $nrow_expect")
    end
    if !isnothing(ref_keyinfo)
        for (k, rv) in ref_keyinfo
            if !haskey(ki, k)
                push!(problems, "missing key $k")
            elseif ki[k].eltype != rv.eltype
                push!(problems, "key $k: eltype $(ki[k].eltype) != reference $(rv.eltype)")
            elseif ki[k].shape[1:end-1] != rv.shape[1:end-1]
                push!(problems, "key $k: leading shape $(ki[k].shape[1:end-1]) != reference $(rv.shape[1:end-1])")
            end
        end
        for k in keys(ki)
            haskey(ref_keyinfo, k) || push!(problems, "unexpected key $k")
        end
    end
    return problems
end

end # module RowContract
