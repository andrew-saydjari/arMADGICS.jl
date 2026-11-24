# compare_h5.jl -- bit-level comparison of a rebuilt release table against a reference.
#
# For every dataset present in either file: presence, dtype, shape, chunk shape, filter pipeline,
# and the VALUES compared bitwise (raw bytes; NaN payloads and -0.0 count), optionally on a row
# slice. Mismatching float columns report the count and the maximum difference in ulps. Every
# dataset attribute and root attribute is compared by value AND by HDF5 string flavour (fixed vs
# variable length, size, padding, character set). Exit status 0 only if everything matches
# (attributes named with --ignore-attr are reported but not counted).
#
# CHAIN POSITION: check (step 20); usable on any pair of release HDF5 files.
# USAGE:   julia -t 4 --project=scripts/release_tables compare_h5.jl REFERENCE CANDIDATE
#              [--rows START:STOP] [--ignore-attr NAME ...] [--only-datasets a,b,...] [--json OUT]
#              [--no-chunk-shortcut]
# INPUTS:  two HDF5 files (read-only). OUTPUT: report on stdout (same lines as compare_h5.py),
#          optional JSON (same structure as compare_h5.py's), exit status 0 / 1.
#
# Julia port of legacy_python/compare_h5.py (2026-11-24). Same checks, same verdicts, same stdout lines, same JSON
# fields; checked against the .py on an equal pair and on a differing pair (allExposure 0.4.0 vs
# the legacy rebuild). Deliberate differences, all of them stricter or faster, never looser:
#   * dtype equality is HDF5's H5Tequal (also sees byte order / string pad / cset), not numpy's.
#   * layout equality compares the chunk shape and the full filter pipeline (filter ids and their
#     parameters, in order; SHUFFLE's automatic element-size parameter excluded), not only h5py's (chunks, compression, compression_opts, shuffle);
#     a difference h5py cannot see is reported in an extra "filters" field.
#   * chunk shortcut: when both datasets have the same dtype, shape, chunking and filters, a chunk
#     whose stored (compressed) bytes and filter mask are identical in both files holds identical
#     values, so it is not decompressed. Only chunks that differ (or are partly outside --rows,
#     or unallocated) are decoded and compared row by row exactly as the .py does, so the counts,
#     first mismatch row and max ulp are the same. --no-chunk-shortcut decodes everything.
#   * datasets are compared in parallel (one task per dataset; HDF5 calls themselves are serialized
#     by HDF5.jl's library lock, the chunk-byte reads and comparisons are not).
#   * variable-length (e.g. vlen string) datasets, on which the .py's byte view fails, are compared
#     by value; scalar datasets count as one row.
# Values are read in blocks of STEP = 4,000,000 rows (first dimension), as in the .py.
using HDF5, ArgParse, JSON
include(joinpath(@__DIR__, "..", "lib", "h5attrs.jl"))

const STEP = 4_000_000
const A = HDF5.API
const Rec = JSON.Object{String,Any}          # insertion-ordered record (Python dict order)

# ----------------------------------------------------------------------------------------------
# Python repr, for stdout lines identical to compare_h5.py's
pyrepr(x::Bool) = x ? "True" : "False"
pyrepr(x::Integer) = string(x)
pyrepr(::Nothing) = "None"
function pyrepr(s::AbstractString)
    q = (occursin('\'', s) && !occursin('"', s)) ? '"' : '\''
    b = IOBuffer(); print(b, q)
    for c in s
        if c == q || c == '\\'
            print(b, '\\', c)
        elseif c == '\n'
            print(b, "\\n")
        else
            print(b, c)
        end
    end
    print(b, q); String(take!(b))
end
pyrepr(v::AbstractVector) = "[" * join((pyrepr(x) for x in v), ", ") * "]"
pyrepr(t::Tuple) = length(t) == 1 ? "(" * pyrepr(t[1]) * ",)" : "(" * join((pyrepr(x) for x in t), ", ") * ")"
pyrepr(d::Union{AbstractDict,AbstractVector{<:Pair}}) =
    "{" * join((pyrepr(k) * ": " * pyrepr(v) for (k, v) in d), ", ") * "}"

# ----------------------------------------------------------------------------------------------
# dataset metadata

"numpy's str(dtype) for an HDF5 type (for reports only; equality uses H5Tequal)."
function np_dtype_str(tid)
    cls = A.h5t_get_class(tid); sz = Int(A.h5t_get_size(tid))
    if cls == A.H5T_INTEGER
        u = A.h5t_get_sign(tid) == A.H5T_SGN_NONE
        sz == 1 && return u ? "uint8" : "int8"
        A.h5t_get_order(tid) == A.H5T_ORDER_BE && return ">" * (u ? "u" : "i") * string(sz)
        return (u ? "uint" : "int") * string(8sz)
    elseif cls == A.H5T_FLOAT
        A.h5t_get_order(tid) == A.H5T_ORDER_BE && return ">f" * string(sz)
        return "float" * string(8sz)
    elseif cls == A.H5T_STRING
        return A.h5t_is_variable_str(tid) ? "object" : "|S" * string(sz)
    elseif cls == A.H5T_ENUM
        return "enum" * string(sz)
    end
    return "class" * string(Int(cls)) * "_size" * string(sz)
end

struct DsMeta
    tid::A.hid_t                     # file datatype (closed by close_meta)
    shape::Vector{Int}               # C order
    chunks::Union{Nothing,Vector{Int}}
    filters::Vector{Tuple{Int,Vector{Int}}}   # (filter id, cd_values) in pipeline order
end

function ds_meta(d::HDF5.Dataset)
    tid = A.h5d_get_type(d)
    sp = A.h5d_get_space(d)
    shape = A.h5s_get_simple_extent_type(sp) == A.H5S_SCALAR ? Int[] : Int.(A.h5s_get_simple_extent_dims(sp, nothing))
    A.h5s_close(sp)
    pl = A.h5d_get_create_plist(d)
    chunks = nothing
    if A.h5p_get_layout(pl) == A.H5D_CHUNKED
        buf = zeros(A.hsize_t, length(shape))
        A.h5p_get_chunk(pl, length(shape), buf)
        chunks = Int.(buf)
    end
    filters = Tuple{Int,Vector{Int}}[]
    for i in 0:(A.h5p_get_nfilters(pl) - 1)
        flags = Ref{Cuint}(0); nel = Ref{Csize_t}(16); cd = zeros(Cuint, 16)
        name = zeros(Cchar, 256); cfg = Ref{Cuint}(0)
        fid = A.h5p_get_filter(pl, i, flags, nel, cd, 256, name, cfg)
        # SHUFFLE's parameter is the element size, filled in by HDF5 from the dtype: not a layout choice
        push!(filters, (Int(fid), fid == 2 ? Int[] : Int.(cd[1:min(Int(nel[]), 16)])))
    end
    A.h5p_close(pl)
    return DsMeta(tid, shape, chunks, filters)
end
close_meta(m::DsMeta) = A.h5t_close(m.tid)

"h5py's (chunks, compression, compression_opts, shuffle), as compare_h5.py prints it."
function h5py_layout_str(m::DsMeta)
    comp = nothing; opts = nothing; shuffle = false
    for (id, cd) in m.filters
        if id == 1
            comp = "gzip"; opts = isempty(cd) ? nothing : cd[1]
        elseif id == 32000
            comp = "lzf"
        elseif id == 4
            comp = "szip"
        elseif id == 2
            shuffle = true
        end
    end
    ch = m.chunks === nothing ? "None" : pyrepr(Tuple(m.chunks))
    return "(" * ch * ", " * pyrepr(comp) * ", " * (comp == "szip" ? "(szip)" : pyrepr(opts)) * ", " * pyrepr(shuffle) * ")"
end

# ----------------------------------------------------------------------------------------------
# raw row-block reads (first dimension = rows, C order)

"Memory type for raw reads: native for numbers, the file type itself otherwise. Caller closes."
memtype_for(tid) = (A.h5t_get_class(tid) == A.H5T_INTEGER || A.h5t_get_class(tid) == A.H5T_FLOAT) ?
                   A.h5t_get_native_type(tid, A.H5T_DIR_ASCEND) : A.h5t_copy(tid)

"Read rows [s, e) (0-based, half-open) of `d` into `buf` as raw bytes of `mt`."
function read_rows!(buf::Vector{UInt8}, d::HDF5.Dataset, mt, shape::Vector{Int}, s::Int, e::Int, rowbytes::Int)
    nb = (e - s) * rowbytes
    length(buf) < nb && resize!(buf, nb)
    if isempty(shape)                                  # scalar dataset
        A.h5d_read(d, mt, A.H5S_ALL, A.H5S_ALL, A.H5P_DEFAULT, buf)
        return buf
    end
    fsp = A.h5d_get_space(d)
    start = zeros(A.hsize_t, length(shape)); start[1] = s
    cnt = A.hsize_t.(shape); cnt[1] = e - s
    A.h5s_select_hyperslab(fsp, A.H5S_SELECT_SET, start, C_NULL, cnt, C_NULL)
    msp = A.h5s_create_simple(length(shape), cnt, C_NULL)
    A.h5d_read(d, mt, msp, fsp, A.H5P_DEFAULT, buf)
    A.h5s_close(msp); A.h5s_close(fsp)
    return buf
end

memeq(a::Vector{UInt8}, ia::Int, b::Vector{UInt8}, ib::Int, n::Int) =
    GC.@preserve a b ccall(:memcmp, Cint, (Ptr{UInt8}, Ptr{UInt8}, Csize_t),
                           pointer(a, ia), pointer(b, ib), n) == 0

# compare_h5.py's ulps(): int view, negatives mapped to INT64_MIN - i, |ia - ib| (wrapping int64)
@inline ulpkey(x::Float64) = (i = reinterpret(Int64, x); i < 0 ? typemin(Int64) - i : i)
@inline ulpkey(x::Float32) = (i = Int64(reinterpret(Int32, x)); i < 0 ? typemin(Int64) - i : i)
@inline ulpkey(x::Float16) = (i = Int64(reinterpret(Int16, x)); i < 0 ? typemin(Int64) - i : i)
@inline ulpdist(a::T, b::T) where {T<:AbstractFloat} = abs(ulpkey(a) - ulpkey(b))

mutable struct Tally
    nbad::Int
    first::Int          # -1 = none
    maxulp::Int
end

"Compare `nrows` rows held in ba/bb (row r starts at byte (r-1)*rowbytes + 1); row0 = file row of r=1."
function tally_block!(t::Tally, ba::Vector{UInt8}, bb::Vector{UInt8}, nrows::Int, rowbytes::Int,
                      row0::Int, ::Type{T}) where {T}
    for r in 1:nrows
        o = (r - 1) * rowbytes + 1
        memeq(ba, o, bb, o, rowbytes) && continue
        t.nbad += 1
        t.first < 0 && (t.first = row0 + r - 1)
        if T <: AbstractFloat
            ne = rowbytes ÷ sizeof(T)
            GC.@preserve ba bb begin
                pa = Ptr{T}(pointer(ba, o)); pb = Ptr{T}(pointer(bb, o))
                for j in 1:ne
                    x = unsafe_load(pa, j); y = unsafe_load(pb, j)
                    (isnan(x) && isnan(y)) && continue
                    u = ulpdist(x, y)
                    u > t.maxulp && (t.maxulp = u)
                end
            end
        end
    end
    return t
end

function float_type(mt)
    A.h5t_get_class(mt) == A.H5T_FLOAT || return Nothing
    sz = A.h5t_get_size(mt)
    return sz == 8 ? Float64 : sz == 4 ? Float32 : sz == 2 ? Float16 : Nothing
end

"Decode and compare rows [s, e) in blocks of STEP rows."
function decode_compare!(t::Tally, da, db, mt, shape, s::Int, e::Int, rowbytes::Int, ::Type{T},
                         ba::Vector{UInt8}, bb::Vector{UInt8}) where {T}
    for bs in s:STEP:(e - 1)
        be = min(bs + STEP, e)
        read_rows!(ba, da, mt, shape, bs, be, rowbytes)
        read_rows!(bb, db, mt, shape, bs, be, rowbytes)
        tally_block!(t, ba, bb, be - bs, rowbytes, bs, T)
    end
    return t
end

"Variable-length datasets: compare values (1-D) through HDF5.jl."
function vlen_compare!(t::Tally, da, db, lo::Int, hi::Int)
    for bs in lo:STEP:(hi - 1)
        be = min(bs + STEP, hi)
        a = da[(bs + 1):be]; b = db[(bs + 1):be]
        for i in eachindex(a)
            isequal(a[i], b[i]) && continue
            t.nbad += 1; t.first < 0 && (t.first = bs + i - 1)
        end
    end
    return t
end

# ----------------------------------------------------------------------------------------------
# chunk shortcut

"Row ranges (0-based, half-open, ascending) inside [lo, hi) that must be decoded."
function rows_to_decode(da, db, fa::String, fb::String, m::DsMeta, lo::Int, hi::Int)
    n = m.shape[1]; cr = m.chunks[1]
    ia = HDF5.get_chunk_info_all(da); ib = HDF5.get_chunk_info_all(db)
    # group chunks by their row offset (Julia-order offsets: the C first dim is the last entry)
    ga = Dict{Int,Vector{Any}}(); gb = Dict{Int,Dict{Any,Any}}()
    for c in ia; push!(get!(() -> Any[], ga, c.offset[end]), c); end
    for c in ib; get!(() -> Dict{Any,Any}(), gb, c.offset[end])[c.offset] = c; end
    nother = isempty(m.shape) ? 1 : prod(cld.(m.shape[2:end], m.chunks[2:end]); init = 1)
    out = Tuple{Int,Int}[]
    bufa = UInt8[]; bufb = UInt8[]
    open(fa, "r") do ioa
        open(fb, "r") do iob
            r0 = (lo ÷ cr) * cr
            for r in r0:cr:(hi - 1)
                rs = max(r, lo); re = min(r + cr, n, hi)
                rs >= re && continue
                full = rs == r && re == min(r + cr, n)
                eq = full && _group_equal(get(ga, r, Any[]), get(gb, r, nothing), nother, ioa, iob, bufa, bufb)
                eq && continue
                if !isempty(out) && out[end][2] == rs
                    out[end] = (out[end][1], re)
                else
                    push!(out, (rs, re))
                end
            end
        end
    end
    return out
end

function _group_equal(ca::Vector{Any}, cb, nother::Int, ioa, iob, bufa, bufb)
    (cb === nothing || length(ca) != nother || length(cb) != nother) && return false
    for c in ca
        d = get(cb, c.offset, nothing)
        d === nothing && return false
        (c.size == d.size && c.filter_mask == d.filter_mask) || return false
        sz = Int(c.size)
        length(bufa) < sz && (resize!(bufa, sz); resize!(bufb, sz))
        GC.@preserve bufa bufb begin
            seek(ioa, c.addr); unsafe_read(ioa, pointer(bufa), sz)
            seek(iob, d.addr); unsafe_read(iob, pointer(bufb), sz)
        end
        memeq(bufa, 1, bufb, 1, sz) || return false
    end
    return true
end

# ----------------------------------------------------------------------------------------------
# one dataset

function compare_dataset(da::HDF5.Dataset, db::HDF5.Dataset, fa::String, fb::String,
                         rows::Union{Nothing,Tuple{Int,Int}}, shortcut::Bool)
    ma = ds_meta(da); mb = ds_meta(db)
    try
        r = Rec()
        r["dtype"] = [np_dtype_str(ma.tid), np_dtype_str(mb.tid)]
        r["shape"] = [ma.shape, mb.shape]
        ok = A.h5t_equal(ma.tid, mb.tid) && ma.shape == mb.shape
        la, lb = h5py_layout_str(ma), h5py_layout_str(mb)
        leq = ma.chunks == mb.chunks && ma.filters == mb.filters
        r["layout_equal"] = leq
        if !leq
            r["layout"] = [la, lb]
            la == lb && (r["filters"] = [string(ma.filters), string(mb.filters)])
        end
        if !ok
            r["values_equal"] = false
            return r
        end
        n = isempty(ma.shape) ? 1 : ma.shape[1]
        lo, hi = rows === nothing ? (0, n) : rows
        hi = min(hi, n)
        isempty(ma.shape) && ((lo, hi) = (0, 1))
        t = Tally(0, -1, 0)
        isvlen = A.h5t_get_class(ma.tid) == A.H5T_VLEN ||
                 (A.h5t_get_class(ma.tid) == A.H5T_STRING && A.h5t_is_variable_str(ma.tid))
        if isvlen
            vlen_compare!(t, da, db, lo, hi)
        else
            mt = memtype_for(ma.tid)
            try
                rowbytes = Int(A.h5t_get_size(mt)) * prod(ma.shape[2:end]; init = 1)
                T = float_type(mt)
                ranges = (shortcut && leq && ma.chunks !== nothing && !isempty(ma.shape)) ?
                         rows_to_decode(da, db, fa, fb, ma, lo, hi) : [(lo, hi)]
                ba = UInt8[]; bb = UInt8[]
                for (s, e) in ranges
                    s < e && decode_compare!(t, da, db, mt, ma.shape, s, e, rowbytes, T, ba, bb)
                end
            finally
                A.h5t_close(mt)
            end
        end
        r["rows_compared"] = hi - lo
        r["values_equal"] = t.nbad == 0
        if t.nbad > 0
            r["n_mismatch"] = t.nbad
            r["first_mismatch_row"] = t.first
            A.h5t_get_class(ma.tid) == A.H5T_FLOAT && (r["max_ulp"] = t.maxulp)
        end
        return r
    finally
        close_meta(ma); close_meta(mb)
    end
end

# ----------------------------------------------------------------------------------------------
# attributes

"Attribute value with the equality h5py + np.array_equal give it in compare_h5.py."
struct AttrVal
    kind::Symbol            # :vstr (str), :fstr (bytes, trailing NULs stripped), :num, :raw, :null
    shape::Vector{Int}
    vals::Vector{Any}
end

function read_attr_value(obj, name::AbstractString)
    aid = A.h5a_open(obj, name, A.H5P_DEFAULT)
    tid = A.h5a_get_type(aid); sid = A.h5a_get_space(aid)
    try
        st = A.h5s_get_simple_extent_type(sid)
        st == A.H5S_NULL && return AttrVal(:null, Int[], Any[np_dtype_str(tid)])
        shape = Int.(A.h5s_get_simple_extent_dims(sid, nothing))
        n = prod(shape; init = 1)
        cls = A.h5t_get_class(tid)
        if cls == A.H5T_STRING && A.h5t_is_variable_str(tid)
            mt = A.h5t_copy(tid)
            ptrs = fill(Ptr{Cchar}(C_NULL), n)
            A.h5a_read(aid, mt, ptrs)
            vals = Any[p == C_NULL ? "" : unsafe_string(p) for p in ptrs]
            A.h5d_vlen_reclaim(mt, sid, A.H5P_DEFAULT, ptrs)
            A.h5t_close(mt)
            return AttrVal(:vstr, shape, vals)
        elseif cls == A.H5T_STRING
            sz = Int(A.h5t_get_size(tid)); buf = zeros(UInt8, n * sz)
            A.h5a_read(aid, tid, buf)
            vals = Any[]
            for i in 1:n
                b = buf[(i - 1) * sz + 1:i * sz]
                while !isempty(b) && b[end] == 0x00; pop!(b); end
                push!(vals, b)
            end
            return AttrVal(:fstr, shape, vals)
        elseif cls == A.H5T_INTEGER || cls == A.H5T_FLOAT
            mt = A.h5t_get_native_type(tid, A.H5T_DIR_ASCEND)
            sz = Int(A.h5t_get_size(mt))
            T = cls == A.H5T_FLOAT ? (sz == 8 ? Float64 : sz == 4 ? Float32 : Float16) :
                A.h5t_get_sign(mt) == A.H5T_SGN_NONE ?
                    (sz == 8 ? UInt64 : sz == 4 ? UInt32 : sz == 2 ? UInt16 : UInt8) :
                    (sz == 8 ? Int64 : sz == 4 ? Int32 : sz == 2 ? Int16 : Int8)
            buf = Vector{T}(undef, n)
            A.h5a_read(aid, mt, buf)
            A.h5t_close(mt)
            return AttrVal(:num, shape, Any[buf...])
        else
            sz = Int(A.h5t_get_size(tid)); buf = zeros(UInt8, n * sz)
            A.h5a_read(aid, tid, buf)
            return AttrVal(:raw, shape, Any[np_dtype_str(tid), buf])
        end
    finally
        A.h5s_close(sid); A.h5t_close(tid); A.h5a_close(aid)
    end
end

attrs_equal(a::AttrVal, b::AttrVal) =
    a.kind == b.kind && a.shape == b.shape && length(a.vals) == length(b.vals) &&
    all(a.vals[i] == b.vals[i] for i in eachindex(a.vals))

function compare_attrs(oa, ob, where::String, ignore::Set{String})
    out = Rec[]
    ka = Set(keys(attributes(oa))); kb = Set(keys(attributes(ob)))
    for k in sort!(collect(union(ka, kb)))
        rec = Rec("where" => where, "attr" => k)
        if !(k in ka) || !(k in kb)
            rec["problem"] = "missing in " * (k in ka ? "candidate" : "reference")
        elseif !attrs_equal(read_attr_value(oa, k), read_attr_value(ob, k))
            rec["problem"] = "value differs"
        else
            fa_, fb_ = attr_flavour(oa, k), attr_flavour(ob, k)
            fa_ != fb_ && (rec["problem"] = "type differs " * pyrepr(fa_) * " vs " * pyrepr(fb_))
        end
        if haskey(rec, "problem")
            rec["ignored"] = k in ignore
            push!(out, rec)
        end
    end
    return out
end

# ----------------------------------------------------------------------------------------------
# driver

function walk(g, prefix::String = "")
    out = String[]
    for k in keys(g)
        o = g[k]
        if o isa HDF5.Group
            append!(out, walk(o, prefix * k * "/"))
        else
            push!(out, prefix * k)
        end
        close(o)
    end
    return out
end

function parse_cli(args)
    s = ArgParseSettings(description = "bit-level comparison of two HDF5 release tables")
    @add_arg_table! s begin
        "reference"
            required = true
        "candidate"
            required = true
        "--rows"
            default = ""
        "--ignore-attr"
            nargs = '*'; default = String[]
        "--only-datasets"
            default = ""
        "--json"
            default = ""
        "--no-chunk-shortcut"
            action = :store_true
            help = "decode every chunk (the .py's algorithm) instead of skipping byte-identical chunks"
    end
    return parse_args(args, s)
end

function compare_files(ref::String, cand::String; rows = nothing, ignore = Set{String}(),
                       only = nothing, shortcut::Bool = true)
    res = Rec("reference" => ref, "candidate" => cand,
              "rows" => rows === nothing ? "all" : string(rows[1], ":", rows[2]),
              "datasets" => Rec(), "attr_problems" => Rec[])
    allok = true
    h5open(ref, "r") do fa
        h5open(cand, "r") do fb
            ka, kb = Set(walk(fa)), Set(walk(fb))
            if only !== nothing
                ka, kb = intersect(ka, only), intersect(kb, only)
            end
            names = sort!(collect(union(ka, kb)))
            both = [k for k in names if k in ka && k in kb]
            results = Dict{String,Tuple{Rec,Vector{Rec}}}()
            lk = ReentrantLock()
            jobs = Channel{String}(length(both))
            foreach(k -> put!(jobs, k), both); close(jobs)
            ntask = min(Threads.nthreads(), max(1, length(both)))
            @sync for _ in 1:ntask
                Threads.@spawn for k in jobs
                    da = fa[k]; db = fb[k]
                    r = compare_dataset(da, db, ref, cand, rows, shortcut)
                    ap = compare_attrs(da, db, k, ignore)
                    close(da); close(db)
                    lock(lk) do; results[k] = (r, ap); end
                end
            end
            for k in names
                if !(k in ka) || !(k in kb)
                    res["datasets"][k] = Rec("problem" => "missing in " * (k in ka ? "candidate" : "reference"))
                    allok = false
                    continue
                end
                r, ap = results[k]
                res["datasets"][k] = r
                allok &= r["values_equal"] && r["layout_equal"]
                append!(res["attr_problems"], ap)
            end
            only === nothing && append!(res["attr_problems"], compare_attrs(fa, fb, "/", ignore))
        end
    end
    allok &= !any(p -> !p["ignored"], res["attr_problems"])
    res["all_equal"] = allok
    return res
end

const EXTRA_KEYS = ("n_mismatch", "max_ulp", "first_mismatch_row", "problem", "layout", "filters", "dtype")

function report(io::IO, res)
    for (k, r) in res["datasets"]
        good = get(r, "values_equal", false) && get(r, "layout_equal", false)
        extra = [kk => v for (kk, v) in r if kk in EXTRA_KEYS && !(kk == "dtype" && get(r, "values_equal", false))]
        println(io, (good ? "OK " : "BAD"), " ", rpad(k, 28), " ", isempty(extra) ? "" : pyrepr(extra))
    end
    for p in res["attr_problems"]
        println(io, (p["ignored"] ? "IGN " : "ATTR"), " ", p["where"], " :: ", p["attr"], " : ", p["problem"])
    end
    println(io, res["all_equal"] ? "ALL EQUAL" : "DIFFERENCES FOUND")
end

function main(args)
    a = parse_cli(args)
    rows = isempty(a["rows"]) ? nothing : Tuple(parse.(Int, split(a["rows"], ":")))
    rows !== nothing && length(rows) != 2 && error("--rows must be START:STOP")
    only = isempty(a["only-datasets"]) ? nothing : Set(String.(split(a["only-datasets"], ",")))
    t = @elapsed res = compare_files(a["reference"], a["candidate"]; rows = rows,
                                     ignore = Set{String}(a["ignore-attr"]), only = only,
                                     shortcut = !a["no-chunk-shortcut"])
    report(stdout, res)
    if !isempty(a["json"])
        open(io -> JSON.json(io, res; pretty = 1), a["json"], "w")
    end
    println(stderr, "work time: ", round(t; digits = 2), " s (", Threads.nthreads(), " threads)")
    return res["all_equal"] ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
