# almanac.jl -- stage S_ALMANAC (and the almanac half of S_OFFTARGET).
#
# Parallel read, ordered apply:
#   * read_night (on Distributed workers) reads ONE telescope-night's exposure table and only the
#     fiber tables that the night's workup rows can reference, into small typed tables. Workers
#     open the almanac read-only; nothing is written in parallel. Transfer per night is tiny
#     (exposure table + <= a few 300-row fiber tables), never per-row data.
#   * apply_night! (on the head process) expands those tables onto the night's rows exactly as
#     the 0.4.0 s6 join did, in sorted (tele, mjd) order so the categorical code tables get the
#     same first-seen order as 0.4.0.

const EXP_INT = ("config_id", "plate_id", "design_id", "field_id", "cart_id", "map_id",
                 "n_read", "fps", "flagged_bad", "chip_flags")
const EXP_FLT = ("seeing", "focus", "collpist", "colpitch", "dithered_pixels")
const FIB_INT = ("sdss_id", "catalogid")
const FIB_FLT = ("ra", "dec", "lambda_eff", "x_focal", "y_focal")

"numpy-assignment casting semantics of the 0.4.0 build (unsafe casts, as numpy does)."
_conv(::Type{T}, x::Bool) where {T<:Integer} = T(x)
_conv(::Type{T}, x::Integer) where {T<:Integer} = x % T
_conv(::Type{T}, x::AbstractFloat) where {T<:Integer} = unsafe_trunc(T, x)
_conv(::Type{T}, x::Real) where {T<:AbstractFloat} = T(x)

"Fixed-length almanac strings as numpy's 'S' dtype yields them (trailing NULs dropped)."
_strs(v) = String[rstrip(String(s), '\0') for s in v]

_readnum(g, c) = haskey(g, c) ? read(g[c]) : nothing
_readstr(g, c) = haskey(g, c) ? _strs(read(g[c])) : nothing

struct FiberTable
    n::Int                                   # length of the "ra" column (0 if absent)
    num::Dict{String,Any}                    # FIB_INT/FIB_FLT name => raw vector (absent => missing key)
    cat::Dict{String,Vector{String}}         # category / hole_type
end

struct NightData
    t::Int
    night::Int
    exposure::Vector{Int64}
    fps::Vector{Int64}
    config_id::Vector{Int64}
    plate_id::Vector{Int64}
    expnum_cols::Dict{String,Any}            # EXP_INT/EXP_FLT raw vectors (absent => missing key)
    image_type::Union{Nothing,Vector{String}}
    has_fibers::Bool
    fiber_keys::Set{String}                  # all keys present under fibers/
    fibers::Dict{String,FiberTable}          # only the keys rows can reference
    fps_night::Bool                          # offtarget: first fibers/ key has on_target (0.4.0 s1)
    fps_flags::Dict{Int32,Vector{Int8}}      # offtarget: config_id => assigned & !on_target per index
end

"Read one telescope-night (worker side). `expnums` = distinct workup expnums of the night."
function read_night(A::HDF5.File, t::Int, night::Int, expnums::Vector{Int32})
    g = A["raw"][telename(t)][string(night)]
    ex = g["exposures"]
    e_exp = Int64.(read(ex["exposure"])); e_fps = Int64.(read(ex["fps"]))
    e_cfg = Int64.(read(ex["config_id"])); e_plt = Int64.(read(ex["plate_id"]))
    cols = Dict{String,Any}()
    for c in (EXP_INT..., EXP_FLT...)
        haskey(ex, c) && (cols[c] = read(ex[c]))
    end
    imt = _readstr(ex, "image_type")
    pos = Dict{Int64,Int}(); for (i, v) in enumerate(e_exp); pos[v] = i; end

    has_f = haskey(g, "fibers")
    fkeys = has_f ? Set(keys(g["fibers"])) : Set{String}()
    want = Set{String}()
    for e in expnums
        i = get(pos, Int64(e), 0); i == 0 && continue
        push!(want, e_fps[i] == 1 ? string(e_cfg[i]) : string(e_plt[i]))
        push!(want, string(e_cfg[i]))        # offtarget looks up by config_id regardless of fps
    end
    fibers = Dict{String,FiberTable}()
    fps_night = false
    flags = Dict{Int32,Vector{Int8}}()
    if has_f
        fg = g["fibers"]
        ks = collect(keys(fg))
        fps_night = !isempty(ks) && haskey(fg[ks[1]], "on_target")
        for k in want
            k in fkeys || continue
            grp = fg[k]
            num = Dict{String,Any}()
            for c in (FIB_INT..., FIB_FLT...)
                haskey(grp, c) && (num[c] = read(grp[c]))
            end
            cat = Dict{String,Vector{String}}()
            for c in CATEGORICAL[1:2]
                haskey(grp, c) && (cat[c] = _strs(read(grp[c])))
            end
            n = haskey(grp, "ra") ? length(grp["ra"]) : 0
            fibers[k] = FiberTable(n, num, cat)
            # offtarget table (0.4.0 offtarget_sidecar/s1 predicate and length check)
            if fps_night && haskey(grp, "on_target")
                asg = read(grp["assigned"]); ont = read(grp["on_target"])
                val = read(grp["valid"]); dis = read(grp["disabled"]); cc = read(grp["category"])
                m = length(asg)
                if length(ont) == m && length(val) == m && length(dis) == m && length(cc) == m
                    flags[Int32(parse(Int, k))] = Int8[(Bool(asg[i]) && !Bool(ont[i])) ? 1 : 0 for i in 1:m]
                end
            end
        end
    end
    return NightData(t, night, e_exp, e_fps, e_cfg, e_plt, cols, imt, has_f, fkeys, fibers, fps_night, flags)
end

"Worker entry: read a batch of nights with one file handle."
function read_nights(path::String, batch::Vector{Tuple{Int,Int,Vector{Int32}}})
    h5open(path, "r") do A
        [read_night(A, t, n, e) for (t, n, e) in batch]
    end
end

"First-seen code table, as the 0.4.0 `cat_maps` dicts (\"\" => 0 preset)."
struct CodeTable
    map::Dict{String,UInt8}
    order::Vector{String}
end
CodeTable() = CodeTable(Dict("" => 0x00), [""])
function code!(ct::CodeTable, vals::AbstractVector{String})
    for v in sort(unique(vals))                       # np.unique order
        haskey(ct.map, v) && continue
        length(ct.order) < 256 || error("categorical column overflows uint8")
        ct.map[v] = UInt8(length(ct.order)); push!(ct.order, v)
    end
    return [ct.map[v] for v in vals]
end

"All columns produced by the almanac stage, preallocated with the 0.4.0 fill values."
struct AlmanacCols
    i64::Dict{String,Vector{Int64}}
    f64::Dict{String,Vector{Float64}}
    i32::Dict{String,Vector{Int32}}
    f32::Dict{String,Vector{Float32}}
    cats::Dict{String,Vector{UInt8}}
    codes::Dict{String,CodeTable}
    join_status::Vector{UInt8}
    era::Vector{UInt8}
    sdss_id_agrees::Vector{UInt8}
    offtarget::Vector{Int8}
end
function AlmanacCols(N::Int)
    AlmanacCols(Dict(c => fill(Int64(-1), N) for c in FIB_INT),
                Dict(c => fill(NaN, N) for c in FIB_FLT),
                Dict(c => fill(Int32(-1), N) for c in EXP_INT),
                Dict(c => fill(NaN32, N) for c in EXP_FLT),
                Dict(c => zeros(UInt8, N) for c in CATEGORICAL),
                Dict(c => CodeTable() for c in CATEGORICAL),
                zeros(UInt8, N), zeros(UInt8, N), zeros(UInt8, N), zeros(Int8, N))
end

"""
Codes for `vals[idx]` with the 0.4.0 first-seen rule: the DISTINCT values among the selected
entries are registered in sorted order (np.unique), then every entry is coded. Works on table
indices so each string is looked up once per table, not once per row.
"""
function code_indexed!(ct::CodeTable, vals::Vector{String}, idx::AbstractVector{<:Integer})
    used = sort!(unique(idx))
    code!(ct, vals[used])                    # registers sorted distinct values (np.unique order)
    tab = Dict{Int,UInt8}(i => ct.map[vals[i]] for i in used)
    return UInt8[tab[i] for i in idx]
end

"Expand one night onto its rows (head side). 0.4.0 s6 join, same order of effects."
function apply_night!(C::AlmanacCols, nd::NightData, rows::Vector{Int}, id::Identity, fia::Vector{Int16})
    pos = Dict{Int64,Int}(); for (i, v) in enumerate(nd.exposure); pos[v] = i; end
    r_ok = Int[]; er = Int[]
    for r in rows
        e = get(pos, Int64(id.expnum[r]), 0)
        if e == 0
            C.join_status[r] = 0x02
        else
            push!(r_ok, r); push!(er, e)
        end
    end
    isempty(r_ok) && return nothing
    for (r, e) in zip(r_ok, er)
        C.era[r] = nd.fps[e] == 1 ? 0x02 : 0x01
    end
    for c in EXP_INT
        haskey(nd.expnum_cols, c) || continue
        src = nd.expnum_cols[c]; dst = C.i32[c]
        for (r, e) in zip(r_ok, er); dst[r] = _conv(Int32, src[e]); end
    end
    for c in EXP_FLT
        haskey(nd.expnum_cols, c) || continue
        src = nd.expnum_cols[c]; dst = C.f32[c]
        for (r, e) in zip(r_ok, er); dst[r] = _conv(Float32, src[e]); end
    end
    if nd.image_type !== nothing
        C.cats["image_type"][r_ok] = code_indexed!(C.codes["image_type"], nd.image_type, er)
    end

    # fiber-table key per exposure (config_id in the FPS era, plate_id otherwise), rows grouped
    # by key in row order; keys visited in sorted STRING order, as np.unique(keystr) did
    ekey = Dict{Int,String}()
    groups = Dict{String,Vector{Int}}()
    for (k, e) in enumerate(er)
        kk = get!(() -> nd.fps[e] == 1 ? string(nd.config_id[e]) : string(nd.plate_id[e]), ekey, e)
        push!(get!(groups, kk, Int[]), k)
    end
    for kk in sort!(collect(keys(groups)))
        rr = r_ok[groups[kk]]
        if !nd.has_fibers || !(kk in nd.fiber_keys)
            C.join_status[rr] .= 0x03
            continue
        end
        ft = nd.fibers[kk]
        if ft.n != 300
            C.join_status[rr] .= 0x04
            continue
        end
        ii = [Int(fia[r]) - 1 for r in rr]             # 0-based almanac fiber index
        bad = [(i < 0) | (i > 299) for i in ii]
        C.join_status[rr[bad]] .= 0x05
        gd = .!bad
        any(gd) || continue
        rr2 = rr[gd]; ii2 = ii[gd] .+ 1
        for c in FIB_INT
            haskey(ft.num, c) || continue
            src = ft.num[c]; dst = C.i64[c]
            for (r, i) in zip(rr2, ii2); dst[r] = _conv(Int64, src[i]); end
        end
        for c in FIB_FLT
            haskey(ft.num, c) || continue
            src = ft.num[c]; dst = C.f64[c]
            for (r, i) in zip(rr2, ii2); dst[r] = _conv(Float64, src[i]); end
        end
        for c in CATEGORICAL[1:2]
            haskey(ft.cat, c) || continue
            C.cats[c][rr2] = code_indexed!(C.codes[c], ft.cat[c], ii2)
        end
        C.join_status[rr2] .= 0x01
        if haskey(ft.num, "sdss_id")
            src = ft.num["sdss_id"]
            for (r, i) in zip(rr2, ii2)
                a = _conv(Int64, src[i]); w = id.sdss_id[r]
                C.sdss_id_agrees[r] = !((a > 0) & (w > 0)) ? 0x03 : (a == w ? 0x01 : 0x02)
            end
        end
    end
    return nothing
end

"Offtarget for one night's rows (0.4.0 offtarget s2): flag = assigned & !on_target at the
row's (config_id, fiber_index_almanac); rows without an FPS table stay 0. `nf` = the survey-wide
maximum fiber-table length (the 0.4.0 matrix height)."
function apply_offtarget!(C::AlmanacCols, nd::NightData, rows::Vector{Int}, fia::Vector{Int16}, nf::Int)
    nd.fps_night || return nothing
    cfg = C.i32["config_id"]
    for r in rows
        v = get(nd.fps_flags, cfg[r], nothing)
        v === nothing && continue
        k = Int(fia[r])
        (1 <= k <= nf) || continue
        C.offtarget[r] = k <= length(v) ? v[k] : Int8(0)
    end
    return nothing
end

"Run the almanac stage: parallel read over workers, ordered apply on the head."
function join_almanac(id::Identity, almanac::AbstractString; nbatch::Int = 0)
    N = length(id)
    fia = fiber_index_almanac.(id.tele, id.adjfiberindx)
    okey = Int64.(id.tele) .* 10^6 .+ Int64.(id.mjd)
    o = sortperm(okey; alg = Base.Sort.DEFAULT_STABLE)
    nights = Tuple{Int,Int}[]; rowsets = Vector{Vector{Int}}()
    b = 1
    while b <= N
        e = b
        while e < N && okey[o[e+1]] == okey[o[b]]; e += 1; end
        k = okey[o[b]]
        push!(nights, (Int(k ÷ 10^6), Int(k % 10^6))); push!(rowsets, o[b:e])
        b = e + 1
    end
    tasks = [(t, n, unique(id.expnum[rs])) for ((t, n), rs) in zip(nights, rowsets)]
    nb = nbatch > 0 ? nbatch : max(1, 8 * max(nworkers(), 1))
    batches = [tasks[i:nb:end] for i in 1:min(nb, length(tasks))]   # interleaved: balances eras
    path = String(almanac)
    t0 = time()
    parts = pmap(bt -> read_nights(path, bt), batches)
    t_read = time() - t0
    nd_by = Dict{Tuple{Int,Int},NightData}()
    for p in parts, nd in p; nd_by[(nd.t, nd.night)] = nd; end

    t0 = time()
    C = AlmanacCols(N)
    for ((t, n), rs) in zip(nights, rowsets)
        apply_night!(C, nd_by[(t, n)], rs, id, fia)
    end
    @info @sprintf("almanac: parallel read %.1f s (%d workers), ordered apply %.1f s", t_read, nworkers(), time() - t0)
    # offtarget needs the survey-wide maximum FPS table length, as the 0.4.0 matrix did
    nf = maximum((length(v) for nd in values(nd_by) for v in values(nd.fps_flags)); init = 0)
    for ((t, n), rs) in zip(nights, rowsets)
        apply_offtarget!(C, nd_by[(t, n)], rs, fia, nf)
    end
    return C, fia, length(nights)
end
