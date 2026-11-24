#!/usr/bin/env julia
#=
build_allTargets.jl (0.4.1) — the allTargets companion table, one row per target (sdss_id).

Copied from 2026_11_19/carton_flags/build_allTargets.jl (the 0.4.0 builder) and changed so that
EVERY dataset of the input `/meta` group is carried into the output automatically. There is no
hard-coded list of columns to read: a column that `almanac add metadata` adds later (e.g. the
`--extended` groups, or anything after them) appears in allTargets without editing this file.

Rows: one per sdss_id > 0 present in the allExposure table, sorted by sdss_id, as in 0.4.0; plus
the solar-system sdss_id(s) (almanac's fixed row, sdss_id 0 = Vesta, `lead` = "solar_system") when
that id is present both in allExposure and in /meta (--include-solar-system, default true). The
0.4.0 inputs contain no such row, so 0.4.0 is reproduced unchanged.

Columns written (output names):
  * sdss_id                       from allExposure (the row key)
  * every /meta dataset, verbatim  ("pass-through"), except those in REPLACED_META below
  * sdss5_target_flags            REBUILT from the (sdss_id, carton_pk) pairs (see below); the /meta
                                  column of the same name is NOT copied (its width is recorded in
                                  the file attribute meta_sdss5_target_flags_nbytes). With
                                  --keep-meta-target-flags true it is also copied, renamed
                                  meta_sdss5_target_flags.
  * n_cartons, n_exposures, metadata_status   derived, exactly as in 0.4.0
  * groups carton_pks/ and carton_table/       exactly as in 0.4.0
If two sources would produce the same output name (a /meta column called n_cartons, carton_pks,
...), the build fails. A 0.4.0 column missing from /meta also fails the build (0.4.1 must be a
superset of 0.4.0).

COLUMN ORDER RULE (HDF5 creation order, tracked in the file, and file attribute `column_order`):
  1. the 44 columns of 0.4.0 in the order of the 0.4.0 product registry (`column_order` in
     allTargets_pass1_dr21.toml) — ORDER_040 below; ORDER_040 fixes ORDER only, it selects nothing;
  2. then every other /meta column, in the order it appears in the /meta file (the /meta group's
     `columns` attribute if present, else its HDF5 creation order);
  3. then meta_sdss5_target_flags if requested; then the groups carton_pks/, carton_table/.
(h5ls and HDF5.jl `keys` list names alphabetically; h5py lists creation order.)

Per-column conventions (same as 0.4.0 for its columns, applied to all):
  * numeric: same element type as in /meta; rows of allExposure ids absent from /meta are filled
    with NaN (floats), -1 (signed ints), 0 (unsigned), "" (strings). /meta's own sentinels
    (NaN, -1, "", 0) are copied unchanged.
  * h5py booleans (int8 enum) are read and written as Int8 0/1 (as 0.4.0).
  * strings: /meta stores fixed-length null-padded strings whose byte width varies between builds;
    they are written as variable-length UTF-8 (as 0.4.0). New columns carry attribute `meta_dtype`
    with the original width.
  * multi-dimensional columns: kept as they are (same per-row shape; only the row axis is
    re-indexed); new columns carry attribute `shape_per_row`. The only one today is
    sdss5_target_flags, which is rebuilt (above).
  * attributes: every attribute of a /meta dataset (today only `description`) is copied to new
    columns, plus `units` parsed from a trailing "[...]" in the description (as 0.4.0), plus
    `source` = "/meta/<name>". The 0.4.0 columns keep exactly the 0.4.0 attributes.

sdss5_target_flags (unchanged from 0.4.0, 2026-11-19 carton_flags): rebuilt from the raw
(sdss_id, carton_pk) pairs with a semaphore mapping CSV (--semaphore-mapping). Per carton_pk:
  1. exact:  the pk has a bit in the mapping -> that bit;
  2. name:   otherwise the bit of a same-name targetdb carton pk that has one (prefer same
             program+mapper, then the highest carton_pk); no new bits are invented;
  3. none:   otherwise not encoded (metadata_status UNKNOWN_CARTON).
Pass --name-fallback false for exact-pk only (astra rule).

Usage (all paths are arguments; nothing points at 0.4.0):
  julia --project=. build_allTargets.jl --release 0.4.1 \
        --allexposure <allExposure.h5> --meta <targets_or_meta.h5> \
        --carton-table <targetdb_carton.csv> \
        --semaphore-mapping <sdss5_target_4_with_groups.csv> --semaphore-version 0.2.9 \
        --out <allTargets.h5>
=#
using HDF5, ArgParse, JSON, Dates, SHA

const STATUS_BITS = [
    (0, "NOT_QUERIED",    "sdss_id not yet queried (partial build only)"),
    (1, "NOT_FOUND",      "queried; no row in catalogdb.sdss_id_to_catalog"),
    (2, "NO_GAIA_DR3",    "no Gaia DR3 source matched"),
    (3, "NO_2MASS",       "no 2MASS PSC source matched"),
    (4, "NO_CARTON",      "no SDSS-V carton assignment found (sdss_id_flat -> target -> carton_to_target)"),
    (5, "UNKNOWN_CARTON", "at least one carton_pk not encoded in sdss5_target_flags (neither the pk nor its carton name has a bit in the sdss_semaphore mapping); see carton_pks/bit_source == 0"),
    (6, "CARTON_BY_NAME", "at least one carton_pk encoded via the bit of a same-name carton_pk from another targeting plan (name-level decoding is right, pk-level decoding is not); see carton_pks/bit_source == 2"),
]
bitval(name) = UInt8(1) << first(b for (b, n, _) in STATUS_BITS if n == name)

# The 0.4.0 output columns, in the order of the 0.4.0 product registry (allTargets_pass1_dr21.toml
# `column_order`). ORDER ONLY: which /meta columns are copied is decided by the input file.
const ORDER_040 = [
    "sdss_id", "catalogid", "version_id", "lead",
    "gaia_dr3_source", "twomass_psc", "twomass_designation", "tic_v8", "allstar_dr17_synspec_rev1",
    "gaia_ra", "gaia_ra_error", "gaia_dec", "gaia_dec_error",
    "gaia_parallax", "gaia_parallax_error", "gaia_pmra", "gaia_pmra_error", "gaia_pmdec", "gaia_pmdec_error",
    "gaia_ruwe", "gaia_duplicated_source",
    "gaia_phot_g_mean_mag", "gaia_phot_bp_mean_mag", "gaia_phot_rp_mean_mag", "gaia_phot_bp_rp_excess_factor",
    "gaia_radial_velocity", "gaia_radial_velocity_error",
    "twomass_j_m", "twomass_j_msigcom", "twomass_h_m", "twomass_h_msigcom", "twomass_k_m", "twomass_k_msigcom",
    "twomass_ph_qual", "twomass_rd_flg", "twomass_bl_flg", "twomass_cc_flg",
    "sdss5_target_flags", "n_cartons", "too", "too_id", "too_program",
    "n_exposures", "metadata_status",
]
# Output columns computed here (not copied from /meta), and output groups.
const DERIVED = ["sdss_id", "sdss5_target_flags", "n_cartons", "n_exposures", "metadata_status"]
const GROUPS = ["carton_pks", "carton_table"]
# /meta columns not passed through under their own name, and why.
const REPLACED_META = Dict("sdss5_target_flags" =>
    "rebuilt from the carton pairs with the --semaphore-mapping (see flags_rule); the /meta column was built by almanac with its own sdss-semaphore default mapping")
const SOLAR_SYSTEM_LEAD = "solar_system"

function parse_cli()
    s = ArgParseSettings(description = "Build allTargets from almanac /meta, passing through every /meta column")
    @add_arg_table! s begin
        "--release"; required = true; help = "release label written into descriptions/attributes, e.g. 0.4.1"
        "--allexposure"; required = true
        "--allexposure-label"; default = ""; help = "name of the allExposure table in descriptions (default: basename of --allexposure without .h5)"
        "--meta"; default = ""
        "--checkpoint"; default = ""
        "--out"; required = true
        "--title"; default = ""; help = "file attribute `title` (default: '<out basename>: target metadata for the <release> APOGEE reduction')"
        "--carton-table"; required = true; help = "CSV dump of targetdb.carton (carton_pk,carton,program,mapper,plan,tag) to embed"
        "--semaphore-mapping"; required = true; help = "sdss_semaphore carton-to-bit CSV (sdss5_target_<V>_with_groups.csv)"
        "--semaphore-version"; required = true; help = "sdss-semaphore release the mapping CSV comes from (recorded only)"
        "--name-fallback"; arg_type = Bool; default = true; help = "map pks without a bit via a same-name pk that has one"
        "--include-solar-system"; arg_type = Bool; default = true; help = "also keep sdss_id <= 0 rows that /meta marks lead == 'solar_system' (almanac's fixed Vesta row) when present in allExposure"
        "--keep-meta-target-flags"; arg_type = Bool; default = false; help = "also copy /meta sdss5_target_flags (almanac's own mapping) as meta_sdss5_target_flags"
    end
    a = parse_args(s)
    xor(isempty(a["meta"]), isempty(a["checkpoint"])) || error("give exactly one of --meta / --checkpoint")
    return a
end

"Flatten the targeting-flags array to K x N (Julia order; h5py sees N x K)."
function flags_matrix(x)
    ndims(x) == 2 && return x
    ndims(x) == 3 && return reshape(x, size(x, 1), size(x, 3))   # (K,1,N) -> (K,N)
    error("unexpected sdss5_target_flags shape $(size(x))")
end

"Describe an on-disk HDF5 datatype (recorded for new columns)."
function dtype_string(d)
    dt = datatype(d)
    cls = HDF5.API.h5t_get_class(dt)
    if cls == HDF5.API.H5T_STRING
        HDF5.API.h5t_is_variable_str(dt) && return "variable-length string"
        return "fixed-length string, $(HDF5.API.h5t_get_size(dt)) bytes; written here as variable-length UTF-8"
    elseif cls == HDF5.API.H5T_ENUM
        return "enum (h5py bool), $(sizeof(dt)) byte; written here as Int8 0/1"
    end
    return string(eltype(d))
end

dataset_attrs(d) = Dict{String,Any}(k => read(attributes(d)[k]) for k in keys(attributes(d)))

"""
Source = a finished /meta file. Columns are read one at a time (low memory). `names` is the /meta
column order: the group's `columns` attribute when present (checked against the datasets), else
the HDF5 iteration order of the group.
"""
function open_meta_file(path)
    f = h5open(path, "r")
    g = f["meta"]
    dsnames = [k for k in keys(g) if g[k] isa HDF5.Dataset]
    attrs = Dict{String,Any}(k => read(attributes(g)[k]) for k in keys(attributes(g)))
    names = if haskey(attrs, "columns")
        c = String.(collect(attrs["columns"]))
        Set(c) == Set(dsnames) && length(unique(c)) == length(c) ||
            error("/meta `columns` attribute disagrees with the datasets present: missing $(setdiff(c, dsnames)), extra $(setdiff(dsnames, c))")
        c
    else
        dsnames
    end
    getcol(c) = read(g[c])
    getdesc(c) = haskey(attributes(g[c]), "description") ? read(attributes(g[c])["description"]) : ""
    getattrs(c) = dataset_attrs(g[c])
    getdtype(c) = dtype_string(g[c])
    notfound = haskey(f, "meta_not_found/sdss_id") ? read(f["meta_not_found/sdss_id"]) : Int64[]
    pid = haskey(f, "meta_carton_pks/sdss_id") ? read(f["meta_carton_pks/sdss_id"]) : Int64[]
    ppk = haskey(f, "meta_carton_pks/carton_pk") ? read(f["meta_carton_pks/carton_pk"]) : Int64[]
    has_pairs = haskey(f, "meta_carton_pks/sdss_id")
    unknown = haskey(attrs, "unknown_carton_pks") ? Int64.(collect(attrs["unknown_carton_pks"])) : Int64[]
    return (; file = f, names, getcol, getdesc, getattrs, getdtype, notfound, pid, ppk, has_pairs,
            unknown, attrs, queried = nothing)
end

"""
Source = a running build's checkpoint (batch_*.h5, group `columns`): a PARTIAL preview. Every
column found in the batch files is carried. NOT exercised for 0.4.1 (format of newer almanac
checkpoints not checked).
"""
function read_checkpoint(dir)
    files = sort(filter(f -> occursin(r"^batch_\d+\.h5$", f), readdir(dir)))
    parts = Dict{String,Vector{Any}}(); order = String[]
    desc = Dict{String,String}(); dts = Dict{String,String}(); cattrs = Dict{String,Dict{String,Any}}()
    queried = Int64[]; pid = Int64[]; ppk = Int64[]; unknown = Set{Int64}()
    for fn in files
        h5open(joinpath(dir, fn), "r") do f
            read(attributes(f)["complete"]) == 1 || return
            append!(queried, read(f["queried_sdss_id"]))
            append!(pid, read(f["carton_sdss_id"])); append!(ppk, read(f["carton_pk"]))
            union!(unknown, Int64.(collect(read(attributes(f)["unknown_carton_pks"]))))
            g = f["columns"]
            haskey(g, "sdss_id") || return
            for c in keys(g)
                g[c] isa HDF5.Dataset || continue
                c in order || push!(order, c)
                v = c == "sdss5_target_flags" ? flags_matrix(read(g[c])) : read(g[c])
                push!(get!(parts, c, Any[]), v)
                desc[c] = haskey(attributes(g[c]), "description") ? read(attributes(g[c])["description"]) : ""
                dts[c] = dtype_string(g[c]); cattrs[c] = dataset_attrs(g[c])
            end
        end
    end
    cols = Dict{String,Any}()
    for (c, vs) in parts
        if c == "sdss5_target_flags"
            K = maximum(size(v, 1) for v in vs)
            cols[c] = reduce(hcat, [size(v, 1) == K ? v : vcat(v, zeros(eltype(v), K - size(v, 1), size(v, 2))) for v in vs])
        else
            cols[c] = ndims(vs[1]) == 1 ? reduce(vcat, vs) : cat(vs...; dims = ndims(vs[1]))
        end
    end
    found = get(cols, "sdss_id", Int64[])
    notfound = setdiff(queried, found)
    manifest = JSON.parsefile(joinpath(dir, "manifest.json"))
    attrs = Dict{String,Any}("checkpoint_dir" => abspath(dir), "n_batch_files" => length(files),
        "almanac_input_path" => get(manifest, "input_path", ""),
        "query_start_utc" => something(get(manifest, "query_start_utc", ""), ""))
    return (; file = nothing, names = order, getcol = c -> cols[c], getdesc = c -> get(desc, c, ""),
            getattrs = c -> cattrs[c], getdtype = c -> dts[c], notfound, pid, ppk, has_pairs = true,
            unknown = collect(unknown), attrs, queried = sort(unique(queried)))
end

fillvalue(::Type{T}) where {T<:AbstractFloat} = T(NaN)
fillvalue(::Type{T}) where {T<:Signed} = T(-1)
fillvalue(::Type{T}) where {T<:Unsigned} = zero(T)
fillvalue(::Type{Bool}) = false
fillvalue(::Type{String}) = ""
fillvalue(::Type{T}) where {T} = error("no fill value for element type $T (add a fillvalue method)")

"Re-index a column (row axis = last) to the output rows; rowidx[i] == 0 -> fill value."
function reindex(v, rowidx)
    N = length(rowidx); T = eltype(v)
    if ndims(v) == 1
        w = fill(fillvalue(T), N)
        for i in 1:N; rowidx[i] > 0 && (w[i] = v[rowidx[i]]); end
        return T == Bool ? UInt8.(w) : w
    end
    w = fill(fillvalue(T), size(v)[1:end-1]..., N)
    V = reshape(v, :, size(v, ndims(v))); W = reshape(w, :, N)
    for i in 1:N; rowidx[i] > 0 && (W[:, i] .= view(V, :, rowidx[i])); end
    return T == Bool ? UInt8.(w) : w
end

"Split one CSV line (double-quoted fields allowed)."
function parse_csv_line(l)
    out = String[]; buf = IOBuffer(); q = false
    for c in l
        if c == '"'; q = !q
        elseif c == ',' && !q; push!(out, String(take!(buf)))
        else; write(buf, c); end
    end
    push!(out, String(take!(buf)))
    return out
end

"Write `v`; arrays with >= 65,536 rows are chunked (65,536 rows) with shuffle + gzip level 4."
function putdata(parent, name, v)
    haskey(parent, name) && error("output name collision: $name written twice")
    n = size(v, ndims(v))
    if n >= 65_536
        ch = ndims(v) == 1 ? (65_536,) : (size(v)[1:end-1]..., 65_536)
        d = create_dataset(parent, name, datatype(v), dataspace(v); chunk = ch, shuffle = true, deflate = 4)
        write(d, v)
    else
        write(parent, name, v)
    end
    return parent[name]
end

read_csv_table(path) = (lines = filter(!isempty, readlines(path));
    hdr = String.(parse_csv_line(lines[1])); rows = [parse_csv_line(l) for l in lines[2:end]];
    all(length(r) == length(hdr) for r in rows) || error("ragged CSV $path");
    Dict(h => [r[j] for r in rows] for (j, h) in enumerate(hdr)))

"""
Resolve every targetdb carton_pk to (bit, source): source 1 = exact pk in the mapping,
2 = via a same-name pk that has a bit, 0 = none (bit -1). Returns Dicts and a notes vector.
"""
function resolve_bits(mapping_csv, carton, name_fallback)
    m = read_csv_table(mapping_csv)
    mpk = parse.(Int64, m["carton_pk"]); mbit = parse.(Int, m["bit"])
    length(unique(mpk)) == length(mpk) || error("duplicate carton_pk in mapping")
    vers = unique(m["SDSSC2BV"]); length(vers) == 1 || error("mapping has several SDSSC2BV: $vers")
    pk2bit = Dict(p => b for (p, b) in zip(mpk, mbit) if p >= 0)
    cpk = parse.(Int64, carton["carton_pk"])
    name = Dict(zip(cpk, carton["carton"])); prog = Dict(zip(cpk, zip(carton["program"], carton["mapper"])))
    for (p, n) in zip(mpk, m["name"])
        p >= 0 && haskey(name, p) && name[p] != n && error("pk $p: mapping name $n != targetdb $(name[p])")
    end
    siblings = Dict{String,Vector{Int64}}()
    for p in cpk; haskey(pk2bit, p) && push!(get!(siblings, name[p], Int64[]), p); end
    bit = Dict{Int64,Int}(); src = Dict{Int64,UInt8}(); via = Dict{Int64,Int64}()
    for p in cpk
        if haskey(pk2bit, p)
            bit[p] = pk2bit[p]; src[p] = 0x01
        elseif name_fallback && haskey(siblings, name[p])
            sib = siblings[name[p]]
            same = filter(q -> prog[q] == prog[p], sib)
            q = maximum(isempty(same) ? sib : same)
            bit[p] = pk2bit[q]; src[p] = 0x02; via[p] = q
        else
            bit[p] = -1; src[p] = 0x00
        end
    end
    return (; bit, src, via, version = only(vers), maxbit = maximum(mbit), nrows = length(mpk))
end

"Units = a trailing \"[...]\" in the description (0.4.0 rule), except an abundance label such as [Fe/H] or [M/H]."
function units_from(desc)
    m = match(r"\[([^\]]+)\]\s*$", desc)
    m === nothing && return ""
    occursin(r"^[A-Z][a-z]?/[A-Z][a-z]?$", m.captures[1]) && return ""   # e.g. zgr_fe_h_confidence "... [Fe/H]"
    return m.captures[1]
end

function main()
    a = parse_cli()
    t0 = time()
    rel = a["release"]
    exlabel = isempty(a["allexposure-label"]) ? replace(basename(a["allexposure"]), r"\.h5$" => "") : a["allexposure-label"]

    # 1. source (opened first: the solar-system ids are taken from it)
    partial = !isempty(a["checkpoint"])
    src = partial ? read_checkpoint(a["checkpoint"]) : open_meta_file(a["meta"])
    "sdss_id" in src.names || error("source has no sdss_id column")
    msid0 = src.getcol("sdss_id")

    # 2. column plan (fail loudly on any name clash)
    passthrough = [c for c in src.names if c != "sdss_id" && !haskey(REPLACED_META, c)]
    produced = vcat(DERIVED, GROUPS, a["keep-meta-target-flags"] ? ["meta_sdss5_target_flags"] : String[])
    clash = intersect(passthrough, produced)
    isempty(clash) || error("name collision: /meta column(s) $(clash) would collide with columns/groups computed by this builder")
    from_meta_040 = [c for c in ORDER_040 if !(c in DERIVED)]
    missing040 = setdiff(from_meta_040, passthrough)
    isempty(missing040) || error("0.4.0 column(s) missing from the /meta input: $(missing040)")
    extra = [c for c in passthrough if !(c in ORDER_040)]                    # /meta order
    column_order = vcat(ORDER_040, extra, a["keep-meta-target-flags"] ? ["meta_sdss5_target_flags"] : String[])
    length(unique(column_order)) == length(column_order) || error("duplicate output column names")
    println("source columns: $(length(src.names)); passed through: $(length(passthrough)) ",
            "($(length(from_meta_040)) carried by 0.4.0, $(length(extra)) added); replaced: ",
            intersect(collect(keys(REPLACED_META)), src.names))

    # 3. solar-system ids (fixed almanac row, lead == "solar_system"; sdss_id <= 0)
    lead_all = "lead" in src.names ? src.getcol("lead") : fill("", length(msid0))
    ss_ids = a["include-solar-system"] ? sort(unique(msid0[(msid0 .<= 0) .& (strip.(lead_all) .== SOLAR_SYSTEM_LEAD)])) : Int64[]
    lead_all = nothing

    # 4. allExposure ids and per-id exposure counts
    exp_ids = h5open(f -> read(f["sdss_id"]), a["allexposure"], "r")
    pos = sort!(filter(s -> s > 0 || insorted(s, ss_ids), exp_ids))
    ids = unique(pos)
    counts = zeros(Int32, length(ids))
    let j = 1
        for x in pos
            while ids[j] != x; j += 1; end
            counts[j] += 1
        end
    end
    N = length(ids)
    ss_in = filter(s -> s <= 0, ids)
    println("allExposure: $(length(exp_ids)) rows, $(N) unique sdss_id kept ($(length(ss_in)) solar-system: $(ss_in))")
    exp_ids = nothing; pos = nothing

    # 5. row map
    perm = sortperm(msid0); msid = msid0[perm]
    length(unique(msid)) == length(msid) || error("duplicate sdss_id in source")
    rowidx = zeros(Int, N)
    for (i, s) in enumerate(ids)
        r = searchsortedfirst(msid, s)
        if r <= length(msid) && msid[r] == s
            rowidx[i] = perm[r]
        end
    end
    n_matched = count(>(0), rowidx)
    println("source rows: $(length(msid)); matched to allExposure ids: $(n_matched)")
    msid = nothing; msid0 = nothing; perm = nothing

    # 6. carton pairs restricted to the table's ids; sdss5_target_flags rebuilt from them
    keep = [insorted(s, ids) for s in src.pid]
    pid, ppk = src.pid[keep], src.ppk[keep]
    keep = nothing
    o = sortperm(collect(zip(pid, ppk))); pid, ppk = pid[o], ppk[o]; o = nothing
    carton = read_csv_table(a["carton-table"])
    res = resolve_bits(a["semaphore-mapping"], carton, a["name-fallback"])
    K = cld(res.maxbit + 1, 8)
    flags = zeros(UInt8, K, N)
    pbit = Vector{Int16}(undef, length(ppk)); psrc = Vector{UInt8}(undef, length(ppk))
    n_cartons = zeros(Int16, N); has_unknown = falses(N); has_name = falses(N)
    for (k, (s, p)) in enumerate(zip(pid, ppk))
        i = searchsortedfirst(ids, s)
        n_cartons[i] += 1
        haskey(res.src, p) || error("carton_pk $p not in --carton-table")
        b = res.bit[p]; pbit[k] = b; psrc[k] = res.src[p]
        if b >= 0
            flags[(b >> 3) + 1, i] |= UInt8(1) << (b & 7)
            psrc[k] == 0x02 && (has_name[i] = true)
        else
            has_unknown[i] = true
        end
    end
    pair_counts = Dict(n => count(==(v), psrc) for (n, v) in (("none", 0x00), ("exact", 0x01), ("name", 0x02)))
    pk_counts = Dict(n => length(unique(ppk[psrc .== v])) for (n, v) in (("none", 0x00), ("exact", 0x01), ("name", 0x02)))
    println("pairs by bit source: ", pair_counts, "; distinct pks: ", pk_counts)

    # 7. status
    notfound = Set(src.notfound)
    queried = src.queried === nothing ? nothing : Set(src.queried)
    status = zeros(UInt8, N)
    gaia = reindex(src.getcol("gaia_dr3_source"), rowidx)
    tm = reindex(src.getcol("twomass_designation"), rowidx)
    for i in 1:N
        s = ids[i]
        if rowidx[i] == 0
            if queried !== nothing && !(s in queried)
                status[i] |= bitval("NOT_QUERIED")
            else
                status[i] |= bitval("NOT_FOUND")
            end
            continue
        end
        gaia[i] <= 0 && (status[i] |= bitval("NO_GAIA_DR3"))
        isempty(strip(tm[i])) && (status[i] |= bitval("NO_2MASS"))
        n_cartons[i] == 0 && src.has_pairs && (status[i] |= bitval("NO_CARTON"))
        has_unknown[i] && (status[i] |= bitval("UNKNOWN_CARTON"))
        has_name[i] && (status[i] |= bitval("CARTON_BY_NAME"))
    end
    gaia = nothing; tm = nothing

    # 8. write (temp file, then rename); creation order tracked = column_order
    meta_flags_nbytes = "sdss5_target_flags" in src.names ? size(flags_matrix(src.getcol("sdss5_target_flags")), 1) : 0
    shapes = Dict{String,String}()
    tmp = a["out"] * ".tmp-$(getpid())"
    h5open(tmp, "w"; fcpl = HDF5.FileCreateProperties(track_order = true)) do f
        function put(name, v, d)
            putdata(f, name, v)
            attributes(f[name])["description"] = d
            u = units_from(d); isempty(u) || (attributes(f[name])["units"] = u)
        end
        for c in column_order
            if c == "sdss_id"
                d = isempty(ss_in) ? "SDSS identifier (one row per sdss_id in the $rel allExposure table)" :
                    "SDSS identifier (one row per sdss_id in the $rel allExposure table; sdss_id $(join(ss_in, ", ")) is a solar-system body, lead '$(SOLAR_SYSTEM_LEAD)')"
                put("sdss_id", ids, d)
            elseif c == "sdss5_target_flags"
                put(c, flags, "SDSS-V targeting bits, sdss_semaphore TargetingFlags with carton-to-bit mapping version $(res.version) (sdss-semaphore $(a["semaphore-version"])); $K bytes per sdss_id. Bits set by exact carton_pk, or by carton name where the pk has no bit (metadata_status bit 6); see file attrs flags_rule")
            elseif c == "n_cartons"
                put(c, n_cartons, "number of SDSS-V carton assignments (rows in carton_pks/ for this sdss_id)")
            elseif c == "n_exposures"
                put(c, counts, "number of rows (exposure x fiber) for this sdss_id in $exlabel")
            elseif c == "metadata_status"
                put(c, status, "bitmask, 0 = ok; bits: " * join(["$b $n ($d)" for (b, n, d) in STATUS_BITS], "; "))
            elseif c in from_meta_040
                # 0.4.0 column: exactly the 0.4.0 treatment (description, 0.4.0 overrides, units)
                d = src.getdesc(c)
                c == "gaia_duplicated_source" && (d *= " (0/1)")
                c == "too" && (d *= " (0/1)")
                c == "catalogid" && (d = "catalogid of one sdss_id_to_catalog row for this sdss_id (version not chosen; see version_id)")
                put(c, reindex(src.getcol(c), rowidx), d)
            else
                # pass-through column added since 0.4.0 (or meta_sdss5_target_flags)
                mc = c == "meta_sdss5_target_flags" ? "sdss5_target_flags" : c
                v = src.getcol(mc)
                mc == "sdss5_target_flags" && (v = flags_matrix(v))
                w = reindex(v, rowidx); v = nothing
                ds = putdata(f, c, w)
                at = src.getattrs(mc)
                d = get(at, "description", "")
                c == "meta_sdss5_target_flags" && (d = "sdss5_target_flags exactly as in /meta (almanac's sdss-semaphore $(get(src.attrs, "sdss_semaphore_version", "?")) default mapping, exact carton_pk only); $(size(w, 1)) bytes per sdss_id. Superseded by sdss5_target_flags")
                attributes(ds)["description"] = d
                for (k, val) in at
                    k == "description" && continue
                    attributes(ds)[k] = val
                end
                u = units_from(d); (isempty(u) || haskey(at, "units")) || (attributes(ds)["units"] = u)
                attributes(ds)["source"] = "/meta/" * mc
                attributes(ds)["meta_dtype"] = src.getdtype(mc)
                if ndims(w) > 1
                    shapes[c] = join(size(w)[1:end-1], "x")
                    attributes(ds)["shape_per_row"] = collect(Int64, size(w)[1:end-1])
                end
                w = nothing
            end
        end
        g = create_group(f, "carton_pks")
        putdata(g, "sdss_id", pid); putdata(g, "carton_pk", ppk)
        putdata(g, "bit", pbit); putdata(g, "bit_source", psrc)
        attributes(g)["description"] = "every (sdss_id, carton_pk) pair returned by targetdb; join carton_pk to carton_table/carton_pk (targetdb.carton.pk) for names"
        attributes(g["bit"])["description"] = "sdss5_target_flags bit set for this pair (-1 = none)"
        attributes(g["bit_source"])["description"] = "0 = no bit (UNKNOWN_CARTON), 1 = exact carton_pk in the semaphore mapping, 2 = bit of a same-name carton_pk (CARTON_BY_NAME)"
        ct = create_group(f, "carton_table")
        cpk = parse.(Int64, carton["carton_pk"])
        for h in ("carton_pk", "carton", "program", "mapper", "plan", "tag")
            putdata(ct, h, h == "carton_pk" ? cpk : String.(carton[h]))
        end
        putdata(ct, "semaphore_bit", Int16[res.bit[p] for p in cpk])
        putdata(ct, "semaphore_bit_source", UInt8[res.src[p] for p in cpk])
        putdata(ct, "semaphore_bit_via_pk", Int64[res.src[p] == 0x02 ? res.via[p] : (res.src[p] == 0x01 ? p : -1) for p in cpk])
        attributes(ct)["description"] = "targetdb.carton joined to targetdb.mapper and targetdb.version (plan, tag), dumped once from targetdb; semaphore_bit / _source / _via_pk = resolved sdss5_target_flags bit, how (0 none, 1 exact, 2 name) and the carton_pk whose bit is used"
        attributes(ct)["source"] = abspath(a["carton-table"])

        fa = attributes(f)
        fa["title"] = isempty(a["title"]) ? "$(replace(basename(a["out"]), r"\.h5$" => "")): target metadata for the $rel APOGEE reduction" : a["title"]
        fa["row_definition"] = "one row per sdss_id > 0 in $exlabel.h5; join on sdss_id" *
            (isempty(ss_in) ? "" : "; plus sdss_id $(join(ss_in, ", ")) (solar-system body: /meta lead '$(SOLAR_SYSTEM_LEAD)', no catalog values)")
        fa["partial"] = UInt8(partial)
        fa["created_utc"] = string(now(UTC)) * "Z"
        fa["builder"] = abspath(PROGRAM_FILE)
        fa["builder_sha256"] = bytes2hex(open(sha256, abspath(PROGRAM_FILE)))
        fa["allexposure_path"] = abspath(a["allexposure"])
        fa["meta_source"] = partial ? abspath(a["checkpoint"]) : abspath(a["meta"])
        for (k, v) in src.attrs
            v isa Union{AbstractString,Number} || continue
            fa["meta_" * k] = v
        end
        fa["sdss_semaphore_version"] = a["semaphore-version"]
        fa["sdss_semaphore_mapping_version"] = parse(Int, res.version)
        fa["sdss_semaphore_mapping_file"] = abspath(a["semaphore-mapping"])
        fa["sdss_semaphore_mapping_sha256"] = bytes2hex(open(sha256, a["semaphore-mapping"]))
        fa["sdss5_target_flags_nbytes"] = K
        fa["flags_rule"] = "sdss5_target_flags rebuilt from carton_pks/ pairs (not copied from /meta). Per carton_pk: (1) exact: bit of the pk in sdss_semaphore mapping version $(res.version); " *
            (a["name-fallback"] ? "(2) name: else the bit of a pk with the same targetdb carton name that has one (prefer same program+mapper, then highest carton_pk); " : "(2) name fallback disabled; ") *
            "(3) none: else not encoded, metadata_status bit 5 UNKNOWN_CARTON. No new bits are assigned. Decode with sdss_semaphore TargetingFlags(array, sdssc2bv=$(res.version)) or SDSSC2BV=$(res.version); name-level decoding is exact, pk-level decoding is not for rows with metadata_status bit 6"
        fa["n_pairs"] = length(ppk)
        fa["n_pairs_bit_exact_pk"] = pair_counts["exact"]
        fa["n_pairs_bit_by_name"] = pair_counts["name"]
        fa["n_pairs_unmapped"] = pair_counts["none"]
        fa["n_carton_pks_exact"] = pk_counts["exact"]
        fa["n_carton_pks_by_name"] = pk_counts["name"]
        fa["n_carton_pks_unmapped"] = pk_counts["none"]
        fa["unmapped_carton_pks"] = sort(unique(ppk[psrc .== 0x00]))
        fa["n_rows"] = N
        fa["n_with_meta"] = n_matched
        for (b, n, _) in STATUS_BITS
            fa["n_status_$(n)"] = count(x -> x & (UInt8(1) << b) != 0, status)
        end
        fa["status_bits"] = JSON.json([Dict("bit" => b, "name" => n, "description" => d) for (b, n, d) in STATUS_BITS])

        # 0.4.1 additions: column provenance and input description
        function newattr(k, v)
            haskey(fa, k) && error("file attribute collision: $k")
            fa[k] = v
        end
        newattr("release", rel)
        newattr("column_order", column_order)
        newattr("column_order_rule", "the 0.4.0 columns in the 0.4.0 registry order, then every other /meta column in /meta order (pass-through, nothing hard-coded), then meta_sdss5_target_flags if requested; groups carton_pks/, carton_table/ last. The column set follows the /meta input: readers must not assume a fixed set of columns.")
        newattr("columns_from_040", ORDER_040)
        newattr("meta_passthrough_columns", [c for c in column_order if c in passthrough])
        newattr("n_meta_passthrough_columns", length(passthrough))
        newattr("columns_added_since_040", [c for c in column_order if !(c in ORDER_040)])
        newattr("meta_columns_not_passed_through", JSON.json(Dict(c => r for (c, r) in REPLACED_META if c in src.names)))
        newattr("meta_n_columns", length(src.names))
        newattr("meta_sdss5_target_flags_nbytes", meta_flags_nbytes)
        haskey(fa, "meta_extended") || newattr("meta_extended", 0)
        newattr("meta_extended_source", haskey(src.attrs, "extended") ? "/meta attribute `extended`" :
            "absent from /meta attributes (almanac without --extended support or default mode): recorded as 0")
        if haskey(src.attrs, "extended_column_groups")
            newattr("meta_extended_column_groups", String.(collect(src.attrs["extended_column_groups"])))
        end
        newattr("multidim_columns", JSON.json(merge(Dict("sdss5_target_flags" => string(K)), shapes)))
        newattr("solar_system_sdss_id", collect(Int64, ss_in))
        newattr("n_solar_system_rows", length(ss_in))
        newattr("solar_system_rule", "sdss_id <= 0 rows are kept only if /meta marks them lead == '$(SOLAR_SYSTEM_LEAD)' (almanac's fixed, never-queried row; all catalog values NaN / -1 / '') and they occur in allExposure; n_exposures counts every allExposure row with that sdss_id" *
            (a["include-solar-system"] ? "" : " (DISABLED by --include-solar-system false)"))
    end
    src.file === nothing || close(src.file)
    mv(tmp, a["out"]; force = true)
    println("wrote $(a["out"]): $N rows, $(n_matched) with metadata, $(length(column_order)) columns ",
            "($(length(extra)) added since 0.4.0), $(length(pid)) carton pairs, ",
            "status counts ", Dict(n => count(x -> x & (UInt8(1) << b) != 0, status) for (b, n, _) in STATUS_BITS),
            "; $(round(time() - t0, digits = 1)) s")
end

main()
