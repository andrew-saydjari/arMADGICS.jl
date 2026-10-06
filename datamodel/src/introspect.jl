## Introspection: read the files and build ProductSchema values. Read-only.

const JULIA_TO_DTYPE = Dict(
    "Int8" => "int8", "Int16" => "int16", "Int32" => "int32", "Int64" => "int64",
    "UInt8" => "uint8", "UInt16" => "uint16", "UInt32" => "uint32", "UInt64" => "uint64",
    "Float16" => "float16", "Float32" => "float32", "Float64" => "float64", "Bool" => "bool",
    "String" => "string",
)

"""
    normalize_dtype(T) -> String

Map a Julia element type (from HDF5.jl or CSV.jl) to a short, numpy-like name.
Fixed-length HDF5 strings become `S<n>`; CSV.jl's inline string widths (String7,
String15, ...) all become `string` so that a change in the longest value does not show
up as a schema change. `Union{Missing,T}` becomes `<T>?`.
"""
function normalize_dtype(@nospecialize(T))
    if T isa Union
        nonmiss = filter(t -> t !== Missing, Base.uniontypes(T))
        if length(nonmiss) == 1 && length(Base.uniontypes(T)) == 2
            return normalize_dtype(nonmiss[1]) * "?"
        end
        return join(sort(normalize_dtype.(nonmiss)), "|")
    end
    T === Missing && return "missing"
    s = string(T)
    haskey(JULIA_TO_DTYPE, s) && return JULIA_TO_DTYPE[s]
    m = match(r"FixedString\{(\d+)", s)
    m !== nothing && return "S" * m.captures[1]
    occursin(r"^(InlineStrings\.)?String\d+$", s) && return "string"
    occursin("VariableArray", s) && return "vlen"
    occursin("String", s) && return "string"
    return lowercase(s)
end

"""
    shape_string(dims, n_rows) -> String

Julia-order shape with the row axis (the last dimension, when it equals `n_rows`)
written as `N`. Scalars are `()`.
"""
function shape_string(dims::Tuple, n_rows::Integer)
    isempty(dims) && return "()"
    parts = [string(d) for d in dims]
    if n_rows > 0 && dims[end] == n_rows
        parts[end] = "N"
    end
    return length(parts) == 1 ? "(" * parts[1] * ",)" : "(" * join(parts, ", ") * ")"
end

"Reverse a Julia-order shape string to the h5py / numpy order."
function h5py_shape(shape::AbstractString)
    shape == "()" && return "()"
    inner = strip(shape, ['(', ')'])
    parts = filter(!isempty, strip.(split(inner, ",")))
    length(parts) == 1 && return "(" * parts[1] * ",)"
    return "(" * join(reverse(parts), ", ") * ")"
end

function plain_value(v)
    v isa AbstractString && return String(v)
    v isa AbstractArray && return [plain_value(x) for x in v]
    v isa Number && return v
    return string(v)
end

function read_attrs(obj)
    out = Dict{String,Any}()
    for k in keys(attrs(obj))
        out[k] = try
            plain_value(attrs(obj)[k])
        catch e
            "<unreadable attribute: $(sprint(showerror, e))>"
        end
    end
    return out
end

function dataset_schema(d::HDF5.Dataset, name, n_rows)
    T = eltype(d)
    dims = size(d)
    return ColumnSchema(String(name), normalize_dtype(T), shape_string(dims, n_rows),
                        read_attrs(d), Dict{String,Any}())
end

"Guess the row count of an HDF5 table: the `n_rows` root attribute, else the most common last dimension."
function infer_rows(f::HDF5.File)
    a = attrs(f)
    if haskey(a, "n_rows")
        return Int(a["n_rows"])
    end
    counts = Dict{Int,Int}()
    for k in keys(f)
        d = f[k]
        d isa HDF5.Dataset || continue
        s = size(d)
        isempty(s) && continue
        counts[s[end]] = get(counts, s[end], 0) + 1
    end
    isempty(counts) && return 0
    return first(sort(collect(counts), by = x -> (-x[2], -x[1])))[1]
end

"""
    introspect_hdf5_table(path; name) -> ProductSchema

One HDF5 file whose root datasets are the columns (allExposure, allVisitRV).
"""
function introspect_hdf5_table(path::AbstractString; name::AbstractString)
    h5open(path, "r") do f
        n = infer_rows(f)
        cols = ColumnSchema[]
        for k in keys(f)
            obj = f[k]
            if obj isa HDF5.Dataset
                push!(cols, dataset_schema(obj, k, n))
            else
                push!(cols, ColumnSchema(String(k), "group", "()", Dict{String,Any}(),
                                         Dict{String,Any}("note" => "HDF5 group")))
            end
        end
        return ProductSchema(String(name), "hdf5_table", n, filesize(path), cols, read_attrs(f),
                             Dict{String,Any}("files" => [basename(path)]))
    end
end

"""
    introspect_hdf5_perkey(source_dir; name, pattern, release_dir="") -> ProductSchema

A family of files `<prefix><key>.h5`, each holding one dataset named `<key>` and a
scalar `hdr` dataset whose attributes carry the pipeline provenance. Each key becomes a
column. When `release_dir` is given, each column records whether its file is present there.
"""
function introspect_hdf5_perkey(source_dir::AbstractString; name::AbstractString,
                                pattern::Regex, release_dir::AbstractString = "")
    files = sort(filter(fn -> occursin(pattern, fn), readdir(source_dir)))
    cols = ColumnSchema[]
    file_attrs = Dict{String,Any}()
    hdr_attrs = Dict{String,Any}()
    hdr_mismatch = String[]
    n_rows = 0
    total = 0
    for fn in files
        key = match(pattern, fn).captures[1]
        path = joinpath(source_dir, fn)
        h5open(path, "r") do f
            fa = read_attrs(f)
            if isempty(file_attrs)
                merge!(file_attrs, fa)
                n_rows = haskey(fa, "workup_nsamp") ? Int(fa["workup_nsamp"]) : infer_rows(f)
            end
            if haskey(f, "hdr")
                ha = read_attrs(f["hdr"])
                if isempty(hdr_attrs)
                    merge!(hdr_attrs, ha)
                elseif ha != hdr_attrs
                    push!(hdr_mismatch, key)
                end
            end
            if haskey(f, key)
                c = dataset_schema(f[key], key, n_rows)
            else
                c = ColumnSchema(String(key), "missing", "()", Dict{String,Any}(), Dict{String,Any}())
                c.extra["note"] = "file has no dataset named $(key)"
            end
            sz = filesize(path)
            total += sz
            c.extra["file"] = fn
            c.extra["size_bytes"] = sz
            if !isempty(release_dir)
                c.extra["in_release_dir"] = isfile(joinpath(release_dir, fn))
            end
            push!(cols, c)
        end
    end
    extra = Dict{String,Any}("hdr_attrs" => hdr_attrs, "hdr_mismatch" => hdr_mismatch,
                             "files" => files)
    return ProductSchema(String(name), "hdf5_perkey", n_rows, total, cols, file_attrs, extra)
end

"""
    introspect_csv(path; name) -> ProductSchema

A CSV file. Lines starting with `#` are kept as the file comment.
"""
function introspect_csv(path::AbstractString; name::AbstractString)
    comment = String[]
    for ln in eachline(path)
        startswith(ln, "#") || break
        push!(comment, strip(lstrip(ln, '#')))
    end
    f = CSV.File(path; comment = "#", ntasks = 1)
    names = String.(f.names)
    types = f.types
    cols = [ColumnSchema(names[i], normalize_dtype(types[i]), "(N,)") for i in eachindex(names)]
    fa = Dict{String,Any}()
    isempty(comment) || (fa["comment"] = join(comment, "\n"))
    return ProductSchema(String(name), "csv", length(f), filesize(path), cols, fa,
                         Dict{String,Any}("files" => [basename(path)]))
end

"""
    introspect_product(pcfg::Dict) -> Union{ProductSchema,Nothing}

Dispatch on the product `kind` from the release config. `reference` products have no
file to read and return `nothing`.
"""
function introspect_product(pcfg::AbstractDict)
    kind = pcfg["kind"]
    name = pcfg["name"]
    if kind == "hdf5_table"
        return introspect_hdf5_table(pcfg["source"]; name = name)
    elseif kind == "hdf5_perkey"
        return introspect_hdf5_perkey(pcfg["source_dir"]; name = name,
                                      pattern = Regex(pcfg["pattern"]),
                                      release_dir = get(pcfg, "release_dir", ""))
    elseif kind == "csv"
        return introspect_csv(pcfg["source"]; name = name)
    elseif kind == "reference"
        return nothing
    else
        error("unknown product kind '$kind' for product $name")
    end
end
