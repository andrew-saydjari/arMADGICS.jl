## Curated registry + release config loading, and per-column description resolution

"""
    load_release(path) -> Dict

Read a release config (`releases/<version>.toml`). Relative `registry` entries are
resolved against `<datamodel>/registry/`.
"""
function load_release(path::AbstractString)
    cfg = TOML.parsefile(path)
    haskey(cfg, "release") || error("release config $path has no `release` key")
    haskey(cfg, "products") || error("release config $path has no [[products]]")
    return cfg
end

"""
    load_registry(path) -> Dict

A registry file has a `[product]` table (page-level text) and a `[columns.<name>]` table
per column. Missing file -> empty registry (every column then falls back to file
attributes and harvested docs, and is flagged UNDOCUMENTED if those are empty too).
"""
function load_registry(path::AbstractString)
    isfile(path) || return Dict{String,Any}("product" => Dict{String,Any}(),
                                             "columns" => Dict{String,Any}(), "_missing" => true)
    reg = TOML.parsefile(path)
    get!(reg, "product", Dict{String,Any}())
    get!(reg, "columns", Dict{String,Any}())
    return reg
end

nonempty(x) = x isa AbstractString ? !isempty(strip(x)) : x !== nothing

"""
    resolve_column(name, reg, col, harvested) -> NamedTuple

Description precedence: curated registry -> dataset `description` attribute -> text
harvested from the release Markdown docs. `source` records which one won (or
`:undocumented`). Bit tables come from the registry; enumerations from the registry or,
failing that, from a `categories` dataset attribute (code = position - 1).
"""
function resolve_column(name::AbstractString, reg::AbstractDict, col, harvested::AbstractDict)
    rc = get(get(reg, "columns", Dict()), name, Dict{String,Any}())
    attrs_ = col === nothing ? Dict{String,Any}() : col.attrs
    desc, source = "", :undocumented
    if nonempty(get(rc, "description", ""))
        desc, source = rc["description"], :registry
    elseif nonempty(get(attrs_, "description", ""))
        desc, source = string(attrs_["description"]), :attribute
    elseif nonempty(get(harvested, name, ""))
        desc, source = harvested[name], :harvested
    end
    values = get(rc, "values", Any[])
    if isempty(values) && haskey(attrs_, "categories")
        cats = attrs_["categories"]
        values = [Dict{String,Any}("value" => i - 1, "meaning" => isempty(c) ? "(empty string: unset)" : c)
                  for (i, c) in enumerate(cats)]
    end
    return (description = desc, source = source,
            units = get(rc, "units", ""),
            null = get(rc, "null", ""),
            bits = get(rc, "bits", Any[]),
            bits_note = get(rc, "bits_note", ""),
            values = values,
            attr_bits = get(attrs_, "bits", ""))
end
