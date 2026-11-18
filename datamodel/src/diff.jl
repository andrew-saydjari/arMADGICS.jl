## Snapshots and diffs

"""
    snapshot_dict(products, release) -> Dict

The per-release schema snapshot: product -> {kind, n_rows, columns -> {dtype, shape}}.
Row counts are recorded for information but are not part of the schema diff (they are
expected to change between releases).
"""
function snapshot_dict(products::AbstractVector{ProductSchema}, release::AbstractString)
    prods = Dict{String,Any}()
    for p in products
        cols = Dict{String,Any}()
        for c in p.columns
            cols[c.name] = Dict{String,Any}("dtype" => c.dtype, "shape" => c.shape)
        end
        prods[p.name] = Dict{String,Any}("kind" => p.kind, "n_rows" => p.n_rows, "columns" => cols)
    end
    return Dict{String,Any}("release" => String(release),
                            "generated" => string(now(UTC)) * "Z",
                            "products" => prods)
end

function write_snapshot(path::AbstractString, snap::AbstractDict)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Schema snapshot written by datamodel/generate.jl. Do not edit by hand;")
        println(io, "# it is the reference the next release is diffed against.")
        TOML.print(io, snap; sorted = true)
    end
    return path
end

read_snapshot(path::AbstractString) = TOML.parsefile(path)

"""
    find_previous_snapshot(dir, release) -> Union{String,Nothing}

The snapshot in `dir` with the largest version number strictly below `release`.
Files that do not parse as `<version>.toml` are ignored.
"""
function find_previous_snapshot(dir::AbstractString, release::AbstractString)
    isdir(dir) || return nothing
    cur = VersionNumber(release)
    best = nothing
    bestv = v"0.0.0-"
    for fn in readdir(dir)
        endswith(fn, ".toml") || continue
        v = tryparse(VersionNumber, fn[1:end-5])
        v === nothing && continue
        if v < cur && v > bestv
            best, bestv = joinpath(dir, fn), v
        end
    end
    return best
end

"""
    diff_snapshots(prev, cur) -> Vector{Change}

Schema changes from `prev` to `cur`: removed/new products, removed/new columns, and
dtype or shape changes of columns present in both. Reference products are not in
snapshots and are never diffed.
"""
function diff_snapshots(prev::AbstractDict, cur::AbstractDict)
    out = Change[]
    pp = get(prev, "products", Dict())
    cp = get(cur, "products", Dict())
    for pname in sort(collect(keys(pp)))
        if !haskey(cp, pname)
            push!(out, Change(:removed_product, pname, ""))
            for col in sort(collect(keys(pp[pname]["columns"])))
                push!(out, Change(:removed_column, pname, col; old = pp[pname]["columns"][col]["dtype"]))
            end
        end
    end
    for pname in sort(collect(keys(cp)))
        ccols = cp[pname]["columns"]
        if !haskey(pp, pname)
            push!(out, Change(:new_product, pname, ""))
            for col in sort(collect(keys(ccols)))
                push!(out, Change(:new_column, pname, col; new = ccols[col]["dtype"]))
            end
            continue
        end
        pcols = pp[pname]["columns"]
        for col in sort(collect(keys(pcols)))
            haskey(ccols, col) || push!(out, Change(:removed_column, pname, col; old = pcols[col]["dtype"]))
        end
        for col in sort(collect(keys(ccols)))
            if !haskey(pcols, col)
                push!(out, Change(:new_column, pname, col; new = ccols[col]["dtype"]))
                continue
            end
            a, b = pcols[col], ccols[col]
            a["dtype"] != b["dtype"] &&
                push!(out, Change(:dtype_changed, pname, col; old = a["dtype"], new = b["dtype"]))
            a["shape"] != b["shape"] &&
                push!(out, Change(:shape_changed, pname, col; old = a["shape"], new = b["shape"]))
        end
    end
    return out
end

"""
    registry_check(product, reg, harvested; registry_file) -> Vector{Change}

Columns in the file with no description from any source (`:undocumented`), and registry
entries naming a column that is not in the file (`:registry_orphan`).
"""
function registry_check(p::ProductSchema, reg::AbstractDict, harvested::AbstractDict;
                        registry_file::AbstractString = "")
    out = Change[]
    for c in p.columns
        r = resolve_column(c.name, reg, c, harvested)
        r.source === :undocumented &&
            push!(out, Change(:undocumented, p.name, c.name; new = c.dtype, note = registry_file))
    end
    present = Set(column_names(p))
    for name in sort(collect(keys(get(reg, "columns", Dict()))))
        name in present || push!(out, Change(:registry_orphan, p.name, name; note = registry_file))
    end
    return out
end

## Loud reporting

const BAR = "!"^96

function banner_line(c::Change)
    pc = isempty(c.column) ? c.product : "$(c.product).$(c.column)"
    if c.kind === :removed_column
        return "COLUMN REMOVED: $pc — are we sure we don't want this column?"
    elseif c.kind === :new_column
        reg = isempty(c.note) ? "the registry" : c.note
        c.old == "documented" &&
            return "NEW COLUMN: $pc ($(c.new)) — already described in $reg; confirm the entry is current"
        return "NEW COLUMN: $pc ($(c.new)) — needs a description in $reg"
    elseif c.kind === :dtype_changed
        return "DTYPE CHANGED: $pc  $(c.old) -> $(c.new) — intended? downstream readers may break"
    elseif c.kind === :shape_changed
        return "SHAPE CHANGED: $pc  $(c.old) -> $(c.new) — intended? downstream readers may break"
    elseif c.kind === :removed_product
        return "PRODUCT REMOVED: $(c.product) — are we sure this file should not ship?"
    elseif c.kind === :new_product
        return "NEW PRODUCT: $(c.product) — needs a registry file and a page description"
    elseif c.kind === :undocumented
        return "UNDOCUMENTED COLUMN: $pc ($(c.new)) — add a description to $(isempty(c.note) ? "the registry" : c.note)"
    elseif c.kind === :registry_orphan
        return "REGISTRY ENTRY WITHOUT COLUMN: $pc is documented in $(c.note) but is not in the files — stale entry or a dropped column?"
    elseif c.kind === :initial
        return "INITIAL: $(c.note)"
    end
    return "$(c.kind): $pc"
end

"""
    report_changes(io, changes; title)

Print every change as its own banner block, so none can scroll past unnoticed.
"""
function report_changes(io::IO, changes::AbstractVector{Change}; title::AbstractString = "")
    isempty(changes) && return
    println(io)
    println(io, BAR)
    isempty(title) || println(io, "!!!  ", title)
    println(io, BAR)
    for c in changes
        println(io, "!!!  ", banner_line(c))
    end
    println(io, BAR)
    println(io)
end
