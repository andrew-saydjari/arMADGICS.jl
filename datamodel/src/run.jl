## One full generator run: introspect -> snapshot -> diff -> report -> render

const DATAMODEL_ROOT = normpath(joinpath(@__DIR__, ".."))

function git_describe(dir)
    try
        sha = readchomp(`git -C $dir rev-parse --short HEAD`)
        dirty = !isempty(readchomp(`git -C $dir status --porcelain -- . ":(exclude)snapshots"`))
        return dirty ? "$sha (uncommitted changes)" : sha
    catch
        return "unknown"
    end
end

function releases_with_product(snapshot_dir, current_release, current_names)
    out = Dict{String,Vector{String}}()
    cur = VersionNumber(current_release)
    if isdir(snapshot_dir)
        for fn in readdir(snapshot_dir)
            endswith(fn, ".toml") || continue
            v = tryparse(VersionNumber, fn[1:end-5])
            (v === nothing || v >= cur) && continue
            snap = try
                read_snapshot(joinpath(snapshot_dir, fn))
            catch
                continue
            end
            for p in keys(get(snap, "products", Dict()))
                push!(get!(out, p, String[]), string(v))
            end
        end
    end
    for n in current_names
        push!(get!(out, n, String[]), current_release)
    end
    for v in values(out)
        sort!(unique!(v), by = VersionNumber)
    end
    return out
end

"""
    run_release(config_path; outdir, snapshot_dir, previous, accept_changes, allow_undocumented,
                write_snapshot, io) -> (exitcode, changes)

Exit codes: 0 clean (or every finding acknowledged); 2 unacknowledged schema changes
against the previous snapshot; 3 undocumented columns or stale registry entries (and no
unacknowledged schema change). `--accept-changes` acknowledges the first,
`--allow-undocumented` the second. The snapshot for this release is written either way,
so a re-run after review sees the same diff.
"""
function run_release(config_path::AbstractString;
                     outdir::AbstractString,
                     snapshot_dir::AbstractString = joinpath(DATAMODEL_ROOT, "snapshots"),
                     registry_dir::AbstractString = joinpath(DATAMODEL_ROOT, "registry"),
                     previous::Union{Nothing,AbstractString} = nothing,
                     accept_changes::Bool = false,
                     allow_undocumented::Bool = false,
                     write_snapshot_file::Bool = true,
                     io::IO = stdout)
    cfg = load_release(config_path)
    release = string(cfg["release"])
    println(io, "datamodel: release $release  (config $(config_path))")

    pcfgs = cfg["products"]
    schemas = Dict{String,ProductSchema}()
    regs = Dict{String,Any}()
    harvested = Dict{String,Dict{String,String}}()
    for pcfg in pcfgs
        name = pcfg["name"]
        regpath = joinpath(registry_dir, get(pcfg, "registry", name * ".toml"))
        regs[name] = load_registry(regpath)
        harvested[name] = harvest_docs(get(pcfg, "harvest", String[]))
        t0 = time()
        p = introspect_product(pcfg)
        if p === nothing
            p = reference_schema(pcfg, regs[name])
        end
        schemas[name] = p
        @printf(io, "  %-28s %-12s %4d columns  %12s rows  (%.1fs)\n", name, p.kind, length(p.columns),
                p.n_rows > 0 ? commas(p.n_rows) : "-", time() - t0)
    end

    data_products = [schemas[p["name"]] for p in pcfgs if p["kind"] != "reference"]
    snap = snapshot_dict(data_products, release)
    snap_path = joinpath(snapshot_dir, release * ".toml")
    if write_snapshot_file
        write_snapshot(snap_path, snap)
        println(io, "  snapshot written: ", snap_path)
    end

    prev_path = previous === nothing ? find_previous_snapshot(snapshot_dir, release) : previous
    changes = Change[]
    previous_label = ""
    regfile(name) = "datamodel/registry/" * get(only(filter(p -> p["name"] == name, pcfgs)), "registry", name * ".toml")
    if prev_path === nothing
        println(io, "  no previous snapshot in $(snapshot_dir): recording the initial schema")
    else
        prev = read_snapshot(prev_path)
        previous_label = string(get(prev, "release", basename(prev_path)))
        println(io, "  comparing against previous snapshot: ", prev_path, " (release ", previous_label, ")")
        for c in diff_snapshots(prev, snap)
            if c.kind === :new_column && haskey(schemas, c.product)
                c = Change(c.kind, c.product, c.column, c.old, c.new, regfile(c.product))
            end
            push!(changes, c)
        end
    end
    for pcfg in pcfgs
        pcfg["kind"] == "reference" && continue
        name = pcfg["name"]
        append!(changes, registry_check(schemas[name], regs[name], harvested[name]; registry_file = regfile(name)))
    end

    schema_changes = filter(is_schema_change, changes)
    doc_problems = filter(is_doc_problem, changes)
    if !isempty(schema_changes)
        report_changes(io, schema_changes;
                       title = "SCHEMA CHANGES SINCE $(previous_label): $(length(schema_changes))" *
                               (accept_changes ? "  (acknowledged with --accept-changes)" : "  (NOT acknowledged)"))
    end
    if !isempty(doc_problems)
        report_changes(io, doc_problems;
                       title = "DOCUMENTATION GAPS: $(length(doc_problems))" *
                               (allow_undocumented ? "  (allowed with --allow-undocumented)" : ""))
    end

    # render
    path_map = Dict{String,String}(string(k) => string(v) for (k, v) in get(cfg, "path_map", Dict()))
    footer = "Generated $(Dates.format(now(UTC), "yyyy-mm-dd HH:MM")) UTC by <code>datamodel/generate.jl</code> " *
             "(arMADGICS.jl $(esc(git_describe(DATAMODEL_ROOT)))). Column descriptions come from the curated registry " *
             "in <code>datamodel/registry/</code>, then from dataset attributes, then from the release documents."
    ctx = (release = release, title = get(cfg, "title", release), public_prefix = get(cfg, "public_prefix", ""),
           path_map = path_map, changes = changes, accepted = accept_changes, previous_label = previous_label,
           footer = footer, registry_rel = pcfg -> "datamodel/registry/" * get(pcfg, "registry", pcfg["name"] * ".toml"),
           releases_by_product = releases_with_product(snapshot_dir, release, [p.name for p in data_products]))
    mkpath(outdir)
    rows = []
    for pcfg in pcfgs
        name = pcfg["name"]
        p = schemas[name]
        pagefile = name * ".html"
        open(joinpath(outdir, pagefile), "w") do f
            render_product(f, pcfg, p, regs[name], harvested[name], ctx)
        end
        preg = get(regs[name], "product", Dict())
        files = p.kind == "hdf5_perkey" ? get(pcfg, "public_path", "") :
                p.kind == "reference" ? "" : basename(get(pcfg, "public_path", ""))
        push!(rows, (name = name, kind = p.kind, page = pagefile, files = files,
                     ftype = Dict("hdf5_table" => "HDF5", "hdf5_perkey" => "HDF5", "csv" => "CSV", "reference" => "-")[p.kind],
                     n_rows = p.n_rows, ncols = length(p.columns), size = p.size_bytes,
                     nund = count(c -> c.kind === :undocumented && c.product == name, changes),
                     summary = get(preg, "summary", "")))
    end
    open(joinpath(outdir, "index.html"), "w") do f
        render_index(f, cfg, rows, ctx)
    end
    println(io, "  HTML written: ", joinpath(outdir, "index.html"), " (+", length(rows), " product pages)")

    nund = count(c -> c.kind === :undocumented, changes)
    norph = count(c -> c.kind === :registry_orphan, changes)
    println(io, "  summary: ", length(schema_changes), " schema change(s), ", nund, " undocumented column(s), ",
            norph, " stale registry entr(y/ies)")
    code = 0
    if !isempty(schema_changes) && !accept_changes
        println(io, "\nFAIL: unacknowledged schema changes. Review each banner above, update the registry,",
                " then re-run with --accept-changes.")
        code = 2
    elseif !isempty(doc_problems) && !allow_undocumented
        println(io, "\nFAIL: columns without documentation (or stale registry entries). Fill in the registry,",
                " or re-run with --allow-undocumented to publish with visible UNDOCUMENTED badges.")
        code = 3
    else
        println(io, "\nOK")
    end
    return code, changes
end
