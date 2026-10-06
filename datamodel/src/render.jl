## Static HTML rendering. Self-contained CSS, no JavaScript.

const CSS = raw"""
:root { --bg:#ffffff; --fg:#1d1f21; --muted:#5f6670; --line:#d6d9de; --head:#eef1f5;
        --stripe:#f8f9fb; --link:#0b5cad; --code:#f2f3f5; --red:#b3141b; --redbg:#fdecec;
        --amber:#8a5a00; --amberbg:#fff6e0; --green:#1d6b2f; --greenbg:#e9f6ec; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#16181b; --fg:#e4e6e9; --muted:#9aa1ab; --line:#353a41; --head:#22262b;
          --stripe:#1b1e22; --link:#7db7f5; --code:#262a30; --red:#ff7b7b; --redbg:#3a1618;
          --amber:#f2c66d; --amberbg:#33290f; --green:#7fd18f; --greenbg:#142a1a; }
}
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--fg);
       font: 15px/1.5 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; }
main { max-width: 1180px; margin: 0 auto; padding: 18px 20px 60px; }
a { color: var(--link); }
h1 { font-size: 1.7em; margin: .4em 0 .3em; }
h2 { font-size: 1.3em; margin-top: 1.8em; border-bottom: 1px solid var(--line); padding-bottom: .2em; }
h3 { font-size: 1.08em; margin-top: 1.4em; }
code, pre { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: .9em;
            background: var(--code); border-radius: 3px; }
code { padding: 0 .25em; }
pre { padding: .7em .9em; overflow-x: auto; line-height: 1.4; }
.crumb { color: var(--muted); font-size: .92em; }
.crumb a { text-decoration: none; }
dl.facts { display: grid; grid-template-columns: max-content 1fr; gap: .25em 1.2em; margin: .6em 0; }
dl.facts dt { font-weight: 600; }
dl.facts dd { margin: 0; }
.tablewrap { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; margin: .6em 0 1em; font-size: .93em; }
th, td { border: 1px solid var(--line); padding: .35em .55em; vertical-align: top; text-align: left; }
th { background: var(--head); }
tbody tr:nth-child(even) > td { background: var(--stripe); }
td.name { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-weight: 600; white-space: nowrap; }
td.type, td.shape, td.num { white-space: nowrap; font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: .92em; }
td.desc p { margin: 0 0 .45em; }
td.desc p:last-child { margin-bottom: 0; }
table.bits { width: auto; margin: .4em 0 .2em; font-size: .92em; }
table.bits td, table.bits th { padding: .2em .5em; background: transparent !important; }
.null { color: var(--muted); }
.src { color: var(--muted); font-size: .82em; }
.badge { display: inline-block; font-weight: 700; font-size: .8em; padding: .05em .45em; border-radius: 3px; }
.badge.undoc { color: #fff; background: #c4161c; }
.badge.yes { color: var(--green); background: var(--greenbg); }
.badge.no { color: var(--muted); background: var(--code); }
.alert { border: 2px solid; border-radius: 6px; padding: .6em 1em; margin: 1em 0; }
.alert h2, .alert h3 { margin-top: .2em; border: 0; }
.alert.red { border-color: var(--red); background: var(--redbg); }
.alert.red h2, .alert.red h3 { color: var(--red); }
.alert.amber { border-color: var(--amber); background: var(--amberbg); }
.alert.green { border-color: var(--green); background: var(--greenbg); }
.alert ul { margin: .3em 0; padding-left: 1.3em; }
.alert li { margin: .15em 0; }
footer { color: var(--muted); font-size: .85em; margin-top: 3em; border-top: 1px solid var(--line); padding-top: .6em; }
ul.toc { columns: 2; }
"""

esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

"""
    sanitize(text, path_map) -> String

Map cluster paths to their public form (`path_map`: cluster prefix => public prefix) and
replace any other absolute cluster path with `[internal path]`, so no page shows a
filesystem location a reader cannot use.
"""
function sanitize(text::AbstractString, path_map::AbstractDict)
    s = String(text)
    for (k, v) in sort(collect(path_map), by = x -> -length(x[1]))
        s = replace(s, k => v)
    end
    s = replace(s, r"/mnt/[^\s;,)'\"]+" => "[internal path]")
    return s
end

function inline_md(s)
    t = replace(esc(s), r"`([^`]+)`" => s"<code>\1</code>")
    return replace(t, r"\*\*([^*]+)\*\*" => s"<b>\1</b>")
end

"""
    md(text; path_map) -> String

A deliberately small Markdown subset: blank-line paragraphs, `- ` bullet lists,
fenced code blocks and inline backticks. Everything else is escaped.
"""
function md(text; path_map = Dict{String,String}())
    t = sanitize(string(text), path_map)
    out = IOBuffer()
    lines = split(t, '\n')
    i = 1
    para = String[]
    function flushpara()
        isempty(para) && return
        print(out, "<p>", inline_md(join(strip.(para), " ")), "</p>")
        empty!(para)
    end
    while i <= length(lines)
        ln = lines[i]
        if startswith(strip(ln), "```")
            flushpara()
            j = i + 1
            code = String[]
            while j <= length(lines) && !startswith(strip(lines[j]), "```")
                push!(code, lines[j]); j += 1
            end
            print(out, "<pre><code>", esc(join(code, "\n")), "</code></pre>")
            i = j + 1
            continue
        elseif startswith(lstrip(ln), "- ")
            flushpara()
            print(out, "<ul>")
            while i <= length(lines) && startswith(lstrip(lines[i]), "- ")
                item = [String(lstrip(lines[i])[3:end])]
                i += 1
                while i <= length(lines) && !isempty(strip(lines[i])) &&
                      !startswith(lstrip(lines[i]), "- ") && startswith(lines[i], "  ")
                    push!(item, strip(lines[i])); i += 1
                end
                print(out, "<li>", inline_md(join(item, " ")), "</li>")
            end
            print(out, "</ul>")
            continue
        elseif isempty(strip(ln))
            flushpara()
        else
            push!(para, ln)
        end
        i += 1
    end
    flushpara()
    return String(take!(out))
end

function human_size(n::Integer)
    n <= 0 && return "—"
    units = ("B", "KB", "MB", "GB", "TB", "PB")
    x = float(n); k = 1
    while x >= 1000 && k < length(units)
        x /= 1000; k += 1
    end
    return k == 1 ? "$(n) B" : @sprintf("%.3g %s", x, units[k])
end

commas(n::Integer) = replace(string(n), r"(?<=\d)(?=(\d{3})+$)" => ",")

function page(io, title, body_fn; crumbs = Pair{String,String}[], footer = "")
    println(io, "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">")
    println(io, "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">")
    println(io, "<title>", esc(title), "</title>\n<style>", CSS, "</style>\n</head>\n<body>\n<main>")
    if !isempty(crumbs)
        print(io, "<div class=\"crumb\">Index of ")
        for (k, (label, href)) in enumerate(crumbs)
            print(io, " / ")
            isempty(href) ? print(io, esc(label)) : print(io, "<a href=\"", esc(href), "\">", esc(label), "</a>")
        end
        println(io, "</div>")
    end
    body_fn(io)
    println(io, "<footer>", footer, "</footer>\n</main>\n</body>\n</html>")
end

function change_li(c::Change)
    pc = isempty(c.column) ? c.product : "$(c.product).$(c.column)"
    code(x) = "<code>" * esc(x) * "</code>"
    if c.kind === :removed_column
        return "<b>COLUMN REMOVED:</b> $(code(pc)) — are we sure we don't want this column?"
    elseif c.kind === :new_column
        return "<b>NEW COLUMN:</b> $(code(pc)) ($(esc(c.new))) — needs a description in $(code(isempty(c.note) ? "the registry" : c.note))"
    elseif c.kind === :dtype_changed
        return "<b>DTYPE CHANGED:</b> $(code(pc)) $(code(c.old)) → $(code(c.new))"
    elseif c.kind === :shape_changed
        return "<b>SHAPE CHANGED:</b> $(code(pc)) $(code(c.old)) → $(code(c.new))"
    elseif c.kind === :removed_product
        return "<b>PRODUCT REMOVED:</b> $(code(c.product)) — are we sure this file should not ship?"
    elseif c.kind === :new_product
        return "<b>NEW PRODUCT:</b> $(code(c.product)) — needs a registry file and a page description"
    elseif c.kind === :undocumented
        return "<span class=\"badge undoc\">UNDOCUMENTED</span> $(code(pc)) ($(esc(c.new))) — add a description to $(code(c.note))"
    elseif c.kind === :registry_orphan
        return "<b>REGISTRY ENTRY WITHOUT COLUMN:</b> $(code(pc)) is in $(code(c.note)) but not in the files"
    end
    return esc(banner_line(c))
end

function alert_boxes(io, ctx; product = nothing)
    sel(cs) = product === nothing ? cs : filter(c -> c.product == product, cs)
    schema = sel(filter(is_schema_change, ctx.changes))
    docs = sel(filter(is_doc_problem, ctx.changes))
    if !isempty(schema)
        status = ctx.accepted ? "acknowledged with <code>--accept-changes</code>" :
                 "<b>not acknowledged</b>: the generator exits non-zero until these are reviewed"
        println(io, "<div class=\"alert red\"><h2>Schema changes since $(esc(ctx.previous_label))</h2>")
        println(io, "<p>", length(schema), " change(s), ", status, ".</p><ul>")
        for c in schema
            println(io, "<li>", change_li(c), "</li>")
        end
        println(io, "</ul></div>")
    end
    if !isempty(docs)
        println(io, "<div class=\"alert red\"><h2>Columns without documentation</h2><ul>")
        for c in docs
            println(io, "<li>", change_li(c), "</li>")
        end
        println(io, "</ul></div>")
    end
    if product === nothing && isempty(schema)
        if ctx.previous_label == ""
            println(io, "<div class=\"alert amber\"><b>No previous snapshot.</b> This run recorded the initial schema for ",
                    esc(ctx.release), "; the next release is checked against it.</div>")
        else
            println(io, "<div class=\"alert green\">No schema changes since ", esc(ctx.previous_label), ".</div>")
        end
    end
end

function bits_table(io, bits; note = "")
    println(io, "<table class=\"bits\"><thead><tr><th>Bit</th><th>Value</th><th>Name</th><th>Meaning</th></tr></thead><tbody>")
    for b in bits
        bit = get(b, "bit", nothing)
        val = bit === nothing ? get(b, "value", "") : 2^bit
        println(io, "<tr><td class=\"num\">", bit === nothing ? "" : bit, "</td><td class=\"num\">", val,
                "</td><td>", isempty(get(b, "name", "")) ? "" : "<code>" * esc(b["name"]) * "</code>", "</td><td>", inline_md(get(b, "meaning", "")), "</td></tr>")
    end
    println(io, "</tbody></table>")
    isempty(note) || println(io, md(note))
end

function values_table(io, values)
    println(io, "<table class=\"bits\"><thead><tr><th>Value</th><th>Meaning</th></tr></thead><tbody>")
    for v in values
        println(io, "<tr><td class=\"num\">", esc(get(v, "value", "")), "</td><td>", inline_md(get(v, "meaning", "")), "</td></tr>")
    end
    println(io, "</tbody></table>")
end

const SOURCE_LABEL = Dict(:attribute => "from the dataset's description attribute",
                          :harvested => "from the release documentation")

function desc_cell(io, r, ctx)
    print(io, "<td class=\"desc\">")
    if r.source === :undocumented
        print(io, "<span class=\"badge undoc\">UNDOCUMENTED</span>")
    else
        print(io, md(r.description; path_map = ctx.path_map))
    end
    isempty(r.null) || print(io, "<p class=\"null\"><b>Missing / sentinel:</b> ", inline_md(sanitize(r.null, ctx.path_map)), "</p>")
    isempty(r.bits) || bits_table(io, r.bits; note = r.bits_note)
    isempty(r.values) || values_table(io, r.values)
    haskey(SOURCE_LABEL, r.source) && print(io, "<p class=\"src\">(", SOURCE_LABEL[r.source], ")</p>")
    println(io, "</td>")
end

function ordered_columns(p::ProductSchema, reg)
    order = get(get(reg, "product", Dict()), "column_order", String[])
    seen = Set{String}()
    out = ColumnSchema[]
    for n in order
        c = getcolumn(p, n)
        c === nothing && continue
        push!(out, c); push!(seen, n)
    end
    for c in p.columns
        c.name in seen || push!(out, c)
    end
    return out
end

"Build a ProductSchema for a `reference` product from its registry entries."
function reference_schema(pcfg, reg)
    cols = ColumnSchema[]
    for (name, rc) in get(reg, "columns", Dict())
        push!(cols, ColumnSchema(name, get(rc, "dtype", ""), get(rc, "shape", ""),
                                 Dict{String,Any}(), Dict{String,Any}()))
    end
    order = get(get(reg, "product", Dict()), "column_order", String[])
    sort!(cols, by = c -> (something(findfirst(==(c.name), order), length(order) + 1), c.name))
    return ProductSchema(pcfg["name"], "reference", 0, 0, cols, Dict{String,Any}(), Dict{String,Any}())
end

function attr_value_html(v, ctx)
    s = v isa AbstractVector ? "[" * join(string.(v), ", ") * "]" : string(v)
    return inline_md(sanitize(s, ctx.path_map))
end

function render_product(io, pcfg, p::ProductSchema, reg, harvested, ctx)
    preg = get(reg, "product", Dict())
    name = p.name
    title = get(preg, "title", name)
    crumbs = ["datamodel" => "index.html", ctx.release => "index.html", name => ""]
    page(io, "Data Model: $name", io -> begin
        println(io, "<h1>Data Model: ", esc(name), "</h1>")
        title == name || println(io, "<p><i>", esc(title), "</i></p>")
        alert_boxes(io, ctx; product = name)
        get(reg, "_missing", false) && println(io, "<div class=\"alert red\"><h3>No registry file</h3><p>Expected <code>",
                                                esc(ctx.registry_rel(pcfg)), "</code>. Every column below relies on file attributes or harvested docs.</p></div>")
        println(io, "<h2 id=\"general\">General Description</h2>")
        gd = get(preg, "general_description", "")
        println(io, isempty(gd) ? "<p><span class=\"badge undoc\">UNDOCUMENTED</span></p>" : md(gd; path_map = ctx.path_map))
        println(io, "<dl class=\"facts\">")
        nc = get(preg, "naming_convention", "")
        isempty(nc) || println(io, "<dt>Naming Convention</dt><dd>", md(nc), "</dd>")
        loc = get(pcfg, "public_path", "")
        isempty(loc) || println(io, "<dt>Location</dt><dd><code>", esc(ctx.public_prefix * "/" * loc), "</code></dd>")
        if p.kind != "reference"
            println(io, "<dt>Approximate Size</dt><dd>", human_size(p.size_bytes), "</dd>")
        end
        ft = get(preg, "file_type", Dict("hdf5_table" => "HDF5", "hdf5_perkey" => "HDF5", "csv" => "CSV",
                                         "reference" => "reference (no file in this release)")[p.kind])
        println(io, "<dt>File Type</dt><dd>", esc(ft), "</dd>")
        p.n_rows > 0 && println(io, "<dt>Rows</dt><dd>", commas(p.n_rows), get(preg, "row_definition", "") == "" ? "" :
                                     " — " * inline_md(preg["row_definition"]), "</dd>")
        p.kind == "reference" || println(io, "<dt>Columns</dt><dd>", length(p.columns), "</dd>")
        for (k, label) in (("written_by", "Written by"), ("read_by", "Read by"))
            v = get(preg, k, String[])
            isempty(v) && continue
            println(io, "<dt>", label, "</dt><dd>", join(inline_md.(v), "<br>"), "</dd>")
        end
        rels = get(ctx.releases_by_product, name, String[])
        isempty(rels) || println(io, "<dt>Releases</dt><dd>", esc(join(rels, ", ")), "</dd>")
        println(io, "</dl>")

        println(io, "<h2>Sections</h2><ul>")
        has_attrs = !isempty(p.file_attrs) || !isempty(get(p.extra, "hdr_attrs", Dict()))
        has_attrs && println(io, "<li><a href=\"#attrs\">File attributes</a></li>")
        println(io, "<li><a href=\"#columns\">", p.kind == "reference" ? "Contents" : "Columns", "</a></li>")
        haskey(preg, "notes") && println(io, "<li><a href=\"#notes\">Notes</a></li>")
        haskey(preg, "examples") && println(io, "<li><a href=\"#examples\">Examples</a></li>")
        println(io, "</ul>")

        if has_attrs
            println(io, "<h2 id=\"attrs\">File attributes</h2>")
            ad = get(preg, "attrs_description", "")
            isempty(ad) || println(io, md(ad))
            attr_notes = get(preg, "attr_notes", Dict())
            if !isempty(p.file_attrs)
                println(io, "<div class=\"tablewrap\"><table><thead><tr><th>Key</th><th>Value</th><th>Comment</th></tr></thead><tbody>")
                for k in sort(collect(keys(p.file_attrs)))
                    println(io, "<tr><td class=\"name\">", esc(k), "</td><td>", attr_value_html(p.file_attrs[k], ctx),
                            "</td><td>", inline_md(get(attr_notes, k, "")), "</td></tr>")
                end
                println(io, "</tbody></table></div>")
            end
            hdr = get(p.extra, "hdr_attrs", Dict())
            if !isempty(hdr)
                println(io, "<h3>Provenance: attributes of the <code>hdr</code> dataset</h3>")
                println(io, "<div class=\"tablewrap\"><table><thead><tr><th>Key</th><th>Value</th><th>Comment</th></tr></thead><tbody>")
                for k in sort(collect(keys(hdr)))
                    println(io, "<tr><td class=\"name\">", esc(k), "</td><td>", attr_value_html(hdr[k], ctx),
                            "</td><td>", inline_md(get(attr_notes, "hdr." * k, "")), "</td></tr>")
                end
                println(io, "</tbody></table></div>")
                mm = get(p.extra, "hdr_mismatch", String[])
                if !isempty(mm)
                    println(io, "<div class=\"alert red\"><b>Provenance differs between files:</b> ", esc(join(mm, ", ")), "</div>")
                else
                    println(io, "<p class=\"src\">Identical in all ", length(p.columns), " files.</p>")
                end
            end
        end

        println(io, "<h2 id=\"columns\">", p.kind == "reference" ? "Contents" : "Columns", "</h2>")
        cd_ = get(preg, "columns_description", "")
        isempty(cd_) || println(io, md(cd_))
        perkey = p.kind == "hdf5_perkey"
        hdf = p.kind in ("hdf5_table", "hdf5_perkey")
        println(io, "<div class=\"tablewrap\"><table><thead><tr><th>Name</th><th>Type</th>")
        if hdf
            println(io, "<th>Shape<br><span class=\"src\">Julia · h5py</span></th>")
        elseif p.kind == "reference"
            println(io, "<th>Shape</th>")
        end
        println(io, "<th>Units</th><th>Description</th>")
        perkey && println(io, "<th>File size</th><th>In ", esc(ctx.release), " release dir</th>")
        println(io, "</tr></thead><tbody>")
        for c in ordered_columns(p, reg)
            r = resolve_column(c.name, reg, c, harvested)
            println(io, "<tr id=\"col-", esc(c.name), "\"><td class=\"name\">", esc(c.name), "</td><td class=\"type\">", esc(c.dtype), "</td>")
            if hdf
                println(io, "<td class=\"shape\">", esc(c.shape), "<br>", esc(h5py_shape(c.shape)), "</td>")
            elseif p.kind == "reference"
                println(io, "<td class=\"shape\">", esc(c.shape), "</td>")
            end
            println(io, "<td>", inline_md(r.units), "</td>")
            desc_cell(io, r, ctx)
            if perkey
                println(io, "<td class=\"num\">", human_size(get(c.extra, "size_bytes", 0)), "</td>")
                inrel = get(c.extra, "in_release_dir", false)
                println(io, "<td>", inrel ? "<span class=\"badge yes\">yes</span>" : "<span class=\"badge no\">no</span>", "</td>")
            end
            println(io, "</tr>")
        end
        println(io, "</tbody></table></div>")
        if haskey(preg, "notes")
            println(io, "<h2 id=\"notes\">Notes</h2>", md(preg["notes"]; path_map = ctx.path_map))
        end
        if haskey(preg, "examples")
            println(io, "<h2 id=\"examples\">Examples</h2>", md(preg["examples"]; path_map = ctx.path_map))
        end
    end; crumbs = crumbs, footer = ctx.footer)
end

function render_index(io, cfg, rows, ctx)
    page(io, "$(ctx.title) data model", io -> begin
        println(io, "<h1>", esc(ctx.title), " — data model</h1>")
        intro = get(cfg, "description", "")
        isempty(intro) || println(io, md(intro; path_map = ctx.path_map))
        alert_boxes(io, ctx)
        ncol = sum(r -> r.ncols, filter(r -> r.kind != "reference", rows); init = 0)
        nund = count(c -> c.kind === :undocumented, ctx.changes)
        println(io, "<dl class=\"facts\"><dt>Release</dt><dd>", esc(ctx.release), "</dd>")
        println(io, "<dt>Location</dt><dd><code>", esc(ctx.public_prefix), "/</code></dd>")
        println(io, "<dt>Products</dt><dd>", count(r -> r.kind != "reference", rows), " (", ncol, " columns)</dd>")
        println(io, "<dt>Undocumented columns</dt><dd>", nund == 0 ? "0" : "<span class=\"badge undoc\">$(nund)</span>", "</dd>")
        println(io, "<dt>Compared against</dt><dd>", isempty(ctx.previous_label) ? "nothing (first snapshot)" : esc(ctx.previous_label), "</dd></dl>")

        println(io, "<h2>Products</h2><div class=\"tablewrap\"><table><thead><tr><th>Product</th><th>File(s)</th><th>Type</th><th>Rows</th><th>Columns</th><th>Size</th><th>Undocumented</th><th>Summary</th></tr></thead><tbody>")
        for r in filter(r -> r.kind != "reference", rows)
            println(io, "<tr><td class=\"name\"><a href=\"", esc(r.page), "\">", esc(r.name), "</a></td><td><code>", esc(r.files),
                    "</code></td><td>", esc(r.ftype), "</td><td class=\"num\">", r.n_rows > 0 ? commas(r.n_rows) : "—",
                    "</td><td class=\"num\">", r.ncols, "</td><td class=\"num\">", human_size(r.size),
                    "</td><td class=\"num\">", r.nund == 0 ? "0" : "<span class=\"badge undoc\">$(r.nund)</span>",
                    "</td><td>", inline_md(r.summary), "</td></tr>")
        end
        println(io, "</tbody></table></div>")
        refs = filter(r -> r.kind == "reference", rows)
        if !isempty(refs)
            println(io, "<h2>Reference pages</h2><ul>")
            for r in refs
                println(io, "<li><a href=\"", esc(r.page), "\">", esc(r.name), "</a> — ", inline_md(r.summary), "</li>")
            end
            println(io, "</ul>")
        end
        docs = get(cfg, "documents", Any[])
        if !isempty(docs)
            println(io, "<h2>Accompanying documents</h2><div class=\"tablewrap\"><table><thead><tr><th>File</th><th>Contents</th></tr></thead><tbody>")
            for d in docs
                println(io, "<tr><td><code>", esc(ctx.public_prefix * "/" * d["public_path"]), "</code></td><td>", inline_md(get(d, "summary", "")), "</td></tr>")
            end
            println(io, "</tbody></table></div>")
        end
        conv = get(cfg, "conventions", "")
        isempty(conv) || println(io, "<h2>Conventions shared by all products</h2>", md(conv; path_map = ctx.path_map))
    end; crumbs = ["datamodel" => "index.html", ctx.release => ""], footer = ctx.footer)
end
