## Harvest column descriptions from the Markdown documents shipped with a release.
# Lowest-priority description source: used only when neither the registry nor a dataset
# `description` attribute documents a column.

"""
    harvest_markdown_tables(text) -> Dict{String,String}

Rows of Markdown tables whose first cell names one or more backticked columns, e.g.

    | `v_rad_mean` | Float64 | unweighted mean of the member `v_rad`, km/s |
    | `ra`, `dec`  | Float64 | ICRS position, deg |

map each named column to the last cell of its row.
"""
function harvest_markdown_tables(text::AbstractString)
    out = Dict{String,String}()
    for ln in split(text, '\n')
        s = strip(ln)
        startswith(s, "|") || continue
        cells = strip.(split(strip(s, '|'), '|'))
        length(cells) >= 2 || continue
        names = [m.captures[1] for m in eachmatch(r"`([A-Za-z_][A-Za-z0-9_]*)`", cells[1])]
        isempty(names) && continue
        # only accept a first cell that is nothing but backticked names and separators
        occursin(r"^[`A-Za-z0-9_,\s]+$", cells[1]) || continue
        desc = String(cells[end])
        for n in names
            haskey(out, n) || (out[n] = desc)
        end
    end
    return out
end

"""
    harvest_columns_block(text) -> Dict{String,String}

A plain-text `COLUMNS` block (as in FIBER_OFFSETS_TABLE.md): after a line that is just
`COLUMNS`, lines indented two spaces start an entry `name[, name]  description`, and
more deeply indented lines continue the previous entry. A blank line ends the block.
"""
function harvest_columns_block(text::AbstractString)
    out = Dict{String,String}()
    lines = split(text, '\n')
    i = findfirst(l -> strip(l) == "COLUMNS", lines)
    i === nothing && return out
    current = String[]
    buf = String[]
    function flush!()
        isempty(current) && return
        d = join(buf, " ")
        for n in current
            out[n] = d
        end
        empty!(current); empty!(buf)
    end
    for ln in lines[i+1:end]
        isempty(strip(ln)) && (flush!(); break)
        startswith(strip(ln), "```") && (flush!(); break)
        m = match(r"^  ([A-Za-z_][A-Za-z0-9_]*(?:, [A-Za-z_][A-Za-z0-9_]*)*)\s{2,}(.*)$", ln)
        if m !== nothing
            flush!()
            append!(current, strip.(split(m.captures[1], ",")))
            push!(buf, strip(m.captures[2]))
        else
            m2 = match(r"^  ([A-Za-z_][A-Za-z0-9_]*)\s*$", ln)
            if m2 !== nothing
                flush!()
                push!(current, m2.captures[1])
            else
                push!(buf, strip(ln))
            end
        end
    end
    flush!()
    return out
end

"""
    harvest_docs(paths) -> Dict{String,String}

Union of both harvesters over every existing file in `paths` (first hit wins).
"""
function harvest_docs(paths)
    out = Dict{String,String}()
    for p in paths
        isfile(p) || continue
        text = read(p, String)
        for d in (harvest_markdown_tables(text), harvest_columns_block(text))
            for (k, v) in d
                haskey(out, k) || (out[k] = v)
            end
        end
    end
    return out
end
