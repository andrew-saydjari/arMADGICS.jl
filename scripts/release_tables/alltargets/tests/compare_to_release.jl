# Regression: compare a released allTargets (0.4.0) with a file built by the 0.4.1 builder from the
# same inputs. For every dataset of the release (top level and inside groups): present, same HDF5
# datatype (H5Tequal: class, size, sign, string padding/cset/vlen), same shape, same chunking and
# filters, same values (NaN == NaN), same dataset attributes. Then lists datasets only in the new
# file, and file attributes that differ / are new. One dataset in memory at a time. Read-only.
# Usage: julia --project=. tests/compare_to_release.jl RELEASE.h5 NEW.h5
using HDF5
old, new = ARGS
same(x, y) = isequal(x, y) || (x isa AbstractFloat && y isa AbstractFloat && isnan(x) && isnan(y))
function walk(g, prefix = "")
    out = String[]
    for k in keys(g)
        o = g[k]; p = prefix * k
        o isa HDF5.Group ? append!(out, walk(o, p * "/")) : push!(out, p)
    end
    out
end
layout(d) = (dc = HDF5.get_create_properties(d);
    (try HDF5.get_chunk(d) catch; () end, string(dc.filters)))
dattrs(d) = Dict(k => read(attributes(d)[k]) for k in keys(attributes(d)))
h5open(old, "r") do fo; h5open(new, "r") do fn
    po = walk(fo); pn = walk(fn)
    println("release datasets: $(length(po)); new datasets: $(length(pn))")
    nok = 0; probs = String[]
    for p in po
        if !haskey(fn, p); push!(probs, "$p: MISSING in new"); continue; end
        a, b = fo[p], fn[p]
        HDF5.API.h5t_equal(datatype(a), datatype(b)) > 0 || push!(probs, "$p: datatype differs ($(eltype(a)) vs $(eltype(b)))")
        size(a) == size(b) || (push!(probs, "$p: shape $(size(a)) vs $(size(b))"); continue)
        layout(a) == layout(b) || push!(probs, "$p: layout $(layout(a)) vs $(layout(b))")
        va, vb = read(a), read(b)
        nd = count(.!same.(va, vb)); nd == 0 || push!(probs, "$p: $nd values differ")
        da, db = dattrs(a), dattrs(b)
        da == db || push!(probs, "$p: dataset attributes differ: release $(da) vs new $(db)")
        nd == 0 && da == db && (nok += 1)
        va = vb = nothing; GC.gc()
    end
    println("release datasets identical in dtype, shape, values and attributes: $nok / $(length(po))")
    println("problems: $(length(probs))"); foreach(x -> println("  ", x), probs)
    extra = setdiff(pn, po)
    println("datasets only in new ($(length(extra))): ", join(extra, ", "))
    # group attributes
    for g in ("carton_pks", "carton_table")
        dattrs(fo[g]) == dattrs(fn[g]) || println("group $g attributes differ: $(dattrs(fo[g])) vs $(dattrs(fn[g]))")
    end
    ao, an = dattrs(fo), dattrs(fn)
    println("file attributes: release $(length(ao)), new $(length(an))")
    for k in sort(collect(keys(ao)))
        haskey(an, k) || (println("  attr $k MISSING in new"); continue)
        isequal(ao[k], an[k]) || println("  attr $k differs: release $(repr(ao[k])[1:min(end,140)]) | new $(repr(an[k])[1:min(end,140)])")
    end
    println("  new-only attrs: ", join(sort(setdiff(collect(keys(an)), collect(keys(ao)))), ", "))
end; end
