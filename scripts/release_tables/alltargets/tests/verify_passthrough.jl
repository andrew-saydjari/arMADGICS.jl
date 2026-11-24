# Check that every /meta column reached an allTargets file built by build_allTargets.jl (0.4.1).
# For each /meta dataset: present in the output (sdss_id as the key; sdss5_target_flags rebuilt,
# compared only for width), same element type (strings: same text), same per-row shape, same
# values row-by-row after joining on sdss_id (NaN == NaN), same description. Prints the
# solar-system (sdss_id <= 0) rows. Read-only on both files.
# Usage: julia --project=. tests/verify_passthrough.jl META.h5 ALLTARGETS.h5
using HDF5, JSON
meta, out = ARGS
same(x, y) = isequal(x, y) || (x isa AbstractFloat && y isa AbstractFloat && isnan(x) && isnan(y))
h5open(meta, "r") do fm; h5open(out, "r") do fo
    g = fm["meta"]
    mnames = haskey(attributes(g), "columns") ? String.(read(attributes(g)["columns"])) : collect(keys(g))
    msid = read(g["sdss_id"]); osid = read(fo["sdss_id"])
    rows = [searchsortedfirst(osid, s) for s in msid]          # output row of each /meta row (or past end)
    inout = [r <= length(osid) && osid[r] == s for (r, s) in zip(rows, msid)]
    println("meta rows $(length(msid)), output rows $(length(osid)), meta rows present in output $(count(inout))")
    sel = findall(inout); orow = rows[sel]
    bad = String[]; nchk = 0
    for c in mnames
        c == "sdss_id" && continue
        if c == "sdss5_target_flags"
            println("  sdss5_target_flags: /meta $(size(g[c])) (rebuilt in output: $(size(fo[c])), file attr meta_sdss5_target_flags_nbytes = $(read(attributes(fo)["meta_sdss5_target_flags_nbytes"])))")
            continue
        end
        haskey(fo, c) || (push!(bad, "$c: MISSING"); continue)
        mv = read(g[c]); ov = read(fo[c])
        eltype(mv) == eltype(ov) || (mv isa Vector{String} && ov isa Vector{String}) || push!(bad, "$c: eltype $(eltype(mv)) -> $(eltype(ov))")
        size(mv)[1:end-1] == size(ov)[1:end-1] || push!(bad, "$c: per-row shape $(size(mv)) -> $(size(ov))")
        M = reshape(mv, :, size(mv, ndims(mv)))[:, sel]; O = reshape(ov, :, size(ov, ndims(ov)))[:, orow]
        all(same.(M, O)) || push!(bad, "$c: VALUES differ in $(count(.!same.(M, O))) entries")
        dm = haskey(attributes(g[c]), "description") ? read(attributes(g[c])["description"]) : ""
        do_ = read(attributes(fo[c])["description"])
        startswith(do_, dm) || c == "catalogid" || push!(bad, "$c: description '$dm' -> '$do_'")
        nchk += 1
    end
    println("checked $nchk pass-through columns of $(length(mnames)) /meta columns; problems: $(length(bad))")
    foreach(b -> println("  ", b), bad)
    fa = attributes(fo)
    for k in ("release", "meta_extended", "meta_extended_source", "meta_almanac_version", "meta_almanac_commit",
              "meta_source", "n_meta_passthrough_columns", "meta_n_columns", "multidim_columns",
              "solar_system_sdss_id", "n_solar_system_rows", "row_definition")
        println("  attr $k = ", haskey(fa, k) ? read(fa[k]) : "(absent)")
    end
    println("  attr meta_extended_column_groups = ", haskey(fa, "meta_extended_column_groups") ? read(fa["meta_extended_column_groups"]) : "(absent)")
    co = read(fa["column_order"]); println("  output columns: $(length(co)); added since 0.4.0: $(length(read(fa["columns_added_since_040"])))")
    for i in findall(<=(0), osid)
        println("solar-system row sdss_id=$(osid[i]):")
        nonfill = String[]
        for c in co
            d = fo[c]; v = ndims(d) == 1 ? d[i] : d[:, i]
            isfill = v isa AbstractFloat ? isnan(v) : v isa AbstractString ? isempty(v) : v isa AbstractVector ? all(iszero, v) : (v == -1 || v == 0)
            isfill || push!(nonfill, "$c=$(repr(v))")
        end
        println("  non-fill values: ", join(nonfill, ", "))
        n = length(read(fo["carton_pks/sdss_id"]) |> x -> filter(==(osid[i]), x))
        println("  carton pairs: $n; metadata_status = $(fo["metadata_status"][i]) (bits: ",
            join([d["name"] for d in JSON.parse(read(fa["status_bits"])) if fo["metadata_status"][i] & (UInt8(1) << d["bit"]) != 0], ", "), ")")
    end
end; end
