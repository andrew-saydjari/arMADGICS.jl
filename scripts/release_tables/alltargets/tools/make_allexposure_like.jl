# TEST HELPER ONLY (not a release tool). Build a small allExposure-like file from an almanac file:
# one row per (object exposure x fiber of that exposure's fiber table), with sdss_id, observatory,
# mjd, fiber_map_key and category. build_allTargets.jl reads only `sdss_id` from allExposure, so
# this is enough to exercise it; n_exposures from it is NOT comparable to a real allExposure
# (which applies its own row selection, s6_build_allexposure.py).
# Usage: julia --project=. tools/make_allexposure_like.jl ALMANAC.h5 OUT.h5
using HDF5
alm, out = ARGS
sid = Int64[]; obs = String[]; mjd = Int64[]; key = String[]; cat = String[]
h5open(alm, "r") do f
    for o in keys(f["raw"]), m in keys(f["raw/$o"])
        g = f["raw/$o/$m"]
        haskey(g, "exposures") && haskey(g, "fibers") || continue
        it = read(g["exposures/image_type"]); fk = read(g["exposures/fiber_map_key"])
        for (t, k) in zip(it, fk)
            t == "object" && haskey(g["fibers"], k) || continue
            s = read(g["fibers/$k/sdss_id"]); c = read(g["fibers/$k/category"])
            append!(sid, s); append!(cat, c)
            append!(obs, fill(o, length(s))); append!(mjd, fill(parse(Int, m), length(s))); append!(key, fill(k, length(s)))
        end
    end
end
h5open(out, "w") do f
    write(f, "sdss_id", sid); write(f, "observatory", obs); write(f, "mjd", mjd)
    write(f, "fiber_map_key", key); write(f, "category", cat)
    attributes(f)["description"] = "TEST allExposure-like table: object exposures x fibers from $alm (made by tools/make_allexposure_like.jl)"
end
println("wrote $out: $(length(sid)) rows, $(count(>(0), sid)) with sdss_id > 0, $(count(==(0), sid)) with sdss_id 0, $(length(unique(filter(>(0), sid)))) unique sdss_id > 0")
