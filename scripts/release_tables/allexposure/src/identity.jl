# identity.jl -- the workup row identity (stage S_IDENTITY) and the per-exposure time harvest
# (stage S_TIME). Row order of the whole table = row order of full_list_info.h5.

"Row identity of the workup, struct-of-arrays. tele: 1 = apo, 2 = lco."
struct Identity
    tele::Vector{Int8}
    mjd::Vector{Int32}
    expnum::Vector{Int32}
    adjfiberindx::Vector{Int16}
    sdss_id::Vector{Int64}
end
Base.length(id::Identity) = length(id.tele)

const TELENAME = (1 => "apo", 2 => "lco")
telename(t::Integer) = t == 1 ? "apo" : "lco"

"Read full_list_info.h5 (0.4.0 s1_build_identity.jl semantics: tele == \"apo\" -> 1, else 2)."
function read_full_list(path::AbstractString)
    h5open(path, "r") do h
        adj = Int16.(read(h["adjfiberindx"]))
        exn = Int32.(read(h["expnum"]))
        sid = Int64.(read(h["sdss_id"]))
        mjd_s = read(h["mjd"])
        tele_s = read(h["tele"])
        mjd = Vector{Int32}(undef, length(mjd_s))
        Threads.@threads for i in eachindex(mjd_s)
            mjd[i] = parse(Int32, mjd_s[i])
        end
        tele = Vector{Int8}(undef, length(tele_s))
        Threads.@threads for i in eachindex(tele_s)
            tele[i] = tele_s[i] == "apo" ? Int8(1) : Int8(2)
        end
        Identity(tele, mjd, exn, adj, sid)
    end
end

"Identity cached in the 0.4.0 identity.h5 format (interchangeable with the legacy chain)."
function load_identity(full_list::AbstractString, cache::AbstractString)
    if !isempty(cache) && isfile(cache)
        return h5open(cache, "r") do f
            Identity(Int8.(read(f["tele"])), Int32.(read(f["mjd"])), Int32.(read(f["expnum"])),
                     Int16.(read(f["adjfiberindx"])), Int64.(read(f["sdss_id"])))
        end
    end
    id = read_full_list(full_list)
    if !isempty(cache)
        h5open(cache * ".tmp", "w") do f
            f["adjfiberindx"] = id.adjfiberindx; f["expnum"] = id.expnum
            f["sdss_id"] = id.sdss_id; f["mjd"] = id.mjd; f["tele"] = id.tele
        end
        mv(cache * ".tmp", cache; force = true)
    end
    return id
end

"Exposure key used everywhere the 0.4.0 code joined on (tele, mjd, expnum)."
@inline expkey(t, m, e) = Int64(t) * 10^12 + Int64(m) * 10^6 + Int64(e)

fiber_index_almanac(t::Int8, adj::Int16) = Int16(301 - (adj - Int16(300) * (Int16(t) - Int16(1))))

# ------------------------------------------------------------------------------ time harvest

"Harvested mid-exposure times, one row per distinct workup exposure, sorted by (tele, mjd, expnum)."
struct ExposureTimes
    tele::Vector{Int8}
    mjd::Vector{Int32}
    expnum::Vector{Int32}
    mjd_mid::Vector{Float64}
    mjd_mid_precise::Vector{Float64}
    nread_total::Vector{Int32}
    source::Vector{String}
    status::Vector{String}
end

const CHIPS = ("a", "b", "c", "R", "G", "B")

"""
Harvest one (tele, night): ar1Dunical first, then ar1Dcal per chip as fallback (0.4.0 s4
semantics, including keeping partially read values when a later read throws). Runs on a worker.
"""
function harvest_night(apred::String, t::Int, night::Int, exps::Vector{Int})
    tn = telename(t)
    d = joinpath(apred, string(night))
    res = Vector{Tuple{Int,Int,Int,Float64,Float64,Int,String,String}}()
    if !isdir(d)
        for e in exps; push!(res, (t, night, e, NaN, NaN, 0, "", "no_night_dir")); end
        return res
    end
    for e in exps
        mid = NaN; midp = NaN; nread = 0; src = ""; status = "no_object_file"
        e4 = lpad(e, 4, '0')
        cands = [(joinpath(d, "ar1Dunical_$(tn)_$(night)_$(e4)_object.h5"), "unical")]
        for c in CHIPS
            push!(cands, (joinpath(d, "ar1Dcal_$(tn)_$(night)_$(e4)_$(c)_object.h5"), "cal_$c"))
        end
        for (path, tag) in cands
            isfile(path) || continue
            try
                h5open(path, "r") do h
                    g = h["metadata"]
                    haskey(g, "mjd_mid_exposure") && (mid = Float64(read(g["mjd_mid_exposure"])))
                    haskey(g, "mjd_mid_exposure_precise") && (midp = Float64(read(g["mjd_mid_exposure_precise"])))
                    haskey(g, "nread_total") && (nread = Int(read(g["nread_total"])))
                end
                src, status = tag, "ok"
            catch
                status = "read_error"
            end
            status == "ok" && break
        end
        push!(res, (t, night, e, mid, midp, nread, src, status))
    end
    return res
end

"Harvest all distinct exposures of `id` in parallel over Distributed workers (one task per night)."
function harvest_times(id::Identity, apred::AbstractString)
    keys_ = unique(sort(expkey.(id.tele, id.mjd, id.expnum)))
    tasks = Dict{Tuple{Int,Int},Vector{Int}}()
    for k in keys_
        t = Int(k ÷ 10^12); m = Int((k ÷ 10^6) % 10^6); e = Int(k % 10^6)
        push!(get!(tasks, (t, m), Int[]), e)
    end
    tl = [(t, n, sort(v)) for ((t, n), v) in sort(collect(tasks))]
    ap = String(apred)
    parts = pmap(x -> harvest_night(ap, x[1], x[2], x[3]), tl; batch_size = 4)
    rows = reduce(vcat, parts)
    sort!(rows, by = r -> (r[1], r[2], r[3]))
    length(rows) == length(keys_) || error("harvest returned $(length(rows)) rows for $(length(keys_)) exposures")
    ExposureTimes(Int8[r[1] for r in rows], Int32[r[2] for r in rows], Int32[r[3] for r in rows],
                  Float64[r[4] for r in rows], Float64[r[5] for r in rows], Int32[r[6] for r in rows],
                  String[r[7] for r in rows], String[r[8] for r in rows])
end

"Write/read the 0.4.0 exposure_times.h5 format (interchangeable with the legacy chain)."
function save_times(path, et::ExposureTimes)
    h5open(path * ".tmp", "w") do f
        f["tele"] = et.tele; f["mjd"] = et.mjd; f["expnum"] = et.expnum
        f["mjd_mid_exposure"] = et.mjd_mid; f["mjd_mid_exposure_precise"] = et.mjd_mid_precise
        f["nread_total"] = et.nread_total
        f["source"] = et.source; f["status"] = et.status
    end
    mv(path * ".tmp", path; force = true)
end
function load_times(path)
    h5open(path, "r") do f
        str(x) = [rstrip(String(s), '\0') for s in x]
        ExposureTimes(Int8.(read(f["tele"])), Int32.(read(f["mjd"])), Int32.(read(f["expnum"])),
                      read(f["mjd_mid_exposure"]), read(f["mjd_mid_exposure_precise"]),
                      Int32.(read(f["nread_total"])), str(read(f["source"])), str(read(f["status"])))
    end
end
