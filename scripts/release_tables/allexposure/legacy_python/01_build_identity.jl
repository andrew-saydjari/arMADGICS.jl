# 01_build_identity.jl -- cache the row-identity table of an arMADGICS workup.
#
# CHAIN POSITION: allExposure step 1 of 5 (see ../README.md). First script of the whole chain.
# PURPOSE: read the string-valued row identity from the workup's full_list_info.h5 and store it as
#          compact integers (tele 1=apo 2=lco, mjd, expnum, adjfiberindx, sdss_id). Every later
#          allExposure step and the exclusion-list builder key on this file.
# INPUTS:  --full-list   <raw run dir>/full_list_info.h5   (arMADGICS raw run, e.g. raw_57618_61230)
# OUTPUTS: --out         identity.h5  (datasets adjfiberindx, expnum, sdss_id, mjd, tele)
# SOURCE:  work/2026_09_13/validation_anomaly_hunt/s1_build_identity.jl (0.4.0). Logic unchanged;
#          the two hard-coded paths became the arguments above.
#
# Usage: julia --project=scripts/release_tables -t 8 01_build_identity.jl --full-list F --out identity.h5
using HDF5, ArgParse

function parse_cli()
    s = ArgParseSettings(description = "Cache the workup row identity as integers")
    @add_arg_table! s begin
        "--full-list"; required = true; help = "full_list_info.h5 of the arMADGICS raw run"
        "--out"; required = true; help = "output identity.h5"
    end
    return parse_args(s)
end

const A = parse_cli()
const FULL = A["full-list"]
const OUT = A["out"]

println("reading full_list_info.h5 ...");
flush(stdout)
t0 = time()
h = h5open(FULL, "r")
adjfiberindx = Int16.(read(h["adjfiberindx"]))
expnum = Int32.(read(h["expnum"]))
sdss_id = read(h["sdss_id"])
println("ints done ", round(time() - t0, digits = 1), "s");
flush(stdout)
mjd_s = read(h["mjd"])
println("mjd strings done ", round(time() - t0, digits = 1), "s");
flush(stdout)
tele_s = read(h["tele"])
close(h)
println("tele strings done ", round(time() - t0, digits = 1), "s");
flush(stdout)

mjd = Vector{Int32}(undef, length(mjd_s))
Threads.@threads for i in eachindex(mjd_s)
    mjd[i] = parse(Int32, mjd_s[i])
end
tele = Vector{Int8}(undef, length(tele_s))
Threads.@threads for i in eachindex(tele_s)
    tele[i] = tele_s[i] == "apo" ? Int8(1) : Int8(2)
end
println("parsed ", round(time() - t0, digits = 1), "s")

h5open(OUT, "w") do f
    f["adjfiberindx"] = adjfiberindx
    f["expnum"] = expnum
    f["sdss_id"] = sdss_id
    f["mjd"] = mjd
    f["tele"] = tele  # 1=apo 2=lco
end
println("wrote ", OUT, " nrows=", length(mjd), " total ", round(time() - t0, digits = 1), "s")
