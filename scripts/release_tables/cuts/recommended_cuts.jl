# recommended_cuts.jl -- the default quality mask of the released allExposure table (Julia).
# Same definition as recommended_cuts.py / RECOMMENDED_CUTS.md:
#   qa_bitmask == 0, ingestBit == 0, RV_flag == 0, bary_flag == 0, isfinite(v_rad), offtarget == 0
# `include_offtarget = false` gives the 5-cut mask used by allVisitRV and the fiber-offset fit.
# Usage: julia --project=scripts/release_tables recommended_cuts.jl ALLEXPOSURE.h5
using HDF5

function recommended_mask(path::AbstractString; include_offtarget::Bool = true)
    h5open(path, "r") do f
        keep = @. (read(f["qa_bitmask"]) == 0) & (read(f["ingestBit"]) == 0) &
                  (read(f["RV_flag"]) == 0) & (read(f["bary_flag"]) == 0) & isfinite(read(f["v_rad"]))
        include_offtarget && (keep .&= read(f["offtarget"]) .== 0)
        keep
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    p = ARGS[1]
    println("kept without offtarget cut: ", count(recommended_mask(p; include_offtarget = false)))
    println("kept with offtarget cut   : ", count(recommended_mask(p)))
end
