# recommended_cuts.jl -- the default quality mask of the released allExposure table (Julia).
#
# Definition = data/work/mwm/0.4.0/public/RECOMMENDED_CUTS.md ("The recipe"):
#   qa_bitmask == 0, ingestBit == 0, RV_flag == 0, bary_flag == 0, isfinite(v_rad), offtarget == 0
# `include_offtarget = false` gives the 5-cut mask that allVisitRV membership and the fiber-offset
# fit used (0.4.0 built those before the offtarget column existed).
#
# CHAIN POSITION: check (step 20), on the final allExposure.
# USAGE:   julia --project=scripts/release_tables recommended_cuts.jl ALLEXPOSURE.h5 [--expect N5 N6]
#          prints the kept counts; with --expect, exits non-zero unless they equal N5 (5 cuts) and
#          N6 (6 cuts). 0.4.0: N5 = 25516821, N6 = 25207025.
# INPUT:   allExposure h5 (qa_bitmask, ingestBit, RV_flag, bary_flag, v_rad, offtarget)
# OUTPUT:  stdout only (the work time goes to stderr).
#
# Julia port of legacy_python/recommended_cuts.py (2026-11-24); reproduces it exactly (same counts on 0.4.0, same
# stdout, same --expect exit status). Columns are read one at a time and folded into the mask.
# (The earlier version of this file broadcast `read` over the column elements via `@.` and failed
# with a MethodError; it never ran.)
using HDF5, ArgParse

"AND `keep` with `pred(x)` elementwise over column `name` (read one column at a time)."
function and_column!(keep::BitVector, f::HDF5.File, name::AbstractString, pred::F) where {F}
    col = read(f[name])                  # element type known only at run time: function barrier
    length(col) == length(keep) || error("$name: length $(length(col)) != $(length(keep))")
    return _and!(keep, col, pred)
end

function _and!(keep::BitVector, col::Vector{T}, pred::F) where {T,F}
    @inbounds for i in eachindex(col)
        keep[i] = keep[i] & pred(col[i])
    end
    return keep
end

"Boolean keep-mask (true = keep), one entry per allExposure row."
function recommended_mask(path::AbstractString; include_offtarget::Bool = true)
    h5open(path, "r") do f
        keep = trues(length(f["v_rad"]))
        and_column!(keep, f, "qa_bitmask", iszero)
        and_column!(keep, f, "ingestBit", iszero)
        and_column!(keep, f, "RV_flag", iszero)
        and_column!(keep, f, "bary_flag", iszero)
        and_column!(keep, f, "v_rad", isfinite)
        include_offtarget && and_column!(keep, f, "offtarget", iszero)
        keep
    end
end

function parse_cli(args)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "allexposure"
            required = true
        "--expect"
            nargs = 2; arg_type = Int
    end
    return parse_args(args, s)
end

function main(args)
    a = parse_cli(args)
    k5 = k6 = 0
    t = @elapsed begin
        k5 = count(recommended_mask(a["allexposure"]; include_offtarget = false))
        k6 = count(recommended_mask(a["allexposure"]; include_offtarget = true))
    end
    println("kept without offtarget cut: $k5\nkept with offtarget cut   : $k6")
    println(stderr, "work time: ", round(t; digits = 2), " s")
    ex = a["expect"]
    if !isempty(ex) && (k5, k6) != (ex[1], ex[2])
        println("MISMATCH: expected [$(ex[1]), $(ex[2])]")
        return 1
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
