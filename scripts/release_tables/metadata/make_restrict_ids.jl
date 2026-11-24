# make_restrict_ids.jl -- the sdss_id set that `almanac add metadata --restrict-to` queries.
#
# CHAIN POSITION: metadata step 1 of 2 (after the allExposure build, before build_metadata.sh).
# USAGE:   julia --project=scripts/release_tables make_restrict_ids.jl ALLEXPOSURE.h5 OUT.npy
# INPUT:   ARGS[1]  allExposure h5 (reads `sdss_id`)
# OUTPUT:  ARGS[2]  .npy of the unique sdss_id > 0, sorted ascending, int64 (numpy format 1.0,
#          little-endian '<i8', C order), byte-identical to what np.save writes.
# For 0.4.0 this was 4,019,650 ids. NOTE for 0.4.1: sdss_id 0 (the solar-system target, Vesta) is
# excluded by `> 0` here.
#
# Julia port of legacy_python/make_restrict_ids.py (2026-11-24); reproduces it exactly (the .npy is byte-identical
# to the Python output and to the 0.4.0 reference file). The .npy writer below is numpy's
# format.write_array for a 1-D little-endian int64 array: header dict with numpy's trailing
# growth-axis spaces, padded with spaces + '\n' so that the data start on a 64-byte boundary.
using HDF5

"Unique positive sdss_id, sorted ascending (numpy: np.unique(s[s > 0]) after astype(int64))."
function restrict_ids(path::AbstractString)
    s = h5open(f -> read(f["sdss_id"]), path, "r")
    return sort!(unique!(Int64[Int64(x) for x in s if x > 0]))
end

"numpy's .npy header (format 1.0) for a 1-D C-order little-endian int64 array of length n."
function npy_header(n::Integer)
    shape = "($(n),)"
    h = "{'descr': '<i8', 'fortran_order': False, 'shape': $(shape), }"
    h *= " "^(21 - length(string(n)))            # numpy GROWTH_AXIS_MAX_DIGITS = 21
    hlen = length(h) + 1                          # + the final '\n'
    padlen = 64 - ((6 + 2 + 2 + hlen) % 64)       # magic(6) + version(2) + uint16 length(2); 64 if aligned
    hlen + padlen > typemax(UInt16) && error("header too long for npy format 1.0")
    io = IOBuffer()
    write(io, UInt8[0x93], b"NUMPY", UInt8[0x01, 0x00], htol(UInt16(hlen + padlen)))
    write(io, h, " "^padlen, "\n")
    return take!(io)
end

function write_npy(path::AbstractString, v::Vector{Int64})
    endswith(path, ".npy") || (path *= ".npy")   # np.save appends the extension
    open(path, "w") do io
        write(io, npy_header(length(v)))
        write(io, htol.(v))
    end
    return nothing
end

function main(args)
    length(args) == 2 || error("usage: make_restrict_ids.jl ALLEXPOSURE.h5 OUT.npy")
    t = @elapsed begin
        ids = restrict_ids(args[1])
        write_npy(args[2], ids)
    end
    println("$(length(ids)) unique sdss_id > 0 -> $(args[2])")
    println("work time: ", round(t; digits = 2), " s")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
