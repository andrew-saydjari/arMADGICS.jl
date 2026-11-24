"""
    AllExposure

Builder of the APOGEE `allExposure` release table: one row per arMADGICS workup row
(exposure x fiber), joining almanac metadata, harvested mid-exposure times, barycentric and
radial velocities, quality flags, Galactic coordinates, fiber zero points and the FPS off-target
flag in ONE pass. Pure Julia by default; `Legacy040Numerics` reproduces the 0.4.0 table bit for bit
(see README).

Entry point: `build(BuildConfig(...))`, or `bin/build_allexposure.jl` from the command line.
"""
module AllExposure

using HDF5, CSV, JSON, Printf, Dates, Distributed
import ERFA_jll, SHA

export BuildConfig, build, JuliaNumerics, Legacy040Numerics, SCHEMA

include("schema.jl")
include("h5util.jl")
include("identity.jl")
include("almanac.jl")
include("derived.jl")
include("bary.jl")
include("build.jl")

end
