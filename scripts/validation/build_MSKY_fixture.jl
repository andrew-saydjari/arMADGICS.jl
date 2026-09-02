# Build the M-SKY validation fixtures from the M1/M2/M3 fixture (build_M123_fixture.jl):
# apo 58588 exposure 0011 with poisoned SKY fibers, to demonstrate/regress the
# sky-prior poisoning path (getSkyRough/combine_sky_fibers).
#
# Variants (each a pipeline-shaped mini reduxBase):
#   poison_partialnan/ : sky fiberindx 102 gets NaN flux on every other pixel.
#                        Pre-guard this fiber PASSES the scale z-cut and its NaNs
#                        enter VLocSkyLines -> every fiber's RV on the exposure is
#                        poisoned. Post-guard it is excluded + flagged.
#   poison_allnan/     : sky fiberindx 102 fully NaN with good mask (the A4 shape).
#                        Pre-guard this fiber is dropped only by the ACCIDENT that
#                        its NaN scale fails the z-cut comparison; post-guard it is
#                        excluded explicitly + flagged. Pre-guard output on this
#                        fixture = the 34-sky-fiber no-poison reference that the
#                        post-guard poison_partialnan run must match bit-identically.
#   poison_lowsky/     : all sky fibers except two fully NaN -> post-guard the
#                        surviving count falls below SKY_MIN_FIBERS and the sky-line
#                        component is skipped (flagged), the batch still completes.
#
# Usage: julia --project=. scripts/validation/build_MSKY_fixture.jl <m123_fixture_dir> <out_dir>
# Reads (read-only): <m123_fixture_dir>/apred/58588/ar1Dunical_apo_58588_0011_object.h5

using HDF5

length(ARGS) >= 2 || error("usage: build_MSKY_fixture.jl <m123_fixture_dir> <out_dir>")
m123_dir = abspath(ARGS[1])
out_dir = abspath(ARGS[2])

const RELPATH = joinpath("apred", "58588", "ar1Dunical_apo_58588_0011_object.h5")
src_file = joinpath(m123_dir, RELPATH)

# sky fiberindxs for apo 58588 exposure 11 (from the almanac fibers table via
# get_fibTargDict; category "sky*"): verified 2026-09-02 against
# allobs_57600_61160.h5 (35 fibers)
const SKYFIBS = [1, 6, 7, 11, 12, 19, 25, 31, 42, 55, 56, 59, 60, 61, 72, 78, 91, 102,
    114, 120, 127, 132, 144, 175, 193, 199, 217, 223, 234, 235, 241, 247, 253, 270, 295]
const POISON_FIB = 102 # the sky fiber adjacent to the M123 run fibers
const LOWSKY_KEEP = [1, 6] # survivors in the low-count variant

flux0, ivar0, msk0 = h5open(src_file) do f
    read(f["flux_1d"]), read(f["ivar_1d"]), read(f["mask_1d"])
end
println("source: $src_file  flux size ", size(flux0))

function write_variant(name, flux, ivar, msk)
    dst = joinpath(out_dir, name, RELPATH)
    mkpath(dirname(dst))
    rm(dst; force=true)
    h5open(dst, "w") do f
        f["flux_1d"] = flux
        f["ivar_1d"] = ivar
        f["mask_1d"] = msk
    end
    println("wrote $dst")
end

# partial-NaN sky fiber (poisons pre-guard)
flux = copy(flux0)
flux[1:2:end, POISON_FIB] .= NaN
write_variant("poison_partialnan", flux, ivar0, msk0)

# all-NaN sky fiber, A4 shape (silently dropped pre-guard; flagged post-guard)
flux = copy(flux0)
msk = copy(msk0)
flux[:, POISON_FIB] .= NaN
msk[:, POISON_FIB] .= true
write_variant("poison_allnan", flux, ivar0, msk)

# all but two sky fibers dead -> post-guard sky-line component skipped + flagged
flux = copy(flux0)
msk = copy(msk0)
for fib in setdiff(SKYFIBS, LOWSKY_KEEP)
    flux[:, fib] .= NaN
    msk[:, fib] .= true
end
write_variant("poison_lowsky", flux, ivar0, msk)
