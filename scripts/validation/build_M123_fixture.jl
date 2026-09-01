# Build the M1/M2/M3 validation fixture from real 2026_05_01 outputs.
#
# Extracts apo 58588 exposure 0011 (dead fiber_id 211 -> fiberindx/adjfib 90,
# domeflat relthrpt 1.4e-4) into a pipeline-shaped mini reduxBase, and injects
# two synthetic pathological fibers to exercise the M2 guards:
#   fiberindx 100 -> all-masked spectrum
#   fiberindx 101 -> NaN flux presented with an all-good mask (the A4 shape)
# Also writes a compact ~10-fiber x 8700 extract (fixture_small.h5).
#
# Usage: julia --project=. scripts/validation/build_M123_fixture.jl <fixture_dir>
# Reads (read-only): /mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/apred/58588/

using HDF5

src_file = "/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/apred/58588/ar1Dunical_apo_58588_0011_object.h5"
fixture_dir = length(ARGS) >= 1 ? ARGS[1] : error("usage: build_M123_fixture.jl <fixture_dir>")

# fibers the driver runs (fiberindx = adjfiberindx for apo):
# 85-89, 92, 93 healthy sci; 90 dead (fiber_id 211); 100, 101, 103 injected pathologies
const RUN_FIBERS = [85, 86, 87, 88, 89, 90, 92, 93, 100, 101, 103]
const ALLMASKED_FIBER = 100
const NANFLUX_FIBER = 101
const PARTIALNAN_FIBER = 103 # NaN flux on half the good pixels reaches the solve pre-fix (NOTE: 102 is a sky fiber - do not poison it)

flux, ivar, msk = h5open(src_file) do f
    read(f["flux_1d"]), read(f["ivar_1d"]), read(f["mask_1d"])
end
println("source: $src_file  flux size ", size(flux))

# inject the two synthetic bad fibers (overwriting two healthy sci fibers)
msk_out = copy(msk)
flux_out = copy(flux)
msk_out[:, ALLMASKED_FIBER] .= false
flux_out[:, NANFLUX_FIBER] .= NaN
msk_out[:, NANFLUX_FIBER] .= true # A4 shape: NaN flux, good mask
flux_out[1:2:end, PARTIALNAN_FIBER] .= NaN # partial A4 shape: NaN flux on good pixels

outpath = joinpath(fixture_dir, "apred", "58588")
mkpath(outpath)
outfile = joinpath(outpath, "ar1Dunical_apo_58588_0011_object.h5")
rm(outfile; force=true)
h5open(outfile, "w") do f
    f["flux_1d"] = flux_out
    f["ivar_1d"] = ivar
    f["mask_1d"] = msk_out
end
println("wrote pipeline-shaped fixture: $outfile")

# compact extract for archival/unit-style use
smallfile = joinpath(fixture_dir, "fixture_small.h5")
rm(smallfile; force=true)
h5open(smallfile, "w") do f
    f["fiberindx"] = RUN_FIBERS
    f["flux_1d"] = flux_out[:, RUN_FIBERS]
    f["ivar_1d"] = ivar[:, RUN_FIBERS]
    f["mask_1d"] = collect(msk_out[:, RUN_FIBERS])
    attrs(f)["source"] = src_file
    attrs(f)["note"] = "fiber 90 = dead fiber_id 211; fiber 100 = injected all-masked; fiber 101 = injected NaN-flux/good-mask (A4 shape); fiber 103 = injected partial NaN flux on good pixels"
end
println("wrote compact fixture: $smallfile")
