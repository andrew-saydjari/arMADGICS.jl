# Tests for the per-exposure sky-bundle cache (src/skyCache.jl).
#
# These exercise the three hazards the cache has to survive in the pass-1 bulk run:
#   1. COLLISION  - the same (mjd, expnum) exists at BOTH observatories (31,307 such
#                   pairs measured over 57618-61230). Two telescopes' bundles must never
#                   alias, and the proof here is content-level, not just key-level.
#   2. RACE       - many workers hit the same exposure concurrently. A reader must never
#                   observe a partially written entry.
#   3. STALENESS  - a changed input must NOT be served from an old entry.
# plus the fail-open contract: every cache failure degrades to recomputation.

import Distributed
using Distributed: addprocs, rmprocs, pmap

# ---------------------------------------------------------------------------------
# synthetic fixture: a minimal reduxBase + almanac that get_sky_fiber_indices accepts
# ---------------------------------------------------------------------------------

# the production wavelength grid (pipeline.jl defines this global; the test suite does not)
if !@isdefined(logUniWaveAPOGEE)
    logUniWaveAPOGEE = 10 .^ range((start = 4.179 - 125 * 6.0e-6), step = 6.0e-6, length = 8575 + 125)
end
const SKYCACHE_NPIX = length(logUniWaveAPOGEE)
# only the sky columns are ever read, so the fixture's ar1Dunical needs just a few
const SKYCACHE_NCOL = 20

"""
Write `raw/<tele>/<mjd>/exposures` and `raw/<tele>/<mjd>/fibers/<plate_id>` for a
plate-era night (mjd below the FPS divide, so `plate_id` is the config column and no
FPI-guide overwrite happens). Row i of the exposures table carries `exposure == i`,
which is what `get_fibTargDict` indexes by.
"""
function write_fixture_almanac!(path, teles, mjd, nexp; plate_id = 9999, nsky = 12)
    h5open(path, "w") do f
        for tele in teles
            g = create_group(f, "raw/$(tele)/$(mjd)/exposures")
            g["exposure"] = collect(1:nexp)
            g["image_type"] = fill("object", nexp)
            g["plate_id"] = fill(plate_id, nexp)
            g["config_id"] = fill(plate_id, nexp)
            g["n_read"] = fill(47, nexp)
            g["chip_flags"] = fill(7, nexp)
            g["flagged_bad"] = zeros(Int, nexp)
            gf = create_group(f, "raw/$(tele)/$(mjd)/fibers/$(plate_id)")
            # fiberIndx = 301 - fiberID; make fiber INDICES 1:nsky the sky fibers
            gf["fiber_id"] = [301 - i for i in 1:300]
            gf["category"] = [i <= nsky ? "sky" : "science" for i in 1:300]
            gf["sdss_id"] = collect(1:300)
            gf["fiber_type"] = fill("APOGEE", 300)
        end
    end
    return path
end

"""Write a fixture ar1Dunical file whose columns carry a recognisable `tag`."""
function write_fixture_1d!(reduxBase, tele, mjd, expnum, tag; seed = 1)
    mkpath(joinpath(reduxBase, "apred", mjd))
    p = get_1Duni_name(reduxBase, tele, mjd, expnum)
    rng = Random.MersenneTwister(seed)
    h5open(p, "w") do f
        f["flux_1d"] = 100.0 .* tag .+ randn(rng, SKYCACHE_NPIX, SKYCACHE_NCOL)
        f["ivar_1d"] = fill(1.0, SKYCACHE_NPIX, SKYCACHE_NCOL)
        f["mask_1d"] = collect(trues(SKYCACHE_NPIX, SKYCACHE_NCOL))  # Bool, as AR.jl writes it
    end
    return p
end

function make_fixture(root; teles = ["apo"], mjd = "58588", expnum = 11, nexp = 12,
        tags = Dict("apo" => 1.0))
    reduxBase = joinpath(root, "redux")
    mkpath(reduxBase)
    alm = write_fixture_almanac!(joinpath(root, "almanac.h5"), teles, mjd, nexp)
    for tele in teles
        write_fixture_1d!(reduxBase, tele, mjd, expnum, tags[tele]; seed = hash(tele) % 1000 + 1)
    end
    return (reduxBase = reduxBase, almanacFile = alm, mjd = mjd, expnum = expnum)
end

bundles_equal(a, b) = (a.skyfibIndxs == b.skyfibIndxs && a.skyFibBits == b.skyFibBits &&
                       a.mskSky == b.mskSky && a.skyBit == b.skyBit &&
                       a.nSkyFibers == b.nSkyFibers && a.survspec == b.survspec &&
                       a.survivar == b.survivar && a.survmsk == b.survmsk)

fixture_key(fx, tele) = sky_cache_key(tele, fx.mjd, fx.expnum,
    get_1Duni_name(fx.reduxBase, tele, fx.mjd, fx.expnum), fx.almanacFile; skyZcut = 10)

@testset "sky cache: disabled by default; a hit reproduces a miss exactly" begin
    mktempdir() do root
        fx = make_fixture(root)
        withenv("ARM_SKY_CACHE_DIR" => nothing) do
            @test isnothing(sky_cache_root())
        end
        withenv("ARM_SKY_CACHE_DIR" => "") do
            @test isnothing(sky_cache_root())
        end

        nocache = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = nothing)
        @test nocache.nSkyFibers == 12
        @test size(nocache.survspec) == (SKYCACHE_NPIX, 12)

        croot = joinpath(root, "cache")
        miss = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = croot)
        @test bundles_equal(miss, nocache)
        key = fixture_key(fx, "apo")
        path = sky_cache_path(croot, "apo", fx.mjd, fx.expnum, key)
        @test isfile(path)                      # the miss published an entry
        hit = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = croot)
        @test bundles_equal(hit, nocache)       # the hit is bit-identical to the miss

        # prove the value really comes FROM the cache: poison the entry with a valid
        # but different bundle and check that it is what comes back
        poisoned = merge(nocache, (survspec = nocache.survspec .+ 1.0,))
        @test write_sky_bundle!(path, key, poisoned)
        @test get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile;
            cache_root = croot).survspec == poisoned.survspec
    end
end

@testset "sky cache: COLLISION - same mjd+expnum at both telescopes never aliases" begin
    mktempdir() do root
        # the measured hazard: 31,307 (mjd, expnum) pairs occur at BOTH apo and lco.
        fx = make_fixture(root; teles = ["apo", "lco"], tags = Dict("apo" => 1.0, "lco" => 7.0))
        croot = joinpath(root, "cache")
        b_apo = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = croot)
        b_lco = get_sky_bundle(fx.reduxBase, "lco", fx.mjd, fx.expnum, fx.almanacFile; cache_root = croot)

        k_apo, k_lco = fixture_key(fx, "apo"), fixture_key(fx, "lco")
        @test k_apo != k_lco
        p_apo = sky_cache_path(croot, "apo", fx.mjd, fx.expnum, k_apo)
        p_lco = sky_cache_path(croot, "lco", fx.mjd, fx.expnum, k_lco)
        @test p_apo != p_lco
        @test isfile(p_apo) && isfile(p_lco)
        # telescope is in BOTH the directory and the file name
        @test occursin(Base.Filesystem.path_separator * "apo" * Base.Filesystem.path_separator, p_apo)
        @test occursin(Base.Filesystem.path_separator * "lco" * Base.Filesystem.path_separator, p_lco)
        @test startswith(basename(p_apo), "skybundle_apo_")
        @test startswith(basename(p_lco), "skybundle_lco_")

        # content-level proof: the two telescopes' data differ and each entry serves its own
        @test b_apo.survspec != b_lco.survspec
        @test read_sky_bundle(p_apo, k_apo).survspec == b_apo.survspec
        @test read_sky_bundle(p_lco, k_lco).survspec == b_lco.survspec
        # a cross-read is refused: the key is stored INSIDE the file, not just in its name
        @test isnothing(read_sky_bundle(p_apo, k_lco))
        @test isnothing(read_sky_bundle(p_lco, k_apo))

        # and the end-to-end path never mixes them up
        @test get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile;
            cache_root = croot).survspec == b_apo.survspec
        @test get_sky_bundle(fx.reduxBase, "lco", fx.mjd, fx.expnum, fx.almanacFile;
            cache_root = croot).survspec == b_lco.survspec
    end
end

@testset "sky cache: INVALIDATION - changed inputs are never served stale" begin
    mktempdir() do root
        fx = make_fixture(root)
        croot = joinpath(root, "cache")
        b1 = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = croot)
        k1 = fixture_key(fx, "apo")

        # (a) the ar1Dunical is re-reduced -> different key -> the NEW value is returned.
        #     This is the exact bug class that has bitten this project before.
        sleep(0.05)
        write_fixture_1d!(fx.reduxBase, "apo", fx.mjd, fx.expnum, 5.0; seed = 99)
        k2 = fixture_key(fx, "apo")
        @test k2 != k1
        b2 = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = croot)
        @test b2.survspec != b1.survspec
        @test isfile(sky_cache_path(croot, "apo", fx.mjd, fx.expnum, k1))  # old entry lingers
        @test isfile(sky_cache_path(croot, "apo", fx.mjd, fx.expnum, k2))  # but is never used

        # (b) the guard tunable changes -> different key
        ar1D = get_1Duni_name(fx.reduxBase, "apo", fx.mjd, fx.expnum)
        @test sky_cache_key("apo", fx.mjd, fx.expnum, ar1D, fx.almanacFile; skyZcut = 5) != k2

        # (c) the almanac is rebuilt -> different key, and the new sky set is used
        sleep(0.05)
        write_fixture_almanac!(fx.almanacFile, ["apo"], fx.mjd, 12; nsky = 14)
        k3 = fixture_key(fx, "apo")
        @test k3 != k2
        b3 = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = croot)
        @test b3.nSkyFibers == 14

        # (d) a schema bump moves the whole tree, so no v1 entry can be reached
        @test occursin("v$(SKY_CACHE_SCHEMA)", sky_cache_path(croot, "apo", fx.mjd, fx.expnum, k3))

        # (e) a missing input forces a miss rather than a guess
        @test isnothing(sky_cache_key("apo", fx.mjd, fx.expnum,
            joinpath(root, "nope.h5"), fx.almanacFile; skyZcut = 10))
    end
end

@testset "sky cache: RACE (tasks) - a reader never observes a partial entry" begin
    mktempdir() do root
        fx = make_fixture(root)
        croot = joinpath(root, "cache")
        truth = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = nothing)
        key = fixture_key(fx, "apo")
        path = sky_cache_path(croot, "apo", fx.mjd, fx.expnum, key)

        bad = Threads.Atomic{Int}(0)
        seen = Threads.Atomic{Int}(0)
        stop = Ref(false)
        reader = @async begin
            while !stop[]
                got = read_sky_bundle(path, key)
                if !isnothing(got)
                    Threads.atomic_add!(seen, 1)
                    bundles_equal(got, truth) || Threads.atomic_add!(bad, 1)
                end
                yield()
            end
        end
        @sync for _ in 1:16
            @async write_sky_bundle!(path, key, truth)
        end
        stop[] = true
        wait(reader)
        @test bad[] == 0                    # every observed entry was complete
        @test isfile(path)
        @test bundles_equal(read_sky_bundle(path, key), truth)
        # every publish either renamed into place or cleaned up after itself
        @test readdir(dirname(path)) == [basename(path)]
    end
end

@testset "sky cache: RACE (processes) - concurrent publishes stay consistent" begin
    mktempdir() do root
        fx = make_fixture(root)
        croot = joinpath(root, "cache")
        truth = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = nothing)
        key = fixture_key(fx, "apo")
        path = sky_cache_path(croot, "apo", fx.mjd, fx.expnum, key)

        # Real OS processes and a real rename(2): this is the bulk-run hazard, where
        # workers share nothing but the filesystem.
        procs = addprocs(4, exeflags = ["--project=$(dirname(Base.active_project()))"])
        try
            skysrc = abspath(joinpath(@__DIR__, "..", "src", "skyCache.jl"))
            # remotecall_eval rather than @everywhere: @everywhere expands to a
            # :toplevel expression and cannot be used inside a testset/closure
            Distributed.remotecall_eval(Main, procs, quote
                using HDF5, Random
                logUniWaveAPOGEE = 10 .^ range((start = 4.179 - 125 * 6.0e-6), step = 6.0e-6, length = 8575 + 125)
                include($skysrc)
            end)
            payload = (skyfibIndxs = truth.skyfibIndxs, skyFibBits = truth.skyFibBits,
                mskSky = truth.mskSky, skyBit = truth.skyBit, nSkyFibers = truth.nSkyFibers,
                survspec = truth.survspec, survivar = truth.survivar, survmsk = truth.survmsk)
            out = pmap(Distributed.WorkerPool(procs), 1:8) do _
                ok = write_sky_bundle!(path, key, payload)
                got = read_sky_bundle(path, key)
                (ok, isnothing(got) ? nothing : (got.nSkyFibers, sum(got.survspec), count(got.survmsk)))
            end
            ref = (truth.nSkyFibers, sum(truth.survspec), count(truth.survmsk))
            @test all(first, out)                       # every publish succeeded
            @test all(o -> o[2] == ref, out)            # every read-back was complete
            @test bundles_equal(read_sky_bundle(path, key), truth)
            @test readdir(dirname(path)) == [basename(path)]   # no orphan temp files
        finally
            rmprocs(procs)
        end
    end
end

@testset "sky cache: FAIL-OPEN - corrupt, truncated and unwritable all recompute" begin
    mktempdir() do root
        fx = make_fixture(root)
        croot = joinpath(root, "cache")
        truth = get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = nothing)
        key = fixture_key(fx, "apo")
        path = sky_cache_path(croot, "apo", fx.mjd, fx.expnum, key)
        get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile; cache_root = croot)
        @test isfile(path)

        # (a) garbage in place of an HDF5 file
        write(path, "not an hdf5 file at all")
        @test isnothing(read_sky_bundle(path, key))
        @test bundles_equal(get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum,
            fx.almanacFile; cache_root = croot), truth)

        # (b) valid HDF5 with inconsistent shapes
        h5open(path, "w") do f
            f["schema"] = SKY_CACHE_SCHEMA
            f["key"] = key
            f["skyfibIndxs"] = collect(1:12)
            f["skyFibBits"] = zeros(Int, 12)
            f["mskSky"] = UInt8.(trues(12))
            f["skyBit"] = 0
            f["nSkyFibers"] = 12
            f["survspec"] = zeros(3, 3)
            f["survivar"] = zeros(3, 3)
            f["survmsk"] = UInt8.(trues(3, 3))
        end
        @test isnothing(read_sky_bundle(path, key))
        @test bundles_equal(get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum,
            fx.almanacFile; cache_root = croot), truth)

        # (c) a stale schema or a foreign key is a miss, never a silent reuse
        h5open(path, "w") do f
            f["schema"] = SKY_CACHE_SCHEMA - 1
            f["key"] = key
        end
        @test isnothing(read_sky_bundle(path, key))
        @test isnothing(read_sky_bundle(joinpath(root, "does_not_exist.h5"), key))

        # (d) an unwritable cache root does not crash the run, and stops retrying
        ro = joinpath(root, "readonly")
        mkpath(ro)
        chmod(ro, 0o500)
        try
            SKY_CACHE_WRITE_DISABLED[] = false
            @test bundles_equal(get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum,
                fx.almanacFile; cache_root = ro), truth)
            @test SKY_CACHE_WRITE_DISABLED[]
            @test bundles_equal(get_sky_bundle(fx.reduxBase, "apo", fx.mjd, fx.expnum,
                fx.almanacFile; cache_root = ro), truth)
        finally
            chmod(ro, 0o700)
            SKY_CACHE_WRITE_DISABLED[] = false
        end
    end
end

@testset "sky cache: getSky4visit is unchanged with the cache on" begin
    mktempdir() do root
        fx = make_fixture(root)
        croot = joinpath(root, "cache")
        rng = Random.MersenneTwister(20260907)
        npix = SKYCACHE_NPIX
        x = range(-1, 1, length = npix)
        V_skycont = 50 .* hcat(ones(npix), x, x .^ 2, x .^ 3)
        V_skyline_faint = 10 .* randn(rng, npix, 6)
        skymsk = trues(npix)
        skymsk[1:50] .= false

        off = getSky4visit(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile,
            skymsk, V_skyline_faint, V_skycont; cache_root = nothing)
        on_cold = getSky4visit(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile,
            skymsk, V_skyline_faint, V_skycont; cache_root = croot)
        on_warm = getSky4visit(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile,
            skymsk, V_skyline_faint, V_skycont; cache_root = croot)
        @test length(off) == 7
        for i in 1:length(off)
            @test isequal(off[i], on_cold[i])
            @test isequal(off[i], on_warm[i])
        end
        # the warm call really did hit: a poisoned entry changes the answer
        key = fixture_key(fx, "apo")
        path = sky_cache_path(croot, "apo", fx.mjd, fx.expnum, key)
        b = read_sky_bundle(path, key)
        write_sky_bundle!(path, key, merge(b, (survspec = b.survspec .+ 50.0,)))
        poisoned = getSky4visit(fx.reduxBase, "apo", fx.mjd, fx.expnum, fx.almanacFile,
            skymsk, V_skyline_faint, V_skycont; cache_root = croot)
        @test !isequal(poisoned[2], off[2])
    end
end
