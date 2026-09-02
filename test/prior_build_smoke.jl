# M4 smoke test: run sample_sky.jl's main body (sample_sky_main) on 1 fiber x 1 night
# of real data, for BOTH telescopes. The LCO leg is the regression for the silent-empty
# LCO sky-sampling bug (adjusted fiber index 301-600 passed where the native 1-300
# index was expected -> every getSkyRough call returned nothing).
#
# Requires the 2026_05_01 redux products + 2026_04 priors on ceph (read-only);
# SKIPPED (green) when those inputs are unavailable, e.g. on CI.

using Distributed, JLD2, SortFilters, DataFrames

include("../scripts/prior_build/prior_utils.jl")
include("../scripts/prior_build/sample_sky_defs.jl")

@testset "M4: sample_sky smoke (1 fiber x 1 night, apo+lco)" begin
    reduxBase = "/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/"
    almanacFile = joinpath(reduxBase, "almanac", "allobs_57600_61160.h5")
    prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
    prior_dict = Dict{String,String}()
    prior_dict["starCont"] = prior_dir * "2026_04_26/star_priors/APOGEE_starcont_svd_60_f"
    prior_dict["StarContChipGapMsk"] = prior_dir * "2026_04_25/StarContChipGapMsk.h5"

    if !(isfile(almanacFile) && isfile(prior_dict["StarContChipGapMsk"]))
        @info "M4 sample_sky smoke test SKIPPED: real-data inputs not available"
    else
        mjd = "60000"
        for tele in ("apo", "lco")
            run_lst = get_telemjd_runlist_from_almanac(almanacFile, tele, mjd, accepted_fibtypes=["sky"])
            @test !isempty(run_lst)
            # most-observed sky fiber that night
            adjfib = mode([r.adjfiberindx for r in run_lst])
            if tele == "lco"
                @test adjfib > 300 # LCO runlist carries adjusted (301-600) indices
            end

            out_dir = joinpath(mktempdir(), "sky_prior_disk")
            iterlist = sample_sky_main(reduxBase, almanacFile, adjfib:adjfib;
                prior_dict=prior_dict, out_dir=out_dir,
                tele_mjd_pairs=[(tele, mjd)], loc_parallel=false)
            @test length(iterlist) == 1
            @test iterlist[1][1] == adjfib

            fibtag = lpad(adjfib, 3, "0")
            fluxf = joinpath(out_dir, "skyflux_$(fibtag).h5")
            contf = joinpath(out_dir, "skycont_$(fibtag).h5")
            linesf = joinpath(out_dir, "skyline_$(fibtag).h5")
            @test isfile(fluxf)
            @test isfile(contf)
            @test isfile(linesf)

            skyflux = h5read(fluxf, "skyflux")
            skycont = h5read(contf, "skycont")
            skylines = h5read(linesf, "skyline")
            nsamp = size(skyflux, 2)
            # THE M4/LCO regression: sky sampling must be NONEMPTY (was silently empty for LCO)
            @test nsamp > 0
            @test size(skycont) == (length(wavetarg), nsamp)
            @test size(skylines) == (length(wavetarg), nsamp)
            # the decomposition products, not stale/garbage bindings
            @test eltype(skycont) <: AbstractFloat
            @test eltype(skylines) <: AbstractFloat
            @test all(isfinite, skycont)
            @test all(isfinite, skylines)
            @test any(!=(0), skycont)
            @test any(!=(0), skylines)
            @test skycont != skyflux
            println("M4 smoke [$tele $mjd]: adjfibindx=$adjfib, n_sky_samples=$nsamp")

            # checkpoint-restart path: flux/ivar/msk files exist, cont/lines deleted ->
            # stage 2 must run standalone (count_sky is recovered from disk; it used to
            # be an UndefVarError here)
            rm(contf); rm(linesf)
            get_sky_samples((adjfib, iterlist[1][2]); reduxBase=reduxBase,
                almanacFile=almanacFile, prior_dict=prior_dict, out_dir=out_dir)
            @test isfile(contf)
            @test isfile(linesf)
            @test h5read(contf, "skycont") == skycont
        end
    end
end
