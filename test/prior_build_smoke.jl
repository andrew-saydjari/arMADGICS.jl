# M4 smoke test: run sample_sky.jl's main body (sample_sky_main) on 1 fiber x 1 night
# of real data, for BOTH telescopes. The LCO leg is the regression for the silent-empty
# LCO sky-sampling bug (adjusted fiber index 301-600 passed where the native 1-300
# index was expected -> every getSkyRough call returned nothing).
#
# E5 extension: feed those fresh samples straight into BOTH prior builders
# (build_skyCont, build_skyLines from build_sky_defs.jl) at the smallest viable
# config and assert finite, sanely-shaped prior outputs — the end-to-end
# producer/consumer contract (the .jdat-era builders could not read the sampler's
# HDF5 output at all).
#
# Requires the 2026_05_01 redux products + 2026_04 priors on ceph (read-only);
# SKIPPED (green) when those inputs are unavailable, e.g. on CI.

using Distributed, JLD2, SortFilters, DataFrames

include("../scripts/prior_build/prior_utils.jl")
include("../scripts/prior_build/sample_sky_defs.jl")
isdefined(Main, :gspice) || include("../scripts/prior_build/gspice.jl")
# already included by sky_prior_contract.jl when running under runtests.jl
isdefined(Main, :build_skyCont) || include("../scripts/prior_build/build_sky_defs.jl")

@testset "M4/E5: sample_sky + prior builders smoke (1 fiber x 1 night, apo+lco)" begin
    reduxBase = "/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/"
    almanacFile = joinpath(reduxBase, "almanac", "allobs_57600_61160.h5")
    prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
    prior_dict = Dict{String,String}()
    prior_dict["starCont"] = prior_dir * "2026_04_26/star_priors/APOGEE_starcont_svd_60_f"
    prior_dict["StarContChipGapMsk"] = prior_dir * "2026_04_25/StarContChipGapMsk.h5"

    if !(isfile(almanacFile) && isfile(prior_dict["StarContChipGapMsk"]))
        @info "M4/E5 sample_sky smoke test SKIPPED: real-data inputs not available"
    else
        mjd = "60000"
        sample_runs = Tuple{String,Int,String}[] # (tele, adjfib, sample out_dir) for E5
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

            push!(sample_runs, (tele, adjfib, out_dir))
        end

        # ---- E5: prior builders consume the fresh sampler output directly ----
        chipgap_msk_path = prior_dict["StarContChipGapMsk"]
        @testset "E5: build_skyCont on fresh samples ($(tele))" for (tele, adjfib, out_dir) in sample_runs
            prior_out = joinpath(mktempdir(), "sky_priors")
            # smallest viable config: the SVD dimension is fixed by the pixel count,
            # but only nsub columns are kept — 5 is plenty for a shape/finiteness smoke
            nsub = 5
            fnameC = build_skyCont(adjfib; sample_dir=out_dir,
                chipgap_msk_path=chipgap_msk_path, out_dir=prior_out, nsub=nsub)
            @test isfile(fnameC)
            Vmat = h5read(fnameC, "Vmat")
            λv = h5read(fnameC, "λv")
            @test size(Vmat) == (length(wavetarg), nsub)
            @test all(isfinite, Vmat)
            @test any(!=(0), Vmat)
            @test length(λv) == nsub
            @test all(isfinite, λv)
            @test all(>=(0), λv)
            @test issorted(λv, rev=true) # SVD spectrum comes ordered
            chipgapmsk = read_chipgap_msk(chipgap_msk_path, adjfib)
            @test all(Vmat[.!chipgapmsk, :] .== 0) # prior lives only on unmasked pixels
            println("E5 smoke [$tele]: skyCont prior ", size(Vmat), ", λv[1]=", λv[1])
        end

        @testset "E5: build_skyLines on fresh samples ($(tele))" for (tele, adjfib, out_dir) in sample_runs
            prior_out = joinpath(mktempdir(), "sky_priors")
            # smallest viable config: few SVD columns, a single gspice masking pass,
            # and usamp_factor large enough that ~50 one-night samples are accepted
            # against ~8000 modeled pixels (production: 7). maxbadpix is loosened the
            # same way so single-night sample masks pass, and reg_eps is raised from
            # the production 1e-3: with a rank-~50 covariance over ~7000 pixels the
            # gspice predicted variance sits at the rounding floor and can come out
            # ~-1e-5 (sqrt DomainError); the ridge keeps it positive
            # (predcovar >= eps*|w|^2 in exact arithmetic).
            nsub = 5
            fnameF, fnameFG = build_skyLines(adjfib; sample_dir=out_dir,
                chipgap_msk_path=chipgap_msk_path, out_dir=prior_out,
                nsub_faint=nsub, nsub_bright=nsub, nsigma_schedule=[20],
                usamp_factor=1000, maxbadpix=4000, reg_eps=1.0)
            for (tag, fn) in (("faint", fnameF), ("faintGSPICE", fnameFG))
                @test isfile(fn)
                Vmat = h5read(fn, "Vmat")
                λv = h5read(fn, "λv")
                submsk = Bool.(h5read(fn, "submsk"))
                @test size(Vmat) == (length(wavetarg), nsub)
                @test all(isfinite, Vmat)
                @test any(!=(0), Vmat)
                @test length(λv) == nsub
                @test all(isfinite, λv)
                @test all(>=(0), λv)
                @test length(submsk) == length(wavetarg)
                @test 0 < count(submsk) < length(wavetarg)
                @test all(Vmat[.!submsk, :] .== 0) # prior lives only on the faint submask
                println("E5 smoke [$tele]: skyLines $tag prior ", size(Vmat),
                    ", npix_submsk=", count(submsk), ", λv[1]=", λv[1])
            end
        end
    end
end
