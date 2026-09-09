# E5 unit test: pin the sky-prior producer/consumer contract on synthetic data.
#
# The sampler (sample_sky_defs.jl) writes skyivar_* = INVERSE variance; gspice's
# gspice_covar_iter_mask consumes INVERSE variance (gspice_standard_scale scales each
# spectrum by refscale = sqrt(mean ivar), i.e. into S/N units). The .jdat-era builder
# read a VARIANCE skyvar_* file and inverted it; the E5 wiring
# (gspice_ivar_from_skyivar) passes skyivar through UNINVERTED. This test makes a
# wrong inversion visibly break an ivar-weighted-mean recovery, and shows the
# masked-pixel (skyivar = 0) guard: direct use is finite zero-weight, inversion
# injects Infs.

isdefined(Main, :gspice) || include("../scripts/prior_build/gspice.jl")
include("../scripts/prior_build/build_sky_defs.jl")

@testset "E5: skyivar -> gspice ivar contract (synthetic)" begin
    rng = Random.MersenneTwister(20260903)
    npix, nspec = 20, 400
    truemean = 5.0
    bias = 10.0
    # half quiet, trustworthy spectra (sigma=0.05); half noisy spectra (sigma=20)
    # that also carry a systematic bias — correct ivar weighting must suppress them
    sigma = vcat(fill(0.05, nspec ÷ 2), fill(20.0, nspec ÷ 2))
    offset = vcat(fill(0.0, nspec ÷ 2), fill(bias, nspec ÷ 2))
    flux = truemean .+ offset .+ sigma .* randn(rng, nspec, npix)

    # sampler convention: (nwave, nsamp) matrix of inverse variances
    skyivar_sub = collect(((1 ./ sigma .^ 2) .* ones(nspec, npix))')

    ivar = gspice_ivar_from_skyivar(skyivar_sub)
    @test size(ivar) == (nspec, npix)          # gspice's (nspec, npix) convention
    @test ivar == collect(skyivar_sub')        # passes through UNINVERTED
    @test all(isfinite, ivar)

    # gspice_standard_scale weights each spectrum by refscale^2 = mean ivar; with the
    # correct wiring the quiet spectra dominate, with variance-in-place-of-ivar the
    # noisy spectra get the (hugely) wrong weight
    _, refscale, _ = gspice.gspice_standard_scale(flux, ivar)
    w = vec(refscale) .^ 2
    @test w[1] / w[end] ≈ (sigma[end] / sigma[1])^2  # = 1.6e5: quiet >> noisy
    specmean = vec(sum(flux, dims=2)) ./ npix
    wmean = sum(w .* specmean) / sum(w)
    @test abs(wmean - truemean) < 0.02 # biased/noisy half correctly downweighted

    # THE WRONG WIRING (treating skyivar as variance and inverting, or equivalently
    # treating a variance file as ivar) visibly breaks the weighted-mean recovery
    ivar_wrong = 1 ./ ivar
    _, refscale_w, _ = gspice.gspice_standard_scale(flux, ivar_wrong)
    ww = vec(refscale_w) .^ 2
    @test ww[end] / ww[1] ≈ (sigma[end] / sigma[1])^2  # weights flipped to the noisy half
    wmean_wrong = sum(ww .* specmean) / sum(ww)
    @test abs(wmean_wrong - truemean) > bias / 2 # dragged to the biased half's mean

    # masked pixels: sampler writes skyivar = 0 there; direct use is valid zero-weight
    # gspice input (finite), while inverting would inject Infs
    skyivar0 = copy(skyivar_sub)
    skyivar0[1:3, 1] .= 0.0
    ivar0 = gspice_ivar_from_skyivar(skyivar0)
    @test all(isfinite, ivar0)
    @test count(iszero, ivar0) == 3
    @test any(isinf, 1 ./ skyivar0)  # what the old 1 ./ fluxvar path would have produced

    # guard rails: non-finite or negative "ivar" (e.g. a variance file with Inf at
    # masked pixels fed in by mistake) is rejected loudly
    skyivar_bad = copy(skyivar_sub); skyivar_bad[1, 1] = Inf
    @test_throws ErrorException gspice_ivar_from_skyivar(skyivar_bad)
    skyivar_neg = copy(skyivar_sub); skyivar_neg[1, 1] = -1.0
    @test_throws ErrorException gspice_ivar_from_skyivar(skyivar_neg)
end

@testset "Sky-sample exposure-guard stamp (round trip)" begin
    # The sampler stamps `exposure_guard` on skyflux_NNN.h5 and the builder's
    # read_sky_sample reports it, so a prior built from unguarded samples is
    # visibly labelled as such. This test round-trips a real HDF5 file because
    # the reader is easy to get subtly wrong: HDF5.jl's `attrs(f)[name]` returns
    # the attribute VALUE, and wrapping it in read() throws — which, swallowed
    # by the reader's catch, silently downgrades a correctly-stamped sample dir
    # to "no stamp". That is exactly the class of silent failure this guard
    # exists to prevent, so it is pinned here.
    mktempdir() do dir
        # (a) properly stamped, guarded
        h5open(joinpath(dir, "skyflux_295.h5"), "w") do f
            f["skyflux"] = zeros(3, 2)
            attrs(f)["exposure_guard"] = "EXPFLAG_NO_SCIENCE"
            attrs(f)["exposure_guard_mask"] = EXPFLAG_NO_SCIENCE
            attrs(f)["exposure_guard_almanac"] = "/some/decorated.h5"
            attrs(f)["exposure_guard_unguarded_env"] = "0"
        end
        empty!(_SKY_GUARD_REPORTED)
        @test_logs min_level = Base.CoreLogging.Warn begin
            check_sky_sample_guard_stamp(dir, 295)   # guarded -> no warning
        end

        # (b) stamped but explicitly unguarded -> loud warning
        dir2 = joinpath(dir, "unguarded"); mkpath(dir2)
        h5open(joinpath(dir2, "skyflux_295.h5"), "w") do f
            f["skyflux"] = zeros(3, 2)
            attrs(f)["exposure_guard"] = "EXPFLAG_NO_SCIENCE"
            attrs(f)["exposure_guard_almanac"] = "/some/undecorated.h5"
            attrs(f)["exposure_guard_unguarded_env"] = "1"
        end
        empty!(_SKY_GUARD_REPORTED)
        @test_logs (:warn, r"EXPLICITLY UNGUARDED") match_mode = :any begin
            check_sky_sample_guard_stamp(dir2, 295)
        end

        # (c) no stamp at all (pre-2026-09-08 samples) -> loud warning
        dir3 = joinpath(dir, "unstamped"); mkpath(dir3)
        h5open(joinpath(dir3, "skyflux_295.h5"), "w") do f
            f["skyflux"] = zeros(3, 2)
        end
        empty!(_SKY_GUARD_REPORTED)
        @test_logs (:warn, r"UNGUARDED SKY SAMPLES") match_mode = :any begin
            check_sky_sample_guard_stamp(dir3, 295)
        end

        # reported once per sample dir, not once per fiber
        empty!(_SKY_GUARD_REPORTED)
        check_sky_sample_guard_stamp(dir3, 295)
        @test_logs min_level = Base.CoreLogging.Warn begin
            check_sky_sample_guard_stamp(dir3, 295)
        end
    end
end
