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
