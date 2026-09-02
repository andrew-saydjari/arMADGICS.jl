@testset "M1: snr metric + Inf-safe nanzero*" begin
    # Inf-safety of the nanzero* family: +/-Inf must be treated like NaN/0 (no data)
    @test nanzeromedian([1.0, Inf, 2.0, 3.0]) == 2.0
    @test nanzeromedian([1.0, -Inf, 2.0, 3.0]) == 2.0
    @test isnan(nanzeromedian([Inf, -Inf, 0.0, NaN]))
    @test nanzeromean([1.0, Inf, 3.0]) == 2.0
    @test nanzerosum([1.0, -Inf, 3.0]) == 4.0
    @test nanzeroiqr([1.0, 2.0, Inf, 3.0]) == nanzeroiqr([1.0, 2.0, 3.0])
    # unchanged existing semantics: NaN and 0 filtered
    @test nanzeromedian([1.0, NaN, 0.0, 3.0]) == 2.0

    # median_snr: known flux/ivar
    npix = 100
    fspec = fill(10.0, npix)
    fivar = fill(4.0, npix)  # sigma = 0.5 -> snr = 20
    msk = trues(npix)
    @test median_snr(fspec, fivar, msk) == 20.0

    # the pre-M1 (inverted) formula would have given median(flux*sigma) = 5.0
    @test nanzeromedian(fspec ./ sqrt.(fivar)) == 5.0
    @test median_snr(fspec, fivar, msk) != nanzeromedian(fspec ./ sqrt.(fivar))

    # ivar = 0 (masked/no-data) pixels: must not pollute the masked snr
    fivar2 = copy(fivar)
    msk2 = trues(npix)
    fivar2[1:50] .= 0.0
    msk2[1:50] .= false
    @test median_snr(fspec, fivar2, msk2) == 20.0
    # even if a zero-ivar pixel sneaks inside the mask, the zero product is filtered
    msk3 = trues(npix)
    @test median_snr(fspec, fivar2, msk3) == 20.0

    # Inf ivar and NaN flux inside the mask are filtered, not propagated
    fspec4 = fill(10.0, npix)
    fivar4 = fill(4.0, npix)
    fivar4[7] = Inf
    fspec4[13] = NaN
    @test median_snr(fspec4, fivar4, trues(npix)) == 20.0

    # fully-masked spectrum -> NaN, not an error
    @test isnan(median_snr(fspec, fivar, falses(npix)))

    # aborted-exposure regression shape: tiny flux with matched ivar must give snr ~ flux*sqrt(ivar), not its inverse
    fspec5 = fill(0.1, npix)
    fivar5 = fill(100.0, npix) # sigma = 0.1 -> snr = 1
    @test median_snr(fspec5, fivar5, trues(npix)) == 1.0
    @test nanzeromedian(fspec5 ./ sqrt.(fivar5)) == 0.01 # the old, wrong answer
end

@testset "M2: validate_exposure ingest guards" begin
    npix = 8700

    # clean spectrum: no flags, mask untouched
    fspec = fill(10.0, npix)
    fivar = fill(4.0, npix)
    fmsk = trues(npix)
    fmsk[1:100] .= false
    msk, starscale0, ingestBit = validate_exposure(fspec, fivar, fmsk)
    @test ingestBit == 0
    @test !ingest_fatal(ingestBit)
    @test msk == fmsk
    @test starscale0 == 10.0

    # all-masked spectrum: fatal (too few good pixels), completes without throwing
    msk, starscale0, ingestBit = validate_exposure(fspec, fivar, falses(npix))
    @test (ingestBit & 2^2) != 0
    @test ingest_fatal(ingestBit)
    @test count(msk) == 0

    # A4 fixture shape: NaN flux presented with an all-good mask
    fspec_nan = fill(NaN, npix)
    msk, starscale0, ingestBit = validate_exposure(fspec_nan, fivar, trues(npix))
    @test (ingestBit & 2^1) != 0 # flux all NaN/zero
    @test (ingestBit & 2^3) != 0 # non-finite flux inside good mask
    @test (ingestBit & 2^2) != 0 # nothing left after cleaning
    @test (ingestBit & 2^6) != 0 # starscale undefined
    @test ingest_fatal(ingestBit)
    @test count(msk) == 0
    @test isnan(starscale0)

    # Inf/zero/negative ivar pixels inside the good mask: masked + flagged, not fatal
    fivar_bad = fill(4.0, npix)
    fivar_bad[11] = Inf
    fivar_bad[12] = 0.0
    fivar_bad[13] = -1.0
    fivar_bad[14] = NaN
    msk, starscale0, ingestBit = validate_exposure(fspec, fivar_bad, trues(npix))
    @test (ingestBit & 2^4) != 0
    @test !ingest_fatal(ingestBit)
    @test all(.!msk[11:14])
    @test count(msk) == npix - 4

    # tiny-but-nonzero ivar (M3): masked + flagged, not fatal
    fivar_tiny = fill(4.0, npix)
    fivar_tiny[21:25] .= 1e-19 # AR good_pix allows ivar > 1e-20
    msk, starscale0, ingestBit = validate_exposure(fspec, fivar_tiny, trues(npix))
    @test (ingestBit & 2^5) != 0
    @test !ingest_fatal(ingestBit)
    @test all(.!msk[21:25])
    @test count(msk) == npix - 5

    # non-positive median flux (dead/garbage fiber): fatal starscale flag
    fspec_neg = fill(-17.0, npix)
    msk, starscale0, ingestBit = validate_exposure(fspec_neg, fivar, trues(npix))
    @test (ingestBit & 2^6) != 0
    @test ingest_fatal(ingestBit)
    @test starscale0 == -17.0

    # too few good pixels even though each is individually fine
    fmsk_few = falses(npix)
    fmsk_few[1:(INGEST_MIN_GOODPIX-1)] .= true
    msk, starscale0, ingestBit = validate_exposure(fspec, fivar, fmsk_few)
    @test (ingestBit & 2^2) != 0
    @test ingest_fatal(ingestBit)
end

@testset "M3: apVisit continuum floor" begin
    npix = 200
    fspec = fill(100.0, npix)
    skyLines = fill(5.0, npix)
    skyCont = fill(15.0, npix)
    starCont = fill(80.0, npix)
    finalmsk = trues(npix)
    finalmsk[1:10] .= false
    cont_scale = 80.0

    apV = apVisit_from_components(fspec, skyLines, skyCont, starCont, finalmsk, cont_scale)
    @test length(apV) == npix
    @test all(isnan.(apV[1:10]))                      # outside finalmsk
    @test all(apV[11:end] .== (100.0 - 20.0) / 80.0)  # (flux - sky)/cont

    # near-zero and negative-tiny continuum pixels must come back NaN, not +/-huge
    starCont_bad = copy(starCont)
    starCont_bad[20] = 1e-8
    starCont_bad[21] = -1e-8
    starCont_bad[22] = 0.0
    starCont_bad[23] = NaN
    starCont_bad[24] = Inf
    apV = apVisit_from_components(fspec, skyLines, skyCont, starCont_bad, finalmsk, cont_scale)
    @test all(isnan.(apV[20:24]))
    @test all(isfinite.(apV[25:end]))
    @test apV[30] == 1.0

    # pixels below the relative floor (1e-3 * |cont_scale|) are masked
    starCont_faint = copy(starCont)
    starCont_faint[40] = 1e-3 * cont_scale * 0.5 # below floor
    starCont_faint[41] = 1e-3 * cont_scale * 2.0 # above floor
    apV = apVisit_from_components(fspec, skyLines, skyCont, starCont_faint, finalmsk, cont_scale)
    @test isnan(apV[40])
    @test isfinite(apV[41])

    # non-finite cont_scale (e.g. starscale1 = NaN) masks everything rather than erroring
    apV = apVisit_from_components(fspec, skyLines, skyCont, starCont, finalmsk, NaN)
    @test all(isnan.(apV))
end

@testset "M2: failed_pipeline_out placeholder shapes" begin
    npix = 8700
    nstarcoef = 50
    lvllens = [281, 81, 61]
    simplemsk = falses(npix)
    fspec = fill(NaN, npix)
    fivar = fill(NaN, npix)
    ingestBit = 2^1 | 2^2
    skyBit = 2^0
    out = failed_pipeline_out(simplemsk, NaN, NaN, fspec, fivar, 0, NaN, ingestBit, skyBit, nstarcoef, lvllens)

    @test length(out) == 5
    # meta block: same accessors multi_spectra_batch uses
    meta = out[1]
    @test meta[1] == 0                   # data_pix_cnt
    @test isnan(meta[2]) && isnan(meta[3])
    @test length(meta[4]) == npix && length(meta[5]) == npix
    @test meta[8] === simplemsk
    @test meta[11] == ingestBit          # ingestBit column
    @test meta[12] == skyBit             # M-SKY skyBit column
    # RV block
    @test isnan(Float64.(out[2][1][1]))  # RV_pixoff_final
    @test out[2][1][6] == INGEST_FAIL_RV_FLAG
    @test isnan(out[2][1][7])            # RV_pix_var
    for (i, n) in enumerate(lvllens)
        @test length(out[2][2][i][3]) == n # RV_p5delchi2_lvl*
    end
    # chi2 block
    @test out[3][3] == 0                 # final_pix_cnt
    # component block: 11 entries with the shapes the extractor expects
    @test length(out[4]) == 11
    for i in [1, 2, 3, 4, 5, 6, 9, 11]
        @test length(out[4][i]) == npix
    end
    @test length(out[4][7]) == nstarcoef # starLines coefficients
    @test out[4][8] isa Float64          # tot_p5chi2 scalar
    @test out[4][10] == falses(npix)     # final mask
    @test length(out[5]) == npix         # dflux_starlines
end

@testset "M-SKY: validate_sky_fiber" begin
    npix = 8700
    flux = fill(100.0, npix)
    flux[1:200] .= 0.0 # chip-gap-like zeros are fine
    ivar = fill(1.0, npix)
    msk = trues(npix)
    msk[1:200] .= false

    # healthy sky fiber: usable
    @test validate_sky_fiber(flux, ivar, msk) == 0

    # a single non-finite pixel anywhere is enough to poison VLocSkyLines -> excluded
    flux_nan = copy(flux); flux_nan[4000] = NaN
    @test (validate_sky_fiber(flux_nan, ivar, msk) & SKYFIB_NONFINITE_BIT) != 0
    flux_inf = copy(flux); flux_inf[4000] = Inf
    @test (validate_sky_fiber(flux_inf, ivar, msk) & SKYFIB_NONFINITE_BIT) != 0
    # even at a pixel that is masked-bad for the fiber itself (the prior is built unmasked)
    flux_nan2 = copy(flux); flux_nan2[10] = NaN # inside the msk=false region
    @test (validate_sky_fiber(flux_nan2, ivar, msk) & SKYFIB_NONFINITE_BIT) != 0

    # the A4 upstream shape: all-NaN flux presented with an all-good mask
    bit = validate_sky_fiber(fill(NaN, npix), fill(NaN, npix), trues(npix))
    @test (bit & SKYFIB_NONFINITE_BIT) != 0
    @test (bit & SKYFIB_ALLNANZERO_BIT) != 0
    @test (bit & SKYFIB_LOWGOODPIX_BIT) != 0
    @test (bit & SKYFIB_EXCLUDE_BITS) != 0

    # all-zero flux (finite, but no data)
    bit = validate_sky_fiber(zeros(npix), ivar, msk)
    @test (bit & SKYFIB_ALLNANZERO_BIT) != 0
    @test (bit & SKYFIB_NONFINITE_BIT) == 0
    @test (bit & SKYFIB_EXCLUDE_BITS) != 0

    # too few good pixels (mask, bad ivar both shrink the count)
    msk_few = falses(npix); msk_few[1:(SKY_MIN_GOODPIX-1)] .= true
    @test (validate_sky_fiber(flux, ivar, msk_few) & SKYFIB_LOWGOODPIX_BIT) != 0
    ivar_bad = zeros(npix)
    @test (validate_sky_fiber(flux, ivar_bad, msk) & SKYFIB_LOWGOODPIX_BIT) != 0

    # non-positive (finite) median flux: recorded, but NOT an exclusion bit
    # (census: ~3% of exposures carry finite faint sky fibers with slightly negative
    # medians -- harmless, kept, AKS to decide any future exclusion)
    bit = validate_sky_fiber(fill(-5.0, npix), ivar, msk)
    @test (bit & SKYFIB_NEGSCALE_BIT) != 0
    @test (bit & SKYFIB_EXCLUDE_BITS) == 0
end

@testset "M-SKY: combine_sky_fibers" begin
    npix = 2000
    nsky = 8
    rng = Random.MersenneTwister(1234)
    base = 50.0 .+ 10.0 .* rand(rng, npix)
    skyspec = base .* (0.8 .+ 0.4 .* rand(rng, 1, nsky)) .+ randn(rng, npix, nsky)
    skyspec[1:20, :] .= 0.0 # chip-gap-like pixels: zero in every fiber
    skyivar = ones(npix, nsky)
    skymsk = trues(npix, nsky)

    # healthy: outputs are bit-identical to the pre-guard formulas
    nSky, meanL, V, mskloc, skyBit, mskSky, bits = combine_sky_fibers(skyspec, skyivar, skymsk)
    @test nSky == nsky
    @test skyBit == 0
    @test all(mskSky)
    @test all(bits .== 0)
    mean_ref = dropdims(nanzeromean(skyspec, 2), dims=2)
    finrows = isfinite.(mean_ref)
    @test mskloc == finrows
    @test all(.!finrows[1:20]) # all-zero pixels carry no sky info
    @test meanL[finrows] == mean_ref[finrows]
    @test V[finrows, :] == ((skyspec .- mean_ref) ./ sqrt(nsky))[finrows, :]
    # and the guard zeroes what the pre-guard code left as NaN
    @test all(meanL[.!finrows] .== 0)
    @test all(V[.!finrows, :] .== 0)
    @test all(isfinite, V) && all(isfinite, meanL)

    # one partially-NaN sky fiber: excluded; survivors match the healthy-subset construction
    skyspec_p = copy(skyspec)
    skyspec_p[100:2:end, 3] .= NaN
    nSky, meanL, V, mskloc, skyBit, mskSky, bits = combine_sky_fibers(skyspec_p, skyivar, skymsk)
    @test nSky == nsky - 1
    @test (skyBit & SKY_EXCLUDED_FIBER_BIT) != 0
    @test !mskSky[3] && (bits[3] & SKYFIB_NONFINITE_BIT) != 0
    @test all(isfinite, V) && all(isfinite, meanL)
    keep = setdiff(1:nsky, 3)
    nSky2, meanL2, V2, mskloc2, skyBit2, _, _ = combine_sky_fibers(skyspec[:, keep], skyivar[:, keep], skymsk[:, keep])
    @test skyBit2 == 0
    @test meanL == meanL2 && V == V2 && mskloc == mskloc2

    # all-NaN (A4) sky fiber: excluded and flagged (pre-guard it was only dropped by
    # the accident that its NaN scale fails the z-cut comparison)
    skyspec_a = copy(skyspec)
    skyspec_a[:, 5] .= NaN
    nSky, meanL, V, mskloc, skyBit, mskSky, bits = combine_sky_fibers(skyspec_a, skyivar, skymsk)
    @test nSky == nsky - 1
    @test (skyBit & SKY_EXCLUDED_FIBER_BIT) != 0
    @test (bits[5] & SKYFIB_ALLNANZERO_BIT) != 0
    @test all(isfinite, V)

    # scale z-cut outlier: excluded via the pre-existing cut, recorded
    skyspec_z = copy(skyspec)
    skyspec_z[:, 2] .*= 1e6
    nSky, meanL, V, mskloc, skyBit, mskSky, bits = combine_sky_fibers(skyspec_z, skyivar, skymsk)
    @test nSky == nsky - 1
    @test (skyBit & SKY_ZCUT_FIBER_BIT) != 0
    @test (bits[2] & SKYFIB_ZCUT_BIT) != 0
    @test (skyBit & SKY_EXCLUDED_FIBER_BIT) == 0

    # negative-median (finite) sky fiber within the z-cut: KEPT, recorded in skyBit
    scales_neg = [-5.0, 10.0, 25.0, 40.0, 55.0, 70.0, 85.0, 100.0]
    skyspec_n = scales_neg' .+ 0.1 .* randn(rng, npix, nsky)
    nSky, meanL, V, mskloc, skyBit, mskSky, bits = combine_sky_fibers(skyspec_n, skyivar, skymsk)
    @test nSky == nsky # nothing excluded
    @test all(mskSky)
    @test (bits[1] & SKYFIB_NEGSCALE_BIT) != 0
    @test (skyBit & SKY_NEGSCALE_FIBER_BIT) != 0
    @test (skyBit & SKY_EXCLUDED_FIBER_BIT) == 0

    # too few surviving fibers: flagged no-op sky component (zeros), never NaN
    nSky, meanL, V, mskloc, skyBit, mskSky, bits = combine_sky_fibers(skyspec[:, 1:2], skyivar[:, 1:2], skymsk[:, 1:2])
    @test nSky == 2
    @test (skyBit & SKY_TOO_FEW_FIBERS_BIT) != 0
    @test size(V) == (npix, 1)
    @test all(V .== 0) && all(meanL .== 0)
    @test all(mskloc)

    # every candidate bad: fallback, both flags
    allbad = fill(NaN, npix, 3)
    nSky, meanL, V, mskloc, skyBit, mskSky, bits = combine_sky_fibers(allbad, ones(npix, 3), trues(npix, 3))
    @test nSky == 0
    @test (skyBit & SKY_TOO_FEW_FIBERS_BIT) != 0
    @test (skyBit & SKY_EXCLUDED_FIBER_BIT) != 0
    @test all(V .== 0) && all(meanL .== 0)

    # degenerate spread (identical scales -> IQR = 0): keep all fibers instead of the
    # pre-guard behavior (0/0 z-scores silently discarded every fiber)
    ident = repeat(base, 1, 4)
    nSky, meanL, V, mskloc, skyBit, mskSky, bits = combine_sky_fibers(ident, ones(npix, 4), trues(npix, 4))
    @test nSky == 4
    @test skyBit == 0
end
