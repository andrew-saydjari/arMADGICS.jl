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

@testset "M2: failed_pipeline_out placeholder shapes" begin
    npix = 8700
    nstarcoef = 50
    lvllens = [281, 81, 61]
    simplemsk = falses(npix)
    fspec = fill(NaN, npix)
    fivar = fill(NaN, npix)
    ingestBit = 2^1 | 2^2
    out = failed_pipeline_out(simplemsk, NaN, NaN, fspec, fivar, 0, NaN, ingestBit, nstarcoef, lvllens)

    @test length(out) == 5
    # meta block: same accessors multi_spectra_batch uses
    meta = out[1]
    @test meta[1] == 0                   # data_pix_cnt
    @test isnan(meta[2]) && isnan(meta[3])
    @test length(meta[4]) == npix && length(meta[5]) == npix
    @test meta[8] === simplemsk
    @test meta[11] == ingestBit          # new ingestBit column
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
