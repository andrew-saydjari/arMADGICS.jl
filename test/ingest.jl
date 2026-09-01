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
