# Unit tests for the E5 sky runtime wiring (audit items 2/4):
# select_sky_fibers extraction, the restored obs-count threshold, and the
# prior-based continuum/line decomposition (sky_decomp).

@testset "sky wiring: select_sky_fibers == combine_sky_fibers guard verdicts" begin
    rng = Random.MersenneTwister(20260905)
    npix, ncand = 2000, 8
    skyspec = 100 .+ randn(rng, npix, ncand)
    skyivar = ones(npix, ncand)
    skymsk = trues(npix, ncand)
    skyspec[:, 3] .= NaN                    # excluded by validate (non-finite)
    skyspec[:, 5] .= 1e5 .+ randn(rng, npix) # excluded by the scale z-cut

    mskSky, nSky, skyBit, bits = select_sky_fibers(skyspec, skyivar, skymsk)
    nSky2, meanL, V, mskloc, skyBit2, mskSky2, bits2 = combine_sky_fibers(skyspec, skyivar, skymsk)
    # the refactor must be behavior-neutral: identical verdicts through both paths
    @test mskSky == mskSky2
    @test nSky == nSky2 == 6
    @test skyBit == skyBit2
    @test bits == bits2
    @test !mskSky[3] && !mskSky[5]
    @test (skyBit & SKY_EXCLUDED_FIBER_BIT) != 0
    @test (skyBit & SKY_ZCUT_FIBER_BIT) != 0
end

@testset "sky wiring: sky_obs_count_mask restores the DR17 obs-count threshold" begin
    npix = 10
    # pixel i has i sky-fiber observations (non-nanzero line residuals)
    outLines = zeros(npix, 8)
    for i in 1:npix, j in 1:min(i, 8)
        outLines[i, j] = 1.0
    end
    m = sky_obs_count_mask(outLines, 5)
    # DR17 semantics: strictly MORE THAN sky_obs_thresh observations required
    @test m == (min.(1:npix, 8) .> 5)
    @test !m[5] && m[6]
    # NaN counts as no observation
    outLines[7, :] .= NaN
    @test !sky_obs_count_mask(outLines, 5)[7]
end

@testset "sky wiring: sky_decomp separates continuum from lines (synthetic)" begin
    rng = Random.MersenneTwister(20260906)
    npix = 400
    x = range(-1, 1, length=npix)
    # smooth continuum basis (poly-like) and sparse line basis
    V_skycont = 50 .* hcat(ones(npix), x, x .^ 2, x .^ 3)
    linepix = 25:25:375
    V_skyline = zeros(npix, length(linepix))
    for (k, p) in enumerate(linepix)
        V_skyline[p-1:p+1, k] .= [30.0, 100.0, 30.0]
    end

    truecont = 20 .+ 10 .* x .^ 2
    truelines = zeros(npix)
    for p in linepix[1:2:end]
        truelines[p-1:p+1] .+= [15.0, 50.0, 15.0]
    end
    noise = 0.1 .* randn(rng, npix)
    spec = truecont .+ truelines .+ noise
    outvar = fill(0.01, npix)
    msk = trues(npix)
    msk[1:10] .= false # some masked pixels

    contvec = sky_decomp(spec, outvar, msk, V_skyline, V_skycont)
    @test length(contvec) == npix
    @test all(isfinite, contvec)
    # continuum recovered away from lines and masked edges
    offline = trues(npix)
    for p in linepix
        offline[p-3:p+3] .= false
    end
    offline[1:14] .= false
    @test maximum(abs.((contvec.-truecont)[offline])) < 1.5
    # the line flux stays out of the continuum component at line pixels
    @test maximum(abs.((contvec.-truecont)[.!offline .& msk])) < 10.0
end
