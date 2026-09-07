## Ingest Module
# minor rewrite on the horizon

using AstroTime
import ApogeeReduction: get_fibTargDict, fiberID2fiberIndx, read_almanac_exp_df

## M-SKY sky-prior guards
# The exposure-level sky-line prior (VLocSkyLines) is an empirical covariance built
# directly from the sky-fiber spectra of the exposure. A single NaN/Inf pixel in any
# included sky fiber poisons the Woodbury solves for EVERY fiber on the exposure
# (all RVs on the exposure degrade/turn to garbage), so sky fibers are validated
# before they may enter the construction, mirroring the M2 ingest philosophy.

# Minimum good pixels for a sky fiber to enter the sky-line prior (matches INGEST_MIN_GOODPIX)
const SKY_MIN_GOODPIX = 500
# Minimum surviving sky fibers required to build the empirical sky-line prior. Below this,
# the sky-line component is SKIPPED (zeros, exactly a no-op in the solve) and flagged via
# skyBit — we do not invent a sky model. Threshold flagged for AKS review: with < 3 fibers
# the mean-subtracted empirical covariance is essentially meaningless (1 fiber -> V == 0,
# 2 fibers -> a single +/- pair).
const SKY_MIN_FIBERS = 3

# per-sky-fiber flag bits (validate_sky_fiber; logged when the guard fires):
#   2^0 non-finite flux anywhere in the spectrum (any such pixel can poison VLocSkyLines:
#       the raw spectrum enters the prior unmasked; healthy ar1Dunical fibers are finite
#       everywhere, with 0 at chip gaps, so this is not over-exclusive)   [EXCLUDES]
#   2^1 flux entirely NaN/zero (the A4 upstream NaN/good-mask shape)      [EXCLUDES]
#   2^2 fewer than SKY_MIN_GOODPIX pixels with good mask, finite flux,
#       finite ivar > 0                                                   [EXCLUDES]
#   2^3 scale = nanzeromedian(flux) non-positive. RECORDED ONLY, does NOT exclude:
#       the 2026_05_01 census found ~3% of object exposures (all APO) carry faint sky
#       fibers with slightly negative medians (background-subtraction noise, full
#       good-pixel counts) — finite, so not poisonous; excluding them would change
#       results on those healthy exposures without a poisoning justification. Whether
#       to exclude them is an AKS decision (they remain subject to the z-cut).
#   2^4 excluded by the pre-existing scale z-cut (recorded by combine_sky_fibers)
const SKYFIB_NONFINITE_BIT = 2^0
const SKYFIB_ALLNANZERO_BIT = 2^1
const SKYFIB_LOWGOODPIX_BIT = 2^2
const SKYFIB_NEGSCALE_BIT = 2^3
const SKYFIB_ZCUT_BIT = 2^4
# bits that remove a fiber from the sky construction (2^3 is informational only)
const SKYFIB_EXCLUDE_BITS = SKYFIB_NONFINITE_BIT | SKYFIB_ALLNANZERO_BIT | SKYFIB_LOWGOODPIX_BIT

# exposure-level skyBit codes (recorded per spectrum in the batch output):
#   2^0 >= 1 candidate sky fiber excluded by validate_sky_fiber (poison guard fired)
#   2^1 >= 1 validated sky fiber excluded by the scale z-cut (pre-existing cut, now recorded)
#   2^2 surviving sky-fiber count < SKY_MIN_FIBERS -> sky-line component skipped (zeros)
#   2^3 no sky fibers in the fiber configuration at all -> sky-line component skipped
#   2^4 >= 1 used sky fiber has non-positive median scale (recorded only, fiber kept)
const SKY_EXCLUDED_FIBER_BIT = 2^0
const SKY_ZCUT_FIBER_BIT = 2^1
const SKY_TOO_FEW_FIBERS_BIT = 2^2
const SKY_NO_FIBERS_BIT = 2^3
const SKY_NEGSCALE_FIBER_BIT = 2^4

"""
    validate_sky_fiber(flux, ivar, msk; min_goodpix=SKY_MIN_GOODPIX)

Per-sky-fiber sanity checks before the fiber may enter the exposure-level sky-line
prior. Returns a bit flag (0 = fiber usable); bit codes documented above.
"""
function validate_sky_fiber(flux, ivar, msk; min_goodpix::Int=SKY_MIN_GOODPIX)
    bit = 0
    if any(!isfinite, flux)
        bit |= SKYFIB_NONFINITE_BIT
    end
    if all(isnanorzero, flux)
        bit |= SKYFIB_ALLNANZERO_BIT
    end
    good = Bool.(msk) .& isfinite.(flux) .& isfinite.(ivar) .& (ivar .> 0)
    if count(good) < min_goodpix
        bit |= SKYFIB_LOWGOODPIX_BIT
    end
    scale = nanzeromedian(flux)
    # non-finite scale only happens when the flux is entirely NaN/zero/Inf, which the
    # exclusion bits above already catch; a finite non-positive scale is recorded only
    if isfinite(scale) && (scale <= 0)
        bit |= SKYFIB_NEGSCALE_BIT
    end
    return bit
end

"""
    combine_sky_fibers(skyspec, skyivar, skymsk; skyZcut=10, min_fibers=SKY_MIN_FIBERS)

Build the exposure-level local sky-line prior from candidate sky-fiber spectra
(columns), excluding fibers whose `validate_sky_fiber` flag intersects
`SKYFIB_EXCLUDE_BITS` and fibers cut by the pre-existing scale z-cut. Returns
`(nSkyFibers, meanLocSkyLines, VLocSkyLines, msk_local_skyLines, skyBit, mskSky, skyFibBits)`.

Guarantees post-guard: `VLocSkyLines` and `meanLocSkyLines` are finite everywhere.
Pixels where every surviving sky fiber is NaN/zero (chip gaps, globally masked pixels)
would give a NaN mean/prior column; they are zeroed AND removed via
`msk_local_skyLines` (these pixels carry no sky information). On a healthy exposure
those pixels are already excluded from every solve by the chip-gap/fiber masks, so
the healthy-path outputs are bit-identical to the pre-guard code.

If fewer than `min_fibers` fibers survive, the sky-line component is skipped:
mean = zeros and V = a single zero column (exactly a no-op in the Woodbury solves),
flagged via `skyBit` — no sky model is invented.
"""
function combine_sky_fibers(skyspec, skyivar, skymsk; skyZcut=10, min_fibers::Int=SKY_MIN_FIBERS)
    npix, ncand = size(skyspec)
    skyFibBits = [validate_sky_fiber(view(skyspec, :, j), view(skyivar, :, j), view(skymsk, :, j)) for j in 1:ncand]
    mskValid = ((skyFibBits .& SKYFIB_EXCLUDE_BITS) .== 0)

    # scale z-cut (pre-existing), computed over validated fibers only (identical to the
    # old cut when all fibers validate: nanzero* already filtered NaN/Inf/0 scales)
    skyScale = dropdims(nanzeromedian(skyspec, 1), dims=1)
    skyMed = nanzeromedian(skyScale[mskValid])
    skyIQR = nanzeroiqr(skyScale[mskValid])
    mskZ = if isfinite(skyIQR) && (skyIQR > 0)
        skyZ = (skyScale .- skyMed) ./ skyIQR
        (abs.(skyZ) .< skyZcut)
    else
        # degenerate spread (identical/too-few scales): no outlier information -> keep all
        trues(ncand)
    end
    for j in findall(mskValid .& .!mskZ)
        skyFibBits[j] |= SKYFIB_ZCUT_BIT
    end
    mskSky = mskValid .& mskZ
    nSkyFibers = count(mskSky)

    skyBit = 0
    if any(.!mskValid)
        skyBit |= SKY_EXCLUDED_FIBER_BIT
    end
    if any(mskValid .& .!mskZ)
        skyBit |= SKY_ZCUT_FIBER_BIT
    end
    if any(j -> mskSky[j] && ((skyFibBits[j] & SKYFIB_NEGSCALE_BIT) != 0), 1:ncand)
        skyBit |= SKY_NEGSCALE_FIBER_BIT
    end

    if nSkyFibers < min_fibers
        skyBit |= SKY_TOO_FEW_FIBERS_BIT
        return nSkyFibers, zeros(npix), zeros(npix, 1), ones(Bool, npix), skyBit, mskSky, skyFibBits
    end

    meanLocSkyLines = dropdims(nanzeromean(skyspec[:, mskSky], 2), dims=2)
    VLocSkyLines = (skyspec[:, mskSky] .- meanLocSkyLines) ./ sqrt(nSkyFibers)
    # pixels with no sky information from any surviving fiber (all NaN/zero -> NaN mean)
    msk_local_skyLines = isfinite.(meanLocSkyLines)
    meanLocSkyLines[.!msk_local_skyLines] .= 0
    VLocSkyLines[.!msk_local_skyLines, :] .= 0
    return nSkyFibers, meanLocSkyLines, VLocSkyLines, msk_local_skyLines, skyBit, mskSky, skyFibBits
end

function getSkyRough(reduxBase, tele, mjd, expnum, almanacFile; skyZcut=10, sky_obs_thresh=5, fibindxoi=nothing)
    # hacks
    f = h5open(almanacFile)
    fibtargDict, fiber_sdss_id_Dict = get_fibTargDict(f, tele, mjd, expnum)
    close(f)

    fibtypelist = map(x -> fibtargDict[x], 1:300)
    skyfibIndxs = findall(map(x->x[1:3] == "sky", fibtypelist)) # allows for skyB fibers

    if !isnothing(fibindxoi)
        if (length(skyfibIndxs) == 0)
            return nothing
        elseif !(fibindxoi in skyfibIndxs)
            return nothing
        end
    end

    if length(skyfibIndxs) == 0
        # M-SKY: was NaN*ones for mean/V, which NaN-poisoned every fiber's solve on the
        # exposure; now a flagged no-op sky component (see combine_sky_fibers docstring)
        npix = length(logUniWaveAPOGEE)
        return 0, zeros(npix), zeros(npix), zeros(npix, 1), ones(Bool, npix), SKY_NO_FIBERS_BIT | SKY_TOO_FEW_FIBERS_BIT
    end

    #get ar1Dname
    ar1Dfname = get_1Duni_name(reduxBase, tele, mjd, expnum)

    # could speed up by only reading the columns we need
    f = jldopen(ar1Dfname)
    skyspec = f["flux_1d"][:, skyfibIndxs]
    skyivar = f["ivar_1d"][:, skyfibIndxs]
    skymsk = f["mask_1d"][:, skyfibIndxs]
    close(f)

    nSkyFibers, meanLocSkyLines, VLocSkyLines, msk_local_skyLines, skyBit, mskSky, skyFibBits =
        combine_sky_fibers(skyspec, skyivar, skymsk; skyZcut=skyZcut)
    if skyBit != 0
        flagged = [(skyfibIndxs[j], skyFibBits[j]) for j in findall(skyFibBits .!= 0)]
        println("getSkyRough: sky guard flagged tele=$tele, mjd=$mjd, expnum=$expnum: skyBit=$skyBit, (fiberindx, skyFibBit)=$flagged")
        flush(stdout)
    end

    if !isnothing(fibindxoi)
        if !(fibindxoi in skyfibIndxs[mskSky])
            return nothing
        else
            localIndx = findfirst(skyfibIndxs.==fibindxoi)
            return skyspec[:,localIndx], skyivar[:,localIndx], skymsk[:,localIndx]
        end
    end

    # msk_local_skyLines = dropdims(sum(.!isnanorzero.(skyspec[:, mskSky]), dims=2), dims=2) .> sky_obs_thresh
    meanLocSky = zero(meanLocSkyLines) # hack and ignores VLocSky
    return nSkyFibers, meanLocSky, meanLocSkyLines, VLocSkyLines, msk_local_skyLines, skyBit
end

## M2/M3 ingest guards
# Minimum number of good pixels (out of 8700) required to attempt the MADGICS solve
const INGEST_MIN_GOODPIX = 500
# Pixels with ivar below this fraction of the median good-pixel ivar are masked (M3):
# AR's good_pix only requires ivar > 1e-20, so tiny-but-nonzero ivar pixels otherwise
# enter Diagonal(1 ./ fivar) with ~1e20 variances and dominate the deblend.
const INGEST_TINY_IVAR_RELFAC = 1e-6

# ingestBit failure/flag codes (recorded per spectrum in the batch output):
#   2^0 runtime error caught during the solve (spectrum skipped)
#   2^1 flux is entirely NaN/zero                              (fatal -> skip)
#   2^2 fewer than INGEST_MIN_GOODPIX good pixels after checks (fatal -> skip)
#   2^3 non-finite flux inside the good mask (pixels masked)
#   2^4 non-finite or non-positive ivar inside the good mask (pixels masked)
#   2^5 tiny-ivar pixels masked (below INGEST_TINY_IVAR_RELFAC * median good ivar)
#   2^6 starscale0 = nanzeromedian(flux) non-finite or <= 0    (fatal -> skip)
const INGEST_RUNTIME_ERROR_BIT = 2^0
const INGEST_FATAL_BITS = 2^0 | 2^1 | 2^2 | 2^6

ingest_fatal(ingestBit::Int) = (ingestBit & INGEST_FATAL_BITS) != 0

"""
    validate_exposure(fspec, fivar, fmsk; min_goodpix, tiny_ivar_relfac)

Sanity-check a 1D uni-cal spectrum before it enters the MADGICS solve (M2/M3).
Returns `(msk, starscale0, ingestBit)` where `msk` is the good-pixel mask with
non-finite flux, non-finite/non-positive ivar, and tiny-ivar pixels removed,
`starscale0 = nanzeromedian(fspec)`, and `ingestBit` encodes the flag bits above.
Callers should skip the spectrum (but not crash) when `ingest_fatal(ingestBit)`.
"""
function validate_exposure(fspec, fivar, fmsk;
        min_goodpix::Int = INGEST_MIN_GOODPIX,
        tiny_ivar_relfac::Float64 = INGEST_TINY_IVAR_RELFAC)
    ingestBit = 0
    msk = Bool.(fmsk)

    if all(isnanorzero, fspec)
        ingestBit |= 2^1 # flux is literally all NaNs/zeros (cf. A4 upstream bug)
    end

    badflux = msk .& .!isfinite.(fspec)
    if any(badflux)
        ingestBit |= 2^3
        msk = msk .& .!badflux
    end

    badivar = msk .& (.!isfinite.(fivar) .| (fivar .<= 0))
    if any(badivar)
        ingestBit |= 2^4
        msk = msk .& .!badivar
    end

    if any(msk)
        medivar = median(fivar[msk])
        tinyivar = msk .& (fivar .< tiny_ivar_relfac * medivar)
        if any(tinyivar)
            ingestBit |= 2^5
            msk = msk .& .!tinyivar
        end
    end

    starscale0 = nanzeromedian(fspec)
    if !(isfinite(starscale0) && (starscale0 > 0))
        ingestBit |= 2^6
    end

    if count(msk) < min_goodpix
        ingestBit |= 2^2
    end

    return msk, starscale0, ingestBit
end

function getExposure(reduxBase, tele, mjd, expnum, adjfiberindx)
    fiberindx = adjfiberindx2fiberindx(adjfiberindx)
    ar1Dfname = get_1Duni_name(reduxBase, tele, mjd, expnum)
    f = h5open(ar1Dfname)
    fspec = f["flux_1d"][:, fiberindx]
    fivar = f["ivar_1d"][:, fiberindx]
    fmsk = f["mask_1d"][:, fiberindx]
    close(f)
    metaexport = []
    return fspec, fivar, fmsk, metaexport
end

function get_telemjd_runlist_from_almanac(
        almanacFile, tele, mjd;
        accepted_fibtypes::Vector{String} = ["sci", "tel"]
    )
    teleind = (tele[1:3] == "lco") ? 2 : 1
    f = h5open(almanacFile)
    df_exp = read_almanac_exp_df(f, tele, mjd)
    msk_obj = (df_exp.image_type .== "object")
    msk_obj .&= (df_exp.n_read .> 3) .& (df_exp.chip_flags .== 7) .& (df_exp.flagged_bad .== 0)
    row_exp = df_exp[msk_obj, :].exposure
    run_lsts = []
    for expnum in row_exp
        fibtargDict, fiber_sdss_id_Dict = get_fibTargDict(f, tele, mjd, expnum)
        fibtypelist = map(x -> fibtargDict[x], 1:300)
        fiber_sdss_id_list = map(x -> fiber_sdss_id_Dict[x], 1:300)
        msk_fiberok = map(x -> x[1:3] in accepted_fibtypes, fibtypelist)
        targfibIndxs = findall(msk_fiberok)
        adjtargfibIndxs = targfibIndxs .+ (teleind - 1) * 300

        iterexp = Iterators.zip(
            Iterators.repeated(tele), Iterators.repeated(mjd),
            Iterators.repeated(expnum), adjtargfibIndxs, fiber_sdss_id_list[targfibIndxs]
        )
        iterexp_named = map(
            ((t, m, e, f, s),) -> (tele=t, mjd=m, expnum=e, adjfiberindx=f, sdss_id=s),
            iterexp
        )
        push!(run_lsts, collect(iterexp_named))
    end
    run_lst = vcat(run_lsts...)
    close(f)
    return run_lst
end

function sky_decomp(outvec, outvar, simplemsk, V_skyline_bright, V_skyline_faint, V_skycont)
    ## Select data for use (might want to handle mean more generally)
    Xd_obs = outvec[simplemsk]

    ## Set up residuals prior
    A = Diagonal(outvar[simplemsk])
    Ainv = Diagonal(1 ./ outvar[simplemsk])

    ## Set up priors
    V_skyline_bright_c = V_skyline_bright
    V_skyline_bright_r = V_skyline_bright_c[simplemsk, :]
    V_skyline_faint_c = V_skyline_faint
    V_skyline_faint_r = V_skyline_faint_c[simplemsk, :]
    V_skycont_c = V_skycont
    V_skycont_r = V_skycont_c[simplemsk, :]

    # Compute sky line/continuum separation
    Vcomb = hcat(V_skyline_bright_r, V_skyline_faint_r, V_skycont_r)
    Ctotinv = LowRankMultMat([Ainv, Vcomb], wood_precomp_mult, wood_fxn_mult)
    x_comp_lst = deblend_components_all_asym(Ctotinv, Xd_obs, (V_skycont_r,), (V_skycont_c,))

    return x_comp_lst[1]
end
