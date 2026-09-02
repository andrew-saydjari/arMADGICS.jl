## Core per-spectrum pipeline (moved verbatim from pipeline.jl for testability)
# Relies on globals defined by the including script: reduxBase, almanacFile,
# logUniWaveAPOGEE, slvl_tuple, minRVres, RV_err_step, refine_iters, cache_dir, inject_cache_dir.

# gross global with reduxBase
function pipeline_single_spectra(argtup, prior_vec; caching=true, sky_caching=false, skyCont_off=false, skyLines_off=false, rv_split=true, ddstaronly=false, cache_dir=cache_dir, inject_cache_dir=inject_cache_dir)
    tele = argtup.tele
    mjd = argtup.mjd
    expnum = argtup.expnum
    adjfiberindx = argtup.adjfiberindx

    # V_skycont,chebmsk_exp,V_skyline_bright,V_skyline_faint,skymsk_bright,skymsk_faint,skymsk,V_starcont,V_starlines_refLSF, V_starlines, msk_starCor, V_dib_lst, V_dib_soft_lst, V_dib_noLSF_soft_lst = prior_vec
    chebmsk_exp, skymsk_bright, skymsk_faint, skymsk, V_starcont, V_starlines_refLSF, V_starlines, msk_starCor = prior_vec

    # M2: defaults so a failure record with correct shapes can be built even if
    # ingest itself throws (missing file, malformed almanac entry, ...)
    npix = length(chebmsk_exp)
    nstarcoef = size(V_starlines, 2)
    ingestBit = 0
    nSkyFibers = 0
    skyscale0 = NaN
    starscale0 = NaN
    snr = NaN
    simplemsk = falses(npix)
    fspec = fill(NaN, npix)
    fivar = fill(NaN, npix)

    try
        out = []

        # This could/should shift to a per night preprocessing
        # Get Sky Prior
        # ingest bit should start here to encode cases with not enough sky
        nSkyFibers, meanLocSky, meanLocSkyLines, VLocSkyLines, msk_local_skyLines = getSkyRough(reduxBase, tele, mjd, expnum, almanacFile)
        skyscale0 = nanzeromedian(meanLocSky)

        # Get the Exposure (Visit) Spectrum
        fspec, fivar, fmsk, metaexport = getExposure(reduxBase, tele, mjd, expnum, adjfiberindx)

        # M2/M3: sanity-check the spectrum before it can poison or kill the solve;
        # validate_exposure also strips non-finite flux/ivar and tiny-ivar pixels
        fmsk_clean, starscale0, ingestBit = validate_exposure(fspec, fivar, fmsk)

        simplemsk = fmsk_clean .& skymsk .& msk_local_skyLines
        # M1 fix: per-pixel snr is flux/sigma = flux.*sqrt.(ivar); the old
        # fspec./sqrt.(fivar) was flux*sigma (inverted) and, being computed unmasked,
        # let ivar=0 pixels inject +/-Inf into the median (nanzero* are now Inf-safe too).
        snr = median_snr(fspec, fivar, simplemsk)

        # combining with the sky/cheb masks can only shrink the good-pixel set
        if count(simplemsk) < INGEST_MIN_GOODPIX
            ingestBit |= 2^2
        end
        if ingest_fatal(ingestBit)
            println("Skipping spectrum (ingestBit=$ingestBit) for tele=$tele, mjd=$mjd, expnum=$expnum, adjfiberindx=$adjfiberindx")
            flush(stdout)
            return failed_pipeline_out(simplemsk, starscale0, skyscale0, fspec, fivar, nSkyFibers, snr, ingestBit, nstarcoef, collect(length.(slvl_tuple)))
        end

        push!(out, (count(simplemsk), starscale0, skyscale0, nanify(fspec[simplemsk], simplemsk), nanify(fivar[simplemsk], simplemsk), count(isnan.(fspec[simplemsk])), count(isnan.(fivar[simplemsk])), simplemsk, nSkyFibers, snr, ingestBit)) # 1

        if skyCont_off
            meanLocSky .= 0
            VLocSky .= 0
        end
        if skyLines_off
            meanLocSkyLines .= 0
            VLocSkyLines .= 0
            V_skyline_bright .= 0
            V_skyline_faint .= 0
        end

        # Make RV Mask for DD Model (else leave as simplemsk)
        rvmsk = copy(simplemsk)
        if ddstaronly
            lshift, rshift = extrema(vcat(map(x -> [extrema(x)...], slvl_tuple)...))
            lshift = floor(Int, lshift)
            rshift = ceil(Int, rshift)
            rvmsk .&= ShiftedArrays.circshift(msk_starCor, rshift)
            rvmsk .&= ShiftedArrays.circshift(msk_starCor, lshift)
        end

        ## Select data for use (might want to handle mean more generally)
        ## Mask full RV scan range
        Xd_obs = (fspec.-meanLocSky.-meanLocSkyLines)[rvmsk]
        wave_obs = logUniWaveAPOGEE[rvmsk]

        ## Set up residuals prior
        A = Diagonal(1 ./ fivar[rvmsk])
        Ainv = Diagonal(fivar[rvmsk])

        ## Set up priors
        # V_skyline_bright_c = V_skyline_bright
        # V_skyline_bright_r = V_skyline_bright_c[simplemsk,:]
        V_skyline_faint_c = VLocSkyLines
        V_skyline_faint_r = V_skyline_faint_c[rvmsk, :]
        # V_skyline_tot_c = hcat(V_skyline_bright_c,V_skyline_faint_c)
        # V_skyline_tot_r = hcat(V_skyline_bright_r,V_skyline_faint_r)
        V_skyline_tot_c = V_skyline_faint_c
        V_skyline_tot_r = V_skyline_faint_r
        # V_locSky_c = VLocSky
        # V_locSky_r = V_locSky_c[rvmsk,:]
        V_starCont_c = abs(starscale0) * V_starcont
        V_starCont_r = V_starCont_c[rvmsk, :]

        ## Solve RV of Star
        # compute stellar continuum to modify stellar line prior
        # Vcomb_skylines = hcat(V_skyline_tot_r,V_locSky_r,V_starCont_r);
        Vcomb_skylines = hcat(V_skyline_tot_r, V_starCont_r)
        Ctotinv_skylines = LowRankMultMatIP([Ainv, Vcomb_skylines], wood_precomp_mult_mat([Ainv, Vcomb_skylines], (size(Ainv, 1), size(V_starlines, 2))), wood_fxn_mult, wood_fxn_mult_mat!)
        x_comp_lst = deblend_components_all_asym(Ctotinv_skylines, Xd_obs, (V_starCont_r,), (V_starCont_c,))

        # starCont_Mscale_ref = x_comp_lst[1]
        starCont_Mscale = x_comp_lst[1][rvmsk]
        # Xd_obs = (fspec.-meanLocSky.-meanLocSkyLines.-starCont_Mscale_ref)[rvmsk];

        # ## Adjust the starContinuum covariance to be 1% of the "starScale"
        # starscalep5 = nanzeromedian(starCont_Mscale)
        # V_starCont_c = starCont_var*abs(starscalep5)*V_starcont
        # V_starCont_r = V_starCont_c[rvmsk,:]

        # now take out the skylines to be included in the scanning
        # Vcomb_cur = hcat(V_locSky_r,V_starCont_r);
        Vcomb_cur = V_starCont_r
        Ctotinv_cur = LowRankMultMatIP([Ainv, Vcomb_cur], wood_precomp_mult_mat([Ainv, Vcomb_cur], (size(Ainv, 1), size(V_starlines, 2))), wood_fxn_mult, wood_fxn_mult_mat!)

        # compute delta chi2 for adding skylines (helps normalize the joint chi2 below with starLines)
        chi2skyoffset = woodbury_update_inv_tst(
            LowRankMultMatIP([Ainv, Vcomb_cur], wood_precomp_mult_mat([Ainv, Vcomb_cur], (size(Ainv, 1), size(V_skyline_tot_r, 2))), wood_fxn_mult, wood_fxn_mult_mat!),
            Xd_obs,
            V_skyline_tot_r
        )

        pre_Vslice = zeros(count(rvmsk), size(V_starlines, 2))
        chi2_wrapper_partial = if rv_split
            AinvV1 = Ctotinv_cur * V_skyline_tot_r
            XdAinvV1 = reshape(Xd_obs, 1, :) * AinvV1
            V1TAinvV1 = V_skyline_tot_r' * AinvV1
            Base.Fix2(chi2_wrapper_split, (rvmsk, Ctotinv_cur, Xd_obs, starCont_Mscale, V_starlines, pre_Vslice, AinvV1, XdAinvV1, V1TAinvV1, chi2skyoffset))
        else
            error("rv_split must be true")
        end
        lout = sampler_1d_hierarchy_var(chi2_wrapper_partial, slvl_tuple, minres=minRVres, stepx=RV_err_step)
        svalc = lout[1][3]
        push!(out, lout) # 2

        # re-estiamte starScale before re-creating the priors with the new finalRV msk
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r, V_starlines_ru = update_Ctotinv_Vstarstarlines_asym(svalc, Ctotinv_skylines.matList[1], rvmsk, starCont_Mscale, Vcomb_skylines, V_starlines, V_starlines_refLSF)
        x_comp_lst = deblend_components_all_asym(Ctotinv_fut, Xd_obs, (V_starCont_r,), (V_starCont_c,))

        # Change data mask based on final inferred RV
        finalmsk = copy(simplemsk)
        if ddstaronly
            rvshift = sign(lout[1][3]) * ceil(abs(lout[1][3]))
            finalmsk .&= ShiftedArrays.circshift(msk_starCor, rvshift)
        end

        starCont_Mscale_ref = x_comp_lst[1]
        starCont_Mscale = starCont_Mscale_ref[finalmsk]
        Xd_obs = (fspec.-meanLocSky.-meanLocSkyLines)[finalmsk]
        wave_obs = logUniWaveAPOGEE[finalmsk]

        starscale1 = nanzeromedian(starCont_Mscale)

        ## Set up residuals prior
        A = Diagonal(1 ./ fivar[finalmsk])
        Ainv = Diagonal(fivar[finalmsk])

        ## Set up priors
        # V_skyline_bright_r = V_skyline_bright_c[simplemsk,:]
        V_skyline_faint_r = V_skyline_faint_c[finalmsk, :]
        V_skyline_tot_r = V_skyline_faint_r
        # V_locSky_r = V_locSky_c[finalmsk,:]
        V_starCont_c = abs(starscale1) * V_starcont
        V_starCont_r = V_starCont_c[finalmsk, :]

        # Vcomb_skylines = hcat(V_skyline_tot_r,V_locSky_r,V_starCont_r);
        Vcomb_skylines = hcat(V_skyline_tot_r, V_starCont_r)
        Ctotinv_skylines = LowRankMultMatIP([Ainv, Vcomb_skylines], wood_precomp_mult_mat([Ainv, Vcomb_skylines], (size(Ainv, 1), size(V_starlines, 2))), wood_fxn_mult, wood_fxn_mult_mat!)

        x_comp_lst = deblend_components_all(Ctotinv_skylines, Xd_obs, (V_starCont_r,))
        starCont_Mscale = x_comp_lst[1]

        # update the Ctotinv to include the stellar line component (iterate to refine starCont_Mscale)
        for i = 1:refine_iters
            Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r, V_starlines_ru = update_Ctotinv_Vstarstarlines_asym(svalc, Ctotinv_skylines.matList[1], finalmsk, starCont_Mscale, Vcomb_skylines, V_starlines, V_starlines_refLSF)
            x_comp_lst = deblend_components_all(Ctotinv_fut, Xd_obs, (V_starCont_r,))
            starCont_Mscale = x_comp_lst[1]
        end
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r, V_starlines_ru = update_Ctotinv_Vstarstarlines_asym(svalc, Ctotinv_skylines.matList[1], finalmsk, starCont_Mscale, Vcomb_skylines, V_starlines, V_starlines_refLSF)

        # do a component save without the 15273 DIB
        # the extra Vstarlines_r is duplicated work if a pure dd model, but helps compare flux conservation in both cases
        # x_comp_lst = deblend_components_all_asym_tot(Ctotinv_fut, Xd_obs, 
        #     (A, V_skyline_faint_r, V_locSky_r, V_starCont_r, V_starlines_r, V_starlines_r, V_starlines_r),
        #     (A, V_skyline_faint_r, V_locSky_c, V_starCont_c, V_starlines_ru, V_starlines_c, I),
        # )
        x_comp_lst = deblend_components_all_asym_tot(Ctotinv_fut, Xd_obs,
            (A, V_skyline_faint_r, V_skyline_faint_r, V_starCont_r, V_starlines_r, V_starlines_r, V_starlines_r),
            (A, V_skyline_faint_r, V_skyline_faint_c, V_starCont_c, V_starlines_ru, V_starlines_c, I),
        )



        x_comp_out = []
        push!(x_comp_out, nanify(x_comp_lst[1] .* sqrt.(fivar[finalmsk]), finalmsk)) #z-scored residuals
        push!(x_comp_out, nanify(x_comp_lst[1], finalmsk)) #residuals
        # push!(x_comp_out,nanify(x_comp_lst[2][skymsk_bright[finalmsk]],finalmsk .& skymsk_bright)) #bright sky lines
        push!(x_comp_out, nanify(x_comp_lst[2][skymsk_faint[finalmsk]] .+ meanLocSkyLines[finalmsk.&skymsk_faint], finalmsk .& skymsk_faint)) #faint sky lines
        push!(x_comp_out, nanify(0 .* x_comp_lst[3][chebmsk_exp] .+ meanLocSky[chebmsk_exp], chebmsk_exp)) #sky continuum #hacked to skylines times zero
        push!(x_comp_out, nanify(x_comp_lst[4][chebmsk_exp], chebmsk_exp)) #star continuum
        push!(x_comp_out, x_comp_lst[6:end]...) # starLines, starlines coefficients, and totchi2
        push!(x_comp_out, apVisit_from_components(fspec, x_comp_out[3], x_comp_out[4], x_comp_out[5], finalmsk, starscale1)) #apVisit analog (M3: continuum-floored)
        push!(x_comp_out, finalmsk) # final mask
        push!(x_comp_out, V_starlines_refLSF[:, :, 6] * x_comp_lst[7]) # Restframe StarLine component with reference LSF

        skyscale1 = nanzeromedian(x_comp_out[4])
        dvec = (fspec .- (x_comp_out[2] .+ x_comp_out[3] .+ x_comp_out[4] .+ x_comp_out[5] .* (1 .+ nanify(x_comp_lst[5], finalmsk)))) ./ fspec
        chi2res = x_comp_lst[1]' * (Ainv * x_comp_lst[1])
        push!(out, (chi2res, nanzeroiqr(dvec), count(finalmsk), starscale1, skyscale1)) # 3
        push!(out, x_comp_out) # 4
        dflux_starlines = sqrt_nan.(get_diag_posterior_from_prior_asym(Ctotinv_fut, V_starlines_c, V_starlines_r))
        push!(out, dflux_starlines) # 5

        # # prepare multiplicative factors for DIB prior
        # x_comp_lst = deblend_components_all(Ctotinv_fut, Xd_obs, (V_starCont_r,V_starlines_r))
        # starCont_Mscale = x_comp_lst[1]
        # starFull_Mscale = starCont_Mscale.+x_comp_lst[2]

        # Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r, V_starlines_ru = update_Ctotinv_Vstarstarlines_asym(svalc,Ctotinv_skylines.matList[1],finalmsk,starCont_Mscale,Vcomb_skylines,V_starlines,V_starlines_refLSF)
        # Ctotinv_cur, Ctotinv_fut = Ctotinv_fut, Ctotinv_cur; Vcomb_cur, Vcomb_fut = Vcomb_fut, Vcomb_cur # swap to updated covariance finally

        # # currently, this is modeling each DIB seperately... I think we want to change this later, just easier parallel structure
        # for dib_ind = 1:length(dib_center_lst) # eventually need to decide if these are cumulative or not
        #     V_dib = V_dib_lst[dib_ind_prior[dib_ind]]
        #     V_dib_soft = V_dib_soft_lst[dib_ind_prior[dib_ind]]
        #     # V_dib_noLSF = V_dib_noLSF_lst[dib_ind_prior[dib_ind]]
        #     V_dib_noLSF_soft = V_dib_noLSF_soft_lst[dib_ind_prior[dib_ind]]

        #     pre_Vslice = zeros(count(finalmsk),size(V_dib,2))
        #     lvltuple_dib = lvltuple_lst[dib_ind]
        #     dib_center = dib_center_lst[dib_ind]
        #     # scan_offset deprecated with individual DIB priors
        #     scan_offset = 0 #findmin(abs.(logUniWaveAPOGEE.-dib_center_lst[dib_ind]))[2].-findmin(abs.(logUniWaveAPOGEE.-dib_center_lst[1]))[2]

        #     ## Solve DIB parameters for just a single DIB
        #     # one of the main questions is how many times to compute components and where
        #     chi2_wrapper_partial = Base.Fix2(chi2_wrapper2d,(finalmsk,Ctotinv_cur,Xd_obs,wave_obs,starFull_Mscale,Vcomb_cur,V_dib,pre_Vslice,dib_center,scan_offset))
        #     lout = sampler_2d_hierarchy_var(chi2_wrapper_partial,lvltuple_dib,step1=DIB_pix_err_step,step2=DIB_sig_err_step,minres1=minDIBvelres,minres2=minDIBsigres)
        #     opt_tup = lout[1][3]
        #     push!(out,lout) # 6

        #     ## Shift the marginalization sampling (should this be wrapped inside the function?)
        #     # especially because we need to do bounds handling
        #     svalMarg = svalMarg0 .+ opt_tup[1]
        #     sigMarg = shift_trim_range(sigMarg0,opt_tup[2]; minv=minSigval, maxv=maxSigval)
        #     samp_lst = Iterators.product(svalMarg,sigMarg)

        #     intupf = (finalmsk,Ctotinv_cur,Xd_obs,wave_obs,starFull_Mscale,Vcomb_cur,V_dib,pre_Vslice,dib_center,scan_offset)
        #     chi2lst, fluxlst, dfluxlst = sample_chi2_flux_dflux(samp_lst,intupf) #shouldn't this take chi2_wrapper_partial as an argument?
        #     refchi2val = minimum(chi2lst) #this should just be set to the min found at the 2d step
        #     lout = marginalize_flux_err(chi2lst, fluxlst, dfluxlst, refchi2val)
        #     push!(out,lout) # 7

        #     # Compute some final components for export (still need to implement DIB iterative refinement)
        #     Ctotinv_fut, Vcomb_fut, V_dibc, V_dibr = update_Ctotinv_Vdib_asym(
        #         opt_tup,Ctotinv_cur.matList[1],finalmsk,starFull_Mscale,Vcomb_cur,V_dib_soft,V_dib_noLSF_soft,scan_offset)

        #     x_comp_lst = deblend_components_all_asym_tot(Ctotinv_fut, Xd_obs, 
        #         (A, V_skyline_faint_r, V_locSky_r, V_starCont_r, V_dibr, V_starlines_r, V_dibr),
        #         (A, V_skyline_faint_r, V_locSky_c, V_starCont_c, V_dibr, V_starlines_c, V_dibc),
        #     )

        #     x_comp_out = []
        #     push!(x_comp_out,nanify(x_comp_lst[1].*sqrt.(fivar[finalmsk]),finalmsk)) # z-scored residuals
        #     push!(x_comp_out,nanify(x_comp_lst[1],finalmsk)) # residuals
        #     # push!(x_comp_out,nanify(x_comp_lst[2][skymsk_bright[finalmsk]],finalmsk .& skymsk_bright)) #bright sky lines
        #     push!(x_comp_out,nanify(x_comp_lst[2][skymsk_faint[finalmsk]].+meanLocSkyLines[finalmsk .& skymsk_faint],finalmsk .& skymsk_faint)) # faint sky lines
        #     push!(x_comp_out,nanify(x_comp_lst[3][chebmsk_exp].+meanLocSky[chebmsk_exp],chebmsk_exp)) #sky continuum
        #     push!(x_comp_out,nanify(x_comp_lst[4][chebmsk_exp],chebmsk_exp)) #star continuum
        #     push!(x_comp_out,nanify(x_comp_lst[5],finalmsk)) # dib flux 
        #     push!(x_comp_out,x_comp_lst[6:end]...) # starLines, dib, and totchi2

        #     chi2res = x_comp_lst[1]'*(Ainv*x_comp_lst[1])
        #     push!(out,(chi2res,)) # 8

        #     push!(out,x_comp_out) # 9
        # end

        return out
    catch e
        # M2: per-spectrum failures are non-fatal — record a failure code in the
        # batch output instead of killing the whole pmap (the old code rethrew here)
        ingestBit |= INGEST_RUNTIME_ERROR_BIT
        println("Error in pipeline_single_spectra for tele=$tele, mjd=$mjd, expnum=$expnum, adjfiberindx=$adjfiberindx (recorded ingestBit=$ingestBit): ", sprint(showerror, e))
        flush(stdout)
        return failed_pipeline_out(simplemsk, starscale0, skyscale0, fspec, fivar, nSkyFibers, snr, ingestBit, nstarcoef, collect(length.(slvl_tuple)))
    end
end

# M3: apVisit analog = (flux - sky) / starContinuum. Pixels whose inferred continuum
# is non-finite or below rel_floor*|cont_scale| (cont_scale = starscale1, the median
# inferred continuum) are set to NaN instead of being divided into +/-huge/Inf values
# (broken or very faint fibers have near-zero continuum). A non-finite cont_scale
# masks the whole vector.
function apVisit_from_components(fspec, skyLines, skyCont, starCont, finalmsk, cont_scale; rel_floor=1e-3)
    cont_floor = rel_floor * abs(cont_scale)
    ok = finalmsk .& isfinite.(starCont) .& (abs.(starCont) .> cont_floor)
    return nanify(((fspec .- (skyLines .+ skyCont)) ./ starCont)[ok], ok)
end

# RV_flag value marking a spectrum that never entered the RV scan (see ingestBit for why)
const INGEST_FAIL_RV_FLAG = 2^6

"""
    failed_pipeline_out(simplemsk, starscale0, skyscale0, fspec, fivar, nSkyFibers, snr, ingestBit, nstarcoef, lvllens)

Build a placeholder `out` with EXACTLY the same nesting/shapes as a successful
`pipeline_single_spectra` return, so `extractor`/`multi_spectra_batch` can save
mixed success/failure batches. All science quantities are NaN (counts 0, masks
false); the failure reason is carried in the `ingestBit` column (bit codes in
src/ingest.jl) and `RV_flag` is set to `INGEST_FAIL_RV_FLAG`.
"""
function failed_pipeline_out(simplemsk, starscale0, skyscale0, fspec, fivar, nSkyFibers, snr, ingestBit, nstarcoef, lvllens)
    npix = length(simplemsk)
    out = []
    # 1: meta block (mirrors the success push, incl. the new ingestBit column)
    push!(out, (count(simplemsk), starscale0, skyscale0,
        nanify(fspec[simplemsk], simplemsk), nanify(fivar[simplemsk], simplemsk),
        count(isnan.(fspec[simplemsk])), count(isnan.(fivar[simplemsk])),
        simplemsk, nSkyFibers, snr, ingestBit))
    # 2: RV block (sampler_1d_hierarchy_var shape)
    lvlouts = [((NaN, NaN, NaN, NaN, 1, INGEST_FAIL_RV_FLAG), fill(NaN, n), fill(NaN, n)) for n in lvllens]
    push!(out, ((NaN, NaN, NaN, NaN, 1, INGEST_FAIL_RV_FLAG, NaN), lvlouts))
    # 3: (chi2res, avg_flux_conservation, final_pix_cnt, starscale1, skyscale1)
    push!(out, (NaN, NaN, 0, NaN, NaN))
    # 4: component block (11 entries mirroring x_comp_out)
    x_comp_out = []
    for i = 1:6
        push!(x_comp_out, fill(NaN, npix)) # z-resid, resid, skyLines, skyCont, starCont, starLines
    end
    push!(x_comp_out, fill(NaN, nstarcoef)) # starLines coefficients
    push!(x_comp_out, NaN) # tot_p5chi2
    push!(x_comp_out, fill(NaN, npix)) # apVisit analog
    push!(x_comp_out, falses(npix)) # final mask
    push!(x_comp_out, fill(NaN, npix)) # restframe starLines
    push!(out, x_comp_out)
    # 5: dflux_starlines
    push!(out, fill(NaN, npix))
    return out
end
