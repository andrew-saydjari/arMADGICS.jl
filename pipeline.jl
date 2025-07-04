## This is the main pipeline that will batch over APOGEE files
# Author - Andrew Saydjari, CfA

import Pkg; using Dates; t0 = now(); t_then = t0;
using InteractiveUtils; versioninfo()
Pkg.instantiate(); Pkg.precompile()
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Package activation took $dt"); t_then = t_now; flush(stdout)
using BLISBLAS
using Distributed, SlurmClusterManager, Suppressor, DataFrames, DelimitedFiles
addprocs(SlurmManager(),exeflags=["--project=./"])
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker allocation took $dt"); t_then = t_now; flush(stdout)
println("Running Main on ", gethostname()); flush(stdout)

@everywhere begin
    using BLISBLAS
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays, JLD2, FileIO
    using Interpolations, SparseArrays, ParallelDataTransfer, AstroTime, Suppressor
    using ThreadPinning, ApogeeReduction, DataFrames

    prior_dir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/"
    src_dir = "./"
    include(src_dir*"src/utils.jl")
    include(src_dir*"src/gridSearch.jl")
    include(src_dir*"src/componentAndPosteriors.jl")
    include(src_dir*"src/fileNameHandling.jl")
    include(src_dir*"src/ingest.jl")
    include(src_dir*"src/lowRankPrescription.jl")
    include(src_dir*"src/marginalizeEW.jl")
    include(src_dir*"src/spectraInterpolation.jl")
    include(src_dir*"src/chi2Wrappers.jl")
    
    using StatsBase, ProgressMeter
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

println(BLAS.get_config()); flush(stdout)
using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

# These global allocations for the injest are messy... but we plan on changing the ingest
# relatively soon... so don't worry for now.
@everywhere begin
    refine_iters = 5
    ddstaronly = true
    runlist_range = 1:600 # 295, 245, 335, 101
    batchsize = 10 #100

    # Step Size for Chi2 Surface Error Bars
    RV_err_step = 4
    # DIB_pix_err_step = 3 # consider increasing to 4 (self consistency + LSF test)
    # DIB_sig_err_step = 3
    # Flux marginalize region
    # sigMarg0 = -50//100:10//100:50//100
    # svalMarg0 = -0//10:1//10:0//10;

    cache_dir = "../local_cache/"
    out_dir="../outdir/"
    inject_cache_dir = replace(cache_dir,"local_cache/"=>"inject_local_cache/")

    reduxBase = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/sandbox51/airflow-ApogeeReduction.jl/daily/outdir"
    almanacFile = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/sandbox51/airflow-ApogeeReduction.jl/daily/outdir/almanac/objects_60834.h5"
    tele = "lco"
    mjd = "60834"

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # Input List (not really a prior, but an input file we search for stars conditioned on)
    # prior_dict["runlists"] = prior_dir*"2024_03_11/inject_15273only_295_real/injection_input_lst_"
    prior_dict["runlists"] = prior_dir*"2024_03_15/outlists/star/dr17_dr17_star_input_lst_msked_" # repackaged for cross platform/version from 2024_03_05

    # Sky Priors
    prior_dict["skycont"] = prior_dir*"2024_02_21/arMADGICS.jl/src/prior_build/sky_priors/APOGEE_skycont_svd_30_f"
    prior_dict["skyLines_bright"] = prior_dir*"2024_02_21/arMADGICS.jl/src/prior_build/sky_priors/APOGEE_skyline_bright_GSPICE_svd_120_f"
    prior_dict["skyLines_faint"] = prior_dir*"2024_02_21/arMADGICS.jl/src/prior_build/sky_priors/APOGEE_skyline_faint_GSPICE_svd_120_f"

    # Star Priors
    # prior_dict["starCont"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/star_priors/APOGEE_starcont_svd_60_f"
    prior_dict["chebmsk"] = prior_dir*"2025_06_16/chebmsk_exp.h5"
    prior_dict["starCont"] = prior_dir*"2025_06_16/APOGEE_starcont_svd_60_rough.h5"
    prior_dict["starLines_refLSF"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/starLine_priors_norm94/APOGEE_stellar_kry_50_subpix_th_22500.h5"
    # prior_dict["starLines_LSF"] = prior_dir*"2024_03_16/arMADGICS.jl/src/prior_build/starLine_priors_norm94_dd/APOGEE_starCor_svd_50_subpix_f" # DD Version
    # prior_dict["starLines_LSF"] = prior_dir*"2024_02_21/arMADGICS.jl/src/prior_build/starLine_priors_norm94/APOGEE_stellar_kry_50_subpix_f" # TH Version

    # # DIB Priors
    # dib_waves = [15273, 15672]
    # for dib in dib_waves
    #     # prior_dict["DIB_noLSF_$(dib)"] = prior_dir*"2024_03_05/arMADGICS.jl/src/prior_build/dib_priors/precomp_dust_1_$(dib)_analyticDeriv_stiff.h5"
    #     scratch1/working/
    #     prior_dict["DIB_noLSF_soft_$(dib)"] = prior_dir*"2024_03_05/arMADGICS.jl/src/prior_build/dib_priors/precomp_dust_3_$(dib)_analyticDeriv_soft.h5"
    #     prior_dict["DIB_LSF_$(dib)"] = prior_dir*"2024_03_05/arMADGICS.jl/src/prior_build/dib_priors/precomp_dust_1_$(dib)_analyticDerivLSF_stiff_"
    #     prior_dict["DIB_LSF_soft_$(dib)"] = prior_dir*"2024_03_05/arMADGICS.jl/src/prior_build/dib_priors/precomp_dust_3_$(dib)_analyticDerivLSF_soft_"
    # end

    # # maps dib scans to dib_waves
    # dib_ind_prior = Dict{Int,Int}()
    # dib_ind_prior[1] = 1
    # dib_ind_prior[2] = 1
    # dib_ind_prior[3] = 2
    # dib_ind_prior[4] = 2
end

# it would be great to move this into a parameter file that is read for each run
# I should revisit the error bars in the context of chi2 versus frame number trends
@everywhere begin
    # minimum step sizes of the priors (might want to store in the file structures of those priors and read in)
    minRVres = 1//10
    # minDIBvelres = 1//10
    # minDIBsigres = 1//100

    ## sigrng is somewhat counterintuitively defined in the chi2Wrappers.jl file
    minSigval, maxSigval = extrema(sigrng)

    # Star Wave
    lvl1 = -70:1//2:70
    lvl2 = -8:2//10:8
    lvl3 = -3:1//10:3
    slvl_tuple = (lvl1,lvl2,lvl3)
    # tuple1dprint(slvl_tuple)

    # # (Wave, Sig) DIB
    # dib_center_lst = map(x->dib_waves[dib_ind_prior[x]],1:length(dib_ind_prior)) # not clear we need this anymore since we don't scanOffset
    # lvl1d_15273wide = ((-137:4:150),(18//10:18//10))
    # lvl1d_15672wide = ((-54:4:150),(18//10:18//10))
    # lvl1d_narrow = ((-54:4:54),(18//10:18//10))
    # lvl1d = ((-150:4:150),(18//10:18//10))
    # lvl2d = ((0:0), (-7//5:4//100:11//5))
    # lvl3d = ((-18:2//10:18), (0:0))
    # lvl4d = ((0:0), (-90//100:2//100:90//100))
    # lvl5d = ((-1:2//10:1), (0:0))
    # lvl6d = ((0:0), (-10//100:2//100:10//100))
    # lvl7d = ((-6//10:2//10:6//10), (0:0))
    # lvl8d = ((0:0), (-6//100:2//100:6//100))
    # lvl9d = ((-4//10:1//10:4//10), (-4//100:1//100:4//100));
    # lvltuple = (lvl1d, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    # lvltuple_15273wide = (lvl1d_15273wide, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    # lvltuple_15672wide = (lvl1d_15672wide, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    # lvltuple_narrow = (lvl1d_narrow, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    # lvltuple_lst = [lvltuple_15273wide, lvltuple_narrow, lvltuple_15672wide, lvltuple_narrow]
    # # tuple2dprint(lvltuple)
end

@everywhere begin
    logUniWaveAPOGEE = 10 .^range((start=4.179-125*6.0e-6),step=6.0e-6,length=8575+125)
    minw, maxw = extrema(logUniWaveAPOGEE)
    
    c = 299792.458; # in km/s
    delLog = 6e-6; 

    Xd_stack = zeros(3*2048)
    Xd_std_stack = zeros(3*2048)
    waveobs_stack = zeros(3*2048)
    waveobs_stack_old = zeros(3*2048)
    pixmsk_stack = zeros(Int,3*2048)
    fullBit = zeros(Int,3*2048);
    outvec = zeros(length(logUniWaveAPOGEE))
    outvar = zeros(length(logUniWaveAPOGEE))
    cntvec = zeros(Int,length(logUniWaveAPOGEE));
end

@everywhere begin
    function pipeline_single_spectra(argtup, prior_vec; caching=true, sky_caching=false, skyCont_off=false, skyLines_off=false, rv_split=true, ddstaronly=false, cache_dir=cache_dir, inject_cache_dir=inject_cache_dir)
        reduxBase, tele, mjd, expnum, adjfiberindx = argtup[2:end]
        # V_skycont,chebmsk_exp,V_skyline_bright,V_skyline_faint,skymsk_bright,skymsk_faint,skymsk,V_starcont,V_starlines_refLSF, V_starlines, msk_starCor, V_dib_lst, V_dib_soft_lst, V_dib_noLSF_soft_lst = prior_vec
        chebmsk_exp,skymsk_bright,skymsk_faint,skymsk,V_starcont,V_starlines_refLSF,V_starlines,msk_starCor = prior_vec
        out = []

        # This could/should shift to a per night preprocessing
        # Get Sky Prior
        meanLocSky, meanLocSkyLines, VLocSkyLines, msk_local_skyLines = getSkyRough(reduxBase, tele, mjd, expnum)
        skyscale0 = nanzeromedian(meanLocSky)

        # Get the Exposure (Visit) Spectrum
        fspec, fivar, fmsk, metaexport = getExposure(reduxBase, tele, mjd, expnum, adjfiberindx)
        starscale0 = nanzeromedian(fspec)
        
        simplemsk = fmsk .& skymsk .& msk_local_skyLines;
        
        push!(out,(count(simplemsk), starscale0, skyscale0, nanify(fspec[simplemsk],simplemsk), nanify(fivar[simplemsk],simplemsk),count(isnan.(fspec[simplemsk])),count(isnan.(fivar[simplemsk])),simplemsk)) # 1

        if skyCont_off
            meanLocSky.=0
            VLocSky.=0
        end
        if skyLines_off
            meanLocSkyLines.=0
            VLocSkyLines.=0
            V_skyline_bright.=0
            V_skyline_faint.=0
        end

        # Make RV Mask for DD Model (else leave as simplemsk)
        rvmsk = copy(simplemsk)
        if ddstaronly
            lshift, rshift = extrema(vcat(map(x->[extrema(x)...],slvl_tuple)...)); lshift = floor(Int,lshift); rshift = ceil(Int,rshift)
            rvmsk .&= ShiftedArrays.circshift(msk_starCor,rshift)
            rvmsk .&= ShiftedArrays.circshift(msk_starCor,lshift)
        end

        ## Select data for use (might want to handle mean more generally)
        ## Mask full RV scan range
        Xd_obs = (fspec.-meanLocSky.-meanLocSkyLines)[rvmsk];
        wave_obs = logUniWaveAPOGEE[rvmsk]

        ## Set up residuals prior
        A = Diagonal(1 ./fivar[rvmsk]);
        Ainv = Diagonal(fivar[rvmsk]);

        ## Set up priors
        # V_skyline_bright_c = V_skyline_bright
        # V_skyline_bright_r = V_skyline_bright_c[simplemsk,:]
        V_skyline_faint_c = VLocSkyLines
        V_skyline_faint_r = V_skyline_faint_c[rvmsk,:]
        # V_skyline_tot_c = hcat(V_skyline_bright_c,V_skyline_faint_c)
        # V_skyline_tot_r = hcat(V_skyline_bright_r,V_skyline_faint_r)
        V_skyline_tot_c = V_skyline_faint_c
        V_skyline_tot_r = V_skyline_faint_r
        # V_locSky_c = VLocSky
        # V_locSky_r = V_locSky_c[rvmsk,:]
        V_starCont_c = abs(starscale0)*V_starcont
        V_starCont_r = V_starCont_c[rvmsk,:]

        ## Solve RV of Star
        # compute stellar continuum to modify stellar line prior
        # Vcomb_skylines = hcat(V_skyline_tot_r,V_locSky_r,V_starCont_r);
        Vcomb_skylines = hcat(V_skyline_tot_r,V_starCont_r);
        Ctotinv_skylines = LowRankMultMatIP([Ainv,Vcomb_skylines],wood_precomp_mult_mat([Ainv,Vcomb_skylines],(size(Ainv,1),size(V_starlines,2))),wood_fxn_mult,wood_fxn_mult_mat!);
        x_comp_lst = deblend_components_all_asym(Ctotinv_skylines, Xd_obs, (V_starCont_r, ), (V_starCont_c, ))

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
        Ctotinv_cur = LowRankMultMatIP([Ainv,Vcomb_cur],wood_precomp_mult_mat([Ainv,Vcomb_cur],(size(Ainv,1),size(V_starlines,2))),wood_fxn_mult,wood_fxn_mult_mat!);

        # compute delta chi2 for adding skylines (helps normalize the joint chi2 below with starLines)
        chi2skyoffset = woodbury_update_inv_tst(
            LowRankMultMatIP([Ainv,Vcomb_cur],wood_precomp_mult_mat([Ainv,Vcomb_cur],(size(Ainv,1),size(V_skyline_tot_r,2))),wood_fxn_mult,wood_fxn_mult_mat!),
            Xd_obs,
            V_skyline_tot_r
        )
  
        pre_Vslice = zeros(count(rvmsk),size(V_starlines,2))
        chi2_wrapper_partial = if rv_split
            AinvV1 = Ctotinv_cur*V_skyline_tot_r
            XdAinvV1 = reshape(Xd_obs,1,:)*AinvV1
            V1TAinvV1 = V_skyline_tot_r'*AinvV1
            Base.Fix2(chi2_wrapper_split,(rvmsk,Ctotinv_cur,Xd_obs,starCont_Mscale,V_starlines,pre_Vslice,AinvV1,XdAinvV1,V1TAinvV1,chi2skyoffset))
        else
            error("rv_split must be true")
        end
        lout = sampler_1d_hierarchy_var(chi2_wrapper_partial,slvl_tuple,minres=minRVres,stepx=RV_err_step)
        svalc = lout[1][3]
        push!(out,lout) # 2

        # re-estiamte starScale before re-creating the priors with the new finalRV msk
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r, V_starlines_ru = update_Ctotinv_Vstarstarlines_asym(svalc,Ctotinv_skylines.matList[1],rvmsk,starCont_Mscale,Vcomb_skylines,V_starlines,V_starlines_refLSF)
        x_comp_lst = deblend_components_all_asym(Ctotinv_fut, Xd_obs, (V_starCont_r, ), (V_starCont_c, ))

        # Change data mask based on final inferred RV
        finalmsk = copy(simplemsk)
        if ddstaronly
            rvshift = sign(lout[1][3])*ceil(abs(lout[1][3]))
            finalmsk .&= ShiftedArrays.circshift(msk_starCor,rvshift)
        end

        starCont_Mscale_ref = x_comp_lst[1]
        starCont_Mscale = starCont_Mscale_ref[finalmsk]
        Xd_obs = (fspec.-meanLocSky.-meanLocSkyLines)[finalmsk]
        wave_obs = logUniWaveAPOGEE[finalmsk]

        starscale1 = nanzeromedian(starCont_Mscale)

        ## Set up residuals prior
        A = Diagonal(1 ./fivar[finalmsk]);
        Ainv = Diagonal(fivar[finalmsk]);

        ## Set up priors
        # V_skyline_bright_r = V_skyline_bright_c[simplemsk,:]
        V_skyline_faint_r = V_skyline_faint_c[finalmsk,:]
        V_skyline_tot_r = V_skyline_faint_r
        # V_locSky_r = V_locSky_c[finalmsk,:]
        V_starCont_c = abs(starscale1)*V_starcont
        V_starCont_r = V_starCont_c[finalmsk,:]

        # Vcomb_skylines = hcat(V_skyline_tot_r,V_locSky_r,V_starCont_r);
        Vcomb_skylines = hcat(V_skyline_tot_r,V_starCont_r);
        Ctotinv_skylines = LowRankMultMatIP([Ainv,Vcomb_skylines],wood_precomp_mult_mat([Ainv,Vcomb_skylines],(size(Ainv,1),size(V_starlines,2))),wood_fxn_mult,wood_fxn_mult_mat!);

        x_comp_lst = deblend_components_all(Ctotinv_skylines, Xd_obs, (V_starCont_r,))
        starCont_Mscale = x_comp_lst[1]

        # update the Ctotinv to include the stellar line component (iterate to refine starCont_Mscale)
        for i=1:refine_iters
            Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r, V_starlines_ru = update_Ctotinv_Vstarstarlines_asym(svalc,Ctotinv_skylines.matList[1],finalmsk,starCont_Mscale,Vcomb_skylines,V_starlines,V_starlines_refLSF)
            x_comp_lst = deblend_components_all(Ctotinv_fut, Xd_obs, (V_starCont_r, ))
            starCont_Mscale = x_comp_lst[1]
        end
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r, V_starlines_ru = update_Ctotinv_Vstarstarlines_asym(svalc,Ctotinv_skylines.matList[1],finalmsk,starCont_Mscale,Vcomb_skylines,V_starlines,V_starlines_refLSF)
        
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
        push!(x_comp_out,nanify(x_comp_lst[1].*sqrt.(fivar[finalmsk]),finalmsk)) #z-scored residuals
        push!(x_comp_out,nanify(x_comp_lst[1],finalmsk)) #residuals
        # push!(x_comp_out,nanify(x_comp_lst[2][skymsk_bright[finalmsk]],finalmsk .& skymsk_bright)) #bright sky lines
        push!(x_comp_out,nanify(x_comp_lst[2][skymsk_faint[finalmsk]].+meanLocSkyLines[finalmsk .& skymsk_faint],finalmsk .& skymsk_faint)) #faint sky lines
        push!(x_comp_out,nanify(0 .*x_comp_lst[3][chebmsk_exp].+meanLocSky[chebmsk_exp],chebmsk_exp)) #sky continuum #hacked to skylines times zero
        push!(x_comp_out,nanify(x_comp_lst[4][chebmsk_exp],chebmsk_exp)) #star continuum
        push!(x_comp_out,x_comp_lst[6:end]...) # starLines, starlines coefficients, and totchi2
        push!(x_comp_out,nanify(((fspec.-(x_comp_out[3].+x_comp_out[4]))./x_comp_out[5])[finalmsk],finalmsk)) #apVisit analog
        push!(x_comp_out,finalmsk) # final mask
        push!(x_comp_out,V_starlines_refLSF[:,:,6]*x_comp_lst[7]) # Restframe StarLine component with reference LSF

        skyscale1 = nanzeromedian(x_comp_out[4])
        dvec = (fspec .-(x_comp_out[2].+x_comp_out[3].+x_comp_out[4].+x_comp_out[5].*(1 .+ nanify(x_comp_lst[5],finalmsk))))./fspec;
        chi2res = x_comp_lst[1]'*(Ainv*x_comp_lst[1])
        push!(out,(chi2res,nanzeroiqr(dvec),count(finalmsk),starscale1,skyscale1)) # 3
        push!(out,x_comp_out) # 4
        dflux_starlines = sqrt_nan.(get_diag_posterior_from_prior_asym(Ctotinv_fut, V_starlines_c, V_starlines_r))
        push!(out,dflux_starlines) # 5
                
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
    end
end

@everywhere begin
    function multi_spectra_batch(indsubset; out_dir=out_dir, ddstaronly=ddstaronly)
        ### Set up
        out = []
        startind = indsubset[1][1]
        tele = indsubset[1][4]
        fiberindx = indsubset[1][end]
        teleind = (tele[1:6] == "lco25m") ? 2 : 1
        adjfibindx = (teleind-1)*300 + fiberindx

        ### Save and cache restart handling
        savename = join([out_dir,lpad(adjfibindx,3,"0"),"arMADGICS_fiber_"*lpad(adjfibindx,3,"0")*"_batch_"*lpad(startind,7,"0")*".h5"],"/")
        dirName = splitdir(savename)[1]
        if !ispath(dirName)
            mkpath(dirName)
        end
        if !isfile(savename)
            # We are loading the priors EVERY time, so there is no benefit to ordering
            # This is not optimal, but reduces scope confusion
            # These first two could be globals, but load them here for consistency
            # V_dib_noLSF_soft_lst = []
            # for dib in dib_waves
            #     local f = h5open(prior_dict["DIB_noLSF_soft_$(dib)"])
            #     push!(V_dib_noLSF_soft_lst, read(f["Vmat"]))
            #     close(f)
            # end
        
            f = h5open(prior_dict["starLines_refLSF"])
            V_starlines_refLSF = read(f["Vmat"])
            close(f)

            ### Need to load the priors here
            # f = h5open(prior_dict["skycont"]*lpad(adjfibindx,3,"0")*".h5")
            # V_skycont = read(f["Vmat"])
            # chebmsk_exp = convert.(Bool,read(f["chebmsk_exp"]))
            # close(f)

            # f = h5open(prior_dict["skyLines_bright"]*lpad(adjfibindx,3,"0")*".h5")
            # V_skyline_bright = read(f["Vmat"])
            # submsk_bright = convert.(Bool,read(f["submsk"]))
            # close(f)

            # f = h5open(prior_dict["skyLines_faint"]*lpad(adjfibindx,3,"0")*".h5")
            # V_skyline_faint = read(f["Vmat"])
            # submsk_faint = convert.(Bool,read(f["submsk"]))
            # close(f)

            # skymsk_bright = chebmsk_exp .& submsk_bright #
            # skymsk_faint = chebmsk_exp .& submsk_faint
            # # global skymsk = chebmsk_exp .& (submsk_bright .| submsk_faint)
            # skymsk = chebmsk_exp .& submsk_faint # completely masking all bright lines b/c detector response is nonlinear;

            f = h5open(prior_dict["chebmsk"])
            chebmsk_exp = read(f["chebmsk_exp"])
            close(f)
            skymsk_bright = chebmsk_exp #.& submsk_bright #
            skymsk_faint = chebmsk_exp #.& submsk_faint
            skymsk = chebmsk_exp #.& submsk_faint # completely masking all bright lines b/c detector response is nonlinear;

            # f = h5open(prior_dict["starCont"]*lpad(adjfibindx,3,"0")*".h5")
            f = h5open(prior_dict["starCont"])
            V_starcont = read(f["Vmat"])
            close(f)

            V_starlines = V_starlines_refLSF #hack
            # f = h5open(prior_dict["starLines_LSF"]*lpad(adjfibindx,3,"0")*".h5")
            # V_starlines = read(f["Vmat"])
            if ddstaronly
                V_starlines_refLSF = V_starlines
                msk_starCor = convert.(Bool,read(f["msk_starCor"]))
            else
                msk_starCor = ones(Bool,length(chebmsk_exp))
            end
            # close(f)

            # V_dib_lst = []
            # for dib in dib_waves
            #     local f = h5open(prior_dict["DIB_LSF_$(dib)"]*lpad(adjfibindx,3,"0")*".h5")
            #     push!(V_dib_lst,read(f["Vmat"]))
            #     close(f)
            # end

            # V_dib_soft_lst = []
            # for dib in dib_waves
            #     local f = h5open(prior_dict["DIB_LSF_soft_$(dib)"]*lpad(adjfibindx,3,"0")*".h5")
            #     push!(V_dib_soft_lst,read(f["Vmat"]))
            #     close(f)
            # end
            
            ### Single spectrum loop
            # prior_vec = (V_skycont,chebmsk_exp,V_skyline_bright,V_skyline_faint,skymsk_bright,skymsk_faint,skymsk,V_starcont,V_starlines_refLSF,V_starlines,msk_starCor,V_dib_lst, V_dib_soft_lst,V_dib_noLSF_soft_lst)
            prior_vec = (chebmsk_exp,V_skyline_bright,V_skyline_faint,skymsk_bright,skymsk_faint,skymsk,V_starcont,V_starlines_refLSF,V_starlines,msk_starCor)
            pipeline_single_spectra_bind(argtup) = pipeline_single_spectra(argtup, prior_vec; ddstaronly=ddstaronly)
            for (ind,indval) in enumerate(indsubset)
                push!(out,pipeline_single_spectra_bind(indval))
            end

            ### Save Exporting
            metai = 1
            RVind, RVchi, RVcom, strpo = 2, 3, 4, 5
            # DIBind, EWind, DIBchi, DIBcom = 6, 7, 8, 9
            # dibsavesz = 4


            framecnts, chipmidtimes, a_relFlux, b_relFlux, c_relFlux, cartVisit, ingest_bit
            
            ## RV Block
            RVextract = [
                # meta info
                (x->x[metai][1],                        "data_pix_cnt"),
                (x->x[metai][2],                        "starscale"),
                (x->x[metai][3],                        "skyscale0"),
                # (x->x[metai][4],                        "frame_counts"),
                # (x->x[metai][5],                        "chip_midtimes"),
                # (x->x[metai][6],                        "a_relFlux"),
                # (x->x[metai][7],                        "b_relFlux"),
                # (x->x[metai][8],                        "c_relFlux"),
                # (x->x[metai][9],                        "cartVisit"),
                # (x->x[metai][10],                       "ingestBit"),
                (x->x[metai][4],                        "flux"),
                (x->x[metai][5],                        "fluxivar"),
                (x->x[metai][6],                        "flux_nans"),
                (x->x[metai][7],                        "fluxerr2_nans"),
                (x->convert(Vector{Int},x[metai][8]),   "simplemsk"),
                (x->adjfibindx,                         "adjfiberindx"),

                (x->Float64.(x[RVind][1][1]),           "RV_pixoff_final"),
                (x->Float64.(x[RVind][1][3]),           "RV_pixoff_disc_final"),
                (x->x[RVind][1][2],                     "RV_minchi2_final"),
                (x->x[RVind][1][6],                     "RV_flag"),
                (x->x[RVind][1][7],                     "RV_pix_var"),
                                    
                (x->x[RVchi][1],                        "RVchi2_residuals"),
                # (x->x[RVchi][2],                        "RVchi2_residuals_flux_scaled"),
                (x->x[RVchi][2],                        "avg_flux_conservation"),
                (x->x[RVchi][3],                        "final_pix_cnt"),
                (x->x[RVchi][4],                        "starscale1"),
                (x->x[RVchi][5],                        "skyscale1"),
                                    
                (x->x[RVind][2][1][3],                  "RV_p5delchi2_lvl1"),
                (x->x[RVind][2][2][3],                  "RV_p5delchi2_lvl2"),
                (x->x[RVind][2][3][3],                  "RV_p5delchi2_lvl3"),

                (x->x[RVcom][1],                        "x_residuals_z_v0"),
                (x->x[RVcom][2],                        "x_residuals_v0"),
                # (x->x[RVcom][3],                        "x_skyLines_bright_v0"),
                (x->x[RVcom][3],                        "x_skyLines_faint_v0"),
                (x->x[RVcom][4],                        "x_skyContinuum_v0"),
                (x->x[RVcom][5],                        "x_starContinuum_v0"),
                # (x->x[RVcom][6],                        "x_starLineCor_v0"),
                (x->x[RVcom][6],                        "x_starLines_v0"),
                (x->x[RVcom][7],                        "x_starLineCof_v0"),
                (x->x[RVcom][8],                        "tot_p5chi2_v0"), 
                (x->x[RVcom][9],                        "apVisit_v0"),       
                (x->Int.(x[RVcom][10]),                 "finalmsk"),
                (x->x[RVcom][11],                       "x_starLines_restFrame_v0"),      
                                    
                (x->x[strpo],                           "x_starLines_err_v0"),    
            ]
            # ## DIB Block
            # DIBextract = []
            # for (dibindx,dibw) in enumerate(dib_center_lst)
            #     dib = string(round(Int,dibw))
            #     dibind = string(dibindx)
            #     push!(DIBextract,[
            #     # Further chi2 refinement does not have fixed sizing because can hit grid edge
            #     (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][1][1]),        "DIB_pixoff_final_$(dibind)_$(dib)"),
            #     (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][1][2]),        "DIB_sigval_final_$(dibind)_$(dib)"),
            #     (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][3][1]),        "DIB_pixoff_disc_final_$(dibind)_$(dib)"),
            #     (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][3][2]),        "DIB_sigval_disc_final_$(dibind)_$(dib)"),
            #     (x->x[DIBind+dibsavesz*(dibindx-1)][1][2],                     "DIB_minchi2_final_$(dibind)_$(dib)"),
            #     (x->x[DIBind+dibsavesz*(dibindx-1)][1][6],                     "DIB_flag_$(dibind)_$(dib)"),
            #     (x->[x[DIBind+dibsavesz*(dibindx-1)][1][7:11]...],             "DIB_hess_var_$(dibind)_$(dib)"),
                                    
            #     (x->x[DIBind+dibsavesz*(dibindx-1)][2][1][3],                  "DIB_p5delchi2_lvl1_$(dibind)_$(dib)"),
            #     (x->x[DIBind+dibsavesz*(dibindx-1)][2][2][3],                  "DIB_p5delchi2_lvl2_$(dibind)_$(dib)"),
            #     (x->x[DIBind+dibsavesz*(dibindx-1)][2][3][3],                  "DIB_p5delchi2_lvl3_$(dibind)_$(dib)"),

            #     (x->x[EWind+dibsavesz*(dibindx-1)][1],                         "EW_dib_$(dibind)_$(dib)"),
            #     (x->x[EWind+dibsavesz*(dibindx-1)][2],                         "EW_dib_err_$(dibind)_$(dib)"),
                                    
            #     (x->x[DIBchi+dibsavesz*(dibindx-1)][1],                        "DIBchi2_residuals_$(dibind)_$(dib)"),
            #     # (x->x[DIBchi+dibsavesz*(dibindx-1)][2],                        "DIBchi2_residuals_flux_scaled_$(dibind)_$(dib)"),

            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][1],                        "x_residuals_z_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][2],                        "x_residuals_v1_$(dibind)_$(dib)"),
            #     # (x->x[DIBcom+dibsavesz*(dibindx-1)][3],                        "x_skyLines_bright_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][3],                        "x_skyLines_faint_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][4],                        "x_skyContinuum_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][5],                        "x_starContinuum_v1_$(dibind)_$(dib)"),
            #     # (x->x[DIBcom+dibsavesz*(dibindx-1)][6],                        "x_starLineCor_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][6],                        "x_dib_flux_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][7],                        "x_starLines_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][8],                        "x_dib_v1_$(dibind)_$(dib)"),
            #     (x->x[DIBcom+dibsavesz*(dibindx-1)][9],                        "tot_p5chi2_v1_$(dibind)_$(dib)"),
            #     ])
            # end
            # extractlst = vcat(RVextract...,DIBextract...)
            extractlst = vcat(RVextract...)
                
            hdr_dict = Dict(   
                    "pipeline"=>"arMADGICS.jl",
                    "git_branch"=>git_branch,   
                    "git_commit"=>git_commit,
            )           
            h5write(savename,"hdr","This is only a header")
            h5writeattr(savename,"hdr",hdr_dict)        
            for elelst in extractlst
                extractor(out,elelst[1],elelst[2],savename)
            end
        end
        return 0
    end

    function extractor(x,elemap,elename,savename)
        len = length(x)
        exobj = elemap(x[1])
        outmat = zeros(eltype(exobj),size(exobj)...,len)
        for i=1:len
            flush(stdout)
            outmat[.. ,i] .= elemap(x[i])
        end
        h5write(savename,elename,outmat)
    end
end

# want to extend this to handle multiple tele and mjd
run_lst = get_telemjd_runlist_from_almanac(almanacFile, tele, mjd)

iterlst = []
Base.length(f::Iterators.Flatten) = sum(length, f.it)
for adjfibindx in 1:600
    subiter = filter(x->x[end].==adjfibindx,run_lst)
    subiterpart = Iterators.partition(subiter,batchsize)
    push!(iterlst,subiterpart)
end
ittot = Iterators.flatten(iterlst)
lenargs = length(ittot)
nwork = length(workers())
println("Batches to Do: $lenargs, number of workers: $nwork")
flush(stdout)

pout = @showprogress pmap(multi_spectra_batch,ittot)
# pout = @showprogress pmap(multi_spectra_batch,ittot,on_error=ex->2)
writedlm(out_dir*"pout_arMADGICS.txt",pout)
rmprocs(workers())

t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t0)); println("Total script runtime: $dt"); t_then = t_now; flush(stdout)