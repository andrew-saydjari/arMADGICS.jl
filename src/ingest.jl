## Ingest Module
# minor rewrite on the horizon

using AstroTime

function getSkyRough(reduxBase, tele, mjd, expnum; skyZcut=10, sky_obs_thresh=5)
    # hacks 
    almanacFile = get_almanac_file(reduxBase, mjd)
    f = h5open(almanacFile)
    fibtargDict = get_fibTargDict(f, tele, mjd, expnum)
    close(f)

    fibtypelist = map(x->fibtargDict[x],1:300)
    skyfibIDs = findall(fibtypelist.=="sky");
    skyfibIndxs = fiberID2fiberIndx.(skyfibIDs);

    #get ar1Dname
    ar1Dfname = get_1Duni_name(reduxBase, tele, mjd, expnum)

    # could speed up by only reading the columns we need
    f = jldopen(ar1Dfname) 
    skyspec = f["flux_1d"][:,skyfibIndxs]
    skyivar = f["ivar_1d"][:,skyfibIndxs];
    skymsk = f["mask_1d"][:,skyfibIndxs];
    close(f)

    skyScale = dropdims(nanzeromedian(skyspec,1),dims=1);
    skyMed = nanzeromedian(skyScale)
    skyIQR = nanzeroiqr(skyScale)
    skyZ = (skyScale.-skyMed)./skyIQR;
    mskSky = (abs.(skyZ).<skyZcut)

    msk_local_skyLines = dropdims(sum(.!isnanorzero.(skyspec[:,mskSky]),dims=2),dims=2).>sky_obs_thresh
    meanLocSkyLines = dropdims(nanzeromean(skyspec[:,mskSky],2),dims=2);
    VLocSkyLines = (skyspec[:,mskSky].-meanLocSkyLines)./sqrt(count(mskSky));
    meanLocSky = zero(meanLocSkyLines) # hack and ignores VLocSky
    return meanLocSky, meanLocSkyLines, VLocSkyLines, msk_local_skyLines
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

function get_telemjd_runlist_from_almanac(almanacFile, tele, mjd)
    teleind = (tele[1:3] == "lco") ? 2 : 1
    f = h5open(almanacFile)
    df_exp = read_almanac_exp_df(f, tele, mjd)
    msk_obj = (df_exp.exptype .== "OBJECT")
    row_exp = map(x->last(x,4),df_exp[msk_obj,:].exposure_str)
    run_lsts = []
    for expnum in row_exp
        expnumInt = parse(Int, expnum)
        fibtargDict = get_fibTargDict(f, tele, mjd, expnum)
        fibtypelist = map(x->fibtargDict[x],1:300)
        # should we be sky subtracting the sky fibers (seems like yes, but in Bayesian context?)
        targfibIDs = findall((fibtypelist.=="sci") .| (fibtypelist.=="tel"));
        targfibIndxs = fiberID2fiberIndx.(targfibIDs) .+ (teleind-1)*300;
        iterexp = Iterators.zip(Iterators.repeated(tele),Iterators.repeated(mjd),Iterators.repeated(expnumInt),targfibIndxs)
        push!(run_lsts, collect(iterexp))
    end
    run_lst = vcat(run_lsts...)
    close(f)
    return run_lst
end

function sky_decomp(outvec,outvar,simplemsk,V_skyline_bright,V_skyline_faint,V_skycont)   
    ## Select data for use (might want to handle mean more generally)
    Xd_obs = outvec[simplemsk];

    ## Set up residuals prior
    A = Diagonal(outvar[simplemsk]);
    Ainv = Diagonal(1 ./outvar[simplemsk]);

    ## Set up priors
    V_skyline_bright_c = V_skyline_bright
    V_skyline_bright_r = V_skyline_bright_c[simplemsk,:]
    V_skyline_faint_c = V_skyline_faint
    V_skyline_faint_r = V_skyline_faint_c[simplemsk,:]
    V_skycont_c = V_skycont
    V_skycont_r = V_skycont_c[simplemsk,:]
    
    # Compute sky line/continuum separation
    Vcomb = hcat(V_skyline_bright_r,V_skyline_faint_r,V_skycont_r);
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    x_comp_lst = deblend_components_all_asym(Ctotinv, Xd_obs, (V_skycont_r, ), (V_skycont_c, ))

    return x_comp_lst[1]
end

function stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberindx; telluric_div=false, cache_dir="../local_cache")

    plateFile = build_platepath(release_dir,redux_ver,tele,field,plate,mjd,"a")
    frame_lst = getFramesFromPlate(plateFile)
    
    #make a dictionary of values for chip a,b,c from fluxing fits
    thrptDict = Dict{String,Float64}()
    f = FITS(cache_fluxname(tele,field,plate,mjd; cache_dir=cache_dir))
    for chip in ["a","b","c"]
        thrpt = read(f[chip],fiberindx)
        thrptDict[chip] = thrpt
    end
    cartVisit = parse(Int,read_header(f[1])["CARTID"])
    close(f)

    ingest_bit = 0
    fill!(outvec,0)
    fill!(outvar,0)
    fill!(cntvec,0)
    if telluric_div
        fill!(telvec,0)
    end
    time_lsts = [[],[],[]]
    for imid in frame_lst
        fill!(Xd_stack,0)
        fill!(Xd_std_stack,0)
        fill!(waveobs_stack,0)
        fill!(pixmsk_stack,0)
        if telluric_div
            fill!(telluric_stack,0)
        end
        fill!(fullBit,0)
        for (chipind,chip) in enumerate(["c","b","a"]) #needs to be c,b,a for chip ind to be right
            fname = build_framepath(release_dir,redux_ver,tele,mjd,imid,chip)
            f = FITS(fname)
            hdr = read_header(f[1])
            midtime = modified_julian(TAIEpoch(hdr["DATE-OBS"]))+(hdr["EXPTIME"]/2/3600/24)days #TAI or UTC?
            push!(time_lsts[chipind],AstroTime.value(midtime))
            Xd = read(f[2],:,fiberindx)
            Xd_stack[(1:2048).+(chipind-1)*2048] .= Xd[end:-1:1]; # ./thrptDict[chip] has already been done at the 2D -> 1D level
            Xd_std = read(f[3],:,fiberindx)
            Xd_std_stack[(1:2048).+(chipind-1)*2048] .= Xd_std[end:-1:1].*err_factor.(Xd[end:-1:1],Ref(err_correct_Dict[join([tele[1:6],chip],"_")])); # ./thrptDict[chip] has already been done at the 2D -> 1D level
            pixmsk = read(f[4],:,fiberindx);
            pixmsk_stack[(1:2048).+(chipind-1)*2048] .= pixmsk[end:-1:1]
            waveobsa = read(f[5],:,fiberindx);
            waveobs_stack[(1:2048).+(chipind-1)*2048] .= waveobsa[end:-1:1]
            fullBit[(1:2048).+(chipind-1)*2048] .+= 2^chipind
            close(f)
            if telluric_div
                vpath = build_visitpath(release_dir,redux_ver,tele,field,plate,mjd,fiberindx)
                cpath = visit2cframe(vpath,tele,imid,chip)
                f = FITS(cpath)
                telluric = read(f[8]);
                tellmsk = dropdims(sum(telluric,dims=1).!=0,dims=1)
                tellindx = find_nearest_nz(tellmsk,fiberindx)
                telluric_stack[(1:2048).+(chipind-1)*2048] .= telluric[end:-1:1,tellindx]
                close(f)
            end
        end
        fullBit[((pixmsk_stack .& 2^0).!=0)] .+= 2^4 # call pixmask bit 0 bad
        fullBit[fullBit.==0] .+= 2^4 # call chip gaps bad for alt space

        goodpix = ((pixmsk_stack .& 2^0).==0) .& ((fullBit .& 2^4).==0) .& (.!isnan.(Xd_std_stack)) .& (Xd_std_stack.< (10^10))
        if telluric_div
            Xd_stack./= telluric_stack
            Xd_std_stack./= telluric_stack
        end

        obsBit = fullBit[goodpix]
        Xd_obs = Xd_stack[goodpix]
        Xd_std_obs = Xd_std_stack[goodpix];
        waveobs = waveobs_stack[goodpix];
        pixindx = (1:length(waveobs_stack))[goodpix]

        Rinv = generateInterpMatrix_sparse_inv(waveobs,obsBit,wavetarg,pixindx)
        normvec = dropdims(sum(Rinv,dims=2),dims=2)
        msk_inter = (normvec.!=0)

        fullvec = Rinv*Xd_obs
        fullvec[.!msk_inter] .= 0
        varvec =  (Rinv.^2)*(Xd_std_obs.^2)
        varvec[.!msk_inter] .= 0

        outvec .+= fullvec
        outvar .+= varvec
        cntvec .+= msk_inter

        if telluric_div
            Rinv = generateInterpMatrix_sparse_inv(waveobs_stack,ones(Int,length(fullBit)).*2^3,wavetarg,(1:length(waveobs_stack)));
            telvec .+= Rinv*telluric_stack
        end

        if all(isnanorzero.(Xd_stack)) && ((ingest_bit & 2^1)==0)
            ingest_bit += 2^1 # ap1D flux is literally NaNs (for at least one of the exposures)
        elseif all(.!((fullBit .& 2^4).==0)) && ((ingest_bit & 2^2)==0)
            ingest_bit += 2^2 # all masked by DRP pixmask  (for at least one of the exposures)
        elseif all(isnanorzero.(Xd_std_stack)) && ((ingest_bit & 2^3)==0)
            ingest_bit += 2^3 # either upstream std NaNs or err_factor NaNed  (for at least one of the exposures)
        end
    end
    framecnts = maximum(cntvec) # a little shocked that I throw it away if it is bad in even one frame
    outvec./=framecnts
    outvar./=(framecnts^2)
    if telluric_div
        telvec./=framecnts
    end

    if all(isnanorzero.(outvec))
        ingest_bit += 2^4 # all NaNs or zeros after interp
    elseif (thrptDict["a"]<0) || (thrptDict["b"]<0) || (thrptDict["c"]<0)
        ingest_bit += 2^5 # bad thrpt below thrpt_cut, NaNed by arMADGICS.jl
    elseif isnan(thrptDict["a"]) || isnan(thrptDict["b"]) || isnan(thrptDict["c"])
        ingest_bit += 2^6 # NaNs in apFlux file, however arMADGICS.jl does not depend on these values
    end

    if (thrptDict["a"]<0) || (thrptDict["b"]<0) || (thrptDict["c"]<0)
        outvec.*=NaN
    end
    
    simplemsk = (cntvec.==framecnts)
    starscale = nanzeromedian(outvec[simplemsk])

    goodframeIndx = length.(time_lsts).!=0
    chipmidtimes = zeros(3)
    chipmidtimes[goodframeIndx] .= mean.(time_lsts[goodframeIndx]) #consider making this flux weighted (need to worry about skyline variance driving it)
    chipmidtimes[.!goodframeIndx] .= NaN
    metaexport = (starscale,framecnts,thrptDict["a"],thrptDict["b"],thrptDict["c"],cartVisit,ingest_bit)
    if telluric_div
        return outvec, outvar, cntvec, chipmidtimes, metaexport, telvec
    end
    return outvec, outvar, cntvec, chipmidtimes, metaexport
end