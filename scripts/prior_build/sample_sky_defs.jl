## Core definitions for sample_sky.jl (M4): grab sky-fiber spectra and decompose them
## into continuum and line components, to serve as samples for building the sky prior.
#
# This file holds the script's main body as callable functions (T2.6 smoke-test pattern)
# so that both the production driver (sample_sky.jl, which includes it @everywhere) and
# the smoke test (test/prior_build_smoke.jl) can run it on real inputs.
#
# Requires (from the includer): HDF5, JLD2, LinearAlgebra, LowRankOps, SortFilters,
# Distributed (for the loc_parallel/runlist_parallel paths), plus
# src/utils.jl, src/fileNameHandling.jl, src/ingest.jl, src/lowRankPrescription.jl,
# src/componentAndPosteriors.jl, and scripts/prior_build/prior_utils.jl.
# Author - Andrew Saydjari

wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat

const _TFUN_FILES = Dict{String,HDF5.File}()

function get_tfun_file(path::String)
    if !haskey(_TFUN_FILES, path) || !isopen(_TFUN_FILES[path])
        _TFUN_FILES[path] = h5open(path, "r")
    end
    return _TFUN_FILES[path]
end

function close_tfun_files()
    for f in values(_TFUN_FILES)
        if isopen(f)
            close(f)
        end
    end
    empty!(_TFUN_FILES)
    return nothing
end

function sky_smooth_wrapper(indx,fluxp,ivarp,skymskp,chipgapmsk,Vcontinuum)
    fluxf = get_tfun_file(fluxp)
    ivarf = get_tfun_file(ivarp)
    skymskf = get_tfun_file(skymskp)

    fvec = fluxf["skyflux"][:,indx]
    fivar = ivarf["skyivar"][:,indx]
    skymsk = skymskf["skymsk"][:,indx]

    simplemsk = skymsk .& chipgapmsk;
    fcont, flines = sky_smooth_fit(fvec,fivar,simplemsk,Vcontinuum)
    return fcont, flines
end

function sky_smooth_fit(fvec,fivar,simplemsk,Vcontinuum)
    ## Select data for use (might want to handle mean more generally)
    wave_obs = wavetarg[simplemsk]
    Xd_obs0 = fvec[simplemsk];

    tmsk = ret_qlines(Xd_obs0,wave_obs)
    Xd_obs = fvec[simplemsk][tmsk];

    ## Set up residuals prior
    Ainv = Diagonal(fivar[simplemsk][tmsk]);

    ## Set up priors
    V_smooth_c = Vcontinuum
    V_smooth_r = Vcontinuum[simplemsk,:][tmsk,:]

    # Compute sky line/continuum separation
    Vcomb = V_smooth_r
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    x_comp_lst = deblend_components_all_asym(Ctotinv, Xd_obs, (V_smooth_r, ), (V_smooth_c, ))

    fnew = fvec.-x_comp_lst[1]
    fnew[.!simplemsk].=0
    return x_comp_lst[1], fnew
end

# contScale is a tuning parameter we might want to investigate
function get_sky_samples((adjfibindx, runlist);
        reduxBase, almanacFile, prior_dict,
        out_dir="sky_prior_disk/", contscale=5e2, loc_parallel=false, seed=2023)
    mkpath(out_dir)
    ## ingest all sky spectra, save fluxes, ivars, pixel masks
    savefluxname = joinpath(out_dir, "skyflux_" * lpad(adjfibindx, 3, "0") * ".h5")
    saveivarname = joinpath(out_dir, "skyivar_" * lpad(adjfibindx, 3, "0") * ".h5")
    savemskname = joinpath(out_dir, "skymsk_" * lpad(adjfibindx, 3, "0") * ".h5")
    if !isfile(savefluxname) || !isfile(saveivarname) || !isfile(savemskname)

        # M4: getSkyRough selects sky fibers in the NATIVE (1:300) fiber index space of the
        # ar1Dunical file; passing the adjusted index (301-600 for LCO) made LCO sky sampling
        # silently empty. Convert before passing. The M-SKY guard inside getSkyRough then
        # validates the fiber (returns nothing when it is excluded), so no extra checks here.
        fibindxoi = adjfiberindx2fiberindx(adjfibindx)
        getSkyRough_partial(argtup) = getSkyRough(reduxBase, argtup.tele, argtup.mjd, argtup.expnum, almanacFile, fibindxoi=fibindxoi)
        pout = if loc_parallel
            @showprogress pmap(getSkyRough_partial,runlist);
        else
            map(getSkyRough_partial,runlist);
        end

        mskpout = .!isnothing.(pout)
        count_sky = count(mskpout)
        # M4: this failure mode (all exposures returning nothing) used to write empty
        # sample files silently; fail loudly instead.
        (count_sky > 0) || error("get_sky_samples: 0 usable sky spectra for adjfibindx=$adjfibindx out of $(length(runlist)) runlist entries (tele=$(isempty(runlist) ? "?" : first(runlist).tele)); check the fiber index convention and the sky-fiber guard log.")
        skyflux = zeros(length(wavetarg), count_sky)
        skyivar = zeros(length(wavetarg), count_sky)
        skymsk  = zeros(Bool, length(wavetarg), count_sky)
        for (ind, pindx) in enumerate(findall(mskpout))
            skyflux[:, ind] .= pout[pindx][1]
            skyivar[:, ind] .= pout[pindx][2]
            skymsk[:, ind] .= pout[pindx][3]
        end

        # Write results to HDF5 files
        h5open(savefluxname, "w") do f
            write(f, "skyflux", skyflux)
        end
        h5open(saveivarname, "w") do f
            write(f, "skyivar", skyivar)
        end
        h5open(savemskname, "w") do f
            write(f, "skymsk", skymsk)
        end
        skyflux = nothing
        skyivar = nothing
        skymsk = nothing
    end

    # M4: recover the sample count from the checkpoint file so a restart (flux files
    # exist, cont/lines files do not) works; count_sky was previously undefined here.
    count_sky = h5open(savefluxname, "r") do f
        size(f["skyflux"], 2)
    end
    (count_sky > 0) || error("get_sky_samples: existing checkpoint $savefluxname holds 0 sky samples for adjfibindx=$adjfibindx; delete it and rerun (pre-fix LCO runs wrote empty files).")

    ## fit the sky continuum
    savecontname = joinpath(out_dir, "skycont_" * lpad(adjfibindx, 3, "0") * ".h5")
    savelinesname = joinpath(out_dir, "skyline_" * lpad(adjfibindx, 3, "0") * ".h5")
    if !isfile(savecontname) || !isfile(savelinesname)
        chipgapmsk = h5open(prior_dict["StarContChipGapMsk"], "r") do f
            if adjfibindx > 300
                read(f["lco"])
            else
                read(f["apo"])
            end
        end
        f = h5open(prior_dict["starCont"]*lpad(adjfibindx,3,"0")*".h5")
        Vcontinuum = read(f["Vmat"])
        close(f)

        sky_smooth_wrapper_partial(indx) = sky_smooth_wrapper(indx,savefluxname,saveivarname,savemskname,chipgapmsk,Vcontinuum)

        pout = if loc_parallel
            @showprogress pmap(sky_smooth_wrapper_partial,1:count_sky);
        else
            map(sky_smooth_wrapper_partial,1:count_sky);
        end
        skycont = zeros(length(wavetarg), count_sky)
        skylines = zeros(length(wavetarg), count_sky)
        for (ind, subp) in enumerate(pout)
            skycont[:, ind] .= subp[1]
            skylines[:, ind] .= subp[2]
        end
        # Write results to HDF5 files (M4: previously wrote the stale skyflux/skyivar
        # bindings, i.e. `nothing`, instead of the decomposition just computed)
        h5open(savecontname, "w") do f
            write(f, "skycont", skycont)
        end
        h5open(savelinesname, "w") do f
            write(f, "skyline", skylines)
        end
        skycont = nothing
        skylines = nothing
    end
    close_tfun_files()
end

function collect_tele_mjd_pairs(almanacFile)
    f = h5open(almanacFile)
    tele_mjd_pairs = Tuple{String,String}[]
    root = haskey(f, "raw") ? f["raw"] : f
    for tele in keys(root)
        for mjd in keys(root[tele])
            push!(tele_mjd_pairs, (tele, mjd))
        end
    end
    close(f)
    return tele_mjd_pairs
end

"""
    sample_sky_main(reduxBase, almanacFile, runlist_range; prior_dict, out_dir,
                    tele_mjd_pairs=nothing, loc_parallel=false, runlist_parallel=false)

Main body of sample_sky.jl: build the sky-fiber runlist from the almanac (all
(tele, mjd) pairs unless `tele_mjd_pairs` is given) and generate sky samples for
every adjusted fiber index in `runlist_range`. `loc_parallel` pmaps within a fiber,
`runlist_parallel` pmaps the almanac runlist construction. Returns the per-fiber
iteration list `[(adjfibindx, runlist), ...]`.
"""
function sample_sky_main(reduxBase, almanacFile, runlist_range;
        prior_dict, out_dir="sky_prior_disk/",
        tele_mjd_pairs=nothing, loc_parallel=false, runlist_parallel=false)
    # Collect all (tele, mjd) pairs first
    if isnothing(tele_mjd_pairs)
        tele_mjd_pairs = collect_tele_mjd_pairs(almanacFile)
    end
    # M4: kwarg is accepted_fibtypes (was passed as fibertype= -> MethodError)
    get_runlist_partial(argtup) = get_telemjd_runlist_from_almanac(almanacFile, argtup[1], argtup[2], accepted_fibtypes=["sky"])
    run_lsts = if runlist_parallel
        pmap(get_runlist_partial, tele_mjd_pairs)
    else
        map(get_runlist_partial, tele_mjd_pairs)
    end
    run_lst = vcat(run_lsts...)

    # this takes like 5 min to run (speed up)
    # M4: entries must be (adjfibindx, runlist) tuples to match get_sky_samples
    iterlist_full = []
    for adjfibindx in runlist_range
        subiter = filter(x -> x.adjfiberindx == adjfibindx, run_lst)
        push!(iterlist_full, (adjfibindx, subiter))
    end

    for itertup in iterlist_full
        get_sky_samples(itertup; reduxBase=reduxBase, almanacFile=almanacFile,
            prior_dict=prior_dict, out_dir=out_dir, loc_parallel=loc_parallel)
    end
    return iterlist_full
end
