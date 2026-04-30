## This is script grabs a bunch of sky fiber spectra and decomposes them into continuum and line components, to serve as samples for building the sky prior
# Author - Andrew Saydjari

import Pkg;
using Dates;
t0 = now();
t_then = t0;
using InteractiveUtils;
versioninfo();
Pkg.update("ApogeeReduction");
Pkg.instantiate();
Pkg.precompile();
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Package activation took $dt"); t_then = t_now; flush(stdout);
using BLISBLAS
using Distributed, SlurmClusterManager, Suppressor, DataFrames

proj_path = dirname(Base.active_project()) * "/"
if "SLURM_NTASKS" in keys(ENV)
    using SlurmClusterManager
    addprocs(SlurmManager(), exeflags = ["--project=$proj_path"])
else
    addprocs(30, exeflags = ["--project=$proj_path"]) # change to a workers per node variable
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Worker allocation took $dt"); t_then = t_now; flush(stdout);
println("Running Main on ", gethostname()); flush(stdout);

@everywhere begin
    using BLISBLAS
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays
    using Interpolations, SparseArrays, ParallelDataTransfer, AstroTime, Suppressor
    using StatsBase, ProgressMeter
    using SortFilters, BasisFunctions, Random
end

@passobj 1 workers() proj_path
@everywhere begin
    prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
    src_dir = "$proj_path"
    include(src_dir*"src/utils.jl")
    include(src_dir*"src/gridSearch.jl")
    include(src_dir*"src/componentAndPosteriors.jl")
    include(src_dir*"src/fileNameHandling.jl")
    include(src_dir*"src/ingest.jl")
    include(src_dir*"src/lowRankPrescription.jl")
    include(src_dir*"src/marginalizeEW.jl")
    include(src_dir*"src/spectraInterpolation.jl")
    include(src_dir*"src/chi2Wrappers.jl")
    include(src_dir*"scripts/prior_build/prior_utils.jl")
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Worker loading took $dt"); flush(stdout);
println(BLAS.get_config()); flush(stdout);

using LibGit2;
println(proj_path)
git_branch, git_commit, git_clean = initalize_git(proj_path);
@passobj 1 workers() git_branch;
@passobj 1 workers() git_commit;
@passobj 1 workers() git_clean;

@everywhere begin
    runlist_range = 295:295 #1:600 #295, 245, 335, 101

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    prior_dict["skycont"] = prior_dir*"2026_04_26/star_priors/APOGEE_starcont_svd_60_f"
    prior_dict["StarContChipGapMsk"] = prior_dir*"2026_04_25/StarContChipGapMsk.h5"
end

@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
end

@everywhere begin
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
end

@everywhere begin
    # contScale is a tuning parameter we might want to investigate
    function get_sky_samples((adjfibindx, runlist);contscale=5e2,loc_parallel=false,seed=2023)
        ## ingest all sky spectra, save fluxes, ivars, pixel masks
        savefluxname = "sky_prior_disk/skyflux_" * lpad(adjfibindx, 3, "0") * ".h5"
        saveivarname = "sky_prior_disk/skyivar_" * lpad(adjfibindx, 3, "0") * ".h5"
        savemskname = "sky_prior_disk/skymsk_" * lpad(adjfibindx, 3, "0") * ".h5"
        if !isfile(savefluxname) || !isfile(saveivarname) || !isfile(savemskname)

            getSkyRough_partial(argtup) = getSkyRough(reduxBase, argtup.tele, argtup.mjd, argtup.expnum, almanacFile, fibindxoi=adjfibindx)
            pout = if loc_parallel 
                @showprogress pmap(getSkyRough_partial,runlist);
            else
                map(getSkyRough_partial,runlist);
            end

            mskpout = .!isnothing.(pout)
            count_sky = count(mskpout)
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
   

        ## fit the sky continuum
        savecontname = "sky_prior_disk/skycont_" * lpad(adjfibindx, 3, "0") * ".h5"
        savelinesname = "sky_prior_disk/skyline_" * lpad(adjfibindx, 3, "0") * ".h5"
        if !isfile(savecontname) || !isfile(savelinesname)
            chipgapmsk = h5open(prior_dict["StarContChipGapMsk"], "r") do f
                if adjfibindx > 300
                    read(f["lco"])
                else
                    read(f["apo"])
                end
            end
            f = h5open(prior_dict["starCont"]*lpad(adjfiberindx ,3,"0")*".h5")
            Vcontinuum = read(f["Vmat"])
            close(f)

            sky_smooth_wrapper_patial(indx) = sky_smooth_wrapper(indx,savefluxname,saveivarname,savemskname,chipgapmsk,V_starcont)

            pout = if loc_parallel 
                @showprogress pmap(sky_smooth_wrapper_patial,1:count_sky);
            else
                map(sky_smooth_wrapper_patial,1:count_sky);
            end
            skycont = zeros(length(wavetarg), count_sky)
            skylines  = zeros(Bool, length(wavetarg), count_sky)
            for (ind, subp) in enumerate(pout)
                skycont[:, ind] .= subp[1]
                skylines[:, ind] .= subp[2]
            end
            # Write results to HDF5 files
            h5open(savecontname, "w") do f
                write(f, "skycont", skyflux)
            end
            h5open(savelinesname, "w") do f
                write(f, "skyline", skyivar)
            end
            skyflux = nothing
            skyivar = nothing
        end
        close_tfun_files()
    end
end

# Collect all (tele, mjd) pairs first
f = h5open(almanacFile)
tele_mjd_pairs = []
if haskey(f, "raw")
    for tele in keys(f["raw"])
        for mjd in keys(f["raw"][tele])
            push!(tele_mjd_pairs, (tele, mjd))
        end
    end
else
    for tele in keys(f)
        for mjd in keys(f[tele])
            push!(tele_mjd_pairs, (tele, mjd))
        end
    end
end
close(f)
@everywhere get_telemjd_runlist_from_almanac_partial(argtup) = get_telemjd_runlist_from_almanac(almanacFile, argtup[1], argtup[2], fibertype=["sky"])
run_lsts = pmap(get_telemjd_runlist_from_almanac_partial, tele_mjd_pairs)
run_lst = vcat(run_lsts...)

# this takes like 5 min to run (speed up)
iterlist_full = []
for adjfiberindx in runlist_range
    subiter = filter(x -> x[:adjfiberindx] .== adjfiberindx, run_lst)
    push!(iterlist_full, subiter)
end

get_sky_samples(iterlist_full[1],loc_parallel=true)
# @showprogress pmap(get_sky_samples,iterlist_full)
