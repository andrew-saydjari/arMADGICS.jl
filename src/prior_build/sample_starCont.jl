## This is script grabs a bunch of transfer functions and telluric components fit to domeflats for building the starContinuum prior. (For apMADGICS we used HOT STD stars instead)
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
using Distributed, ArgParse, SlurmClusterManager, Suppressor, DataFrames, DelimitedFiles
using ApogeeReduction: read_almanac_exp_df, get_fibTargDict, check_type_for_jld2

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
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays, JLD2, FileIO
    using Interpolations, SparseArrays, AstroTime, Suppressor
    using ApogeeReduction, DataFrames, ParallelDataTransfer
    using StatsBase, ProgressMeter
    using SortFilters, BasisFunctions, Random, DustExtinction, DelimitedFiles
    using ApogeeReduction: check_type_for_jld2, adjFiberIndx2FiberIndx
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Worker loading took $dt"); flush(stdout);
@passobj 1 workers() proj_path
println(BLAS.get_config()); flush(stdout);

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
    include(src_dir*"src/prior_build/prior_utils.jl")
end

using LibGit2;
println(proj_path)
git_branch, git_commit, git_clean = initalize_git(proj_path);
@passobj 1 workers() git_branch;
@passobj 1 workers() git_commit;
@passobj 1 workers() git_clean;

@everywhere begin
    runlist_range = 295 #1:600 #295, 245, 335, 101

    nsamp = 10_000

    Teff_rng = 4_000:1:10_000
    Av_rng = 0:1e-4:5
    Rv_rng = 2.6:1e-4:3.6

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # Data for Cals (not really a prior, but an input the results depend on in detail)
    prior_dict["LSF_mat_APO"] = prior_dir*"2026_04_25/mat_lsf_out/sp_combolsfmat_norm_6_" # last made 2023_04_01 by AKS
    prior_dict["LSF_mat_LCO"] = prior_dir*"2026_04_26/mat_lsf_out/sp_combolsfmat_norm_6_" # last made 2023_04_07 by AKS
    prior_dict["fracTellSamples_APO"] = prior_dir*"2026_04_25/outsamptell_apo.jdat" # last made 2023_04_03 by AKS
    prior_dict["fracTellSamples_LCO"] = prior_dir*"2026_04_26/outsamptell_lco.jdat" # last made 2023_04_07 by AKS

    # Location of the Tfun samples
    tell_base = "/mnt/home/acasey/scratch/20260220-arjl-domeflats/"
    prior_dict["tfun_samples_APO"] = tell_base*"20260323_apo.h5"
    prior_dict["tfun_samples_LCO"] = tell_base*"20260323_lco.h5"
    prior_dict["tfun_sample_lst_APO"] = prior_dir*"2026_04_25/20260323_apo_tfunlist.jdat"
    prior_dict["tfun_sample_lst_LCO"] = prior_dir*"2026_04_25/20260323_lco_tfunlist.jdat"
end

@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg);
    x_model = 15000:0.01:17000
end

@everywhere begin
    function genModSamp(intup,tellFracSamples,TfunSamplef,Ksp,nvecLSF,Atell,fiberindx)
        Teff,Av,Rv,Tfunindx,Tfracindx = intup
        bbs = blackbody.(Ref(Teff), x_model*1e-8);
        rvec = redden_mult(x_model,Av,Rv);
        Tfunsample = exp.(Atell*TfunSamplef["theta"][:,fiberindx,Tfunindx])
        return tellFracSamples[:,Tfracindx].*Tfunsample./nanzeromedian(Tfunsample).*((Ksp*(rvec.*bbs))./nvecLSF);
    end
end

@everywhere begin
    function gen_starCont_samples(adjfibindx;loc_parallel=false,seed=203)
        fiberindx = adjFiberIndx2FiberIndx(adjfibindx)

        savename = "tell_prior_disk/starCont_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            Ksp = if adjfibindx>300
                deserialize(prior_dict["LSF_mat_LCO"]*lpad(adjfibindx-300,3,"0")*".jdat");
            else
                deserialize(prior_dict["LSF_mat_APO"]*lpad(adjfibindx,3,"0")*".jdat");
            end

            nvecLSF = dropdims(sum(Ksp,dims=2),dims=2); # used only in starCont sample gen

            # Get tell obs to use from disk
            TfunSamplef =  if adjfibindx > 300
                h5open(prior_dict["tfun_samples_LCO"], "r")
            else
                h5open(prior_dict["tfun_samples_APO"], "r")
            end
            Atell = permutedims(read(TfunSamplef["design_matrix"]),[2,1])
            Tfungoodindxlist = if adjfiberindx > 300
                deserialize(prior_dict["tfun_sample_lst_LCO"]);
            else
                deserialize(prior_dict["tfun_sample_lst_APO"]);
            end
            # Get list for which obs are "useable"


            # Load the fraction telluric model samples (10k random from visits, stacked frames)
            tellFracSamples = if adjfibindx > 300
                deserialize(prior_dict["fracTellSamples_LCO"])
            else
                deserialize(prior_dict["fracTellSamples_APO"])
            end

            # Generate StarCont Samples
            rng = MersenneTwister(seed)
            # draw over parameter space (consider marginalizing over the dust exponent)
            Teff_lst = rand(rng,Teff_rng,nsamp)
            Av_lst = rand(rng,Av_rng,nsamp)
            Rv_lst = rand(rng,Rv_rng,nsamp)
            Tfunindx_lst = rand(rng,Tfungoodindxlist[fiberindx],nsamp);
            Tfracindx_lst = rand(rng,1:size(tellFracSamples,2),nsamp);
            itobj = Iterators.zip(Teff_lst,Av_lst,Rv_lst,Tfunindx_lst,Tfracindx_lst)

            genModSamp_bound(itobj) = genModSamp(itobj,tellFracSamples,TfunSamplef,Ksp,nvecLSF,Atell,fiberindx)
            pout = if loc_parallel
                @showprogress pmap(genModSamp_bound,itobj); #not very fast because of passing
            else
                map(genModSamp_bound,itobj);
            end
            outsamp = zeros(length(wavetarg),size(pout,1));
            for i=1:size(pout,1)
                outsamp[:,i].=pout[i]
            end
            serialize(savename,outsamp)
        end
    end
end

gen_starCont_samples(runlist_range,loc_parallel=true)
# @showprogress pmap(gen_starCont_samples,1:600)