## This is script grabs a bunch of transfer functions and telluric components fit to domeflats for building the starContinuum prior. (For apMADGICS we used HOT STD stars instead)
# Author - Andrew Saydjari

import Pkg;
using Dates;
t0 = now();
t_then = t0;
using InteractiveUtils;
versioninfo();
# ARM_SKIP_PKG_OPS=1 skips update/instantiate/precompile for reproducible reruns
# against an already-resolved Manifest (mirrors build_starCont.jl; E4b: Pkg.update
# would move the ApogeeReduction repo-rev=main pin and break bit-comparability
# with the 2026_09_03 E4 corpus).
if get(ENV, "ARM_SKIP_PKG_OPS", "0") != "1"
    Pkg.update("ApogeeReduction");
    Pkg.instantiate();
    Pkg.precompile();
end
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
    # Worker count override (ARM_STARCONT_NWORKERS); default 20 = dedicated-node
    # (ccalin051) etiquette cap. Slurm path above is unaffected.
    addprocs(parse(Int, get(ENV, "ARM_STARCONT_NWORKERS", "20")),
        exeflags = ["--project=$proj_path"])
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
    using ApogeeReduction: check_type_for_jld2, adjFiberIndx2FiberIndx, get_lsf_matrix
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
    runlist_range = 295 #1:600 #295, 245, 335, 101

    nsamp = 10_000

    Teff_rng = 4_000:1:10_000
    Av_rng = 0:1e-4:5
    Rv_rng = 2.6:1e-4:3.6

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # Secured prior inputs (plan v2 P1/E2; provenance + sha256 in prior_inputs/PROVENANCE.md)
    prior_inputs_dir = prior_dir*"2026_08_31/prior_inputs/"

    # Data for Cals (not really a prior, but an input the results depend on in detail)
    prior_dict["LSF_path_APO"] = prior_inputs_dir*"lsf_20260427/fpiLSFparams_REGULARIZED_apo_60861.h5"
    prior_dict["LSF_path_LCO"] = prior_inputs_dir*"lsf_20260427/fpiLSFparams_REGULARIZED_lco_60861.h5"
    prior_dict["fracTellSamples_APO"] = prior_dir*"2026_04_25/outsamptell_apo.jdat" # last made 2023_04_03 by AKS
    prior_dict["fracTellSamples_LCO"] = prior_dir*"2026_04_26/outsamptell_lco.jdat" # last made 2023_04_07 by AKS

    # Tfun samples: E3 bug-fixed telluric transfer-function refit products
    # (T_out=T_init rerun bug fixed; see 2026_09_02/telluric_refit_full/). Override the
    # directory without editing this file via ARM_TFUN_BASE (expects
    # <tell_base>/tellurics_refit_20260902_{apo,lco}.h5) or each file via
    # ARM_TFUN_{APO,LCO}. The pre-E3 delivered products (secured copies of acasey's
    # 20260220-arjl-domeflats 20260323_{apo,lco}.h5) live in
    # <prior_inputs_dir>tellurics_20260220_arjl_domeflats/ if a rollback is needed.
    tell_base = get(ENV, "ARM_TFUN_BASE",
        prior_dir*"2026_09_02/telluric_refit_full/out/")
    prior_dict["tfun_samples_APO"] = get(ENV, "ARM_TFUN_APO", tell_base*"tellurics_refit_20260902_apo.h5")
    prior_dict["tfun_samples_LCO"] = get(ENV, "ARM_TFUN_LCO", tell_base*"tellurics_refit_20260902_lco.h5")
    # Usable-(fiber,exposure) lists REBUILT against the E3 artifacts (row indices differ
    # from the delivered files!) with the 2026-09-03 consumption cuts: medflux>400
    # (2026_04_25-era convention), APO medflux<=10,000 bright cut, chi_sq_fiber<=p99.9
    # per telescope. Build script + BUILD_REPORT.txt sit next to the lists.
    # Override via ARM_TFUN_LIST_BASE or ARM_TFUN_LIST_{APO,LCO}. The old delivered-file
    # lists (2026_04_25/20260323_{apo,lco}_tfunlist.jdat) do NOT apply to the E3 files.
    tfun_list_base = get(ENV, "ARM_TFUN_LIST_BASE",
        prior_dir*"2026_09_03/tfunlists_refit20260902/")
    prior_dict["tfun_sample_lst_APO"] = get(ENV, "ARM_TFUN_LIST_APO", tfun_list_base*"20260902_apo_tfunlist.jdat")
    prior_dict["tfun_sample_lst_LCO"] = get(ENV, "ARM_TFUN_LIST_LCO", tfun_list_base*"20260902_lco_tfunlist.jdat")
end

@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg);
    x_model = 15000:0.01:17000
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

    function genModSamp(intup,tellFracSamples,tfun_path,Ksp,Atell,fiberindx)
        Teff,Av,Rv,Tfunindx,Tfracindx = intup
        bbs = blackbody.(Ref(Teff), x_model*1e-8);
        rvec = redden_mult(x_model,Av,Rv);
        TfunSamplef = get_tfun_file(tfun_path)
        Tfunsample = exp.(Atell*TfunSamplef["theta"][:,fiberindx,Tfunindx])
        return tellFracSamples[:,Tfracindx].*Tfunsample./nanzeromedian(Tfunsample).*(Ksp*(rvec.*bbs));
    end
end

@everywhere begin
    function gen_starCont_samples(adjfibindx;loc_parallel=false,seed=203)
        fiberindx = adjFiberIndx2FiberIndx(adjfibindx)

        savename = "tell_prior_disk/starCont_"*lpad(adjfibindx,3,"0")*".jdat"
        mkpath(dirname(savename))
        if !isfile(savename)
            Ksp = if adjfibindx>300
                get_lsf_matrix(adjfibindx, prior_dict["LSF_path_LCO"]);
            else
                get_lsf_matrix(adjfibindx, prior_dict["LSF_path_APO"]);
            end

            # Get tell obs to use from disk
            tfun_path =  if adjfibindx > 300
                prior_dict["tfun_samples_LCO"]
            else
                prior_dict["tfun_samples_APO"]
            end
            Atell = h5open(tfun_path, "r") do f
                permutedims(read(f["design_matrix"]), [2,1])
            end
            Tfungoodindxlist = if adjfibindx > 300
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

            genModSamp_bound(itobj) = genModSamp(itobj,tellFracSamples,tfun_path,Ksp,Atell,fiberindx)
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
            close_tfun_files()
        end
    end
end

# gen_starCont_samples(runlist_range,loc_parallel=true)
# ARM_STARCONT_RANGE restricts the adjusted fibers (default 1:600). Accepts
# comma-separated tokens, each "a:b" or a single index (e.g. "388,448,459,519").
# E4b used 301:600; the intelligent-policy delta run resamples the 4 changed fibers.
run_range = let s = get(ENV, "ARM_STARCONT_RANGE", "1:600")
    sort(unique(reduce(vcat, [occursin(":", tok) ?
        collect(parse(Int, split(tok, ":")[1]):parse(Int, split(tok, ":")[2])) :
        [parse(Int, tok)] for tok in split(s, ",")])))
end
println("sampling adjusted-fiber range: $run_range"); flush(stdout)
@showprogress pmap(gen_starCont_samples, run_range)