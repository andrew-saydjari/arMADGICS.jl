## This is script grabs uses the samples from sample_starCont.jl to build a covariance matrix prior for the star continuum
# Author - Andrew Saydjari,

import Pkg;
using Dates;
t0 = now();
t_then = t0;
using InteractiveUtils;
versioninfo();
# ARM_SKIP_PKG_OPS: skip the Pkg mutation/instantiate step (e.g. local E6-style
# regression builds against an already-instantiated env). Default: run them.
if !haskey(ENV, "ARM_SKIP_PKG_OPS")
    Pkg.update("ApogeeReduction");
    Pkg.instantiate();
    Pkg.precompile();
end
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
    # ARM_STARCONT_NWORKERS: local worker count (default 30, historical)
    addprocs(parse(Int, get(ENV, "ARM_STARCONT_NWORKERS", "30")), exeflags = ["--project=$proj_path"])
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
    using SortFilters, BasisFunctions, Random, DustExtinction
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
println(proj_path); flush(stdout)
git_branch, git_commit, git_clean = initalize_git(proj_path);
@passobj 1 workers() git_branch;
@passobj 1 workers() git_commit;
@passobj 1 workers() git_clean;

@everywhere begin
    runlist_range = 295 #1:600 #295, 245, 335, 101

    nsub = 60

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # StarCont Samples. Env overrides (all optional, defaults = historical paths):
    #   ARM_STARCONT_SAMPLES    - sample-file prefix (…/starCont_), e.g. the E4
    #                             2026_09_03 tell_prior_disk generation
    #   ARM_STARCONT_CHIPGAPMSK - chip-gap mask h5 (keys apo/lco)
    #   ARM_STARCONT_OUTDIR     - output dir for the built per-fiber h5 priors
    prior_dict["starcont"] = get(ENV, "ARM_STARCONT_SAMPLES",
        prior_dir*"2026_04_26/tell_prior_disk/starCont_")
    prior_dict["StarContChipGapMsk"] = get(ENV, "ARM_STARCONT_CHIPGAPMSK",
        prior_dir*"2026_04_25/StarContChipGapMsk.h5")
    out_dir = get(ENV, "ARM_STARCONT_OUTDIR", "star_priors")

    # ARM_STARCONT_FIBERS: comma list and/or a:b ranges (e.g. "101,245,295,335,350").
    # Default 1:600 (historical full build).
    build_fibers = if haskey(ENV, "ARM_STARCONT_FIBERS")
        sort(unique(reduce(vcat, [occursin(":", tok) ?
            collect(parse(Int, split(tok, ":")[1]):parse(Int, split(tok, ":")[2])) :
            [parse(Int, tok)] for tok in split(ENV["ARM_STARCONT_FIBERS"], ",")])))
    else
        collect(1:600)
    end
end

@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg);
end

@everywhere begin    
    function build_starCont(adjfibindx)
        savename = joinpath(out_dir, "APOGEE_starcont_svd_"*string(nsub)*"_f"*lpad(adjfibindx,3,"0")*".h5")
        mkpath(dirname(savename))
        if !isfile(savename)
            starcontfname = prior_dict["starcont"]*lpad(adjfibindx,3,"0")*".jdat"
            starcont = deserialize(starcontfname)

            # This should probably go back to being fiber dependent like I had in apMADGICS
            # chbmskfname = prior_dict["chebmsk_exp"]*lpad(adjfibindx,3,"0")*".jdat"
            # chebmsk_exp = deserialize(chbmskfname);
            chipgapmsk = h5open(prior_dict["StarContChipGapMsk"], "r") do f
                if adjfibindx > 300
                    read(f["lco"])
                else
                    read(f["apo"])
                end
            end

            specsum = dropdims(sum(starcont,dims=1),dims=1)
            Vred = starcont[chipgapmsk,specsum.>0];
            mnorm = mean(filter(.!iszero,Vred))
            Vred./=mnorm
            # weights = ones(size(Vred,2));
            # Vred .*= reshape(weights,1,:);
            nsamp = size(Vred,2)
            # norm_weights = weights'*weights
            Cexp = Vred*Vred'
            # Csky./=norm_weights
            Cexp./=nsamp

            SF = svd(Cexp);
            EVEC = zeros(length(wavetarg),size(SF.U,2))
            EVEC[chipgapmsk,:].=SF.U;

            h5write(savename,"Vmat",EVEC[:,1:nsub]*Diagonal(sqrt.(SF.S[1:nsub])))
            h5write(savename,"λv",SF.S[1:nsub])
            h5write(savename,"chipgapmsk",chipgapmsk)
        end
    end
end

# build_starCont(runlist_range)
@showprogress pmap(build_starCont,build_fibers)