## This is script grabs a bunch of sky fiber spectra and decomposes them into continuum and line components, to serve as samples for building the sky prior
# The main body lives in sample_sky_defs.jl (sample_sky_main) so tests can call it directly (M4).
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
    using StatsBase, ProgressMeter, DataFrames, JLD2, FileIO
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
    include(src_dir*"scripts/prior_build/sample_sky_defs.jl")
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

    # Input data locations (M4: these were previously undefined in this script)
    reduxBase = "/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/"
    almanacFile = reduxBase * "almanac/allobs_57600_61160.h5"

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # star-continuum SVD basis used as the smooth basis for the sky continuum fit
    # (M4: key was "skycont" while get_sky_samples read "starCont" -> KeyError)
    prior_dict["starCont"] = prior_dir*"2026_04_26/star_priors/APOGEE_starcont_svd_60_f"
    prior_dict["StarContChipGapMsk"] = prior_dir*"2026_04_25/StarContChipGapMsk.h5"
end

sample_sky_main(reduxBase, almanacFile, runlist_range;
    prior_dict=prior_dict, out_dir="sky_prior_disk/",
    loc_parallel=true, runlist_parallel=true)
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Sky sample generation took $dt"); flush(stdout);
