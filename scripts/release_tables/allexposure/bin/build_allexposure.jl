#!/usr/bin/env julia
# build_allexposure.jl -- command-line entry point of the Julia allExposure builder.
#
# Parallelism: Distributed worker processes for the HDF5-bound reads (the mid-exposure-time
# harvest over ~144k small files and the per-night almanac join; HDF5.jl serializes all HDF5 calls
# inside one process, so threads cannot parallelize them), and threads (`julia -t N`) for the
# row-wise stages. Output is a single ordered write by the head process.
# Worker count: --workers, else SLURM_CPUS_ON_NODE when under Slurm, else 8. Never nproc.
#
# Environment (set by ../run_release_tables.sh; required on Ceph, as for the workup):
#   CEPHTWEAKS_LAZYIO=1  HDF5_USE_FILE_LOCKING=FALSE
using ArgParse, Distributed

function cli()
    s = ArgParseSettings(description = "Build the allExposure release table (single pass)")
    @add_arg_table! s begin
        "--full-list"; required = true; help = "workup raw run full_list_info.h5 (row identity)"
        "--almanac"; required = true; help = "almanac allobs_<start>_<end>.h5 used by the workup"
        "--apred"; required = true; help = "apred 1D output root (<mjd>/ar1Dunical_*_object.h5)"
        "--workup"; required = true; help = "arMADGICS workup directory (arMADGICS_out_*.h5)"
        "--exclusion-lists"; required = true; help = "dir with bad_nights.csv, bad_fibers.csv, rv_railed.csv"
        "--fiber-offsets"; default = ""; help = "fiber_offsets.csv; omit to build without the v_rad_corr columns"
        "--no-offtarget"; action = :store_true; help = "omit the offtarget column"
        "--out"; required = true; help = "output allExposure h5"
        "--report"; default = ""; help = "output JSON build report"
        "--cache-dir"; required = true; help = "cache for identity, exposure times, kernel outputs (reused if present)"
        "--workers"; arg_type = Int; default = 0; help = "Distributed workers (0: SLURM_CPUS_ON_NODE or 8)"
        "--python"; default = get(ENV, "RELEASE_TABLES_PYTHON", "python3"); help = "Python with numpy, h5py, astropy (legacy kernels)"
        "--iers-a"; default = ""; help = "pin an IERS-A finals2000A file for the barycentric kernel (default: astropy IERS_Auto, as 0.4.0)"
        "--numerics"; default = "legacy"; help = "legacy (bit-identical to 0.4.0) or native (pure Julia libm)"
    end
    parse_args(s)
end

const A = cli()
nw = A["workers"] > 0 ? A["workers"] : parse(Int, get(ENV, "SLURM_CPUS_ON_NODE", "8"))
proj = dirname(@__DIR__)
nw > 1 && addprocs(nw; exeflags = "--project=$proj")
@everywhere using AllExposure

cfg = BuildConfig(full_list = A["full-list"], almanac = A["almanac"], apred = A["apred"],
    workup = A["workup"], lists = A["exclusion-lists"], fiber_offsets = A["fiber-offsets"],
    offtarget = !A["no-offtarget"], out = A["out"], report = A["report"], cache = A["cache-dir"],
    python = A["python"], iers_a = A["iers-a"],
    numerics = A["numerics"] == "legacy" ? LegacyNumerics() : A["numerics"] == "native" ? NativeNumerics() :
               error("--numerics must be legacy or native"))
build(cfg)
