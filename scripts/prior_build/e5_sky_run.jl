## E5 pass-1 sky-prior retraining driver (DR21 testbed, ccalin051, NO Slurm).
#
# Staged execution (each stage checkpointed; rerun-safe):
#   --stage=lists      pooled per-target sample lists (pooling lives here; light)
#   --stage=packs      per-exposure guard cache (merged combine_sky_fibers once/exposure)
#   --stage=gather     per-SOURCE-fiber sample files from packs
#   --stage=stats      per-sample screening statistics -> screens/e5_sample_stats.h5
#   --stage=assemble   per-TARGET pooled skyflux/skyivar/skymsk (screens applied;
#                      requires screens/e5_screen_decisions.h5 from the screen-design
#                      step unless E5_UNSCREENED=1)
#   --stage=decompose  merged get_sky_samples second stage (target-fiber FRESH E4b
#                      starCont basis) -> skycont/skyline files
#   --stage=build      merged build_skyCont + build_skyLines -> built priors
#
# Config via env (all have E5 defaults for the 2026_09_04 pass-1 layout):
#   E5_OUT        output root (default /mnt/ceph/.../2026_09_04/prior_outputs/sky_pass1)
#   E5_REDUX      testbed reduction base (default .../2026_09_03/testbed_run/)
#   E5_ALMANAC    testbed almanac (default .../2026_09_03/testbed_prep/almanac_testbed_dr21.h5)
#   E5_STARCONT   FRESH starCont prior prefix (default <E5_OUT>/starcont_fresh/APOGEE_starcont_svd_60_f)
#   E5_CHIPGAP    chip-gap mask (default .../2026_04_25/StarContChipGapMsk.h5 — static
#                 instrument mask, same file the E4b starCont builds used)
#   E5_NWORKERS   local workers (default 8)
#   E5_FIBERS     fiber subset "a:b" or comma list (default 1:600)
#   E5_SAMPLE_SUBDIR  sample dir name (default samples; QA drop-variant uses samples_unscreened)
#   E5_UNSCREENED     =1: assemble without screens (drop-variant QA)
# Author - Andrew Saydjari (E5 pass 1)

using Dates
t0 = now(); t_then = t0;
using InteractiveUtils
versioninfo()
using Distributed, DataFrames

stage = let s = filter(a -> startswith(a, "--stage="), ARGS)
    isempty(s) && error("usage: julia e5_sky_run.jl --stage=lists|packs|gather|stats|assemble|decompose|build")
    String(split(s[1], "=")[2])
end

e5_out = get(ENV, "E5_OUT", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1")
reduxBase = get(ENV, "E5_REDUX", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/testbed_run/")
almanacFile = get(ENV, "E5_ALMANAC", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/testbed_prep/almanac_testbed_dr21.h5")
starcont_prefix = get(ENV, "E5_STARCONT", joinpath(e5_out, "starcont_fresh", "APOGEE_starcont_svd_60_f"))
chipgap_msk_path = get(ENV, "E5_CHIPGAP", "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5")
nworkers_req = parse(Int, get(ENV, "E5_NWORKERS", "8"))
sample_subdir = get(ENV, "E5_SAMPLE_SUBDIR", "samples")
unscreened = get(ENV, "E5_UNSCREENED", "0") == "1"

fibers = let s = get(ENV, "E5_FIBERS", "1:600")
    sort(unique(reduce(vcat, [occursin(":", tok) ?
        collect(parse(Int, split(tok, ":")[1]):parse(Int, split(tok, ":")[2])) :
        [parse(Int, tok)] for tok in split(s, ",")])))
end

list_dir = joinpath(e5_out, "outlists")
pack_dir = joinpath(e5_out, "packs")
src_dir = joinpath(e5_out, "src_samples")
sample_dir = joinpath(e5_out, sample_subdir)
screen_dir = joinpath(e5_out, "screens")
built_dir = joinpath(e5_out, unscreened ? "built_unscreened" : "built")
for d in (e5_out, list_dir, screen_dir)
    mkpath(d)
end

println("e5_sky_run stage=$stage out=$e5_out fibers=$(first(fibers)):$(last(fibers)) nworkers=$nworkers_req unscreened=$unscreened")
flush(stdout)

proj_path = dirname(Base.active_project()) * "/"
addprocs(nworkers_req, exeflags=["--project=$proj_path"])
t_now = now(); println("Worker allocation took $(Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then)))"); t_then = t_now; flush(stdout)

@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(parse(Int, get(ENV, "E5_BLAS_THREADS", "1")))
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays
    using Interpolations, SparseArrays, ParallelDataTransfer, AstroTime, Suppressor
    using StatsBase, ProgressMeter, DataFrames, JLD2, FileIO
    using SortFilters, BasisFunctions, Random, Statistics
end
@passobj 1 workers() proj_path
@everywhere begin
    src_code = "$proj_path"
    include(src_code * "src/utils.jl")
    include(src_code * "src/gridSearch.jl")
    include(src_code * "src/componentAndPosteriors.jl")
    include(src_code * "src/fileNameHandling.jl")
    include(src_code * "src/ingest.jl")
    include(src_code * "src/lowRankPrescription.jl")
    include(src_code * "src/marginalizeEW.jl")
    include(src_code * "src/spectraInterpolation.jl")
    include(src_code * "src/chi2Wrappers.jl")
    include(src_code * "scripts/prior_build/prior_utils.jl")
    include(src_code * "scripts/prior_build/sample_sky_defs.jl")
    include(src_code * "scripts/prior_build/build_sky_defs.jl")
    include(src_code * "scripts/prior_build/gspice.jl")
    include(src_code * "scripts/prior_build/e5_sky_pass1_defs.jl")
end
t_now = now(); println("Worker loading took $(Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then)))"); t_then = t_now; flush(stdout)

using LibGit2
git_branch, git_commit, git_clean = initalize_git(proj_path)
println("branch=$git_branch commit=$git_commit clean=$git_clean"); flush(stdout)

@passobj 1 workers() reduxBase
@passobj 1 workers() almanacFile
@passobj 1 workers() list_dir
@passobj 1 workers() pack_dir
@passobj 1 workers() src_dir
@passobj 1 workers() sample_dir
@passobj 1 workers() screen_dir
@passobj 1 workers() built_dir

if stage == "lists"
    counts = e5_make_pooled_lists(almanacFile, list_dir)
    println("pooled candidate counts: min=$(minimum(counts)) med=$(median(counts)) max=$(maximum(counts))")
    for adjfib in (1, 30, 31, 76, 226, 300, 301, 388, 519, 600)
        offset = adjfib > 300 ? 300 : 0
        println("  target $adjfib: pool=$(collect(e5_pool_range(adjfib - offset))) n=$(counts[adjfib])")
    end

elseif stage == "packs"
    expo = e5_unique_exposures(list_dir)
    println("unique exposures: $(length(expo))"); flush(stdout)
    res = @showprogress pmap(x -> begin
        try
            e5_extract_pack(reduxBase, almanacFile, x[1], x[2], x[3]; pack_dir=pack_dir)
        catch err
            println("PACK ERROR $(x): $err"); flush(stdout)
            :error
        end
    end, expo)
    println("packs: done=$(count(==(:done), res)) exists=$(count(==(:exists), res)) nosky=$(count(==(:nosky), res)) error=$(count(==(:error), res))")

elseif stage == "gather"
    mkpath(src_dir)
    res = @showprogress pmap(adjfib -> begin
        try
            (adjfib, e5_gather_source_samples(adjfib, list_dir, pack_dir, src_dir)...)
        catch err
            println("GATHER ERROR fiber $adjfib: $err"); flush(stdout)
            (adjfib, -1, -1)
        end
    end, fibers)
    nc = [r[2] for r in res]; nk = [r[3] for r in res]
    println("gather: kept min=$(minimum(nk)) med=$(median(nk)) max=$(maximum(nk)); guard-dropped total=$(sum(nc)-sum(nk)) of $(sum(nc)) candidates")
    h5open(joinpath(screen_dir, "e5_guard_counts.h5"), "w") do fh
        fh["adjfiberindx"] = Int[r[1] for r in res]
        fh["n_candidates"] = Int.(nc)
        fh["n_kept_guard"] = Int.(nk)
    end

elseif stage == "stats"
    res = @showprogress pmap(adjfib -> (adjfib, e5_source_stats(adjfib, src_dir)), fibers)
    h5open(joinpath(screen_dir, "e5_sample_stats.h5"), "w") do fh
        for (adjfib, st) in res
            g = create_group(fh, lpad(adjfib, 3, "0"))
            g["mjd"] = st.mjd; g["expnum"] = st.expnum
            g["scale"] = st.scale; g["goodfrac"] = st.goodfrac
            g["chi2r_full"] = st.chi2r_full; g["chi2r_cont"] = st.chi2r_cont
        end
    end
    println("stats written for $(length(res)) fibers")

elseif stage == "assemble"
    drop_dict = Dict{Int,BitVector}()
    if !unscreened
        dec_path = joinpath(screen_dir, "e5_screen_decisions.h5")
        isfile(dec_path) || error("assemble: $dec_path missing; run the screen-design step (or set E5_UNSCREENED=1 for the drop-variant)")
        h5open(dec_path, "r") do fh
            for k in keys(fh)
                startswith(k, "drop_") || continue
                drop_dict[parse(Int, k[6:end])] = BitVector(read(fh[k]) .!= 0)
            end
        end
        println("screen decisions loaded for $(length(drop_dict)) source fibers; total drops=$(sum(sum, values(drop_dict)))")
    end
    @passobj 1 workers() drop_dict
    res = @showprogress pmap(adjfib -> begin
        try
            (adjfib, e5_assemble_pooled(adjfib, list_dir, src_dir, sample_dir;
                drop_dict=drop_dict, variant_unscreened=unscreened))
        catch err
            println("ASSEMBLE ERROR fiber $adjfib: $err"); flush(stdout)
            (adjfib, -1)
        end
    end, fibers)
    ns = [r[2] for r in res]
    println("assembled samples/target: min=$(minimum(ns)) med=$(median(ns)) max=$(maximum(ns)); failures=$(count(<(0), ns))")

elseif stage == "decompose"
    prior_dict = Dict{String,String}()
    prior_dict["starCont"] = starcont_prefix
    prior_dict["StarContChipGapMsk"] = chipgap_msk_path
    @passobj 1 workers() prior_dict
    res = @showprogress pmap(adjfib -> begin
        try
            get_sky_samples((adjfib, []); reduxBase=reduxBase, almanacFile=almanacFile,
                prior_dict=prior_dict, out_dir=sample_dir, loc_parallel=false)
            (adjfib, :ok)
        catch err
            println("DECOMPOSE ERROR fiber $adjfib: $err"); flush(stdout)
            (adjfib, :error)
        end
    end, fibers)
    println("decompose: ok=$(count(r -> r[2] == :ok, res)) error=$(count(r -> r[2] == :error, res))")

elseif stage == "build"
    mkpath(built_dir)
    res = @showprogress pmap(adjfib -> begin
        r = Dict{Symbol,Any}(:adjfib => adjfib)
        try
            t1 = time()
            build_skyCont(adjfib; sample_dir=sample_dir, chipgap_msk_path=chipgap_msk_path,
                out_dir=built_dir, nsub=30)
            r[:cont_s] = round(time() - t1, digits=1)
            t1 = time()
            build_skyLines(adjfib; sample_dir=sample_dir, chipgap_msk_path=chipgap_msk_path,
                out_dir=built_dir) # GSPICE knobs untouched: usamp_factor=7 etc. defaults
            r[:lines_s] = round(time() - t1, digits=1)
            r[:status] = :ok
        catch err
            r[:status] = :error
            r[:err] = sprint(showerror, err)[1:min(end, 2000)]
            println("BUILD ERROR fiber $adjfib (STOP-AND-REPORT per hold-until-crash policy): $(r[:err])"); flush(stdout)
        end
        r
    end, fibers)
    nok = count(r -> r[:status] == :ok, res)
    println("build: ok=$nok error=$(length(res) - nok)")
    for r in res
        if r[:status] != :ok
            println("  CRASHED fiber $(r[:adjfib]): $(r[:err])")
        end
    end
else
    error("unknown stage $stage")
end

t_now = now(); println("Stage $stage took $(Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then)))"); flush(stdout)
