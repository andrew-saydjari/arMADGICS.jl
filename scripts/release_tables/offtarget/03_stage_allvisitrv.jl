# 03_stage_allvisitrv.jl -- copy-and-swap: add `offtarget` to allVisitRV.
#
# CHAIN POSITION: offtarget step 3 of 3; produces the final allVisitRV (25 datasets).
# USAGE:   julia --project=scripts/release_tables 03_stage_allvisitrv.jl \
#              --src ALLVISITRV_24.h5 --flag offtarget_visit_flag.h5 --dst ALLVISITRV.h5
# INPUTS:  --src   allVisitRV with 24 columns (read-only)
#          --flag  offtarget_visit_flag.h5 from 02_sidecar_and_visits.jl (`offtarget`, one per visit)
# OUTPUTS: --dst   copy of --src plus exactly one new dataset, `offtarget` (int8, 0/1; chunks
#                  (65536,), SHUFFLE then deflate 9, like the host file's other datasets) with a
#                  `description` attribute (scalar vlen UTF-8 string, as h5py writes a str).
#                  Nothing existing is touched; an existing --dst is replaced.
# SOURCE:  work/2026_09_16/offtarget_sidecar/s4_stage_columns.py (0.4.0), via legacy_python/03_stage_allvisitrv.py.
#
# Julia port of legacy_python/03_stage_allvisitrv.py (2026-11-24); reproduces it exactly (dtype, shape, chunking,
# filter pipeline, values and attribute flavour, checked with compare_h5.py and compare_h5.jl). The file
# bytes are not identical (nor are two runs of the .py: the dataset object header carries creation
# timestamps, as with h5py's defaults; the object header is also laid out differently by a few bytes).
using HDF5, ArgParse
include(joinpath(@__DIR__, "..", "lib", "h5attrs.jl"))

const CHUNKS = (65536,)
const GZIP = 9
const DESC = "FPS off-target flag, constant within a visit (verified).  1 where the almanac fiber " *
    "table for this visit's (telescope, mjd, config_id, fiber) has assigned == 1 and " *
    "on_target == 0, i.e. the robotic positioner was assigned this target but did not " *
    "reach it, so the spectra are of whatever sky the fiber did reach while still " *
    "carrying the assigned target's identity.  Task 94's sdss_id 85604486 is the worked " *
    "example: a field star 9 arcmin away at S/N 96, clean on every other flag.  0 for all " *
    "plate-era visits -- the plate era has no assigned/on_target columns and cannot be " *
    "flagged this way.  Source: pass1_dr21/almanac/allobs_57618_61230.h5."

log(a...) = (println(a...); flush(stdout))

function parse_cli(args)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--src"
            help = "allVisitRV with 24 columns"; required = true
        "--flag"
            help = "offtarget_visit_flag.h5 from step 2"; required = true
        "--dst"
            help = "output allVisitRV (25 columns)"; required = true
    end
    return parse_args(args, s)
end

"Read the visit flag and convert to Int8 exactly as numpy's astype(np.int8) (wrap-around)."
function load_flag(path::AbstractString, n::Integer)
    raw = h5open(f -> read(f["offtarget"]), path, "r")
    flag = raw isa Vector{Int8} ? raw : (raw .% Int8)
    @assert size(flag) == (n,) size(flag)
    @assert all(x -> x == 0 || x == 1, flag)
    return flag
end

function stage(src::AbstractString, flagpath::AbstractString, dst::AbstractString)
    n = h5open(f -> length(f["sdss_id"]), src, "r")
    log("\n=== ", basename(dst), " ===")
    flag = load_flag(flagpath, n)
    log("flag loaded: $(length(flag)) rows, $(sum(Int, flag)) flagged")

    isfile(dst) && rm(dst)
    log("copying...")
    cp(src, dst)
    log("copied, $(filesize(dst)) bytes")

    h5open(dst, "r+") do f
        haskey(f, "offtarget") && error("refusing to clobber an existing offtarget dataset")
        # h5py's create_dataset(shuffle=True, compression="gzip") puts SHUFFLE before DEFLATE
        d = create_dataset(f, "offtarget", datatype(Int8), dataspace((n,));
                           chunk = CHUNKS, filters = [HDF5.Filters.Shuffle(), HDF5.Filters.Deflate(GZIP)])
        write(d, flag)
        put_str_attr!(d, "description", DESC)
        close(d)
    end
    h5open(dst, "r") do f
        d = f["offtarget"]
        log("added: dtype=int8 shape=$(size(d)) chunks=$(HDF5.get_chunk(d)) ",
            "comp=gzip$(GZIP) shuffle=True sum=$(sum(Int, read(d)))")
    end
    log("final size $(filesize(dst)) bytes")
    log("md5: " * first(split(read(`md5sum $dst`, String))))
    return nothing
end

function main(args)
    a = parse_cli(args)
    t = @elapsed stage(a["src"], a["flag"], a["dst"])
    log("work time: ", round(t; digits = 2), " s")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
