## workup_serial.jl — W2 serial-tier CLI: Distributed readers + single
## streaming writer over the W1 row contract. See WorkupSerial.jl for design.
#
# Usage (example, fiber-subset workup of the 2026_05_01 partial corpus):
#   nice -n 10 julia --project=. workup_serial.jl \
#       --rawdir /mnt/ceph/.../outdir/arMADGICS/raw_57600_61160 \
#       --outdir /path/to/output_dir \
#       --fibers 1:3 --nworkers 8
#
# Raw batch files are NEVER deleted by this script (W5: deletion is a
# separate, later step gated on a passing W3 validation report).

import Pkg
using Dates
using Distributed
using ArgParse

function parse_commandline()
    s = ArgParseSettings(description = "W2 serial-tier arMADGICS workup (streaming, contract-driven)")
    @add_arg_table s begin
        "--rawdir"
        help = "arMADGICS raw batch dir (contains NNN/ fiber dirs, batch_info.txt, full_list_info.h5)"
        arg_type = String
        required = true
        "--outdir"
        help = "output directory for arMADGICS_out_<key>.h5 files (created if absent)"
        arg_type = String
        required = true
        "--fibers"
        help = "contiguous adjfiberindx window to work up, as F1:F2 (default all, 1:600)"
        arg_type = String
        default = "1:600"
        "--allow-missing"
        help = "tolerate missing batch files: their rows are sentinel-filled and flagged in a missing_row_mask dataset (default: any missing batch is a hard error)"
        action = :store_true
        "--nworkers"
        help = "distributed reader processes (0 = inline single-process)"
        arg_type = Int
        default = 8
        "--batchsize"
        help = "producer batch size"
        arg_type = Int
        default = 100
        "--resume"
        help = "resume from <outdir>/workup_serial.ckpt (skips batches already written + flushed)"
        action = :store_true
        "--ckpt-every"
        help = "flush outputs + append checkpoint every N batches"
        arg_type = Int
        default = 200
        "--progress-every"
        help = "log progress every N batches"
        arg_type = Int
        default = 100
    end
    return parse_args(s)
end

const parg = parse_commandline()

function parse_fibers(s::AbstractString)
    m = match(r"^(\d+):(\d+)$", strip(s))
    isnothing(m) && error("--fibers must be F1:F2 (got '$s')")
    return parse(Int, m.captures[1]):parse(Int, m.captures[2])
end

nw = parg["nworkers"]
nw > 0 && addprocs(nw; exeflags = "--project=$(Base.active_project())")

const RCPATH = joinpath(@__DIR__, "RowContract.jl")
const WSPATH = joinpath(@__DIR__, "WorkupSerial.jl")
@everywhere begin
    include($RCPATH)
    include($WSPATH)
end
using .RowContract
using .WorkupSerial

WorkupSerial.run_workup(parg["rawdir"], parg["outdir"];
    fibers = parse_fibers(parg["fibers"]),
    allow_missing = parg["allow-missing"],
    batchsize = parg["batchsize"],
    resume = parg["resume"],
    ckpt_every = parg["ckpt-every"],
    progress_every = parg["progress-every"])

nw > 0 && rmprocs(workers())
