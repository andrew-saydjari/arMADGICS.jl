## workup_old_reference.jl — the OLD in-memory workup (repo-root workup.jl at
## origin/main 4ef02ca), minimally adapted so it can run as the W4 regression
## reference on a fiber-subset mini-corpus. The accumulation/writing logic is
## UNCHANGED — byte-for-byte the same key discovery, size-sorted key order,
## sorted-glob positional row concatenation, per-key output files, hdr attr
## copy — so its outputs are exactly what the historical script would have
## produced on that corpus.
#
# ADAPTATIONS vs workup.jl (each is orthogonal to the numerics):
#   1. `using SlurmClusterManager` removed from the import line and the
#      SLURM_NTASKS addprocs branch deleted (the regression runs on a
#      workstation, never under Slurm).
#   2. `addprocs(32)` → `addprocs(8)` (shared-node courtesy cap for the
#      regression session; worker count does not affect output values).
#   3. The final `rm.(deblendf)` REMOVED: the old script deleted every raw
#      batch file after the workup (the M8 hazard closed by W5). The
#      regression must never delete corpus files — and the mini-corpus
#      entries are symlinks into the real 2026_05_01 corpus.
#   4. `Pkg.instantiate(); Pkg.precompile()` dropped (run via an already
#      instantiated scratch project).
#
# The mini-corpus handed to --outdir must contain: fiber subdirectories
# (symlinks are fine), and a batch_info.txt truncated to exactly the rows of
# those fibers (the old script derives nsamp = number of batch_info rows and
# fills rows positionally from the sorted file glob — so the fiber subset
# must be complete, with no missing batches, for its output to be correct).

using Dates;
t0 = now();
t_then = t0;
using InteractiveUtils;
versioninfo();

using Distributed, ArgParse, TimerOutputs
using Glob, DelimitedFiles, HDF5, EllipsisNotation
t_now = now();
dt = Dates.canonicalize(Dates.CompoundPeriod(t_now - t_then));
println("Package activation took $dt");
t_then = t_now;
flush(stdout);

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--outdir"
        required = false
        help = "output directory"
        arg_type = String
        default = "../outdir/arMADGICS/raw/"
    end
    return parse_args(s)
end
parg = parse_commandline()

foldname = parg["outdir"]
savename = replace(parg["outdir"], "raw" => "wu_th") * "arMADGICS_out.h5"
dirName = splitdir(savename)[1]
if !ispath(dirName)
    mkpath(dirName)
end

deblendf = sort(glob("*/*batch*.h5", "$foldname"));
lfile = length(deblendf)

if lfile == 0
    println("No batches found in $foldname. Exiting...")
    exit()
end

batchinfo = readdlm(foldname * "batch_info.txt", ',', comments=true)
nsamp = size(batchinfo, 1)

println("Total spectra processed: $nsamp")
println("Number of files: $lfile")

f = h5open(deblendf[1])
keyslst = keys(f)
hdr = h5readattr(deblendf[1], "hdr")
deleteat!(keyslst, keyslst .== "hdr")
keysizelst = []
keytypelst = []
for (indkey, keyval) in enumerate(keyslst)
    push!(keysizelst, size(f[keyval]))
    push!(keytypelst, eltype(f[keyval]))
end
close(f)

dataset_size = prod.(keysizelst) .* sizeof.(keytypelst);
p = sortperm(dataset_size);

keyslst = keyslst[p]
keysizelst = keysizelst[p]
keytypelst = keytypelst[p];

proj_path = dirname(Base.active_project()) * "/"
addprocs(8, exeflags = ["--project=$proj_path"])

@everywhere begin
    using HDF5, ProgressMeter
    function get_data(intup)
        fnameval, keyval = intup
        try
            return h5read(fnameval, keyval)
        catch
            println("Error reading $keyval from $fnameval")
            return h5read(fnameval, keyval)
        end
    end
end

@showprogress for (indkey, keyval) in enumerate(keyslst)
    savename_sub = chop(savename, tail=3) * "_" * keyval * ".h5"
    if !isfile(savename_sub)
        println("Started accumulating: $keyval")
        flush(stdout)
        if (indkey != 1) & (@isdefined data_out)
            cond = (keysizelst[indkey][1:(end-1)] == keysizelst[indkey-1][1:(end-1)])
            cond &= (keytypelst[indkey] == keytypelst[indkey-1])
            if cond
                fill!(data_out, 0)
            else
                global data_out = zeros(keytypelst[indkey], keysizelst[indkey][1:(end-1)]..., nsamp)
            end
        else
            global data_out = zeros(keytypelst[indkey], keysizelst[indkey][1:(end-1)]..., nsamp)
        end

        itarg = Iterators.product(deblendf, [keyval])
        pout = @showprogress pmap(get_data, itarg)
        println("Started writing: $keyval")
        flush(stdout)
        bind = 1
        for ind in 1:length(pout)
            @views sublen = size(pout[ind])[end]
            @views data_out[.., bind:(bind+sublen-1)] .= pout[ind]
            bind += sublen
        end
        g = h5open(savename_sub, "w")
        write(g, "hdr", "This is only a header")
        write(g, keyval, data_out)
        close(g)
        println("Completed writing: $keyval, Number of Observed Samples: $(bind-1), Number of files: $lfile")
        flush(stdout)
        pout = nothing
        @everywhere GC.gc()
        h5writeattr(savename_sub, "hdr", hdr)
    end
end

rmprocs(workers())
