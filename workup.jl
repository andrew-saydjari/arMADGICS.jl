## This is a workup script for arMADGICS
## This possibly has memory leaks for large batches...
## need to revisit this/robustify it.
import Pkg;
using Dates;
t0 = now();
t_then = t0;
using InteractiveUtils;
versioninfo();
Pkg.instantiate();
Pkg.precompile();

using Distributed, SlurmClusterManager, ArgParse, TimerOutputs
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
if "SLURM_NTASKS" in keys(ENV)
    using SlurmClusterManager
    addprocs(SlurmManager(), exeflags = ["--project=$proj_path"])
else
    addprocs(32, exeflags = ["--project=$proj_path"]) # change to a workers per node variable
end

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

rm.(deblendf)
rmprocs(workers())