# Larger scan of the existing 46k batches. Signature: RV at grid edge AND NaN chi2.
# Classify: finite snr -> clean-own-flux (P3 sky-mean poisoning); NaN snr -> self-NaN target (A4/M2 class).
using HDF5, Random
raw = "/mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/arMADGICS/raw_57600_61160"
rng = MersenneTwister(7)
sampled_files = String[]
for fib in shuffle(rng, lpad.(1:89, 3, '0'))
    fdir = joinpath(raw, fib)
    isdir(fdir) || continue
    fl = filter(endswith(".h5"), readdir(fdir))
    append!(sampled_files, joinpath.(fdir, shuffle(rng, fl)[1:min(14, length(fl))]))
end
println("scanning $(length(sampled_files)) batch files")
nrow = 0; nsig = 0; np3 = 0; nself = 0; nfiles_p3 = 0
for (k, f) in enumerate(sampled_files)
    rv, chi2, snr = try
        h5open(f) do fh
            vec(read(fh["RV_pixoff_final"])), vec(read(fh["RVchi2_residuals"])), vec(read(fh["snr"]))
        end
    catch e
        continue
    end
    global nrow += length(rv)
    sig = (abs.(rv) .>= 70) .& isnan.(chi2)
    global nsig += count(sig)
    p3 = sig .& isfinite.(snr)
    global np3 += count(p3)
    global nself += count(sig .& .!isfinite.(snr))
    any(p3) && (global nfiles_p3 += 1)
    k % 300 == 0 && (println("  $k files, $nrow rows, sig=$nsig p3=$np3"); flush(stdout))
end
p = np3 / nrow
println("rows: $nrow  signature: $nsig  P3(clean-flux): $np3  self-NaN: $nself")
println("P3 per-row rate: $(round(p, sigdigits=3)) +/- $(round(sqrt(np3)/nrow, sigdigits=2))")
println("files (of $(length(sampled_files))) with >=1 P3 row: $nfiles_p3")
