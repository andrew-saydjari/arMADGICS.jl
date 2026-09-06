## Measure what the DEPLOYED (DR17-era) sky priors actually split bright/faint:
# pixel counts of the bright and faint submsk datasets in the 2025_07_31 dump.
using HDF5, Statistics, Printf
const D = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/sky_priors"
fibs = vcat(1:25:300, 301:25:600)
bf = Float64[]; ff = Float64[]; frac = Float64[]
for adjfib in fibs
    n = lpad(adjfib, 3, "0")
    pb = joinpath(D, "APOGEE_skyline_bright_svd_120_f$n.h5")
    pf = joinpath(D, "APOGEE_skyline_faint_svd_120_f$n.h5")
    (isfile(pb) && isfile(pf)) || continue
    nb = count(Bool.(h5read(pb, "submsk")))
    nf = count(Bool.(h5read(pf, "submsk")))
    push!(bf, nb); push!(ff, nf); push!(frac, nb / (nb + nf))
    if adjfib in (1, 151, 301, 451)
        @printf("f%s bright=%d faint=%d bright_frac=%.4f\n", n, nb, nf, nb / (nb + nf))
    end
end
apo = 1:count(<=(300), fibs)
@printf("\nDR17-era deployed priors (%d fibers sampled):\n", length(bf))
@printf("  APO bright pixels: med=%.0f  bright fraction med=%.4f (range %.4f-%.4f)\n",
    median(bf[apo]), median(frac[apo]), minimum(frac[apo]), maximum(frac[apo]))
lco = (length(apo)+1):length(bf)
@printf("  LCO bright pixels: med=%.0f  bright fraction med=%.4f (range %.4f-%.4f)\n",
    median(bf[lco]), median(frac[lco]), minimum(frac[lco]), maximum(frac[lco]))
