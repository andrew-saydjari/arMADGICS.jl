# Predict per-column draw changes for ALL 300 APO fibers under the FINAL
# (AKS-approved 2026-09-05) tfunlist vs the ORIGINAL 20260902 apo list, and
# certify that no removed entry (E2 exposures / F1 freaks / C3-tightened) can
# be drawn. Mirrors 2026_09_05/c2_analysis/rng_predict_intelligent.jl.
using Serialization, Random, Printf, Statistics
tl_orig = deserialize("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/20260902_apo_tfunlist.jdat")
tl_new = deserialize("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_final/20260905_apo_tfunlist.jdat")
frac = deserialize("/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/outsamptell_apo.jdat")
nfrac = size(frac, 2)
println("apo fracTellSamples columns: ", nfrac)
nsamp = 10_000
function draws(tlfib)
    rng = MersenneTwister(203)
    rand(rng, 4_000:1:10_000, nsamp); rand(rng, 0:1e-4:5, nsamp); rand(rng, 2.6:1e-4:3.6, nsamp)
    T = rand(rng, tlfib, nsamp); F = rand(rng, 1:nfrac, nsamp)
    T, F
end
diffcols = Dict{Int,Vector{Int}}()
nch = Int[]
for fb in 1:300
    To, Fo = draws(tl_orig[fb]); Tn, Fn = draws(tl_new[fb])
    removedset = Set(setdiff(tl_orig[fb], tl_new[fb]))
    @assert !any(in(removedset), Tn) "removed entry drawable on fiber $fb"
    cols = findall((To .!= Tn) .| (Fo .!= Fn))
    diffcols[fb] = cols
    push!(nch, length(cols))
end
@printf("predicted changed cols/fiber: min %d med %.0f max %d; fibers with 10000 changed: %d\n",
    minimum(nch), median(nch), maximum(nch), count(==(10000), nch))
for fb in (28, 76, 226)
    @printf("fiberindx %3d: list %d->%d; predicted changed cols %d\n",
        fb, length(tl_orig[fb]), length(tl_new[fb]), length(diffcols[fb]))
end
serialize("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/pass1c_run/rng_predict_final_apo.jdat", diffcols)
println("saved")
