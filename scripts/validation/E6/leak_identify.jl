# E6 step 5a: identify the leaked lco exposures (59160/0018, 60291/0009) in the
# tellurics h5, find which tfunlist entries they contribute, and reproduce the
# per-fiber MersenneTwister(203) sample draws to count affected sample columns.
using HDF5, Serialization, Random, Statistics

tf = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5"
paths = h5open(tf, "r") do f
    read(f["paths"])
end
println("n paths = ", length(paths))
println("example path: ", paths[1])

leak_idx = Int[]
for (i, p) in enumerate(paths)
    if (occursin("59160", p) && occursin("0018", p)) || (occursin("60291", p) && occursin("0009", p))
        push!(leak_idx, i)
        println("LEAK candidate exposure idx $i: ", p)
    end
end

tl = deserialize("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/20260902_lco_tfunlist.jdat")
println("tfunlist: ", length(tl), " fibers, lengths ", extrema(length.(tl)))

# per-fiber: how many leaked entries are in the tfunlist
nleak_fib = [count(in(Set(leak_idx)), tl[fb]) for fb in 1:300]
println("total leaked entries across fibers: ", sum(nleak_fib), " (addendum says ~229)")

# reproduce the sampling draws for every lco fiber (adjfib 301..600 -> fiberindx via
# adjFiberIndx2FiberIndx). Draw order (sample_starCont.jl): Teff, Av, Rv, Tfunindx, Tfracindx.
nsamp = 10_000
Teff_rng = 4_000:1:10_000
Av_rng = 0:1e-4:5
Rv_rng = 2.6:1e-4:3.6
# fracTellSamples lco: need only its column count for the Tfrac draw (comes after Tfunindx,
# so it does not even matter for reproducing Tfunindx_lst; skip loading the 700MB file).

# same mapping as ApogeeReduction.adjFiberIndx2FiberIndx (mod1(adjfib,300))
adjFiberIndx2FiberIndx(adjfibindx) = mod1(adjfibindx, 300)

leakset = Set(leak_idx)
ncols_affected = zeros(Int, 300)
affected_cols = Dict{Int,Vector{Int}}()  # adjfib => sample column indices
for adjfib in 301:600
    fiberindx = adjFiberIndx2FiberIndx(adjfib)
    rng = MersenneTwister(203)
    rand(rng, Teff_rng, nsamp)
    rand(rng, Av_rng, nsamp)
    rand(rng, Rv_rng, nsamp)
    Tfunindx_lst = rand(rng, tl[fiberindx], nsamp)
    cols = findall(in(leakset), Tfunindx_lst)
    ncols_affected[fiberindx] = length(cols)
    if !isempty(cols)
        affected_cols[adjfib] = cols
    end
end
println("fibers (of 300 lco) with >=1 affected sample column: ", count(>(0), ncols_affected))
println("total affected sample columns across lco fibers: ", sum(ncols_affected))
top = sortperm(ncols_affected, rev=true)[1:10]
for fb in top
    println("fiberindx $fb (adjfib $(fb+300)): tfunlist leaked entries=", nleak_fib[fb],
        " of ", length(tl[fb]), "; affected sample cols=", ncols_affected[fb])
end
serialize("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E6_run/leak_affected_cols.jdat",
    (leak_idx=leak_idx, nleak_fib=nleak_fib, ncols_affected=ncols_affected, affected_cols=affected_cols))
println("saved leak_affected_cols.jdat")
