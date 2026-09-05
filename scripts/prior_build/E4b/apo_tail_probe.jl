using HDF5, Statistics, StatsBase, Printf, Serialization
mf = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/tfunlist_audit_apo.h5","r") do f; read(f["medflux"]); end
chif = h5open("/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_apo.h5","r") do f; Float64.(read(f["chi_sq_fiber"])); end
N = size(mf,2); nonfpi = setdiff(1:300,[76,226])
ec = [median(@view chif[nonfpi,i]) for i in 1:N]
eb = [median(@view mf[nonfpi,i]) for i in 1:N]
println("APO chi2stat tail: counts above thresholds:")
for t in (200,500,1000,2000,5000,10000,20000,40000,50000)
    @printf("  >%6d: %4d exposures\n", t, count(>(t), ec))
end
s = sort(ec, rev=true)
println("ranked chi2stat 30..70: ", join(round.(Int, s[30:70]), " "))
bright = findall(>(10_000), eb)
println("bright-median exposures: ", length(bright), "; their chi2stat range: ", round(Int,minimum(ec[bright])), "-", round(Int,maximum(ec[bright])))
hi = findall(>(20_000), ec)
println("chi2stat>20k exposures: ", length(hi), "; overlap with bright set: ", length(intersect(hi,bright)))
# proposed E2 (exposure screen: chi2stat > 20k OR brightness-gap): delta accounting vs current apo policy
p999 = 58328.4680
base = (mf .> 400.0) .& (chif .<= p999)
cur = base .& (mf .<= 10_000.0)               # current APO keep
new = copy(base); for e in hi; new[:,e] .= false; end   # E2 keep (no ceiling)
restored = count(new .& .!cur); removed = count(cur .& .!new)
fibch = count(fb -> (@view new[fb,:]) != (@view cur[fb,:]), 1:300)
@printf("E2(chi2stat>20k) vs current: restored %d entries, newly removed %d, fibers changed %d/300\n", restored, removed, fibch)
r76 = count(new[76,:] .& .!cur[76,:]); r226 = count(new[226,:] .& .!cur[226,:])
println("restored on FPI f76: ", r76, "  f226: ", r226)
pf = [count(@view new[fb,:]) for fb in 1:300]
println("per-fiber kept under E2: min ", minimum(pf), " med ", median(pf))
# what >10k entries remain under E2, and their chi2
rem10 = new .& (mf .> 10_000)
@printf("kept >10k entries under E2: %d, chi2 med %.1f max %.0f\n", count(rem10), median(chif[rem10]), maximum(chif[rem10]))
serialize("apo_E2_probe.jdat", (ec=ec, eb=eb, hi=hi, bright=bright))
