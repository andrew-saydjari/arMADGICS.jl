using HDF5, Statistics, StatsBase, Printf
TELL = Dict("apo"=>"/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_apo.h5",
            "lco"=>"/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5")
AUD = Dict("apo"=>"/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/tfunlist_audit_apo.h5",
           "lco"=>"/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/tfunlists_lcoC2/tfunlist_audit_lco.h5")
FPI = Dict("apo"=>[76,226], "lco"=>[88,219])
for tele in ("apo","lco")
    mf = h5open(AUD[tele],"r") do f; read(f["medflux"]); end
    chif = h5open(TELL[tele],"r") do f; Float64.(read(f["chi_sq_fiber"])); end
    N = size(mf,2); nonfpi = setdiff(1:300, FPI[tele])
    base = (mf .> 400.0)
    if tele == "apo"
        # chi2 by class for >10k entries
        for (nm, fibs) in (("FPI fibers 76/226", FPI[tele]), ("non-FPI fibers", nonfpi))
            sel = falses(300,N); sel[fibs,:] .= true
            s10 = sel .& base .& (mf .> 10_000)
            slo = sel .& base .& (mf .<= 10_000)
            @printf("[apo] %s: n>10k=%d chi2 med %.0f (q25 %.0f q75 %.0f) | <=10k chi2 med %.1f\n",
                nm, count(s10), median(chif[s10]), quantile(chif[s10],.25), quantile(chif[s10],.75), median(chif[slo]))
        end
        # are non-FPI >10k entries concentrated in whole-bright exposures?
        em = [median(@view mf[nonfpi,i]) for i in 1:N]
        brightexp = findall(>(10_000), em)
        s10nf = falses(300,N); s10nf[nonfpi,:] .= true; s10nf .&= base .& (mf .> 10_000)
        inbright = sum(count(@view s10nf[:,i]) for i in brightexp; init=0)
        @printf("[apo] non-FPI >10k entries: %d total; %d (%.0f%%) inside the %d bright-median exposures\n",
            count(s10nf), inbright, 100inbright/count(s10nf), length(brightexp))
    end
    # exposure-level chi2 statistic (non-FPI median chi2 per exposure)
    ec = [median(@view chif[nonfpi,i]) for i in 1:N]
    m, s = median(ec), mad(ec, normalize=true)
    top = sortperm(ec, rev=true)[1:8]
    @printf("[%s] exposure chi2 (non-FPI median): med %.1f MAD %.1f; top-8 rows: %s vals: %s\n",
        tele, m, s, string(top), join(round.(ec[top],digits=0), " "))
    # gap isolation on ec
    ord = sortperm(ec, rev=true); cap = max(1, ceil(Int, 0.005N))
    istar = 0
    for i in 1:cap
        (ec[ord[i]] / ec[ord[i+1]] >= 2.0) && (istar = i)
    end
    @printf("[%s] chi2-exposure gap screen (G=2, cap %d): flags %d exposures %s\n",
        tele, cap, istar, istar>0 ? string(sort(ord[1:istar])) : "")
    # robust threshold view: exposures above med+10MAD
    for k in (10, 20, 50)
        @printf("[%s]   exposures with chi2stat > med+%d*MAD (=%.0f): %d\n", tele, k, m+k*s, count(>(m+k*s), ec))
    end
end
