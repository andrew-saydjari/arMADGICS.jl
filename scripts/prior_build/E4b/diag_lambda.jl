# Diagnostic: WHY do the pass-1c built priors show LOWER subleading variance
# than the pass-1 ceiling builds on f76/f226/f171, when the policy RESTORES
# entries? Decompose each fiber's list delta (restored vs removed) and the
# chi2/medflux character of each class.
using HDF5, Serialization, Statistics, StatsBase, Printf
P = "/mnt/ceph/users/sdssv/work/asaydjari/"
mf, p999f, p999u = h5open(P*"2026_09_05/tfunlists_final/tfunlist_audit_apo.h5","r") do f
    read(f["medflux"]), read(f["chi_p999_full"]), read(f["chi_p999_used"])
end
chif = h5open(P*"2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_apo.h5","r") do f
    Float64.(read(f["chi_sq_fiber"]))
end
new = deserialize(P*"2026_09_05/tfunlists_final/20260905_apo_tfunlist.jdat")
old = deserialize(P*"2026_09_03/tfunlists_refit20260902/20260902_apo_tfunlist.jdat")
lp(p) = h5open(p,"r") do f; read(f["λv"]); end
NB = P*"2026_09_05/prior_outputs/starCont_20260905_final/built"
OB = P*"2026_09_04/prior_outputs/starCont_pass1/built_apo"
@printf("APO C3 gate: full-fleet %.1f -> used(post-E2) %.1f\n", p999f, p999u)
println("fiber | old  new  | restored (chi2 med, mf med) | removed (chi2 med, mf med) | l1r l2r l3r")
lam2 = zeros(300)
for fb in 1:300
    o, n = Set(old[fb]), Set(new[fb])
    res = collect(setdiff(n,o)); rem = collect(setdiff(o,n))
    lo = lp(joinpath(OB,"APOGEE_starcont_svd_60_f"*lpad(fb,3,"0")*".h5"))
    ln = lp(joinpath(NB,"APOGEE_starcont_svd_60_f"*lpad(fb,3,"0")*".h5"))
    lam2[fb] = ln[2]/lo[2]
    if fb in (28,42,76,90,138,171,226,150,300)
        rc = isempty(res) ? NaN : median(chif[fb,res]); rm_ = isempty(res) ? NaN : median(mf[fb,res])
        dc = isempty(rem) ? NaN : median(chif[fb,rem]); dm = isempty(rem) ? NaN : median(mf[fb,rem])
        @printf("%5d | %5d %5d | %5d (%9.1f, %7.0f) | %5d (%9.1f, %7.0f) | %.3f %.3f %.3f\n",
            fb, length(old[fb]), length(new[fb]), length(res), rc, rm_, length(rem), dc, dm,
            ln[1]/lo[1], ln[2]/lo[2], ln[3]/lo[3])
    end
end
@printf("\nfleet lambda2 ratio: min %.3f (fiber %d)  p10 %.3f  med %.4f  p90 %.3f  max %.3f (fiber %d)\n",
    minimum(lam2), argmin(lam2), quantile(lam2,.1), median(lam2), quantile(lam2,.9), maximum(lam2), argmax(lam2))
@printf("fibers with lambda2 ratio < 0.9: %d ; > 1.1: %d\n", count(<(0.9),lam2), count(>(1.1),lam2))
# how many entries did the tightened gate remove per fiber, fleet-wide
ng = [count(j -> chif[fb,j] > p999u, old[fb]) for fb in 1:300]
@printf("entries per fiber removed by the tightened C3 gate: min %d med %.0f max %d (total %d)\n",
    minimum(ng), median(ng), maximum(ng), sum(ng))
# correlation: lambda2 drop vs gate removals
println("corr(lambda2 ratio, gate removals) = ", round(cor(lam2, Float64.(ng)), digits=3))
