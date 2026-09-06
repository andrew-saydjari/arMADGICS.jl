using HDF5, Statistics, Printf
P="/mnt/ceph/users/sdssv/work/asaydjari/"
NB=P*"2026_09_05/prior_outputs/starCont_20260905_final/built"; OB=P*"2026_09_04/prior_outputs/starCont_pass1/built_apo"
lp(p)=h5open(p,"r") do f; read(f["λv"]); end
tot=zeros(300); tail=zeros(300)
for fb in 1:300
    lo=lp(joinpath(OB,"APOGEE_starcont_svd_60_f"*lpad(fb,3,"0")*".h5")); ln=lp(joinpath(NB,"APOGEE_starcont_svd_60_f"*lpad(fb,3,"0")*".h5"))
    tot[fb]=sum(ln)/sum(lo); tail[fb]=median(ln[5:60]./lo[5:60])
    fb in (76,226,171,244,90) && @printf("f%3d: sum(lam) ratio %.4f ; median ratio k=5..60 %.3f ; lam1 %.4f lam2 %.3f\n", fb, tot[fb], tail[fb], ln[1]/lo[1], ln[2]/lo[2])
end
@printf("\nfleet sum(lam) ratio: med %.4f  p10 %.4f  p90 %.4f\n", median(tot), quantile(tot,.1), quantile(tot,.9))
@printf("fleet median-ratio(k=5..60): med %.3f  p10 %.3f  p90 %.3f ; fibers with >1: %d/300\n", median(tail), quantile(tail,.1), quantile(tail,.9), count(>(1.0),tail))
