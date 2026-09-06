using HDF5, LinearAlgebra, Printf, Statistics
B="/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1"
OLD="/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/sky_priors"
ang(V1,V2,k)=acosd(clamp(minimum(svdvals(Matrix(qr(V1[:,1:k]).Q)'*Matrix(qr(V2[:,1:k]).Q))),0,1))
println("=== drop-variant: principal angle vs truncation k, with eigenvalue context ===")
for f in [10, 460]
    n=lpad(f,3,"0")
    p1=joinpath(B,"built","APOGEE_skyline_faint_GSPICE_svd_120_f$n.h5")
    p2=joinpath(B,"built_unscreened","APOGEE_skyline_faint_GSPICE_svd_120_f$n.h5")
    V1=h5read(p1,"Vmat"); V2=h5read(p2,"Vmat"); l1=h5read(p1,"λv"); l2=h5read(p2,"λv")
    @printf("f%s GSPICE: lambda1=%.3e lambda30=%.3e lambda30/lambda1=%.2e  |dlam|/lam max(1..30)=%.2e\n",
        n, l1[1], l1[30], l1[30]/l1[1], maximum(abs.(l1[1:30].-l2[1:30])./l1[1:30]))
    @printf("  lambda ratios near truncation: l29/l30=%.4f l30/l31=%.4f l30/l35=%.4f\n",
        l1[29]/l1[30], l1[30]/l1[31], l1[30]/l1[35])
    print("  worst principal angle vs k: ")
    for k in [1,2,3,5,10,20,30,50]
        @printf("k=%d:%.2f ", k, ang(V1,V2,k))
    end
    println()
    # energy-weighted agreement: how much of prior 1's variance lies in prior 2's span
    for k in [5,30]
        Q2=Matrix(qr(V2[:,1:k]).Q)
        cap=sum(abs2, Q2'*V1[:,1:k])/sum(abs2, V1[:,1:k])
        @printf("  energy of V1(1:%d) captured by span V2(1:%d): %.4f\n", k, k, cap)
    end
end
println("\n=== old-vs-new skyCont: is the support the same? ===")
for f in [10, 460]
    n=lpad(f,3,"0")
    Vn=h5read(joinpath(B,"built","APOGEE_skycont_svd_30_f$n.h5"),"Vmat")
    Vo=h5read(joinpath(OLD,"APOGEE_skycont_svd_30_f$n.h5"),"Vmat")
    sn=vec(sum(abs.(Vn),dims=2)) .> 0; so=vec(sum(abs.(Vo),dims=2)) .> 0
    common = sn .& so
    @printf("f%s support: new=%d old=%d common=%d (new-only=%d old-only=%d)\n",
        n, count(sn), count(so), count(common), count(sn .& .!so), count(so .& .!sn))
    a_full = ang(Vn,Vo,5)
    a_com = ang(Vn[common,:],Vo[common,:],5)
    @printf("  worst principal angle k=5: full grid %.2f deg | restricted to common support %.2f deg\n", a_full, a_com)
    for k in [1,2,3]
        @printf("    k=%d common-support angle %.2f deg\n", k, ang(Vn[common,:],Vo[common,:],k))
    end
end
