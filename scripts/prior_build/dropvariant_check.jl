## E5 drop-variant sensitivity check: screened vs unscreened built priors.
# Principal angles between the leading subspaces + relative eigenvalue drift.
# Usage: julia --project=plots_project dropvariant_check.jl [fiber ...]
using HDF5, LinearAlgebra, Printf

const B = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1"

function principal_angle(V1, V2, k)
    Q1 = Matrix(qr(V1[:, 1:k]).Q)
    Q2 = Matrix(qr(V2[:, 1:k]).Q)
    s = svdvals(transpose(Q1) * Q2)
    return acosd(clamp(minimum(s), 0, 1)), minimum(s)
end

fibers = isempty(ARGS) ? [10, 76, 295, 351, 460, 519] : parse.(Int, ARGS)
cases = [("skyCont k=10", "APOGEE_skycont_svd_30_f", 10),
    ("skyCont k=30", "APOGEE_skycont_svd_30_f", 30),
    ("GSPICE k=5", "APOGEE_skyline_faint_GSPICE_svd_120_f", 5),
    ("GSPICE k=30", "APOGEE_skyline_faint_GSPICE_svd_120_f", 30),
    ("faint k=30", "APOGEE_skyline_faint_svd_120_f", 30)]

for f in fibers
    n = lpad(f, 3, "0")
    for (tag, fn, k) in cases
        p1 = joinpath(B, "built", fn * n * ".h5")
        p2 = joinpath(B, "built_unscreened", fn * n * ".h5")
        if isfile(p1) && isfile(p2)
            V1 = h5read(p1, "Vmat")
            V2 = h5read(p2, "Vmat")
            a, c = principal_angle(V1, V2, k)
            l1 = h5read(p1, "λv")
            l2 = h5read(p2, "λv")
            dl = maximum(abs.(l1[1:k] .- l2[1:k]) ./ max.(l1[1:k], 1e-30))
            @printf("f%s %-13s worst principal angle %8.4f deg (cos=%.7f)  max rel dlambda %.3e\n",
                n, tag, a, c, dl)
        else
            @printf("f%s %-13s PENDING (screened=%s unscreened=%s)\n", n, tag, isfile(p1), isfile(p2))
        end
    end
end
