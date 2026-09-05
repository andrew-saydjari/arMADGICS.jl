## E7 QA: structural + old-vs-new checks on the per-fiber TH starLines priors.
##
## (a) structural: every built file has key "Vmat", Float64 (8700,50,10),
##     finite; column norms sane vs the refLSF prior.
## (b) old-vs-new: for fibers with an old 2024-era norm94 file on disk,
##     per-component cosine similarity (findx 6) and principal subspace angles
##     between the 50-dim old/new bases; expect smooth LSF-shaped differences.
## (LSF cross-check vs the 2026_04_26 materialized matrices is a separate
##  script: qa_lsf_crosscheck.jl)
##
## Usage: julia --project=<E7 env> qa_starlines_perfiber.jl
## Env overrides: ARM_STARLINES_OUTDIR, ARM_E7_OLDREF_DIR, ARM_E7_QA_OUT

using LinearAlgebra, HDF5, Statistics, Printf

prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
new_dir = get(ENV, "ARM_STARLINES_OUTDIR",
    prior_dir * "2026_09_05/prior_outputs/starLines_perfiber/")
old_dir = get(ENV, "ARM_E7_OLDREF_DIR",
    prior_dir * "2026_09_05/E7_run/old_norm94_ref/")
qa_out = get(ENV, "ARM_E7_QA_OUT",
    prior_dir * "2026_09_05/E7_run/qa_starlines_results.h5")
reflsf_path = prior_dir * "2025_07_31/prior_dump/APOGEE_stellar_kry_50_subpix_th_22500.h5"

fname(dir, fib) = joinpath(dir, "APOGEE_stellar_kry_50_subpix_f" * lpad(fib, 3, "0") * ".h5")

Vref = h5read(reflsf_path, "Vmat")
refnorms = [norm(Vref[:, k, 6]) for k in 1:size(Vref, 2)]

## (a) structural sweep over all built files
built = [fib for fib in 1:600 if isfile(fname(new_dir, fib))]
println("built files present: $(length(built))/600")
nbad = 0
colnorm_ratio = fill(NaN, 50, length(built)) # ||new_k|| / ||ref_k|| at findx 6
for (j, fib) in enumerate(built)
    V = h5read(fname(new_dir, fib), "Vmat")
    ok = size(V) == (8700, 50, 10) && eltype(V) == Float64 && all(isfinite, V)
    if !ok
        global nbad += 1
        println("STRUCTURAL FAIL f$fib: size=$(size(V)) eltype=$(eltype(V)) finite=$(all(isfinite, V))")
    end
    colnorm_ratio[:, j] .= [norm(V[:, k, 6]) for k in 1:50] ./ refnorms
end
println("structural failures: $nbad")
if !isempty(built)
    r = vec(colnorm_ratio)
    @printf("colnorm ratio new/ref (findx6): min %.4f  med %.4f  max %.4f\n",
        minimum(r), median(r), maximum(r))
end

## (b) old-vs-new comparisons where old files exist
old_fibs = [fib for fib in built if isfile(fname(old_dir, fib))]
println("old-reference fibers available: ", old_fibs)
cosims = fill(NaN, 50, length(old_fibs))
prin_angles = fill(NaN, 50, length(old_fibs))
for (j, fib) in enumerate(old_fibs)
    Vn = h5read(fname(new_dir, fib), "Vmat")[:, :, 6]
    Vo = h5read(fname(old_dir, fib), "Vmat")[:, :, 6]
    # per-component cosine similarity (components are paired by construction:
    # both are K*Vout of the same Vout)
    for k in 1:50
        cosims[k, j] = dot(Vn[:, k], Vo[:, k]) / (norm(Vn[:, k]) * norm(Vo[:, k]))
    end
    # principal angles between the spanned subspaces
    Qn = Matrix(qr(Vn).Q)
    Qo = Matrix(qr(Vo).Q)
    s = svd(Qn' * Qo).S
    prin_angles[:, j] .= acosd.(clamp.(s, -1, 1))
    @printf("f%03d: cos-sim comp1-5: %s | worst comp: %.4f (k=%d) | max principal angle %.2f deg\n",
        fib, join([@sprintf("%.4f", cosims[k, j]) for k in 1:5], ","),
        minimum(cosims[:, j]), argmin(cosims[:, j]), maximum(prin_angles[:, j]))
end

h5open(qa_out, "w") do fh
    fh["built_fibers"] = built
    fh["colnorm_ratio_findx6"] = colnorm_ratio
    fh["old_fibers"] = old_fibs
    fh["cosine_similarity_findx6"] = cosims
    fh["principal_angles_deg"] = prin_angles
    HDF5.attrs(fh)["new_dir"] = new_dir
    HDF5.attrs(fh)["old_dir"] = old_dir
    HDF5.attrs(fh)["reflsf"] = reflsf_path
end
println("QA results -> $qa_out")
