## E7: build per-fiber LSF-projected TH starLines priors with the new (FPI-fit) LSFs.
## Replaces the 2023-era materialized-LSF convolution step of build_starLines.jl
## (steps 1-2 of that script — the Korg sampling + Krylov compression — are NOT
## rerun: the canonical fullres basis Vout is a pinned input, REQUIRED to be the
## same Vmat that generated the deployed refLSF prior so that the runtime
## asymmetric (fiber-LSF fit / refLSF report) coefficient pairing stays valid).
##
## Per fiber (adjfibindx 1:600), per sub-pixel index findx 1:10:
##   Vsubpix[:,:,findx] = get_lsf_matrix(adjfibindx, lsfpath; fstep=(6-findx)/10) * Vout
## get_lsf_matrix is row-normalized internally (E1 contract), so the historical
## `./ nvecLSF` is subsumed. fstep convention: index 6 = no shift; the OUTPUT
## grid is shifted by +fstep so the spectrum appears shifted by offrng[findx] =
## (findx-6)/10 pix, matching indTenth() at runtime and the refLSF file.
##
## Deterministic — no RNG. Run on ccalin051 under nice; no Slurm.
# Author - Andrew Saydjari (script by Claude, E7)

import Pkg
using Dates
t0 = now()
using InteractiveUtils
versioninfo()
if !haskey(ENV, "ARM_SKIP_PKG_OPS")
    Pkg.update("ApogeeReduction")
    Pkg.instantiate()
    Pkg.precompile()
end
println("Package ops took $(Dates.canonicalize(Dates.CompoundPeriod(now()-t0)))")
flush(stdout)

using Distributed

nworkers_req = parse(Int, get(ENV, "ARM_STARLINES_NWORKERS", "8"))
proj_path = dirname(Base.active_project()) * "/"
addprocs(nworkers_req, exeflags = ["--project=$proj_path"])
println("Running Main on ", gethostname(), " with $(nworkers()) workers")
flush(stdout)

@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using HDF5, SparseArrays, Dates
    using ApogeeReduction: get_lsf_matrix, adjFiberIndx2FiberIndx
end

# ---- configuration (env-overridable, E6 pattern) ----
prior_dir = "/mnt/ceph/users/sdssv/work/asaydjari/"
vout_path = get(ENV, "ARM_STARLINES_VOUT",
    prior_dir * "2026_09_05/prior_inputs/starLines_fullres/APOGEE_stellar_kry_50_fullres.h5")
lsf_path_apo = get(ENV, "ARM_STARLINES_LSF_APO",
    prior_dir * "2026_08_31/prior_inputs/lsf_20260427/fpiLSFparams_REGULARIZED_apo_60861.h5")
lsf_path_lco = get(ENV, "ARM_STARLINES_LSF_LCO",
    prior_dir * "2026_08_31/prior_inputs/lsf_20260427/fpiLSFparams_REGULARIZED_lco_60861.h5")
out_dir = get(ENV, "ARM_STARLINES_OUTDIR",
    prior_dir * "2026_09_05/prior_outputs/starLines_perfiber/")
fiber_range = let s = get(ENV, "ARM_STARLINES_FIBERS", "1:600")
    lo, hi = parse.(Int, split(s, ":"))
    lo:hi
end
mkpath(out_dir)

for p in (vout_path, lsf_path_apo, lsf_path_lco)
    isfile(p) || error("missing input: $p")
end

@everywhere begin
    nsub_out = 50
    len_outwave = 8700
    # sub-pixel convention: findx 6 = no shift (matches indTenth and the
    # deployed refLSF file); output grid shifted by fstep = (6-findx)/10.
    fsteprng = (5 // 10):(-1 // 10):(-4 // 10)
end

@everywhere function build_fiber(adjfibindx, Vout, lsf_path_apo, lsf_path_lco, out_dir)
    lname = joinpath(out_dir,
        "APOGEE_stellar_kry_$(nsub_out)_subpix_f" * lpad(adjfibindx, 3, "0") * ".h5")
    isfile(lname) && return :cached
    lsfpath = adjfibindx > 300 ? lsf_path_lco : lsf_path_apo
    Vsubpix = zeros(len_outwave, nsub_out, length(fsteprng))
    for (findx, fstep) in enumerate(fsteprng)
        Ksp = get_lsf_matrix(adjfibindx, lsfpath; fstep = Float64(fstep))
        size(Ksp) == (len_outwave, size(Vout, 1)) ||
            error("fiber $adjfibindx findx $findx: LSF matrix size $(size(Ksp))")
        Vsubpix[:, :, findx] .= Ksp * Vout
    end
    any(!isfinite, Vsubpix) && error("fiber $adjfibindx: non-finite output")
    tmpname = lname * ".tmp"
    h5write(tmpname, "Vmat", Vsubpix)
    mv(tmpname, lname; force = true)
    return :built
end

Vout = h5read(vout_path, "Vmat")
@assert size(Vout) == (200001, nsub_out)
@everywhere workers() Vout = $Vout

println("Building fibers $fiber_range -> $out_dir")
flush(stdout)
t1 = now()
using ProgressMeter
res = @showprogress pmap(fiber_range) do adjfibindx
    t = @elapsed r = build_fiber(adjfibindx, Vout, lsf_path_apo, lsf_path_lco, out_dir)
    (adjfibindx, r, t)
end
println("Build of $(length(fiber_range)) fibers took $(Dates.canonicalize(Dates.CompoundPeriod(now()-t1)))")
built = count(x -> x[2] == :built, res)
cached = count(x -> x[2] == :cached, res)
times = [x[3] for x in res if x[2] == :built]
println("built=$built cached=$cached; per-fiber s: min/med/max = ",
    isempty(times) ? "n/a" : join(round.((minimum(times), sort(times)[end÷2+1], maximum(times)), digits = 1), "/"))
flush(stdout)
