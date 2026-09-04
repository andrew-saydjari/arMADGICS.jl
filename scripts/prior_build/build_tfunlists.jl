# build_tfunlists.jl — E4 step 1: rebuild tfun_sample_lst against the NEW E3 telluric
# refit artifacts (tellurics_refit_20260902_{apo,lco}.h5), applying the three
# USER-MANDATED consumption conditions:
#   C1: keep (fiber,exposure) entries with median transfer flux > 400
#       (medflux per viz.jl fig06: nzmed over all 8700 pixels of exp.(A*theta[:,fib,i]))
#   C2: exclude entries with median transfer flux above a per-telescope bright cut:
#       APO > 10,000 (fig16 bright-end nonlinearity regime); LCO > 3,000 (AKS 2026-09-04:
#       poisoning guard — the LCO fleet is faint, median ~872 / normal max 1,884, and the
#       only exposures above 3k are two anomalous domeflats at 5.6k/6.5k medflux whose
#       median chi_sq_fiber is 7-8x the fleet median: lco 59160 exp 0018, lco 60291
#       exp 0009. Any threshold in the 1,884-5,631 gap removes exactly those two.
#       The 2026-09-03 E4 run predates this LCO cut ("no LCO cut" mandate) — ~229
#       fiber-entries from those two exposures leaked past C3; assessed not worth an
#       E4 rebuild, to be applied at the pass-2 refit.)
#   C3: exclude entries with chi_sq_fiber above the per-telescope p99.9
#       (percentile recomputed exactly from the file, StatsBase.percentile(vec, 99.9))
# Output format matches what sample_starCont.jl consumes: a serialized Vector (push! onto
# `[]`, as in 2026_04_25 cells 33/39) of 300 per-fiber Vector{Int} exposure-index lists,
# each sorted by DESCENDING medflux (old-convention ordering).
# NOTE (deviation by design, per AKS conditions): the 2026_04_25-era rank-based bright-end
# trim (drop top 99 brightest, cap at N-100) is SUPERSEDED by the absolute C2 cut.
#
# Also: validation pass reproducing the OLD 2026_04_25 notebook logic (plain median,
# p[100:min(N-100,lindx)]) on the DELIVERED files and comparing to the old lists,
# to certify the medflux implementation.
#
# Run: nice -n 10 julia --project=/mnt/home/asaydjari/gitcode/worktrees/arM-E4run -t 16 build_tfunlists.jl

using LinearAlgebra, Statistics, StatsBase, HDF5, Serialization, Printf, Dates
BLAS.set_num_threads(1)   # threading is over exposures; keep gemm single-threaded

const NEWD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out"
const DELD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_08_31/prior_inputs/tellurics_20260220_arjl_domeflats"
const OLDL = "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25"
const OUTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902"
mkpath(OUTD)

report_io = open(joinpath(OUTD, "BUILD_REPORT.txt"), "w")
logmsg(args...) = (println(args...); println(report_io, args...); flush(stdout); flush(report_io))
logmsg("tfunlist build — ", now())
logmsg("julia $(VERSION), nthreads=$(Threads.nthreads())")

nzmed(v) = median(x for x in v if isfinite(x) && x != 0)

# medflux (300,N): per-(fiber,exposure) median of the transfer function exp(A*theta).
# f64 A and theta, per viz.jl fig06 (loadfile promotes both to Float64).
function medflux_of(path; plain_median::Bool=false)
    A, th = h5open(path, "r") do f
        Float64.(permutedims(read(f["design_matrix"]), (2, 1))), read(f["theta"])
    end
    N = size(th, 3)
    mf = zeros(300, N)
    nbadcol = Threads.Atomic{Int}(0)
    Threads.@threads for i in 1:N
        M = A * Float64.(@view th[:, :, i])   # (8700,300)
        map!(exp, M, M)
        for fb in 1:300
            col = @view M[:, fb]
            nb = count(x -> !isfinite(x) || x == 0, col)
            nb > 0 && Threads.atomic_add!(nbadcol, 1)
            mf[fb, i] = (plain_median || nb == 0) ? median(col) : nzmed(col)
        end
    end
    mf, nbadcol[]
end

# ── Validation: reproduce OLD 2026_04_25 lists from DELIVERED files ──────────────────────
for tele in ("lco", "apo")
    t0 = time()
    mf, _ = medflux_of(joinpath(DELD, "20260323_$(tele).h5"); plain_median=true)
    old = deserialize(joinpath(OLDL, "20260323_$(tele)_tfunlist.jdat"))
    N = size(mf, 2)
    nmatch = 0; ndiff_entries = 0
    for fb in 1:300
        p = sortperm(mf[fb, :], rev=true)
        lindx = findlast(mf[fb, p] .> 400)
        lindxf = minimum([N - 100, lindx])
        repro = p[100:lindxf]
        if repro == old[fb]
            nmatch += 1
        else
            ndiff_entries += length(symdiff(Set(repro), Set(old[fb])))
        end
    end
    logmsg(@sprintf("[VALIDATE %s] delivered-file old-convention reproduction: %d/300 fibers exact; symdiff entries across non-matching fibers: %d (%.0fs)",
        tele, nmatch, ndiff_entries, time() - t0))
end

# ── Build NEW lists with the three mandated conditions ───────────────────────────────────
for tele in ("apo", "lco")
    t0 = time()
    path = joinpath(NEWD, "tellurics_refit_20260902_$(tele).h5")
    mf, nbad = medflux_of(path)
    chif = h5open(path, "r") do f
        Float64.(read(f["chi_sq_fiber"]))   # (300,N)
    end
    N = size(mf, 2)
    @assert size(chif) == (300, N)
    nfin = count(isfinite, chif)
    p999 = percentile(vec(chif), 99.9)
    logmsg(@sprintf("\n[%s] %s", tele, path))
    logmsg(@sprintf("[%s] N exposures = %d; entries = %d; medflux cols w/ nonfinite-or-zero pixels: %d", tele, N, 300N, nbad))
    logmsg(@sprintf("[%s] chi_sq_fiber: finite %d/%d, p99.9 = %.4f (recomputed; MERGE_REPORT cross-check: apo 58328.47 / lco 295.91)", tele, nfin, 300N, p999))

    # per-telescope C2 thresholds (see header): APO physically motivated (fig16
    # nonlinearity regime); LCO a poisoning guard sitting in the wide gap between
    # the normal fleet (max 1,884) and the two known-bad bright domeflats (5.6k+).
    bright_cut = tele == "apo" ? 10_000.0 : 3_000.0
    fail1 = mf .<= 400.0                 # C1: not above 400
    fail2 = mf .> bright_cut             # C2: APO bright end
    fail3 = .!(chif .<= p999)            # C3: chi2 screen (NaN chi2 also excluded)
    keep = .!(fail1 .| fail2 .| fail3)

    logmsg(@sprintf("[%s] C1 (medflux<=400)   excludes %7d entries (%.2f%%)", tele, count(fail1), 100count(fail1) / 300N))
    logmsg(@sprintf("[%s] C2 (medflux>%g) excludes %7d entries (%.2f%%)  [only-C2: %d]", tele, bright_cut, count(fail2), 100count(fail2) / 300N, count(fail2 .& .!fail1 .& .!fail3)))
    logmsg(@sprintf("[%s] C3 (chi2>p99.9)     excludes %7d entries (%.2f%%)  [only-C3: %d]", tele, count(fail3), 100count(fail3) / 300N, count(fail3 .& .!fail1 .& .!fail2)))
    logmsg(@sprintf("[%s] KEPT %d / %d entries (%.2f%%)", tele, count(keep), 300N, 100count(keep) / 300N))
    # whole-exposure view of C2 for context (median over fibers)
    expo_med = [median(@view mf[:, i]) for i in 1:N]
    logmsg(@sprintf("[%s] context: exposures whose fiber-median medflux > %g: %d/%d", tele, bright_cut, count(>(bright_cut), expo_med), N))

    Tfunlists = []
    perfib = zeros(Int, 300)
    for fb in 1:300
        idx = findall(@view keep[fb, :])
        sort!(idx, by=i -> mf[fb, i], rev=true)   # old-convention descending-medflux order
        perfib[fb] = length(idx)
        push!(Tfunlists, idx)
    end
    logmsg(@sprintf("[%s] per-fiber kept: min %d (fiber %d), median %.0f, max %d; fibers with <100 usable: %d",
        tele, minimum(perfib), argmin(perfib), median(perfib), maximum(perfib), count(<(100), perfib)))
    @assert minimum(perfib) > 0 "fiber with EMPTY tfun list — sampling would fail"

    outlist = joinpath(OUTD, "20260902_$(tele)_tfunlist.jdat")
    serialize(outlist, Tfunlists)
    logmsg(@sprintf("[%s] wrote %s (%d B)", tele, outlist, filesize(outlist)))

    h5open(joinpath(OUTD, "tfunlist_audit_$(tele).h5"), "w") do f
        f["medflux"] = mf
        f["chi_p999"] = p999
        f["bright_cut"] = bright_cut
        f["faint_cut"] = 400.0
        f["kept_mask"] = UInt8.(keep)
        f["perfiber_kept"] = perfib
    end
    logmsg(@sprintf("[%s] done in %.0fs", tele, time() - t0))
end

logmsg("\nALL-TFUNLISTS-OK  ", now())
close(report_io)
