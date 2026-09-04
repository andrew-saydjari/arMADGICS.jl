# E4b step 1: rebuild the LCO tfunlist ONLY, with the committed C2_LCO = 3000 bright
# cut (AKS 2026-09-04 decision; scripts/prior_build/build_tfunlists.jl @ main c457154).
# The APO list is NOT regenerated — the 2026_09_03 APO list stays canonical, and the
# 2026_09_03 lco list is RETAINED (no deletions); this writes to a NEW directory.
#
# Logic is a verbatim subset of the committed builder (same medflux_of, same C1/C2/C3,
# same descending-medflux ordering, same output format), plus a delta-accounting
# section against the 2026_09_03 lco list.
#
# PRE-REGISTERED EXPECTATION vs REALITY (found in E4b pre-flight from the audit h5):
# the addendum expected the cut to newly exclude ONLY the two leaked domeflats
# (lco 59160 exp 0018 = row 2578, lco 60291 exp 0009 = row 3830; ~229 kept entries).
# But C2 is a PER-ENTRY cut, and per-entry the 1,884-5,631 "gap" does not exist:
# fiberindx 88 and 219 are systematically bright (~2,140 exposures each with
# medflux>3000, q95 ~5-6k) and fiberindx 148/159 have 36 such entries each.
# Predicted kept delta: 4,569 = 229 (leak) + ~4,340 (collateral, 4 fibers).
# This script REPORTS the exact attribution; the committed cut is applied as committed.
#
# Run: cd <outdir-parent> && nice -n 10 julia +1.11.6 \
#   --project=/mnt/home/asaydjari/gitcode/worktrees/arM-E4b -t 16 \
#   scripts/prior_build/E4b/build_tfunlist_lcoC2.jl

using LinearAlgebra, Statistics, StatsBase, HDF5, Serialization, Printf, Dates
BLAS.set_num_threads(1)

const NEWD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out"
const OLDLIST = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/20260902_lco_tfunlist.jdat"
const OUTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/tfunlists_lcoC2"
mkpath(OUTD)

report_io = open(joinpath(OUTD, "BUILD_REPORT.txt"), "a")
logmsg(args...) = (println(args...); println(report_io, args...); flush(stdout); flush(report_io))
logmsg("E4b LCO tfunlist build (C2_LCO=3000) — ", now())
logmsg("julia $(VERSION), nthreads=$(Threads.nthreads())")

nzmed(v) = median(x for x in v if isfinite(x) && x != 0)

# identical to committed build_tfunlists.jl
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

tele = "lco"
t0 = time()
path = joinpath(NEWD, "tellurics_refit_20260902_$(tele).h5")
mf, nbad = medflux_of(path)
chif = h5open(path, "r") do f
    Float64.(read(f["chi_sq_fiber"]))
end
N = size(mf, 2)
@assert size(chif) == (300, N)
nfin = count(isfinite, chif)
p999 = percentile(vec(chif), 99.9)
logmsg(@sprintf("\n[%s] %s", tele, path))
logmsg(@sprintf("[%s] N exposures = %d; entries = %d; medflux cols w/ nonfinite-or-zero pixels: %d", tele, N, 300N, nbad))
logmsg(@sprintf("[%s] chi_sq_fiber: finite %d/%d, p99.9 = %.4f (2026_09_03 report: 295.9077)", tele, nfin, 300N, p999))

bright_cut = 3_000.0   # committed C2_LCO (main @ c457154)
fail1 = mf .<= 400.0
fail2 = mf .> bright_cut
fail3 = .!(chif .<= p999)
keep = .!(fail1 .| fail2 .| fail3)

logmsg(@sprintf("[%s] C1 (medflux<=400)   excludes %7d entries (%.2f%%)", tele, count(fail1), 100count(fail1) / 300N))
logmsg(@sprintf("[%s] C2 (medflux>%g) excludes %7d entries (%.2f%%)  [only-C2: %d]", tele, bright_cut, count(fail2), 100count(fail2) / 300N, count(fail2 .& .!fail1 .& .!fail3)))
logmsg(@sprintf("[%s] C3 (chi2>p99.9)     excludes %7d entries (%.2f%%)  [only-C3: %d]", tele, count(fail3), 100count(fail3) / 300N, count(fail3 .& .!fail1 .& .!fail2)))
logmsg(@sprintf("[%s] KEPT %d / %d entries (%.2f%%)", tele, count(keep), 300N, 100count(keep) / 300N))
expo_med = [median(@view mf[:, i]) for i in 1:N]
logmsg(@sprintf("[%s] context: exposures whose fiber-median medflux > %g: %d/%d", tele, bright_cut, count(>(bright_cut), expo_med), N))

Tfunlists = []
perfib = zeros(Int, 300)
for fb in 1:300
    idx = findall(@view keep[fb, :])
    sort!(idx, by=i -> mf[fb, i], rev=true)
    perfib[fb] = length(idx)
    push!(Tfunlists, idx)
end
logmsg(@sprintf("[%s] per-fiber kept: min %d (fiber %d), median %.0f, max %d; fibers with <100 usable: %d",
    tele, minimum(perfib), argmin(perfib), median(perfib), maximum(perfib), count(<(100), perfib)))
@assert minimum(perfib) > 0 "fiber with EMPTY tfun list — sampling would fail"

outlist = joinpath(OUTD, "20260904_lco_tfunlist.jdat")
serialize(outlist, Tfunlists)
logmsg(@sprintf("[%s] wrote %s (%d B)", tele, outlist, filesize(outlist)))

h5open(joinpath(OUTD, "tfunlist_audit_lco.h5"), "w") do f
    f["medflux"] = mf
    f["chi_p999"] = p999
    f["bright_cut"] = bright_cut
    f["faint_cut"] = 400.0
    f["kept_mask"] = UInt8.(keep)
    f["perfiber_kept"] = perfib
end
logmsg(@sprintf("[%s] build done in %.0fs", tele, time() - t0))

# ── Delta accounting vs the 2026_09_03 lco list ─────────────────────────────────────────
logmsg("\n--- DELTA vs 2026_09_03 20260902_lco_tfunlist.jdat ---")
old = deserialize(OLDLIST)
LEAK = [2578, 3830]   # lco 59160/0018, 60291/0009 (E6-verified row indices)
leakset = Set(LEAK)

function delta_accounting(old, Tfunlists, leakset)
    n_removed = 0; n_added = 0
    removed_leak = 0; removed_other = 0
    removed_other_byfib = Dict{Int,Int}()
    order_ok = 0; fib_changed = Int[]
    for fb in 1:300
        oldfb = old[fb]; newfb = Tfunlists[fb]
        rem = setdiff(Set(oldfb), Set(newfb)); add = setdiff(Set(newfb), Set(oldfb))
        n_removed += length(rem); n_added += length(add)
        for e in rem
            if e in leakset
                removed_leak += 1
            else
                removed_other += 1
                removed_other_byfib[fb] = get(removed_other_byfib, fb, 0) + 1
            end
        end
        if newfb == filter(!in(rem), oldfb)
            order_ok += 1
        end
        (isempty(rem) && isempty(add)) || push!(fib_changed, fb)
    end
    n_removed, n_added, removed_leak, removed_other, removed_other_byfib, order_ok, fib_changed
end
n_removed, n_added, removed_leak, removed_other, removed_other_byfib, order_ok, fib_changed =
    delta_accounting(old, Tfunlists, leakset)
logmsg(@sprintf("entries removed: %d (leak exposures %s: %d; other: %d); entries ADDED: %d (expect 0)",
    n_removed, string(LEAK), removed_leak, removed_other, n_added))
logmsg(@sprintf("fibers changed: %d/300; fibers with order = old-minus-removed (positional stability): %d/300",
    length(fib_changed), order_ok))
for (fb, c) in sort(collect(removed_other_byfib), by=x->-x[2])
    logmsg(@sprintf("  non-leak removals: fiberindx %3d (adjfib %3d): %d entries (kept %d -> %d)",
        fb, fb+300, c, length(old[fb]), perfib[fb]))
end
# leaked-entry retention check
still = sum(count(in(leakset), Tfunlists[fb]) for fb in 1:300)
logmsg(@sprintf("leaked-exposure entries remaining in NEW lists: %d (MUST be 0)", still))
@assert still == 0
@assert n_added == 0

logmsg("\nE4B-LCO-TFUNLIST-OK  ", now())
close(report_io)
