# Intelligent bright-cut tfunlist build (AKS overnight directive 2026-09-05):
# "remove bad exposures, don't penalize bright fibers, only remove the most
# egregious outlier exposures."
# Policy (both telescopes, fresh build from the E3 refit inputs):
#   C1: medflux > 400 (unchanged faint cut)
#   C3: chi_sq_fiber <= per-telescope p99.9 (unchanged, recomputed)
#   E1: GAP-ISOLATED bright-exposure screen — sort exposure fiber-median medflux
#       descending; flag the smallest top prefix (capped at 0.1% of exposures)
#       ending at a multiplicative gap em[i]/em[i+1] >= G=2. Flagged exposures are
#       removed for ALL fibers. Data: LCO top medians 6538, 5631 | gap 2.99x | 1884
#       (fleet continuum) -> flags exactly the two leaked domeflats (rows 2578,
#       3830 = lco 59160/0018, 60291/0009; 34-44 fleet-MADs above the median).
#       APO top medians 15834..15780 (ratio ~1.001, smooth continuum) -> flags
#       NOTHING: bright domeflats are legitimate calibration variety at APO.
#   A1 (APO only): absolute per-entry ceiling medflux <= 10,000 — detector
#       nonlinearity physics (fig16), applies regardless of fiber brightness.
# NO per-fiber trims: a systematically bright fiber (LCO FPI 88/219; APO 28/226)
# is normal for itself. Per-fiber MAD screens were rejected because the FPI
# fibers are BIMODAL — their bright mode is normal, and C3 already screens shape
# anomalies per entry. See 2026_09_05/c2_analysis/C2_POLICY_REPORT.txt.
# Expected deltas: APO list IDENTICAL to 2026_09_03 20260902_apo_tfunlist.jdat;
# LCO = 2026_09_03 original minus the 229 kept leak entries (FPI collateral of
# the superseded C2_LCO=3000 restored).
#
# Run: nice -n 10 julia +1.11.6 --project=<arM-E4b> -t 16 build_tfunlists_intelligent.jl

using LinearAlgebra, Statistics, StatsBase, HDF5, Serialization, Printf, Dates
BLAS.set_num_threads(1)

const NEWD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out"
const OUTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_intelligent"
const REF = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/20260902_apo_tfunlist.jdat",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/20260902_lco_tfunlist.jdat")
const REF_LCOC2 = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/tfunlists_lcoC2/20260904_lco_tfunlist.jdat"
const LEAK = [2578, 3830]
mkpath(OUTD)

report_io = open(joinpath(OUTD, "BUILD_REPORT.txt"), "a")
logmsg(args...) = (println(args...); println(report_io, args...); flush(stdout); flush(report_io))
logmsg("intelligent tfunlist build — ", now())
logmsg("julia $(VERSION), nthreads=$(Threads.nthreads())")

nzmed(v) = median(x for x in v if isfinite(x) && x != 0)
function medflux_of(path)   # identical to committed build_tfunlists.jl
    A, th = h5open(path, "r") do f
        Float64.(permutedims(read(f["design_matrix"]), (2, 1))), read(f["theta"])
    end
    N = size(th, 3)
    mf = zeros(300, N)
    Threads.@threads for i in 1:N
        M = A * Float64.(@view th[:, :, i])
        map!(exp, M, M)
        for fb in 1:300
            col = @view M[:, fb]
            nb = count(x -> !isfinite(x) || x == 0, col)
            mf[fb, i] = nb == 0 ? median(col) : nzmed(col)
        end
    end
    mf
end

# E1: gap-isolation flagging of exposure medians
function gap_flag(em; G=2.0, capfrac=0.001)
    N = length(em)
    ord = sortperm(em, rev=true)
    cap = max(1, ceil(Int, capfrac * N))
    istar = 0
    for i in 1:cap
        if em[ord[i]] / em[ord[i+1]] >= G
            istar = i   # deepest gap within the cap
        end
    end
    istar == 0 ? Int[] : sort(ord[1:istar])
end

for tele in ("apo", "lco")
    t0 = time()
    path = joinpath(NEWD, "tellurics_refit_20260902_$(tele).h5")
    mf = medflux_of(path)
    chif = h5open(path, "r") do f
        Float64.(read(f["chi_sq_fiber"]))
    end
    N = size(mf, 2)
    p999 = percentile(vec(chif), 99.9)
    em = [median(@view mf[:, i]) for i in 1:N]
    m, s = median(em), mad(em, normalize=true)
    flagged = gap_flag(em)
    logmsg(@sprintf("\n[%s] %s  N=%d  chi2 p99.9=%.4f", tele, path, N, p999))
    logmsg(@sprintf("[%s] exposure medians: med %.0f MAD %.0f; E1 gap screen (G=2, cap %d) flags %d exposures: %s",
        tele, m, s, max(1, ceil(Int, 0.001N)), length(flagged), string(flagged)))
    for e in flagged
        logmsg(@sprintf("[%s]   flagged exp row %d: median medflux %.0f (%.0f MADs above fleet median)",
            tele, e, em[e], (em[e] - m) / s))
    end

    fail1 = mf .<= 400.0
    fail3 = .!(chif .<= p999)
    failE = falses(300, N); for e in flagged; failE[:, e] .= true; end
    failA = tele == "apo" ? (mf .> 10_000.0) : falses(300, N)
    keep = .!(fail1 .| fail3 .| failE .| failA)

    logmsg(@sprintf("[%s] C1 (mf<=400)      excludes %7d", tele, count(fail1)))
    logmsg(@sprintf("[%s] C3 (chi2>p99.9)   excludes %7d", tele, count(fail3)))
    logmsg(@sprintf("[%s] E1 (gap exposures) excludes %6d entries (%d whole exposures)", tele, count(failE), length(flagged)))
    logmsg(@sprintf("[%s] A1 (apo mf>10k)   excludes %7d", tele, count(failA)))
    logmsg(@sprintf("[%s] KEPT %d / %d (%.2f%%)", tele, count(keep), 300N, 100count(keep) / 300N))

    Tfunlists = []
    perfib = zeros(Int, 300)
    for fb in 1:300
        idx = findall(@view keep[fb, :])
        sort!(idx, by=i -> mf[fb, i], rev=true)
        perfib[fb] = length(idx)
        push!(Tfunlists, idx)
    end
    logmsg(@sprintf("[%s] per-fiber kept: min %d (fiber %d), median %.0f, max %d", tele,
        minimum(perfib), argmin(perfib), median(perfib), maximum(perfib)))
    @assert minimum(perfib) > 100
    if tele == "lco"
        logmsg(@sprintf("[lco] FPI fibers kept (MUST be ~full lists): f88=%d f219=%d f148=%d f159=%d",
            perfib[88], perfib[219], perfib[148], perfib[159]))
    end

    outlist = joinpath(OUTD, "20260905_$(tele)_tfunlist.jdat")
    serialize(outlist, Tfunlists)
    h5open(joinpath(OUTD, "tfunlist_audit_$(tele).h5"), "w") do f
        f["medflux"] = mf; f["chi_p999"] = p999
        f["expo_median"] = em; f["E1_flagged"] = flagged
        f["kept_mask"] = UInt8.(keep); f["perfiber_kept"] = perfib
    end
    logmsg(@sprintf("[%s] wrote %s (%d B) in %.0fs", tele, outlist, filesize(outlist), time() - t0))

    # delta accounting vs the ORIGINAL 2026_09_03 list
    ref = deserialize(REF[tele])
    nident = count(fb -> Tfunlists[fb] == ref[fb], 1:300)
    changed = [fb for fb in 1:300 if Tfunlists[fb] != ref[fb]]
    nrem = sum(length(setdiff(Set(ref[fb]), Set(Tfunlists[fb]))) for fb in 1:300)
    nadd = sum(length(setdiff(Set(Tfunlists[fb]), Set(ref[fb]))) for fb in 1:300)
    remleak = sum(count(in(Set(LEAK)), collect(setdiff(Set(ref[fb]), Set(Tfunlists[fb])))) for fb in 1:300)
    logmsg(@sprintf("[%s] vs ORIGINAL 20260902 list: identical fibers %d/300; removed %d (leak %d, other %d); ADDED %d (expect 0)",
        tele, nident, nrem, remleak, nrem - remleak, nadd))
    if tele == "lco"
        refc2 = deserialize(REF_LCOC2)
        identc2 = count(fb -> Tfunlists[fb] == refc2[fb], 1:300)
        diffc2 = [fb for fb in 1:300 if Tfunlists[fb] != refc2[fb]]
        logmsg(@sprintf("[lco] vs lcoC2 (C2=3000) list: identical fibers %d/300; differing fibers: %s  => ONLY these need resampling (their lcoC2 samples are invalid under this policy; all other fibers' lcoC2 samples remain exactly valid)",
            identc2, string(diffc2)))
        serialize(joinpath(OUTD, "lco_changed_fibers.jdat"),
            (vs_original=changed, vs_lcoC2=diffc2))
    else
        logmsg(nident == 300 ? "[apo] APO LIST IDENTICAL TO 20260903 — no APO resampling or rebuild needed" :
            "[apo] APO LIST CHANGED — resampling required for $(length(changed)) fibers: $(changed)")
    end
end
logmsg("\nINTELLIGENT-TFUNLISTS-OK  ", now())
close(report_io)
