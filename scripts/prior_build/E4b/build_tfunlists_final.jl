# FINAL tfunlist cut policy (AKS-APPROVED 2026-09-05: "I agree with your choices
# on the fpi, finish it"). Executes the REVISED UNIFIED POLICY of
# 2026_09_05/c2_analysis/C2_POLICY_REPORT_v2.txt with the open choices resolved:
#   K = 20k ("cliff"), F1 ON, APO C3 p99.9 RECOMPUTED on the post-E2 fleet.
# Policy (both telescopes):
#   C1: medflux > 400 (unchanged faint cut)
#   C3: chi_sq_fiber <= per-telescope p99.9. APO gate RECOMPUTED over the
#       post-E2 fleet (the 38 catastrophic exposures inflate the full-fleet
#       percentile); LCO gate stays the full-fleet value (the approved recompute
#       is APO-specific — E2 removes only the 2 leak domeflats at LCO and the
#       approved outcome is LCO byte-identical to the intelligent list).
#   E2: exposure screen on NON-FPI fleet statistics (FPI = fibers with
#       fbright > 5%: APO {76,226}, LCO {88,219}), egregious-only:
#       (i)  brightness gap-isolation (G=2, cap 0.1%) on non-FPI exposure-median
#            medflux  -> LCO {2578,3830} leak domeflats; APO none (smooth).
#       (ii) fit-quality cliff: non-FPI exposure-median chi2 > K=20,000
#            -> APO 38 catastrophic exposures (chi2stat 31k-54.6k vs fleet
#            median 15.5); LCO none.
#   F1: per-fiber gap-isolated entry screen (G=2, cap 0.1% of the fiber's
#       C1&C3-kept entries) -> 3 isolated APO freak entries (f28 @277k,
#       f90 @312k, f226 @80k medflux); at LCO only leak entries (subset of
#       E2(i) columns -> no additional effect, asserted).
#   NO absolute ceilings (the APO 10k ceiling is STRUCK — it suppressed 3,311
#   pristine FPI bright-mode entries, chi2 med 6), NO per-fiber trims.
# F1 is computed on the C1&C3 base (v2-report semantics) and unioned with E2;
# final keep = C1 & C3 & !E2 & !F1.
# Expected: LCO byte-identical to tfunlists_intelligent; APO restores ~3,546
# pristine FPI-fiber entries and removes ~40 (38 exposures + freaks), ~293+
# fibers change (more where the recomputed gate bites).
#
# Run: nice -n 10 julia +1.11.6 --project=<arM-E4b> -t 8 build_tfunlists_final.jl

using LinearAlgebra, Statistics, StatsBase, HDF5, Serialization, Printf, Dates
BLAS.set_num_threads(1)

const OUTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_final"
const TELL = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_apo.h5",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5")
# medflux reused from the intelligent-build audits (identical medflux_of code;
# cross-asserted below against the independent 2026_09_03 / 2026_09_04 audits)
const AUD = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_intelligent/tfunlist_audit_apo.h5",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_intelligent/tfunlist_audit_lco.h5")
const AUD_XCHECK = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/tfunlist_audit_apo.h5",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/tfunlists_lcoC2/tfunlist_audit_lco.h5")
const REF = Dict(   # CURRENT pass-1b consumption lists
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/20260902_apo_tfunlist.jdat",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_intelligent/20260905_lco_tfunlist.jdat")
const INTEL_LCO = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_intelligent/20260905_lco_tfunlist.jdat"
const K_CLIFF = 20_000.0
const EXPECT_FPI = Dict("apo" => [76, 226], "lco" => [88, 219])
mkpath(OUTD)

report_io = open(joinpath(OUTD, "BUILD_REPORT.txt"), "a")
logmsg(args...) = (println(args...); println(report_io, args...); flush(stdout); flush(report_io))
logmsg("FINAL tfunlist build (AKS-approved policy, 2026-09-05) — ", now())
logmsg("julia $(VERSION), nthreads=$(Threads.nthreads())")
logmsg("policy: C1 + C3(APO p99.9 recomputed post-E2; LCO full-fleet) + E2 non-FPI exposure screen [(i) brightness gap G=2 cap 0.1%; (ii) chi2 cliff K=20k] + F1 per-fiber gap screen; NO ceilings, NO per-fiber trims")

function gap_flag(em; G=2.0, capfrac=0.001)
    N = length(em)
    ord = sortperm(em, rev=true)
    cap = max(1, ceil(Int, capfrac * N))
    istar = 0
    for i in 1:cap
        (em[ord[i]] / em[ord[i+1]] >= G) && (istar = i)
    end
    istar == 0 ? Int[] : sort(ord[1:istar])
end

for tele in ("apo", "lco")
    t0 = time()
    mf = h5open(AUD[tele], "r") do f
        read(f["medflux"])
    end
    mfx = h5open(AUD_XCHECK[tele], "r") do f
        read(f["medflux"])
    end
    @assert mf == mfx "medflux mismatch between intelligent audit and $(AUD_XCHECK[tele])"
    chif = h5open(TELL[tele], "r") do f
        Float64.(read(f["chi_sq_fiber"]))
    end
    N = size(mf, 2)
    p999_full = percentile(vec(chif), 99.9)

    # FPI (bimodal bright-fiber) determination: fbright > 5% on the C1&C3(full) base
    basefull = (mf .> 400.0) .& (chif .<= p999_full)
    fbright = zeros(300)
    for fb in 1:300
        v = mf[fb, basefull[fb, :]]
        fbright[fb] = count(>(2 * median(v)), v) / length(v)
    end
    FPIset = findall(>(0.05), fbright)
    @assert FPIset == EXPECT_FPI[tele] "FPI set $(FPIset) != expected $(EXPECT_FPI[tele])"
    nonfpi = setdiff(1:300, FPIset)
    logmsg(@sprintf("\n[%s] N=%d exposures; FPI fibers (fbright>5%%): %s (fbright %s)", tele, N,
        string(FPIset), join([@sprintf("%.1f%%", 100fbright[fb]) for fb in FPIset], " ")))

    # E2 clause (i): brightness gap-isolation on non-FPI exposure medians
    em = [median(@view mf[nonfpi, i]) for i in 1:N]
    E2i = gap_flag(em)
    m_em, s_em = median(em), mad(em, normalize=true)
    logmsg(@sprintf("[%s] E2(i) brightness gap screen (non-FPI medians: med %.0f MAD %.0f): flags %d: %s",
        tele, m_em, s_em, length(E2i), string(E2i)))
    for e in E2i
        logmsg(@sprintf("[%s]   E2(i) row %d: non-FPI median medflux %.0f (%.0f MADs)", tele, e, em[e], (em[e] - m_em) / s_em))
    end

    # E2 clause (ii): fit-quality cliff on non-FPI exposure-median chi2
    ec = [median(@view chif[nonfpi, i]) for i in 1:N]
    E2ii = findall(>(K_CLIFF), ec)
    m_ec, s_ec = median(ec), mad(ec, normalize=true)
    logmsg(@sprintf("[%s] E2(ii) chi2 cliff K=%.0f (non-FPI chi2stat: med %.1f MAD %.1f): flags %d exposures%s",
        tele, K_CLIFF, m_ec, s_ec, length(E2ii),
        isempty(E2ii) ? "" : @sprintf("; chi2stat range %.0f-%.0f", minimum(ec[E2ii]), maximum(ec[E2ii]))))
    E2 = sort(union(E2i, E2ii))
    logmsg(@sprintf("[%s] E2 union: %d exposures", tele, length(E2)))

    # C3 gate: APO recomputed post-E2; LCO full-fleet (approved LCO-unchanged outcome)
    keepcols = setdiff(1:N, E2)
    p999_post = percentile(vec(chif[:, keepcols]), 99.9)
    p999 = tele == "apo" ? p999_post : p999_full
    logmsg(@sprintf("[%s] C3 gate: full-fleet p99.9 = %.4f; post-E2 p99.9 = %.4f; USED = %.4f (%s)",
        tele, p999_full, p999_post, p999,
        tele == "apo" ? "RECOMPUTED post-E2 per approved policy" : "full-fleet; recompute is APO-specific"))

    fail1 = mf .<= 400.0
    fail3 = .!(chif .<= p999)
    base = .!(fail1 .| fail3)
    if tele == "apo"
        dC3 = count(fail3) - count(.!(chif .<= p999_full))
        logmsg(@sprintf("[apo] C3 delta from gate recompute: %d additional entries fail C3 (%d -> %d)",
            dC3, count(.!(chif .<= p999_full)), count(fail3)))
    end

    # F1: per-fiber gap-isolated entry screen on the C1&C3 base (v2 semantics)
    failF = falses(300, N)
    nF1 = 0
    for fb in 1:300
        idx = findall(@view base[fb, :])
        v = mf[fb, idx]
        ord = sortperm(v, rev=true)
        cap = max(1, ceil(Int, 0.001 * length(idx)))
        istar = 0
        for i in 1:cap
            (v[ord[i]] / v[ord[i+1]] >= 2.0) && (istar = i)
        end
        if istar > 0
            for j in idx[ord[1:istar]]
                failF[fb, j] = true
                nF1 += 1
            end
        end
    end
    F1outside = [(fb, e) for fb in 1:300, e in 1:N if failF[fb, e] && !(e in E2)]
    logmsg(@sprintf("[%s] F1 flags %d entries total; %d OUTSIDE E2 exposures (these are the effective F1 removals):",
        tele, nF1, length(F1outside)))
    for (fb, e) in F1outside
        logmsg(@sprintf("[%s]   F1 fiberindx %3d exp row %d: medflux %.0f chi2 %.1f", tele, fb, e, mf[fb, e], chif[fb, e]))
    end
    if tele == "lco"
        @assert isempty(F1outside) "LCO F1 flags entries outside E2 — breaks approved LCO-unchanged outcome"
    end

    failE = falses(300, N)
    for e in E2
        failE[:, e] .= true
    end
    keep = base .& .!failE .& .!failF

    logmsg(@sprintf("[%s] C1 (mf<=400)        excludes %7d", tele, count(fail1)))
    logmsg(@sprintf("[%s] C3 (chi2>p99.9)     excludes %7d", tele, count(fail3)))
    logmsg(@sprintf("[%s] E2 exposures        excludes %6d entries (%d exposures) beyond C1/C3", tele,
        count(base .& failE), length(E2)))
    logmsg(@sprintf("[%s] F1 freak entries    excludes %6d beyond C1/C3/E2", tele, count(base .& .!failE .& failF)))
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
    for fb in FPIset
        logmsg(@sprintf("[%s] FPI fiber %d kept %d (MUST be ~full list)", tele, fb, perfib[fb]))
    end

    outlist = joinpath(OUTD, "20260905_$(tele)_tfunlist.jdat")
    serialize(outlist, Tfunlists)
    h5open(joinpath(OUTD, "tfunlist_audit_$(tele).h5"), "w") do f
        f["medflux"] = mf
        f["chi_p999_full"] = p999_full
        f["chi_p999_postE2"] = p999_post
        f["chi_p999_used"] = p999
        f["fpi_fibers"] = FPIset
        f["expo_median_nonfpi"] = em
        f["expo_chi2_nonfpi"] = ec
        f["E2i_flagged"] = E2i
        f["E2ii_flagged"] = E2ii
        f["F1_fiber"] = Int[fb for (fb, e) in Tuple.(findall(failF))]
        f["F1_expo"] = Int[e for (fb, e) in Tuple.(findall(failF))]
        f["kept_mask"] = UInt8.(keep)
        f["perfiber_kept"] = perfib
    end
    logmsg(@sprintf("[%s] wrote %s (%d B) in %.0fs", tele, outlist, filesize(outlist), time() - t0))

    # ── delta accounting vs the CURRENT pass-1b consumption list ─────────────
    ref = deserialize(REF[tele])
    changed = [fb for fb in 1:300 if Tfunlists[fb] != ref[fb]]
    logmsg(@sprintf("[%s] vs CURRENT pass-1b list: identical fibers %d/300; changed fibers %d",
        tele, 300 - length(changed), length(changed)))
    if tele == "apo"
        # class accounting vs the current APO keep = C1 & C3(full) & 10k ceiling
        cur = basefull .& (mf .<= 10_000.0)
        restored = keep .& .!cur
        removed = cur .& .!keep
        rFPI = sum(count(@view restored[fb, :]) for fb in FPIset)
        r76 = count(@view restored[76, :]); r226 = count(@view restored[226, :])
        remE2 = count(removed .& failE)
        remF1 = count(removed .& failF .& .!failE)
        remC3 = count(removed .& fail3)
        logmsg(@sprintf("[apo] RESTORED %d entries (ceiling struck): f76 %d, f226 %d, other %d; of these %d have mf>10k",
            count(restored), r76, r226, count(restored) - rFPI, count(restored .& (mf .> 10_000.0))))
        logmsg(@sprintf("[apo] REMOVED %d entries: %d in E2 exposures, %d F1 freaks, %d by tightened C3 gate",
            count(removed), remE2, remF1, remC3))
        logmsg(@sprintf("[apo] FPI kept counts: f76 %d -> %d (%+d); f226 %d -> %d (%+d)",
            count(@view cur[76, :]), perfib[76], perfib[76] - count(@view cur[76, :]),
            count(@view cur[226, :]), perfib[226], perfib[226] - count(@view cur[226, :])))
        logmsg("[apo] changed fibers: ", string(changed))
        serialize(joinpath(OUTD, "apo_changed_fibers.jdat"), changed)
    else
        # byte-identity vs the intelligent list (approved outcome)
        ident = read(outlist) == read(INTEL_LCO)
        logmsg("[lco] BYTE-IDENTITY vs tfunlists_intelligent 20260905_lco_tfunlist.jdat: ",
            ident ? "IDENTICAL (approved outcome confirmed)" : "MISMATCH — POLICY VIOLATION")
        @assert ident
        @assert isempty(changed)
    end
end
logmsg("\nFINAL-TFUNLISTS-OK  ", now())
close(report_io)
