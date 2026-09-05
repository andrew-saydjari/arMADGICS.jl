# C2 policy analysis REVISION 2 (AKS correction 2026-09-05):
#  - The APO 10k per-entry ceiling is NOT physics-based (previous framing wrong).
#  - APO has the same FPI/bright-fiber bimodality as LCO (fibers 76/226 are the
#    APO analogs of LCO 88/219); the retained ceiling was again dropping a ton
#    of entries from them.
# Evaluate AKS's designs:
#  (a) exposure-level flagging on NON-FPI fleet statistics (both telescopes);
#  (b) fiberID-dependent thresholds -> implemented as a per-fiber entry-level
#      gap screen (a bright MODE is not gap-isolated, so bright fibers are never
#      penalized; only egregious isolated per-fiber outlier entries flag).
# Produces the revised policy table + APO delta costs. ANALYSIS ONLY.
using HDF5, Statistics, StatsBase, Printf, Serialization, Dates

const OUTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/c2_analysis"
io = open(joinpath(OUTD, "C2_POLICY_REPORT_v2.txt"), "w")
logmsg(args...) = (println(args...); println(io, args...); flush(stdout); flush(io))
logmsg("C2 policy analysis v2 (corrected framing) — ", now())

const TELL = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_apo.h5",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5")
const AUD = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/tfunlist_audit_apo.h5",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/tfunlists_lcoC2/tfunlist_audit_lco.h5")
const LEAK = [2578, 3830]
const REFLIST = Dict(   # CURRENT pass-1b consumption lists
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/20260902_apo_tfunlist.jdat",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/tfunlists_intelligent/20260905_lco_tfunlist.jdat")

data = Dict{String,Any}()
for tele in ("apo", "lco")
    mf, p999 = h5open(AUD[tele], "r") do f
        read(f["medflux"]), read(f["chi_p999"])
    end
    chif = h5open(TELL[tele], "r") do f
        Float64.(read(f["chi_sq_fiber"]))
    end
    base = (mf .> 400.0) .& (chif .<= p999)
    data[tele] = (mf=mf, chif=chif, base=base, N=size(mf, 2))
end

# ── 1. Bimodality characterization (objective metric, both telescopes) ─────────────────
logmsg("\n== Per-fiber bright-mode (bimodality) characterization ==")
logmsg("metric: fbright = fraction of base-kept entries with medflux > 2x the fiber's own median")
FPI = Dict{String,Vector{Int}}()
for tele in ("apo", "lco")
    mf, base = data[tele].mf, data[tele].base
    fbright = zeros(300); q90r = zeros(300)
    for fb in 1:300
        v = mf[fb, base[fb, :]]
        m = median(v)
        fbright[fb] = count(>(2m), v) / length(v)
        q90r[fb] = quantile(v, 0.9) / m
    end
    thr = 0.05   # >5% of entries in a >2x bright mode = bimodal bright fiber
    set = findall(>(thr), fbright)
    FPI[tele] = set
    logmsg(@sprintf("[%s] fibers with fbright > 5%%: %s", tele, string(set)))
    for fb in sort(set, by=fb -> -fbright[fb])
        v = mf[fb, base[fb, :]]
        logmsg(@sprintf("[%s]   fiberindx %3d: fbright %.1f%%  q50 %.0f  q90 %.0f  q99 %.0f  max %.0f  (q90/q50 = %.1f)",
            tele, fb, 100fbright[fb], median(v), quantile(v, .9), quantile(v, .99), maximum(v), q90r[fb]))
    end
    # context: next-highest fbright values (continuum fibers)
    rest = sort(fbright[setdiff(1:300, set)], rev=true)[1:5]
    logmsg(@sprintf("[%s]   next-highest fbright among remaining fibers: %s", tele,
        join([@sprintf("%.2f%%", 100x) for x in rest], " ")))
    serialize(joinpath(OUTD, "fbright_$(tele).jdat"), (fbright=fbright, set=set))
end

# ── 2. E1a: exposure gap screen on FPI-EXCLUDED fleet statistics ───────────────────────
logmsg("\n== E1a: gap-isolated exposure screen, exposure medians over NON-FPI fibers only ==")
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
E1a = Dict{String,Vector{Int}}()
for tele in ("apo", "lco")
    mf = data[tele].mf; N = data[tele].N
    nonfpi = setdiff(1:300, FPI[tele])
    em = [median(@view mf[nonfpi, i]) for i in 1:N]
    flagged = gap_flag(em)
    E1a[tele] = flagged
    m, s = median(em), mad(em, normalize=true)
    top = sort(em, rev=true)[1:6]
    logmsg(@sprintf("[%s] non-FPI exposure medians: med %.0f MAD %.0f; top-6: %s; E1a flags %d: %s",
        tele, m, s, join(round.(Int, top), " "), length(flagged), string(flagged)))
    for e in flagged
        logmsg(@sprintf("[%s]   flagged row %d: non-FPI median %.0f (%.0f MADs)", tele, e, em[e], (em[e] - m) / s))
    end
end

# ── 3. F1: per-fiber entry-level gap screen (fiberID-dependent, design (b)) ────────────
logmsg("\n== F1: per-fiber gap-isolated entry screen (G=2, cap 0.1% of the fiber's kept entries) ==")
F1 = Dict{String,Any}()
for tele in ("apo", "lco")
    mf, base = data[tele].mf, data[tele].base
    flagged = [Int[] for _ in 1:300]
    for fb in 1:300
        idx = findall(@view base[fb, :])
        v = mf[fb, idx]
        ord = sortperm(v, rev=true)
        cap = max(1, ceil(Int, 0.001 * length(idx)))
        istar = 0
        for i in 1:cap
            (v[ord[i]] / v[ord[i+1]] >= 2.0) && (istar = i)
        end
        istar > 0 && (flagged[fb] = idx[ord[1:istar]])
    end
    F1[tele] = flagged
    nf = count(!isempty, flagged); tot = sum(length, flagged)
    logmsg(@sprintf("[%s] F1 flags entries on %d/300 fibers; total %d entries; max per fiber %d",
        tele, nf, tot, maximum(length.(flagged))))
    for fb in findall(!isempty, flagged)
        v = mf[fb, flagged[fb]]
        vb = mf[fb, data[tele].base[fb, :]]
        logmsg(@sprintf("[%s]   fiberindx %3d: %d flagged (medflux %s; fiber q99 %.0f)",
            tele, fb, length(flagged[fb]), join(round.(Int, sort(v, rev=true)), " "), quantile(vb, .99)))
    end
end

# chi2 character of the APO >10k entries (is the old ceiling catching anomalies?)
let (mf, chif, base) = (data["apo"].mf, data["apo"].chif, data["apo"].base)
    sel10k = base .& (mf .> 10_000.0)
    logmsg(@sprintf("\n[apo] chi2 of base-kept >10k entries: median %.1f vs all base-kept %.1f (p99.9 gate %.0f) — %s",
        median(chif[sel10k]), median(chif[base]), 58328.0,
        median(chif[sel10k]) < 2 * median(chif[base]) ? "NOT anomalous as a class" : "elevated"))
end

# ── 4. Policy variants + deltas vs CURRENT pass-1b lists ───────────────────────────────
logmsg("\n== Policy variants (base = C1&C3) ==")
function build_keep(tele; E1flag=Int[], F1flag=nothing, ceil10k=false)
    mf, base = data[tele].mf, data[tele].base
    keep = copy(base)
    for e in E1flag; keep[:, e] .= false; end
    if F1flag !== nothing
        for fb in 1:300, e in F1flag[fb]; keep[fb, e] = false; end
    end
    ceil10k && (keep .&= (mf .<= 10_000.0))
    keep
end
function delta_vs_current(tele, keep)
    mf = data[tele].mf
    ref = deserialize(REFLIST[tele])
    newlists = [sort(findall(@view keep[fb, :]), by=i -> mf[fb, i], rev=true) for fb in 1:300]
    changed = [fb for fb in 1:300 if newlists[fb] != ref[fb]]
    changed, newlists
end
variants = [
    ("status-quo (E1 + APO 10k ceiling)", Dict("apo" => (E1a["apo"], nothing, true), "lco" => ([LEAK...], nothing, false))),
    ("E1a only (ceiling dropped)", Dict("apo" => (E1a["apo"], nothing, false), "lco" => (E1a["lco"], nothing, false))),
    ("REVISED: E1a + F1 (no ceilings)", Dict("apo" => (E1a["apo"], F1["apo"], false), "lco" => (E1a["lco"], F1["lco"], false))),
]
for (nm, cfg) in variants
    logmsg("\n--- $nm ---")
    for tele in ("apo", "lco")
        E1f, F1f, c10 = cfg[tele]
        keep = build_keep(tele; E1flag=E1f, F1flag=F1f, ceil10k=c10)
        mf, base = data[tele].mf, data[tele].base
        removed = count(base) - count(keep)
        perfib = [count(@view keep[fb, :]) for fb in 1:300]
        fpirem = sum(count(@view base[fb, :]) - count(@view keep[fb, :]) for fb in FPI[tele])
        leaksurv = tele == "lco" ? sum(count(keep[:, e]) for e in LEAK) : -1
        changed, _ = delta_vs_current(tele, keep)
        logmsg(@sprintf("[%s] removed %d; FPI/bright-fiber removed %d; kept min %d med %.0f; %s; DELTA vs current list: %d fibers change",
            tele, removed, fpirem, minimum(perfib), median(perfib),
            tele == "lco" ? "leak surviving $(leaksurv)" : "(apo)", length(changed)))
        if tele == "apo" && length(changed) > 0 && length(changed) <= 300
            logmsg(@sprintf("[apo]   resample+rebuild cost if adopted: %d fibers x ~4.5 min sample + ~5 min build (=~%.1f h at 16 workers)",
                length(changed), length(changed) * 9.5 / 16 / 60))
        end
    end
end
serialize(joinpath(OUTD, "c2_policy_v2_results.jdat"), (FPI=FPI, E1a=E1a, F1=F1))
logmsg("\nsaved c2_policy_v2_results.jdat  ", now())
close(io)
