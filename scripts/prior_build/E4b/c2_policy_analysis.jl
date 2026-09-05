# C2 bright-cut policy analysis (AKS 2026-09-05): the committed per-entry
# C2_LCO=3000 removes 4,569 entries (229 leak + 4,340 collateral on the FPI
# fibers) — too big. Evaluate candidate policies on BOTH telescopes against:
#   (a) removes the 2 leaked LCO domeflats' entries (rows 2578/3830);
#   (b) removes/controls the APO medflux>10,000 nonlinearity regime;
#   (c) collateral on the FPI/bright fibers (LCO fiberindx 88/219 [+148/159];
#       APO 28/226/115) — should be SMALL;
#   (d) simplicity/robustness.
# Policies (all on top of base = C1 (mf>400) & C3 (chi2<=p99.9)):
#   const-X  : & mf <= X                     (X swept)
#   rank99   : drop the 99 brightest base-kept entries per fiber
#              (modernized 2026_04_25 convention; the old N-100 faint cap is
#              reported but not part of the candidate)
#   hybrid   : rank99 + safety constant (LCO 8000 / APO 12000)
# Pure analysis — no tfunlist is written.
using HDF5, Statistics, StatsBase, Printf, Serialization, Dates

const OUTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/c2_analysis"
mkpath(OUTD)
io = open(joinpath(OUTD, "C2_POLICY_REPORT.txt"), "w")
logmsg(args...) = (println(args...); println(io, args...); flush(stdout); flush(io))
logmsg("C2 bright-cut policy analysis — ", now())

const TELL = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_apo.h5",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5")
const AUD = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/tfunlist_audit_apo.h5",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/tfunlists_lcoC2/tfunlist_audit_lco.h5")
const LEAK = [2578, 3830]              # lco leaked domeflats (row indices)
const SPECIAL = Dict("lco" => [88, 219, 148, 159], "apo" => [28, 226, 115])
const APO_REGIME = 10_000.0            # physically-motivated nonlinearity threshold

data = Dict{String,Any}()
for tele in ("apo", "lco")
    mf, p999 = h5open(AUD[tele], "r") do f
        read(f["medflux"]), read(f["chi_p999"])
    end
    chif = h5open(TELL[tele], "r") do f
        Float64.(read(f["chi_sq_fiber"]))
    end
    base = (mf .> 400.0) .& (chif .<= p999)     # C1 & C3 only
    data[tele] = (mf=mf, base=base, p999=p999, N=size(mf, 2))
    logmsg(@sprintf("[%s] N=%d exposures; base(C1&C3)-kept %d entries (per-fiber med %.0f)",
        tele, size(mf, 2), count(base), median([count(@view base[fb, :]) for fb in 1:300])))
end

# ── leak-entry brightness anatomy (LCO) ────────────────────────────────────────────────
let (mf, base) = (data["lco"].mf, data["lco"].base)
    for e in LEAK
        kept = findall(@view base[:, e])
        v = mf[kept, e]
        logmsg(@sprintf("[lco] leak exp %d: %d base-kept fibers; medflux min %.0f med %.0f max %.0f",
            e, length(kept), minimum(v), median(v), maximum(v)))
    end
    allleak = vcat([mf[base[:, e], e] for e in LEAK]...)
    logmsg(@sprintf("[lco] ALL base-kept leak entries: n=%d, min=%.0f  => a constant cut must sit BELOW %.0f to remove every leak entry",
        length(allleak), minimum(allleak), minimum(allleak)))
    # FPI fibers' own brightness for contrast
    for fb in SPECIAL["lco"]
        v = mf[fb, base[fb, :]]
        logmsg(@sprintf("[lco] fiberindx %3d (FPI/bright): base-kept %d; q50 %.0f q90 %.0f q99 %.0f max %.0f; n>3000: %d  n>5000: %d  n>6000: %d  n>8000: %d",
            fb, length(v), quantile(v, .5), quantile(v, .9), quantile(v, .99), maximum(v),
            count(>(3000), v), count(>(5000), v), count(>(6000), v), count(>(8000), v)))
    end
end
# APO regime anatomy
let (mf, base) = (data["apo"].mf, data["apo"].base)
    n10k = [count(mf[fb, base[fb, :]] .> APO_REGIME) for fb in 1:300]
    logmsg(@sprintf("[apo] base-kept entries >10k per fiber: min %d med %.0f max %d (fiber %d); fibers with >99: %d",
        minimum(n10k), median(n10k), maximum(n10k), argmax(n10k), count(>(99), n10k)))
    for fb in SPECIAL["apo"]
        v = mf[fb, base[fb, :]]
        logmsg(@sprintf("[apo] fiberindx %3d (bright): base-kept %d; q50 %.0f q90 %.0f q99 %.0f max %.0f; n>10k: %d n>12k: %d",
            fb, length(v), quantile(v, .5), quantile(v, .9), quantile(v, .99), maximum(v),
            count(>(10_000), v), count(>(12_000), v)))
    end
end

# ── policy machinery ───────────────────────────────────────────────────────────────────
# returns keep::BitMatrix given policy
function apply_policy(tele, kind; X=Inf, rank=0)
    mf, base = data[tele].mf, data[tele].base
    keep = copy(base)
    if kind == :const || kind == :hybrid
        keep .&= (mf .<= X)
    end
    if kind == :rank || kind == :hybrid
        for fb in 1:300
            idx = findall(@view keep[fb, :])
            if length(idx) > rank
                bright = partialsortperm(mf[fb, idx], 1:rank, rev=true)
                keep[fb, idx[bright]] .= false
            end
        end
    end
    keep
end

function eval_policy(tele, name, keep)
    mf, base = data[tele].mf, data[tele].base
    removed = base .& .!keep
    perfib_kept = [count(@view keep[fb, :]) for fb in 1:300]
    out = Dict{String,Any}("name" => name, "tele" => tele,
        "removed_total" => count(removed),
        "kept_min" => minimum(perfib_kept), "kept_med" => median(perfib_kept),
        "special_removed" => sum(count(@view removed[fb, :]) for fb in SPECIAL[tele]))
    if tele == "lco"
        leak_surv = sum(count(keep[:, e]) for e in LEAK)
        # per-fiber check: any fiber where a leak entry survives?
        fib_bad = [fb for fb in 1:300 if any(keep[fb, e] for e in LEAK)]
        out["leak_surviving"] = leak_surv
        out["leak_fibers_bad"] = length(fib_bad)
        out["collateral"] = count(removed) - (sum(count(base[:, e]) for e in LEAK) - leak_surv)
    else
        out["regime_surviving"] = count(keep .& (mf .> APO_REGIME))
        out["collateral"] = count(removed .& (mf .<= APO_REGIME))
    end
    out
end

results = []
# LCO menu (const 4000-8000 expected to FAIL criterion (a): min leak medflux 3987)
for X in (3000.0, 3500.0, 3900.0, 4000.0, 5000.0, 6000.0, 8000.0)
    push!(results, eval_policy("lco", "const-$(Int(X))", apply_policy("lco", :const, X=X)))
end
push!(results, eval_policy("lco", "rank99", apply_policy("lco", :rank, rank=99)))
push!(results, eval_policy("lco", "rank99+6000", apply_policy("lco", :hybrid, X=6000.0, rank=99)))
push!(results, eval_policy("lco", "rank99+8000", apply_policy("lco", :hybrid, X=8000.0, rank=99)))
# APO menu
for X in (10_000.0, 12_000.0)
    push!(results, eval_policy("apo", "const-$(Int(X))", apply_policy("apo", :const, X=X)))
end
push!(results, eval_policy("apo", "rank99", apply_policy("apo", :rank, rank=99)))
push!(results, eval_policy("apo", "rank99+12000", apply_policy("apo", :hybrid, X=12_000.0, rank=99)))

logmsg("\n=== POLICY TABLE (LCO) — base = C1&C3; criteria: leak gone / collateral small ===")
logmsg(@sprintf("%-14s %14s %12s %12s %16s %10s %10s", "policy", "leak surviving", "removed", "collateral", "FPI-fib removed", "kept min", "kept med"))
for r in results
    r["tele"] == "lco" || continue
    logmsg(@sprintf("%-14s %8d (%d fib) %12d %12d %16d %10d %10.0f",
        r["name"], r["leak_surviving"], r["leak_fibers_bad"], r["removed_total"],
        r["collateral"], r["special_removed"], r["kept_min"], r["kept_med"]))
end
logmsg("\n=== POLICY TABLE (APO) — criteria: >10k regime suppressed / collateral small ===")
logmsg(@sprintf("%-14s %15s %12s %12s %18s %10s %10s", "policy", ">10k surviving", "removed", "collateral", "brightfib removed", "kept min", "kept med"))
for r in results
    r["tele"] == "apo" || continue
    logmsg(@sprintf("%-14s %15d %12d %12d %18d %10d %10.0f",
        r["name"], r["regime_surviving"], r["removed_total"],
        r["collateral"], r["special_removed"], r["kept_min"], r["kept_med"]))
end

# rank99 per-fiber leak check detail (criterion a for the rank policy)
let keep = apply_policy("lco", :rank, rank=99)
    mf, base = data["lco"].mf, data["lco"].base
    maxrank = 0
    for fb in 1:300, e in LEAK
        base[fb, e] || continue
        idx = findall(@view base[fb, :])
        r = count(mf[fb, idx] .> mf[fb, e]) + 1   # brightness rank of the leak entry
        maxrank = max(maxrank, r)
    end
    logmsg(@sprintf("\n[lco] rank99 leak check: worst brightness rank of any base-kept leak entry = %d (must be <= 99 for rank99 to remove all)", maxrank))
end

# fine LCO constant sweep for the figure
sweepX = collect(3000.0:250.0:9000.0)
sweep = map(sweepX) do X
    r = eval_policy("lco", "s", apply_policy("lco", :const, X=X))
    (X, r["leak_surviving"], r["collateral"], r["special_removed"])
end
serialize(joinpath(OUTD, "c2_policy_results.jdat"),
    (results=results, sweep=sweep, sweepX=sweepX))
logmsg("\nsaved c2_policy_results.jdat  ", now())
close(io)
