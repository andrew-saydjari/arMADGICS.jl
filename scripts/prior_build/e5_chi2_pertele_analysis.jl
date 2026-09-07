## E5: QUANTIFY ONLY — what changes if the chi2 outlier screen's p99.9 thresholds are
## computed PER TELESCOPE instead of pooled across both (AKS 2026-09-07: "should always
## have been per telescope").
#
# Nothing is rebuilt and no decision file is rewritten. This reads the existing
# screens/e5_sample_stats.h5 and reports, for pooled vs per-telescope thresholds:
#   (1) the thresholds themselves;
#   (2) additional samples screened, per telescope, absolute and %, overall and per fiber;
#   (3) the propagated ASSEMBLED count per TARGET fiber (pool = in-MTP-block +/-5, same
#       telescope) vs the GSPICE floor (~1015-1022);
#   (4) concentration diagnostics (MTP block / MJD era) for a qualitative-surprise check.
#
# The chi2 screen fires when chi2r_cont > CHI2CONT_MAX OR chi2r_full > CHI2FULL_MAX; the
# bright screen (robust z of log10 scale within the SOURCE fiber) is already per-fiber and
# so already per-telescope — it is held fixed here and reported for context only.
# Usage: julia --project=<arM-E5b> e5_chi2_pertele_analysis.jl
# Author - Andrew Saydjari (E5 pass 1)
using HDF5, Statistics, Printf, Dates

const OUT = get(ENV, "E5_OUT", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1")
const SC = joinpath(OUT, "screens")
const RESDIR = get(ENV, "E5_RESDIR", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_07/e5_combined")
mkpath(RESDIR)

const SCALE_ZMAX = 6.0
tele_of(adjfib) = adjfib <= 300 ? :apo : :lco
mtp_block(fiberindx::Int) = fld(fiberindx - 1, 30) + 1
native(adjfib) = adjfib <= 300 ? adjfib : adjfib - 300
function pool_range(fiberindx::Int; wind::Int=5)
    blo = 30 * fld(fiberindx - 1, 30) + 1
    return max(fiberindx - wind, blo):min(fiberindx + wind, blo + 29)
end

# ---------------------------------------------------------------- load sample stats
fibs = Int[]
per_fib = Dict{Int,NamedTuple}()
h5open(joinpath(SC, "e5_sample_stats.h5"), "r") do fh
    for k in sort(collect(keys(fh)))
        adjfib = parse(Int, k); g = fh[k]
        scale = read(g["scale"]); c2f = read(g["chi2r_full"]); c2c = read(g["chi2r_cont"])
        mjd = read(g["mjd"]); expnum = read(g["expnum"])
        n = length(scale)
        z = fill(NaN, n); pos = scale .> 0
        if count(pos) >= 10
            lx = log10.(scale[pos]); med = median(lx)
            madv = max(1.4826 * median(abs.(lx .- med)), 1e-3)
            z[pos] .= (lx .- med) ./ madv
        end
        per_fib[adjfib] = (scale=scale, z=z, c2f=c2f, c2c=c2c, mjd=mjd, expnum=expnum)
        push!(fibs, adjfib)
    end
end
ntot = sum(length(per_fib[f].scale) for f in fibs)

# pooled vs per-telescope value pools
gather(sel, field) = reduce(vcat, [filter(isfinite, getfield(per_fib[f], field)) for f in fibs if sel(f)])
all_c2f = gather(_ -> true, :c2f); all_c2c = gather(_ -> true, :c2c)
apo_c2f = gather(f -> f <= 300, :c2f); apo_c2c = gather(f -> f <= 300, :c2c)
lco_c2f = gather(f -> f > 300, :c2f); lco_c2c = gather(f -> f > 300, :c2c)

chi2_threshold(v; p_hi=0.999, p_ref=0.99, floor_fac=3.0) =
    max(quantile(v, p_hi), floor_fac * quantile(v, p_ref))

POOL_FULL = chi2_threshold(all_c2f); POOL_CONT = chi2_threshold(all_c2c)
APO_FULL = chi2_threshold(apo_c2f);  APO_CONT = chi2_threshold(apo_c2c)
LCO_FULL = chi2_threshold(lco_c2f);  LCO_CONT = chi2_threshold(lco_c2c)

# which arm of max() bound each threshold (p99.9 vs the 3x-p99 floor)?
armof(v) = quantile(v, 0.999) >= 3 * quantile(v, 0.99) ? "p99.9" : "3xp99 FLOOR"

io = open(joinpath(RESDIR, "chi2_pertele_report.txt"), "w")
# NB: @printf needs a literal format string, so route through Printf.format at runtime
function emit(fmt, args...)
    s = Printf.format(Printf.Format(fmt), args...)
    print(io, s); print(stdout, s); flush(io); flush(stdout)
end

emit("E5 chi2 screen: POOLED vs PER-TELESCOPE p99.9 thresholds (%s)\n", string(now()))
emit("ANALYSIS ONLY — nothing rebuilt, no decision file rewritten.\n\n")
emit("samples: %d total (%d APO, %d LCO across %d source fibers)\n\n",
    ntot, sum(length(per_fib[f].scale) for f in fibs if f <= 300),
    sum(length(per_fib[f].scale) for f in fibs if f > 300), length(fibs))

emit("## 1. Thresholds (MEASURED)\n\n")
emit("%-28s %14s %14s %14s %14s\n", "quantity", "p99", "p99.9", "3x p99", "THRESHOLD")
for (nm, v, thr) in (("chi2r_cont POOLED", all_c2c, POOL_CONT), ("chi2r_cont APO", apo_c2c, APO_CONT),
    ("chi2r_cont LCO", lco_c2c, LCO_CONT), ("chi2r_full POOLED", all_c2f, POOL_FULL),
    ("chi2r_full APO", apo_c2f, APO_FULL), ("chi2r_full LCO", lco_c2f, LCO_FULL))
    emit("%-28s %14.3f %14.3f %14.3f %14.3f  [%s]\n", nm, quantile(v, 0.99),
        quantile(v, 0.999), 3 * quantile(v, 0.99), thr, armof(v))
end
emit("\nmedian chi2r_cont: APO %.3f  LCO %.3f  (ratio %.2fx)\n",
    median(apo_c2c), median(lco_c2c), median(apo_c2c) / median(lco_c2c))
emit("median chi2r_full: APO %.3f  LCO %.3f  (ratio %.2fx)\n\n",
    median(apo_c2f), median(lco_c2f), median(apo_c2f) / median(lco_c2f))

# ---------------------------------------------------------------- drop decisions
"drop mask under a given (cont, full) threshold pair per telescope; bright screen fixed"
function drops(f, cont_thr, full_thr)
    st = per_fib[f]; n = length(st.scale)
    br = falses(n); c2 = falses(n)
    for k in 1:n
        (isfinite(st.z[k]) && st.z[k] > SCALE_ZMAX) && (br[k] = true)
        ((isfinite(st.c2c[k]) && st.c2c[k] > cont_thr) ||
         (isfinite(st.c2f[k]) && st.c2f[k] > full_thr)) && (c2[k] = true)
    end
    return br, c2
end
thr_pooled(f) = (POOL_CONT, POOL_FULL)
thr_pertele(f) = f <= 300 ? (APO_CONT, APO_FULL) : (LCO_CONT, LCO_FULL)

kept_pooled = Dict{Int,BitVector}(); kept_pertele = Dict{Int,BitVector}()
for f in fibs
    b1, c1 = drops(f, thr_pooled(f)...); b2, c2 = drops(f, thr_pertele(f)...)
    kept_pooled[f] = .!(b1 .| c1); kept_pertele[f] = .!(b2 .| c2)
end

function tallies(keptd, sel)
    ns = sum(length(per_fib[f].scale) for f in fibs if sel(f))
    nd = sum(count(.!keptd[f]) for f in fibs if sel(f))
    return nd, ns, 100nd / ns
end
emit("## 2. Screened fractions (MEASURED)\n\n")
emit("%-12s %14s %14s %14s\n", "telescope", "pooled", "per-telescope", "additional")
for (nm, sel) in (("APO", f -> f <= 300), ("LCO", f -> f > 300), ("BOTH", _ -> true))
    d1, ns, p1 = tallies(kept_pooled, sel); d2, _, p2 = tallies(kept_pertele, sel)
    emit("%-12s %6d (%.3f%%) %6d (%.3f%%) %6d (%+.3f pp)\n", nm, d1, p1, d2, p2, d2 - d1, p2 - p1)
end

# ---------------------------------------------------------------- per-source-fiber
emit("\n## 3. Worst-affected SOURCE fibers (extra samples screened per-tele vs pooled)\n\n")
extra = [(f, count(.!kept_pertele[f]) - count(.!kept_pooled[f]), length(per_fib[f].scale)) for f in fibs]
sort!(extra, by=x -> -x[2])
emit("%-8s %-6s %10s %10s %10s\n", "adjfib", "tele", "n_samp", "extra", "extra%")
for (f, e, n) in extra[1:min(15, end)]
    emit("%-8d %-6s %10d %10d %9.2f%%\n", f, string(tele_of(f)), n, e, n > 0 ? 100e / n : 0.0)
end
srcloss = [(f, n > 0 ? 100 * (count(.!kept_pertele[f]) - count(.!kept_pooled[f])) / n : 0.0)
           for (f, _, n) in extra]
emit("\nworst SOURCE-fiber relative extra loss: %.2f%% (adjfib %d)\n",
    maximum(x -> x[2], srcloss), srcloss[argmax([x[2] for x in srcloss])][1])

# ---------------------------------------------------------------- assembled counts
emit("\n## 4. Propagated ASSEMBLED counts per TARGET fiber (GSPICE floor ~1015-1022)\n\n")
function assembled(keptd)
    a = zeros(Int, 600)
    for t in 1:600
        fi = native(t); off = t <= 300 ? 0 : 300
        a[t] = sum(count(keptd[g + off]) for g in pool_range(fi) if haskey(keptd, g + off))
    end
    return a
end
A1 = assembled(kept_pooled); A2 = assembled(kept_pertele)
emit("pooled screening:      min=%d med=%.1f max=%d\n", minimum(A1), median(A1), maximum(A1))
emit("per-telescope screening: min=%d med=%.1f max=%d\n", minimum(A2), median(A2), maximum(A2))
rel = [A1[t] > 0 ? 100 * (A1[t] - A2[t]) / A1[t] : 0.0 for t in 1:600]
emit("worst TARGET relative loss: %.3f%% (fiber %d: %d -> %d)\n",
    maximum(rel), argmax(rel), A1[argmax(rel)], A2[argmax(rel)])
emit("targets below 3000 assembled: %d   below 1022 (GSPICE floor): %d\n",
    count(<(3000), A2), count(<(1022), A2))
emit("min assembled / GSPICE floor = %.1fx\n", minimum(A2) / 1022)
srt = sortperm(A2)[1:8]
emit("\nlowest-8 assembled targets under per-telescope screening:\n")
for t in srt
    emit("  f%03d (%s): pooled %d -> per-tele %d  (%.3f%% loss, %.1fx floor)\n",
        t, string(tele_of(t)), A1[t], A2[t], rel[t], A2[t] / 1022)
end

# ---------------------------------------------------------------- concentration
emit("\n## 5. Concentration diagnostics for the EXTRA drops (surprise check)\n\n")
xf = Int[]; xm = Int[]
for f in fibs
    d = kept_pooled[f] .& .!kept_pertele[f]     # newly dropped only
    for k in findall(d); push!(xf, f); push!(xm, per_fib[f].mjd[k]); end
end
emit("extra drops: %d total (APO %d, LCO %d)\n", length(xf),
    count(<=(300), xf), count(>(300), xf))
if !isempty(xf)
    blkc = Dict{Tuple{Symbol,Int},Int}()
    for f in xf; k = (tele_of(f), mtp_block(native(f))); blkc[k] = get(blkc, k, 0) + 1; end
    emit("MTP blocks touched: %d of 20 (10 per telescope)\n", length(blkc))
    top = sort(collect(blkc), by=x -> -x[2])[1:min(6, end)]
    emit("top MTP blocks: %s\n", join(["$(t[1][1])-blk$(t[1][2]):$(t[2])" for t in top], "  "))
    emit("  -> largest single block holds %.1f%% of extra drops\n", 100 * top[1][2] / length(xf))
    mjc = Dict{Int,Int}(); for m in xm; mjc[m] = get(mjc, m, 0) + 1; end
    tm = sort(collect(mjc), by=x -> -x[2])[1:min(8, end)]
    emit("distinct MJDs: %d;  top MJDs: %s\n", length(mjc),
        join(["$(t[1]):$(t[2])" for t in tm], "  "))
    emit("  -> largest single MJD holds %.1f%% of extra drops\n", 100 * tm[1][2] / length(xf))
    emit("distinct source fibers touched: %d of %d\n", length(unique(xf)), length(fibs))
end

# ---------------------------------------------------------------- RE-ADMISSIONS
# The change is NOT one-directional: a per-telescope p99.9 is LOOSER on APO (whose chi2
# distribution is inflated) as well as TIGHTER on LCO. Samples currently screened that
# would be let back in are the risk side of this change, because the screen's designed
# target is the E4 light-leak anomaly class, which arrives as whole bad EXPOSURES.
emit("\n## 6. RE-ADMITTED samples (currently screened, kept under per-telescope) \n\n")
rf = Int[]; rm = Int[]; re_ = Int[]
for f in fibs
    d = .!kept_pooled[f] .& kept_pertele[f]
    for k in findall(d)
        push!(rf, f); push!(rm, per_fib[f].mjd[k]); push!(re_, per_fib[f].expnum[k])
    end
end
emit("re-admitted: %d total (APO %d, LCO %d)\n", length(rf), count(<=(300), rf), count(>(300), rf))
if !isempty(rf)
    emit("distinct source fibers touched: %d;  distinct MJDs: %d;  distinct (mjd,exp): %d\n",
        length(unique(rf)), length(unique(rm)), length(unique(collect(zip(rm, re_)))))
    ec = Dict{Tuple{Int,Int},Int}()
    for (m, e) in zip(rm, re_); ec[(m, e)] = get(ec, (m, e), 0) + 1; end
    te = sort(collect(ec), by=x -> -x[2])[1:min(10, end)]
    emit("top (mjd,expnum) by re-admitted fiber count:\n")
    for (k, v) in te
        emit("   mjd %d exp %-5d : %4d fibers re-admitted\n", k[1], k[2], v)
    end
    emit("  -> largest single EXPOSURE holds %.1f%% of re-admissions; top-10 hold %.1f%%\n",
        100 * te[1][2] / length(rf), 100 * sum(x -> x[2], te) / length(rf))
    nbig = count(x -> x[2] >= 50, collect(ec))
    emit("exposures with >=50 fibers re-admitted: %d (that shape = a whole anomalous exposure,\n", nbig)
    emit("   which is exactly the class the chi2 screen was designed to remove)\n")
end

# ---------------------------------------------------------------- variant for AKS
# NOT AKS's instruction and NOT adopted here — reported only, because it is the obvious
# way to get the LCO tightening he asked for without the APO re-admission side effect.
emit("\n## 7. VARIANT (agent suggestion, NOT adopted, NOT AKS's words):\n")
emit("   per-telescope threshold applied TIGHTENING-ONLY, thr = min(pooled, per-telescope)\n\n")
tighten(f) = f <= 300 ? (min(POOL_CONT, APO_CONT), min(POOL_FULL, APO_FULL)) :
             (min(POOL_CONT, LCO_CONT), min(POOL_FULL, LCO_FULL))
kept_tight = Dict{Int,BitVector}()
for f in fibs
    b, c = drops(f, tighten(f)...); kept_tight[f] = .!(b .| c)
end
emit("%-12s %14s %14s %14s\n", "telescope", "pooled", "tighten-only", "additional")
for (nm, sel) in (("APO", f -> f <= 300), ("LCO", f -> f > 300), ("BOTH", _ -> true))
    d1, ns, p1 = tallies(kept_pooled, sel); d3, _, p3 = tallies(kept_tight, sel)
    emit("%-12s %6d (%.3f%%) %6d (%.3f%%) %6d (%+.3f pp)\n", nm, d1, p1, d3, p3, d3 - d1, p3 - p1)
end
nre = sum(count(.!kept_pooled[f] .& kept_tight[f]) for f in fibs)
A3 = assembled(kept_tight)
emit("re-admitted under this variant: %d (by construction 0)\n", nre)
emit("assembled: min=%d med=%.1f max=%d  (pooled min was %d)\n",
    minimum(A3), median(A3), maximum(A3), minimum(A1))
emit("targets below 3000: %d;  below GSPICE floor 1022: %d;  min/floor = %.1fx\n",
    count(<(3000), A3), count(<(1022), A3), minimum(A3) / 1022)
close(io)
println("\nreport -> ", joinpath(RESDIR, "chi2_pertele_report.txt"))
