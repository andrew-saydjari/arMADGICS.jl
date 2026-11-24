# 02_build_visits.jl -- compute the 21 allVisitRV columns (visit membership + reductions).
#
# CHAIN POSITION: allVisitRV step 2 of 4 (after 01_seqmap.jl, before 03_write_allvisitrv.jl).
# INPUTS:  --allexposure  the allExposure table (needs tele, mjd, expnum, sdss_id, mjd_mid_exposure,
#                         dithered_pixels, v_rad, e_v_rad, snr, qa_bitmask, ingestBit, RV_flag,
#                         bary_flag, ra, dec, glon, glat, adjfiberindx, era, field_id, plate_id,
#                         config_id, category)
#          --seqmap       seqmap.jls from step 1
#          --reference    OPTIONAL previously released allVisitRV; if given, every carried column must
#                         equal it row for row or the script refuses to write (the 0.4.0 alignment proof)
# OUTPUTS: --out          v2_cols.jls (Serialization NamedTuple of the column vectors)
# SOURCE:  work/2026_09_16/allvisitrv/s9_jack_build.jl (0.4.0). Logic unchanged apart from: paths
#          became the arguments above; the alignment proof became optional (--reference); the
#          scratch reference CSV (v_rad_scatter_plain_vs_jack.csv.gz) was removed; Jackknife.jl is a
#          project dependency instead of an include() from a user depot.
#
# Rebuild allVisitRV (21 datasets), with `v_rad_scatter_jack`, the bias-corrected JACKKNIFE
# scatter, ALONGSIDE the plain sample standard deviation `v_rad_scatter`, which is kept because
# the jackknife is undefined at n_exp == 2 (decision 2026-09-16: keep both columns).
#
# ---------------------------------------------------------------------------------------------
# THE ESTIMATOR  (AJ 169, 167 -- arXiv:2408.07126, Section 4, Eqs. 7-8-10)
#
#   Eq 10  the statistic that is jackknifed, the ordinary sample standard deviation
#
#             sigma_hat = sqrt( 1/(N-1) * sum_i ( x_i - (1/N) sum_j x_j )^2 )
#
#   Eq 7   mean of that statistic over the N leave-one-out subsamples
#
#             theta_tilde_jack = (1/N) * sum_i f(x_1, ..., x_{i-1}, x_{i+1}, ..., x_N)
#
#   Eq 8   the bias-corrected jackknife estimate -- THIS is the delivered number
#
#             theta_hat_jack = N * theta_hat - (N-1) * theta_tilde_jack
#
# i.e.  v_rad_scatter_jack = N*sigma_hat(x) - (N-1) * mean_i sigma_hat(x without x_i)
#
# The paper's own words (Section 4): "we can calculate the spread sigma(vbary) with the standard
# deviation of vbary using jackknife resampling to help offset biases introduced by the small
# sample sizes ... using a simple standard deviation without jackknifing leads to 17% and 11%
# underestimation of the scatter for sample size N = 3 and N = 4, respectively."
#
# The reference code path is `jack_std` in ApogeeReduction.jl/src/utils.jl:
#
#     function jack_std(x)
#         y = filter(!isnanorzero, x)
#         ly = length(y)
#         if ly > 2
#             return Jackknife.estimate(std, y)
#         else
#             return Inf
#         end
#     end
#
# and Jackknife.jl's `estimate(f, x) = n*f(x) - (n-1)*mean(leaveoneout(f, x))` is literally Eq 8
# with f = Statistics.std (the (N-1)-normalised, mean-centred std of Eq 10).  Paper and code agree.
#
# N < 3 GUARD: the leave-one-out subsamples of a 2-element sample are singletons, whose (N-1)-
# normalised std is 0/0.  The estimator is UNDEFINED at N = 2.  The reference code returns Inf there;
# this file writes NaN instead (same meaning, "undefined", but the sentinel already used by the
# rest of this table for an undefined float).  70.5 % of the visits here have n_exp == 2, so
# 70.5 % of `v_rad_scatter_jack` is NaN -- which is exactly why the plain `v_rad_scatter` is
# retained as a 21st column rather than replaced.  See allVisitRV_README.md.
# ---------------------------------------------------------------------------------------------
#
# Membership is recomputed from scratch with the SAME code path as s2_build.jl / s6_addcols.jl
# (almanac object sequence x sdss_id, release default mask, balanced-dither-only, same ordering)
# and then every carried column is required to equal the live released file row for row before
# anything is written.  READ-ONLY outside this directory.
using HDF5, Printf, Statistics, Serialization, Base.Threads, ArgParse
using Jackknife     # the reference jackknife code path (Jackknife.estimate), used by selftest()

function parse_cli()
    s = ArgParseSettings(description = "compute allVisitRV columns")
    @add_arg_table! s begin
        "--allexposure"; required = true
        "--seqmap"; required = true
        "--reference"; default = ""
        "--out"; required = true
    end
    parse_args(s)
end
const ARGV = parse_cli()
const P    = ARGV["allexposure"]
const LIVE = ARGV["reference"]
say(s) = (println(s); flush(stdout))

# ---- jackknifed scatter, Eq 8 applied to Eq 10 -----------------------------------------------
# Two-pass mean/sum-of-squares per subsample: numerically identical to Statistics.std, and stable
# (no sum-of-squares cancellation).  n <= 44 here, so the O(n^2) form is free.
@inline function subsample_std(y::Vector{Float64}, skip::Int)
    n = length(y)
    m = skip == 0 ? n : n - 1
    mu = 0.0
    @inbounds for j in 1:n
        j == skip && continue
        mu += y[j]
    end
    mu /= m
    ss = 0.0
    @inbounds for j in 1:n
        j == skip && continue
        d = y[j] - mu
        ss += d * d
    end
    return sqrt(ss / (m - 1))
end

function jack_scatter(y::Vector{Float64})
    n = length(y)
    n >= 3 || return NaN                      # the reference `ly > 2` guard (it returns Inf; we use NaN)
    acc = 0.0
    @inbounds for i in 1:n
        acc += subsample_std(y, i)
    end
    return n * subsample_std(y, 0) - (n - 1) * (acc / n)   # Eq 8
end

# ---- prove our implementation IS Jackknife.estimate(std, .) ----------------------------------
function selftest()
    # deterministic, wide dynamic range (no RNG dependency)
    cases = Vector{Vector{Float64}}()
    for n in [3, 4, 5, 6, 8, 13, 44], s in [1e-5, 1e-3, 1.0, 1e2], off in [0.0, -273.15, 412.7]
        push!(cases, [off + s * sin(2.3i + 0.7n) + s * 0.3cos(7.1i) for i in 1:n])
    end
    worst = 0.0
    for y in cases
        a = jack_scatter(y)
        b = Jackknife.estimate(std, y)
        worst = max(worst, abs(a - b) / max(abs(b), eps()))
    end
    say(@sprintf("self-test vs Jackknife.estimate(std, .) on %d samples: worst rel. diff = %.3e", length(cases), worst))
    worst < 1e-10 || error("jackknife implementation does not match Jackknife.jl")
    # n = 2 really is undefined in his convention
    y2 = [1.0, 2.0]
    say("   n=2: ours = $(jack_scatter(y2)) ; Jackknife.estimate(std, [1.0,2.0]) = $(Jackknife.estimate(std, y2)) ; jack_std guard returns Inf")
end

function main()
    t0 = time()
    selftest()

    seqmap = deserialize(ARGV["seqmap"])
    lut = Dict{Tuple{Int8,Int32}, Tuple{Vector{Int32},Vector{Int32},Vector{Int32}}}()
    gid = Int32(0)
    for k in sort(collect(keys(seqmap)))
        v = seqmap[k]
        st = Int32[a for (a, b) in v]; en = Int32[b for (a, b) in v]
        o = sortperm(st); st = st[o]; en = en[o]
        lut[k] = (st, en, Int32[(gid += Int32(1)) for _ in st])
    end

    f = h5open(P, "r")
    tele = read(f["tele"]); mjd = read(f["mjd"]); expnum = read(f["expnum"])
    sid = read(f["sdss_id"]); mjdmid = read(f["mjd_mid_exposure"])
    dith = read(f["dithered_pixels"]); vr = read(f["v_rad"]); ev = read(f["e_v_rad"])
    snr = read(f["snr"])
    qa = read(f["qa_bitmask"]); ing = read(f["ingestBit"]); rvf = read(f["RV_flag"]); bf = read(f["bary_flag"])
    ra = read(f["ra"]); dec = read(f["dec"]); glon = read(f["glon"]); glat = read(f["glat"])
    afi = read(f["adjfiberindx"]); era = read(f["era"]); fld = read(f["field_id"])
    plate = read(f["plate_id"]); config = read(f["config_id"]); cat = read(f["category"])
    close(f)
    N = length(mjd)
    say("allExposure rows = $N ($(round(time()-t0, digits=1)) s)")

    seqid = zeros(Int32, N)
    @threads for i in 1:N
        e = get(lut, (tele[i], mjd[i]), nothing); e === nothing && continue
        st, en, ids = e; j = searchsortedlast(st, expnum[i])
        (j >= 1 && expnum[i] <= en[j]) && (seqid[i] = ids[j])
    end
    mask = (qa .== 0) .& (ing .== 0) .& (rvf .== 0) .& (bf .== 0) .& isfinite.(vr)
    dlab = Vector{Float32}(undef, N)
    @threads for i in 1:N
        dlab[i] = (!isfinite(dith[i]) || dith[i] == 0f0) ? Float32(NaN) : dith[i]
    end

    valid = (seqid .> 0) .& (sid .> 0)
    idx = findall(valid)
    key = [(seqid[i], sid[i]) for i in idx]
    ord = sortperm(key); idx = idx[ord]; key = key[ord]
    NR = length(idx)
    bnd = Int[1]
    for k in 2:NR
        key[k] == key[k-1] || push!(bnd, k)
    end
    push!(bnd, NR + 1)
    NV = length(bnd) - 1
    say("visits (pre-mask) = $NV")

    O_sid = Int64[]; O_tele = Int8[]; O_mjd = Int32[]; O_mid = Float64[]; O_v = Float64[]
    O_n = Int16[]; O_snrm = Float64[]; O_snrq = Float64[]; O_sc = Float64[]; O_jk = Float64[]
    O_ev = Float64[]
    O_ra = Float64[]; O_dec = Float64[]; O_glon = Float64[]; O_glat = Float64[]
    O_fib = Int16[]; O_era = Int8[]; O_fld = Int32[]; O_plate = Int32[]; O_config = Int32[]
    O_cat = UInt8[]
    buf = Int[]; y = Float64[]
    for b in 1:NV
        rg = bnd[b]:(bnd[b+1] - 1)
        empty!(buf)
        for k in rg
            r = idx[k]; mask[r] && push!(buf, r)
        end
        isempty(buf) && continue
        any(r -> isnan(dlab[r]), buf) && continue
        vals = unique(dlab[r] for r in buf)
        length(vals) == 2 || continue
        c1 = count(r -> dlab[r] == vals[1], buf)
        c1 == length(buf) - c1 || continue

        r0 = buf[1]; n = length(buf)
        resize!(y, n)
        @inbounds for t in 1:n
            y[t] = Float64(vr[buf[t]])
        end

        push!(O_sid, sid[r0]); push!(O_tele, tele[r0]); push!(O_mjd, mjd[r0])
        push!(O_mid, mean(mjdmid[r] for r in buf)); push!(O_v, mean(vr[r] for r in buf))
        push!(O_n, Int16(n))
        push!(O_snrm, mean(snr[r] for r in buf))
        push!(O_snrq, sqrt(sum(Float64(snr[r])^2 for r in buf)))
        push!(O_sc, n >= 2 ? std(y) : NaN)
        push!(O_jk, jack_scatter(y))
        push!(O_ev, mean(Float64(ev[r]) for r in buf) / sqrt(n))
        push!(O_ra, ra[r0]); push!(O_dec, dec[r0]); push!(O_glon, glon[r0]); push!(O_glat, glat[r0])
        push!(O_fib, Int16(Int(afi[r0]) - 300 * (Int(tele[r0]) - 1))); push!(O_era, Int8(era[r0]))
        push!(O_fld, fld[r0]); push!(O_plate, plate[r0]); push!(O_config, config[r0])
        push!(O_cat, cat[r0])
    end
    say("balanced visits = $(length(O_sid))  ($(round(time()-t0, digits=1)) s)")

    o = sortperm(collect(zip(O_tele, O_mjd, O_sid, O_mid)))
    O_sid = O_sid[o]; O_tele = O_tele[o]; O_mjd = O_mjd[o]; O_mid = O_mid[o]; O_v = O_v[o]
    O_n = O_n[o]; O_snrm = O_snrm[o]; O_snrq = O_snrq[o]; O_sc = O_sc[o]; O_jk = O_jk[o]
    O_ev = O_ev[o]; O_ra = O_ra[o]; O_dec = O_dec[o]; O_glon = O_glon[o]; O_glat = O_glat[o]
    O_fib = O_fib[o]; O_era = O_era[o]; O_fld = O_fld[o]; O_plate = O_plate[o]
    O_config = O_config[o]; O_cat = O_cat[o]
    telestr = [t == 1 ? "apo" : "lco" for t in O_tele]

    # ---- ALIGNMENT PROOF (optional): every carried column must equal a reference file, row for row
    if !isempty(LIVE)
        h = h5open(LIVE, "r")
        live = Dict{String,Any}()
        for k in keys(h); live[k] = read(h[k]); end
        close(h)
        mine = Dict{String,Any}("sdss_id"=>O_sid, "mjd"=>O_mjd, "mjd_mid"=>O_mid,
            "v_rad_mean"=>O_v, "n_exp"=>O_n, "snr_mean"=>O_snrm, "snr_visit"=>O_snrq,
            "v_rad_scatter"=>O_sc, "e_v_rad_visit"=>O_ev, "ra"=>O_ra, "dec"=>O_dec,
            "glon"=>O_glon, "glat"=>O_glat, "fiber"=>O_fib, "era"=>O_era, "field_id"=>O_fld,
            "plate_id"=>O_plate, "config_id"=>O_config, "category"=>O_cat,
            "tele"=>telestr)
        say("")
        say("--- value-for-value match against the live released file (carried datasets) ---")
        allok = true
        for k in sort([k for k in keys(live) if haskey(mine, k)])   # carried columns only
            a = live[k]; b = mine[k]
            k == "tele" && (a = String.(a))
            ok = length(a) == length(b) && all(a .== b)
            allok &= ok
            @printf("   %-15s n=%d  identical=%s\n", k, length(b), ok)
        end
        say("ALL CARRIED COLUMNS REPRODUCED EXACTLY: $allok")
        allok || error("alignment failed -- refusing to write")
    end

    # ---- statistics of the new column vs the old one
    n2 = count(==(Int16(2)), O_n)
    fin = findall(isfinite, O_jk)
    say("")
    say("v_rad_scatter_jack: rows = $(length(O_jk)), finite = $(length(fin)), NaN = $(length(O_jk)-length(fin))")
    say("   NaN rows are exactly the n_exp == 2 rows: $(length(O_jk) - length(fin) == n2)")
    say("   negative values (jackknife can undershoot below 0): $(count(<(0), O_jk[fin]))")
    qs = [0.01, 0.05, 0.16, 0.25, 0.50, 0.75, 0.84, 0.95, 0.99]
    jk = sort(O_jk[fin]); pl_same = sort(O_sc[fin]); pl_all = sort(filter(isfinite, O_sc))
    say(@sprintf("%-8s %12s %12s %12s", "pctile", "jack(n>=3)", "plain(n>=3)", "plain(all)"))
    for q in qs
        @printf("%-8.0f %12.6g %12.6g %12.6g\n", 100q, quantile(jk, q), quantile(pl_same, q), quantile(pl_all, q))
    end
    say(@sprintf("ratio of medians jack/plain on the same n>=3 rows = %.4f", quantile(jk, 0.5) / quantile(pl_same, 0.5)))
    for nn in [4, 6, 8, 10, 12]
        sel = findall(==(Int16(nn)), O_n)
        isempty(sel) && continue
        a = median(O_jk[sel]); b = median(O_sc[sel])
        @printf("   n_exp=%2d : median jack = %.6g  median plain = %.6g  ratio = %.4f  (%d rows)\n", nn, a, b, a/b, length(sel))
    end

    serialize(ARGV["out"],
        (; O_sid, O_tele, O_mjd, O_mid, O_v, O_n, O_snrm, O_snrq, O_sc, O_jk, O_ev, O_ra, O_dec,
           O_glon, O_glat, O_fib, O_era, O_fld, O_plate, O_config, O_cat, telestr))
    say("wrote $(ARGV["out"])  ($(round(time()-t0, digits=1)) s)")
end

main()
