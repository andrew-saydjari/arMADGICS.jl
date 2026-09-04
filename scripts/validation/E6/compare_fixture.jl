# E6 step 4b: compare the M123 fixture run with the OLD (production rough) vs NEW
# (E6 fiber-295 build) starCont prior. Reads the per-fiber results h5 and the
# batchstyle component dumps.
using HDF5, Statistics, Printf

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E6_run"
fold = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(RUND, "fixture_oldprior.h5")
fnew = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(RUND, "fixture_newprior.h5")
println("baseline: $fold\nswapped:  $fnew")

rd(f) = h5open(f, "r") do fh
    Dict(k => read(fh[k]) for k in keys(fh))
end
o = rd(fold); n = rd(fnew)
fibs = o["fiberindx"]
println("=== E6 runtime swap comparison (old=production rough | new=E6 f295 build) ===")
println(rpad("fib", 5), rpad("status o/n", 12), rpad("chi2res_old", 14), rpad("chi2res_new", 14),
    rpad("dchi2/old", 12), rpad("RV_old", 10), rpad("RV_new", 10), rpad("dRV(pix)", 10),
    rpad("ss0 o/n-1", 12), rpad("flag o/n", 9))
fails = String[]
for (i, fib) in enumerate(fibs)
    so, sn = o["status"][i], n["status"][i]
    so == sn || push!(fails, "f$fib status $so vs $sn")
    co, cn = o["chi2res"][i], n["chi2res"][i]
    rvo, rvn = o["RV_pixoff_final"][i], n["RV_pixoff_final"][i]
    flo, fln = o["RV_flag"][i], n["RV_flag"][i]
    dchi = (cn - co) / abs(co)
    ss = n["starscale0"][i] / o["starscale0"][i] - 1
    @printf("%-5d%-12s%-14.6g%-14.6g%-12.3e%-10.4g%-10.4g%-10.4g%-12.3e%-9s\n",
        fib, "$so/$sn", co, cn, isfinite(dchi) ? dchi : NaN, rvo, rvn, rvn - rvo, ss, "$flo/$fln")
    # prior-independent quantities must be identical
    for k in ("snr_pipe", "ingestBit", "skyBit", "nSkyFibers")
        vo, vn = o[k][i], n[k][i]
        (isequal(vo, vn)) || push!(fails, "f$fib $k differs: $vo vs $vn")
    end
    if so == "ok"
        flo == fln || push!(fails, "f$fib RV_flag differs: $flo vs $fln")
        (!isfinite(dchi) || abs(dchi) < 0.05) || push!(fails, "f$fib chi2res rel change $(dchi)")
        (isnan(rvo) && isnan(rvn)) || abs(rvn - rvo) <= 0.2 || push!(fails, "f$fib RV moved $(rvn-rvo) pix")
    end
end

# component-level comparison from the batchstyle dumps
bo = replace(fold, ".h5" => "_batchstyle.h5")
bn = replace(fnew, ".h5" => "_batchstyle.h5")
if isfile(bo) && isfile(bn)
    co = rd(bo); cn = rd(bn)
    for k in ("x_starContinuum_v0", "x_skyContinuum_v0", "x_starLines_v0", "apVisit_v0")
        A, B = co[k], cn[k]
        d = filter(isfinite, vec(B .- A))
        a = filter(isfinite, vec(A))
        rel = maximum(abs.(d)) / max(1e-30, quantile(abs.(a), 0.99))
        @printf("component %-22s max|Δ|/q99|old| = %.3e   rms Δ = %.3e\n", k, rel, sqrt(mean(d .^ 2)))
        if k == "x_starContinuum_v0"
            # pixel-scale artifact check: high-pass the difference per fiber; a healthy
            # continuum change is smooth, so the pixel-to-pixel diff of Δ should be tiny
            for j in axes(A, 2)
                dd = B[:, j] .- A[:, j]
                ddf = filter(isfinite, diff(dd))
                if !isempty(ddf)
                    hp = sqrt(mean(ddf .^ 2))
                    sm = sqrt(mean(filter(isfinite, dd) .^ 2))
                    if sm > 0 && hp / sm > 0.5 && sm > 1e-3 * quantile(abs.(a), 0.99)
                        push!(fails, "pixel-scale structure in Δ starContinuum fixture col $j (hp/sm=$(hp/sm))")
                    end
                end
            end
        end
    end
else
    println("batchstyle dumps missing — component comparison skipped")
end

println()
if isempty(fails)
    println("RUNTIME-SWAP-PASS")
else
    println("RUNTIME-SWAP-FAIL:")
    foreach(println, fails)
end
