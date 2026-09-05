# C2 bright-cut policy analysis figures (dark theme). Companion to
# c2_policy_analysis.jl; publishes to 2026_09_05/plots/c2_cut_analysis/.
# Run: cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_05 && nice -n 10 \
#   julia +1.11.6 --project=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfun_viz \
#   /mnt/home/asaydjari/gitcode/worktrees/arM-E4b/scripts/prior_build/E4b/make_figures_c2.jl
using CairoMakie, ColorSchemes, Printf, Serialization, HDF5, Statistics, StatsBase

black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")

const OUTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/c2_analysis"
const PLOTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/plots/c2_cut_analysis"
mkpath(PLOTD)

const TELL = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_apo.h5",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/telluric_refit_full/out/tellurics_refit_20260902_lco.h5")
const AUD = Dict(
    "apo" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfunlists_refit20260902/tfunlist_audit_apo.h5",
    "lco" => "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/tfunlists_lcoC2/tfunlist_audit_lco.h5")
const LEAK = [2578, 3830]

function load_tele(tele)
    mf, p999 = h5open(AUD[tele], "r") do f
        read(f["medflux"]), read(f["chi_p999"])
    end
    chif = h5open(TELL[tele], "r") do f
        Float64.(read(f["chi_sq_fiber"]))
    end
    mf, (mf .> 400.0) .& (chif .<= p999)
end
mf_lco, base_lco = load_tele("lco")
mf_apo, base_apo = load_tele("apo")
res = deserialize(joinpath(OUTD, "c2_policy_results.jdat"))

qb(mf, base, q) = [length(v) > 0 ? quantile(v, q) : NaN for v in
    (mf[fb, base[fb, :]] for fb in 1:300)]

## LCO fig 1: per-fiber quantile bands + leak entries
begin
    fig = Figure(size=(1150, 460))
    ax = Axis(fig[1, 1], yscale=log10, xlabel="fiberindx", ylabel="medflux",
        title="LCO per-fiber base-kept medflux (C1&C3): q50/q90/q99/max — FPI fibers 88/219 (+148/159); leak entries as dots")
    for (q, c, lw) in ((0.5, :deepskyblue, 1.2), (0.9, :cyan, 1.0), (0.99, :orange, 1.0))
        lines!(ax, 1:300, qb(mf_lco, base_lco, q), color=c, linewidth=lw, label="q$(Int(100q))")
    end
    mx = [maximum(mf_lco[fb, base_lco[fb, :]]) for fb in 1:300]
    scatter!(ax, 1:300, mx, color=:gray70, markersize=3, label="max")
    for e in LEAK
        fibs = findall(base_lco[:, e])
        scatter!(ax, fibs, mf_lco[fibs, e], color=:orangered, markersize=5,
            label=e == LEAK[1] ? "leak entries" : nothing)
    end
    hlines!(ax, [3000], color=:magenta, linestyle=:dash, label="C2=3000 (committed)")
    hlines!(ax, [3987], color=:springgreen, linestyle=:dot, label="min leak = 3987")
    for fb in (88, 148, 159, 219)
        vlines!(ax, [fb], color=(:white, 0.15))
        text!(ax, fb, 1.2e4, text="$fb", color=:white, fontsize=11)
    end
    axislegend(ax, position=:lt, framevisible=false, nbanks=4)
    save(joinpath(PLOTD, "fig_lco_quantiles.png"), fig, px_per_unit=2)
end

## LCO fig 2: histograms
begin
    fig = Figure(size=(1150, 640))
    picks = ((88, "fiberindx 88 (FPI)"), (219, "fiberindx 219 (FPI)"),
        (148, "fiberindx 148"), (150, "fiberindx 150 (typical, = adjfib 450)"))
    for (i, (fb, ttl)) in enumerate(picks)
        r, c = fldmod1(i, 2)
        ax = Axis(fig[r, c], xscale=log10, xlabel=r == 2 ? "medflux" : "", ylabel="entries",
            title=ttl)
        v = mf_lco[fb, base_lco[fb, :]]
        hist!(ax, v, bins=10 .^ range(log10(401), log10(15000), 80), color=(:deepskyblue, 0.8))
        vlines!(ax, [3000], color=:magenta, linestyle=:dash)
        vlines!(ax, [3987], color=:springgreen, linestyle=:dot)
        for e in LEAK
            base_lco[fb, e] && vlines!(ax, [mf_lco[fb, e]], color=:orangered, linewidth=2)
        end
        r99 = partialsort(v, 99, rev=true)
        vlines!(ax, [r99], color=:yellow, linestyle=:dashdot)
    end
    Label(fig[0, :], "LCO medflux distributions — magenta: C2=3000; green: min-leak 3987; red: this fiber's leak entries; yellow: rank-99 threshold", fontsize=14)
    save(joinpath(PLOTD, "fig_lco_hists.png"), fig, px_per_unit=2)
end

## LCO fig 3: constant sweep
begin
    X = res.sweepX
    leak = [s[2] for s in res.sweep]; fpi = [s[4] for s in res.sweep]
    fig = Figure(size=(1000, 420))
    ax = Axis(fig[1, 1], xlabel="constant cut C2_LCO", ylabel="entries (log)", yscale=log10,
        title="LCO constant-cut sweep: no X both removes all 229 leak entries AND spares the FPI fibers")
    lines!(ax, X, max.(leak, 0.5), color=:orangered, linewidth=2.5, label="leak entries SURVIVING")
    lines!(ax, X, max.(fpi, 0.5), color=:magenta, linewidth=2.5, label="FPI-fiber entries REMOVED (collateral)")
    hlines!(ax, [396], color=:yellow, linestyle=:dashdot, label="rank99 FPI collateral (396)")
    vlines!(ax, [3987], color=:springgreen, linestyle=:dot, label="min leak = 3987")
    axislegend(ax, position=:ct, framevisible=false)
    save(joinpath(PLOTD, "fig_lco_sweep.png"), fig, px_per_unit=2)
end

## APO fig 1: quantile bands
begin
    fig = Figure(size=(1150, 460))
    ax = Axis(fig[1, 1], yscale=log10, xlabel="fiberindx", ylabel="medflux",
        title="APO per-fiber base-kept medflux (C1&C3): q50/q90/q99/max — bright fibers 28/226/115 and 76")
    for (q, c, lw) in ((0.5, :deepskyblue, 1.2), (0.9, :cyan, 1.0), (0.99, :orange, 1.0))
        lines!(ax, 1:300, qb(mf_apo, base_apo, q), color=c, linewidth=lw, label="q$(Int(100q))")
    end
    mx = [maximum(mf_apo[fb, base_apo[fb, :]]) for fb in 1:300]
    scatter!(ax, 1:300, mx, color=:gray70, markersize=3, label="max")
    hlines!(ax, [10_000], color=:magenta, linestyle=:dash, label="C2=10,000 (nonlinearity)")
    for fb in (28, 76, 115, 226)
        vlines!(ax, [fb], color=(:white, 0.15))
        text!(ax, fb, 1.5e5, text="$fb", color=:white, fontsize=11)
    end
    axislegend(ax, position=:lt, framevisible=false, nbanks=3)
    save(joinpath(PLOTD, "fig_apo_quantiles.png"), fig, px_per_unit=2)
end

## APO fig 2: histograms
begin
    fig = Figure(size=(1150, 640))
    picks = ((28, "fiberindx 28 (λ1 outlier; max 277k)"), (226, "fiberindx 226 (1,542 entries >10k)"),
        (76, "fiberindx 76 (1,764 entries >10k)"), (150, "fiberindx 150 (typical)"))
    for (i, (fb, ttl)) in enumerate(picks)
        r, c = fldmod1(i, 2)
        ax = Axis(fig[r, c], xscale=log10, xlabel=r == 2 ? "medflux" : "", ylabel="entries", title=ttl)
        v = mf_apo[fb, base_apo[fb, :]]
        hist!(ax, v, bins=10 .^ range(log10(401), log10(3e5), 90), color=(:deepskyblue, 0.8))
        vlines!(ax, [10_000], color=:magenta, linestyle=:dash)
        r99 = partialsort(v, 99, rev=true)
        vlines!(ax, [r99], color=:yellow, linestyle=:dashdot)
    end
    Label(fig[0, :], "APO medflux distributions — magenta: C2=10,000 nonlinearity cut; yellow: rank-99 threshold (fails to reach the >10k regime on 76/226)", fontsize=14)
    save(joinpath(PLOTD, "fig_apo_hists.png"), fig, px_per_unit=2)
end

## APO fig 3: per-fiber >10k counts vs the 99 rank budget
begin
    n10k = [count(mf_apo[fb, base_apo[fb, :]] .> 10_000) for fb in 1:300]
    fig = Figure(size=(1000, 420))
    ax = Axis(fig[1, 1], xlabel="fiberindx", ylabel="base-kept entries with medflux > 10,000",
        yscale=log10, title="Why rank99 cannot replace the APO constant: fibers 76 & 226 exceed the 99-entry budget")
    stem!(ax, 1:300, max.(n10k, 0.5), color=:deepskyblue, stemcolor=:deepskyblue, markersize=3)
    hlines!(ax, [99], color=:yellow, linestyle=:dashdot, label="rank-99 budget")
    for fb in (76, 226)
        text!(ax, fb, n10k[fb], text="f$fb: $(n10k[fb])", color=:orangered, fontsize=12, offset=(5, 3))
        scatter!(ax, [fb], [n10k[fb]], color=:orangered, markersize=8, marker=:diamond)
    end
    axislegend(ax, position=:rt, framevisible=false)
    save(joinpath(PLOTD, "fig_apo_rank_regime.png"), fig, px_per_unit=2)
end

println("figures written to ", PLOTD)
