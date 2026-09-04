# E4b step 6: overlay figures for the LCO regeneration (dark theme).
# Run with the tfun_viz env:
#   cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run && nice -n 10 \
#   julia +1.11.6 --project=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfun_viz \
#   /mnt/home/asaydjari/gitcode/worktrees/arM-E4b/scripts/prior_build/E4b/make_figures_E4b.jl
using CairoMakie, ColorSchemes, Printf, Serialization, HDF5, Statistics, LinearAlgebra

black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E4b_run"
const BUILT = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/starCont_pass1/built_lco"
const E6D = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E6_run"
const PLOTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/plots/e4b_lco_regen"
mkpath(PLOTD)

wavetarg = 10 .^ range(4.179 - 125 * 6.0e-6, step=6.0e-6, length=8700)
loadprior(path) = h5open(path, "r") do f
    read(f["Vmat"]), read(f["λv"])
end
newp(f) = joinpath(BUILT, "APOGEE_starcont_svd_60_f$(lpad(f,3,"0")).h5")
e6p(f, s="") = joinpath(E6D, "built_new", "APOGEE_starcont_svd_60_f$(lpad(f,3,"0"))$(s).h5")

V450n, l450n = loadprior(newp(450))
V450L, l450L = loadprior(e6p(450))
V450d, l450d = loadprior(e6p(450, "_dropleak"))
V595n, l595n = loadprior(newp(595))
V595L, l595L = loadprior(e6p(595))
V595d, l595d = loadprior(e6p(595, "_dropleak"))

cLEAK = :orangered; cDROP = :cyan; cREGEN = :springgreen

## Fig 1: f450 + f595 lambda spectra, three variants
begin
    fig = Figure(size=(1100, 440))
    for (p, fib, lL, ld, ln_) in ((1, 450, l450L, l450d, l450n), (2, 595, l595L, l595d, l595n))
        ax = Axis(fig[1, p], yscale=log10, xlabel="component k", ylabel=L"\lambda_k",
            title="f$fib (lco)")
        lines!(ax, 1:60, lL, color=cLEAK, linewidth=2, label="leaked build (20260903 samples)")
        lines!(ax, 1:60, ld, color=cDROP, linewidth=2, linestyle=:dash, label="E6 drop-variant (columns removed)")
        lines!(ax, 1:60, ln_, color=cREGEN, linewidth=2, label="E4b regen (C2_LCO=3000 samples)")
        axislegend(ax, position=:rt, framevisible=false)
    end
    text!(fig.content[1], 8, l450L[2], text=@sprintf("λ2 deflation %.2fx", l450L[2] / l450n[2]),
        color=:white, fontsize=14)
    Label(fig[0, :], "E4b: built LCO starCont prior spectra — leaked vs drop-variant vs regenerated", fontsize=16)
    save(joinpath(PLOTD, "fig1_lambda_spectra_f450_f595.png"), fig, px_per_unit=2)
end

## Fig 2: f450 mode-2 imprint gone (full + zoom on 15,250-15,350 A)
begin
    fig = Figure(size=(1200, 640))
    for (row, rng, ttl) in ((1, 1:8700, "full range"), (2, findall(w -> 15150 <= w <= 15450, wavetarg), "leak imprint region"))
        ax = Axis(fig[row, 1], xlabel=row == 2 ? L"\lambda\;(\AA)" : "", ylabel=L"V[:,2]", title="f450 mode 2 — $ttl")
        sL = 1.0
        sd = sign(dot(V450L[:, 2], V450d[:, 2])); sn = sign(dot(V450L[:, 2], V450n[:, 2]))
        lines!(ax, wavetarg[rng], V450L[rng, 2], color=(cLEAK, 0.9), linewidth=1.2, label="leaked build")
        lines!(ax, wavetarg[rng], sd .* V450d[rng, 2], color=(cDROP, 0.9), linewidth=1.2, label="drop-variant")
        lines!(ax, wavetarg[rng], sn .* V450n[rng, 2], color=(cREGEN, 0.9), linewidth=1.2, label="E4b regen")
        row == 1 && axislegend(ax, position=:rt, framevisible=false)
    end
    Label(fig[0, :], "E4b: f450 mode-2 — the ~15,250-15,350 Å leaked-domeflat imprint is gone", fontsize=16)
    save(joinpath(PLOTD, "fig2_f450_mode2_imprint.png"), fig, px_per_unit=2)
end

## Fig 3: lambda ratios vs leaked build (f450, f595): regen and drop-variant
begin
    fig = Figure(size=(1100, 420))
    for (p, fib, lL, ld, ln_) in ((1, 450, l450L, l450d, l450n), (2, 595, l595L, l595d, l595n))
        ax = Axis(fig[1, p], xlabel="component k", ylabel=L"\lambda_k / \lambda_k^{leaked}",
            title="f$fib", yscale=log10)
        lines!(ax, 1:60, ld ./ lL, color=cDROP, linewidth=2, label="drop-variant / leaked")
        lines!(ax, 1:60, ln_ ./ lL, color=cREGEN, linewidth=2, label="E4b regen / leaked")
        hlines!(ax, [1.0], color=:white, linestyle=:dot)
        axislegend(ax, position=:rb, framevisible=false)
    end
    Label(fig[0, :], "E4b: per-mode variance vs the leaked build — regen matches the drop-variant prediction", fontsize=16)
    save(joinpath(PLOTD, "fig3_lambda_ratio_vs_leaked.png"), fig, px_per_unit=2)
end

## Fig 4: corpus overview — per-fiber lambda1, lambda2 for all 300 new builds
begin
    q = deserialize(joinpath(RUND, "qa_built_results.jdat"))
    fig = Figure(size=(1100, 420))
    ax1 = Axis(fig[1, 1], xlabel="adjusted fiber index", ylabel=L"\lambda_1", yscale=log10, title="leading mode")
    ax2 = Axis(fig[1, 2], xlabel="adjusted fiber index", ylabel=L"\lambda_2", yscale=log10, title="second mode")
    scatter!(ax1, 301:600, q.lam1, color=:deepskyblue, markersize=4)
    scatter!(ax2, 301:600, q.lam2, color=:deepskyblue, markersize=4)
    for (ax, v) in ((ax1, q.lam1), (ax2, q.lam2))
        for (fib, c) in ((450, cLEAK), (388, :magenta), (519, :magenta), (448, :yellow), (459, :yellow))
            scatter!(ax, [fib], [v[fib-300]], color=c, markersize=10, marker=:diamond)
            text!(ax, fib, v[fib-300], text="f$fib", color=c, fontsize=10, offset=(4, 4))
        end
    end
    Label(fig[0, :], "E4b: built LCO prior scales across all 300 fibers (diamonds: f450 ex-leak; magenta f388/f519 bright-fiber trim; yellow f448/f459)", fontsize=14)
    save(joinpath(PLOTD, "fig4_corpus_lambda.png"), fig, px_per_unit=2)
end

## Fig 5: sample-level — predicted vs measured changed columns per fiber
begin
    pred = deserialize(joinpath(RUND, "rng_predict_final.jdat"))
    npred = [haskey(pred.diffcols, fb) ? length(pred.diffcols[fb]) : 0 for fb in 1:300]
    fig = Figure(size=(1100, 400))
    ax = Axis(fig[1, 1], xlabel="adjusted fiber index", ylabel="changed sample columns (of 10,000)",
        title="RNG-stream prediction of per-fiber sample changes (QA verified measured == predicted)")
    stem!(ax, 301:600, npred, color=:deepskyblue, stemcolor=:deepskyblue, markersize=3)
    for (fib, c) in ((450, cLEAK), (388, :magenta), (519, :magenta))
        scatter!(ax, [fib], [npred[fib-300]], color=c, markersize=10, marker=:diamond)
    end
    save(joinpath(PLOTD, "fig5_changed_columns.png"), fig, px_per_unit=2)
end

println("figures written to ", PLOTD)
