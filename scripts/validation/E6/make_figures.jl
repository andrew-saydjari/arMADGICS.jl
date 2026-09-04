# E6 step 6: figures for the prior-swap regression (dark theme).
# Run with the tfun_viz env: julia --project=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/tfun_viz make_figures.jl
using CairoMakie, ColorSchemes, Printf, Serialization, HDF5, Statistics, LinearAlgebra

black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")

const RUND = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E6_run"
const PLOTD = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/plots/e6_regression"
mkpath(PLOTD)

res = deserialize(joinpath(RUND, "structural_results.jdat"))
FIBERS = sort(collect(keys(res)))
wavetarg = 10 .^ range(4.179 - 125 * 6.0e-6, step=6.0e-6, length=8700)

fibcolors = ColorSchemes.glasbey_bw_minc_20_maxl_70_n256.colors

## Fig 1: singular value spectra old vs new
begin
    fig = Figure(size=(1000, 420))
    for (p, tele, fset) in ((1, "APO", filter(<=(300), FIBERS)), (2, "LCO", filter(>(300), FIBERS)))
        ax = Axis(fig[1, p], yscale=log10, xlabel="component k", ylabel=L"\lambda_k",
            title="$tele: singular-value spectra")
        for (ci, fib) in enumerate(fset)
            lines!(ax, 1:60, res[fib]["lam_old"], color=(fibcolors[ci], 0.9), linestyle=:dash)
            lines!(ax, 1:60, res[fib]["lam_new"], color=fibcolors[ci], label="f$fib")
        end
        axislegend(ax, position=:rt, framevisible=false)
    end
    Label(fig[0, :], "E6: built starCont prior spectra — dashed = OLD (2026_04_26 samples), solid = NEW (E4 2026_09_03)", fontsize=16)
    save(joinpath(PLOTD, "fig1_lambda_spectra.png"), fig, px_per_unit=2)
end

## Fig 2: lambda ratio new/old
begin
    fig = Figure(size=(900, 400))
    ax = Axis(fig[1, 1], xlabel="component k", ylabel=L"\lambda_k^{new}/\lambda_k^{old}",
        title="Per-mode variance ratio (NEW/OLD)")
    for (ci, fib) in enumerate(FIBERS)
        lines!(ax, 1:60, res[fib]["lam_new"] ./ res[fib]["lam_old"], color=fibcolors[ci],
            label="f$fib" * (fib > 300 ? " (lco)" : " (apo)"))
    end
    hlines!(ax, [1.0], color=:white, linestyle=:dot)
    axislegend(ax, position=:rt, framevisible=false, nbanks=2)
    save(joinpath(PLOTD, "fig2_lambda_ratio.png"), fig, px_per_unit=2)
end

## Fig 3: leading components old vs new for f295 (apo) and f595 (lco)
begin
    fig = Figure(size=(1200, 800))
    for (row, fib) in enumerate((295, 595))
        Vo = res[fib]["V_old_lead"]; Vn = res[fib]["V_new_lead"]
        for k in 1:3
            ax = Axis(fig[row, k], xlabel=k == 2 ? L"\lambda\;(\AA)" : "", ylabel=k == 1 ? "f$fib  V[:,·]" : "",
                title=row == 1 ? "component $k" : "")
            s = sign(dot(Vo[:, k], Vn[:, k]))
            lines!(ax, wavetarg, Vo[:, k], color=(:cyan, 0.7), label="old")
            lines!(ax, wavetarg, s .* Vn[:, k], color=(:orange, 0.7), label="new (sign-aligned)")
            lines!(ax, wavetarg, s .* Vn[:, k] .- Vo[:, k] .- 0.05, color=(:magenta, 0.8), label="Δ  (offset −0.05)")
            k == 1 && row == 1 && axislegend(ax, position=:rb, framevisible=false)
        end
    end
    Label(fig[0, :], "E6: leading starCont components, OLD vs NEW builds (apo f295 / lco f595)", fontsize=16)
    save(joinpath(PLOTD, "fig3_leading_components.png"), fig, px_per_unit=2)
end

## Fig 4: per-fiber alignment: |cor| of matched components and subspace cosines
begin
    fig = Figure(size=(1000, 420))
    ax1 = Axis(fig[1, 1], xlabel="component k", ylabel=L"|\mathrm{cor}(V^{old}_k, V^{new}_k)|",
        title="Matched-component alignment")
    ax2 = Axis(fig[1, 2], xlabel="principal angle index", ylabel=L"\cos\theta_i",
        title="Leading-10 subspace principal angles")
    for (ci, fib) in enumerate(FIBERS)
        lab = "f$fib" * (fib > 300 ? " (lco)" : " (apo)")
        lines!(ax1, 1:10, res[fib]["cors"], color=fibcolors[ci], label=lab)
        scatter!(ax1, 1:10, res[fib]["cors"], color=fibcolors[ci], markersize=5)
        lines!(ax2, 1:10, sort(res[fib]["cos_oldnew"], rev=true), color=fibcolors[ci], label=lab)
    end
    ylims!(ax1, 0, 1.05); ylims!(ax2, 0, 1.05)
    axislegend(ax1, position=:lb, framevisible=false)
    save(joinpath(PLOTD, "fig4_alignment.png"), fig, px_per_unit=2)
end

## Fig 5: leak drop-variant comparison (f450 and f595): lambda ratio + component deltas
begin
    fig = Figure(size=(1100, 800))
    for (row, fib, ncols) in ((1, 450, 2), (2, 595, 3))
        leakfile = joinpath(RUND, "built_new", "APOGEE_starcont_svd_60_f$(fib)_dropleak.h5")
        mainfile = joinpath(RUND, "built_new", "APOGEE_starcont_svd_60_f$(fib).h5")
        (isfile(leakfile) && isfile(mainfile)) || continue
        Vm, lm = h5open(mainfile, "r") do f; read(f["Vmat"]), read(f["λv"]); end
        Vd, ld = h5open(leakfile, "r") do f; read(f["Vmat"]), read(f["λv"]); end
        ax1 = Axis(fig[row, 1], xlabel="component k", ylabel=L"\lambda_k^{drop}/\lambda_k^{full} - 1",
            title="f$fib: drop $ncols leaked columns (of 10,000)")
        lines!(ax1, 1:60, ld ./ lm .- 1, color=:orange)
        scatter!(ax1, 1:60, ld ./ lm .- 1, color=:orange, markersize=4)
        ax2 = Axis(fig[row, 2], xlabel=L"\lambda\;(\AA)", ylabel="Δ V[:,k]",
            title="f$fib: component change (sign-aligned)")
        for k in 1:3
            s = sign(dot(Vm[:, k], Vd[:, k]))
            lines!(ax2, wavetarg, s .* Vd[:, k] .- Vm[:, k], color=fibcolors[k], label="k=$k")
        end
        axislegend(ax2, position=:rb, framevisible=false)
        @printf("leak variant f%d: max |λ ratio - 1| = %.3e ; max |cor-1| k1-3 = %.3e\n",
            fib, maximum(abs.(ld ./ lm .- 1)),
            maximum(1 - abs(dot(Vm[:, k], Vd[:, k]) / (norm(Vm[:, k]) * norm(Vd[:, k]))) for k in 1:3))
    end
    Label(fig[0, :], "E6 leak imprint: NEW lco builds vs rebuilds dropping the 59160/0018 + 60291/0009 sample columns", fontsize=16)
    save(joinpath(PLOTD, "fig5_leak_variant.png"), fig, px_per_unit=2)
end

## Fig 6: runtime swap summary — generation shift AND isolated E4 comparison
fold = joinpath(RUND, "fixture_oldprior.h5")       # production rough (2025 generation)
fctl = joinpath(RUND, "fixture_oldbuildprior.h5")  # OLD-samples f295 build (current builder+mask)
fnew = joinpath(RUND, "fixture_newprior.h5")       # NEW-samples f295 build
if isfile(fold) && isfile(fnew) && isfile(fctl)
    rd(f) = h5open(f, "r") do fh; Dict(k => read(fh[k]) for k in keys(fh)); end
    o = rd(fold); c = rd(fctl); n = rd(fnew)
    fibs = o["fiberindx"]
    fig = Figure(size=(1100, 440))
    ax1 = Axis(fig[1, 1], xlabel="fixture fiber", ylabel=L"\chi^2_{res}\;\mathrm{ratio} - 1",
        title="Residual chi2 change", xticks=(1:length(fibs), string.(fibs)))
    dgen = c["chi2res"] ./ o["chi2res"] .- 1  # generation shift (rough -> current build, OLD samples)
    de4 = n["chi2res"] ./ c["chi2res"] .- 1   # isolated E4 sample swap
    x = collect(1:length(fibs))
    barplot!(ax1, x .- 0.2, replace(dgen, NaN => 0.0), width=0.35, color=:gray70,
        label="generation shift (rough → 2026 build, OLD samples)")
    barplot!(ax1, x .+ 0.2, replace(de4, NaN => 0.0), width=0.35, color=:cyan,
        label="E4 sample swap (OLD → NEW build)")
    axislegend(ax1, position=:rt, framevisible=false)
    ax2 = Axis(fig[1, 2], xlabel="fixture fiber", ylabel="ΔRV (pix)",
        title="RV shift (isolated E4 swap)", xticks=(1:length(fibs), string.(fibs)))
    drv = n["RV_pixoff_final"] .- c["RV_pixoff_final"]
    barplot!(ax2, x, replace(drv, NaN => 0.0), color=:orange)
    Label(fig[0, :], "E6 runtime swap on M123 fixture (f295 builds; rough = production 2025 prior)", fontsize=16)
    save(joinpath(PLOTD, "fig6_runtime_swap.png"), fig, px_per_unit=2)
end

## Fig 7: fixture starContinuum component overlay for one healthy fiber
bo = joinpath(RUND, "fixture_oldprior_batchstyle.h5")
bc = joinpath(RUND, "fixture_oldbuildprior_batchstyle.h5")
bn = joinpath(RUND, "fixture_newprior_batchstyle.h5")
if isfile(bo) && isfile(bn) && isfile(bc)
    xo = h5open(bo, "r") do f; read(f["x_starContinuum_v0"]); end
    xc = h5open(bc, "r") do f; read(f["x_starContinuum_v0"]); end
    xn = h5open(bn, "r") do f; read(f["x_starContinuum_v0"]); end
    j = 1 # first fixture fiber (85, healthy sci)
    fig = Figure(size=(1100, 420))
    ax = Axis(fig[1, 1], xlabel=L"\lambda\;(\AA)", ylabel="flux",
        title="Fixture fiber 85: inferred starContinuum")
    lines!(ax, wavetarg, xo[:, j], color=(:gray70, 0.8), label="production rough")
    lines!(ax, wavetarg, xc[:, j], color=(:cyan, 0.8), label="OLD-samples build")
    lines!(ax, wavetarg, xn[:, j], color=(:orange, 0.8), label="NEW-samples build")
    axislegend(ax, position=:rb, framevisible=false)
    ax2 = Axis(fig[1, 2], xlabel=L"\lambda\;(\AA)", ylabel="Δ flux", title="differences")
    lines!(ax2, wavetarg, xc[:, j] .- xo[:, j], color=(:gray70, 0.9), label="generation (build − rough)")
    lines!(ax2, wavetarg, xn[:, j] .- xc[:, j], color=(:magenta, 0.9), label="E4 swap (new − old build)")
    axislegend(ax2, position=:rb, framevisible=false)
    save(joinpath(PLOTD, "fig7_fixture_continuum.png"), fig, px_per_unit=2)
end

println("figures written to ", PLOTD)
