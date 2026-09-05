# e7_plots.jl — figures for the E7 per-fiber starLines prior build.
# Run: nice -n 12 julia +1.11.0 --project=. e7_plots.jl
using CairoMakie, ColorSchemes, Printf
using HDF5, Statistics, StatsBase, LinearAlgebra

black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")

const PRIOR_DIR = "/mnt/ceph/users/sdssv/work/asaydjari/"
const NEW_DIR = PRIOR_DIR * "2026_09_05/prior_outputs/starLines_perfiber/"
const OLD_DIR = PRIOR_DIR * "2026_09_05/E7_run/old_norm94_ref/"
const QA_H5 = PRIOR_DIR * "2026_09_05/E7_run/qa_starlines_results.h5"
const REF_H5 = PRIOR_DIR * "2025_07_31/prior_dump/APOGEE_stellar_kry_50_subpix_th_22500.h5"
const PD = PRIOR_DIR * "2026_09_05/plots/e7_starlines/"
mkpath(PD)

const CM_SEQ = :linear_bmy_10_95_c71_n256
const CM_DIV = :diverging_bwr_40_95_c42_n256
const COL_NEW = "#00E5FF" # bright cyan
const COL_OLD = "#FFB300" # bright amber
const COL_REF = "#B0BEC5" # grey
const COL_ACC = "#FF4081" # bright pink

fname(dir, fib) = joinpath(dir, "APOGEE_stellar_kry_50_subpix_f" * lpad(fib, 3, "0") * ".h5")
wavetarg = 10 .^ range(4.179 - 125 * 6.0e-6, step = 6.0e-6, length = 8700)

## Fig 1: component overlays old vs new vs refLSF at two fibers ---------------
for (fib, comps) in ((150, (1, 5)), (450, (1, 5)))
    isfile(fname(NEW_DIR, fib)) && isfile(fname(OLD_DIR, fib)) || continue
    Vn = h5read(fname(NEW_DIR, fib), "Vmat")[:, :, 6]
    Vo = h5read(fname(OLD_DIR, fib), "Vmat")[:, :, 6]
    Vr = h5read(REF_H5, "Vmat")[:, :, 6]
    fig = Figure(size = (1400, 700))
    px = 3500:4100 # a line-rich stretch
    for (row, k) in enumerate(comps)
        ax = Axis(fig[row, 1], xlabel = row == length(comps) ? "wavelength (Å)" : "",
            ylabel = "component $k",
            title = row == 1 ?
                    "starLines component, fiber $fib (subpix 6): old 2023-LSF vs new FPI-LSF vs refLSF" :
                    "")
        lines!(ax, wavetarg[px], Vr[px, k], color = COL_REF, linewidth = 1,
            label = "refLSF (R=22500)")
        lines!(ax, wavetarg[px], Vo[px, k], color = COL_OLD, linewidth = 1.2,
            label = "old (2023 LSF)")
        lines!(ax, wavetarg[px], Vn[px, k], color = COL_NEW, linewidth = 1.2,
            label = "new (FPI LSF)")
        row == 1 && axislegend(ax, position = :rt, framevisible = false)
        ax2 = Axis(fig[row, 2], xlabel = row == length(comps) ? "wavelength (Å)" : "",
            ylabel = "new − old")
        lines!(ax2, wavetarg[px], Vn[px, k] .- Vo[px, k], color = COL_ACC, linewidth = 1)
    end
    save(joinpath(PD, "components_f$(lpad(fib,3,"0")).png"), fig, px_per_unit = 2)
    println("fig components f$fib done")
end

## Fig 2: QA summary — cosine similarity + principal angles -------------------
if isfile(QA_H5)
    old_fibs = h5read(QA_H5, "old_fibers")
    cs = h5read(QA_H5, "cosine_similarity_findx6")
    pa = h5read(QA_H5, "principal_angles_deg")
    fig = Figure(size = (1400, 550))
    ax = Axis(fig[1, 1], xlabel = "component", ylabel = "cosine similarity (old·new)",
        title = "per-component old-vs-new similarity (subpix 6)")
    cols = get(colorschemes[CM_SEQ], range(0.15, 0.95, length = length(old_fibs)))
    for (j, fib) in enumerate(old_fibs)
        lines!(ax, 1:50, cs[:, j], color = cols[j], linewidth = 1.5,
            label = "f$(lpad(fib,3,"0"))")
    end
    axislegend(ax, position = :lb, framevisible = false, nbanks = 2)
    ax2 = Axis(fig[1, 2], xlabel = "principal-angle index", ylabel = "angle (deg)",
        title = "old/new 50-dim subspace principal angles")
    for (j, fib) in enumerate(old_fibs)
        lines!(ax2, 1:50, sort(pa[:, j]), color = cols[j], linewidth = 1.5)
    end
    save(joinpath(PD, "qa_similarity.png"), fig, px_per_unit = 2)
    println("fig qa_similarity done")

    ## Fig 3: column-norm ratio map (all built fibers)
    built = h5read(QA_H5, "built_fibers")
    cn = h5read(QA_H5, "colnorm_ratio_findx6")
    fig = Figure(size = (1200, 500))
    ax = Axis(fig[1, 1], xlabel = "adjfiberindx", ylabel = "component",
        title = "column-norm ratio ‖new‖/‖refLSF‖ (subpix 6)")
    hm = heatmap!(ax, built, 1:50, permutedims(cn), colormap = CM_SEQ)
    Colorbar(fig[1, 2], hm)
    save(joinpath(PD, "qa_colnorm_ratio.png"), fig, px_per_unit = 2)
    println("fig qa_colnorm_ratio done")
end

## Fig 4: fixture deltas (hack vs per-fiber) ----------------------------------
fx_a = PRIOR_DIR * "2026_09_05/E7_run/fixture_refLSFhack.h5"
fx_b = PRIOR_DIR * "2026_09_05/E7_run/fixture_perfiber.h5"
if isfile(fx_a) && isfile(fx_b)
    fibs = h5read(fx_a, "fiberindx")
    rva, rvb = h5read(fx_a, "RV_pixoff_final"), h5read(fx_b, "RV_pixoff_final")
    ca, cb = h5read(fx_a, "RV_minchi2_final"), h5read(fx_b, "RV_minchi2_final")
    ok = isfinite.(rva) .& isfinite.(rvb)
    fig = Figure(size = (1250, 500))
    ax = Axis(fig[1, 1], xlabel = "fixture fiber", ylabel = "ΔRV (pix, per-fiber − hack)",
        title = "M123 fixture: RV shift from replacing the refLSF hack",
        xticks = (1:count(ok), string.(fibs[ok])))
    barplot!(ax, 1:count(ok), (rvb .- rva)[ok], color = COL_NEW)
    ax2 = Axis(fig[1, 2], xlabel = "fixture fiber", ylabel = "Δ(RV min χ²) (per-fiber − hack)",
        title = "χ² improvement (negative = better fit)",
        xticks = (1:count(ok), string.(fibs[ok])))
    barplot!(ax2, 1:count(ok), (cb .- ca)[ok], color = COL_ACC)
    save(joinpath(PD, "fixture_deltas.png"), fig, px_per_unit = 2)
    println("fig fixture_deltas done")
    for i in findall(ok)
        @printf("fib %d: RV %.4f -> %.4f (Δ %.4f pix); minchi2 %.2f -> %.2f (Δ %.2f)\n",
            fibs[i], rva[i], rvb[i], rvb[i] - rva[i], ca[i], cb[i], cb[i] - ca[i])
    end
end

println("plots -> $PD")
