## E5 pass-1 sky-prior QA: sample-count accounting, built-prior spectra,
## drop-variant sensitivity, old-vs-new reference comparison, public figures.
# Run AFTER --stage=build (and the drop-variant builds for E5_QA_FIBERS).
# Env: E5_OUT, E5_PLOTDIR, E5_QA_FIBERS (default "10,76,295,351,519,460"),
#      E5_OLD_PRIORS (default 2025_07_31 dump; REFERENCE comparison only).
# Author - Andrew Saydjari (E5 pass 1)

using HDF5, Statistics, StatsBase, Printf, LinearAlgebra, Dates
using CairoMakie, ColorSchemes
black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")

e5_out = get(ENV, "E5_OUT", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1")
plot_dir = get(ENV, "E5_PLOTDIR", "/mnt/home/asaydjari/ceph/working/2026_09_04/plots/e5_sky")
old_prior_dir = get(ENV, "E5_OLD_PRIORS", "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/sky_priors")
qa_fibers = [parse(Int, t) for t in split(get(ENV, "E5_QA_FIBERS", "10,76,295,351,519,460"), ",")]
mkpath(plot_dir)
built = joinpath(e5_out, "built")
built_u = joinpath(e5_out, "built_unscreened")
report = String[]

wavetarg = 10 .^ range((4.179 - 125 * 6.0e-6), step=6.0e-6, length=8575 + 125)

skycont_name(d, f) = joinpath(d, "APOGEE_skycont_svd_30_f" * lpad(f, 3, "0") * ".h5")
skyfaint_name(d, f) = joinpath(d, "APOGEE_skyline_faint_svd_120_f" * lpad(f, 3, "0") * ".h5")
skygspice_name(d, f) = joinpath(d, "APOGEE_skyline_faint_GSPICE_svd_120_f" * lpad(f, 3, "0") * ".h5")

## 1. sample accounting
nsamp_assembled = zeros(Int, 600)
for f in 1:600
    p = joinpath(e5_out, "samples", "skyflux_" * lpad(f, 3, "0") * ".h5")
    if isfile(p)
        nsamp_assembled[f] = h5open(p, "r") do fh; size(fh["skyflux"], 2); end
    end
end
guard = h5open(joinpath(e5_out, "screens", "e5_guard_counts.h5"), "r") do fh
    (adjfib=read(fh["adjfiberindx"]), nc=read(fh["n_candidates"]), nk=read(fh["n_kept_guard"]))
end
pool = h5open(joinpath(e5_out, "outlists", "e5_pool_summary.h5"), "r") do fh
    (n=read(fh["n_pooled_candidates"]),)
end
push!(report, @sprintf("assembled samples/target: min=%d med=%.0f max=%d (floor ~1015-1022 for GSPICE)",
    minimum(nsamp_assembled), median(nsamp_assembled), maximum(nsamp_assembled)))

begin
    fig = Figure(size=(1300, 450))
    ax = Axis(fig[1, 1], xlabel="adjfiberindx", ylabel="samples",
        title="E5 pass-1 pooled sample accounting", yscale=log10)
    lines!(ax, 1:600, max.(pool.n, 1), color=:gray70, label="pooled candidates (list)")
    scatter!(ax, 1:600, max.(nsamp_assembled, 1), markersize=4, color=:cyan,
        label="assembled (guard+screens)")
    hlines!(ax, [1022], color=:orange, linestyle=:dash, label="GSPICE hard floor")
    vlines!(ax, collect(30.5:30:570.5), color=(:gray, 0.3), linewidth=0.5)
    vlines!(ax, [300.5], color=:orange, linewidth=1)
    axislegend(ax, position=:rb)
    save(joinpath(plot_dir, "fig_qa_counts.png"), fig, px_per_unit=2)
end

## 2. built-prior eigenvalue spectra (all fibers heatmap + examples)
function lam_matrix(namer, dir, nsub)
    L = fill(NaN, nsub, 600)
    for f in 1:600
        p = namer(dir, f)
        isfile(p) || continue
        L[:, f] .= h5read(p, "λv")[1:nsub]
    end
    return L
end
Lc = lam_matrix(skycont_name, built, 30)
Lf = lam_matrix(skyfaint_name, built, 120)
Lg = lam_matrix(skygspice_name, built, 120)
nbuilt = count(.!isnan.(Lc[1, :]))
push!(report, @sprintf("built priors present: skyCont %d/600, faint %d/600, faint-GSPICE %d/600",
    nbuilt, count(.!isnan.(Lf[1, :])), count(.!isnan.(Lg[1, :]))))

begin
    fig = Figure(size=(1400, 800))
    for (i, (L, ttl)) in enumerate([(Lc, "skyCont λv (rank 30)"), (Lf, "skyLines faint λv (rank 120)"), (Lg, "skyLines faint GSPICE λv (rank 120)")])
        nbuilt_i = count(.!isnan.(L[1, :]))
        ax = Axis(fig[i, 1], xlabel=i == 3 ? "adjfiberindx" : "", ylabel="mode",
            title=ttl * " — $(nbuilt_i)/600 built")
        Z = log10.(max.(L, 1e-12))
        fin = filter(isfinite, Z)
        if isempty(fin)
            text!(ax, 0.5, 0.5, text="no priors built yet", space=:relative, align=(:center, :center))
            continue
        end
        # explicit colorrange: Makie cannot derive limits when columns are all-NaN
        # (partial builds), which previously threw inside Colorbar/extract_colormap
        crange = extrema(fin)
        crange = crange[1] == crange[2] ? (crange[1] - 1, crange[2] + 1) : crange
        hm = heatmap!(ax, 1:600, 1:size(L, 1), Z', colormap=:CET_L8, colorrange=crange)
        Colorbar(fig[i, 2], hm, label="log10 λ")
    end
    save(joinpath(plot_dir, "fig_qa_lambda_heatmaps.png"), fig, px_per_unit=2)
end

## 3. drop-variant sensitivity (screened vs unscreened) for QA fibers
princ_angles(V1, V2, k) = svdvals(Matrix(qr(V1[:, 1:k]).Q)' * Matrix(qr(V2[:, 1:k]).Q))
for f in qa_fibers
    for (namer, k, tag) in [(skycont_name, 10, "skyCont k=10"), (skygspice_name, 30, "GSPICE k=30")]
        p1, p2 = namer(built, f), namer(built_u, f)
        if isfile(p1) && isfile(p2)
            V1 = h5read(p1, "Vmat"); V2 = h5read(p2, "Vmat")
            s = princ_angles(V1, V2, k)
            worst = acosd(clamp(minimum(s), 0, 1))
            push!(report, @sprintf("drop-variant fiber %03d %s: worst principal angle %.3f deg (cos=%.6f)", f, tag, worst, minimum(s)))
        else
            push!(report, @sprintf("drop-variant fiber %03d %s: MISSING (%s / %s)", f, tag, isfile(p1), isfile(p2)))
        end
    end
end

## 4. old-vs-new reference comparison (QA context only; old = DR17-era corpus)
for f in qa_fibers
    oldc = joinpath(old_prior_dir, "APOGEE_skycont_svd_30_f" * lpad(f, 3, "0") * ".h5")
    newc = skycont_name(built, f)
    if isfile(oldc) && isfile(newc)
        Vo = h5read(oldc, "Vmat"); Vn = h5read(newc, "Vmat")
        s = princ_angles(Vo, Vn, 5)
        push!(report, @sprintf("old-vs-new skyCont fiber %03d: worst principal angle (k=5) %.2f deg", f, acosd(clamp(minimum(s), 0, 1))))
    end
end

## 5. example prior modes + mean sky for a QA fiber pair
begin
    fig = Figure(size=(1500, 900))
    for (row, f) in enumerate(qa_fibers[1:min(3, end)])
        p = skygspice_name(built, f)
        isfile(p) || continue
        V = h5read(p, "Vmat")
        ax = Axis(fig[row, 1], xlabel=row == 3 ? "wavelength [Å]" : "", ylabel="mode amp",
            title="fiber $(lpad(f,3,"0")) skyLines faint GSPICE modes 1-3")
        for m in 1:3
            lines!(ax, wavetarg, V[:, m] .+ 0.35 * (m - 1) * maximum(abs.(V[:, 1])), linewidth=0.5)
        end
    end
    save(joinpath(plot_dir, "fig_qa_modes.png"), fig, px_per_unit=2)
end

open(joinpath(plot_dir, "qa_report.txt"), "w") do io
    println(io, "E5 pass-1 sky-prior QA — $(Dates.now())")
    for r in report
        println(io, r)
    end
end
println(join(report, "\n"))
println("QA figures -> $plot_dir")
