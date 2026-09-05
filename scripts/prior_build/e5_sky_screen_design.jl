## E5 pass-1: outlier-screen DESIGN from the sampled data (AKS-agreed screens,
## E6 SVD-amplifies-rare-outliers lesson). Reads screens/e5_sample_stats.h5
## (written by e5_sky_run.jl --stage=stats), fits per-source-fiber robust
## distributions, chooses thresholds, and writes:
##   screens/e5_screen_decisions.h5   drop_NNN (UInt8 vectors over each source
##                                    fiber's kept samples) + threshold attrs
##   screens/e5_screen_report.txt     rejection statistics
##   <plot_dir>/*.png                 design figures (dark theme)
#
# Screen definitions (design rationale in the E5 MANIFEST):
#   BRIGHT (bit 1): robust z of log10(scale) within the source fiber's own
#     sample distribution, z = (x - med)/(1.4826 MAD), flagged when z > SCALE_ZMAX.
#     Sky brightness has heavy natural tails (moon, twilight, airglow storms);
#     the screen targets samples brighter than the fiber's distribution supports
#     (the E4 light-leak class sat far above the per-fiber locus).
#   CHI2 (bit 2): shape deviation from the fiber's median normalized spectrum,
#     chi2r = mean over good pixels of ivar*(flux - scale*ref)^2, both
#     full-spectrum and continuum-restricted (ref < its 60th pctl). Flagged when
#     chi2r_cont > CHI2CONT_MAX or chi2r_full > CHI2FULL_MAX. The continuum
#     restriction targets continuum-shape anomalies (leaks, scattered light)
#     without penalizing natural airglow-line variability.
# Thresholds are set GLOBALLY from the pooled standardized distributions (see
# code): defaults SCALE_ZMAX=6, CHI2CONT_MAX / CHI2FULL_MAX at the p99.9 break
# of their global distributions, overridable via env for the sensitivity check.
# Author - Andrew Saydjari (E5 pass 1)

using HDF5, Statistics, StatsBase, Printf, Dates
# figures need CairoMakie (separate plots project); E5_SKIP_FIGS=1 runs the
# threshold design + decisions only (rerun later with figs for the report)
const SKIP_FIGS = get(ENV, "E5_SKIP_FIGS", "0") == "1"
if !SKIP_FIGS
    using CairoMakie, ColorSchemes
    black_latexfonts = merge(theme_black(), theme_latexfonts())
    set_theme!(black_latexfonts)
    CairoMakie.disable_mime!("svg", "pdf", "text/html")
end

e5_out = get(ENV, "E5_OUT", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1")
screen_dir = joinpath(e5_out, "screens")
plot_dir = get(ENV, "E5_PLOTDIR", "/mnt/home/asaydjari/ceph/working/2026_09_04/plots/e5_sky")
mkpath(plot_dir)

SCALE_ZMAX = parse(Float64, get(ENV, "E5_SCALE_ZMAX", "6.0"))

stats_path = joinpath(screen_dir, "e5_sample_stats.h5")
fibs = Int[]; all_z = Float64[]; all_c2f = Float64[]; all_c2c = Float64[]
per_fib = Dict{Int,NamedTuple}()
h5open(stats_path, "r") do fh
    for k in sort(collect(keys(fh)))
        adjfib = parse(Int, k)
        g = fh[k]
        scale = read(g["scale"]); c2f = read(g["chi2r_full"]); c2c = read(g["chi2r_cont"])
        mjd = read(g["mjd"]); expnum = read(g["expnum"])
        n = length(scale)
        z = fill(NaN, n)
        pos = scale .> 0
        if count(pos) >= 10
            lx = log10.(scale[pos])
            med = median(lx); madv = max(1.4826 * median(abs.(lx .- med)), 1e-3)
            z[pos] .= (lx .- med) ./ madv
        end
        per_fib[adjfib] = (scale=scale, z=z, c2f=c2f, c2c=c2c, mjd=mjd, expnum=expnum)
        push!(fibs, adjfib)
        append!(all_z, filter(isfinite, z))
        append!(all_c2f, filter(isfinite, c2f))
        append!(all_c2c, filter(isfinite, c2c))
    end
end
ntot = sum(length(per_fib[f].scale) for f in fibs)
println("loaded stats: $(length(fibs)) fibers, $ntot samples")

# chi2 thresholds: p99.9 of the global distribution, with a floor at 3x p99
# (guards against a pathologically tight distribution making the screen twitchy)
function chi2_threshold(v; p_hi=0.999, p_ref=0.99, floor_fac=3.0)
    q_hi = quantile(v, p_hi)
    q_ref = quantile(v, p_ref)
    return max(q_hi, floor_fac * q_ref)
end
CHI2FULL_MAX = parse(Float64, get(ENV, "E5_CHI2FULL_MAX", string(chi2_threshold(all_c2f))))
CHI2CONT_MAX = parse(Float64, get(ENV, "E5_CHI2CONT_MAX", string(chi2_threshold(all_c2c))))
@printf("thresholds: SCALE_ZMAX=%.2f CHI2CONT_MAX=%.3f CHI2FULL_MAX=%.3f\n",
    SCALE_ZMAX, CHI2CONT_MAX, CHI2FULL_MAX)

# apply
ndrop_bright = 0; ndrop_chi2 = 0; ndrop = 0
drop_rows = Tuple{Int,Int,Int,Float64,Float64,Float64,Int}[] # adjfib, mjd, expnum, z, c2c, c2f, reason
h5open(joinpath(screen_dir, "e5_screen_decisions.h5"), "w") do fh
    attrs(fh)["SCALE_ZMAX"] = SCALE_ZMAX
    attrs(fh)["CHI2CONT_MAX"] = CHI2CONT_MAX
    attrs(fh)["CHI2FULL_MAX"] = CHI2FULL_MAX
    attrs(fh)["definition"] = "bit1: robust z of log10(scale) within source fiber > SCALE_ZMAX; bit2: chi2r_cont > CHI2CONT_MAX or chi2r_full > CHI2FULL_MAX (vs per-fiber median normalized spectrum)"
    for f in fibs
        global ndrop, ndrop_bright, ndrop_chi2
        st = per_fib[f]
        n = length(st.scale)
        drop = zeros(UInt8, n)
        for k in 1:n
            r = 0
            if isfinite(st.z[k]) && st.z[k] > SCALE_ZMAX
                r |= 1
            end
            if (isfinite(st.c2c[k]) && st.c2c[k] > CHI2CONT_MAX) ||
               (isfinite(st.c2f[k]) && st.c2f[k] > CHI2FULL_MAX)
                r |= 2
            end
            if r != 0
                drop[k] = UInt8(r)
                ndrop += 1
                (r & 1 != 0) && (ndrop_bright += 1)
                (r & 2 != 0) && (ndrop_chi2 += 1)
                push!(drop_rows, (f, st.mjd[k], st.expnum[k], st.z[k], st.c2c[k], st.c2f[k], r))
            end
        end
        fh["drop_" * lpad(f, 3, "0")] = drop
    end
end

open(joinpath(screen_dir, "e5_screen_report.txt"), "w") do io
    @printf(io, "E5 outlier screens (designed %s)\n", string(Dates.now()))
    @printf(io, "samples: %d total\n", ntot)
    @printf(io, "thresholds: SCALE_ZMAX=%.2f CHI2CONT_MAX=%.3f CHI2FULL_MAX=%.3f\n",
        SCALE_ZMAX, CHI2CONT_MAX, CHI2FULL_MAX)
    @printf(io, "dropped: %d (%.3f%%) | bright %d | chi2 %d (overlap %d)\n",
        ndrop, 100ndrop / ntot, ndrop_bright, ndrop_chi2, ndrop_bright + ndrop_chi2 - ndrop)
    println(io, "\nadjfib mjd expnum z_scale chi2r_cont chi2r_full reason")
    for r in sort(drop_rows)
        @printf(io, "%3d %5d %4d %8.2f %10.3f %10.3f %d\n", r...)
    end
end
println("dropped $ndrop of $ntot ($(round(100ndrop/ntot, digits=3))%) — report at $(joinpath(screen_dir, "e5_screen_report.txt"))")

## design figures
if !SKIP_FIGS
begin
    fig = Figure(size=(1500, 500))
    ax1 = Axis(fig[1, 1], xlabel="robust z of log10(scale)", ylabel="count",
        yscale=log10, title="bright-outlier screen")
    hist!(ax1, clamp.(all_z, -10, 30), bins=200)
    vlines!(ax1, [SCALE_ZMAX], color=:orange, linewidth=2)
    ax2 = Axis(fig[1, 2], xlabel="log10 chi2r (continuum-restricted)", ylabel="count",
        yscale=log10, title="chi2 screen (cont)")
    hist!(ax2, log10.(clamp.(all_c2c, 1e-3, 1e6)), bins=200)
    vlines!(ax2, [log10(CHI2CONT_MAX)], color=:orange, linewidth=2)
    ax3 = Axis(fig[1, 3], xlabel="log10 chi2r (full)", ylabel="count",
        yscale=log10, title="chi2 screen (full)")
    hist!(ax3, log10.(clamp.(all_c2f, 1e-3, 1e6)), bins=200)
    vlines!(ax3, [log10(CHI2FULL_MAX)], color=:orange, linewidth=2)
    save(joinpath(plot_dir, "fig_screens_design.png"), fig, px_per_unit=2)
end
begin
    # per-fiber drop fraction map
    dropfrac = zeros(600)
    for f in fibs
        st = per_fib[f]
        n = length(st.scale)
        nd = count(k -> (isfinite(st.z[k]) && st.z[k] > SCALE_ZMAX) ||
            (isfinite(st.c2c[k]) && st.c2c[k] > CHI2CONT_MAX) ||
            (isfinite(st.c2f[k]) && st.c2f[k] > CHI2FULL_MAX), 1:n)
        dropfrac[f] = n > 0 ? nd / n : NaN
    end
    fig = Figure(size=(1200, 420))
    ax = Axis(fig[1, 1], xlabel="adjfiberindx", ylabel="screen drop fraction",
        title="per-source-fiber screen rejection")
    scatter!(ax, 1:600, dropfrac, markersize=4, color=:cyan)
    vlines!(ax, collect(30.5:30:570.5), color=(:gray, 0.4), linewidth=0.5)
    vlines!(ax, [300.5], color=:orange, linewidth=1)
    save(joinpath(plot_dir, "fig_screens_perfiber.png"), fig, px_per_unit=2)
end
end # !SKIP_FIGS
println(SKIP_FIGS ? "figures SKIPPED (E5_SKIP_FIGS=1)" : "figures -> $plot_dir")
