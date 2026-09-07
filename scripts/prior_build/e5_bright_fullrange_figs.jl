## E5: FULL-WAVELENGTH-RANGE audit figures of the bright-line mask (AKS request,
## 2026-09-06): DR17's bright mask vs the calibrated `linedetect` detector, per fiber,
## across the entire 15075-17000 A grid, for several fibers at EACH telescope.
#
# Every masked and every missed line must be auditable, so resolution is not traded for
# range: each fiber gets its own figure of 5 stacked wavelength bands (~1740 px / ~385 A
# each), which at px_per_unit=2 is ~2 output pixels per spectral pixel.
#
# Mask agreement is drawn as three EXPLICIT, mutually exclusive categories rather than as
# overlapping translucent colours (the earlier figure relied on a cyan+magenta blend, so
# "both" and "DR17-only" were told apart only by hue):
#   grey    = BOTH agree (DR17 bright AND new detector)
#   cyan    = DR17-only  (we miss it)
#   magenta = new-only   (we add it)
# DR17-only spans are annotated with their central wavelength so the misses are
# individually accountable.
#
# Usage: julia --project=<plots env> e5_bright_fullrange_figs.jl [adjfib ...]
using HDF5, Statistics, Printf, Dates
using CairoMakie, ColorSchemes
set_theme!(merge(theme_black(), theme_latexfonts()))
CairoMakie.disable_mime!("svg", "pdf", "text/html")

const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
include(ARM * "scripts/prior_build/e5_bright_line_detect.jl")
const SC = get(ENV, "E5_SCREENS", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1/screens")
const DR = get(ENV, "E5_OLD_PRIORS", "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/sky_priors")
const PD = get(ENV, "E5_PLOTDIR", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/plots/e5_sky")
mkpath(PD)

wavetarg = 10 .^ range((4.179 - 125 * 6.0e-6), step=6.0e-6, length=8575 + 125)

# detector settings = the calibrated `linedetect` policy
const SCALE_WINDOW, KCUT, DILATION = 2001, 90.0, 4

function fiber_data(adjfib)
    n = lpad(adjfib, 3, "0")
    ms = h5read(joinpath(SC, "median_sky_$n.h5"), "median_sky")
    sub = Bool.(h5read(joinpath(SC, "median_sky_$n.h5"), "submsk"))
    x = fill(NaN, length(wavetarg)); x[sub] .= ms
    db = Bool.(h5read(joinpath(DR, "APOGEE_skyline_bright_svd_120_f$n.h5"), "submsk"))
    df = Bool.(h5read(joinpath(DR, "APOGEE_skyline_faint_svd_120_f$n.h5"), "submsk"))
    valid = sub .& (db .| df)
    scale = running_spread_fast(x, SCALE_WINDOW; kind=:mad, stride=50)
    base = [isfinite(x[i]) && isfinite(scale[i]) && scale[i] > 0 && x[i] > KCUT * scale[i]
            for i in eachindex(x)]
    new = dilate_msk(base, DILATION)
    for i in eachindex(new); isfinite(x[i]) || (new[i] = false); end
    st = line_overlap_stats(db, new, valid)
    return (n=n, adjfib=adjfib, x=x, sub=sub, db=db, valid=valid, new=new, scale=scale, st=st)
end

# short note on why each fiber is in the set (printed in the figure subtitle)
const WHY = Dict(
    1 => "detector edge (fiberindx 1)", 10 => "best APO agreement",
    76 => "worst APO agreement AND FPI-guide fiber (pooled-only prior)",
    200 => "lowest APO line recall (red-end DR17 misses)",
    295 => "detector edge (fiberindx 295)",
    388 => "FPI-guide fiber (pooled-only prior)",
    450 => "best LCO agreement", 519 => "WORST fiber overall; FPI-guide, ~zero own samples",
    570 => "worst-tier LCO agreement", 600 => "detector edge (fiberindx 300)")

function make_fig(d; nband=5)
    tele = d.adjfib <= 300 ? "APO" : "LCO"
    both = d.db .& d.new
    d17only = d.db .& .!d.new
    newonly = d.new .& .!d.db
    npix = length(wavetarg)
    edges = round.(Int, range(1, npix + 1, length=nband + 1))

    fig = Figure(size=(2000, 1700))
    why = get(WHY, d.adjfib, "")
    Label(fig[0, 1], @sprintf("%s fiber %s — %s\npixel IoU %.3f | DR17 lines recovered %.1f%% (%d missed) | precision %.3f (%d spurious lines) | bright %.2f%% (DR17 %.2f%%)",
            tele, d.n, why, d.st.pixel_iou, 100d.st.recall, d.st.n_ref_lines - d.st.recovered,
            d.st.precision, d.st.spurious,
            100count(d.new .& d.sub) / count(d.sub), 100d.st.ref_pix / count(d.valid)),
        fontsize=19, color=:white, tellwidth=false)

    for b in 1:nband
        lo, hi = edges[b], edges[b+1] - 1
        ax = Axis(fig[b, 1], xlabel=b == nband ? "wavelength [Å]" : "",
            ylabel="median sky line flux", yscale=Makie.pseudolog10,
            xgridvisible=false, ygridcolor=(:white, 0.12))
        xlims!(ax, wavetarg[lo], wavetarg[hi])

        # mask spans, drawn beneath the spectrum
        for (msk, col) in ((both, (:gray70, 0.30)), (d17only, (:cyan, 0.42)), (newonly, (:magenta, 0.42)))
            for r in mask_runs(msk)
                (last(r) < lo || first(r) > hi) && continue
                a = max(first(r), lo); z = min(last(r), hi)
                vspan!(ax, wavetarg[a], wavetarg[min(z + 1, npix)], color=col)
            end
        end

        xs = wavetarg[lo:hi]
        lines!(ax, xs, replace(d.x[lo:hi], NaN => 0.0), color=(:white, 0.9), linewidth=0.6)
        lines!(ax, xs, KCUT .* replace(d.scale[lo:hi], NaN => 0.0), color=:orange, linewidth=1.3)

        # annotate the DR17-only misses individually
        ymax = maximum(filter(isfinite, d.x[lo:hi]); init=1.0)
        for r in mask_runs(d17only)
            (last(r) < lo || first(r) > hi) && continue
            length(r) < 2 && continue
            λc = wavetarg[clamp((first(r) + last(r)) ÷ 2, 1, npix)]
            text!(ax, λc, ymax, text=@sprintf("%.0f", λc), color=:cyan, fontsize=10,
                rotation=pi / 2, align=(:left, :center))
        end
    end

    # one explicit legend for the whole figure
    els = [PolyElement(color=(:gray70, 0.30)), PolyElement(color=(:cyan, 0.42)),
        PolyElement(color=(:magenta, 0.42)), LineElement(color=:orange),
        LineElement(color=(:white, 0.9))]
    Legend(fig[nband+1, 1],
        els, ["both agree (pixels)", "DR17-only pixels (missed)", "new-only pixels (added, incl. line-wing extensions)",
            "detector cut = $(Int(KCUT)) × running MAD", "median sky line flux"],
        orientation=:horizontal, framevisible=false, labelsize=14, tellheight=true)

    out = joinpath(PD, "fig_bright_fullrange_$(tele)_f$(d.n).png")
    save(out, fig, px_per_unit=2)
    return out
end

fibers = isempty(ARGS) ? [1, 10, 76, 200, 388, 450, 519, 600] : parse.(Int, ARGS)
for a in fibers
    d = fiber_data(a)
    out = make_fig(d)
    @printf("wrote %s  (IoU %.3f, recall %.3f, prec %.3f)\n", out, d.st.pixel_iou,
        d.st.recall, d.st.precision)
    flush(stdout)
end
println("done")
