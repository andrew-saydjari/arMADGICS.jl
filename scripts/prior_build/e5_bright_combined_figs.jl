## E5: VALIDATION FIGURES for the COMBINED bright sky-line mask (AKS 2026-09-07).
#
# AKS: "To validate, please show me the fibers where the result differs the most from the
# final delivered bright mask, with the sign of the difference clear (over or undermasked
# on the per fiber products you made compared to the one made by combining all of the
# fibers)."
#
# Same full-range stacked-band layout, log-y (pseudolog10), dark theme and MUTUALLY
# EXCLUSIVE colour categories as the fig_bright_fullrange_* series, with the reference
# swapped from DR17 to the delivered COMBINED mask. Colour convention is carried over
# unchanged (cyan = reference-only, magenta = candidate-only), which here means:
#
#   grey    = BOTH  (per-fiber AND combined)
#   cyan    = combined-only  -> the per-fiber product UNDER-masked this pixel
#   magenta = per-fiber-only -> the per-fiber product OVER-masked this pixel
#
# Figures produced:
#   fig_bright_combined_summary.png            evidence for combining at all
#   fig_bright_combined_diff_<TELE>_f<NNN>.png worst offenders, BOTH directions
# Usage: julia --project=<plots env> e5_bright_combined_figs.jl [variant]
# Author - Andrew Saydjari (E5 pass 1)
using HDF5, Statistics, Printf, Dates
using CairoMakie, ColorSchemes
set_theme!(merge(theme_black(), theme_latexfonts()))
CairoMakie.disable_mime!("svg", "pdf", "text/html")

const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
include(ARM * "scripts/prior_build/e5_bright_line_detect.jl")
const SC = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1/screens"
const DR = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/sky_priors"
const RESDIR = get(ENV, "E5_RESDIR", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_07/e5_combined")
const PD = get(ENV, "E5_PLOTDIR", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/plots/e5_sky")
mkpath(PD)
const VARIANT = isempty(ARGS) ? "drop12" : ARGS[1]
const SCALE_WINDOW, KCUT, DILATION = 2001, 90.0, 4

wavetarg = 10 .^ range((4.179 - 125 * 6.0e-6), step=6.0e-6, length=8575 + 125)
NPIX = length(wavetarg)

M = Bool.(h5read(joinpath(RESDIR, "e5_bright_perfiber.h5"), "mask"))
S = Bool.(h5read(joinpath(RESDIR, "e5_bright_perfiber.h5"), "submsk"))
nflag = h5read(joinpath(RESDIR, "e5_bright_perfiber.h5"), "nflag")
nvalid = h5read(joinpath(RESDIR, "e5_bright_perfiber.h5"), "nvalid")
CB = Bool.(h5read(joinpath(RESDIR, "e5_bright_combined.h5"), "mask_" * VARIANT))
nflag_apo = h5read(joinpath(RESDIR, "e5_bright_combined.h5"), "nflag_apo")
nflag_lco = h5read(joinpath(RESDIR, "e5_bright_combined.h5"), "nflag_lco")
dr_nflag = h5read(joinpath(RESDIR, "e5_bright_combined.h5"), "dr17_nflag")

## ------------------------------------------------------------------ summary figure
begin
    fig = Figure(size=(1900, 1150))
    Label(fig[0, 1:2], "Combining the bright sky-line mask across all 600 fibers — evidence",
        fontsize=24, color=:white, tellwidth=false)

    # (a) detection count vs wavelength
    ax = Axis(fig[1, 1:2], xlabel="wavelength [Å]", ylabel="fibers flagging this pixel",
        title="per-pixel detection count (of 600 fibers); flat-topped plateaus = lines every fiber agrees on",
        xgridvisible=false, ygridcolor=(:white, 0.12))
    band!(ax, wavetarg, zeros(NPIX), Float64.(nflag), color=(:cyan, 0.35))
    lines!(ax, wavetarg, Float64.(nflag), color=:cyan, linewidth=0.5)
    lines!(ax, wavetarg, Float64.(nvalid), color=(:orange, 0.8), linewidth=1.0)
    hlines!(ax, [300.0], color=:magenta, linewidth=1.5, linestyle=:dash)
    hlines!(ax, [3.0], color=(:white, 0.55), linewidth=1.0, linestyle=:dot)
    axislegend(ax, [LineElement(color=:cyan), LineElement(color=(:orange, 0.8)),
            LineElement(color=:magenta, linestyle=:dash), LineElement(color=(:white, 0.55), linestyle=:dot)],
        ["fibers flagging", "fibers where pixel is valid",
            "MAJORITY cut (nflag=300, recommended)", "drop-1-2 cut (nflag=3)"],
        position=:lt, labelsize=13, framevisible=false)
    xlims!(ax, wavetarg[1], wavetarg[end])

    # (b) histogram of nflag among flagged pixels
    nf = filter(>(0), nflag)
    ax2 = Axis(fig[2, 1], xlabel="fibers flagging the pixel", ylabel="pixels",
        yscale=log10, title="detection count of FLAGGED pixels — strongly bimodal")
    hist!(ax2, Float64.(nf), bins=60, color=(:cyan, 0.7))
    vlines!(ax2, [2.5], color=:magenta, linewidth=2, linestyle=:dash)

    # (c) APO vs LCO detection count
    ax3 = Axis(fig[2, 2], xlabel="fibers flagging (APO, of 300)", ylabel="fibers flagging (LCO, of 300)",
        title="per-telescope agreement: off-diagonal = a line one site sees and the other does not")
    sel = (nflag_apo .+ nflag_lco) .> 0
    scatter!(ax3, Float64.(nflag_apo[sel]), Float64.(nflag_lco[sel]),
        markersize=3, color=(:cyan, 0.35))
    lines!(ax3, [0, 300], [0, 300], color=(:orange, 0.8), linewidth=1)

    # (d) per-fiber over/under vs the combined mask
    over = [count(M[:, f] .& .!CB .& S[:, f]) for f in 1:600]
    under = [count(CB .& .!M[:, f] .& S[:, f]) for f in 1:600]
    ax4 = Axis(fig[3, 1:2], xlabel="adjfiberindx", ylabel="pixels differing from combined mask",
        title="per-fiber difference from the delivered combined mask ($VARIANT), signed")
    scatter!(ax4, 1:600, Float64.(over), markersize=4, color=:magenta,
        label="per-fiber OVER-masked (per-fiber only)")
    scatter!(ax4, 1:600, -Float64.(under), markersize=4, color=:cyan,
        label="per-fiber UNDER-masked (combined only)")
    hlines!(ax4, [0.0], color=(:white, 0.4), linewidth=1)
    vlines!(ax4, [300.5], color=:orange, linewidth=1.5)
    text!(ax4, 150, maximum(over) * 0.9, text="APO", color=:orange, fontsize=15, align=(:center, :center))
    text!(ax4, 450, maximum(over) * 0.9, text="LCO", color=:orange, fontsize=15, align=(:center, :center))
    axislegend(ax4, position=:rt, labelsize=13, framevisible=false)

    save(joinpath(PD, "fig_bright_combined_summary.png"), fig, px_per_unit=2)
    println("wrote fig_bright_combined_summary.png"); flush(stdout)
end

## ------------------------------------------------------- per-fiber difference figures
function fiber_spec(adjfib)
    n = lpad(adjfib, 3, "0")
    ms = h5read(joinpath(SC, "median_sky_$n.h5"), "median_sky")
    sub = Bool.(h5read(joinpath(SC, "median_sky_$n.h5"), "submsk"))
    x = fill(NaN, NPIX); x[sub] .= ms
    scale = running_spread_fast(x, SCALE_WINDOW; kind=:mad, stride=50)
    return x, scale
end

function make_diff_fig(adjfib, why; nband=5)
    tele = adjfib <= 300 ? "APO" : "LCO"
    n = lpad(adjfib, 3, "0")
    x, scale = fiber_spec(adjfib)
    sub = S[:, adjfib]; pf = M[:, adjfib]
    both = pf .& CB .& sub
    overm = pf .& .!CB .& sub          # per-fiber OVER-masked  (magenta = candidate-only)
    underm = CB .& .!pf .& sub         # per-fiber UNDER-masked (cyan   = reference-only)
    nov, nun = count(overm), count(underm)
    npf, ncb = count(pf), count(CB .& sub)
    iou = count(both) / max(count((pf .| CB) .& sub), 1)

    edges = round.(Int, range(1, NPIX + 1, length=nband + 1))
    fig = Figure(size=(2000, 1700))
    Label(fig[0, 1], @sprintf("%s fiber %s vs COMBINED mask (%s) — %s\nper-fiber %d px | combined %d px | OVER-masked %d px (%.2f%% of fiber's valid) | UNDER-masked %d px (%.2f%%) | net %+d | IoU %.4f",
            tele, n, VARIANT, why, npf, ncb, nov, 100nov / count(sub), nun, 100nun / count(sub),
            nov - nun, iou),
        fontsize=19, color=:white, tellwidth=false)

    for b in 1:nband
        lo, hi = edges[b], edges[b+1] - 1
        ax = Axis(fig[b, 1], xlabel=b == nband ? "wavelength [Å]" : "",
            ylabel="median sky line flux", yscale=Makie.pseudolog10,
            xgridvisible=false, ygridcolor=(:white, 0.12))
        xlims!(ax, wavetarg[lo], wavetarg[hi])
        for (msk, col) in ((both, (:gray70, 0.30)), (underm, (:cyan, 0.45)), (overm, (:magenta, 0.45)))
            for r in mask_runs(msk)
                (last(r) < lo || first(r) > hi) && continue
                a = max(first(r), lo); z = min(last(r), hi)
                vspan!(ax, wavetarg[a], wavetarg[min(z + 1, NPIX)], color=col)
            end
        end
        xs = wavetarg[lo:hi]
        lines!(ax, xs, replace(x[lo:hi], NaN => 0.0), color=(:white, 0.9), linewidth=0.6)
        lines!(ax, xs, KCUT .* replace(scale[lo:hi], NaN => 0.0), color=:orange, linewidth=1.3)
        # annotate the disagreement spans with their central wavelength
        ymax = maximum(filter(isfinite, x[lo:hi]); init=1.0)
        for (msk, col) in ((overm, :magenta), (underm, :cyan))
            for r in mask_runs(msk)
                (last(r) < lo || first(r) > hi) && continue
                length(r) < 3 && continue
                λc = wavetarg[clamp((first(r) + last(r)) ÷ 2, 1, NPIX)]
                text!(ax, λc, ymax, text=@sprintf("%.0f", λc), color=col, fontsize=10,
                    rotation=pi / 2, align=(:left, :center))
            end
        end
    end
    els = [PolyElement(color=(:gray70, 0.30)), PolyElement(color=(:magenta, 0.45)),
        PolyElement(color=(:cyan, 0.45)), LineElement(color=:orange), LineElement(color=(:white, 0.9))]
    Legend(fig[nband+1, 1], els,
        ["both agree", "per-fiber OVER-masked (this fiber flagged it, combined does NOT)",
            "per-fiber UNDER-masked (combined flags it, this fiber did NOT)",
            "detector cut = $(Int(KCUT)) × running MAD", "median sky line flux"],
        orientation=:horizontal, framevisible=false, labelsize=14, tellheight=true)
    out = joinpath(PD, "fig_bright_combined_diff_$(tele)_f$(n).png")
    save(out, fig, px_per_unit=2)
    return out, nov, nun
end

# pick the worst offenders in BOTH directions, from the measured per-fiber differences
over = [count(M[:, f] .& .!CB .& S[:, f]) for f in 1:600]
under = [count(CB .& .!M[:, f] .& S[:, f]) for f in 1:600]
picks = Tuple{Int,String}[]
for (lbl, v) in (("WORST OVER-masker", over), ("WORST UNDER-masker", under))
    for tele in ((1:300, "APO"), (301:600, "LCO"))
        rng, tn = tele
        o = sort(collect(rng), by=f -> -v[f])
        for (i, f) in enumerate(o[1:2])
            push!(picks, (f, "$lbl #$i at $tn (over=$(over[f]) px, under=$(under[f]) px)"))
        end
    end
end
seen = Set{Int}()
picks = [p for p in picks if !(p[1] in seen) && (push!(seen, p[1]); true)]
for (f, why) in picks
    out, nov, nun = make_diff_fig(f, why)
    @printf("wrote %s  (over=%d under=%d)\n", out, nov, nun); flush(stdout)
end
println("done")
