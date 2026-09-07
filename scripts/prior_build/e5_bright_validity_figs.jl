## E5: figures for the VALIDITY-CORRECT combined bright mask (AKS 2026-09-07, 2nd pass).
#
#   fig_bright_validity_chipedge.png  — the chip-edge demonstration AKS asked for:
#       a region where LCO fibers HAVE data and APO fibers have NONE, showing that an LCO
#       fiber which has data there now gets its bright pixel, where the all-600 majority
#       rule made that arithmetically impossible.
#   fig_bright_validity_lco.png       — the LCO complaint diagnosed: per-telescope
#       detection FRACTION over ELIGIBLE fibers, with every pixel the old rule lost.
#
# Dark theme, log-y, mutually exclusive colours, same conventions as the other E5 figures.
# Usage: julia --project=<plots env> e5_bright_validity_figs.jl
# Author - Andrew Saydjari (E5 pass 1)
using HDF5, Statistics, Printf
using CairoMakie, ColorSchemes
set_theme!(merge(theme_black(), theme_latexfonts()))
CairoMakie.disable_mime!("svg", "pdf", "text/html")

const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
include(ARM * "scripts/prior_build/e5_bright_line_detect.jl")
const SC = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1/screens"
const RESDIR = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_07/e5_combined"
const PD = get(ENV, "E5_PLOTDIR", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/plots/e5_sky")
mkpath(PD)
const NPIX = 8575 + 125
const SCALE_WINDOW, KCUT = 2001, 90.0
wavetarg = 10 .^ range((4.179 - 125 * 6.0e-6), step=6.0e-6, length=NPIX)

C = joinpath(RESDIR, "e5_bright_combined.h5")
NEW = Bool.(h5read(C, "mask_telemaj_union"))
OLD = Bool.(h5read(C, "mask_majority"))
nvA = h5read(C, "nval_apo"); nvL = h5read(C, "nval_lco")
frA = h5read(C, "frac_apo"); frL = h5read(C, "frac_lco")   # -1 where undefined
S = Bool.(h5read(joinpath(RESDIR, "e5_bright_perfiber.h5"), "submsk"))

function fiber_spec(adjfib)
    n = lpad(adjfib, 3, "0")
    ms = h5read(joinpath(SC, "median_sky_$n.h5"), "median_sky")
    sub = Bool.(h5read(joinpath(SC, "median_sky_$n.h5"), "submsk"))
    x = fill(NaN, NPIX); x[sub] .= ms
    return x, running_spread_fast(x, SCALE_WINDOW; kind=:mad, stride=50), sub
end

## ------------------------------------------------------------- 1. chip-edge demonstration
# the LCO-only-coverage pixels: LCO fibers have data, APO fibers have none
lco_only = [nvL[p] > 0 && nvA[p] == 0 for p in 1:NPIX]
gained = NEW .& .!OLD
# centre the figure on the LCO-only region that the new mask recovers
target = findfirst(p -> lco_only[p] && gained[p], 1:NPIX)
target === nothing && (target = findfirst(lco_only))
lo = max(target - 260, 1); hi = min(target + 260, NPIX)
# an LCO fiber that HAS data at that pixel, and an APO fiber (which cannot)
lcofib = findfirst(f -> S[target, f], 301:600) + 300
apofib = 150
xl, scl, subl = fiber_spec(lcofib)
xa, sca, suba = fiber_spec(apofib)

fig = Figure(size=(2000, 1450))
Label(fig[0, 1], @sprintf("Chip-edge demonstration: coverage must NOT be baked into the bright mask\nregion %.1f-%.1f Å — %d pixels here have LCO data and NO APO data at all",
        wavetarg[lo], wavetarg[hi], count(view(lco_only, lo:hi))),
    fontsize=22, color=:white, tellwidth=false)

# (a) coverage
ax1 = Axis(fig[1, 1], ylabel="fibers with valid data", xgridvisible=false,
    title="(a) eligible-fiber count per telescope — APO coverage ENDS inside this window",
    ygridcolor=(:white, 0.12))
lines!(ax1, wavetarg[lo:hi], Float64.(nvA[lo:hi]), color=:orange, linewidth=2.2, label="APO eligible (of 300)")
lines!(ax1, wavetarg[lo:hi], Float64.(nvL[lo:hi]), color=:cyan, linewidth=2.2, label="LCO eligible (of 300)")
vlines!(ax1, [wavetarg[target]], color=:magenta, linewidth=1.6, linestyle=:dash)
axislegend(ax1, position=:lc, labelsize=13, framevisible=false)
xlims!(ax1, wavetarg[lo], wavetarg[hi]); ylims!(ax1, -15, 330)

# (b) the LCO fiber that HAS data there
ax2 = Axis(fig[2, 1], ylabel="median sky line flux", yscale=Makie.pseudolog10,
    xgridvisible=false, ygridcolor=(:white, 0.12),
    title=@sprintf("(b) LCO fiber %d — it HAS data across the edge, and now gets to use it", lcofib))
for (msk, col) in ((OLD .& NEW, (:gray70, 0.30)), (gained, (:magenta, 0.55)))
    for r in mask_runs(BitVector(msk))
        (last(r) < lo || first(r) > hi) && continue
        vspan!(ax2, wavetarg[max(first(r), lo)], wavetarg[min(last(r) + 1, hi)], color=col)
    end
end
lines!(ax2, wavetarg[lo:hi], replace(xl[lo:hi], NaN => 0.0), color=(:white, 0.9), linewidth=1.0)
lines!(ax2, wavetarg[lo:hi], KCUT .* replace(scl[lo:hi], NaN => 0.0), color=:orange, linewidth=1.4)
vlines!(ax2, [wavetarg[target]], color=:magenta, linewidth=1.6, linestyle=:dash)
xlims!(ax2, wavetarg[lo], wavetarg[hi])

# (c) an APO fiber over the same range — no data past the edge
ax3 = Axis(fig[3, 1], xlabel="wavelength [Å]", ylabel="median sky line flux",
    yscale=Makie.pseudolog10, xgridvisible=false, ygridcolor=(:white, 0.12),
    title=@sprintf("(c) APO fiber %d — its data STOPS at the edge; under the old all-600 rule its absence outvoted LCO", apofib))
# shade the APO-invalid stretches as spans (a band! with +/-1e9 sentinels distorts the
# pseudolog axis by 8 decades and makes the real flux unreadable)
for r in mask_runs(BitVector([!suba[p] for p in 1:NPIX]))
    (last(r) < lo || first(r) > hi) && continue
    vspan!(ax3, wavetarg[max(first(r), lo)], wavetarg[min(last(r) + 1, hi)], color=(:red, 0.18))
end
lines!(ax3, wavetarg[lo:hi], replace(xa[lo:hi], NaN => 0.0), color=(:white, 0.7), linewidth=1.0)
vlines!(ax3, [wavetarg[target]], color=:magenta, linewidth=1.6, linestyle=:dash)
xlims!(ax3, wavetarg[lo], wavetarg[hi])
let fin = filter(isfinite, xa[lo:hi])
    isempty(fin) || ylims!(ax3, -2 * max(abs(minimum(fin)), 1.0), 1.35 * max(maximum(fin), 1.0))
end
text!(ax3, wavetarg[clamp(target - 230, 1, NPIX)], 0.0,
    text="APO has NO valid pixels here (red)", color=(:red, 0.95), fontsize=15,
    align=(:left, :center))

Legend(fig[4, 1],
    [PolyElement(color=(:gray70, 0.30)), PolyElement(color=(:magenta, 0.55)),
        PolyElement(color=(:red, 0.16)), LineElement(color=:orange),
        LineElement(color=:magenta, linestyle=:dash)],
    ["bright in BOTH old and new mask", "RECOVERED by the new per-telescope union",
        "no APO coverage", "detector cut = 90 × running MAD (that fiber)",
        "the recovered LCO-only pixel"],
    orientation=:horizontal, framevisible=false, labelsize=14, tellheight=true)
save(joinpath(PD, "fig_bright_validity_chipedge.png"), fig, px_per_unit=2)
println("wrote fig_bright_validity_chipedge.png  (target ", round(wavetarg[target], digits=2),
    " Å, LCO fiber ", lcofib, ")"); flush(stdout)

## ------------------------------------------------------------------ 2. the LCO diagnosis
fig2 = Figure(size=(2000, 1150))
Label(fig2[0, 1], "Why LCO lines were being lost: detection fraction over ELIGIBLE fibers, per telescope",
    fontsize=22, color=:white, tellwidth=false)

ax = Axis(fig2[1, 1], xlabel="wavelength [Å]", ylabel="fraction of ELIGIBLE fibers flagging",
    xgridvisible=false, ygridcolor=(:white, 0.12),
    title="APO vs LCO detection fraction; magenta ticks = pixels the old all-600 majority lost")
lines!(ax, wavetarg, [frA[p] < 0 ? NaN : frA[p] for p in 1:NPIX], color=(:orange, 0.85), linewidth=0.8)
lines!(ax, wavetarg, [frL[p] < 0 ? NaN : frL[p] for p in 1:NPIX], color=(:cyan, 0.85), linewidth=0.8)
hlines!(ax, [0.5], color=(:white, 0.5), linewidth=1.2, linestyle=:dash)
for p in findall(gained)
    vlines!(ax, [wavetarg[p]], color=(:magenta, 0.85), linewidth=1.2)
end
axislegend(ax, [LineElement(color=(:orange, 0.85)), LineElement(color=(:cyan, 0.85)),
        LineElement(color=(:white, 0.5), linestyle=:dash), LineElement(color=(:magenta, 0.85))],
    ["APO fraction of eligible", "LCO fraction of eligible", "majority cut (0.5)",
        "recovered by the new rule"], position=:lt, labelsize=13, framevisible=false)
xlims!(ax, wavetarg[1], wavetarg[end]); ylims!(ax, -0.03, 1.08)

# scatter: APO vs LCO fraction, recovered pixels highlighted
ax2b = Axis(fig2[2, 1], xlabel="APO detection fraction (of eligible)",
    ylabel="LCO detection fraction (of eligible)",
    title="each flagged pixel; points in the upper-left are lines LCO sees and APO does not")
sel = [(frA[p] >= 0 || frL[p] >= 0) && (max(frA[p], frL[p]) > 0.02) for p in 1:NPIX]
scatter!(ax2b, [frA[p] < 0 ? 0.0 : frA[p] for p in findall(sel)],
    [frL[p] < 0 ? 0.0 : frL[p] for p in findall(sel)], markersize=5, color=(:gray70, 0.5))
scatter!(ax2b, [frA[p] < 0 ? 0.0 : frA[p] for p in findall(gained)],
    [frL[p] < 0 ? 0.0 : frL[p] for p in findall(gained)], markersize=11, color=:magenta)
hlines!(ax2b, [0.5], color=(:white, 0.4), linewidth=1, linestyle=:dash)
vlines!(ax2b, [0.5], color=(:white, 0.4), linewidth=1, linestyle=:dash)
lines!(ax2b, [0, 1], [0, 1], color=(:orange, 0.6), linewidth=1)
save(joinpath(PD, "fig_bright_validity_lco.png"), fig2, px_per_unit=2)
println("wrote fig_bright_validity_lco.png"); flush(stdout)
println("done")
