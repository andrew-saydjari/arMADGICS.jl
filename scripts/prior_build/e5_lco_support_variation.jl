## E5: WHY does the number of pixels removed from the faint support vary, and why more at LCO?
#
# AKS (2026-09-07): "why do the number of pixels removed from the faint support vary at LCO?
# Is that because there is a skyline at the edge of the wavelength coverage, which varies
# slightly across the chip and thus fiber number?"
#
# The combined bright mask is ONE fiber-independent 650-px array. What varies per fiber is
#   removed[f] = |bright ∩ submsk_f|
# so all variation comes from bright pixels that a fiber's own valid support does not cover.
# submsk_f factorises EXACTLY (build_sky_defs.jl:345-352):
#   submsk_f = (obscnt_f .>= 10) .& chipgapmsk[telescope(f)]
# chipgapmsk is per TELESCOPE (one array for all 300 APO fibers, one for all 300 LCO), so
# it can only produce a CONSTANT offset within a telescope. Every per-fiber difference must
# therefore come from obscnt_f -- the number of contributing exposures at that pixel.
#
# This script measures that decomposition and locates the varying pixels in wavelength
# relative to the chip footprints.
#
# Run:
#   julia --project=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E5_run/rebuild_qa \
#         -t 16 scripts/prior_build/e5_lco_support_variation.jl
# Author - Andrew Saydjari (E5)

using HDF5, Statistics, Printf, StatsBase
using CairoMakie, ColorSchemes
black_latexfonts = merge(theme_black(), theme_latexfonts())
set_theme!(black_latexfonts)
CairoMakie.disable_mime!("svg", "pdf", "text/html")

const NPIX = 8700
const APO = 1:300
const LCO = 301:600

base = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1/built"
newd = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1/built_combined_telemaj_union"
combined_path = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_07/e5_combined/e5_bright_combined.h5"
perfiber_path = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_07/e5_combined/e5_bright_perfiber.h5"
chipgap_path = "/mnt/ceph/users/sdssv/work/asaydjari/2026_04_25/StarContChipGapMsk.h5"
outdir = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/plots/e5_sky"
mkpath(outdir)

# same grid every script in this workstream declares (build_sky_defs.jl:21)
wavetarg = 10 .^ range((4.179 - 125 * 6.0e-6), step = 6.0e-6, length = 8575 + 125)
@assert length(wavetarg) == NPIX

# submsk is stored Int64 0/1, NOT Bool -- coerce on read (e5_rebuild_qa.jl:48-52)
readmsk(p) = h5read(p, "submsk") .!= 0
skyfaint_name(d, f) = joinpath(d, "APOGEE_skyline_faint_svd_120_f" * lpad(f, 3, "0") * ".h5")

rep = IOBuffer()
say(args...) = (s = string(args...); println(s); println(rep, s))

say("# E5 faint-support variation: decomposition and chip-edge test")
say("grid: ", @sprintf("%.3f", wavetarg[1]), " - ", @sprintf("%.3f", wavetarg[end]), " A, ", NPIX, " px")

bright = Bool.(h5read(combined_path, "mask_telemaj_union"))
say("bright mask telemaj_union: ", sum(bright), " px")
bidx = findall(bright)

chip_apo = Bool.(h5read(chipgap_path, "apo"))
chip_lco = Bool.(h5read(chipgap_path, "lco"))
say("chipgapmsk: apo ", sum(chip_apo), " px, lco ", sum(chip_lco), " px")

# ---------------------------------------------------------------- load all 600 submsk
S = falses(NPIX, 600)
have = falses(600)
Threads.@threads for f in 1:600
    p = skyfaint_name(base, f)
    isfile(p) || continue
    try
        S[:, f] = readmsk(p)
        have[f] = true
    catch
    end
end
say("baseline submsk loaded for ", count(have), "/600 fibers")
@assert count(have) == 600

# sanity: the baseline (legacy policy) submsk must equal (obscnt>=10) & chipgapmsk, so it
# must be a SUBSET of the telescope chip mask. This is what licenses the factorisation.
viol = sum(count(S[:, f] .& .!(f <= 300 ? chip_apo : chip_lco)) for f in 1:600)
say("submsk pixels outside the telescope chipgapmsk (must be 0): ", viol)

# ---------------------------------------------------------------- removed per fiber
removed = [count(bright .& S[:, f]) for f in 1:600]
say("")
say("## 1. removed = |bright ∩ submsk_f| (MEASURED, 600 fibers)")
say(@sprintf("all      : med %d  min %d  max %d  range %d", median(removed), minimum(removed), maximum(removed), maximum(removed) - minimum(removed)))
for (nm, rng) in (("APO", APO), ("LCO", LCO))
    r = removed[rng]
    say(@sprintf("%-9s: med %d  min %d  max %d  range %d  mean %.2f  sd %.2f  n_distinct %d",
        nm, median(r), minimum(r), maximum(r), maximum(r) - minimum(r), mean(r), std(r), length(unique(r))))
end

# ---------------------------------------------------------------- ceiling per telescope
ceil_apo = count(bright .& chip_apo)
ceil_lco = count(bright .& chip_lco)
say("")
say("## 2. Decomposition: chip footprint (constant within a telescope) vs obscnt (per fiber)")
say(@sprintf("|bright ∩ chipgapmsk_apo| = %d   (ceiling for every APO fiber; %d bright px are off the APO chips)", ceil_apo, sum(bright) - ceil_apo))
say(@sprintf("|bright ∩ chipgapmsk_lco| = %d   (ceiling for every LCO fiber; %d bright px are off the LCO chips)", ceil_lco, sum(bright) - ceil_lco))
for (nm, rng, cl) in (("APO", APO, ceil_apo), ("LCO", LCO, ceil_lco))
    d = cl .- removed[rng]
    say(@sprintf("%s: obscnt-driven shortfall below the ceiling: med %d  min %d  max %d  (mean %.2f)",
        nm, median(d), minimum(d), maximum(d), mean(d)))
end

# ---------------------------------------------------------------- which bright px are missed
say("")
say("## 3. WHICH bright pixels are missed, and by how many fibers")
nmiss_apo = [count(!, S[p, APO]) for p in bidx]     # of 300
nmiss_lco = [count(!, S[p, LCO]) for p in bidx]
for (nm, nm_arr, cg) in (("APO", nmiss_apo, chip_apo), ("LCO", nmiss_lco, chip_lco))
    allmiss = count(nm_arr .== 300)
    part = count(0 .< nm_arr .< 300)
    none = count(nm_arr .== 0)
    offchip = count(.!cg[bidx])
    say(@sprintf("%s: of %d bright px -> %d covered by ALL 300 fibers, %d missed by SOME (1-299), %d missed by ALL 300 (of which %d are off-chip)",
        nm, length(bidx), none, part, allmiss, offchip))
end

# the partially-missed pixels are the entire source of within-telescope variation
for (nm, nm_arr, rng, cg) in (("APO", nmiss_apo, APO, chip_apo), ("LCO", nmiss_lco, LCO, chip_lco))
    sel = findall(0 .< nm_arr .< 300)
    say("")
    say("### $nm partially-covered bright pixels (the ONLY source of per-fiber variation): $(length(sel))")
    if !isempty(sel)
        say(@sprintf("  %-6s %-12s %-8s %-9s %-10s", "pix", "wave[A]", "nmiss", "onchip?", "d_edge[px]"))
        # distance to the nearest chip-footprint boundary for this telescope
        runs = Tuple{Int,Int}[]
        let inrun = false, s0 = 0
            for i in 1:NPIX
                if cg[i] && !inrun
                    inrun = true; s0 = i
                elseif !cg[i] && inrun
                    inrun = false; push!(runs, (s0, i - 1))
                end
            end
            inrun && push!(runs, (s0, NPIX))
        end
        edges = vcat([[a, b] for (a, b) in runs]...)
        for k in sel
            p = bidx[k]
            de = minimum(abs.(p .- edges))
            say(@sprintf("  %-6d %-12.3f %-8d %-9s %-10d", p, wavetarg[p], nm_arr[k], cg[p] ? "yes" : "NO", de))
        end
        des = [minimum(abs.(bidx[k] .- edges)) for k in sel]
        say(@sprintf("  distance to nearest chip-footprint edge: med %d px, min %d, max %d", median(des), minimum(des), maximum(des)))
        # null: the same statistic for ALL bright pixels on-chip
        allon = [p for p in bidx if cg[p]]
        desall = [minimum(abs.(p .- edges)) for p in allon]
        say(@sprintf("  NULL (all %d on-chip bright px): med %d px, min %d, max %d", length(allon), median(desall), minimum(desall), maximum(desall)))
    end
end

# ---------------------------------------------------------------- which bright px are off-chip
say("")
say("## 3b. Bright pixels OFF a telescope's chip footprint (constant loss, not variation)")
for (nm, cg) in (("APO", chip_apo), ("LCO", chip_lco))
    off = [p for p in bidx if !cg[p]]
    if isempty(off)
        say("$nm: none -- all 650 bright px sit on the $nm chip footprint")
    else
        say(@sprintf("%s: %d px, grid %d-%d, %.3f-%.3f A", nm, length(off), minimum(off), maximum(off),
            wavetarg[minimum(off)], wavetarg[maximum(off)]))
        say("  pixels: ", join(off, ","))
    end
end

# --------------------------------------------- DIRECT test: per-fiber chip-edge position
# If AKS is right, the variation is a chip EDGE that sits at a slightly different grid pixel
# for each fiber. Measure that edge directly: within the window around the varying region,
# find each fiber's bluest supported pixel and check the support is a clean truncation
# (all-false then all-true), which is the signature of an edge and not of speckled dropout.
say("")
say("## 3c. DIRECT measurement of the per-fiber chip edge in the varying window")
const WLO, WHI = 6400, 6500
edgepix = fill(0, 600)
clean = falses(600)
for f in 1:600
    w = S[WLO:WHI, f]
    i = findfirst(w)
    if i === nothing
        edgepix[f] = 0
    else
        edgepix[f] = WLO + i - 1
        clean[f] = all(w[i:end])   # contiguous from the edge to the end of the window
    end
end
for (nm, rng) in (("APO", APO), ("LCO", LCO))
    e = edgepix[rng]
    say(@sprintf("%s: bluest supported pixel in %d:%d -> med %d  min %d  max %d  spread %d px (%.2f A)  n_distinct %d  clean truncation %d/300",
        nm, WLO, WHI, median(e), minimum(e), maximum(e), maximum(e) - minimum(e),
        wavetarg[maximum(e)] - wavetarg[minimum(e)], length(unique(e)), count(clean[rng])))
end
let e = edgepix[LCO]
    say(@sprintf("LCO edge vs fiber number: Spearman = %+.3f   Pearson = %+.3f", corspearman(float.(1:300), float.(e)), cor(float.(1:300), float.(e))))
    say(@sprintf("LCO removed vs LCO edge pixel: Spearman = %+.3f (must be ~-1 if the edge alone explains it)",
        corspearman(float.(e), float.(removed[LCO]))))
    say(@sprintf("removed[f] == %d - (edge[f] - %d) exactly for %d/300 LCO fibers",
        ceil_lco, minimum(e), count(removed[LCO] .== ceil_lco .- (e .- minimum(e)))))
end

# ---------------------------------------------------------------- fiber-number trend
say("")
say("## 4. Does removed track fiber number? (a smooth detector trend vs scatter)")
for (nm, rng) in (("APO", APO), ("LCO", LCO))
    r = removed[rng]
    x = collect(1:300)
    sp = corspearman(float.(x), float.(r))
    pe = cor(float.(x), float.(r))
    # neighbour roughness: |r[i+1]-r[i]| vs the spread. Smooth trend => small steps.
    d1 = abs.(diff(r))
    say(@sprintf("%s: Spearman(fiber, removed) = %+.3f   Pearson = %+.3f   median |Δ neighbour| = %.1f px   sd = %.2f px",
        nm, sp, pe, median(d1), std(r)))
end
# structure of the departure: is it a gentle global gradient or confined to one end?
let r = removed[LCO], cl = ceil_lco
    atceil = findall(r .== cl)
    below = findall(r .< cl)
    say(@sprintf("LCO: %d/300 fibers sit exactly at the ceiling (%d px); %d fall short.", length(atceil), cl, length(below)))
    say(@sprintf("     first LCO fiber below the ceiling: LCO#%d (f%03d); last fiber at the ceiling: LCO#%d (f%03d)",
        minimum(below), 300 + minimum(below), maximum(atceil), 300 + maximum(atceil)))
    say(@sprintf("     of the %d short fibers, %d have LCO index > %d -- the shortfall is confined to one end of the fiber array, not spread over it",
        length(below), count(below .> 230), 230))
end

# per-fiber list of the worst (most-missing) fibers
say("")
say("## 5. Fibers with the largest shortfall")
for (nm, rng, cl) in (("APO", APO, ceil_apo), ("LCO", LCO, ceil_lco))
    d = cl .- removed[rng]
    ord = sortperm(d, rev = true)[1:min(10, length(d))]
    say("$nm worst 10 (fiber, shortfall below ceiling, removed): ",
        join([string("f", lpad(first(rng) + i - 1, 3, "0"), "(", d[i], ",", removed[first(rng)+i-1], ")") for i in ord], " "))
end

# ---------------------------------------------------------------- FIGURES
cA = colorant"#4FC3F7"
cL = colorant"#FF7043"

# Fig 1: removed vs fiber number, per telescope, plus the directly measured chip edge
fig1 = Figure(size = (1150, 900))
for (i, (nm, rng, cl, col)) in enumerate((("APO (fibers 1-300)", APO, ceil_apo, cA), ("LCO (fibers 301-600)", LCO, ceil_lco, cL)))
    ax = Axis(fig1[i, 1], xlabel = "",
        ylabel = "bright px removed",
        title = nm * @sprintf("   ceiling |bright ∩ chipmask| = %d px", cl))
    scatter!(ax, collect(rng), removed[rng], color = col, markersize = 5)
    hlines!(ax, [cl], color = :white, linestyle = :dash, linewidth = 1.2)
    text!(ax, first(rng) + 4, cl, text = "chip-footprint ceiling", align = (:left, :bottom), color = :white, fontsize = 11)
    ylims!(ax, minimum(removed) - 3, 653)
    xlims!(ax, first(rng) - 3, last(rng) + 3)
end
let ax = Axis(fig1[3, 1], xlabel = "adjusted fiber index", ylabel = "bluest supported pixel",
        title = "Directly measured chip edge: each fiber's bluest valid pixel in grid window $WLO:$WHI (LCO red-chip blue edge)")
    scatter!(ax, collect(LCO), edgepix[LCO], color = cL, markersize = 5)
    scatter!(ax, collect(APO), edgepix[APO], color = cA, markersize = 5)
    text!(ax, 5, edgepix[1], text = "APO: fixed at $(edgepix[1])", align = (:left, :bottom), color = cA, fontsize = 11)
    xlims!(ax, -3, 603)
end
Label(fig1[0, 1], "Pixels removed from the faint sky-line support, per fiber\ncombined:telemaj_union (650 px, fiber-independent) ∩ each fiber's own submsk",
    fontsize = 15, tellwidth = false)
save(joinpath(outdir, "fig_lco_removed_vs_fiber.png"), fig1, px_per_unit = 2)
say("")
say("wrote fig_lco_removed_vs_fiber.png")

# Fig 2: wavelength location of the partially-covered bright pixels, chip footprints shaded
fig2 = Figure(size = (1250, 760))
for (i, (nm, nm_arr, cg, col)) in enumerate((("APO", nmiss_apo, chip_apo, cA), ("LCO", nmiss_lco, chip_lco, cL)))
    ax = Axis(fig2[i, 1], xlabel = i == 2 ? "wavelength [Å]" : "",
        ylabel = "fibers WITHOUT support (of 300)", title = nm,
        yscale = identity)
    # shade the chip footprints
    runs = Tuple{Int,Int}[]
    let inrun = false, s0 = 0
        for j in 1:NPIX
            if cg[j] && !inrun
                inrun = true; s0 = j
            elseif !cg[j] && inrun
                inrun = false; push!(runs, (s0, j - 1))
            end
        end
        inrun && push!(runs, (s0, NPIX))
    end
    for (a, b) in runs
        vspan!(ax, wavetarg[a], wavetarg[b], color = (:white, 0.07))
        vlines!(ax, [wavetarg[a], wavetarg[b]], color = (:white, 0.55), linestyle = :dot, linewidth = 1.0)
    end
    # all bright pixels at y=0 for context
    scatter!(ax, wavetarg[bidx], fill(0.0, length(bidx)), color = (:grey, 0.45), markersize = 3)
    sel = findall(nm_arr .> 0)
    if !isempty(sel)
        stem!(ax, wavetarg[bidx[sel]], float.(nm_arr[sel]), color = col, markersize = 7,
            stemcolor = col, stemwidth = 1.5)
    end
    ylims!(ax, -12, 320)
    xlims!(ax, wavetarg[1], wavetarg[end])
    # zoom on the varying region, same construction
    az = Axis(fig2[i, 2], xlabel = i == 2 ? "wavelength [Å]" : "", title = nm * " zoom: 16465-16490 Å")
    for (a, b) in runs
        vspan!(az, wavetarg[a], wavetarg[b], color = (:white, 0.07))
        vlines!(az, [wavetarg[a], wavetarg[b]], color = (:white, 0.55), linestyle = :dot, linewidth = 1.2)
    end
    scatter!(az, wavetarg[bidx], fill(0.0, length(bidx)), color = (:grey, 0.6), markersize = 5)
    if !isempty(sel)
        stem!(az, wavetarg[bidx[sel]], float.(nm_arr[sel]), color = col, markersize = 7,
            stemcolor = col, stemwidth = 1.5)
    end
    ylims!(az, -12, 320)
    xlims!(az, 16465.0, 16490.0)
end
colsize!(fig2.layout, 1, Relative(0.66))
Label(fig2[0, 1:2], "Where the 650 bright pixels lose per-fiber support\nshaded = that telescope's chip footprint (StarContChipGapMsk); dotted = footprint edges; grey dots = all bright px",
    fontsize = 14, tellwidth = false)
save(joinpath(outdir, "fig_lco_missing_support_wave.png"), fig2, px_per_unit = 2)
say("wrote fig_lco_missing_support_wave.png")

open(joinpath(outdir, "lco_support_variation.txt"), "w") do io
    write(io, String(take!(rep)))
end
println("report -> ", joinpath(outdir, "lco_support_variation.txt"))
