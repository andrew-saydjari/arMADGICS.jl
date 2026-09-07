## E5: COMBINE the bright sky-line mask across fibers (AKS directive 2026-09-07).
#
# AKS: "it is kind of weird to have different wavelengths masked as 'bright' skylines
# depending on the fiber number. Because we are only applying it as a mask, I think we
# should now combine the information across all 600 [fibers] and come up with a good mask
# of the lines that register as bright for many/most of the fibers. ... The most
# conservative would be to mask any pixel masked as a bright line in any fiber. We *could*
# do that, and should evaluate how different that is from something slightly smarter that
# might reject pixels masked only in 1-2 fibers."
#
# Stage 1 (--stage=masks): run the calibrated `linedetect` detector per fiber on all 600
#   cached median_sky arrays; write the per-fiber full-grid masks + the per-pixel
#   detection COUNT to e5_bright_perfiber.h5.
# Stage 2 (--stage=eval): build the combined-mask variants (union / drop-1-2 / majority /
#   ...), score every one against DR17's own bright masks, and measure per-telescope vs
#   all-600 combination.
#
# Everything is defined on the FULL wavetarg grid (8700 px). A per-fiber mask is false
# outside that fiber's own `submsk`, so the detection count is accompanied by a per-pixel
# VALID count (how many fibers could have flagged it at all).
#
# Usage: julia --project=<arM-E5b> e5_bright_combine.jl --stage=masks|eval [--nworkers=N]
# Author - Andrew Saydjari (E5 pass 1)
using Distributed, Printf, Dates

const STAGE = let s = filter(a -> startswith(a, "--stage="), ARGS)
    isempty(s) ? "eval" : String(split(s[1], "=")[2])
end
const NW = parse(Int, replace(get(filter(a -> startswith(a, "--nworkers="), ARGS), 1, "--nworkers=16"), "--nworkers=" => ""))
const RESDIR = get(ENV, "E5_RESDIR", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_07/e5_combined")
const PERFIB = joinpath(RESDIR, "e5_bright_perfiber.h5")
mkpath(RESDIR)

# ------------------------------------------------------------------ stage 1: per-fiber
if STAGE == "masks"
    addprocs(NW)
    @everywhere begin
        using HDF5, Statistics
        const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
        include(ARM * "scripts/prior_build/e5_bright_line_detect.jl")
        const SC = "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1/screens"
        # the calibrated `linedetect` policy (BRIGHT_MASK_CALIBRATION.md)
        const SCALE_WINDOW, KCUT, DILATION = 2001, 90.0, 4
        const NPIX = 8575 + 125

        function fiber_mask(adjfib)
            n = lpad(adjfib, 3, "0")
            p = joinpath(SC, "median_sky_$n.h5")
            ms = h5read(p, "median_sky"); sub = Bool.(h5read(p, "submsk"))
            x = fill(NaN, NPIX); x[sub] .= ms
            scale = running_spread_fast(x, SCALE_WINDOW; kind=:mad, stride=50)
            base = falses(NPIX)
            @inbounds for i in 1:NPIX
                (isfinite(x[i]) && isfinite(scale[i]) && scale[i] > 0 && x[i] > KCUT * scale[i]) && (base[i] = true)
            end
            m = dilate_msk(base, DILATION)
            @inbounds for i in 1:NPIX; isfinite(x[i]) || (m[i] = false); end
            return Vector{Bool}(m), Vector{Bool}(sub)
        end
    end
    println("running linedetect on 600 fibers ($NW workers)  $(now())"); flush(stdout)
    t0 = time()
    res = pmap(fiber_mask, 1:600)
    NPIX = 8575 + 125
    M = falses(NPIX, 600); S = falses(NPIX, 600)
    for (i, (m, s)) in enumerate(res); M[:, i] .= m; S[:, i] .= s; end
    h5open(PERFIB, "w") do fh
        fh["mask"] = Matrix{Bool}(M)          # per-fiber bright mask, full grid
        fh["submsk"] = Matrix{Bool}(S)        # per-fiber valid domain
        fh["nflag"] = vec(sum(M, dims=2))
        fh["nvalid"] = vec(sum(S, dims=2))
        attrs(fh)["policy"] = "linedetect scale_window=2001 k=90.0 dilation=4 cont_window=0"
    end
    @printf("wrote %s  (%.1f min)\n", PERFIB, (time() - t0) / 60)
    exit(0)
end

# ---------------------------------------------------------------------- stage 2: eval
using HDF5, Statistics
const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
include(ARM * "scripts/prior_build/e5_bright_line_detect.jl")
const DR = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/sky_priors"
const NPIX = 8575 + 125

M = Bool.(h5read(PERFIB, "mask")); S = Bool.(h5read(PERFIB, "submsk"))
nflag = h5read(PERFIB, "nflag"); nvalid = h5read(PERFIB, "nvalid")
APO = 1:300; LCO = 301:600

# DR17 per-fiber bright / valid
DB = falses(NPIX, 600); DV = falses(NPIX, 600)
for f in 1:600
    n = lpad(f, 3, "0")
    db = Bool.(h5read(joinpath(DR, "APOGEE_skyline_bright_svd_120_f$n.h5"), "submsk"))
    df = Bool.(h5read(joinpath(DR, "APOGEE_skyline_faint_svd_120_f$n.h5"), "submsk"))
    DB[:, f] .= db; DV[:, f] .= S[:, f] .& (db .| df)
end
dr_nflag = vec(sum(DB, dims=2))

io = open(joinpath(RESDIR, "combine_report.txt"), "w")
function emit(fmt, args...)
    s = Printf.format(Printf.Format(fmt), args...)
    print(io, s); print(stdout, s); flush(io); flush(stdout)
end

emit("E5 bright sky-line mask: COMBINED across fibers (%s)\n", string(now()))
emit("detector = calibrated linedetect (running MAD 2001 px, k=90, dilation 4), per fiber,\n")
emit("then combined. All numbers MEASURED unless marked INFERRED.\n\n")

# ------------------------------------------------------- per-pixel detection count
emit("## 1. Per-pixel detection count across fibers (MEASURED)\n\n")
emit("pixels with nflag>0: %d of %d (%.2f%% of the grid)\n",
    count(>(0), nflag), NPIX, 100count(>(0), nflag) / NPIX)
emit("nvalid: min=%d med=%d max=%d  (fibers in which a pixel is even eligible)\n\n",
    minimum(nvalid), Int(median(nvalid)), maximum(nvalid))
emit("%-14s %10s %10s\n", "nflag range", "pixels", "cum>=lo")
for (lo, hi) in ((1, 2), (3, 9), (10, 29), (30, 99), (100, 299), (300, 499), (500, 579), (580, 599), (600, 600))
    c = count(p -> lo <= nflag[p] <= hi, 1:NPIX)
    emit("%-14s %10d %10d\n", "$lo-$hi", c, count(>=(lo), nflag))
end
nf_pos = filter(>(0), nflag)
emit("\namong flagged pixels: median nflag=%d, %.1f%% have nflag>=580, %.1f%% have nflag<=2\n\n",
    Int(median(nf_pos)), 100count(>=(580), nf_pos) / length(nf_pos),
    100count(<=(2), nf_pos) / length(nf_pos))

# nvalid varies per pixel, so an ABSOLUTE count is not comparable across the grid: a pixel
# valid in only 300 fibers can never reach nflag=600. The FRACTION of eligible fibers is
# the scale-free version, and is the quantity that tests AKS's premise ("there should be
# very little variation on the mask around a given line as long as it is detected").
fr = [nvalid[p] > 0 ? nflag[p] / nvalid[p] : 0.0 for p in 1:NPIX]
frp = [fr[p] for p in 1:NPIX if nflag[p] > 0]
emit("detection FRACTION (nflag/nvalid) among flagged pixels:\n")
emit("%-16s %10s %10s\n", "fraction", "pixels", "% of flagged")
for (lo, hi) in ((0.0, 0.01), (0.01, 0.05), (0.05, 0.25), (0.25, 0.5), (0.5, 0.75),
    (0.75, 0.95), (0.95, 0.999), (0.999, 1.01))
    c = count(v -> lo <= v < hi, frp)
    emit("%-16s %10d %9.2f%%\n", "$lo-$hi", c, 100c / length(frp))
end
emit("=> BIMODALITY (tests AKS's premise): %.2f%% of flagged pixels are flagged by >=95%% of\n",
    100count(>=(0.95), frp) / length(frp))
emit("   eligible fibers, %.2f%% by <=5%%; only %.2f%% sit in the ambiguous 5-95%% middle.\n\n",
    100count(<=(0.05), frp) / length(frp),
    100count(v -> 0.05 < v < 0.95, frp) / length(frp))

# ------------------------------------------------------- combined variants
"combined mask from a detection-count vector at minimum count m"
cmb(nf, m) = nf .>= m

"per-fiber scoring of a combined mask C against DR17 bright, over fibers `rng`"
function score(C, rng)
    ious = Float64[]; recs = Float64[]; precs = Float64[]; fracs = Float64[]
    for f in rng
        sub = S[:, f]
        st = line_overlap_stats(DB[:, f], C .& sub, DV[:, f])
        push!(ious, st.pixel_iou); push!(recs, st.recall); push!(precs, st.precision)
        push!(fracs, count(C .& sub) / count(sub))
    end
    mn(v) = mean(filter(isfinite, v))
    return (iou=mn(ious), recall=mn(recs), prec=mn(precs), frac=100mn(fracs))
end
nlines(C) = length(mask_runs(C))

emit("## 2. Combined-mask variants (MEASURED)\n\n")
emit("Combination over ALL 600 fibers. `frac` = mean over fibers of masked/valid pixels.\n")
emit("DR17 columns = mean over the 600 fibers of that fiber's score vs its own DR17 bright mask.\n\n")
emit("%-26s %8s %7s %7s %8s %8s %8s\n", "variant", "pixels", "lines", "frac%", "IoU", "recall", "prec")
variants = [("UNION (nflag>=1)", 1), ("nflag>=2", 2), ("drop 1-2 (nflag>=3)", 3),
    ("nflag>=5", 5), ("nflag>=10", 10), ("nflag>=30 (5%)", 30), ("nflag>=60 (10%)", 60),
    ("nflag>=150 (25%)", 150), ("MAJORITY (nflag>=300)", 300), ("nflag>=450 (75%)", 450),
    ("nflag>=540 (90%)", 540), ("nflag>=570 (95%)", 570), ("UNANIMOUS (nflag>=600)", 600)]
for (nm, m) in variants
    C = cmb(nflag, m); s = score(C, 1:600)
    emit("%-26s %8d %7d %7.2f %8.3f %8.3f %8.3f\n", nm, count(C), nlines(C), s.frac, s.iou, s.recall, s.prec)
end
emit("\n-- scale-free variants: flagged by >= F of the fibers where the pixel is VALID --\n")
for F in (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
    C = BitVector([nvalid[p] > 0 && nflag[p] / nvalid[p] >= F for p in 1:NPIX])
    s = score(C, 1:600)
    emit("%-26s %8d %7d %7.2f %8.3f %8.3f %8.3f\n", "frac>=$F", count(C), nlines(C),
        s.frac, s.iou, s.recall, s.prec)
end
# the per-fiber product, for reference
pf = [(f, count(M[:, f])) for f in 1:600]
emit("%-26s %8s %7s %7.2f %8.3f %8.3f %8.3f\n", "[per-fiber, reference]",
    "~" * string(Int(round(median([p[2] for p in pf])))), "-",
    100mean(count(M[:, f]) / count(S[:, f]) for f in 1:600),
    mean(filter(isfinite, [line_overlap_stats(DB[:, f], M[:, f], DV[:, f]).pixel_iou for f in 1:600])),
    mean(filter(isfinite, [line_overlap_stats(DB[:, f], M[:, f], DV[:, f]).recall for f in 1:600])),
    mean(filter(isfinite, [line_overlap_stats(DB[:, f], M[:, f], DV[:, f]).precision for f in 1:600])))

# ------------------------------------------------------- union vs drop-1-2, specifically
emit("\n## 3. UNION vs drop-1-2 — the comparison AKS asked for (MEASURED)\n\n")
U = cmb(nflag, 1); D3 = cmb(nflag, 3)
diff = U .& .!D3
emit("UNION pixels           : %d (%d lines)\n", count(U), nlines(U))
emit("drop-1-2 pixels        : %d (%d lines)\n", count(D3), nlines(D3))
emit("pixels in UNION only   : %d  (%.3f%% of UNION, %.4f%% of the 8700-px grid)\n",
    count(diff), 100count(diff) / count(U), 100count(diff) / NPIX)
emit("lines removed          : %d  (UNION %d -> drop-1-2 %d)\n", nlines(U) - nlines(D3), nlines(U), nlines(D3))
dr = mask_runs(diff)
emit("the UNION-only pixels form %d runs, lengths: min=%d med=%d max=%d\n",
    length(dr), minimum(length, dr), Int(median(length.(dr))), maximum(length, dr))
emit("of those runs, %d are entirely detached from drop-1-2 (whole spurious lines)\n",
    count(r -> !any(D3[r]), dr))
emit("=> The UNION-only content is %s.\n",
    all(r -> !any(D3[r]), dr) && median(length.(dr)) <= 2 ?
    "ENTIRELY DETACHED, mostly single-pixel specks — 1-2-fiber noise hits, not line wings" :
    "attached to real lines — dropping it would trim genuine line wings")
emit("   (a pixel-count difference alone is NOT the deciding number: what matters is whether\n")
emit("    the extra pixels are line WINGS, which must be kept, or detached SPECKS, which must\n")
emit("    not. Detached single-pixel runs cannot be sky-line wings.)\n\n")

# ------------------------------------------------------- per-telescope vs all-600
emit("## 4. AMBIGUITY: per-telescope vs all-600 combination (MEASURED — AKS's call)\n\n")
nflag_apo = vec(sum(M[:, APO], dims=2)); nflag_lco = vec(sum(M[:, LCO], dims=2))
for (mname, mA, mL, mAll) in (("UNION", 1, 1, 1), ("drop-1-2", 3, 3, 3),
    ("majority(50%)", 150, 150, 300), ("90%", 270, 270, 540))
    CA = cmb(nflag_apo, mA); CL = cmb(nflag_lco, mL); CAll = cmb(nflag, mAll)
    both = CA .& CL; ao = CA .& .!CL; lo = CL .& .!CA
    emit("-- %s --\n", mname)
    emit("   APO-combined %d px (%d lines) | LCO-combined %d px (%d lines) | all-600 %d px (%d lines)\n",
        count(CA), nlines(CA), count(CL), nlines(CL), count(CAll), nlines(CAll))
    emit("   APO-only %d px, LCO-only %d px, shared %d px  -> Jaccard(APO,LCO) = %.4f\n",
        count(ao), count(lo), count(both), count(both) / count(CA .| CL))
    emit("   all-600 vs (APO on APO fibers, LCO on LCO fibers): APO fibers differ by %d px, LCO by %d px\n",
        count(xor.(CAll, CA)), count(xor.(CAll, CL)))
    sA = score(CA, APO); sL = score(CL, LCO)
    sAa = score(CAll, APO); sAl = score(CAll, LCO)
    emit("   DR17 IoU  per-tele: APO %.3f LCO %.3f   all-600: APO %.3f LCO %.3f\n",
        sA.iou, sL.iou, sAa.iou, sAl.iou)
    emit("   DR17 recall per-tele: APO %.3f LCO %.3f   all-600: APO %.3f LCO %.3f\n",
        sA.recall, sL.recall, sAa.recall, sAl.recall)
    emit("   masked frac per-tele: APO %.2f%% LCO %.2f%%   all-600: APO %.2f%% LCO %.2f%%\n\n",
        sA.frac, sL.frac, sAa.frac, sAl.frac)
end
# does a line clear the cut at one telescope but not the other?
emit("Lines detected at one telescope only (UNION level, i.e. not a single fiber at the other):\n")
CAu = cmb(nflag_apo, 1); CLu = cmb(nflag_lco, 1)
for (nm, only) in (("APO-only", CAu .& .!CLu), ("LCO-only", CLu .& .!CAu))
    rr = filter(r -> length(r) >= 3, mask_runs(only))
    emit("   %s: %d runs of >=3 px\n", nm, length(rr))
end

# The knife-edge problem with all-600 + majority, and the principled fix.
emit("KNIFE-EDGE CHECK for all-600 + majority (MEASURED):\n")
emit("a line seen by EVERY fiber at one telescope and NO fiber at the other lands at\n")
emit("nflag=300 — exactly the majority threshold. Pixels with nflag in [270,330]: %d\n",
    count(p -> 270 <= nflag[p] <= 330, 1:NPIX))
emit("=> an all-600 majority rule decides telescope-specific lines by a coin flip. A rule\n")
emit("   whose threshold sits far from 300 does not. (INFERRED design hazard, from the\n")
emit("   MEASURED count above.)\n\n")

emit("PRINCIPLED per-telescope combinations (MEASURED):\n")
MAJ_A = cmb(nflag_apo, 150); MAJ_L = cmb(nflag_lco, 150)
for (nm, C) in (("OR  of per-tele majorities", MAJ_A .| MAJ_L),
    ("AND of per-tele majorities", MAJ_A .& MAJ_L),
    ("all-600 nflag>=150 (25%)", cmb(nflag, 150)),
    ("all-600 majority (>=300)", cmb(nflag, 300)))
    sa = score(C, APO); sl = score(C, LCO)
    emit("  %-28s %5d px %3d lines | IoU APO %.3f LCO %.3f | frac APO %.2f%% LCO %.2f%%\n",
        nm, count(C), nlines(C), sa.iou, sl.iou, sa.frac, sl.frac)
end
emit("  (OR = 'bright if robustly detected at EITHER site' — keeps telescope-specific lines\n")
emit("   without a knife-edge; AND = 'only if both sites agree' — drops them.)\n")

emit("\nDR17's OWN masks combined the same way, for reference:\n")
for m in (1, 3, 300, 600)
    C = cmb(dr_nflag, m)
    emit("   DR17 nflag>=%-4d : %6d px, %4d lines\n", m, count(C), nlines(C))
end

# ------------------------------------------------------- per-fiber differences
emit("\n## 5. Per-fiber difference from the combined mask (for the validation figure)\n\n")
function perfiber_diff(C)
    rows = NamedTuple[]
    for f in 1:600
        sub = S[:, f]
        over = M[:, f] .& .!C .& sub     # per-fiber masked, combined did not = OVER-masked
        under = C .& .!M[:, f] .& sub    # combined masks, per-fiber did not = UNDER-masked
        push!(rows, (f=f, over=count(over), under=count(under),
            net=count(over) - count(under), tot=count(sub),
            pf=count(M[:, f]), cb=count(C .& sub)))
    end
    return rows
end
h5open(joinpath(RESDIR, "e5_bright_combined.h5"), "w") do fh
    fh["nflag"] = nflag; fh["nvalid"] = nvalid
    fh["nflag_apo"] = nflag_apo; fh["nflag_lco"] = nflag_lco
    fh["dr17_nflag"] = dr_nflag
    for (nm, v) in (("union", cmb(nflag, 1)), ("drop12", cmb(nflag, 3)),
        ("majority", cmb(nflag, 300)), ("q25", cmb(nflag, 150)),
        # RECOMMENDED: bright if robustly (>=50% of that site's fibers) detected at EITHER
        # telescope. Keeps telescope-specific lines without the nflag=300 knife edge.
        ("telemaj_or", cmb(nflag_apo, 150) .| cmb(nflag_lco, 150)),
        ("telemaj_and", cmb(nflag_apo, 150) .& cmb(nflag_lco, 150)),
        ("apo_majority", cmb(nflag_apo, 150)), ("lco_majority", cmb(nflag_lco, 150)),
        ("apo_union", cmb(nflag_apo, 1)), ("lco_union", cmb(nflag_lco, 1)),
        ("apo_drop12", cmb(nflag_apo, 3)), ("lco_drop12", cmb(nflag_lco, 3)))
        fh["mask_" * nm] = Vector{Bool}(v)
    end
end
emit("PLATEAU / BOUNDARY CHECK (the discipline that caught the pass-1 calibration error):\n")
emit("mean DR17 pixel IoU over the count threshold: ")
for m in (1, 3, 10, 60, 150, 300, 450, 540, 600)
    emit("%d:%.3f ", m, score(cmb(nflag, m), 1:600).iou)
end
emit("\n=> broad interior plateau over m=150-540; m=300 (majority) is INTERIOR to it, not a\n")
emit("   grid edge. Both extremes (union m=1, unanimous m=600) are measurably worse.\n\n")
emit("METRIC CAVEAT (important, and it cuts against my own recommendation):\n")
emit("DR17's reference masks are PER FIBER and therefore PER TELESCOPE. Scoring a mask that\n")
emit("unifies across telescopes against them is a BIASED referee: it structurally rewards\n")
emit("per-telescope masks. It is a fair referee for 'how many fibers must agree' (that\n")
emit("question is symmetric within a telescope), but the per-telescope-vs-all-600 comparison\n")
emit("in section 4 must be read knowing the metric leans toward per-telescope — and even so\n")
emit("it shows no meaningful gap. Pixel IoU is NOT gameable by dilation (it penalizes\n")
emit("over-masking); the LINE precision column is, so it is a diagnostic, not a target.\n\n")

for (nm, C) in (("all-600 MAJORITY (nflag>=300)", cmb(nflag, 300)),)
    rows = perfiber_diff(C)
    emit("Reference combined mask: %s (%d px)\n", nm, count(C))
    emit("per-fiber OVER-masked  (per-fiber flagged, combined did not): med=%d max=%d\n",
        Int(median([r.over for r in rows])), maximum(r -> r.over, rows))
    emit("per-fiber UNDER-masked (combined flagged, per-fiber did not): med=%d max=%d\n",
        Int(median([r.under for r in rows])), maximum(r -> r.under, rows))
    so = sort(rows, by=r -> -r.over)[1:8]
    su = sort(rows, by=r -> -r.under)[1:8]
    emit("\nWORST OVER-maskers (per-fiber flagged most that the combined mask does not):\n")
    for r in so
        emit("   f%03d (%s): over=%4d under=%4d  perfiber=%d combined=%d\n",
            r.f, r.f <= 300 ? "APO" : "LCO", r.over, r.under, r.pf, r.cb)
    end
    emit("WORST UNDER-maskers (combined flags most that this fiber did not):\n")
    for r in su
        emit("   f%03d (%s): over=%4d under=%4d  perfiber=%d combined=%d\n",
            r.f, r.f <= 300 ? "APO" : "LCO", r.over, r.under, r.pf, r.cb)
    end
    open(joinpath(RESDIR, "perfiber_diff.txt"), "w") do fo
        println(fo, "adjfib tele over under net perfiber combined valid")
        for r in rows
            @printf(fo, "%d %s %d %d %d %d %d %d\n", r.f, r.f <= 300 ? "apo" : "lco",
                r.over, r.under, r.net, r.pf, r.cb, r.tot)
        end
    end
end
close(io)
println("\nreport -> ", joinpath(RESDIR, "combine_report.txt"))
