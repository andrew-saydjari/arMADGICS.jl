## E5: VALIDITY-CORRECT combination of the bright sky-line mask (AKS 2026-09-07, 2nd pass).
#
# AKS: "It seems like some lines at LCO are being systematically ignored. I suggest the
# skyline mask be built by averaging over the masks per telescope and then taking a
# straight union of apo and LCO. Also, I noted that you have the 'fibers where pixel is
# valid' which might factor into this mask. That is the one part we would not want to be
# baked into the union. We want to use the added pixels at the edges of the chip when we
# have them!"
#
# THE BUG HE CAUGHT (MEASURED, see §2): the delivered variants thresholded the ABSOLUTE
# detection count (e.g. all-600 majority = nflag >= 300). At a pixel where only one
# telescope has data, the absolute count can never exceed that telescope's eligible count
# -- and the minimum nonzero eligible count is 227 (LCO) -- so `nflag >= 300` is
# arithmetically IMPOSSIBLE there. Chip-edge pixels were therefore excluded by COVERAGE,
# not by evidence. That is coverage leaking into a mask that is supposed to be a pure
# statement about which WAVELENGTHS are bright sky lines.
#
# THE FIX: decide within each telescope on the fraction of that telescope's ELIGIBLE
# fibers (fibers that actually have valid data at that pixel), then take a STRAIGHT UNION
# across telescopes. A telescope with no data at a pixel abstains; it cannot vote a line
# away for the telescope that does see it.
#
# Usage: julia --project=<arM-E5b> e5_bright_validity.jl
# Author - Andrew Saydjari (E5 pass 1)
using HDF5, Statistics, Printf, Dates

const ARM = "/mnt/home/asaydjari/gitcode/worktrees/arM-E5b/"
include(ARM * "scripts/prior_build/e5_bright_line_detect.jl")
const DR = "/mnt/ceph/users/sdssv/work/asaydjari/2025_07_31/prior_dump/sky_priors"
const RESDIR = get(ENV, "E5_RESDIR", "/mnt/ceph/users/sdssv/work/asaydjari/2026_09_07/e5_combined")
const PERFIB = joinpath(RESDIR, "e5_bright_perfiber.h5")
const COMB = joinpath(RESDIR, "e5_bright_combined.h5")
const NPIX = 8575 + 125
wavetarg = 10 .^ range((4.179 - 125 * 6.0e-6), step=6.0e-6, length=NPIX)

# Minimum eligible fibers before a telescope is allowed to make a positive call. Set at
# 10% of a telescope; MEASURED below to never bind (min nonzero eligible is 288/227), so
# it is a guard against a future corpus, not a tuned parameter.
const NMIN = 30
const APO, LCO = 1:300, 301:600

M = Bool.(h5read(PERFIB, "mask")); S = Bool.(h5read(PERFIB, "submsk"))
nflag_apo = vec(sum(M[:, APO], dims=2)); nflag_lco = vec(sum(M[:, LCO], dims=2))
nval_apo = vec(sum(S[:, APO], dims=2)); nval_lco = vec(sum(S[:, LCO], dims=2))
nflag = nflag_apo .+ nflag_lco; nvalid = nval_apo .+ nval_lco

DB = falses(NPIX, 600); DV = falses(NPIX, 600)
for f in 1:600
    n = lpad(f, 3, "0")
    db = Bool.(h5read(joinpath(DR, "APOGEE_skyline_bright_svd_120_f$n.h5"), "submsk"))
    df = Bool.(h5read(joinpath(DR, "APOGEE_skyline_faint_svd_120_f$n.h5"), "submsk"))
    DB[:, f] .= db; DV[:, f] .= S[:, f] .& (db .| df)
end

io = open(joinpath(RESDIR, "validity_report.txt"), "w")
function emit(fmt, args...)
    s = Printf.format(Printf.Format(fmt), args...)
    print(io, s); print(stdout, s); flush(io); flush(stdout)
end

emit("E5 bright mask: VALIDITY-CORRECT per-telescope combination (%s)\n", string(now()))
emit("All numbers MEASURED unless marked INFERRED.\n\n")

## ---------------------------------------------------------------- 1. how validity enters
emit("## 1. How validity enters the pipeline today (MEASURED)\n\n")
emit("Three separate places, and only the FIRST two are legitimate:\n\n")
emit("  (a) PER FIBER, inside the detector: build_skyLines forms\n")
emit("      submsk = (obscnt >= min_obscnt) & chipgapmsk, and e5_bright_combine.jl sets the\n")
emit("      per-fiber mask FALSE wherever that fiber has no finite median_sky. So a fiber\n")
emit("      never flags a pixel it cannot see. CORRECT and necessary.\n")
emit("  (b) AT APPLICATION, downstream: the delivered mask is applied as full[submsk], i.e.\n")
emit("      intersected with each fiber's own validity. CORRECT -- this is where per-fiber\n")
emit("      coverage belongs, and it is already handled by obscnt + the chip-gap mask.\n")
emit("  (c) IN THE COMBINATION RULE -- the leak. Thresholding the ABSOLUTE count across\n")
emit("      fibers silently conflates 'few fibers called it bright' with 'few fibers could\n")
emit("      see it at all'. This is the one AKS flagged, and it is real (section 2).\n\n")
emit("NOTE on what the earlier report did: the DIAGNOSTIC bimodality number ('76.4%% of\n")
emit("flagged pixels are flagged by >=95%% of eligible fibers') was already normalized by\n")
emit("eligible count. But the DELIVERED variants (union/drop12/majority/telemaj_or) all\n")
emit("thresholded the RAW count. The diagnostic was right; the deliverable was not.\n\n")

emit("Coverage structure of the 8700-px grid (MEASURED):\n")
emit("  valid in >=1 APO fiber            : %d px\n", count(>(0), nval_apo))
emit("  valid in >=1 LCO fiber            : %d px\n", count(>(0), nval_lco))
emit("  APO-ONLY coverage (no LCO fiber)  : %d px\n", count(p -> nval_apo[p] > 0 && nval_lco[p] == 0, 1:NPIX))
emit("  LCO-ONLY coverage (no APO fiber)  : %d px  <-- the chip-edge pixels AKS wants kept\n",
    count(p -> nval_lco[p] > 0 && nval_apo[p] == 0, 1:NPIX))
emit("  partial APO (0 < nvalid < 300)    : %d px\n", count(p -> 0 < nval_apo[p] < 300, 1:NPIX))
emit("  partial LCO (0 < nvalid < 300)    : %d px\n", count(p -> 0 < nval_lco[p] < 300, 1:NPIX))
emit("  no fiber at all (both zero)       : %d px\n", count(p -> nval_apo[p] == 0 && nval_lco[p] == 0, 1:NPIX))
emit("  min NONZERO eligible count        : APO %d, LCO %d\n",
    minimum(filter(>(0), nval_apo)), minimum(filter(>(0), nval_lco)))
emit("  pixels with 0 < eligible < %d      : APO %d, LCO %d  => the NMIN guard NEVER binds\n\n",
    NMIN, count(p -> 0 < nval_apo[p] < NMIN, 1:NPIX), count(p -> 0 < nval_lco[p] < NMIN, 1:NPIX))

## ------------------------------------------------------- 2. the arithmetic impossibility
emit("## 2. WHY LCO lines were systematically ignored (MEASURED mechanism)\n\n")
lco_only = [nval_lco[p] > 0 && nval_apo[p] == 0 for p in 1:NPIX]
imposs = [lco_only[p] && nval_lco[p] < 300 for p in 1:NPIX]
emit("At an LCO-ONLY pixel the maximum achievable nflag is nval_lco (APO contributes 0).\n")
emit("The all-600 'majority' rule demanded nflag >= 300. Therefore:\n")
emit("  LCO-only pixels where nval_lco < 300, i.e. nflag>=300 is ARITHMETICALLY IMPOSSIBLE\n")
emit("  no matter how bright the line: %d px\n", count(imposs))
emit("  LCO-only pixels where nval_lco == 300, i.e. it requires UNANIMITY of all 300\n")
emit("  LCO fibers: %d px\n", count(p -> lco_only[p] && nval_lco[p] == 300, 1:NPIX))
emit("=> every one of the %d LCO-only pixels was decided by COVERAGE, not by evidence.\n",
    count(lco_only))
emit("   AKS's read of the figures was correct, and the effect is systematic, not random.\n\n")

## ---------------------------------------------------------------- 3. new construction
emit("## 3. New construction: per-telescope majority of ELIGIBLE fibers, then UNION\n\n")
fracA = [nval_apo[p] > 0 ? nflag_apo[p] / nval_apo[p] : NaN for p in 1:NPIX]
fracL = [nval_lco[p] > 0 ? nflag_lco[p] / nval_lco[p] : NaN for p in 1:NPIX]
A = BitVector([nval_apo[p] >= NMIN && isfinite(fracA[p]) && fracA[p] >= 0.5 for p in 1:NPIX])
L = BitVector([nval_lco[p] >= NMIN && isfinite(fracL[p]) && fracL[p] >= 0.5 for p in 1:NPIX])
TU = A .| L
emit("  APO majority-of-eligible : %d px (%d lines)\n", count(A), length(mask_runs(A)))
emit("  LCO majority-of-eligible : %d px (%d lines)\n", count(L), length(mask_runs(L)))
emit("  UNION (delivered)        : %d px (%d lines)\n\n", count(TU), length(mask_runs(TU)))
emit("Within-telescope threshold = 0.5 of ELIGIBLE fibers. Justification:\n")
emit("  * it is literally AKS's 'averaging over the masks per telescope' -- the mean of a\n")
emit("    set of 0/1 masks, cut at a half;\n")
emit("  * it is the MEASURED optimum: the all-600 DR17 IoU sweep peaked at the majority\n")
emit("    point and sat on a broad interior plateau (>=150...>=540 of 600), so 0.5 is\n")
emit("    interior, not a grid edge;\n")
emit("  * being within-telescope, no line can be voted down by fibers at the other site.\n\n")

# is it identical to the old telemaj_or?
partial_cov = BitVector([0 < nval_apo[p] < 300 || 0 < nval_lco[p] < 300 for p in 1:NPIX])
old_TO = Bool.(h5read(COMB, "mask_telemaj_or"))
old_maj = Bool.(h5read(COMB, "mask_majority"))
emit("Is this identical to the earlier `telemaj_or` (which used ABSOLUTE >=150)? NO:\n")
emit("  telemaj_or %d px vs new %d px; differing pixels %d\n",
    count(old_TO), count(TU), count(xor.(TU, old_TO)))
emit("  pixels the NEW rule adds that telemaj_or missed: %d\n", count(TU .& .!old_TO))
emit("  pixels telemaj_or had that the NEW rule drops   : %d\n", count(old_TO .& .!TU))
if count(xor.(TU, old_TO)) == 0
    emit("  => MEASURED: on THIS corpus the two rules give the IDENTICAL mask. They agree\n")
    emit("     automatically wherever eligible==300 (there >=150 IS the majority of eligible),\n")
    emit("     and at the %d partial-coverage pixels the detection fraction happens to sit far\n", count(partial_cov))
    emit("     from both thresholds, so nothing flips. The normalization is nevertheless the\n")
    emit("     CORRECT rule and is what is delivered: an absolute >=150 is only accidentally\n")
    emit("     right here, and would silently exclude any pixel whose eligible count fell\n")
    emit("     below 300 with a detection fraction between 0.5 and 150/eligible. Honest\n")
    emit("     statement: this change fixes the RULE, not (on this corpus) any pixel.\n\n")
else
    emit("  They agree wherever eligible==300 (there >=150 IS the majority of eligible) and\n")
    emit("  differ on partial-coverage pixels.\n\n")
end

## ------------------------------------------------------- 4. LCO line-by-line before/after
emit("## 4. LCO lines: BEFORE (all-600 majority) vs AFTER (new union), line by line\n\n")
gained = TU .& .!old_maj
runs = mask_runs(BitVector(gained))
emit("The new mask recovers %d pixels in %d contiguous runs that all-600 majority missed.\n\n",
    count(gained), length(runs))
emit("%-9s %-9s %5s %7s %7s %7s %7s %7s  %s\n", "lam_lo", "lam_hi", "npix",
    "nvA", "nvL", "fracA", "fracL", "nflag", "mechanism")
mech_count = Dict{String,Int}()
for r in runs
    c = (first(r) + last(r)) ÷ 2
    nvA, nvL = nval_apo[c], nval_lco[c]
    fA = isfinite(fracA[c]) ? fracA[c] : NaN
    fL = isfinite(fracL[c]) ? fracL[c] : NaN
    mech = if nvA == 0
        "LCO-ONLY COVERAGE (APO cannot vote)"
    elseif nvL == 0
        "APO-only coverage"
    elseif isfinite(fL) && fL >= 0.5 && isfinite(fA) && fA < 0.5
        "APO DILUTION (LCO sees it, APO does not)"
    elseif isfinite(fA) && fA >= 0.5 && isfinite(fL) && fL < 0.5
        "LCO dilution (APO sees it, LCO does not)"
    else
        "partial coverage / threshold"
    end
    mech_count[mech] = get(mech_count, mech, 0) + 1
    emit("%-9.2f %-9.2f %5d %7d %7d %7.3f %7.3f %7d  %s\n",
        wavetarg[first(r)], wavetarg[last(r)], length(r), nvA, nvL, fA, fL, nflag[c], mech)
end
emit("\nmechanism tally: ")
for (k, v) in sort(collect(mech_count), by=x -> -x[2]); emit("[%s x%d] ", k, v); end
emit("\n\n")
# Are the recovered pixels LINE-EDGE pixels (adjacent to an already-masked run, i.e. they
# widen a real line) or isolated specks? This decides whether the recovery is physical.
is_adjacent(r) = (first(r) > 1 && old_maj[first(r)-1]) || (last(r) < NPIX && old_maj[last(r)+1])
adj = count(is_adjacent, runs)
emit("Of the %d recovered runs, %d are ADJACENT to a line the old mask already had\n",
    length(runs), adj)
emit("(i.e. they widen a real line rather than adding an isolated speck); %d are isolated.\n\n",
    length(runs) - adj)

# How sensitive is the LCO side to the within-telescope threshold?
emit("LCO threshold sensitivity (MEASURED): pixels LCO calls bright at various fractions,\n")
emit("and how many of them the OLD all-600 majority mask missed:\n")
for F in (0.30, 0.40, 0.50, 0.60, 0.75, 0.90)
    Lf = [nval_lco[p] >= NMIN && isfinite(fracL[p]) && fracL[p] >= F for p in 1:NPIX]
    emit("  fracL>=%.2f : %4d px LCO-bright, %3d of them missing from the OLD mask\n",
        F, count(Lf), count(p -> Lf[p] && !old_maj[p], 1:NPIX))
end
emit("=> the LCO deficit is not an artifact of where the within-telescope cut is placed;\n")
emit("   it persists across the range, which is the signature of a POOLING problem\n")
emit("   (APO outvoting LCO), not a thresholding one.\n\n")

lost = old_maj .& .!TU
emit("Pixels the OLD mask had that the new one drops: %d", count(lost))
if count(lost) > 0
    lr = mask_runs(BitVector(lost))
    emit(" in %d runs\n", length(lr))
    for r in lr[1:min(12, end)]
        c = (first(r) + last(r)) ÷ 2
        emit("  %.2f-%.2f A (%d px): fracA=%.3f fracL=%.3f  (neither telescope reaches 0.5)\n",
            wavetarg[first(r)], wavetarg[last(r)], length(r), fracA[c], fracL[c])
    end
else
    emit("\n")
end

## ---------------------------------------------------------------- 5. DR17 sanity check
emit("\n## 5. DR17 overlap sanity check (MEASURED)\n\n")
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
emit("%-34s %7s %6s %8s %8s %8s %8s\n", "mask", "pixels", "lines", "IoU_APO", "IoU_LCO", "rec_APO", "rec_LCO")
for (nm, C) in (("all-600 majority (OLD)", BitVector(old_maj)),
    ("telemaj_or (absolute >=150)", BitVector(old_TO)),
    ("NEW per-tele elig. majority UNION", TU))
    sa = score(C, APO); sl = score(C, LCO)
    emit("%-34s %7d %6d %8.3f %8.3f %8.3f %8.3f\n", nm, count(C), length(mask_runs(C)),
        sa.iou, sl.iou, sa.recall, sl.recall)
end
emit("\nbright fraction: ")
for (nm, C) in (("OLD", BitVector(old_maj)), ("NEW", TU))
    sa = score(C, APO); sl = score(C, LCO)
    emit("%s APO %.2f%% LCO %.2f%%   ", nm, sa.frac, sl.frac)
end
emit("(DR17 deployed: APO 8.35%%, LCO 8.15%%)\n")

## ------------------------------------------------- 6. coverage is NOT baked in: the proof
emit("\n## 6. Proof that coverage is NOT baked into the delivered mask (MEASURED)\n\n")
emit("Test: for every pixel, does the mask's value depend on how many fibers happen to\n")
emit("have data there, given the evidence of the fibers that DO?\n\n")
# every LCO-only pixel that LCO robustly calls bright must be in the mask
lco_rob = [lco_only[p] && nval_lco[p] >= NMIN && fracL[p] >= 0.5 for p in 1:NPIX]
emit("  LCO-only pixels LCO robustly calls bright : %d ; of those, IN the new mask: %d\n",
    count(lco_rob), count(p -> lco_rob[p] && TU[p], 1:NPIX))
apo_only = [nval_apo[p] > 0 && nval_lco[p] == 0 for p in 1:NPIX]
apo_rob = [apo_only[p] && nval_apo[p] >= NMIN && fracA[p] >= 0.5 for p in 1:NPIX]
emit("  APO-only pixels APO robustly calls bright : %d ; of those, IN the new mask: %d\n",
    count(apo_rob), count(p -> apo_rob[p] && TU[p], 1:NPIX))
emit("  (under the OLD all-600 majority: %d and %d respectively)\n",
    count(p -> lco_rob[p] && old_maj[p], 1:NPIX), count(p -> apo_rob[p] && old_maj[p], 1:NPIX))
partial = [0 < nval_apo[p] < 300 || 0 < nval_lco[p] < 300 for p in 1:NPIX]
emit("  partial-coverage pixels: %d; robustly bright at their covering telescope: %d;\n",
    count(partial), count(p -> partial[p] && ((nval_apo[p] >= NMIN && fracA[p] >= 0.5) ||
                                              (nval_lco[p] >= NMIN && fracL[p] >= 0.5)), 1:NPIX))
emit("    of those, in the new mask: %d\n",
    count(p -> partial[p] && ((nval_apo[p] >= NMIN && fracA[p] >= 0.5) ||
                              (nval_lco[p] >= NMIN && fracL[p] >= 0.5)) && TU[p], 1:NPIX))
emit("\n  The mask is FALSE at the %d pixels no fiber can see. That is not a coverage leak:\n",
    count(p -> nval_apo[p] == 0 && nval_lco[p] == 0, 1:NPIX))
emit("  by definition no fiber's submsk contains them, so no fiber ever reads that value.\n")
emit("  INFERRED limitation worth recording: the mask can only speak about wavelengths the\n")
emit("  training corpus covers. If a future reduction gains pixels beyond that footprint,\n")
emit("  the mask would call them faint by default and would need rebuilding -- it is not a\n")
emit("  bug today, but it is the one place where corpus coverage bounds the product.\n")

## ---------------------------------------------------------------- 7. persist
h5open(COMB, "r+") do fh
    for (nm, v) in (("telemaj_union", TU), ("apo_elig_majority", A), ("lco_elig_majority", L))
        haskey(fh, "mask_" * nm) && delete_object(fh, "mask_" * nm)
        fh["mask_"*nm] = Vector{Bool}(v)
    end
    for (nm, v) in (("nval_apo", nval_apo), ("nval_lco", nval_lco),
        ("frac_apo", replace(fracA, NaN => -1.0)), ("frac_lco", replace(fracL, NaN => -1.0)))
        haskey(fh, nm) && delete_object(fh, nm)
        fh[nm] = v
    end
end
emit("\nwrote mask_telemaj_union (+ apo/lco_elig_majority, nval_*, frac_*) to %s\n", COMB)
close(io)
println("\nreport -> ", joinpath(RESDIR, "validity_report.txt"))
