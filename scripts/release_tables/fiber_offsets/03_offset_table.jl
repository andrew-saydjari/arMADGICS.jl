# 03_offset_table.jl
#
# CHAIN POSITION: fiber offsets step 3 of 5 (see ../README.md).
# INPUTS:  ARGS[1] WORKDIR (reads offsets_raw.csv), ARGS[2] bad_fibers.csv (exclusion lists)
# OUTPUTS: WORKDIR/fiber_offsets.csv (review header; step 4 rewrites it), WORKDIR/fiber_offsets_bare.csv, WORKDIR/FIBER_OFFSETS_TABLE.md
# SOURCE:  work/2026_09_16/fiber_offsets/s3_table.jl (0.4.0). Logic unchanged; hard-coded paths became
#          ARGS (WORKDIR = the fiber-offset working directory shared by steps 1-5).
#
# s3_table.jl -- turn the raw solve into the shipped fiber_offsets.csv, decide which blocks are
# CALIBRATABLE, and cross-check the non-calibratable set against the release chronic bad-fiber
# list (the maintainer's explicit callout: DR17 had fibers that could not be calibrated and had to be
# dropped; are they already excluded by our cuts, or not?).
#
# CALIBRATABLE  :=  a solution exists  AND  n_stars_used >= 100  AND  err <= 15 m/s.
#   * 15 m/s is 3x the median bootstrap error (5.0 m/s) and sits in the clean gap between the
#     99th percentile of the error distribution (10.4 m/s) and the next value up (19.7 m/s).
#   * It also caps the noise the correction itself injects: (15/100)^2 = 2 % of the per-fiber
#     variance being removed, so a calibratable fiber's correction is >= 98 % signal.
#   * n_stars_used is the count of stars that actually INFORM the block -- stars seen on >= 2
#     blocks.  Stars confined to one block contribute exactly nothing (see s2 header).
using HDF5, CSV, DataFrames, Statistics, Printf, Dates

const DIR  = ARGS[1]
const BADF = ARGS[2]
const ERRMAX = 0.015     # km/s
const NSMIN  = 100
say(s) = (println(s); flush(stdout))

d = CSV.read(joinpath(DIR, "offsets_raw.csv"), DataFrame)
d.calibratable = isfinite.(d.offset_kms) .& (d.n_stars_used .>= NSMIN) .& (d.err_kms .<= ERRMAX) .& isfinite.(d.err_kms)

# reason code for every non-calibratable block
d.reason = map(eachrow(d)) do r
    r.calibratable                  ? "" :
    r.n_exp == 0                    ? "NO_DATA_AFTER_MASK" :
    !isfinite(r.offset_kms)         ? "DISCONNECTED" :
    r.n_stars_used < NSMIN          ? "TOO_FEW_LINKING_STARS" :
                                      "ERR_TOO_LARGE"
end

# offsets of non-calibratable blocks are NaN in the shipped table (the measured value is kept
# in offset_meas_kms so the diagnosis is not lost)
d.offset_meas_kms = copy(d.offset_kms)
d.err_kms_orig    = copy(d.err_kms)
d.offset_kms = [r.calibratable ? r.offset_kms : NaN for r in eachrow(d)]
d.err_kms    = [r.calibratable ? r.err_kms    : NaN for r in eachrow(d)]

say(@sprintf("calibratable blocks: %d / 1200  (non-calibratable %d)", count(d.calibratable), count(.!d.calibratable)))
for g in groupby(sort(d[.!d.calibratable, :], [:reason, :tele, :fiber]), :reason)
    say(@sprintf("  %-24s  %d blocks", g.reason[1], nrow(g)))
end

# ---------------------------------------------------------------- chronic cross-check
bf = CSV.read(BADF, DataFrame)
chronic = Set((String(r.tele), Int(r.fiber)) for r in eachrow(bf) if String(r.scope) == "chronic")
say("\nrelease chronic (tele, fiber) list: $(length(chronic)) fibers")

nc = d[.!d.calibratable, :]
ncf = unique([(r.tele, r.fiber) for r in eachrow(nc)])
inb  = [p for p in ncf if p in chronic]
notb = [p for p in ncf if !(p in chronic)]
say(@sprintf("non-calibratable fibers (any era): %d ; already on the chronic list: %d ; NOT excluded: %d",
    length(ncf), length(inb), length(notb)))
say("\n  --- NOT already excluded by the release cuts (the maintainer's callout) ---")
for (tl, fbr) in sort(notb)
    sub = nc[(nc.tele .== tl) .& (nc.fiber .== fbr), :]
    fid = 301 - fbr   # APOGEE FIBERID convention, = 301 - adjfiberindx_local
    for r in eachrow(sub)
        @printf("    %s fiber %3d (FIBERID %3d) era %d : measured %+8.1f +- %6.1f m/s, n_stars %4d, n_exp %6d  [%s]\n",
            tl, fbr, fid, r.era, 1000r.offset_meas_kms, 1000r.err_kms_orig, r.n_stars_used, r.n_exp, r.reason)
    end
end
# chronic fibers that never appear as non-calibratable (they have no data at all, so they cannot)
say("\n  --- chronic fibers, and whether they show up here ---")
for (tl, fbr) in sort(collect(chronic))
    sub = d[(d.tele .== tl) .& (d.fiber .== fbr), :]
    @printf("    %s fiber %3d : n_exp %d + %d (era 1 + 2) -> %s\n", tl, fbr,
        sub.n_exp[sub.era .== 1][1], sub.n_exp[sub.era .== 2][1],
        all(sub.n_exp .== 0) ? "no surviving rows (already fully excluded by qa_bitmask)" : "HAS SURVIVING ROWS")
end

# ---------------------------------------------------------------- ship
out = DataFrame(tele = d.tele, fiber = d.fiber, era = d.era,
                offset_kms = d.offset_kms, err_kms = d.err_kms,
                n_stars_used = d.n_stars_used, calibratable = d.calibratable,
                offset_meas_kms = d.offset_meas_kms, err_meas_kms = d.err_kms_orig,
                err_curv_kms = d.err_curv_kms, n_stars_any = d.n_stars_any,
                n_exp = d.n_exp, reason = d.reason)
sort!(out, [:tele, :era, :fiber])

hdr = """
# fiber_offsets.csv -- per (telescope, fiber, era) radial-velocity zero points, pass-1 DR21 (0.4.0)
# Written $(Dates.format(now(), "yyyy-mm-dd HH:MM")) by 2026_09_16/fiber_offsets/s2_solve.jl + s3_table.jl.
# STAGED FOR REVIEW.  Nothing here has been released.
#
# WHAT THESE ARE.  EMPIRICAL CALIBRATIONS, not physical corrections.  Each number is the
# constant velocity zero point of one (telescope, fiber, era) block, measured by the
# ubercalibration least-squares of Saydjari et al. 2025 (AJ 169, 167), Sec. 4: every exposure
# is modelled as y = mu_star + f_block, and mu and f are solved JOINTLY (stars seen on a single
# block contribute exactly zero, so the estimator is undiluted).  The block index absorbs the
# per-fiber term, the APO/LCO telescope term and the plate->FPS era step at once, on ONE common
# velocity frame.
#
# WE DO NOT KNOW WHY THE TELESCOPES DIFFER.  The implied cross-telescope term (LCO - APO,
# -253 m/s in the FPS era) has no established physical cause; it is open task #77.  Applying
# this table removes the offset; it does not explain it.  The same holds for the era steps
# (+24.6 / -61.7 m/s APO/LCO), which reproduce the independent 2026_09_13 measurements
# (+24.4 +- 1.3 / -60.9 +- 3.2 m/s) to better than 1 m/s but are equally unexplained.
#
# GAUGE.  Offsets have zero unweighted mean over the solved blocks, so v_rad - offset carries
# no net velocity shift relative to v_rad.  Only DIFFERENCES between blocks are measured; the
# absolute zero point is unconstrained by construction.
#
# TIME DEPENDENCE IS NOT INCLUDED.  One constant per block per era.  2026_09_13 measured a
# median 21.6 m/s epoch-to-epoch drift that this table does not model; it is a floor on how
# well the correction can work and is a stated limitation, not a claim of absence.
#
# COLUMNS
#   tele, fiber      telescope and fiber index 1..300 (= adjfiberindx - 300*(tele-1)).
#                    APOGEE FIBERID = 301 - fiber.
#   era              1 = plate, 2 = FPS
#   offset_kms       the zero point, km/s.  v_rad_corr = v_rad - offset_kms.  NaN where the
#                    block is not calibratable.
#   err_kms          1 sigma, from a 64-replicate bootstrap OVER STARS (carries the within-star
#                    and within-visit correlation).  NaN where not calibratable.
#   n_stars_used     stars that inform this block, i.e. seen on >= 2 blocks
#   calibratable     true / false.  false => offset_kms and err_kms are NaN
#   offset_meas_kms, err_meas_kms   the measured values, kept for non-calibratable blocks too
#                    so the diagnosis is auditable.  DO NOT APPLY THESE.
#   err_curv_kms     the analytic inverse-curvature error (paper Fig. 4, bottom right).  It
#                    assumes independent exposures and runs ~2.7x smaller than err_kms.
#   n_stars_any      stars with any exposure on this block
#   n_exp            exposures on this block passing the release default mask + science + S/N 15
#   reason           why a block is not calibratable
#
# SAMPLE.  Release default mask (qa_bitmask/ingestBit/RV_flag/bary_flag == 0, finite v_rad)
# AND category == science AND snr >= 15 AND an RV-stability cut on the star (robust scatter
# <= 1.0 km/s, max |v - median| <= 2.0 km/s), then 4-sigma residual clipping, 4 iterations.
# 21,032,672 exposures entered; 17,033,038 were in the final fit; 319,724 stars informed it.
"""
open(joinpath(DIR, "fiber_offsets.csv"), "w") do io
    print(io, hdr)
    CSV.write(io, out; append = true, writeheader = true)
end
chmod(joinpath(DIR, "fiber_offsets.csv"), 0o644)
say("\nwrote fiber_offsets.csv ($(nrow(out)) rows)")
