# 04_release_header.jl
#
# CHAIN POSITION: fiber offsets step 4 of 5 (see ../README.md).
# INPUTS:  ARGS[1] WORKDIR (reads fiber_offsets.csv from step 3)
# OUTPUTS: WORKDIR/fiber_offsets.csv (release header), WORKDIR/fiber_offsets_bare.csv (THE RELEASED fiber_offsets.csv: 0.4.0 shipped the bare file under that name), WORKDIR/FIBER_OFFSETS_TABLE.md (released)
# SOURCE:  work/2026_09_16/fiber_offsets/s9_csv_header.jl (0.4.0). Logic unchanged; hard-coded paths became
#          ARGS (WORKDIR = the fiber-offset working directory shared by steps 1-5).
#
# s9_csv_header.jl -- rewrite fiber_offsets.csv with a RELEASE-READY, self-contained header.
#
# The first header was written for internal review: it cited script paths, an internal task
# number and a dated working directory, and said "STAGED FOR REVIEW".  None of that means
# anything to an outside reader and the last would be false once the file ships.  This rewrites
# the preamble so the file explains itself, and leaves every data row byte-identical.
using CSV, DataFrames

const DIR = ARGS[1]
d = CSV.read(joinpath(DIR, "fiber_offsets.csv"), DataFrame; comment = "#")

hdr = """
# fiber_offsets.csv -- per-fiber radial-velocity zero points, pass-1 DR21 (0.4.0)
#
# WHAT THIS IS.  Each APOGEE fiber has its own radial-velocity zero point.  This file gives that
# zero point for every (telescope, fiber, era) combination, in km/s.  To use it:
#
#     v_rad_corr = v_rad - offset_kms          matched on (tele, fiber, era)
#
# HOW IT WAS DERIVED.  Stars observed on more than one fiber, telescope or era tie the fibers
# together; a joint least-squares fit solves for the star velocities and the fiber zero points at
# the same time, so a star seen on only one fiber contributes nothing and cannot dilute the
# result.  The method is the one in Saydjari et al. 2025, AJ 169, 167 (arXiv:2408.07126),
# Section 4.  The fit used science targets passing the release default mask at S/N >= 15, with
# binaries and variables held out by a cut on each star's own velocity scatter.
#
# ONE SET OF OFFSETS COVERS EVERYTHING.  Because the telescopes and eras share one fit, the
# APO/LCO difference (about 250 m/s) and the plate->FPS step (+25 m/s at APO, -61 m/s at LCO) are
# absorbed into these numbers.  Do not apply a separate telescope or era correction on top.
#
# THESE ARE EMPIRICAL CALIBRATIONS, NOT EXPLANATIONS.  The reason the two telescopes disagree by
# 250 m/s is not understood, and neither are the era steps.  Subtracting this table removes those
# differences; it does not account for them.  Treat the offsets as a measured instrument
# property, not as a physical model.  Only differences between fibers are measured -- the
# absolute zero point is not.  The offsets average to zero, so subtracting them shifts no
# velocity overall.
#
# NO TIME DEPENDENCE.  One constant per fiber per era.  Drift within an era is not modelled.
#
# NOT CALIBRATABLE.  A fiber is calibratable only if enough stars link it to other fibers and its
# zero point comes out to better than 15 m/s.  Where calibratable is false, offset_kms and
# err_kms are NaN and the fiber must not be corrected.  Four LCO fibers fail: 121, 240 and 270
# (FIBERID 180, 61 and 31), which are the same three fibers Saydjari et al. 2025 could not
# calibrate and dropped, and 182 (FIBERID 119), which simply has too little data in the FPS era.
# Thirteen further fibers appear here with no data at all: they are excluded across the whole survey by
# qa_bitmask, so no row of the release needs an offset for them.  See bad_fibers.csv.
#
# COLUMNS
#   tele, fiber      telescope, and fiber index 1..300 (= adjfiberindx - 300*(tele-1)).
#                    APOGEE FIBERID = 301 - fiber.
#   era              1 = plate, 2 = FPS
#   offset_kms       the zero point, km/s.  Subtract it from v_rad.  NaN if not calibratable.
#   err_kms          1 sigma, from a bootstrap over stars.  Median 5.0 m/s.  NaN if not
#                    calibratable.
#   n_stars_used     stars linking this fiber to at least one other, i.e. the stars that
#                    actually determined this number
#   calibratable     true / false.  false means offset_kms and err_kms are NaN
#   offset_meas_kms  the value the fit returned even where it is not trustworthy, kept so the
#   err_meas_kms     flag can be audited.  DO NOT APPLY THESE.
#   err_curv_kms     the formal least-squares uncertainty.  It assumes exposures are independent
#                    and runs about 2.7x smaller than err_kms; err_kms is the one to use
#   n_stars_any      stars with any exposure on this fiber
#   n_exp            exposures on this fiber that entered the fit
#   reason           why a fiber is not calibratable; empty when it is
"""
open(joinpath(DIR, "fiber_offsets.csv"), "w") do io
    print(io, hdr)
    CSV.write(io, d; append = true, writeheader = true)
end
chmod(joinpath(DIR, "fiber_offsets.csv"), 0o644)
println("rewrote fiber_offsets.csv")

# ---------------------------------------------------------------------------------------------
# The three CSVs already in public/ (bad_fibers, bad_nights, rv_railed) carry NO comment header;
# their prose lives in RECOMMENDED_CUTS.md.  A 52-line '#' preamble therefore breaks that
# directory's convention and breaks pandas.read_csv() unless the caller passes comment='#'.
# Rather than decide for the maintainer, emit the alternative packaging too, so the choice is one line:
#   A) fiber_offsets.csv          -- self-contained, header inline (written above)
#   B) fiber_offsets_bare.csv + FIBER_OFFSETS_TABLE.md -- matches the existing convention
CSV.write(joinpath(DIR, "fiber_offsets_bare.csv"), d)
chmod(joinpath(DIR, "fiber_offsets_bare.csv"), 0o644)

md = replace(hdr, r"^# ?"m => "")
open(joinpath(DIR, "FIBER_OFFSETS_TABLE.md"), "w") do io
    println(io, "# fiber_offsets.csv\n")
    println(io, "```")
    print(io, md)
    println(io, "```")
end
chmod(joinpath(DIR, "FIBER_OFFSETS_TABLE.md"), 0o644)
println("also wrote fiber_offsets_bare.csv + FIBER_OFFSETS_TABLE.md (option B packaging)")
