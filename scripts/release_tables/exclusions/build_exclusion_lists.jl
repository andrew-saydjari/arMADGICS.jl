# build_exclusion_lists.jl -- build the public exclusion CSVs for the pass-1 DR21 data share.
#
# CHAIN POSITION: runs after allExposure identity (needs identity.h5) and BEFORE the qa_bitmask
#   stage of the allExposure builder, which reads the three CSVs it writes.
# Writes ONLY into --out ; every other path is read-only.
#   bad_nights.csv   -- (tele, mjd) exclusion aids, union of 4 sources
#   bad_fibers.csv   -- chronic / fiber_night / fiber_mjd_range granularities
#   rv_railed.csv    -- full RV-grid-railed population (no S/N cut)
#
# INPUTS (all read-only):
#   --wavecal-dir   wavecal census folder holding bad_wavesoln_census.csv, serious_bits.csv,
#                   bad_fibers_census.csv (0.4.0: work/2026_09_11/wavecal_compare)
#   --trace-csv     trace_night_jump_outliers.csv (0.4.0: work/2026_09_13/validation_anomaly_hunt)
#   --identity      identity.h5 (allExposure step 1)
#   --workup        arMADGICS workup directory (RV_pixoff_final, RV_pixoff_disc_final, snr, RV_flag,
#                   ingestBit)
# OUTPUTS: --out    directory receiving the three CSVs
# The upstream census CSVs are research products of separate analyses and are inputs here; the
# SRC_* strings below are copied verbatim into the `source` column of the released CSVs, so they
# name the 0.4.0 provenance and must not be edited if the 0.4.0 lists are to be reproduced.
# SOURCE: work/2026_09_16/public_exclusions/build_public_lists.jl (0.4.0). Logic unchanged; the
#         hard-coded paths became the arguments above.
#
# Usage: julia --project=scripts/release_tables build_exclusion_lists.jl --wavecal-dir D \
#            --trace-csv F --identity identity.h5 --workup W --out OUTDIR
using CSV, DataFrames, HDF5, Printf, Statistics, ArgParse

function parse_cli()
    s = ArgParseSettings(description = "Build bad_nights.csv, bad_fibers.csv, rv_railed.csv")
    @add_arg_table! s begin
        "--wavecal-dir"; required = true
        "--trace-csv"; required = true
        "--identity"; required = true
        "--workup"; required = true
        "--out"; required = true
    end
    return parse_args(s)
end
const ARGV = parse_cli()
const OUT = ARGV["out"]
const WCDIR = ARGV["wavecal-dir"]
const TRACECSV = ARGV["trace-csv"]
const IDENTITY = ARGV["identity"]
const WK = ARGV["workup"]
mkpath(OUT)

const SRC_CENSUS = "2026_09_11/wavecal_compare/bad_wavesoln_census.csv"
const SRC_SERIOUS = "2026_09_11/wavecal_compare/serious_bits.csv"
const SRC_ARCLAMP = "2026_09_14/arclamp_wavecal/check_offdetector.jl"
const SRC_TRACE = "2026_09_13/validation_anomaly_hunt/trace_night_jump_outliers.csv"
const SRC_FIBCEN = "2026_09_11/wavecal_compare/bad_fibers_census.csv"
const SRC_CHI2 = "2026_09_13/validation_anomaly_hunt/ANOMALY_HUNT.md sec 2.4"
const SRC_TRACEBLK = "2026_09_13/validation_anomaly_hunt/ANOMALY_HUNT.md sec 2.1 (T1)"
const SRC_RAIL = "2026_09_13/validation_anomaly_hunt/s17_rvgrid_rail.jl"

say(s) = (println(s); flush(stdout))

# ============================================================ 1. bad_nights.csv
# accumulate per (tele,mjd): ordered class list, reason list, source list
mutable struct NightRec
    classes::Vector{String}
    reasons::Vector{String}
    sources::Vector{String}
end
nights = Dict{Tuple{String, Int}, NightRec}()
function addnight!(tele, mjd, class, reason, source)
    r = get!(nights, (tele, mjd), NightRec(String[], String[], String[]))
    for c in split(class, ";")
        c in r.classes || push!(r.classes, String(c))
    end
    reason in r.reasons || push!(r.reasons, reason)
    source in r.sources || push!(r.sources, source)
    return
end

# -- source A: the 84 RV-flag tele-nights of the wavecal census
cen = CSV.read(joinpath(WCDIR, "bad_wavesoln_census.csv"), DataFrame; types = Dict(:tele => String))
say("census rows: $(nrow(cen))")
nA = 0
for r in eachrow(cen)
    r.primary_class == "UNARBITRATED" && continue   # no verdict possible (weather); not an RV flag
    r.primary_class == "NO_NIGHTAVE" && continue
    global nA += 1
    fx(v, u) = (v === missing || (v isa AbstractFloat && isnan(v))) ? "n/a" : string(v) * u
    ms(v) = v === missing ? "n/a" : string(v)
    reason = "wavecal census: primary=$(r.primary_class) severity=$(ms(r.severity)) " *
             "bwt=$(r.bwt) chips=$(ms(r.chips_impacted)); " *
             "non-chronic bad fibers R/G/B=$(r.n_fib_bad_R)/$(r.n_fib_bad_G)/$(r.n_fib_bad_B); " *
             "max ambiguous fibers/chip=$(r.n_fib_ambig_max); " *
             "max |D|>=2 junk fibers/chip=$(r.n_fib_junk_max); " *
             "sky exposures=$(r.nexp_sky) (bad-well $(r.nexp_badwell)); " *
             "max within-exposure residual MAD=$(fx(r.mad_within_exp_max, " A")); " *
             "worst chip median offset=$(fx(r.med_off_ang_worst, " A")); " *
             "min arbitration coverage=$(r.arb_cov_min) fibers"
    # list the census primary class first, then any secondary classes it carries
    allcls = String[String(r.primary_class)]
    for c in split(String(r.classes), ";")
        String(c) in allcls || push!(allcls, String(c))
    end
    addnight!(r.tele, r.mjd, join(allcls, ";"), reason, SRC_CENSUS)
end
say("source A (wavecal census RV-flag nights): $nA")

# -- source B: the 23 serious_bits nights
ser = CSV.read(joinpath(WCDIR, "serious_bits.csv"), DataFrame; types = Dict(:tele => String))
nB = 0
for r in eachrow(ser)
    r.serious === true || continue
    global nB += 1
    reason = if r.serious_reason == "nofile+no_calframes"
        "data-quality triage: zero same-night calibration frames of any kind (no darks, " *
        "arclamps, quartzflats or domeflats in the almanac) AND fluxing NOFILE; this is the " *
        "57800 pattern -- the night's $(r.n_obj) object exposures are fake/placeholder " *
        "commissioning data and any science derived from them is meaningless, which is why " *
        "the night is segregated rather than investigated"
    else
        "data-quality triage: $(r.n_obj_flagged) of $(r.n_obj) object exposures carry almanac " *
        "exposure_flags (bit 0 predicted_bad, bit 1 engineering, bit 2 notrun); " *
        "$(r.n_calframes) same-night calibration frames; fluxing NOFILE=$(r.nofile)"
    end
    addnight!(r.tele, r.mjd, "COMMISSIONING_FLAGGED", reason, SRC_SERIOUS)
end
say("source B (serious_bits nights): $nB")

# -- source C: apo 58014, found after the census closed
addnight!("apo", 58014, "WAVECAL_BROKEN",
    "arclamp verification: mapping the night's own arclamp peaks through its adopted pass-1 " *
    "wavelength solution puts 6,331 of 33,541 peaks (18.9%) outside the physical detector " *
    "(14000-18000 A), across 26 fibers; three control nights (apo 58383, apo 57708, " *
    "apo 59616) give 0.00%. The solution is unphysical, not merely imprecise, and nothing " *
    "in pass 1 flags it. Found after the wavecal census closed",
    SRC_ARCLAMP)
say("source C (apo 58014): 1")

# -- source D: trace-jump outlier tele-nights
tr = CSV.read(TRACECSV, DataFrame; types = Dict(:tele => String, :chip => String))
say("trace outlier tele-chip-nights: $(nrow(tr))")
nD = 0
for g in groupby(tr, [:tele, :mjd])
    global nD += 1
    tele = g.tele[1]
    mjd = g.mjd[1]
    chips = join(sort(unique(String.(g.chip))), "")
    r2(x) = round(x, digits = 2)
    r3(x) = round(x, digits = 3)
    detail = "night-to-night trace jump vs previous night $(g.prev_mjd[1]) (gap $(g.gap[1]) d): " *
             "chips $chips; median y-center jump $(r2(minimum(g.med_jump)))-" *
             "$(r2(maximum(g.med_jump))) px (max $(r2(maximum(g.max_jump))) px); up to " *
             "$(maximum(g.n_fib_jump_gt0p5)) of 300 fibers jump >0.5 px; median " *
             "neighbor-spacing jump $(r3(median(g.med_spacejump))) px " *
             "(max $(r2(maximum(g.max_spacejump))) px)"
    extra = if tele == "apo" && mjd == 57969
        ". SPECIAL: this is a total trace-solution failure for the night -- 299 of 300 fibers " *
        "move together on B+G and the trace Gaussian width also jumps (median deviation 0.50, " *
        "up to 6x), so the extraction geometry itself is wrong, not just offset"
    elseif tele == "lco" && mjd == 60702
        ". SPECIAL: this is a rigid whole-detector shift of 3.3-3.6 px on all 3 chips and all " *
        "300 fibers with neighbor spacing unchanged (median dspacing 0.01 px) -- a real " *
        "instrument/detector move rather than a fit failure; it is the single largest jump in " *
        "the corpus and only partially recovers at 60705 (1.7 px)"
    else
        ""
    end
    addnight!(tele, mjd, "TRACE_JUMP", detail * extra, SRC_TRACE)
end
say("source D (trace-jump tele-nights): $nD")

# -- assemble
ks = sort(collect(keys(nights)); by = k -> (k[1], k[2]))
bad_nights = DataFrame(
    tele = [k[1] for k in ks],
    mjd = [k[2] for k in ks],
    class = [join(nights[k].classes, ";") for k in ks],
    reason = [join(nights[k].reasons, " | ") for k in ks],
    source = [join(nights[k].sources, ";") for k in ks])
CSV.write(joinpath(OUT, "bad_nights.csv"), bad_nights)
say("bad_nights.csv rows: $(nrow(bad_nights))")

# per-class report (a row can carry several classes)
clscount = Dict{String, Int}()
for k in ks, c in nights[k].classes
    clscount[c] = get(clscount, c, 0) + 1
end
say("--- bad_nights class counts (tele-nights carrying the class; rows may carry several)")
for (c, v) in sort(collect(clscount); by = x -> -x[2])
    @printf("  %-22s %4d\n", c, v)
end
say("--- bad_nights primary (first-listed) class counts")
pc = Dict{String, Int}()
for k in ks
    pc[nights[k].classes[1]] = get(pc, nights[k].classes[1], 0) + 1
end
for (c, v) in sort(collect(pc); by = x -> -x[2])
    @printf("  %-22s %4d\n", c, v)
end

# ============================================================ 2. bad_fibers.csv
fibc = CSV.read(joinpath(WCDIR, "bad_fibers_census.csv"), DataFrame;
    types = Dict(:tele => String, :chip => String, :bwt => String))
say("\nbad_fibers_census rows: $(nrow(fibc))")

# chronic set, reproducing classify_census.jl (>=50 appearances at nlines>=3)
const CHRONIC_N = 50
appear = Dict{Tuple{String, Int}, Int}()
slipn = Dict{Tuple{String, Int}, Int}()
for r in eachrow(fibc)
    (r.nlines === missing || r.nlines < 3) && continue
    (r.D === missing || r.fib_med_ang === missing) && continue
    k = (r.tele, r.fiber)
    appear[k] = get(appear, k, 0) + 1
    if r.bwt == "fpi" && abs(r.D) == 1 && r.margin !== missing && r.margin <= 0.2
        slipn[k] = get(slipn, k, 0) + 1
    end
end
chronic = sort([k for (k, v) in appear if v >= CHRONIC_N]; by = k -> (k[1], -appear[k]))
say("chronic fibers (>=$CHRONIC_N appearances): $(length(chronic)) -> $chronic")
@assert length(chronic) == 11
chronicset = Set(chronic)

fib_rows = DataFrame(tele = String[], fiber = Int[], scope = String[], mjd = Union{Int, Missing}[],
    mjd_start = Union{Int, Missing}[], mjd_end = Union{Int, Missing}[],
    class = String[], reason = String[], source = String[])

# -- scope = chronic : 11 wavecal-chronic fibers
for k in chronic
    push!(fib_rows, (k[1], k[2], "chronic", missing, missing, missing, "WAVECAL_CHRONIC",
        "chronically wavelength-unreliable: $(appear[k]) bad-night appearances in the wavecal " *
        "census (nlines>=3), of which $(get(slipn, k, 0)) are confident +-1 FPI-mode slips; " *
        "typical signature |D|>=2 with 7-33 A offsets or multi-thousand-A garbage -- a dead or " *
        "damaged fiber. These 11 fibers supply 4,983 of the 13,999 census detail rows and are " *
        "excluded from all per-night counts so they do not spuriously flag ~500 nights. Treat " *
        "this (tele, fiber) wavelength solution as untrusted corpus-wide",
        SRC_CENSUS))
end
# -- scope = chronic : the 2 chi2-pathological APO fibers
push!(fib_rows, ("apo", 202, "chronic", missing, missing, missing, "CHI2_PATHOLOGY",
    "median chi2/dof = 37 with 94.6% of rows above chi2/dof = 2; 64% of rows have negative " *
    "skyscale1 (the #1 contributor to the negative-skyscale1 census, 35,502 rows); RV " *
    "zero-point -84 m/s, the largest in the corpus -- while relthrpt (1.00) and S/N (70) look " *
    "entirely normal, so no existing pass-1 flag fires. Uniformly bad, not intermittent " *
    "(chi2 bimodality index 0.0). Cause unresolved", SRC_CHI2))
push!(fib_rows, ("apo", 201, "chronic", missing, missing, missing, "CHI2_PATHOLOGY",
    "median chi2/dof = 18.6 with 91.9% of rows above chi2/dof = 2, while relthrpt and S/N look " *
    "normal, so no existing pass-1 flag fires. Uniformly bad, not intermittent (chi2 " *
    "bimodality index 0.0). Cause unresolved", SRC_CHI2))

# -- scope = fiber_night : the 3,139 flagged (tele, mjd, fiber) pairs
# gate reproduces classify_census.jl: nlines>=3, non-chronic, and either
#   a confident FPI mode slip (bwt=fpi, |D|==1, margin<=0.2) or
#   a sky per-fiber offset (bwt=sky, |fib_med_ang| > 1 A)
fn = Dict{Tuple{String, Int, Int}, Vector{Tuple{String, String, Float64, Float64, Int, Float64}}}()
for r in eachrow(fibc)
    (r.nlines === missing || r.nlines < 3) && continue
    (r.D === missing || r.fib_med_ang === missing) && continue
    (r.tele, r.fiber) in chronicset && continue
    isslip = r.bwt == "fpi" && abs(r.D) == 1 && r.margin !== missing && r.margin <= 0.2
    isoff = r.bwt == "sky" && abs(r.fib_med_ang) > 1.0
    (isslip || isoff) || continue
    k = (r.tele, r.mjd, r.fiber)
    push!(get!(fn, k, Vector{Tuple{String, String, Float64, Float64, Int, Float64}}()),
        (isslip ? "FPI_MODE_SLIP" : "SKY_FIBER_OFFSET", r.chip, Float64(r.D),
            r.margin === missing ? NaN : Float64(r.margin), r.nlines, Float64(r.fib_med_ang)))
end
say("fiber_night distinct (tele,mjd,fiber): $(length(fn))")
@assert length(fn) == 3139
for k in sort(collect(keys(fn)); by = x -> (x[1], x[2], x[3]))
    v = fn[k]
    cls = String[]
    for e in v
        e[1] in cls || push!(cls, e[1])
    end
    det = join([@sprintf("%s: D=%+g margin=%s nlines=%d fiber median offset=%.2f A",
                    e[2], e[3], isnan(e[4]) ? "n/a" : @sprintf("%.3f", e[4]), e[5], e[6])
                for e in sort(v; by = x -> x[2])], "; ")
    push!(fib_rows, (k[1], k[3], "fiber_night", k[2], missing, missing, join(cls, ";"),
        "per-fiber wavecal flag on this night -- " * det, SRC_FIBCEN))
end

# -- scope = fiber_mjd_range : the APO chip-B trace block
for f in 151:165
    push!(fib_rows, ("apo", f, "fiber_mjd_range", missing, 57969, 58015, "TRACE_BLOCK",
        "APO chip B trace block: fibers 151-165 sit 6.3-6.5 px off the corpus-median trace " *
        "y-center with neighbor spacing off by 5.7 px for ~24 consecutive nights, while the " *
        "other 285 fibers are normal. Fifteen fibers is exactly half a 30-fiber MTP bundle, so " *
        "this reads as a block-level trace mis-assignment rather than a drift. Only MJD 58015 " *
        "appears in the bad-wavecal census; 57969-58014 carry no pass-1 flag of any kind",
        SRC_TRACEBLK))
end

CSV.write(joinpath(OUT, "bad_fibers.csv"), fib_rows)
say("bad_fibers.csv rows: $(nrow(fib_rows))")
for s in ("chronic", "fiber_night", "fiber_mjd_range")
    @printf("  scope=%-16s %5d\n", s, count(==(s), fib_rows.scope))
end

# ============================================================ 3. rv_railed.csv
say("\nreading identity + workup for the RV rail population ...")
id = h5open(IDENTITY, "r")
tele = read(id["tele"])
mjd = read(id["mjd"])
expnum = read(id["expnum"])
fib = read(id["adjfiberindx"])
sdss_id = read(id["sdss_id"])
close(id)
po = h5read(joinpath(WK, "arMADGICS_out_RV_pixoff_final.h5"), "RV_pixoff_final")
pd = h5read(joinpath(WK, "arMADGICS_out_RV_pixoff_disc_final.h5"), "RV_pixoff_disc_final")
snr = h5read(joinpath(WK, "arMADGICS_out_snr.h5"), "snr")
RV_flag = h5read(joinpath(WK, "arMADGICS_out_RV_flag.h5"), "RV_flag")
ingestBit = h5read(joinpath(WK, "arMADGICS_out_ingestBit.h5"), "ingestBit")

# criterion reused verbatim from s17_rvgrid_rail.jl
ingest_ok = (ingestBit .& (2^0 | 2^1 | 2^2 | 2^6)) .== 0
ok = ingest_ok .& (RV_flag .== 0) .& .!isnan.(po)
say("usable rows: $(count(ok))")
lo_d = minimum(filter(!isnan, pd[ok]))
hi_d = maximum(filter(!isnan, pd[ok]))
say(@sprintf("coarse-grid (disc) bounds: %.5f .. %.5f pix", lo_d, hi_d))
@assert isapprox(abs(lo_d), 80.6; atol = 1e-6) && isapprox(hi_d, 80.6; atol = 1e-6)
prail = (pd .<= lo_d + 1e-9) .| (pd .>= hi_d - 1e-9)
m = ok .& prail
say("railed rows (full population, no S/N cut): $(count(m))")
say("  of those snr>20: $(count(m .& (snr .> 20)))")
# honest re-check: does any railed row carry a nonzero RV_flag if we drop that gate?
m_noflag = ingest_ok .& .!isnan.(po) .& prail
say("railed rows dropping the RV_flag==0 gate: $(count(m_noflag)) " *
    "(nonzero RV_flag among them: $(count(m_noflag .& (RV_flag .!= 0))))")

const c_kms = 299792.458
const delLog = 6e-6
z2v(z) = c_kms * ((1 + z)^2 - 1) / ((1 + z)^2 + 1)
sel = findall(m)
sel = sel[sortperm([(tele[i], mjd[i], expnum[i], fib[i]) for i in sel])]
rail = DataFrame(
    tele = [tele[i] == 1 ? "apo" : "lco" for i in sel],
    mjd = [string(mjd[i]) for i in sel],
    expnum = [string(expnum[i]) for i in sel],
    adjfiberindx = [string(fib[i]) for i in sel],
    sdss_id = [string(sdss_id[i]) for i in sel],
    snr = snr[sel],
    pixoff = po[sel],
    RV_kms = [z2v(10.0^(po[i] * delLog) - 1) for i in sel],
    RV_flag = [Int(RV_flag[i]) for i in sel])
CSV.write(joinpath(OUT, "rv_railed.csv"), rail)
say("rv_railed.csv rows: $(nrow(rail))  (all RV_flag==0: $(all(rail.RV_flag .== 0)))")
say("DONE")
