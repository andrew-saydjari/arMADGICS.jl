# schema.jl -- the ONE place that defines the allExposure table: every column's name, element
# type, fill value, the build stage that produces it, its on-disk layout, and any attributes.
# The writer (write.jl) and the tests iterate this table; nothing else hard-codes column lists.
#
# Layout classes reproduce the 0.4.0 file exactly (dtype, chunk shape, filter pipeline):
#   :gz      chunked, shuffle + gzip 1, chunk = 13045 rows for 8- and 4-byte types, 26090 for
#            2-byte, 52179 for 1-byte (what h5py chose for the 0.4.0 build)
#   :gz_nosh chunked as :gz but no shuffle (categorical uint8 codes and bary_flag in 0.4.0)
#   :gz13045 chunked 13045 rows, shuffle + gzip 1, regardless of element size (v_rad_corr_flag,
#            offtarget: added by later 0.4.0 patch steps that copied the host convention)
#   :contig  contiguous, no filters (glon, glat, qa_bitmask in 0.4.0)

"Where a column comes from. One stage = one former 0.4.0 script step (see README)."
@enum Stage begin
    S_IDENTITY      # workup row identity (full_list_info.h5)            0.4.0 s1 + s6
    S_ALMANAC       # almanac per-exposure / per-fiber metadata join     0.4.0 s6
    S_TIME          # harvested mid-exposure time                        0.4.0 s4 + s6
    S_BARY          # barycentric correction                             0.4.0 s6
    S_WORKUP        # workup scalar passthrough                          0.4.0 s6
    S_VELOCITY      # v_rel, v_rad, e_v_rad                              0.4.0 s6
    S_BARYFLAG      # bary_flag                                          0.4.0 s9
    S_GALACTIC      # glon, glat                                         0.4.0 s10
    S_QA            # qa_bitmask from the exclusion lists                0.4.0 s10
    S_FIBEROFFSET   # v_rad_corr, v_rad_offset_applied, v_rad_corr_flag  0.4.0 fiber_offsets/s10
    S_OFFTARGET     # offtarget                                          0.4.0 offtarget_sidecar s1-s4
end

struct ColumnSpec
    name::String
    T::DataType
    stage::Stage
    layout::Symbol
    attrs::Vector{Pair{String,Any}}   # dataset attributes (see AttrValue in h5util.jl)
end
ColumnSpec(name, T, stage, layout) = ColumnSpec(name, T, stage, layout, Pair{String,Any}[])

# ---------------------------------------------------------------------------- attribute texts
const BARY_FLAG_BITS = "1=NO_BARY (no barycentric correction available); " *
    "2=BAD_TIME (mid-exposure time implausible for the night; " *
    "LCO 57800 exp 1-8, correction evaluated at the wrong epoch, " *
    "wrong by up to the full +/-30 km/s); " *
    "4=ID_MISMATCH (almanac vs workup sdss_id disagree, RA/Dec " *
    "provenance uncertain, error bounded by plate radius ~0.8 km/s); " *
    "8=NO_RV (RV_pixoff_final not finite)"

const DESC_V_RAD_CORR = "Barycentric radial velocity with the per-fiber zero-point offset removed, km/s: " *
    "v_rad_corr = v_rad - v_rad_offset_applied. The offsets are in public/fiber_offsets.csv, one per " *
    "(telescope, fiber, era), matched on (tele, adjfiberindx - 300*(tele-1), era). They are an EMPIRICAL " *
    "CALIBRATION measured from stars observed on more than one fiber, telescope or era, by a joint " *
    "least-squares fit for the star velocities and the fiber zero points at the same time (method of " *
    "Saydjari et al. 2025, AJ 169, 167, Section 4). Because one fit covers both telescopes and both eras, " *
    "the APO/LCO difference (about 250 m/s) and the plate->FPS step are absorbed here too: do not apply a " *
    "separate telescope or era correction on top. The reason the telescopes disagree is not understood; " *
    "this removes the difference without explaining it. Only differences between fibers are measured, so " *
    "the absolute zero point is unchanged - the offsets average to zero. No time dependence is modelled. " *
    "NaN wherever v_rad_corr_flag is non-zero."
const DESC_V_RAD_OFFSET = "The (telescope, fiber, era) zero point subtracted from v_rad to give " *
    "v_rad_corr, km/s, taken from the offset_kms column of public/fiber_offsets.csv. Range -0.354 to " *
    "+0.275 km/s. NaN wherever v_rad_corr_flag is non-zero."
const DESC_V_RAD_CORR_FLAG = "Why v_rad_corr is or is not defined on this row. 0 = corrected. " *
    "1 = the row's fiber is not calibratable, so no offset exists and v_rad_corr is NaN; these are LCO " *
    "fibers 121, 240, 270 and 182 (FIBERID 180, 61, 31 and 119), where calibratable is false in " *
    "public/fiber_offsets.csv - the first three are the same fibers Saydjari et al. 2025 could not " *
    "calibrate and dropped. 2 = the row does not pass the release default mask " *
    "(qa_bitmask == ingestBit == RV_flag == bary_flag == 0 and finite v_rad; see RECOMMENDED_CUTS.md), " *
    "so no correction was attempted and v_rad_corr is NaN. The uncorrected v_rad is unaffected in every " *
    "case."
const DESC_OFFTARGET = "FPS off-target flag.  1 where the almanac fiber table for this row's " *
    "(telescope, mjd, config_id, fiber) has assigned == 1 and on_target == 0, i.e. the " *
    "robotic positioner was assigned this target but did not reach it, so the spectrum is " *
    "of whatever sky the fiber did reach while still carrying the assigned target's " *
    "identity (sdss_id, ra, dec).  Task 94's sdss_id 85604486 is the worked example: a " *
    "field star 9 arcmin away at S/N 96, clean on every other flag.  0 for all plate-era " *
    "rows -- the plate era has no assigned/on_target columns and cannot be flagged this " *
    "way.  Fiber index convention: fiber_index_almanac = 301 - (adjfiberindx - " *
    "300*(tele-1)).  Source: pass1_dr21/almanac/allobs_57618_61230.h5."

# Attribute value wrappers: the 0.4.0 file mixes fixed-length and variable-length UTF-8 strings
# (a consequence of which tool wrote which step); bit-identity of the HDF5 types needs both.
struct FixedStr; s::String; end
struct VlenStr;  s::String; end
struct VlenStrVec; v::Vector{String}; end

# Categorical columns get their `categories` attribute at write time (codes are data-driven).
const CATEGORICAL = ("category", "hole_type", "image_type")

"The 0.4.0 allExposure schema, in the order the 0.4.0 builder created the datasets."
const SCHEMA = ColumnSpec[
    ColumnSpec("tele",                Int8,    S_IDENTITY, :gz),
    ColumnSpec("mjd",                 Int32,   S_IDENTITY, :gz),
    ColumnSpec("expnum",              Int32,   S_IDENTITY, :gz),
    ColumnSpec("adjfiberindx",        Int16,   S_IDENTITY, :gz),
    ColumnSpec("fiber_index_almanac", Int16,   S_IDENTITY, :gz),
    ColumnSpec("sdss_id_workup",      Int64,   S_IDENTITY, :gz),
    ColumnSpec("sdss_id",             Int64,   S_ALMANAC,  :gz),
    ColumnSpec("catalogid",           Int64,   S_ALMANAC,  :gz),
    ColumnSpec("ra",                  Float64, S_ALMANAC,  :gz),
    ColumnSpec("dec",                 Float64, S_ALMANAC,  :gz),
    ColumnSpec("lambda_eff",          Float64, S_ALMANAC,  :gz),
    ColumnSpec("x_focal",             Float64, S_ALMANAC,  :gz),
    ColumnSpec("y_focal",             Float64, S_ALMANAC,  :gz),
    ColumnSpec("config_id",           Int32,   S_ALMANAC,  :gz),
    ColumnSpec("plate_id",            Int32,   S_ALMANAC,  :gz),
    ColumnSpec("design_id",           Int32,   S_ALMANAC,  :gz),
    ColumnSpec("field_id",            Int32,   S_ALMANAC,  :gz),
    ColumnSpec("cart_id",             Int32,   S_ALMANAC,  :gz),
    ColumnSpec("map_id",              Int32,   S_ALMANAC,  :gz),
    ColumnSpec("n_read",              Int32,   S_ALMANAC,  :gz),
    ColumnSpec("fps",                 Int32,   S_ALMANAC,  :gz),
    ColumnSpec("flagged_bad",         Int32,   S_ALMANAC,  :gz),
    ColumnSpec("chip_flags",          Int32,   S_ALMANAC,  :gz),
    ColumnSpec("seeing",              Float32, S_ALMANAC,  :gz),
    ColumnSpec("focus",               Float32, S_ALMANAC,  :gz),
    ColumnSpec("collpist",            Float32, S_ALMANAC,  :gz),
    ColumnSpec("colpitch",            Float32, S_ALMANAC,  :gz),
    ColumnSpec("dithered_pixels",     Float32, S_ALMANAC,  :gz),
    ColumnSpec("join_status",         UInt8,   S_ALMANAC,  :gz),
    ColumnSpec("era",                 UInt8,   S_ALMANAC,  :gz),
    ColumnSpec("sdss_id_agrees",      UInt8,   S_ALMANAC,  :gz),
    ColumnSpec("mjd_mid_exposure",    Float64, S_TIME,     :gz),
    ColumnSpec("v_barycentric_correction", Float64, S_BARY, :gz),
    ColumnSpec("v_rel",               Float64, S_VELOCITY, :gz),
    ColumnSpec("v_rad",               Float64, S_VELOCITY, :gz),
    ColumnSpec("e_v_rad",             Float64, S_VELOCITY, :gz),
    ColumnSpec("RV_pixoff_final",     Float64, S_WORKUP,   :gz),
    ColumnSpec("RV_pix_var",          Float64, S_WORKUP,   :gz),
    ColumnSpec("RV_flag",             Int64,   S_WORKUP,   :gz),
    ColumnSpec("ingestBit",           Int64,   S_WORKUP,   :gz),
    ColumnSpec("snr",                 Float64, S_WORKUP,   :gz),
    ColumnSpec("category",            UInt8,   S_ALMANAC,  :gz_nosh),
    ColumnSpec("hole_type",           UInt8,   S_ALMANAC,  :gz_nosh),
    ColumnSpec("image_type",          UInt8,   S_ALMANAC,  :gz_nosh),
    ColumnSpec("bary_flag",           UInt8,   S_BARYFLAG, :gz_nosh, ["bits" => VlenStr(BARY_FLAG_BITS)]),
    ColumnSpec("glon",                Float64, S_GALACTIC, :contig),
    ColumnSpec("glat",                Float64, S_GALACTIC, :contig),
    ColumnSpec("qa_bitmask",          Int32,   S_QA,       :contig),
    ColumnSpec("v_rad_corr",          Float64, S_FIBEROFFSET, :gz13045, ["description" => FixedStr(DESC_V_RAD_CORR)]),
    ColumnSpec("v_rad_offset_applied", Float64, S_FIBEROFFSET, :gz13045, ["description" => FixedStr(DESC_V_RAD_OFFSET)]),
    ColumnSpec("v_rad_corr_flag",     UInt8,   S_FIBEROFFSET, :gz13045, ["description" => FixedStr(DESC_V_RAD_CORR_FLAG)]),
    ColumnSpec("offtarget",           Int8,    S_OFFTARGET, :gz13045, ["description" => VlenStr(DESC_OFFTARGET)]),
]

"Root attributes of the 0.4.0 file. `astropy_version` is filled from the barycentric kernel."
function root_attrs(astropy_version::AbstractString; with_fiber_offsets::Bool)
    a = Pair{String,Any}[
        "description" => FixedStr("APOGEE DR21 pass-1 allExposure summary: almanac per-exposure and per-fiber " *
            "metadata joined to the arMADGICS workup rows, with both the observed " *
            "(topocentric) and barycentric velocity frames."),
        "n_rows" => nothing,   # filled by the writer
        "bary_method" => FixedStr("astropy SkyCoord.radial_velocity_correction(kind='barycentric'), evaluated at " *
            "mjd_mid_exposure_precise (UTC) with EarthLocation.of_site; computed via an exact " *
            "linear factorization in the target unit vector (6 basis directions per exposure)."),
        "astropy_version" => FixedStr(astropy_version),
        "velocity_convention" => FixedStr("v_rel = observed/topocentric (NO barycentric correction applied); " *
            "v_barycentric_correction = astropy barycentric correction, km/s; " *
            "v_rad = barycentric-corrected, composed in redshift as (1+z_tot)=(1+z)(1+z_corr), " *
            "matching almanac postprocess/radial_velocity.py finalize_radial_velocities()."),
        "fiber_index_convention" => FixedStr("fiber_index_almanac = 301 - (adjfiberindx - 300*(tele-1)); REVERSED relative to " *
            "adjfiberindx. Verified against almanac sdss_id on 100.00% of testable rows."),
        "join_status_codes" => FixedStr("1=ok 2=no_exposure_row 3=no_fiber_table 4=wrong_length 5=idx_out_of_range"),
        "era_codes" => FixedStr("1=plate 2=fps"),
        "sdss_id_agrees_codes" => FixedStr("1=agree 2=disagree 3=untestable(one side <=0)"),
        "tele_codes" => FixedStr("1=apo 2=lco"),
        "recommended_cut" => FixedStr("bary_flag & 0b1111 == 0 for barycentric science"),
        "galactic_convention" => FixedStr("glon/glat are Galactic longitude/latitude in degrees, from ICRS (ra, dec) via " *
            "astropy SkyCoord(...).galactic; NaN where ra or dec is not finite"),
        "qa_bitmask_codes" => FixedStr("bit0(1)=BAD_NIGHT (tele,mjd) in public/bad_nights.csv; " *
            "bit1(2)=CHRONIC_FIBER (tele,fiber) scope=chronic in public/bad_fibers.csv; " *
            "bit2(4)=BAD_FIBER_NIGHT (tele,mjd,fiber) scope=fiber_night; " *
            "bit3(8)=BAD_FIBER_RANGE (tele,fiber) scope=fiber_mjd_range with " *
            "mjd_start<=mjd<=mjd_end; bit4(16)=RV_RAILED (tele,mjd,expnum,adjfiberindx) in " *
            "public/rv_railed.csv. fiber = adjfiberindx - 300*(tele-1) (product-array " *
            "1-based fiber index, NOT fiber_index_almanac)."),
        "qa_bitmask_source" => FixedStr("pass1_dr21/arMADGICS/public/ exclusion lists, 2026-09-16"),
    ]
    with_fiber_offsets && push!(a, "fiber_offset_convention" => VlenStr(
        "v_rad_corr = v_rad - offset, with one offset per (fiber, telescope, era) from fiber_offsets.csv in " *
        "this directory. The offsets are an empirical calibration derived from stars observed on more than " *
        "one fiber, telescope, or era; the telescope and era differences are absorbed, so no further " *
        "correction should be applied on top. v_rad_corr is NaN on four uncalibratable LCO fibers; see " *
        "v_rad_corr_flag."))
    return a
end

"""
Root attributes of a 0.4.1 (default-numerics) build: the 0.4.0 set with the barycentric method
updated to say what computed it, plus the pinned IERS table (path, sha256, ranges, warnings).
Content text otherwise unchanged (the 0.4.1 content edits listed in the PR are separate decisions).
"""
function root_attrs_041(; with_fiber_offsets::Bool, bary::Symbol, astropy_version::AbstractString,
                        erfa::AbstractString, iers::AbstractDict)
    a = root_attrs(bary === :astropy ? astropy_version : "8.0.1"; with_fiber_offsets)
    method = bary === :astropy ?
        "astropy SkyCoord.radial_velocity_correction(kind='barycentric'), evaluated at " *
        "mjd_mid_exposure_precise (UTC) with EarthLocation.of_site and the pinned IERS table below; computed " *
        "via an exact linear factorization in the target unit vector (6 basis directions per exposure)." :
        "astropy 8.0.1 SkyCoord.radial_velocity_correction(kind='barycentric') reproduced in Julia with the " *
        "ERFA C library (builtin ephemeris: epv00, plan94, moon98), evaluated at mjd_mid_exposure_precise (UTC) " *
        "at the APO/LCO coordinates of astropy EarthLocation.of_site with the pinned IERS table below; " *
        "computed via an exact linear factorization in the target unit vector (6 basis directions per " *
        "exposure). Bit-identical to astropy 8.0.1 for all 144,381 exposures of the 0.4.0 table."
    for (i, p) in enumerate(a)
        p.first == "bary_method" && (a[i] = "bary_method" => FixedStr(method))
        p.first == "astropy_version" && bary === :julia &&
            (a[i] = "astropy_version" => FixedStr("8.0.1 (model reproduced in Julia/ERFA; astropy not run)"))
    end
    push!(a, "erfa_version" => FixedStr(erfa))
    push!(a, "iers_table" => FixedStr(string(iers["recorded_path"])))
    push!(a, "iers_table_sha256" => FixedStr(string(iers["sha256"])))
    push!(a, "iers_b_table" => FixedStr(string(iers["recorded_path_b"])))
    push!(a, "iers_b_table_sha256" => FixedStr(string(iers["sha256_b"])))
    push!(a, "iers_tables_combination" => FixedStr("astropy IERS_Auto: IERS-A (finals2000A) with its Bulletin B columns replaced by the IERS-B (EOP 20 C04) values; B where present, else A"))
    push!(a, "iers_table_mjd_range" => FixedStr(@sprintf("%.0f-%.0f; measured (non-predicted) values before MJD %.0f; final (Bulletin B / C04) values to MJD %.0f",
        iers["first_mjd"], iers["last_mjd"], iers["predictive_mjd"], iers["last_bulletin_b_mjd"])))
    w = iers["warnings"]
    push!(a, "iers_table_warning" => FixedStr(isempty(w) ? "none: every exposure time is inside the measured part of the table" :
                                              "WARNING: " * join(w, "; ")))
    push!(a, "numerics" => FixedStr("10^x: Julia Base.exp10; cos/sin: C libm; no Python in the build"))
    return a
end
