# Exposure-level science verdict (0.4.1): exposure_verdict and its users.
#
# The first two testsets need only HDF5 + Test, so this file also runs on its own:
#     julia --project=<env with HDF5> -e 'using HDF5, Test;
#         include("src/exposureGuard.jl"); include("test/exposure_guard.jl")'
# The runlist testset needs the full arMADGICS stack (ApogeeReduction) and runs
# from test/runtests.jl.

@isdefined(EXPSTATE_EXCLUDED) ||
    include(joinpath(@__DIR__, "..", "src", "exposureGuard.jl"))

@testset "exposure_verdict: each bit, NOTRUN, absent field, flagged_bad" begin
    @test EXPFLAG_NOTRUN == 0x04
    @test (EXPFLAG_NO_SCIENCE & EXPFLAG_NOTRUN) == 0x00
    @test (EXPSTATE_OK, EXPSTATE_EXCLUDED, EXPSTATE_UNKNOWN) == (Int8(0), Int8(1), Int8(2))

    # judged and fine
    @test exposure_verdict(0x00, 0) == (EXPSTATE_OK, "ok")
    # each exclusion bit alone, and both
    @test exposure_verdict(EXPFLAG_PREDICTED_BAD, 0) == (EXPSTATE_EXCLUDED, "predicted_bad")
    @test exposure_verdict(EXPFLAG_ENGINEERING, 0) == (EXPSTATE_EXCLUDED, "engineering")
    @test exposure_verdict(0x03, 0) == (EXPSTATE_EXCLUDED, "predicted_bad+engineering")
    # almanac flagged_bad alone, with AR bits, and with the decoration absent
    @test exposure_verdict(0x00, 1) == (EXPSTATE_EXCLUDED, "almanac_flagged_bad")
    @test exposure_verdict(EXPFLAG_ENGINEERING, 1) ==
          (EXPSTATE_EXCLUDED, "almanac_flagged_bad+engineering")
    @test exposure_verdict(nothing, 1) == (EXPSTATE_EXCLUDED, "almanac_flagged_bad")
    @test exposure_verdict(EXPFLAG_NOTRUN, 1) == (EXPSTATE_EXCLUDED, "almanac_flagged_bad")
    # never judged: UNKNOWN, not OK and not EXCLUDED
    @test exposure_verdict(EXPFLAG_NOTRUN, 0) == (EXPSTATE_UNKNOWN, "notrun")
    @test exposure_verdict(nothing, 0) == (EXPSTATE_UNKNOWN, "flags_absent")
    # engineering on an unjudged exposure (2^1|2^2 is a legal byte): EXCLUDED wins
    @test exposure_verdict(EXPFLAG_ENGINEERING | EXPFLAG_NOTRUN, 0) ==
          (EXPSTATE_EXCLUDED, "engineering")
    # object never judged by the object-frame rule: OK -> UNKNOWN, never hides a bit
    @test exposure_verdict(0x00, 0; object_unjudged = "object_rule_absent") ==
          (EXPSTATE_UNKNOWN, "object_rule_absent")
    @test exposure_verdict(EXPFLAG_PREDICTED_BAD, 0; object_unjudged = "object_rule_absent") ==
          (EXPSTATE_EXCLUDED, "predicted_bad")
    @test exposure_verdict(EXPFLAG_NOTRUN, 0; object_unjudged = "object_rule_absent") ==
          (EXPSTATE_UNKNOWN, "notrun")
    # any integer type is accepted (HDF5 hands back Int8/Int32/UInt8)
    @test exposure_state(Int32(1), Int8(0)) == EXPSTATE_EXCLUDED
    @test exposure_state(UInt8(0), Int32(0)) == EXPSTATE_OK

    # the combined-product cut: only EXCLUDED fails
    @test exposure_usable_for_combined(EXPSTATE_OK)
    @test exposure_usable_for_combined(EXPSTATE_UNKNOWN)
    @test !exposure_usable_for_combined(EXPSTATE_EXCLUDED)

    # object_rule_unjudged
    @test object_rule_unjudged("dark", false, "", NaN, NaN, NaN) === nothing
    @test object_rule_unjudged("object", false, "", NaN, NaN, NaN) == "object_rule_absent"
    @test object_rule_unjudged("Object ", true, "none", NaN, NaN, NaN) == "object_rule_no_evidence"
    @test object_rule_unjudged("object", true, "none", 0.99, NaN, NaN) === nothing
    @test object_rule_unjudged("object", true, "guard_exempt", NaN, NaN, NaN) === nothing
    @test object_rule_unjudged("object", true, "forest", 0.01, -1.0, 1.0) === nothing

    # policies
    @test_throws ErrorException check_exposure_flag_policy("off")
    for pol in EXPOSURE_FLAG_POLICIES
        @test exposure_fit_action(EXPSTATE_OK, "ok", pol) == "fit"
        @test exposure_fit_action(EXPSTATE_UNKNOWN, "notrun", pol) == "fit"
        @test exposure_fit_action(EXPSTATE_EXCLUDED, "almanac_flagged_bad", pol) ==
              "dropped_flagged_bad"
        @test exposure_fit_action(EXPSTATE_EXCLUDED, "almanac_flagged_bad+predicted_bad", pol) ==
              "dropped_flagged_bad"
    end
    @test exposure_fit_action(EXPSTATE_EXCLUDED, "predicted_bad", "flag") == "fit"
    @test exposure_fit_action(EXPSTATE_EXCLUDED, "predicted_bad", "ignore") == "fit"
    @test exposure_fit_action(EXPSTATE_EXCLUDED, "predicted_bad", "drop") == "dropped_policy"
    @test exposure_fit_action(EXPSTATE_EXCLUDED, "engineering", "drop") == "dropped_policy"

    c = count_exposure_states(Int8[0, 1, 1, 2, 0, 0])
    @test (c.n_ok, c.n_excluded, c.n_unknown) == (3, 2, 1)
end

# A night of 8 object exposures plus one dark, decorated in the AR #435 layout.
const EG_TELE = "apo"
const EG_MJD = "58588"   # plate era, same fiber-table layout as test/ingest.jl
function write_eg_almanac(path; decorate = true, object_rule = true, drop_row = 0)
    nexp = 9
    flags = UInt8[0x00, 0x01, 0x02, 0x03, 0x04, 0x00, 0x00, 0x06, 0x00]
    fbad = [0, 0, 0, 0, 0, 1, 0, 0, 0]
    imtype = vcat(fill("object", 8), ["dark"])
    rule = ["none", "forest", "none", "no_sky_light", "n/a", "none", "none", "n/a", "n/a"]
    psky = [0.99, 0.01, 0.98, 0.5, NaN, 0.99, NaN, NaN, NaN]
    h5open(path, "w") do f
        g = create_group(f, "raw/$(EG_TELE)/$(EG_MJD)/exposures")
        g["exposure"] = collect(1:nexp)
        g["image_type"] = imtype
        g["plate_id"] = fill(9999, nexp)
        g["config_id"] = fill(9999, nexp)
        g["n_read"] = fill(47, nexp)
        g["chip_flags"] = fill(7, nexp)
        g["flagged_bad"] = fbad
        gf = create_group(f, "raw/$(EG_TELE)/$(EG_MJD)/fibers/9999")
        gf["fiber_id"] = [301 - i for i in 1:300]
        gf["category"] = [i <= 12 ? "sky" : "science" for i in 1:300]
        gf["sdss_id"] = collect(1:300)
        gf["fiber_type"] = fill("APOGEE", 300)
        if decorate
            keep = [i for i in 1:nexp if i != drop_row]
            gc = create_group(f, "exposure_class/$(EG_TELE)/$(EG_MJD)")
            gc["exposure"] = collect(1:nexp)[keep]
            gc["exposure_flags"] = flags[keep]
            gc["exposure_class_pred"] = fill("sky", nexp)[keep]
            gc["exposure_class_prob"] = fill(0.9, nexp)[keep]
            gc["exposure_class_status"] = fill("ok", nexp)[keep]
            gc["exposure_class_p_sky"] = psky[keep]
            gc["exposure_class_oh_rel"] = fill(NaN, nexp)[keep]
            gc["exposure_class_line_z"] = fill(NaN, nexp)[keep]
            gc["exposure_class_tr_dyn"] = fill(NaN, nexp)[keep]
            object_rule && (gc["exposure_class_object_rule"] = rule[keep])
            gc["engineering"] = UInt8.((flags[keep] .& 0x02) .!= 0)
            gc["engineering_frac"] = fill(NaN, nexp)[keep]
            gc["engineering_carton"] = [((x & 0x02) != 0) ? "manual_fps_position_stars" : ""
                                        for x in flags[keep]]
            gc["engineering_basis"] = fill("none", nexp)[keep]
        end
    end
    return path
end

@testset "exposure_verdicts_for_night: decoration reader and exposure table" begin
    mktempdir() do dir
        alm = write_eg_almanac(joinpath(dir, "alm.h5"))
        imtype = vcat(fill("object", 8), ["dark"])
        fbad = [0, 0, 0, 0, 0, 1, 0, 0, 0]
        v = h5open(alm, "r") do f
            exposure_verdicts_for_night(f, EG_TELE, EG_MJD, collect(1:9), imtype, fbad)
        end
        @test v.reason == ["ok", "predicted_bad", "engineering", "predicted_bad+engineering",
            "notrun", "almanac_flagged_bad", "object_rule_no_evidence", "engineering", "ok"]
        @test v.state == Int8[0, 1, 1, 1, 2, 1, 2, 1, 0]
        @test all(v.present)
        @test v.decoration.has_object_rule
        @test v.decoration.prov["exposure_class_object_rule"][2] == "forest"
        @test v.decoration.prov["engineering_carton"][3] == "manual_fps_position_stars"

        # decoration written before AR #435 (no object rule): objects with a zero byte
        # are UNKNOWN, the dark is OK, exclusion bits still exclude
        alm2 = write_eg_almanac(joinpath(dir, "pre435.h5"); object_rule = false)
        v2 = h5open(alm2, "r") do f
            exposure_verdicts_for_night(f, EG_TELE, EG_MJD, collect(1:9), imtype, fbad)
        end
        @test v2.reason[1] == "object_rule_absent"
        @test v2.reason[9] == "ok"
        @test v2.state[2:4] == Int8[1, 1, 1]
        @test v2.decoration.prov["exposure_class_object_rule"] == fill("", 9)

        # an exposure with no decoration row: UNKNOWN (flags_absent), and the legacy
        # reader still returns 0 for it (unchanged prior-build contract)
        alm3 = write_eg_almanac(joinpath(dir, "hole.h5"); drop_row = 2)
        h5open(alm3, "r") do f
            v3 = exposure_verdicts_for_night(f, EG_TELE, EG_MJD, collect(1:9), imtype, fbad)
            @test v3.reason[2] == "flags_absent"
            @test v3.state[2] == EXPSTATE_UNKNOWN
            @test !v3.present[2]
            @test read_exposure_science_flags(f, EG_TELE, EG_MJD, collect(1:9))[2] == 0x00
        end

        # undecorated almanac: throws, unless the escape hatch is set; then UNKNOWN
        alm4 = write_eg_almanac(joinpath(dir, "undecorated.h5"); decorate = false)
        h5open(alm4, "r") do f
            @test_throws ErrorException exposure_verdicts_for_night(f, EG_TELE, EG_MJD,
                collect(1:9), imtype, fbad)
            withenv(ALMANAC_UNDECORATED_ENV => "1") do
                v4 = exposure_verdicts_for_night(f, EG_TELE, EG_MJD, collect(1:9), imtype,
                    fbad)
                @test v4.state == Int8[2, 2, 2, 2, 2, 1, 2, 2, 2]
                @test count(==("flags_absent"), v4.reason) == 8
            end
        end

        # the exposure table round trip
        rows = [(tele = "apo", mjd = EG_MJD, exposure = i, exposure_flags = v.flags[i],
                    flags_present = v.present[i], flagged_bad = fbad[i],
                    exposure_state = v.state[i], state_reason = v.reason[i],
                    action = exposure_fit_action(v.state[i], v.reason[i], "flag"),
                    n_fibers_queued = 0) for i in 1:8]
        p = write_exposure_guard_table(joinpath(dir, "sub", "exposure_guard.h5"), rows;
            policy = "flag", almanac_file = alm)
        h5open(p, "r") do f
            @test read(f["exposure_state"]) == v.state[1:8]
            @test read(f["state_reason"]) == v.reason[1:8]
            @test read(f["action"])[6] == "dropped_flagged_bad"
            @test read(f["flags_present"]) == fill(0x01, 8)
            @test attrs(f)["exposure_flag_policy"] == "flag"
            @test attrs(f)["combined_product_cut"] == "exposure_state != 1"
        end
        c = log_exposure_guard_summary(rows, "flag")
        @test (c.n_ok, c.n_excluded, c.n_unknown) == (1, 5, 2)
    end
end

if @isdefined(get_telemjd_runlist_and_guard)
    @testset "runlist under each exposure_flag_policy" begin
        mktempdir() do dir
            alm = write_eg_almanac(joinpath(dir, "alm.h5"))
            queued(rl) = sort(unique([r.expnum for r in rl]))

            rl, g = get_telemjd_runlist_and_guard(alm, EG_TELE, EG_MJD;
                accepted_fibtypes = ["sci"], exposure_flag_policy = "flag")
            # 6 is flagged_bad (never reduced); everything else is fitted and labeled
            @test queued(rl) == [1, 2, 3, 4, 5, 7, 8]
            st = Dict(r.expnum => r.exposure_state for r in rl)
            @test st[1] == EXPSTATE_OK && st[2] == EXPSTATE_EXCLUDED && st[5] == EXPSTATE_UNKNOWN
            @test all(r -> r.exposure_flags isa UInt8 && r.exposure_state isa Int8, rl)
            @test [r.exposure for r in g] == collect(1:8)       # the dark is not a candidate
            @test [r.action for r in g][6] == "dropped_flagged_bad"
            @test all(r -> r.n_fibers_queued == (r.exposure == 6 ? 0 : 288), g)
            @test g[2].exposure_class_object_rule == "forest"

            rl_d, g_d = get_telemjd_runlist_and_guard(alm, EG_TELE, EG_MJD;
                accepted_fibtypes = ["sci"], exposure_flag_policy = "drop")
            @test queued(rl_d) == [1, 5, 7]                    # UNKNOWN kept, EXCLUDED not
            @test count(r -> r.action == "dropped_policy", g_d) == 4

            rl_i, _ = get_telemjd_runlist_and_guard(alm, EG_TELE, EG_MJD;
                accepted_fibtypes = ["sci"], exposure_flag_policy = "ignore")
            @test queued(rl_i) == queued(rl)

            # the legacy entry point keeps the prior-build guard and the old tuple shape
            rl_old = get_telemjd_runlist_from_almanac(alm, EG_TELE, EG_MJD;
                accepted_fibtypes = ["sci"])
            @test queued(rl_old) == [1, 5, 7]
            @test keys(first(rl_old)) == (:tele, :mjd, :expnum, :adjfiberindx, :sdss_id)

            @test_throws ErrorException get_telemjd_runlist_and_guard(alm, EG_TELE, EG_MJD;
                exposure_flag_policy = "off")
        end
    end
end
