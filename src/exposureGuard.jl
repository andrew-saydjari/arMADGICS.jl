## Exposure-level science-safety guard (standalone: needs only HDF5)
#
# Kept in its own file so that lightweight prior-build scripts (e.g.
# scripts/prior_build/build_tfunlists.jl) can include ONLY this, without pulling in
# the full ingest/solve stack. src/ingest.jl includes it, so every existing
# includer keeps working unchanged.

## ---------------------------------------------------------------------------
## Exposure-level science-safety guard
##
## AKS 2026-09-08: "we want to be sure no prior builds after any first pass use
## data for which the flagged column is bad or engineering."
##
## Two independent exposure-level verdicts must be honoured before an exposure
## may contribute to a prior (or to any other science sample):
##
##   1. `raw/<tele>/<mjd>/exposures/flagged_bad` — the observers'/almanac flag.
##      Already applied (below) and unchanged.
##   2. `exposure_class/<tele>/<mjd>/exposure_flags` — ApogeeReduction's
##      exposure-level bitmask, written by
##      ApogeeReduction/scripts/cal/decorate_almanac_exptype.jl. See the
##      "Exposure-Level Flag Bits" table in the AR README. Bit 0 is the image
##      classifier's predicted-bad verdict; bit 1 is ENGINEERING (the
##      configuration's science fibers carry an engineering carton, e.g.
##      `manual_fps_position_stars*` FPS positioning frames).
##
## These bit values are mirrored from AR deliberately: arM must not gain a
## compile-time dependency on an AR version just to read two bits, and the
## almanac carries `exposure_flags_bits` as an attribute so a mismatch is
## detectable. Keep them in sync with AR's `src/exposureClassifier.jl`.
##
## The guard is LOUD by design. A silent filter is worse than none: nobody can
## then tell whether it ran. Every path below either logs its exclusions or
## throws.
## ---------------------------------------------------------------------------
const EXPFLAG_PREDICTED_BAD = 0x01
const EXPFLAG_ENGINEERING = 0x02
const EXPFLAG_NO_SCIENCE = EXPFLAG_PREDICTED_BAD | EXPFLAG_ENGINEERING
# 2^2: AR formed NO verdict for this exposure (classifier disabled, its table
# missing, or it threw). Mutually exclusive with 2^0; may co-occur with 2^1.
# Deliberately NOT in EXPFLAG_NO_SCIENCE: "never judged" is unknown, not bad.
const EXPFLAG_NOTRUN = 0x04

## ---------------------------------------------------------------------------
## Per-exposure science state (0.4.1): ONE verdict for every selection point
##
## Three sources say "do not do science with this exposure":
##   - EXPFLAG_PREDICTED_BAD (AR image classifier; for commanded objects the
##     2D object-frame rule, AR PR #435),
##   - EXPFLAG_ENGINEERING (AR engineering-carton check),
##   - almanac `flagged_bad` (observers' lists and the almanac override file,
##     e.g. the `not_science` "robots not moved" entries).
## `exposure_verdict` folds them into one state. Every place that selects
## exposures for a fit, a table or a combined product calls it; nothing else
## re-derives the answer from raw bits.
## ---------------------------------------------------------------------------
"Judged, and no exclusion source fired."
const EXPSTATE_OK = Int8(0)
"An exclusion source fired: must not contribute to any combined product."
const EXPSTATE_EXCLUDED = Int8(1)
"""
No verdict exists (AR 2^2 notrun, no decoration row, or a commanded object the
object-frame rule never judged). Kept, counted and labeled; never silently
treated as good or as bad.
"""
const EXPSTATE_UNKNOWN = Int8(2)

"""
    exposure_verdict(flags, flagged_bad; object_unjudged = nothing) -> (state, reason)

The single exposure-level science verdict.

- `flags`: the AR `exposure_flags` byte for this exposure, or `nothing` when the
  almanac decoration has no row for it (absent ⇒ UNKNOWN, never OK).
- `flagged_bad`: the almanac `raw/<tele>/<mjd>/exposures/flagged_bad` value.
- `object_unjudged`: `nothing`, or the reason a commanded object was never judged
  by the object-frame rule (`"object_rule_absent"`: the decoration predates AR
  #435; `"object_rule_no_evidence"`: every rule input was NaN). It can only turn
  OK into UNKNOWN; it never hides an exclusion bit.

Precedence: EXCLUDED > UNKNOWN > OK. EXCLUDED reasons join the sources that
fired with `+`, in the order `almanac_flagged_bad`, `predicted_bad`,
`engineering`. UNKNOWN reasons: `"flags_absent"`, `"notrun"`, or the
`object_unjudged` reason. Otherwise `"ok"`.
"""
function exposure_verdict(flags::Union{Nothing, Integer}, flagged_bad::Integer;
        object_unjudged::Union{Nothing, AbstractString} = nothing)
    parts = String[]
    flagged_bad != 0 && push!(parts, "almanac_flagged_bad")
    if !isnothing(flags)
        f = UInt8(flags)
        (f & EXPFLAG_PREDICTED_BAD) != 0x00 && push!(parts, "predicted_bad")
        (f & EXPFLAG_ENGINEERING) != 0x00 && push!(parts, "engineering")
    end
    isempty(parts) || return (EXPSTATE_EXCLUDED, join(parts, "+"))
    isnothing(flags) && return (EXPSTATE_UNKNOWN, "flags_absent")
    (UInt8(flags) & EXPFLAG_NOTRUN) != 0x00 && return (EXPSTATE_UNKNOWN, "notrun")
    isnothing(object_unjudged) || return (EXPSTATE_UNKNOWN, String(object_unjudged))
    return (EXPSTATE_OK, "ok")
end

"""
    exposure_state(flags, flagged_bad; object_unjudged = nothing) -> Int8

`first(exposure_verdict(...))`.
"""
exposure_state(flags, flagged_bad; object_unjudged = nothing) =
    first(exposure_verdict(flags, flagged_bad; object_unjudged = object_unjudged))

"""
    exposure_usable_for_combined(state) -> Bool

The cut that every COMBINED product (visit and star RVs and parameters, coadds,
star-level counts) and the recommended cuts apply to a per-exposure row:
anything except EXCLUDED. UNKNOWN passes on purpose (unknown is not bad); select
`state == EXPSTATE_OK` explicitly to keep only positively cleared exposures.
"""
exposure_usable_for_combined(state::Integer) = Int8(state) != EXPSTATE_EXCLUDED

"""
Values of `pipeline.jl --exposure_flag_policy`:

- `"flag"`: EXCLUDED exposures are fitted and kept in the per-exposure outputs
  with `exposure_state = 1`; combined products cut on
  `exposure_usable_for_combined`.
- `"drop"`: EXCLUDED exposures are not fitted (the arM #40 behavior, and the
  default of `get_telemjd_runlist_from_almanac`, so prior builds never see them).
- `"ignore"`: nothing is dropped for the AR bits (the behavior before any
  exposure exclusion); `exposure_state` is still computed and written.

Under every policy almanac `flagged_bad` exposures stay out of the fit, because
AR never reduces them to 1D; they are still listed in the exposure table.
"""
const EXPOSURE_FLAG_POLICIES = ("flag", "drop", "ignore")

function check_exposure_flag_policy(policy::AbstractString)
    policy in EXPOSURE_FLAG_POLICIES ||
        error("exposure_flag_policy must be one of $(EXPOSURE_FLAG_POLICIES), got \"$(policy)\"")
    return String(policy)
end

"""
    exposure_fit_action(state, reason, policy) -> String

What the runlist does with one exposure: `"fit"`, `"dropped_flagged_bad"` (no 1D
product exists) or `"dropped_policy"` (EXCLUDED under policy `"drop"`).
"""
function exposure_fit_action(state::Integer, reason::AbstractString, policy::AbstractString)
    check_exposure_flag_policy(policy)
    occursin("almanac_flagged_bad", reason) && return "dropped_flagged_bad"
    (Int8(state) == EXPSTATE_EXCLUDED && policy == "drop") && return "dropped_policy"
    return "fit"
end

"`(n_ok, n_excluded, n_unknown)` of a vector of `EXPSTATE_*` codes, for logging."
count_exposure_states(states) = (n_ok = count(==(EXPSTATE_OK), states),
    n_excluded = count(==(EXPSTATE_EXCLUDED), states),
    n_unknown = count(==(EXPSTATE_UNKNOWN), states))

"""
Escape hatch for running against an almanac that has not been decorated with the
exposure-level flags (`exposure_class/.../exposure_flags`). Set
`ENV["ARM_ALLOW_UNDECORATED_ALMANAC"] = "1"` to downgrade the hard error to a
loud warning. Intended ONLY for smoke tests and for reproducing pre-2026-09-08
runs — never for a prior build whose products will be used for science.
"""
const ALMANAC_UNDECORATED_ENV = "ARM_ALLOW_UNDECORATED_ALMANAC"

almanac_undecorated_allowed() = get(ENV, ALMANAC_UNDECORATED_ENV, "0") in
                                ("1", "true", "TRUE", "yes")

"""
    read_exposure_science_flags(f, tele, mjd, expnum)

Return the `EXPFLAG_*` bitmask for each exposure number in `expnum`, read from
the almanac's `exposure_class/<tele>/<mjd>` group (open file `f`).

Throws unless the group exists AND carries the `exposure_flags` dataset, because
without it there is no way to honour AKS's requirement that engineering
exposures never reach a prior build. `ARM_ALLOW_UNDECORATED_ALMANAC=1` downgrades
the throw to a loud `@warn` and returns all-zero flags (i.e. no exclusions) —
use only for smoke tests.

An `exposure_class` group that predates the engineering bit (it has
`predicted_bad` but no `exposure_flags`) is treated as undecorated: the
engineering verdict is genuinely absent from such a file and must not be
silently assumed to be "clean". Re-run the AR decoration step.
"""
function read_exposure_science_flags(f, tele, mjd, expnum)
    d = read_exposure_class_decoration(f, tele, mjd, expnum)
    nmiss = count(.!d.present)
    if d.decorated && nmiss > 0
        @warn "Exposure-level science guard: $(nmiss)/$(length(expnum)) exposures in " *
              "$(tele)/$(mjd) are absent from `exposure_class/$(tele)/$(mjd)`; they are " *
              "treated as unflagged here (exposure_state marks them UNKNOWN). " *
              "Re-decorate the almanac if this is unexpected."
    end
    return d.flags
end

# string/float provenance columns copied from the decoration when present
const DECORATION_STRING_COLS = ("exposure_class_pred", "exposure_class_status",
    "exposure_class_object_rule", "engineering_carton", "engineering_basis")
const DECORATION_FLOAT_COLS = ("exposure_class_prob", "exposure_class_p_sky",
    "exposure_class_oh_rel", "exposure_class_line_z", "exposure_class_tr_dyn",
    "engineering_frac")

"""
    read_exposure_class_decoration(f, tele, mjd, expnum) -> NamedTuple

Read the almanac decoration `exposure_class/<tele>/<mjd>` for the exposures
`expnum`, aligned with `expnum`:

- `flags::Vector{UInt8}`: the AR `exposure_flags` byte (0 where absent);
- `present::BitVector`: the decoration has a row for this exposure;
- `decorated::Bool`: the group and its `exposure_flags` dataset exist;
- `has_object_rule::Bool`: the decoration carries `exposure_class_object_rule`
  (written by AR #435 and later);
- `prov::Dict{String, Vector}`: the provenance columns in
  `DECORATION_STRING_COLS` / `DECORATION_FLOAT_COLS` ("" / NaN where absent).

Same contract as before for an undecorated almanac: throws, unless
`ARM_ALLOW_UNDECORATED_ALMANAC=1`, in which case it warns and returns
`present = false` everywhere (which `exposure_verdict` reads as UNKNOWN).
"""
function read_exposure_class_decoration(f, tele, mjd, expnum)
    n = length(expnum)
    prov = Dict{String, Vector}()
    for k in DECORATION_STRING_COLS
        prov[k] = fill("", n)
    end
    for k in DECORATION_FLOAT_COLS
        prov[k] = fill(NaN, n)
    end
    grp = "exposure_class/$(tele)/$(mjd)"
    problem = if !haskey(f, "exposure_class")
        "the almanac has no `exposure_class` group at all"
    elseif !haskey(f, grp)
        "the almanac has no `$(grp)` group"
    elseif !haskey(f[grp], "exposure_flags")
        haskey(f[grp], "predicted_bad") ?
        "`$(grp)` predates the engineering bit (has `predicted_bad` but no " *
        "`exposure_flags`), so the ENGINEERING verdict is simply absent" :
        "`$(grp)` has no `exposure_flags` dataset"
    else
        nothing
    end
    if !isnothing(problem)
        msg = "Exposure-level science guard cannot run for $(tele)/$(mjd): $(problem). " *
              "Decorate the almanac with " *
              "ApogeeReduction/scripts/cal/decorate_almanac_exptype.jl before building " *
              "priors, or set ENV[\"$(ALMANAC_UNDECORATED_ENV)\"]=\"1\" to proceed " *
              "UNGUARDED (bad/engineering exposures WILL enter the sample)."
        if almanac_undecorated_allowed()
            @warn "UNGUARDED PRIOR SAMPLE: " * msg
            return (flags = zeros(UInt8, n), present = falses(n), decorated = false,
                has_object_rule = false, prov = prov)
        end
        error(msg)
    end
    g = f[grp]
    ecexp = read(g["exposure"])
    ecflag = UInt8.(read(g["exposure_flags"]))
    row = Dict(zip(ecexp, eachindex(ecexp)))
    cols = Dict{String, Vector}()
    for k in (DECORATION_STRING_COLS..., DECORATION_FLOAT_COLS...)
        haskey(g, k) && (cols[k] = read(g[k]))
    end
    # an exposure missing from the decoration is a hole in the guard, not a pass
    flags = zeros(UInt8, n)
    present = falses(n)
    for (i, e) in enumerate(expnum)
        j = get(row, e, nothing)
        isnothing(j) && continue
        present[i] = true
        flags[i] = ecflag[j]
        for (k, v) in cols
            prov[k][i] = k in DECORATION_STRING_COLS ? String(v[j]) : Float64(v[j])
        end
    end
    return (flags = flags, present = present, decorated = true,
        has_object_rule = haskey(g, "exposure_class_object_rule"), prov = prov)
end

"""
    object_rule_unjudged(image_type, has_object_rule, object_rule, p_sky, oh_rel, line_z)

`nothing` when the object-frame rule judged this exposure (or it is not a
commanded object); otherwise the UNKNOWN reason passed to `exposure_verdict`:
`"object_rule_absent"` (the decoration predates AR #435, so commanded objects were
exempt rather than judged) or `"object_rule_no_evidence"` (the rule reported
`none` but every input that can fire it was NaN). `tr_dyn` is the guard, which
cannot fire the rule, so it is not consulted.
"""
function object_rule_unjudged(image_type, has_object_rule::Bool, object_rule,
        p_sky, oh_rel, line_z)
    lowercase(strip(String(image_type))) == "object" || return nothing
    has_object_rule || return "object_rule_absent"
    (String(object_rule) == "none" && isnan(p_sky) && isnan(oh_rel) && isnan(line_z)) &&
        return "object_rule_no_evidence"
    return nothing
end

"""
    exposure_verdicts_for_night(f, tele, mjd, expnum, image_type, flagged_bad)

Apply `exposure_verdict` to every exposure of one night (open almanac `f`).
Returns `(state::Vector{Int8}, reason::Vector{String}, flags::Vector{UInt8},
present::BitVector, decoration)`.
"""
function exposure_verdicts_for_night(f, tele, mjd, expnum, image_type, flagged_bad)
    d = read_exposure_class_decoration(f, tele, mjd, expnum)
    n = length(expnum)
    state = Vector{Int8}(undef, n)
    reason = Vector{String}(undef, n)
    for i in 1:n
        unj = d.present[i] ?
              object_rule_unjudged(image_type[i], d.has_object_rule,
            d.prov["exposure_class_object_rule"][i], d.prov["exposure_class_p_sky"][i],
            d.prov["exposure_class_oh_rel"][i], d.prov["exposure_class_line_z"][i]) :
              nothing
        state[i], reason[i] = exposure_verdict(d.present[i] ? d.flags[i] : nothing,
            flagged_bad[i]; object_unjudged = unj)
    end
    return (state = state, reason = reason, flags = d.flags, present = d.present,
        decoration = d)
end

"""
    exposure_science_exclusion_reasons(flags)

Break an `EXPFLAG_*` vector down into a `(n_predicted_bad, n_engineering,
n_both, n_excluded)` tally for logging.
"""
function exposure_science_exclusion_reasons(flags)
    nbad = count(x -> (x & EXPFLAG_PREDICTED_BAD) != 0 && (x & EXPFLAG_ENGINEERING) == 0,
        flags)
    neng = count(x -> (x & EXPFLAG_ENGINEERING) != 0 && (x & EXPFLAG_PREDICTED_BAD) == 0,
        flags)
    nboth = count(x -> (x & EXPFLAG_NO_SCIENCE) == EXPFLAG_NO_SCIENCE, flags)
    (n_predicted_bad = nbad, n_engineering = neng, n_both = nboth,
        n_excluded = nbad + neng + nboth)
end

"""
    almanac_exposure_science_flags(almanacFile, expkeys)

Look up the exposure-level flags for an arbitrary list of `(tele, mjd, expnum)`
keys (`mjd` as String or Int), for consumers that identify exposures by file path
rather than by walking the almanac night by night (e.g. the telluric
transfer-function lists, whose source exposures are recorded as ar1Dunical
paths).

Returns `(flags::Vector{UInt8}, flagged_bad::Vector{Bool}, missing_keys::Int)`.
Keys absent from the almanac get `flags = 0` and are counted in `missing_keys`;
the caller must decide what to do with those and must say so out loud.
"""
function almanac_exposure_science_flags(almanacFile, expkeys)
    flags = zeros(UInt8, length(expkeys))
    fbad = falses(length(expkeys))
    nmiss = 0
    cache_flags = Dict{Tuple{String, String}, Dict{Int, UInt8}}()
    cache_fbad = Dict{Tuple{String, String}, Dict{Int, Bool}}()
    h5open(almanacFile, "r") do f
        for (i, k) in enumerate(expkeys)
            tele, mjd, expnum = String(k[1]), string(k[2]), Int(k[3])
            tm = (tele, mjd)
            if !haskey(cache_flags, tm)
                grp = "raw/$(tele)/$(mjd)/exposures"
                if !haskey(f, grp)
                    cache_flags[tm] = Dict{Int, UInt8}()
                    cache_fbad[tm] = Dict{Int, Bool}()
                else
                    eg = f[grp]
                    en = Int.(read(eg["exposure"]))
                    cache_flags[tm] = Dict(zip(en,
                        read_exposure_science_flags(f, tele, mjd, en)))
                    cache_fbad[tm] = Dict(zip(en, read(eg["flagged_bad"]) .!= 0))
                end
            end
            fl = get(cache_flags[tm], expnum, nothing)
            if isnothing(fl)
                nmiss += 1
            else
                flags[i] = fl
                fbad[i] = cache_fbad[tm][expnum]
            end
        end
    end
    return flags, fbad, nmiss
end

"""
    almanac_science_exposure_census(almanacFile, tele_mjd_pairs; label)

Corpus-wide, up-front accounting of the exposure-level science guard, to be
called ONCE by every prior-build entry point before it assembles a sample.
Prints a banner naming the almanac, the number of candidate object exposures,
and how many are excluded for which reason, then returns the tally.

This exists so that "the guard ran and removed N exposures" is a visible fact in
every prior-build log, not something a reader has to infer from the absence of
warnings. It also fails fast (via `read_exposure_science_flags`) when the
almanac has not been decorated, before hours of compute are spent.
"""
function almanac_science_exposure_census(almanacFile, tele_mjd_pairs;
        label::AbstractString = "prior build")
    ntot = 0; nobj = 0
    nbad = 0; neng = 0; nboth = 0
    excluded = Tuple{String, String, Int, UInt8}[]
    h5open(almanacFile, "r") do f
        for (tele, mjd) in tele_mjd_pairs
            grp = "raw/$(tele)/$(mjd)/exposures"
            haskey(f, grp) || continue
            eg = f[grp]
            # read the columns directly: this file must stay free of the
            # ApogeeReduction/DataFrames stack so lightweight scripts can include it alone
            expnum = read(eg["exposure"])
            msk_obj = read(eg["image_type"]) .== "object"
            msk_obj .&= (read(eg["n_read"]) .> 3) .& (read(eg["chip_flags"]) .== 7) .&
                        (read(eg["flagged_bad"]) .== 0)
            flags = read_exposure_science_flags(f, tele, mjd, expnum)
            ntot += length(expnum)
            nobj += count(msk_obj)
            drop = msk_obj .& ((flags .& EXPFLAG_NO_SCIENCE) .!= 0x00)
            r = exposure_science_exclusion_reasons(flags[drop])
            nbad += r.n_predicted_bad; neng += r.n_engineering; nboth += r.n_both
            for i in findall(drop)
                push!(excluded, (String(tele), String(mjd), Int(expnum[i]), flags[i]))
            end
        end
    end
    nexcl = nbad + neng + nboth
    println("=" ^ 78)
    println("EXPOSURE-LEVEL SCIENCE GUARD ($(label))")
    println("  almanac: $(almanacFile)")
    println("  exposures in corpus:                $(ntot)")
    println("  object exposures passing flagged_bad/n_read/chip_flags: $(nobj)")
    println("  EXCLUDED by exposure_flags:         $(nexcl)")
    println("      predicted_bad only:             $(nbad)")
    println("      engineering only:               $(neng)")
    println("      both:                           $(nboth)")
    println("  remaining for the sample:           $(nobj - nexcl)")
    if !isempty(excluded)
        println("  excluded (tele mjd exposure flags):")
        for (t, m, e, fl) in excluded
            println("      $(t) $(m) $(e) 0x$(string(fl, base = 16, pad = 2))")
        end
    end
    println("=" ^ 78)
    flush(stdout)
    return (n_exposures = ntot, n_object = nobj, n_excluded = nexcl,
        n_predicted_bad = nbad, n_engineering = neng, n_both = nboth,
        excluded = excluded)
end

"""
    write_exposure_guard_table(path, guard_rows; policy, almanac_file)

Write the per-run exposure table (`exposure_guard.h5`): one row per candidate
object exposure, from `get_telemjd_runlist_and_guard`. One dataset per column;
file attributes record the policy, the almanac, and the meaning of the codes.
Rows whose `action` is not `"fit"` are the exposures that appear in no
per-(exposure, fiber) output; this table is their only record.
"""
function write_exposure_guard_table(path, guard_rows; policy::AbstractString,
        almanac_file::AbstractString = "")
    mkpath(dirname(path))
    cols = isempty(guard_rows) ? Symbol[] : collect(keys(first(guard_rows)))
    h5open(path, "w") do f
        for c in cols
            v = [getfield(r, c) for r in guard_rows]
            f[String(c)] = c == :flags_present ? UInt8.(v) : identity.(v)
        end
        a = attrs(f)
        a["exposure_flag_policy"] = String(policy)
        a["almanac_file"] = String(almanac_file)
        a["n_rows"] = length(guard_rows)
        a["exposure_state_codes"] = "0=ok 1=excluded 2=unknown (src/exposureGuard.jl)"
        a["exposure_flags_bits"] = "AR byte as read from the almanac decoration: " *
                                   "2^0=predicted_bad 2^1=engineering 2^2=notrun; 0 where flags_present=0"
        a["action_codes"] = "fit | dropped_policy (excluded, policy=drop) | " *
                            "dropped_flagged_bad (almanac flagged_bad; AR does not reduce these)"
        a["combined_product_cut"] = "exposure_state != 1"
    end
    return path
end

"""
    log_exposure_guard_summary(guard_rows, policy)

Print the run-level count of candidate exposures by state, reason and action.
"""
function log_exposure_guard_summary(guard_rows, policy::AbstractString)
    st = [r.exposure_state for r in guard_rows]
    c = count_exposure_states(st)
    tally(f) = sort(collect(_exposure_guard_tally(map(f, guard_rows))); by = x -> -x[2])
    println("=" ^ 78)
    println("EXPOSURE FLAGS (arMADGICS run, policy=$(policy))")
    println("  candidate object exposures: $(length(guard_rows))")
    println("  ok: $(c.n_ok)   excluded: $(c.n_excluded)   unknown: $(c.n_unknown)")
    println("  by reason:")
    for (k, n) in tally(r -> r.state_reason)
        println("      $(rpad(k, 40)) $(n)")
    end
    println("  by action:")
    for (k, n) in tally(r -> r.action)
        println("      $(rpad(k, 40)) $(n)")
    end
    println("  fitted but EXCLUDED (exposure_state=1, must be cut from combined products): " *
            "$(count(r -> r.exposure_state == EXPSTATE_EXCLUDED && r.action == "fit", guard_rows))")
    println("=" ^ 78)
    flush(stdout)
    return c
end

# dependency-free countmap (this file must stay includable with HDF5 alone)
function _exposure_guard_tally(v)
    d = Dict{eltype(v), Int}()
    for x in v
        d[x] = get(d, x, 0) + 1
    end
    return d
end
