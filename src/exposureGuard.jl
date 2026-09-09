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
            return zeros(UInt8, length(expnum))
        end
        error(msg)
    end
    ecexp = read(f["$(grp)/exposure"])
    ecflag = UInt8.(read(f["$(grp)/exposure_flags"]))
    lut = Dict(zip(ecexp, ecflag))
    # an exposure missing from the decoration is a hole in the guard, not a pass
    out = zeros(UInt8, length(expnum))
    nmiss = 0
    for (i, e) in enumerate(expnum)
        v = get(lut, e, nothing)
        if isnothing(v)
            nmiss += 1
        else
            out[i] = v
        end
    end
    if nmiss > 0
        @warn "Exposure-level science guard: $(nmiss)/$(length(expnum)) exposures in " *
              "$(tele)/$(mjd) are absent from `$(grp)`; they are treated as unflagged. " *
              "Re-decorate the almanac if this is unexpected."
    end
    return out
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
