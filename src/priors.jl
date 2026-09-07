## Prior path construction and per-fiber prior loading (pass-1 runtime integration).
# Shared by pipeline.jl and scripts/validation/run_M123_fixture.jl so the production
# loader and the fixture driver cannot drift (audit 2026_09_05 PER_FIBER_AUDIT.md,
# items 1-5).
#
# Path roots are env-overridable (for prior-swap regression runs and future
# generations); defaults point at the current pass-1 built sets:
#   ARM_STARCONT_PRIOR_DIR  starCont pass-1 per-fiber priors (audit item 1)
#                           default <prior_dir>/2026_09_05/prior_outputs/starCont_pass1c
#                           (pass-1c = AKS-approved final cut-policy regeneration,
#                           supersedes pass1b; identical schema)
#   ARM_SKY_PRIOR_DIR       E5 per-fiber sky priors (audit items 2/4)
#                           default <prior_dir>/2026_09_04/prior_outputs/sky_pass1/built_combined_telemaj_union
#                           (rebuild under the combined bright mask; the older
#                           `built/` set is the pre-mask baseline, retained for
#                           comparison only — its leading faint mode is ~99.8%
#                           a BRIGHT-pixel mode, so do not run pass-1 on it)
#   ARM_CHIPGAP_MSK         per-telescope chip-gap/cheb mask file (audit item 3)
#                           default <prior_dir>/2026_04_25/StarContChipGapMsk.h5
#   ARM_STARLINES_PRIOR_DIR E7 per-fiber TH starLines priors (audit item 2)
#                           default <prior_dir>/2026_09_05/prior_outputs/starLines_perfiber
#   ARM_STARLINES_REFLSF_HACK=1  fall back to the pre-E7 refLSF hack
#                           (V_starlines = V_starlines_refLSF) — regression use only
#   ARM_PRIOR_SUPPORT_RELTOL     relative-power threshold for the prior-support guard
#                           in load_fiber_priors (default 1e-2; set to 0 to DISABLE
#                           the guard for before/after regression comparisons)
#
# One further ARM_* variable lives outside this file; listed here so this comment stays
# the single index of runtime environment overrides:
#   ARM_SKY_CACHE_DIR       root of the per-exposure sky-bundle cache (src/skyCache.jl).
#                           UNSET (the default) = no caching, byte-identical behaviour.
#                           Nothing that cache stores depends on any prior set: every
#                           cached quantity is computed BEFORE the priors are touched.
#                           So repointing any variable above — including
#                           ARM_PRIOR_SUPPORT_RELTOL, whose guard narrows chebmsk_exp
#                           (and hence skymsk) on LCO fibers — cannot stale a cache
#                           entry, and the guard needs no cache invalidation.

"""
    build_prior_dict(prior_dir)

Build the prior path dictionary used by `load_fiber_priors`. `prior_dir` is the
prior root (production: /mnt/ceph/users/sdssv/work/asaydjari/); per-set roots can be
overridden with the ARM_* environment variables documented at the top of
src/priors.jl. Per-fiber entries are PREFIXES completed as
`prefix * lpad(adjfiberindx, 3, "0") * ".h5"`.
"""
function build_prior_dict(prior_dir)
    prior_dict = Dict{String,String}()

    # Star continuum (audit item 1): per-fiber pass-1c builds, split by telescope
    starcont_root = get(ENV, "ARM_STARCONT_PRIOR_DIR",
        joinpath(prior_dir, "2026_09_05/prior_outputs/starCont_pass1c"))
    prior_dict["starCont_apo"] = joinpath(starcont_root, "built_apo", "APOGEE_starcont_svd_60_f")
    prior_dict["starCont_lco"] = joinpath(starcont_root, "built_lco", "APOGEE_starcont_svd_60_f")

    # Sky priors (audit items 2/4): E5 per-fiber skycont + faint GSPICE skyline
    # priors. E5-output-contract notes: (i) no bright skyline priors are produced
    # ("not making bright priors for now, unrecoverable" — build_sky_defs.jl) and no
    # submsk_bright; bright-line pixels are excluded via the faint file's submsk,
    # exactly the DR17 consumption (skymsk = chebmsk & submsk_faint). (ii) the E5
    # skycont files carry no chebmsk_exp dataset (unlike DR17-era files); the
    # chip-gap mask is a separate input.
    sky_root = get(ENV, "ARM_SKY_PRIOR_DIR",
        joinpath(prior_dir, "2026_09_04/prior_outputs/sky_pass1/built_combined_telemaj_union"))
    prior_dict["skycont"] = joinpath(sky_root, "APOGEE_skycont_svd_30_f")
    prior_dict["skyLines_faint"] = joinpath(sky_root, "APOGEE_skyline_faint_GSPICE_svd_120_f")

    # Chip-gap/cheb mask (audit item 3): per-telescope 2026_04_25 masks (apo 7742 /
    # lco 7833 good px) — the same file every E4/E5 prior build trained against, so
    # runtime and priors are mask-consistent. Replaces the single 2025 global mask
    # (7783 px), which MASK_PROVENANCE.md (task #30, 2026_09_05/mask_revisit/) shows
    # descends from LCO FIBER 150's chip coverage applied to both telescopes; with
    # per-fiber priors wired, 102 APO pixels of that mask would be fit with NO
    # starCont prior support. The trace verified the 2026_04_25 masks are
    # bit-reproducible from current corpora and their excluded edge pixels are dead
    # in current data (live-fiber fraction 0.000).
    prior_dict["chebmsk"] = get(ENV, "ARM_CHIPGAP_MSK",
        joinpath(prior_dir, "2026_04_25/StarContChipGapMsk.h5"))

    # starLines (audit item 2, E7 wiring per scripts/validation/E7/REPOINTING_SPEC.md
    # on run/E7-starlines-perfiber):
    # - starLines_refLSF MUST stay pointed at the 2025_07_31 file: the per-fiber
    #   priors are K_fiber * Vout of the SAME fullres Vout (md5 0dcf0ea9...) that
    #   produced it, which keeps the fit-at-fiber-LSF / report-at-refLSF coefficient
    #   pairing in update_Ctotinv_Vstarstarlines_asym valid. Never regenerate or
    #   renormalize either side independently.
    # - starLines_LSF: E7 per-fiber TH priors built with the new FPI LSFs
    #   (get_lsf_matrix, MJD 60861 params); Vmat Float64 (8700, 50, 10),
    #   subpix index 6 = zero shift, amplitude-scaled (no renormalization).
    prior_dict["starLines_refLSF"] = joinpath(prior_dir,
        "2025_07_31/prior_dump/APOGEE_stellar_kry_50_subpix_th_22500.h5")
    starlines_root = get(ENV, "ARM_STARLINES_PRIOR_DIR",
        joinpath(prior_dir, "2026_09_05/prior_outputs/starLines_perfiber"))
    prior_dict["starLines_LSF"] = joinpath(starlines_root, "APOGEE_stellar_kry_50_subpix_f")

    return prior_dict
end

"""
    per_fiber_prior_file(prefix, adjfiberindx)

Complete a per-fiber prior path prefix and fail LOUDLY if the file does not exist
(so a missing/unbuilt prior is an explicit error, not an obscure downstream crash).
"""
function per_fiber_prior_file(prefix, adjfiberindx)
    fname = prefix * lpad(adjfiberindx, 3, "0") * ".h5"
    isfile(fname) || error("per-fiber prior file not found: $fname " *
        "(adjfiberindx=$adjfiberindx). Check the ARM_* environment overrides and " *
        "whether that fiber's prior has been built.")
    return fname
end

"""
    prior_support_mask(V, msk; reltol=prior_support_reltol())

Pixels of `V` (npix x ncomp, or npix x ncomp x nsub) that carry non-negligible prior
power, i.e. whose row norm is at least `reltol` times the median row norm over `msk`.

`sum(k, V[i,k]^2) == (V*V')[i,i]` is the prior VARIANCE at pixel i, so the row norm is
the prior's 1-sigma amplitude there; `reltol=1e-2` therefore rejects pixels whose prior
variance is below 1e-4 of typical. Such a pixel is not "unconstrained" — it is a
near-delta prior pinning the component to ~0, which is far worse.

Returns a Bool vector over all npix (pixels outside `msk` are judged on the same
criterion, and are zero there by construction, so the result is always a subset of the
built support).
"""
function prior_support_mask(V::AbstractArray, msk::AbstractVector{Bool};
        reltol=prior_support_reltol())
    npix = size(V, 1)
    length(msk) == npix || error("prior_support_mask: mask length $(length(msk)) != " *
        "prior npix $npix")
    any(msk) || error("prior_support_mask: mask selects no pixels, so there is no " *
        "reference level to judge prior power against")
    Vm = reshape(V, npix, :)
    rownorm = vec(sqrt.(sum(abs2, Vm, dims=2)))
    ref = median(view(rownorm, msk))
    isfinite(ref) && ref > 0 || error("prior_support_mask: median in-mask row norm " *
        "is $ref (prior carries no power on its own mask)")
    return rownorm .>= reltol * ref
end

"prior-support guard threshold; 0 disables the guard (see ARM_PRIOR_SUPPORT_RELTOL)."
prior_support_reltol() = parse(Float64, get(ENV, "ARM_PRIOR_SUPPORT_RELTOL", "1e-2"))

"""
    load_fiber_priors(prior_dict, adjfiberindx; ddstaronly=false)

Load all per-fiber priors for one adjusted fiber index (1-300 apo, 301-600 lco) and
return the `prior_vec` tuple consumed by `pipeline_single_spectra`:

    (chebmsk_exp, skymsk_bright, skymsk_faint, skymsk, V_starcont,
     V_starlines_refLSF, V_starlines, msk_starCor, V_skycont, V_skyline_faint)

Masking layout (DR17 consumption pattern, apMADGICS pipeline.jl):
- `submsk_faint` (from the E5 faint GSPICE skyline file) encodes
  obs>=min_obscnt & chipgap & faint-line-region; `skymsk = chebmsk_exp & submsk_faint`
  is the solve mask (bright sky lines excluded — nonlinear detector response).
  KNOWN E5 CALIBRATION GAP (2026-09-05, reported): the DR17-era bright/faint
  threshold is preserved in DR17 flux units and currently flags zero pixels in
  ar1Dunical units, so submsk == chip-gap mask until E5 recalibrates; the wiring
  here consumes submsk unchanged so a rebuild flows through.
- `skymsk_bright` is retained in the tuple for layout stability but equals
  `chebmsk_exp` (no per-fiber bright submask exists; the bright component is
  neither modeled nor exported).
- PRIOR-SUPPORT GUARD: the returned `chebmsk_exp` is the per-telescope mask
  INTERSECTED with the pixels where the delivered starCont and skyCont priors
  actually carry power (`prior_support_mask`), so `skymsk` — and hence the
  pipeline's `simplemsk` — cannot fit a pixel whose prior pins the component to
  ~0. Drops 67 px on every LCO fiber and 0 px on every APO fiber (APO is
  bit-identical to the raw mask); see the block comment at the call site for the
  measurement and `ARM_PRIOR_SUPPORT_RELTOL=0` to disable.

`ddstaronly=true` is refused loudly: the per-fiber DD starLines priors are a
pass-2 deliverable (the E7 TH files carry no `msk_starCor`; on pre-integration
main this flag silently read `msk_starCor` from a CLOSED file handle — latent
crash, audit item 2 note).
"""
function load_fiber_priors(prior_dict, adjfiberindx; ddstaronly=false)
    if ddstaronly
        error("ddstaronly=true requires per-fiber DD starLines priors (pass-2: DD " *
            "structurally needs a completed arM run as training data; the E7 TH " *
            "per-fiber files carry no msk_starCor). Refusing to run rather than " *
            "crashing on a missing dataset (on pre-integration main this path read " *
            "msk_starCor from an unopened file).")
    end
    (1 <= adjfiberindx <= 600) || error("adjfiberindx=$adjfiberindx outside 1:600")
    tele_key = adjfiberindx > 300 ? "lco" : "apo"

    # per-telescope chip-gap/cheb mask (audit item 3)
    chebmsk_exp = h5open(prior_dict["chebmsk"]) do f
        convert.(Bool, read(f[tele_key]))
    end

    # starLines (E7 landed): fit basis = per-fiber TH prior with the new FPI LSFs;
    # report/restframe basis = the refLSF prior (V_starlines_refLSF[:,:,6] is the
    # restframe-export basis; coefficient pairing requires the same parent Vout —
    # see build_prior_dict). ARM_STARLINES_REFLSF_HACK=1 restores the pre-E7
    # refLSF-for-all-fibers hack byte-identically (regression comparisons only;
    # E7 measured the hack costs +0.056±0.017 pix RV systematic on the M123 fixture).
    f = h5open(prior_dict["starLines_refLSF"])
    V_starlines_refLSF = read(f["Vmat"])
    close(f)
    V_starlines = if get(ENV, "ARM_STARLINES_REFLSF_HACK", "0") == "1"
        V_starlines_refLSF
    else
        fname = per_fiber_prior_file(prior_dict["starLines_LSF"], adjfiberindx)
        h5open(fname) do f
            read(f["Vmat"])
        end
    end
    # NOTE (deviation from REPOINTING_SPEC.md): the spec sketches an inert
    # ddstaronly branch reading msk_starCor from the per-fiber file; the E7 TH
    # files carry no msk_starCor (DD priors are pass-2), so ddstaronly stays a
    # LOUD refusal above instead of becoming a latent KeyError here.
    msk_starCor = ones(Bool, length(chebmsk_exp))

    # starCont (audit item 1): per-fiber pass-1 prior (files carry Vmat, λv,
    # chipgapmsk — NOT the old rough file's cont_msk / the DR17 chebmsk_exp). The
    # stored chipgapmsk must equal the runtime per-telescope mask: the build trained
    # against the same file, so a mismatch means mixed prior/mask generations.
    fname = per_fiber_prior_file(prior_dict["starCont_"*tele_key], adjfiberindx)
    V_starcont, starcont_chipgapmsk = h5open(fname) do f
        read(f["Vmat"]), convert.(Bool, read(f["chipgapmsk"]))
    end
    starcont_chipgapmsk == chebmsk_exp || error("starCont prior chipgapmsk in $fname " *
        "does not match the runtime per-telescope chebmsk " *
        "($(count(starcont_chipgapmsk)) vs $(count(chebmsk_exp)) good px): prior " *
        "generation and runtime mask are inconsistent (check ARM_CHIPGAP_MSK / " *
        "ARM_STARCONT_PRIOR_DIR).")

    # E5 per-fiber sky priors (audit items 2/4)
    fname = per_fiber_prior_file(prior_dict["skycont"], adjfiberindx)
    V_skycont = h5open(fname) do f
        read(f["Vmat"])
    end
    fname = per_fiber_prior_file(prior_dict["skyLines_faint"], adjfiberindx)
    V_skyline_faint, submsk_faint = h5open(fname) do f
        read(f["Vmat"]), convert.(Bool, read(f["submsk"]))
    end

    # PRIOR-SUPPORT GUARD (2026_09_07; evidence in <prior_dir>2026_09_07/prior_edge_check/).
    #
    # SYMPTOM: 67 in-mask pixels — IDENTICAL on all 300 LCO fibers, 0 on all 300 APO
    # fibers — carry ~1e-3 of the median starCont/skyCont row norm, i.e. ~1e-6 of the
    # prior VARIANCE, with a hard ~700-900x STEP (not a rolloff) at px 321/3656/8510:
    #     276-320 (45 px, blue chip blue edge), 3638-3655 (18 px, green chip blue
    #     edge), 8511-8514 (4 px, red chip red edge).
    #
    # ROOT CAUSE — the mask rule and the signal boundary key off DIFFERENT inputs.
    # A starCont sample is a product of three factors (sample_starCont.jl:148):
    #     tellFracSamples[:,i] .* Tfun./median(Tfun) .* (Ksp*(rvec.*bbs))
    # MEASURED, per factor, at the 67 px:
    #   * Ksp (LSF matrix): row sums are EXACTLY 1.0 at all 8700 px (lsf.jl:88-97
    #     renormalizes). Ruled out — it cannot suppress anything.
    #   * tellFracSamples: value ~0.98 there, so it does NOT suppress. But its
    #     EXACT-ZERO fraction goes 1.00 (px 265) -> 0.49 (270) -> 0.010 (275) ->
    #     0.0000 (276). The chipgapmsk rule is "any exact zero among the samples kills
    #     the pixel", so the mask edge IS this file's exact-zero edge — MEASURED to
    #     match on BOTH telescopes (lco 276==276, apo 341==341). These are the 2023-era
    #     files 2026_04_26/outsamptell_lco.jdat and 2026_04_25/outsamptell_apo.jdat.
    #   * Tfun = exp(Atell*theta): THE SUPPRESSOR. It collapses to ~0.002-0.004 of its
    #     median blueward of px 321 and steps ~400x to 0.73 at px 321 — and MEASURED,
    #     px 321 IS THE SAME BOUNDARY ON BOTH TELESCOPES. The telluric transfer
    #     function simply carries no information below it.
    # So the two boundaries are set by two DIFFERENT input files, and APO vs LCO is
    # decided purely by which side of px 321 the mask edge lands on:
    #     APO mask starts 341 = 20 px REDWARD of 321 -> safely inside support, 0 bad px
    #     LCO mask starts 276 = 45 px BLUEWARD of 321 -> 45 unsupported px (same at the
    #                                                    green/red chip edges: 18 + 4)
    # It is NOT that the LCO detector's usable range differs — the model support edge
    # is identical. NOT an E3 regression either: the pre-E3 delivered telluric products
    # (prior_inputs/tellurics_20260220_arjl_domeflats/20260323_lco.h5) give the SAME
    # px-321 boundary and the same 45 px gap, so this predates the E3 refit.
    # skyCont is NOT independent corroboration: sky_smooth_fit reconstructs each sample
    # as V_smooth_c*coeffs with Vcontinuum = the starCont Vmat (sample_sky_defs.jl:61-67,
    # 144-145), so it inherits this hole verbatim.
    #
    # WHY apMADGICS DIDN'T HAVE IT (build_starCont.jl:92 "This should probably go back
    # to being fiber dependent like I had in apMADGICS"): apMADGICS built chebmsk_exp
    # PER FIBER from measured per-fiber chip spans (medframes; prior_utils.jl:4-42,
    # sample_sky.jl:210-213) AND explicitly zeroed its telluric samples outside that
    # footprint (Vout[msknall,:] .= 0). arM's exp(Atell*theta) is strictly positive
    # everywhere, so that implicit per-fiber support — and the exact zeros the mask rule
    # was designed to detect — are gone. generate_poly_prior still exists at
    # scripts/prior_build/prior_utils.jl:4-42 but is DEAD CODE; no medframes file is
    # referenced anywhere in this repo.
    #
    # WHY THE RUNTIME DOESN'T CATCH IT: the data really is there. MEASURED on
    # ar1Dunical lco/61127 and 61130: frac(ivar>0)=0.997 and flux 0.63-0.81x the
    # interior median at these pixels. obscnt likewise counts real resampler COVERAGE
    # (cntvec .== framecnts, ApogeeReduction ar1D.jl:643-690), so 63 of the 67 sit
    # inside submsk_faint and skymsk/simplemsk keep them: MEASURED 36-48 entering the
    # solve in every one of 34 real (fiber, exposure) cases. The prior then pins
    # starCont/skyCont to ~0 while starLines (rel 1.1) and skyLine_faint (rel 3.2) keep
    # FULL power, so real continuum is forced into the LINE components. px 8511-8514
    # additionally sit under a bright sky line (3.6-13.7x the interior median).
    # The pre-integration prior set did not have this (min in-mask relative power 0.43),
    # so it is a pass-1 regression, not inherited.
    #
    # EFFECT OF THE FIX (lco/61127 exp 11+12, 34 cases; apo/61123+61130, 10 cases):
    #   SNR>=20: residual chi2 falls on 20 of 20 cases, by 26-95% (median -77%).
    #   SNR< 20: median +2.7%, i.e. a wash. The 8 cases that rise (all SNR<7, +2.3 to
    #            +15.3%) are fibers whose RV also jumps 8-66 px between the two runs —
    #            the RV is not determined at that SNR either way.
    #   APO:     10 of 10 cases BIT-IDENTICAL (same npix, drv = dchi2 = 0, apVisit equal).
    #
    # SELF-RETIRING: the threshold is derived from the delivered Vmat, so this becomes
    # a no-op the moment the masks/priors are rebuilt over real per-fiber support
    # (pass-2 task #30). It is deliberately runtime-only — no prior file is modified.
    # THRESHOLD: the population is bimodal with a ~60x-wide empty band (no pixel in any
    # of the 600 fibers has relative power in 1e-2..1e-1; APO min is 0.20), so any
    # reltol in ~3e-3..0.19 yields the identical mask.
    reltol = prior_support_reltol()
    if reltol > 0
        prior_supp = prior_support_mask(V_starcont, chebmsk_exp; reltol=reltol) .&
                     prior_support_mask(V_skycont, chebmsk_exp; reltol=reltol)
        ndrop = count(chebmsk_exp .& .!prior_supp)
        # NEVER silent: silent masking is how this class of bug hid in the first place.
        if ndrop > 0
            @info "load_fiber_priors: prior-support guard dropped $ndrop px " *
                  "(adjfiberindx=$adjfiberindx, $tele_key; $(count(chebmsk_exp)) -> " *
                  "$(count(chebmsk_exp .& prior_supp)) good px, reltol=$reltol)"
        end
        chebmsk_exp = chebmsk_exp .& prior_supp
    end

    skymsk_bright = chebmsk_exp # no bright submask exists (see docstring)
    skymsk_faint = chebmsk_exp .& submsk_faint
    # completely masking all bright lines b/c detector response is nonlinear
    skymsk = chebmsk_exp .& submsk_faint

    return (chebmsk_exp, skymsk_bright, skymsk_faint, skymsk, V_starcont,
        V_starlines_refLSF, V_starlines, msk_starCor, V_skycont, V_skyline_faint)
end
