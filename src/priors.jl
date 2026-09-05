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
#                           default <prior_dir>/2026_09_04/prior_outputs/sky_pass1/built
#   ARM_CHIPGAP_MSK         per-telescope chip-gap/cheb mask file (audit item 3)
#                           default <prior_dir>/2026_04_25/StarContChipGapMsk.h5

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

    # Star continuum (audit item 1): per-fiber pass-1b builds, split by telescope
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
        joinpath(prior_dir, "2026_09_04/prior_outputs/sky_pass1/built"))
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

    # Reference-LSF TH starLines prior (E7 will add the per-fiber starLines set;
    # the refLSF file remains the restframe-export basis even after E7)
    prior_dict["starLines_refLSF"] = joinpath(prior_dir,
        "2025_07_31/prior_dump/APOGEE_stellar_kry_50_subpix_th_22500.h5")

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

`ddstaronly=true` is refused loudly: the per-fiber DD starLines priors are an E7/
pass-2 deliverable (on pre-integration main this flag silently read `msk_starCor`
from a CLOSED file handle — latent crash, audit item 2 note).
"""
function load_fiber_priors(prior_dict, adjfiberindx; ddstaronly=false)
    if ddstaronly
        error("ddstaronly=true requires per-fiber DD starLines priors (E7 deliverable, " *
            "not yet wired; DD structurally needs a completed arM run and rides to " *
            "pass-2). Refusing to run rather than reading msk_starCor from an " *
            "unopened file (latent bug on pre-integration main).")
    end
    (1 <= adjfiberindx <= 600) || error("adjfiberindx=$adjfiberindx outside 1:600")
    tele_key = adjfiberindx > 300 ? "lco" : "apo"

    # per-telescope chip-gap/cheb mask (audit item 3)
    chebmsk_exp = h5open(prior_dict["chebmsk"]) do f
        convert.(Bool, read(f[tele_key]))
    end

    # starLines: reference-LSF TH prior for all fibers.
    # TODO(E7): replace with the per-fiber TH-with-new-LSF priors once branch
    # run/E7-starlines-perfiber delivers them (audit item 2; the
    # "V_starlines = V_starlines_refLSF #hack"); a follow-up commit lands it here.
    # The refLSF prior stays loaded regardless: V_starlines_refLSF[:,:,6] is the
    # restframe-export basis.
    f = h5open(prior_dict["starLines_refLSF"])
    V_starlines_refLSF = read(f["Vmat"])
    close(f)
    V_starlines = V_starlines_refLSF # E7 pending (see TODO above)
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

    skymsk_bright = chebmsk_exp # no bright submask exists (see docstring)
    skymsk_faint = chebmsk_exp .& submsk_faint
    # completely masking all bright lines b/c detector response is nonlinear
    skymsk = chebmsk_exp .& submsk_faint

    return (chebmsk_exp, skymsk_bright, skymsk_faint, skymsk, V_starcont,
        V_starlines_refLSF, V_starlines, msk_starCor, V_skycont, V_skyline_faint)
end
