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

    # Chip-gap/cheb mask: still the single global 2025 file (audit item 3 lands in a
    # later commit of this branch; kept here so the swap is isolated per stage)
    prior_dict["chebmsk"] = joinpath(prior_dir, "2025_07_31/prior_dump/chebmsk_exp.h5")

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
     V_starlines_refLSF, V_starlines, msk_starCor)

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

    # global chip-gap/cheb mask (per-telescope switch is audit item 3, later commit)
    f = h5open(prior_dict["chebmsk"])
    chebmsk_exp = convert.(Bool, read(f["chebmsk_exp"]))
    close(f)

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

    # starCont (audit item 1): per-fiber pass-1b prior (files carry Vmat, λv,
    # chipgapmsk — NOT the old rough file's cont_msk / the DR17 chebmsk_exp; the
    # chipgapmsk dataset becomes a runtime consistency check when the per-telescope
    # mask switch lands)
    fname = per_fiber_prior_file(prior_dict["starCont_"*tele_key], adjfiberindx)
    V_starcont = h5open(fname) do f
        read(f["Vmat"])
    end

    # sky masks: E5 per-fiber submasks land in the sky-wiring commit of this branch;
    # until then this preserves main's behavior (no submask)
    skymsk_bright = chebmsk_exp
    skymsk_faint = chebmsk_exp
    skymsk = chebmsk_exp

    return (chebmsk_exp, skymsk_bright, skymsk_faint, skymsk, V_starcont,
        V_starlines_refLSF, V_starlines, msk_starCor)
end
