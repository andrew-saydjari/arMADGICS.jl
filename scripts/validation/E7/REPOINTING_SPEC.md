# E7 runtime repointing SPEC (pipeline.jl on main — do NOT apply until pass-1 swap is approved)

Target: `arMADGICS.jl/pipeline.jl` @ main 4ef02ca (line numbers from that commit).

## Edit 1 — prior dictionary (lines 140–142)

Current:
```julia
prior_dict["starLines_refLSF"] = prior_dir * "2025_07_31/prior_dump/APOGEE_stellar_kry_50_subpix_th_22500.h5"
# prior_dict["starLines_LSF"] = prior_dir*"2024_03_16/arMADGICS.jl/src/prior_build/starLine_priors_norm94_dd/APOGEE_starCor_svd_50_subpix_f" # DD Version
# prior_dict["starLines_LSF"] = prior_dir*"2024_02_21/arMADGICS.jl/src/prior_build/starLine_priors_norm94/APOGEE_stellar_kry_50_subpix_f" # TH Version
```

New:
```julia
prior_dict["starLines_refLSF"] = prior_dir * "2025_07_31/prior_dump/APOGEE_stellar_kry_50_subpix_th_22500.h5"
prior_dict["starLines_LSF"] = prior_dir * "2026_09_05/prior_outputs/starLines_perfiber/APOGEE_stellar_kry_50_subpix_f" # TH Version, E7 (new FPI LSFs)
```

`starLines_refLSF` MUST stay pointed at the existing 2025_07_31 file: the new
per-fiber priors are `K_fiber · Vout` of the same fullres `Vout`
(md5 0dcf0ea9072b0391741b84814d70c08e) that produced that refLSF file, which is
what keeps the asymmetric fit-at-fiber-LSF / report-at-refLSF coefficient
pairing valid. Do not regenerate or renormalize either side independently.

## Edit 2 — per-batch prior load (lines 547–549)

Current:
```julia
V_starlines = V_starlines_refLSF #hack
# f = h5open(prior_dict["starLines_LSF"]*lpad(adjfiberindx ,3,"0")*".h5")
# V_starlines = read(f["Vmat"])
```

New (un-hack; note the `close(f)` at line 556 is currently commented too):
```julia
f = h5open(prior_dict["starLines_LSF"] * lpad(adjfiberindx, 3, "0") * ".h5")
V_starlines = read(f["Vmat"])
```
and at line 556 restore `close(f)` — BUT only in the `else` branch of the
`ddstaronly` conditional, or restructure: with `ddstaronly = false` (pass-1)
the `msk_starCor` read at line 552 never runs, so the minimal safe form is:

```julia
V_starlines, msk_starCor = if ddstaronly
    f = h5open(prior_dict["starLines_LSF"] * lpad(adjfiberindx, 3, "0") * ".h5")
    V, m = read(f["Vmat"]), convert.(Bool, read(f["msk_starCor"]))
    close(f)
    V_starlines_refLSF = V # dd basis serves both sides
    V, m
else
    f = h5open(prior_dict["starLines_LSF"] * lpad(adjfiberindx, 3, "0") * ".h5")
    V = read(f["Vmat"])
    close(f)
    V, ones(Bool, length(chebmsk_exp))
end
```
(The dd branch stays inert for pass-1; it is written out only so the swap does
not silently break `ddstaronly = true` runs. The dd per-fiber files carry
`msk_starCor`; the TH files do not.)

## No other edits
The hack exists only in `pipeline.jl` (the inject scripts consume past-run
component outputs, not the prior files, and `run_M123_fixture.jl` mirrors of
the hack live on validation branches, not main).

## File contract (what the new files guarantee)
- path: `/mnt/ceph/users/sdssv/work/asaydjari/2026_09_05/prior_outputs/starLines_perfiber/APOGEE_stellar_kry_50_subpix_f<NNN>.h5`, NNN = adjfiberindx 001–600
- key `Vmat`: Float64 (8700, 50, 10); sub-pixel index 6 = zero shift,
  index i holds the basis for spectra shifted by (i−6)/10 pix (same convention
  as the refLSF file and `indTenth`).
- Components are NOT orthonormal (columns scaled by prior amplitude), exactly
  like the old norm94 files; consume with no further normalization.
