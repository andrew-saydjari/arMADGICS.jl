# telluric_fit — dome-flat continuum + sky fit (telluric prior input)

Produces the `tellurics_*.h5` files that `sample_starCont.jl` reads as `Atell`
and `theta`.  Python/JAX; the rest of `prior_build/` is Julia, so it carries its
own `pyproject.toml` / `uv.lock` and is run with `uv run`.

Provenance: this is the successor to `20260323.py`, which lived only in dated
scratch directories (`2026_09_02/telluric_refit_full/`, and a pre-E3 copy).  It
is versioned here because `sample_starCont.jl` in this same directory consumes
its output and the two contracts have to move together.

## The contract with consumers

The output declares its own pixel support:

| dataset / attr | meaning |
|---|---|
| `design_matrix` `(N_PIXELS, n_fourier)` | continuum basis; rows outside the support are **exact zeros** |
| `support` `(N_PIXELS,)` or `(N_PIXELS, N_FIBERS)` bool | **where the model is defined** |
| `attrs.canonical_chip_intervals` | the frozen intervals the basis is built on |
| `attrs.support_bounds` | per-chip `[start, stop)` actually fit |
| `attrs.support_{min_live_frac,min_exposure_frac,edge_buffer}` | derivation thresholds |
| `attrs.support_input_list{,_sha256}`, `attrs.support_n_exposures_scanned` | what it was derived from |
| `attrs.off_support_fill` | `"zero"` (default) or `"nan"` |

`support == any(design_matrix != 0, axis=1)` is asserted at write time, in
`fit_domeflats.py` and again in `merge_shards.py`.

**A consumer MUST apply `support` in linear space.**  The adoption patch for
`sample_starCont.jl` (NOT applied here — it lands with the refit, in one
adoption commit alongside the mask rederivation):

```julia
# alongside the existing Atell read
Atell, tellsupport = h5open(tfun_path, "r") do f
    A = permutedims(read(f["design_matrix"]), [2, 1])          # (8700, n_fourier)
    s = haskey(f, "support") ? read(f["support"]) : trues(size(A, 1))
    # 1-D support reads back as (8700,); a per-fiber one as (n_fiber, 8700)
    A, ndims(s) == 2 ? permutedims(s, [2, 1]) : s
end

# in genModSamp, replacing  Tfunsample = exp.(Atell*theta)
Tfunsample = exp.(Atell * TfunSamplef["theta"][:, fiberindx, Tfunindx])
Tfunsample[.!(ndims(tellsupport) == 2 ? tellsupport[:, fiberindx] : tellsupport)] .= 0.0
```

The zeroing must happen **before** the `./nanzeromedian(Tfunsample)` on the next
line.

`exp()` cannot represent "no support": `exp(0) = 1` exactly, and a consumer that
skips this line reads a placeholder as a throughput of `1/median(Tfun)` (~0.003
for LCO).  That is the bug this rewrite exists to close.  Zeroing in linear
space restores the apMADGICS semantics that `nanzeromedian` (filters zeros) and
`build_starCont.jl`'s `mean(filter(.!iszero, Vred))` are already written for.

`--dead-row-fill nan` makes an ignoring consumer fail loudly instead, but NaN
propagates into `sum(starcont, dims=1)` in `build_starCont.jl`, making
`specsum .> 0` false for every sample column.  Audit runs only.

## Basis vs support

These are separate, and keeping them separate is the point.

* **Basis** — built once on `telluric_support.CANONICAL_CHIP_INTERVALS`, a
  FROZEN constant.  Mode `k`'s period and node positions are absolute pixel
  quantities, identical for every fiber, telescope and epoch, so `theta[:, :, k]`
  is directly comparable across all of them.  Changing these constants
  invalidates that comparison and must never be done to accommodate data.
* **Support** — derived per telescope (optionally per fiber) from the dome-flat
  stack.  It selects which **rows** participate, never which columns exist.
  Widening or narrowing it cannot move a node; `tests/` asserts this.

`20260323.py` conflated the two (`fourier_design_matrix_1d(ei - si, N_MODES)`
built the basis *on* the window), which is why a per-telescope window could not
be adopted there without silently redefining every mode.

## Running

```bash
# 0. once per telescope: derive the support from the dome-flat list
uv run python telluric_support.py --list list_lco_full.txt --telescope lco \
    --out support_lco.json --npz live_lco.npz --n-sample 60

# 1. the fit (one shard; see make_shards.py / run_shard.sh for the array)
uv run python fit_domeflats.py --support support_lco.json \
    --input shards/list_shard_5.txt --output out/shard_5_lco.h5 \
    --no-figures --s2 --s2-iters 50 --lambda-t 1e4 --update-t --s-pixels 1000

# 2. merge + gate
uv run python merge_shards.py --base . --date-tag 20260907
uv run python validate_output.py out/tellurics_refit_20260907_lco.h5 --live live_lco.npz
```

All shards of one telescope must be fit against the **same** support JSON;
`fit_domeflats.py` refuses to resume into an output with a different support and
`merge_shards.py` refuses to merge shards that disagree.

Frozen inputs read from the working directory (not versioned here; sha256 in
`input_sha256.txt`): `T_init.pkl`, `tellurics_init.pkl`, `nmf_sky_lines.pkl`.
`CH4.fits` / `CO2.fits` / `H2O.fits` are no longer read — the import-time
`get_telluric_models()` / `instrument_lsf_sparse_matrix()` call in `20260323.py`
built `K_lsf`, which that file never used.

## Tests

```bash
python tests/test_telluric_support.py      # numpy only, no data needed
```

Covers G1 (support == nonzero rows), G2 (six boundaries vs the measured live
region, with the old literal as a failing regression baseline), G4 (the exp(0)
placeholder canary), and mode anchoring — including a check that the *old*
construction fails the anchoring test, so the test has teeth.

`validate_output.py` runs G1/G2/G4 against a produced h5.
