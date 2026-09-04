# E6 PLAN — starContinuum prior-swap regression (written BEFORE any builds/runs)

Date: 2026-09-04. Node: ccalin051 (nice -n 10, no Slurm). Author: Claude (E6 agent) for AKS.
Code: arMADGICS.jl worktree /mnt/home/asaydjari/gitcode/worktrees/arM-E6,
branch run/E6-prior-swap off origin/main @ 8d3e8a0 (post PR #29 + #30).

## Chain under test
E4 samples (2026_09_03 tell_prior_disk, 600 x Float64 8700x10000)
  -> scripts/prior_build/build_starCont.jl (mask chip gaps, drop zero-sum sample
     columns, grand-mean normalize, C = V V'/nsamp, SVD, keep nsub=60,
     Vmat = U[:,1:60]*sqrt(S[1:60]) zero-padded back to 8700 pixels)
  -> HDF5 {Vmat (8700,60) Float64, λv (60), chipgapmsk (8700)}
  -> pipeline.jl / run_M123_fixture.jl read ONLY f["Vmat"] as V_starcont
     (production path: 2025_07_31/prior_dump/APOGEE_starcont_svd_60_rough.h5,
     a single file used for ALL fibers).

## Planned steps
0. Provenance control: rebuild from the OLD 2026_04_26 samples for candidate fiber
   295 (the default `runlist_range` in both sampler and builder scripts) and compare
   to the production rough file. If (near-)identical, "rough" = old fiber-295 build,
   and the correct swap counterpart is the NEW fiber-295 build. Comparison tolerance:
   sign-aligned column agreement; small deviations acceptable if builder inputs
   (chipgapmsk file) changed since 2025_07_31 — will report what is found.
1. Build OLD and NEW priors for fibers {101, 245, 295, 335} (apo) and {350, 450, 595}
   (lco) via the committed builder with env-override edits (no logic changes).
2. Structural regression NEW vs OLD built per-fiber priors + vs production rough.
3. Runtime swap: run scripts/validation/run_M123_fixture.jl on the existing
   2026_08_31/m123_fixture with (a) production rough, (b) NEW fiber-295 build,
   same code, same fixture. Compare chi2res, RV, starscale, continuum components.
4. Leak imprint: reproduce the per-fiber MersenneTwister(203) draw stream (cheap, no
   spectra), find sample columns whose Tfunindx maps to lco 59160/0018 or 60291/0009
   (exposure-index lookup via tellurics_refit_20260902_lco.h5 metadata /
   tfunlist_audit_lco.h5), and rebuild the most-affected lco fiber's prior with those
   columns DROPPED (builder-level column subset — no resampling needed). Compare.

## Expected-diff statements (pre-registered)

## Pre-build finding (inspection only, before any builds — 2026-09-04)
The production rough file keys are {Vmat, cont_msk, λv} with count(cont_msk)=7855,
which matches NEITHER the current 2026_04_25 StarContChipGapMsk.h5 apo mask (7742)
nor lco (7833). The rough prior therefore predates the current chip-gap mask and
cannot be bit-reproduced by the current builder from any sample generation; the
old-samples control build isolates builder+mask effects, and the OLD-vs-NEW
per-fiber comparison isolates the E4 sample change. Statement S2 is amended:
new builds' chipgapmsk must equal the 2026_04_25 per-telescope masks (7742/7833),
and their Vmat support (nonzero rows) will differ from rough's by construction.
Verified on rough: S3-S5 hold there (offdiag < 9e-16, zero-padding exact).

## Pre-build leak finding (RNG reproduction, no builds)
Exposure indices: 2578 (lco 59160/0018), 3830 (lco 60291/0009); 229 tfunlist
entries across 300 lco fibers (matches addendum). Reproducing the MersenneTwister(203)
draw stream: 187/300 lco fibers have >=1 affected sample column, max 3 per fiber,
415 affected columns total (of 3M lco sample columns). Drop-variant target:
adjfib 595, affected sample columns [668, 7106, 8208] -> rebuild with those
3 columns removed (builder-level column subset; no resampling). Expected effect
bounded at the ~3/10000 covariance-weight level (statement L1).

MUST NOT change (structural invariants; any violation = FAIL):
- S1. Output schema: HDF5 keys Vmat, λv, chipgapmsk; Vmat size (8700, 60),
  eltype Float64; λv length 60 positive & descending; chipgapmsk length 8700.
- S2. chipgapmsk identical to the one stored in the production rough file for the
  same telescope (input 2026_04_25/StarContChipGapMsk.h5 is unchanged by E4).
- S3. Rows of Vmat outside chipgapmsk exactly zero (builder zero-pads).
- S4. No NaN/Inf anywhere in Vmat or λv.
- S5. Column structure: Vmat = U*sqrt(S) with U orthonormal on the masked pixels,
  i.e. (Vmat'Vmat)[masked] = Diagonal(λv) to numerical precision — check
  offdiag(Vmat'Vmat)/sqrt(λi λj) < 1e-8.
- S6. Scale: grand-mean normalization makes Vmat O(1); leading sqrt(λ1) expected
  same order of magnitude as old (QA (c) showed sample medians agree to ~0.1-1%,
  correlations 0.99998+). λ1_new/λ1_old within a factor ~2 (expect ~1 +- tens of %).

EXPECTED to change (physically sensible; absence of change would also be notable):
- E1. Leading component shapes: smooth continuum x telluric-band structure; small
  coherent differences concentrated in telluric absorption bands (~15,850-16,000 and
  ~16,850-17,000 A regions) and near chip edges — the E3 telluric bug fix and the E1
  row-normalized on-the-fly LSF change exactly those imprints. Leading-component
  vector correlation |dot| expected > ~0.99 for k=1-3, decreasing with k
  (trailing components are noise-dominated and can rotate freely).
- E2. Subspace angles: principal angles between old/new leading-10 subspaces small
  (cos > ~0.95 for the well-separated leading modes); larger angles allowed at k
  where the old spectrum has near-degenerate λ.
- E3. Singular-value spectrum: same decay slope; per-mode changes of order the
  sample-level differences (sub-% to few %); NO new outlier modes, no flattening.
- E4. Sign of columns is arbitrary (SVD) — all comparisons sign-aligned.

Runtime swap (fixture, same code, only V_starcont swapped):
- R1. All 11 fixture fibers keep identical status (ok), identical ingestBit/skyBit,
  identical nSkyFibers, snr (prior-independent quantities).
- R2. chi2res changes by small relative amounts (expect <~1% level; the starCont
  component is a broad continuum absorbed largely by starscale) — pixel-scale
  artifacts in x_starContinuum_v0 or order-unity chi2 jumps = FAIL, flag loudly.
- R3. RV_pixoff_final: changes small (<~ the RV grid fine step 0.1 pix) for healthy
  fibers; RV_flag unchanged.
- R4. apVisit_absmax comparable; no new floored pixels beyond +-small count.
- R5. Pathological fibers (90 dead, 100 all-masked, 101 NaN-flux, 103 partial-NaN)
  behave identically (guards are prior-independent).

Leak imprint expectation:
- L1. ~229 leaked (fiber,exposure) entries / 1.55M lco entries; per fiber ~0-3 of
  ~5,175 tfunlist entries -> expected ~0-6 of 10,000 sample columns per fiber.
  Dropping them should perturb Vmat leading components at the < 1e-3 relative level
  (columns enter covariance with weight n/nsamp; brightness is normalized per-sample
  by nanzeromedian, so a bright domeflat's leverage is bounded by its SHAPE anomaly,
  not its flux). Anything visible at the component level = flag.

## Pacing
uptime checked before each heavy phase; skip/wait if 1-min load > 20.
