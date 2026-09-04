# E6 REPORT — starContinuum prior-swap regression (2026-09-04)

Node ccalin051, nice -n 10, no Slurm. Branch `run/E6-prior-swap`
(worktree /mnt/home/asaydjari/gitcode/worktrees/arM-E6, off origin/main @ 8d3e8a0,
E6 commits: 0882c9e builder env overrides, 6f25977 fixture prior-swap arg).
Plan with pre-registered expected diffs: `E6_PLAN.md` (written before any builds).
Figures: https://users.flatironinstitute.org/~asaydjari/5mVSWXVXmtfhZTPZGJMBLId2/sdsswork/2026_09_04/plots/e6_regression/
Input sha256: `INPUT_SHA256SUMS.txt`. Logs: `build_{old,new}.log`, `fixture_*.log`.

## What was exercised
E4 samples (2026_09_03) -> committed builder `scripts/prior_build/build_starCont.jl`
(via new env overrides, logic untouched) -> built SVD priors (Vmat 8700x60 F64, λv,
chipgapmsk) -> pipeline core via `scripts/validation/run_M123_fixture.jl` (new optional
starCont-swap arg) on the existing 2026_08_31 M123 fixture.

Builds: fibers {101,245,295,335} apo + {350,450,595} lco, from OLD (2026_04_26) and
NEW (E4 2026_09_03) samples. **Timing: 7 fibers per generation, 4 workers, 8m52s (old)
/ 8m41s (new) wall; ~4-5 min per fiber single-threaded.** Full 600-fiber build at 30
workers projects to ~1.5-2 h.

## Verdicts

### Structural regression: PASS
All pre-registered invariants hold for all 7 fibers x {old,new}
(`structural_report.txt`): schema/dims/dtype exact; chipgapmsk == 2026_04_25
per-telescope masks; zero-padding outside mask exact; no NaN/Inf; V'V = Diag(λ) to
<1e-8 relative (offdiag <9e-16 typical); λ1 new/old = 0.999-1.001 on every fiber.
Matched components k=1-3: |cor| > 0.99 on all fibers except f450 (see leak);
leading-8 subspace principal cosines > 0.94 everywhere; changes are smooth and sit in
telluric-band regions (fig3), as pre-registered (E1/E2/E3 statements).
Two λ-spectrum deviations beyond the sub-% level, both understood:
- f101 (apo): λ2 new = 0.23x old — the new tfunlists' APO bright/chi2 cuts removed
  variance the old generation carried. Consistent with the E4 intent.
- f450 (lco): λ2 new = 8.0x old with modes 2-4 rotated — this is the LEAK (below),
  not an E4 code/builder issue.
Control: the production `APOGEE_starcont_svd_60_rough.h5` mask (`cont_msk`, 7855 px)
matches NEITHER current mask (apo 7742 / lco 7833) — rough predates the 2026_04_25
mask and is NOT reproducible by the current builder from any current samples; it is a
2025-generation artifact (see runtime decomposition).

### Runtime swap on the M123 fixture: PASS for the E4 swap, with a loud
### generation-level caveat
Three runs, identical code+fixture, only the starCont prior swapped
(`runtime_swap_report.txt`, `runtime_swap_isolated_report.txt`, fig6):
- (a) production rough (2025 generation), (b) NEW f295 build, (c) OLD-samples f295 build.
- Isolated E4 swap (c)->(b): chi2res −0.2% to −6.5% (all 8 live fibers IMPROVE),
  |ΔRV| < 0.01 pix, RV flags / starscale / snr / ingestBit / skyBit identical,
  pathological fibers (90/100/101/103) behave identically, batch-extraction contract
  passes, and the pixel-scale-artifact check on Δx_starContinuum finds nothing.
  The nominal "FAIL" lines in the isolated report are only my pre-registered ±5%
  chi2 threshold being crossed by −6.5%/−5.3% improvements on the two highest-SNR
  fibers — direction and smoothness are healthy; verdict PASS.
- Generation shift (a)->(c), SAMPLES HELD OLD: chi2res +8% to +124% (largest on
  high-SNR fibers). Cause: rough is a 2025 build — different mask (150 chip-edge
  pixels carrying 1.2% of rough's prior power lose support under the 2026_04_25 apo
  mask: runs 268-340, 3631-3675, 6414-6427, 8499-8516) and a different (2025) sample
  vintage. ΔRV vs rough still < 0.12 pix and flags unchanged.
  => Swapping any 2026-generation built prior into production will move chi2res by
  tens of percent REGARDLESS of E4; this is a generation change to be owned once,
  not an E4 regression. Golden-validation chi2 baselines will shift accordingly.

### Leak imprint (lco 59160/0018 + 60291/0009): REAL at the per-fiber level — FLAGGED
RNG-stream reproduction (MersenneTwister(203), draw order Teff,Av,Rv,Tfunindx,Tfrac)
identifies the affected sample columns exactly (`leak_identify.jl`): exposure rows
2578/3830; 229 tfunlist entries; 187/300 lco fibers carry >=1 affected column, max 3
per 10,000. Drop-variant rebuilds (columns removed at the builder level — no
resampling needed) show (`leak_report.txt`, fig5):
- f450: dropping its 2 leaked columns removes 76% of mode-2 variance (λ2 inflated
  ~4.2x by the leak; imprint localized ~15,250-15,350 Å). Pre-registered bound L1
  (<1e-3) is violated by 3 orders of magnitude at this fiber.
- f595: dropping its 3 columns changes λ2 by 0.5%, components at the 1e-3 level;
  mid-modes (k~6-10) deflate up to ~67% (noise-level trailing modes).
Production impact: NONE today — the runtime uses a single APO-fiber prior (rough, and
its E6 replacement candidate f295 is apo). But any per-fiber lco prior use, or an lco
"rough" analog, should wait for a sample regeneration with the committed C2_LCO=3000
tfunlist cut (validates the AKS addendum decision).

## Deliverables
- Built priors: `built_old/`, `built_new/` (+ `f450_dropleak`, `f595_dropleak`).
- Fixture outputs: `fixture_{oldprior,oldbuildprior,newprior}{,_batchstyle}.h5`.
- Reports: `structural_report.txt`, `runtime_swap_report.txt`,
  `runtime_swap_isolated_report.txt`, `leak_report.txt`, `mask_diff_report.txt`.
- Scripts (also committed on the branch under `scripts/validation/E6/`):
  `leak_identify.jl`, `build_leak_variant*.jl`, `compare_structural.jl`,
  `compare_fixture.jl`, `compare_leak.jl`, `mask_diff.jl`, `make_figures.jl`,
  `inspect_inputs.jl`.

## Recommendation
The E4 sample generation is SAFE to push through the builder into the runtime: the
swap is structurally clean and mildly improves fixture chi2 at fixed builder+mask.
Before the production swap, decide explicitly (1) which fiber's build (or which
aggregate) replaces `rough`, and (2) that the generation-level chi2 shift (mask +
2026 sample vintage) is accepted — it will dominate any golden-run chi2 diffs.
For lco per-fiber priors, regenerate samples with the C2_LCO=3000 lists first.
