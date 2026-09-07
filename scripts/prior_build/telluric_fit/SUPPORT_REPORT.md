# Support derivation — measurements behind the frozen constants

2026-09-07.  Every number here is MEASURED unless labelled INFERRED.

Sample: 60 dome flats per telescope, evenly spaced through
`list_<tele>_full.txt` (sorted by mjd, exposure); LCO MJD 57801-61160 of 5457
exposures, APO MJD 57643-61158 of 11716.  300 fibers each.
"live" = `mask_1d != 0 & ivar_1d > 0` — the same condition the fit itself uses
to build weights, so the support cannot claim pixels the fit would discard.
Rule: a pixel is live in an exposure when >= `min_live_frac` (0.5) of fibers are
live; a pixel enters the support when it is live in >= `min_exposure_frac` (0.5)
of exposures.  Same thresholds as FINDINGS.md section 1d.

All pixel indices below are **0-based, half-open** `[start, stop)`, matching the
Python side.  Julia pixel *n* is Python index *n-1*.

## 1. Raw live runs, before erosion

| chip | LCO | APO |
|---|---|---|
| 0 | `[276, 3407)` | `[327, 3438)` |
| 1 | `[3637, 6237)` | `[3665, 6247)` |
| 2 | `[6422, 8515)` | `[6421, 8500)` |

The blue chip starts 51 px bluer at LCO than at APO.  The old literal
`CHIP_INDICES = [(320, 3450), (3655, 6260), (6410, 8510)]` is a single set of
numbers for both.

## 2. Derived support vs the old literal (thresholds 0.5/0.5, `edge_buffer=2`)

| chip | new LCO | new APO | old (both) | LCO change | APO change |
|---|---|---|---|---|---|
| 0 | `[278, 3405)` | `[329, 3436)` | `[320, 3450)` | blue **+42**, red −45 | blue −9, red −14 |
| 1 | `[3639, 6235)` | `[3667, 6245)` | `[3655, 6260)` | blue **+16**, red −25 | blue −12, red −15 |
| 2 | `[6424, 8513)` | `[6423, 8498)` | `[6410, 8510)` | blue −14, red **+3** | blue −13, red −12 |
| total px | 7812 | 7760 | 7835 | | |

Positive = pixels gained (data the old window refused); negative = pixels given
back (grid the old window claimed with no data behind it).  Every APO edge moves
inward: the old `+/- 20` pad was applied as **dilation**, and APO is the
telescope those nominal numbers were closest to, so at APO the pad was pure
over-extension.

Wavelength spans of the new support:
LCO 15132.8-15800.6 / 15852.0-16430.6 / 16473.8-16955.9 A;
APO 15143.4-15807.4 / 15858.1-16432.9 / 16473.6-16952.4 A.

### The 67 dead LCO prior pixels
FINDINGS.md's three runs, 1-based 276-320 / 3638-3655 / 8511-8514 (0-based
275-319 / 3637-3654 / 8510-8513), all sat outside the old window.  Of the 67:

* `edge_buffer=0` recovers **66**.  0-based px 275 is live in 0.28 of exposures,
  below the 0.5 threshold, so it is genuinely marginal rather than lost.
* `edge_buffer=2` (the default) recovers **61** — 2 px at each blue edge and 1
  at the red edge traded for boundary stability.

Recovering 61 vs 66 is a science-range question, not a correctness one: under
the new contract the remaining pixels are *declared* off-support and zeroed in
linear space, which is what apMADGICS did.  The pathology (a placeholder read as
a throughput) is gone either way.

## 3. `edge_buffer` sensitivity — why the default is 2, not 20

Total support pixels, thresholds fixed at 0.5/0.5:

| edge_buffer | 0 | 1 | **2** | 3 | 5 | 10 | 20 |
|---|---|---|---|---|---|---|---|
| LCO | 7824 | 7818 | **7812** | 7806 | 7794 | 7764 | 7704 |
| APO | 7772 | 7766 | **7760** | 7754 | 7742 | 7712 | 7652 |

INFERRED (mine, not AKS's): the two thresholds are what decide whether a pixel
has usable data; the buffer only absorbs the +/-1 px quantisation of the
threshold crossing, so it should be as small as that job allows.  The measured
live-fiber fraction at LCO's red edge is 0.90 / 0.88 / 0.66 / 0.33 at px 8509 /
8511 / 8513 / 8515 — `edge_buffer=5` would discard px 8510-8512 at ~88% live and
recover **none** of that dead run, and `edge_buffer=20` would cost 12 x 20 = 240
px of grid, 3.5x what the whole fix recovers.  One flag (`--edge-buffer`)
changes it; the gate in `validate_output.py` reads the value from the file's
own attrs, so it stays self-consistent whatever is chosen.

## 4. `min_exposure_frac` sensitivity, and an epoch shift at APO

Support bounds at `edge_buffer=2`:

| thr | LCO | APO |
|---|---|---|
| 0.20 | `[277,3412) [3637,6241) [6422,8518)` | `[297,3439) [3640,6248) [6401,8500)` |
| 0.35 | `[278,3406) [3638,6236) [6423,8513)` | `[299,3438) [3641,6246) [6402,8499)` |
| **0.50** | `[278,3405) [3639,6235) [6424,8513)` | `[329,3436) [3667,6245) [6423,8498)` |
| 0.65 | `[280,3404) [3640,6234) [6425,8511)` | `[332,3409) [3669,6223) [6425,8481)` |
| 0.80 | `[287,3403) [3646,6233) [6430,8511)` | `[334,3408) [3671,6222) [6426,8480)` |

LCO is threshold-insensitive.  APO jumps by ~30 px between 0.35 and 0.5, and by
another ~27 px between 0.5 and 0.65.  MEASURED cause — split the same 60 APO
exposures into MJD tertiles:

| tertile | MJD | chip 0 | chip 1 | chip 2 |
|---|---|---|---|---|
| 0 | 57643-59188 | `[329, 3440)` | `[3666, 6249)` | `[6422, 8501)` |
| 1 | 59220-60272 | `[295, 3410)` | `[3638, 6224)` | `[6399, 8482)` |
| 2 | 60312-61158 | `[298, 3413)` | `[3640, 6226)` | `[6401, 8484)` |

**The whole APO footprint shifts ~30 px blueward at MJD ~59200 and stays there**
— both edges of all three chips move together by the same amount, which is the
signature of a wavelength-solution / dispersion change, not a coverage change.
LCO drifts too, by ~+8 px in its last tertile (`[285, 3414)` vs `[276, 3407)`).

This answers FINDINGS.md section 8 item 3 (open: "per-fiber vs per-epoch
drift").  It has a consequence a single global support cannot escape: with
threshold 0.5 the APO support is the *older* epoch's footprint, so the ~40% of
APO exposures on the post-59200 solution have ~30 px of real blue data outside
it — the same class of defect as the LCO bug, one order of magnitude smaller and
epoch-dependent rather than permanent.  Threshold 0.35 takes the union instead
and gives those pixels a smooth ~30 px continuum extrapolation (bounded by the
1000-px Matern prior) for the exposures that lack data there — qualitatively
different from `exp(0)`, but still a declaration that is generous for those
exposures.

**This is AKS's call, not mine.**  The default ships at 0.5 because that is what
FINDINGS.md's gate G2 is written against, and 0.35 is one flag away.  The
structural fix is a per-epoch (or per-exposure) support, which the current
one-`design_matrix`-per-product layout cannot express; it belongs with the
per-fiber work in pass 2.

## 5. Per-fiber spread (context for the pass-2 per-fiber support)

Fibers whose own liveness gives three clean chip runs: 251/300 LCO, 217/300 APO.
Per-fiber chip bounds, p5 / p50 / p95:

| | LCO start | LCO stop | APO start | APO stop |
|---|---|---|---|---|
| chip 0 | 274 / 276 / 279 | 3405 / 3407 / 3408 | 326 / 327 / 334 | 3437 / 3438 / 3441 |
| chip 1 | 3636 / 3637 / 3642 | 6203 / 6237 / 6239 | 3663 / 3664 / 3670 | 6227 / 6246 / 6249 |
| chip 2 | 6419 / 6422 / 6471 | 8493 / 8515 / 8516 | 6420 / 6421 / 6422 | 8464 / 8500 / 8502 |

Mostly a few px, with long tails (LCO chip 2 start p95 = 6471, 49 px past the
median).  `--per-fiber` derives and emits an `(N_PIXELS, N_FIBERS)` support and
applies it through the weights (`cinv[f, p] = 0`, exact row selection at no
memory cost, columns untouched); it is **off by default** — per-fiber support is
FINDINGS.md section 4 item 4, orthogonal to this fix.  The machinery and its
tests are here so pass 2 does not have to reopen the basis question.

MEASURED, and the concrete reason it is not ready to ship: a per-fiber LCO fit
passes G1 and the mode-anchoring tests, but **fails G2 for 49 of 300 fibers**.
Those are exactly the fibers whose own liveness gives *six* runs rather than
three — an interior gap inside a chip — so `bounds_from_live` refuses them and
they fall back to the global bounds, which then claim pixels they do not have.
Per-fiber support therefore needs a **non-contiguous** support per fiber, i.e.
apMADGICS's `msknall` from `generate_poly_prior` (which survives in arM as dead
code), not a per-fiber `[start, stop)` per chip.  That is the pass-2 work item,
and it is now a measured requirement rather than a guess.

## 6. The frozen canonical intervals

`CANONICAL_CHIP_INTERVALS = ((240, 3500), (3600, 6300), (6350, 8560))`

Extremes they must contain (over telescope x epoch-tertile x fiber):

| chip | bluest start | reddest stop | canonical | slack |
|---|---|---|---|---|
| 0 | 274 (LCO per-fiber) | 3443 (APO per-fiber) | `[240, 3500)` | 34 / 57 |
| 1 | 3635 (LCO) | 6251 (APO per-fiber) | `[3600, 6300)` | 35 / 49 |
| 2 | 6399 (APO MJD>59200) | 8520 (LCO) | `[6350, 8560)` | 49 / 40 |

Disjoint and ordered; every chip gap falls entirely between two canonical
intervals.  `bounds_from_live` raises if a future epoch escapes them rather than
clipping silently.  Lengths 3260 / 2700 / 2210 vs the old windows' 3130 / 2605 /
2100, so the Matern length-scale conversion `s_pixels / N` moves by <5%.
