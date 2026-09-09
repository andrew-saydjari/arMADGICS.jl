# arMADGICS <img src="docs/src/assets/logo.png" alt="MADGICS Logo" width="100" align="right"/>

[![Build Status](https://github.com/andrew-saydjari/arMADGICS.jl/actions/workflows/CI_lite.yml/badge.svg?branch=main)](https://github.com/andrew-saydjari/arMADGICS.jl/actions/workflows/CI_lite.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/andrew-saydjari/arMADGICS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/andrew-saydjari/arMADGICS.jl)


Pipeline for APOGEE spectra using Marginalized Analytic Dataspace Gaussian Inference for Component Separation (MADGICS).

Newer version of [apMADGICS.jl](https://github.com/andrew-saydjari/apMADGICS.jl) in order to support improved reductions from [ApogeeReduction.jl](https://github.com/andrew-saydjari/ApogeeReduction.jl).

## Installation

This is a pipeline. Not a package. It is not really meant to be installed. The pipeline.jl script is meant to be run. This repo documents the development and versions of the pipeline for transparency and reproducibility.

If you wish to download the code and have the dependencies required to run it installed, you can install directly from GitHub. 

```julia
import Pkg
Pkg.add(url="https://github.com/andrew-saydjari/arMADGICS.jl")
```

## Inputs

The pipeline runs on a list of tuples. Those tuples have the form `(runindx, release_dir, redux_ver, tele, field, plate, mjd, fiberindx)`. 
- `runindx` is arbitrary and is a linear index for file name incrementing when batching. 
- `release_dir` and `redux_ver` specify which version of the data to point at (see examples below). 
- `tele` is either `apo25m` or `lco25m` for the APO or LCO 2.5 meter locations for APOGEE (North and South, respectively). 
- `field` is a pointing identifier, `plate` is a fiber configuration identifier, and `mjd` is a 5-digit rough date index.
- `fiberindex` is the fiberindex (1-300) of the spectra on the CCD. Note that this index is reversed `(301-fibernum)` from the `fibernum` reported in the allVisit and allStar summary files. However, this index is the real index on the CCD... which is by far the more reasonable index to work in.

One example is 
```
(273, "sdsswork/mwm", "daily", "apo25m", "101689", "6455", "59817", 295)
```
which runs an observation from APOGEE-North that has been reduced from 3D -> 1D by the daily pipeline.

Another example is
```
(10, "dr17", "dr17", "apo25m0000010i", "180+60", "9667", "58126", 295)
```
which shows how injection tests are handled. This is an injection test into sky observations from DR17 taken on APOGEE-North. The trailing "i" indicates that it is an injection and the number after "apo25m" is a dummy index to prevent collisions of multiple injections into the same sky observation.

## gridSearch Module Flag Bits (surfaced as `RV_flag`)

There is still a (much smaller dimensional) space that MADGICS needs to sample over (e.g. radial velocity). We have a custom grid-sampler module to implement that sampling. The flag bits from that module are below. These are written to the output column **`RV_flag`** (`src/gridSearch.jl`; see `pipeline.jl`).

| Value         | Bit         | Meaning     |
| ----------- | ----------- | ----------- |
| 0     | -     | No problems       |
| 1     | 0     | Interpolated minimum not less than minimum (should not occur) |
| 2     | 1     | Minimum index at edge of grid for dimension 1 |
| 4     | 2     | Minimum index at edge of grid for dimension 2 |
| 8     | 3     | Finite difference Hessian beyond grid edge |
| 16    | 4     | Bad curvature of chi2 surface (can't invert full 2d Hessian)|
| 32    | 5     | Very bad curvature of chi2 surface (can't invert diagonal entries)|
| 64    | 6     | Spectrum never entered the RV scan; see `ingestBit` for why (`INGEST_FAIL_RV_FLAG`, set in `src/pipelineCore.jl`, not by gridSearch) |

Observed on the DR21 200-MJD testbed (1,622,474 spectra): 0 = 1,619,707;
1 = 4; 8 = 95; 10 = 878; 32 = 82; 64 = 1,708. Value 10 is bits 1+3 (grid edge
plus Hessian off grid). Value 64 is exactly the `ingestBit != 0` set.

Value 1 is documented above as "should not occur", and it **occurred 4 times**
(2 APO, 2 LCO). That is rare enough to have gone unnoticed and is worth running
down rather than assuming benign.

## Ingest Module Flag Bits (`ingestBit`)

Set by `validate_exposure` (`src/ingest.jl`) when a 1D uni-cal spectrum is checked
before entering the MADGICS solve. Bits marked **fatal** cause the spectrum to be
skipped (`ingest_fatal`); the others record that pixels were masked but the
spectrum was still fitted.

| Value | Bit | Meaning | Fatal |
| --- | --- | --- | --- |
| 0   | -   | No problems | |
| 1   | 0   | Runtime error while ingesting (`INGEST_RUNTIME_ERROR_BIT`) | **yes** |
| 2   | 1   | Flux is entirely NaN/zero | **yes** |
| 4   | 2   | Fewer than `INGEST_MIN_GOODPIX` good pixels after checks | **yes** |
| 8   | 3   | Non-finite flux inside the good mask (those pixels masked) | |
| 16  | 4   | Non-finite or non-positive ivar inside the good mask (those pixels masked) | |
| 32  | 5   | Tiny-ivar pixels masked (below `INGEST_TINY_IVAR_RELFAC` x median good ivar) | |
| 64  | 6   | `starscale0 = nanzeromedian(flux)` non-finite or <= 0 | **yes** |
| 128 | 7   | ApogeeReduction flagged this FIBER's relative throughput unusable on this exposure (`INGEST_RELTHRPT_BROKEN_BIT`) | opt-in |
| 256 | 8   | AR's per-fiber throughput flag was ABSENT from the `ar1Duni` file (`INGEST_RELTHRPT_UNKNOWN_BIT`) | |

`INGEST_FATAL_BITS = 2^0 | 2^1 | 2^2 | 2^6`, plus `2^7` when
`INGEST_RELTHRPT_FATAL` is on.

A skipped spectrum is written out with NaN products and `RV_flag = 64`
(`INGEST_FAIL_RV_FLAG`, `src/pipelineCore.jl`). Before bits 7 and 8 existed,
**`ingestBit != 0` and `RV_flag == 64` were equivalent** — verified exactly on
the DR21 200-MJD testbed, 1,708 spectra of 1,622,474, both directions. Bits 7
and 8 are informational by default and therefore **break that equivalence**: the
correct statement is now `ingest_fatal(ingestBit) <=> RV_flag == 64`.

### Per-fiber throughput (`relthrpt`, `bitmsk_relthrpt`, `ingestBit` bit 7)

ApogeeReduction measures each fiber's relative throughput off a dome flat, per
fiber and per exposure, and records a quality bitmask. **A fiber AR calls broken
is deliberately left UNSCALED by AR** — not zeroed, not masked, not NaN'd — so
its flux sits on an arbitrary scale and any chi2 computed from it is nonsense.
Both fields ride in the `ar1Duni*` file arM already reads, and arM ignored them
entirely until now.

`getExposure` (`src/ingest.jl`) now reads them and `pipeline_single_spectra`
writes two new per-spectrum columns:

| Column | Type | Meaning |
| --- | --- | --- |
| `relthrpt` | Float64 | AR's relative throughput for this fiber on this exposure; the worst finite value across chips. `NaN` if unavailable. |
| `bitmsk_relthrpt` | Int | AR's per-fiber quality bitmask, OR-ed across chips. **`-1` means the field was absent from the reduction: UNKNOWN, not good.** |

AR's bit table (mirrored as `AR_RELTHRPT_*` in `src/ingest.jl`): 1 = low-throughput
warn, 2 = broken (`relthrpt < 0.07`), 4 = no fluxing file (`relthrpt` forced to 1),
8 = `relthrpt` non-finite, 16 = too few good pixels to measure throughput.
**`AR_RELTHRPT_UNUSABLE_BITS = 2 | 8 | 16`** is the aggressive cut, and is exactly
the set AR refuses to flux-scale.

**How to mask chi2 analysis.** Cut on `(bitmsk_relthrpt & 26) != 0`, or
equivalently on `(ingestBit & 128) != 0`. Do NOT cut on bit 1 (warn) — those
fibers are fluxed normally and are fine. Treat `bitmsk_relthrpt < 0` as UNKNOWN
and report it separately rather than folding it into either bucket.

**What this does NOT change.** By default bit 7 is informational: the spectrum is
still solved, still written, and no exposure or fiber is dropped from the
reduction. Setting `ARM_RELTHRPT_FATAL=1` adds bit 7 to `INGEST_FATAL_BITS`, which
skips the solve for such fibers and writes NaN products with `RV_flag = 64`. That
is a real science-behaviour change and is offered as a switch, not assumed.

### Two notions of a broken fiber, and which is authoritative

There are two independent signals and they do not agree. They are kept separate
on purpose.

| | AR `bitmsk_relthrpt` | arM `ingestBit` bit 6 / `skyBit` bit 4 |
| --- | --- | --- |
| Measures | dome-flat throughput, per fiber, per exposure | the science spectrum's own median flux |
| Sees | a fiber that stops delivering lamp light | a fiber whose extracted flux is non-positive, for any reason |
| Blind to | a fiber dead on R or G but healthy on B (chip-B-only fluxing, see AR README); a partial, chip-localised loss that keeps the median above 0.07 | any fiber whose flux stays positive despite being badly wrong |

**AR's `bitmsk_relthrpt` is authoritative for "was this fiber flux-calibrated",**
because it is the flag AR itself acts on: it is the literal record of which fibers
were scaled and which were left alone. Any analysis asking "is this spectrum's
flux scale meaningful, and therefore is its chi2 meaningful" must use it.

**arM's `starscale0`/sky guards remain authoritative for "can this spectrum be
solved at all",** which is a different and narrower question, and they stay in
place as a last-resort net for whatever the dome-flat cut misses.

They should NOT be merged into a single flag. Merging would lose the distinction
between "AR knew and declined to flux it" and "arM found it unusable at solve
time", which are different failure modes with different fixes. Recording both,
per spectrum, is what lets the disagreement be measured.

## Sky Module Flag Bits (`skyBit`)

Set in `src/ingest.jl` while building the per-exposure sky model. This is an
**exposure-level** flag: every spectrum from an exposure carries the same value,
not just the sky fibers that triggered it.

| Value | Bit | Meaning |
| --- | --- | --- |
| 0   | -   | No problems |
| 1   | 0   | A sky fiber was excluded upstream (`SKY_EXCLUDED_FIBER_BIT`) |
| 2   | 1   | At least one sky fiber failed the z-cut (`SKY_ZCUT_FIBER_BIT`) |
| 4   | 2   | Too few sky fibers survived (`SKY_TOO_FEW_FIBERS_BIT`) |
| 8   | 3   | No sky fibers at all (`SKY_NO_FIBERS_BIT`) |
| 16  | 4   | At least one KEPT fiber has a non-positive median (`SKY_NEGSCALE_FIBER_BIT`) |
| 32  | 5   | The sky decomposition went non-finite (`SKY_NONFINITE_DECOMP_BIT`) |
| 64  | 6   | At least one candidate sky fiber excluded by AR's per-fiber throughput flag, before the z-cut (`SKY_RELTHRPT_FIBER_BIT`) |

Bits 4 and 8 are set together when an exposure has no usable sky fibers, so
`skyBit = 12` means the sky model was skipped entirely.

Note bit 16 is informational: those fibers are **kept**. It flags that a fiber
with a non-positive median entered the sky model, which is not by itself an
error but is worth screening on. Bit 64 is not informational: those fibers are
removed.

### The sky-guard verdict is NOT logged. Read it from the products.

`getSky4visit` used to print one line per TARGET FIBER for a verdict that is
EXPOSURE-level: measured on job 7001233, **479,570 lines covering 3,676 distinct
exposures** — a ~130x duplication that was 98% of the job log and made naive line
counts overstate the problem by the same factor. That print is now **suppressed
entirely**.

No information is lost. `skyBit` *is* the verdict and it is a per-spectrum column
of every batch product, so the exposure-level truth is recoverable exactly, per
exposure. Unlike a log it cannot be rotated, truncated, duplicated by an
append-mode resume, or silently invalidated by someone rewording a `println`.

To census it, use ApogeeReduction `test/regression/arm_census.jl` (invoked by
`arm_census.sh`), which reads `ingestBit` and `skyBit` out of the batch products.
**Report the unique-exposure count, not the per-spectrum count**: every spectrum
of an exposure carries the same `skyBit`, so the per-spectrum number is inflated
by the fiber multiplicity — exactly the distortion the old log duplication caused.

### Observed distribution (DR21 200-MJD testbed, 1,622,474 spectra)

Useful as a sanity baseline, not as a statement that these rates are acceptable.

| `skyBit` | APO | LCO |
| --- | ---: | ---: |
| 0  | 797,750 | 344,854 |
| 2  | 175,563 | 273,328 |
| 12 | 0       | 300 |
| 16 | 22,438  | 2,175 |
| 18 | 1,853   | 4,213 |

`skyBit = 2` alone accounts for 27.6% of all spectra and its SNR distribution is
indistinguishable from unflagged data (median 49.1 vs 43.4) — it records that the
sky screen did its job, not that the spectrum is suspect. All 300 `skyBit = 12`
spectra are a single exposure (LCO MJD 57802, exposure 215).

