# AllExposure: the single-pass allExposure builder (Julia)

```bash
julia --project=scripts/release_tables/allexposure -t 8 scripts/release_tables/allexposure/bin/build_allexposure.jl \
    --full-list <raw run>/full_list_info.h5 --almanac <allobs_*.h5> --apred <apred root> \
    --workup <workup dir> --exclusion-lists <dir with bad_nights/bad_fibers/rv_railed.csv> \
    [--fiber-offsets fiber_offsets.csv] [--iers finals2000A.all --iers-b eopc04.1962-now] \
    --cache-dir <scratch> --workers 8 \
    --out allExposure_<name>.h5 --report allExposure_<name>.report.json
#   0.4.0 bit-regression only:  --numerics legacy040 --python <venv>/bin/python   (AVX-512 node)
#   astropy cross-check only:   --bary astropy --python <venv>/bin/python
julia --project=scripts/release_tables/allexposure scripts/release_tables/allexposure/test/runtests.jl
```

One row per arMADGICS workup row, in `full_list_info.h5` order. It replaces the 0.4.0 chain
`s1` (identity) + `s4` (times) + `s6` (build) + `s9` (bary_flag, in place) + `s10` (Galactic + QA, in
place) + `fiber_offsets/s10_final_build.jl` (copy-and-swap) + `offtarget_sidecar/s1,s2,s4` (copy-and-swap)
with one program that writes the final table once. With `--numerics legacy040` it reproduces the released 0.4.0
table bit for bit: all 52 datasets, values, dtypes, chunking, filters, attributes and their HDF5 string
types (2026-11-24, `validation/compare_h5.jl`: ALL EQUAL). The default (0.4.1) build runs no Python.

## Layout

| file | contents |
|---|---|
| `src/schema.jl` | **the schema**: every column's name, element type, producing stage, on-disk layout class and attributes; the root attributes. The writer and the tests iterate it; nothing else lists columns |
| `src/identity.jl` | stage IDENTITY (full_list_info.h5) and TIME (mid-exposure-time harvest, Distributed); caches in the 0.4.0 file formats |
| `src/almanac.jl` | stage ALMANAC and the almanac half of OFFTARGET: parallel per-night reads, ordered apply |
| `src/derived.jl` | row-wise stages: BARY, VELOCITY, BARYFLAG, GALACTIC, QA, FIBEROFFSET; the `Numerics` types |
| `src/bary.jl` | the barycentric correction: astropy 8.0.1 `radial_velocity_correction` ported to Julia on ERFA_jll; IERS table reader |
| `iers/` | the pinned Earth-orientation tables (`finals2000A.all`, `eopc04.1962-now`), `SHA256SUMS`, `fetch_iers.sh` (provenance) |
| `src/build.jl` | `BuildConfig`, `build` (stage orchestration, timing, report), `write_table` |
| `src/h5util.jl` | layout-faithful dataset writes; fixed vs variable-length string attributes |
| `python/legacy_kernels.py` | numpy `10**x` (`--numerics legacy040` only) and the astropy barycentric basis (`--bary astropy` only) |
| `bin/build_allexposure.jl` | CLI |
| `test/runtests.jl` | unit tests (+ a full-table bitwise regression when `ALLEXPOSURE_REFERENCE`/`ALLEXPOSURE_CANDIDATE` are set); `test/fixture_rows_040.csv` = 400 released rows, `test/fixture_bary_040.csv` = astropy's barycentric basis for 240 exposures, as exact hex floats |
| `legacy_python/` | the 0.4.0 reference scripts (Python + the identity script), used by the regression only |

Stages (each one former 0.4.0 step, in `@enum Stage`): IDENTITY, ALMANAC, TIME, BARY, WORKUP,
VELOCITY, BARYFLAG, GALACTIC, QA, FIBEROFFSET (only with `--fiber-offsets`), OFFTARGET.

## Parallel design

Measured on the 0.4.0 inputs, the cost is in two places: the mid-exposure-time harvest (144,381
small `ar1Dunical` files: open latency on Ceph) and the almanac join (3,657 per-night groups of a
5.7 GB HDF5 file). Both are HDF5 reads, and HDF5.jl serializes every HDF5 call inside a process, so
threads cannot overlap them. The builder therefore uses **Distributed worker processes for the
reads** (`pmap` over nights; each worker opens files read-only and returns small typed tables, never
per-row data), **threads for the row-wise stages**, and **one ordered write by the head process**.

The arMADGICS workup's MPI pattern (`workup/workup_mpi.jl`, many ranks writing disjoint row ranges of
the shared output with `CEPHTWEAKS_LAZYIO=1`, `HDF5_USE_FILE_LOCKING=FALSE`) was considered and not
used: here the write is not the bottleneck (about 85 s for 2.7 GB, single writer, gzip 1), while
parallel writes of the chunked gzip+shuffle layout (which bit-identity of the file layout requires)
would need collective MPI-IO through an MPI-enabled HDF5 build and an MPI launch for a job that
otherwise fits one process. A single writer also makes row order trivially deterministic. The two Ceph variables are still set by the driver (concurrent readers of the same
almanac file). Workers: `--workers`, else `SLURM_CPUS_ON_NODE` under Slurm, else 8 (never `nproc`).
Multi-node runs are not supported (one node is enough: 26.7 M rows fit in ~10 GB).

Possible next step if the write ever matters: compress chunks in parallel threads and write them
with `H5Dwrite_chunk` (direct chunk I/O), keeping a single writer.

## Numerics and the barycentric correction

| quantity | 0.4.0 implementation | default (0.4.1) | `--numerics legacy040` (regression only) |
|---|---|---|---|
| `10^(RV_pixoff_final*6e-6)` | numpy `**`: Intel SVML on AVX-512 (1-ulp differences from correctly rounded pow in ~7 % of rows) | `Base.exp10` | numpy, via `legacy_kernels.py pow10` |
| barycentric basis U (4 numbers per exposure) | astropy `radial_velocity_correction` for 6 directions + `pinv`, IERS table chosen by astropy at run time | `src/bary.jl`: the same computation in Julia on ERFA's C library, pinned IERS tables | same as default |
| `cos`/`sin` of the target direction | numpy float64 (= glibc libm, verified) | `ccall` to libm | same |
| `glon`, `glat` | astropy ICRS -> FK5(J2000) -> Galactic | ERFA_jll with astropy's matrices (exact hex floats), scalings and wrap | same |
| everything else | numpy element-wise arithmetic | same operation order in Julia (no FMA contraction, `x^2` = `x*x`, IEEE `sqrt`) | same |

**Default vs 0.4.0.** Only the `10^x` kernel differs, so only `v_rel`, `v_rad`, `e_v_rad` and
`v_rad_corr` change (numbers in `../README.md`, "0.4.1 defaults"). The default gives the same bits on
every CPU: `Base.exp10` over all 26,715,540 rows and the barycentric basis were bit-identical with
Julia code generated for the native CPU (Sapphire Rapids), `haswell`, `znver2` (rome) and generic
`x86-64` (no FMA) (`julia -C <target>`, 2026-11-24). glibc's `cos`/`sin` take the same code path on all
Rusty CPU types (all have FMA). `legacy040` reproduces 0.4.0 only on AVX-512 nodes.

**Barycentric correction in Julia (`src/bary.jl`).** A port of astropy 8.0.1
`SkyCoord.radial_velocity_correction(kind="barycentric")` for a target at infinity: astropy `Time`
scale conversions (`day_frac`, UTC->TAI->TT->TDB with `eraDtdb` at the site, UTC->UT1), IERS
interpolation (`IERS._interpolate`), the observer's GCRS position and velocity (`eraC2tcio`,
`eraC2i06a`, `eraPom00`, `eraSp00`, `eraEra00`), Earth's barycentric velocity (`eraEpv00`), the
gravitational redshift of Earth, Moon, Jupiter and Sun (`eraMoon98`, `eraPlan94`, `eraEpv00`), astropy's
unit-conversion factors as exact hex floats, and the 0.4.0 6-direction fit with numpy's `pinv(A)`.
Matrix products accumulate with `fma` in k order, as numpy's `@` does through OpenBLAS on this
hardware; Julia's `fma` is correctly rounded on every CPU. Checked against astropy with the same
inputs on all 144,381 exposures of 0.4.0: every intermediate (TT/TDB/UT1 two-part JDs, UT1-UTC, polar
motion, GCRS position and velocity, Earth velocity, gravitational redshift), all 866,286 per-direction
astropy values and all 577,524 basis coefficients are **bit-identical**. Time: 15 s on 1 thread, 2.6 s
on 8 threads, against 140 s for astropy. The direct analytic coefficients (`bary_coeffs_direct`,
without the 6-direction fit) differ from the fit by up to 2.0e-10 km/s; they are not used.
Caveat: ERFA_jll 2.0.0 carries its own leap-second table (pyerfa 2.0.1 gets astropy's at run time);
both are current through the 2017 leap second, which covers every exposure, but a future leap
second would need an ERFA_jll update.

**IERS tables (`iers/`).** The barycentric correction depends on UT1-UTC and polar motion. 0.4.0 used
astropy's `IERS_Auto` and recorded no table. Searching on 2026-11-24, the tables of
`astropy-iers-data` 0.2026.9.7.0.56.14 (IERS-A `finals2000A.all` with IERS-B `eopc04.1962-now`
substituted, as `IERS_Auto` does) reproduce every 0.4.0 `v_barycentric_correction` bit; the
installed 0.2026.8.31 differs in 1,387 rows of the last 7 nights (the "one difference" of the first
version of this PR). Both files are vendored (3.76 MB + 5.17 MB; text inputs, not products) so a
build needs no network and the pin travels with the code; `fetch_iers.sh` re-derives them from the
PyPI wheel and checks `SHA256SUMS`. The builder never downloads anything. It refuses a vendored file
whose sha256 differs, records both paths and sha256 in the build report, the log and (default
mode) the file's root attributes (`iers_table`, `iers_table_sha256`, `iers_b_table`,
`iers_b_table_sha256`, `iers_table_mjd_range`, `iers_table_warning`), and logs a `Warning` and sets
`iers_table_warning` if any exposure time lies outside the table or in its predicted part (MJD >=
61280 for the pinned table; the 0.4.0 exposures end at 61230.43). A newer release window needs a
newer table: pass `--iers/--iers-b` (or update the vendored pair and `SHA256SUMS`).

## Speed (0.4.0 inputs, 26,715,540 rows, head node, 8 threads + 8 workers)

Both builds ran on the head node, 2026-11-24, at 8 processes. The load average was about 9-10 from other agents' jobs. The Python harvest also overlapped one of my Julia test builds, and the Julia cold build overlapped the allTargets regression. Treat differences under about 10 % as noise.

| stage | 0.4.0 Python chain, re-run today | Julia builder, cold cache |
|---|---|---|
| row identity | 48.5 s (`01_build_identity.jl`) | 14.8 s |
| mid-exposure times (144,381 files) | 1,497 s (process pool, 8) | 1,370 s (Distributed, 8 workers) |
| almanac join, barycentric, velocities, write (`s6`) | 482 s (join 232 s, single process) | 264 s: almanac 132 s (parallel read 55 s, ordered apply 76 s), wait for astropy kernel 21 s, workup read 30 s, row stages 5 s, write 76 s |
| `bary_flag` (`s9`) | 2.7 s | 0.4 s |
| Galactic + `qa_bitmask` (`s10`) | 18.7 s | 2.5 s |
| `v_rad_corr`, `offtarget` (0.4.0 copy-and-swap steps) | not re-run (0.4.0 log: 38 s for `v_rad_corr`; offtarget steps not logged) | 0.4 s, and inside the join |
| **total** | **2,049 s** (without the copy-and-swap steps) | **1,649 s**, all 52 columns |
| peak memory | 9.2 GB (`s6`) | 9.3 GB head + 8 x ~0.5 GB workers |
| pass B rebuild (cache warm) | n/a | 284 s |

For reference, the 0.4.0 logs show s1 91 s, s4 998 s with 12 processes, and s6 549 s.

- The harvest is bound by Ceph open latency (a few % CPU per process) in both implementations. It speeds up with more processes, not more cores, so `--workers` can exceed the core count there.
- Without the harvest, Julia is about 2x faster (279 s against 552 s).
- The serial part left is the ordered apply of the almanac join, 76 s. It can be split by telescope, or the categorical codes can be built in a second pass, if it ever matters.

**0.4.1 default (2026-11-24, warm cache, 6 threads + 4 workers, head node load ~2):** 167 s total
(almanac join 89 s, barycentric basis in Julia 3.5 s, row stages 4 s, write 67 s), 9.5 GB peak; two runs
170 s and 167 s. The astropy kernel it replaces took 133-148 s (it ran concurrently with the join, so
the saving shows up as the removed "wait for astropy kernel" and freed cores). `--numerics legacy040`:
165 s.

Per-stage wall times of any build are in its `*.report.json` (`timing_s`).

## 0.4.1 content changes (not implemented; separate decisions, file:line in the PR description)

`exposure_flags` / `exposure_state` from arMADGICS PR #45, synthetic `plate_id >= 1e9` with
`original_plate_id`, sdss_id 0 handling, `sdss5_target_flags` width (allTargets), the excluded-exposure
cut in the recommended cuts.
