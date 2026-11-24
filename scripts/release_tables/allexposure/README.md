# AllExposure: the single-pass allExposure builder (Julia)

```bash
julia --project=scripts/release_tables/allexposure -t 8 scripts/release_tables/allexposure/bin/build_allexposure.jl \
    --full-list <raw run>/full_list_info.h5 --almanac <allobs_*.h5> --apred <apred root> \
    --workup <workup dir> --exclusion-lists <dir with bad_nights/bad_fibers/rv_railed.csv> \
    [--fiber-offsets fiber_offsets.csv] [--iers-a finals2000A.all] \
    --cache-dir <scratch> --workers 8 --python <venv>/bin/python \
    --out allExposure_<name>.h5 --report allExposure_<name>.report.json
julia --project=scripts/release_tables/allexposure scripts/release_tables/allexposure/test/runtests.jl
```

One row per arMADGICS workup row, in `full_list_info.h5` order. It replaces the 0.4.0 chain
`s1` (identity) + `s4` (times) + `s6` (build) + `s9` (bary_flag, in place) + `s10` (Galactic + QA, in
place) + `fiber_offsets/s10_final_build.jl` (copy-and-swap) + `offtarget_sidecar/s1,s2,s4` (copy-and-swap)
with one program that writes the final table once. It reproduces the released 0.4.0 table bit for bit
(values, dtypes, chunking, filters, attributes and their HDF5 string types), except for the
IERS-dependent rows explained in `../README.md`.

## Layout

| file | contents |
|---|---|
| `src/schema.jl` | **the schema**: every column's name, element type, producing stage, on-disk layout class and attributes; the root attributes. The writer and the tests iterate it; nothing else lists columns |
| `src/identity.jl` | stage IDENTITY (full_list_info.h5) and TIME (mid-exposure-time harvest, Distributed); caches in the 0.4.0 file formats |
| `src/almanac.jl` | stage ALMANAC and the almanac half of OFFTARGET: parallel per-night reads, ordered apply |
| `src/derived.jl` | row-wise stages: BARY, VELOCITY, BARYFLAG, GALACTIC, QA, FIBEROFFSET; the `Numerics` types |
| `src/build.jl` | `BuildConfig`, `build` (stage orchestration, timing, report), `write_table` |
| `src/h5util.jl` | layout-faithful dataset writes; fixed vs variable-length string attributes |
| `python/legacy_kernels.py` | the two numpy/astropy kernels (see Numerics) |
| `bin/build_allexposure.jl` | CLI |
| `test/runtests.jl` | unit tests (+ a full-table bitwise regression when `ALLEXPOSURE_REFERENCE`/`ALLEXPOSURE_CANDIDATE` are set); `test/fixture_rows_040.csv` = 400 released rows as exact hex floats |
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

## Numerics (what bit-identity with 0.4.0 required)

`--numerics legacy` (default) reproduces 0.4.0; `--numerics native` is pure Julia and CPU-independent.

| quantity | 0.4.0 implementation | legacy mode here | native mode |
|---|---|---|---|
| `cos`/`sin` of the target direction | numpy float64 (= glibc libm, verified) | `ccall` to libm | `Base.cos/sin` |
| `10^(RV_pixoff_final*6e-6)` | numpy `**`: Intel SVML on AVX-512 (1-ulp differences from correctly rounded pow in ~7 % of rows) | taken from `legacy_kernels.py pow10` (numpy) | `Base.exp10` |
| barycentric basis (6 directions per exposure) and `U = pinv(A) @ Y` | astropy `radial_velocity_correction` + numpy/BLAS | `legacy_kernels.py bary-basis` (same code) | same (astropy) |
| `glon`, `glat` | astropy ICRS -> FK5(J2000) -> Galactic: ERFA `s2c`, `rxp`, `c2s`, astropy unit scalings and longitude wrap | ERFA's C library (ERFA_jll) with the same matrices (exact hex floats), scalings and wrap | same |
| everything else | numpy element-wise arithmetic | same operation order in Julia (no FMA contraction, `x^2` = `x*x`, IEEE `sqrt` giving NaN for negative input) | same |

The barycentric basis stays in Python in both modes: it is astropy's model (ERFA plus astropy's
IERS interpolation) by definition, and costs ~150 s for 144k exposures. It runs concurrently with the
almanac join. Both kernels cache their output in `--cache-dir`. Native mode is not bit-identical to
0.4.0. Measured on all 26.7 M rows against the legacy build: `v_barycentric_correction` differs in 98
rows (<= 6.7e-11 km/s), `v_rel` in 1,789,469 rows (<= 1.0e-10 km/s), `v_rad` in 1,681,546 rows
(<= 1.3e-10 km/s), `e_v_rad` in 1,256,229 rows (<= 8 ulp). All far below any physical meaning;
consider native mode for 0.4.1 if CPU independence matters more than continuity with 0.4.0 bits.

## Speed (0.4.0 inputs, 26,715,540 rows, head node, 8 threads + 8 workers)

See the timing table in the PR / build reports (`*.report.json` has per-stage wall times).

## 0.4.1 changes (not implemented; file:line in the PR description)

`exposure_flags` / `exposure_state` from arMADGICS PR #45, synthetic `plate_id >= 1e9` with
`original_plate_id`, sdss_id 0 handling, `sdss5_target_flags` width (allTargets), the excluded-exposure
cut in the recommended cuts.
