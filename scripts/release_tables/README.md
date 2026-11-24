# Release summary tables: allExposure, allVisitRV, allTargets

The scripts that build the APOGEE release summary tables and their side products
(`data/work/mwm/<release>/public/`): `allExposure_<name>.h5`, `allVisitRV_<name>.h5`,
`allTargets_<name>.h5`, `bad_nights.csv`, `bad_fibers.csv`, `rv_railed.csv`, `fiber_offsets.csv`,
`FIBER_OFFSETS_TABLE.md`, `offtarget_exposures.csv`.

**The tables are products. Never commit them.** Every script writes to a path given as an argument,
and `.gitignore` here refuses h5/fits/csv/npy/parquet/jls outputs. The only committed data are inputs:
four small fixtures/reference inputs (< 100 kB each, whitelisted in `.gitignore`) and the two pinned
IERS Earth-orientation tables in `allexposure/iers/` (3.8 + 5.2 MB of text).

These scripts were moved here from the dated work folders that built 0.4.0. Except for allExposure,
the logic is unchanged: only hard-coded paths and release names became arguments, scratch code was
removed, and each script got a header (purpose, inputs, outputs, place in the chain, 0.4.0 source).
**allExposure was ported to Julia** (`allexposure/`, a single-pass builder; see
[allexposure/README.md](allexposure/README.md)), including the barycentric correction, and so were the
small Python helpers. **A default run of the chain uses no Python.** The 0.4.0 Python scripts are kept,
frozen, as regression references: `allexposure/legacy_python/` (allExposure chain) and
`legacy_python/` (helpers).

## 0.4.1 defaults vs. legacy-for-regression

| | default (0.4.1) | legacy, for the 0.4.0 bit-regression only |
|---|---|---|
| `10^x` in `v_rel`, `v_rad`, `e_v_rad` | Julia `Base.exp10`: same bits on every CPU type | `--numerics legacy040` / `NUMERICS=legacy040`: numpy `10**x` via `allexposure/python/legacy_kernels.py pow10`; reproduces 0.4.0 only on AVX-512 nodes |
| barycentric basis | Julia port of astropy 8.0.1 on ERFA_jll (`allexposure/src/bary.jl`), bit-identical to astropy | same (also `--bary astropy`: the astropy kernel, as a cross-check) |
| IERS tables | pinned, vendored: `allexposure/iers/finals2000A.all` + `eopc04.1962-now` (astropy-iers-data 0.2026.9.7.0.56.14), sha256-checked, recorded in the file and report; never downloaded | same |
| allVisitRV attributes, offtarget staging, restrict ids, cuts, comparison | Julia (`*.jl`) | frozen Python in `legacy_python/` |
| root attributes of allExposure | 0.4.0 set, `bary_method`/`astropy_version` updated, plus `erfa_version`, `iers_*`, `numerics` | exactly the 0.4.0 set |

**Python that remains:** none in a default run. Python is used only by (1) `--numerics legacy040`
(numpy `10**x`, to reproduce the 0.4.0 bits), (2) the optional `--bary astropy` cross-check, and (3)
the frozen reference scripts in `legacy_python/` and `allexposure/legacy_python/`, plus the database
step `metadata/build_metadata.sh`, which runs the `almanac` package (a Python tool outside this
repository).

**Default vs 0.4.0** (allExposure, all 26,715,540 rows, 2026-11-24): `v_barycentric_correction`
identical; the `10^x` change alone moves

| column | rows differing | max abs difference | max difference in ulp |
|---|---|---|---|
| `v_rel` | 1,789,469 (6.70 %) | 9.98e-11 km/s | 1 ulp of `10^x`, carried |
| `v_rad` | 1,681,470 (6.29 %) | 1.33e-10 km/s | (as above) |
| `e_v_rad` | 1,256,229 (4.70 %) | 2.91e-11 km/s (max relative 9.7e-16) | 8 |
| `v_rad_corr` | 1,603,833 (6.00 %) | 1.33e-10 km/s | (as above) |

The max relative difference of `v_rel`/`v_rad` (up to 1.9e-6) occurs where the velocity is near zero;
the absolute differences are the meaningful measure. No finite/NaN pattern changes, so every flag,
mask and count is unchanged.

## The chain

Inputs for a release: the almanac file the workup used, the arMADGICS raw run (`full_list_info.h5`)
and workup outputs (`arMADGICS_out_*.h5`), the apred 1D outputs (mid-exposure times), the wavecal /
trace census CSVs behind the exclusion lists, and the target-metadata file written by
`almanac add metadata`.

| # | step | script | inputs | outputs | lang |
|---|---|---|---|---|---|
| 1 | row identity | `allexposure/legacy_python/01_build_identity.jl` (also done inside the builder) | full_list_info.h5 | identity.h5 | Julia |
| 2 | exclusion lists | `exclusions/build_exclusion_lists.jl` | wavecal census CSVs, trace-jump CSV, identity.h5, workup (RV_pixoff_final, RV_pixoff_disc_final, snr, RV_flag, ingestBit) | bad_nights.csv, bad_fibers.csv, rv_railed.csv | Julia |
| 3 | allExposure, pass A | `allexposure/bin/build_allexposure.jl` (no `--fiber-offsets`) | full_list_info.h5, almanac, apred, workup, exclusion lists, pinned IERS table | allExposure without the 3 v_rad_corr columns | Julia |
| 4 | visit sequences | `allvisitrv/01_seqmap.jl` | almanac (`sequences/objects`) | seqmap.jls | Julia |
| 5 | visits | `allvisitrv/02_build_visits.jl` | allExposure (pass A), seqmap.jls | v2_cols.jls | Julia |
| 6 | write allVisitRV | `allvisitrv/03_write_allvisitrv.jl` | v2_cols.jls | allVisitRV (21 datasets) | Julia |
| 7 | allVisitRV attributes | `allvisitrv/04_attrs_allvisitrv.jl` | allExposure, allVisitRV | (in place) descriptions, `categories` | Julia |
| 8 | fiber offsets: working tables | `fiber_offsets/01_build_working_tables.jl` | allExposure (pass A), allVisitRV (21), seqmap.jls | build.h5 | Julia |
| 9 | fiber offsets: solve | `fiber_offsets/02_solve_offsets.jl` | build.h5 | offsets_raw.csv, solve.h5 | Julia |
| 10 | fiber offsets: table | `fiber_offsets/03_offset_table.jl` | offsets_raw.csv, bad_fibers.csv | fiber_offsets.csv (review header) | Julia |
| 11 | fiber offsets: release files | `fiber_offsets/04_release_header.jl` | fiber_offsets.csv | **fiber_offsets_bare.csv (= released fiber_offsets.csv)**, FIBER_OFFSETS_TABLE.md | Julia |
| 12 | allVisitRV corrected velocities | `fiber_offsets/05_apply_allvisitrv.jl` | fiber_offsets.csv, build.h5, allVisitRV (21) | allVisitRV (24) | Julia |
| 13 | allExposure, final | `allexposure/bin/build_allexposure.jl --fiber-offsets` | as step 3 + fiber_offsets.csv (cache reused) | **allExposure** (52 datasets) | Julia |
| 14 | offtarget: FPS scan | `offtarget/01_almanac_scan.jl` | almanac | almanac_fps_flags.h5 | Julia |
| 15 | offtarget: sidecar + visits | `offtarget/02_sidecar_and_visits.jl` | allExposure, allVisitRV (24), almanac_fps_flags.h5 | **offtarget_exposures.csv**, offtarget_visit_flag.h5 | Julia |
| 16 | allVisitRV offtarget | `offtarget/03_stage_allvisitrv.jl` | allVisitRV (24), offtarget_visit_flag.h5 | **allVisitRV** (25 datasets) | Julia |
| 17 | metadata ids | `metadata/make_restrict_ids.jl` | allExposure | sdss_id .npy | Julia |
| 18 | target metadata | `metadata/build_metadata.sh` (`almanac add metadata`, needs the SDSS database) | almanac copy, ids | decorated almanac (`/meta`, ...) | shell + almanac |
| 19 | allTargets | `alltargets/build_allTargets.jl` (0.4.1 version: passes through every `/meta` column) | allExposure, decorated almanac, `alltargets/reference/*` | **allTargets** | Julia |
| 20 | checks | `cuts/recommended_cuts.jl`, `validation/compare_h5.jl` | products | counts / bitwise comparison | Julia |

Why two allExposure passes: the fiber zero points are fitted on allVisitRV, which is built from
allExposure. Pass A builds everything except `v_rad_corr*`; pass B reuses the cached identity, exposure
times and numerics-kernel outputs and only re-joins and re-writes (about 4-5 minutes). In 0.4.0 the
same dependency was handled by appending columns to the released file (copy-and-swap).

`run_release_tables.sh` runs steps 1-17, 19, 20 with the release name and every path as environment
variables (see its header). Step 18 needs the database and is run separately
(`metadata/build_metadata.sh`; the 0.4.0 run record is in its header).

```bash
export RELEASE=0.4.1 NAME=pass1_dr21 OUTDIR=<products dir> WORKDIR=<scratch dir> \
       FULL_LIST=<raw run>/full_list_info.h5 ALMANAC=<allobs_*.h5> APRED=<apred root> WORKUP=<workup dir> \
       WAVECAL_DIR=<wavecal census dir> TRACE_CSV=<trace_night_jump_outliers.csv> \
       META_FILE=<decorated almanac>     # IERS defaults to the vendored allexposure/iers/finals2000A.all
scripts/release_tables/run_release_tables.sh
```

## Environments

- **Julia** (1.11, `julia --project=...`, never the global environment). Three projects:
  - `scripts/release_tables/Project.toml`: every Julia script except allExposure and allTargets
    (HDF5, CSV, DataFrames, ArgParse, JSON, Jackknife, ...). **CSV.jl is pinned to 0.10.17**, the
    version that wrote the 0.4.0 `fiber_offsets.csv`: CSV.jl 1.x quotes empty strings (`""`), which
    changes the released file's bytes. A separate small project, not the repo's root project,
    because the root project pulls in the whole pipeline (Korg, MPI, ...) and these scripts need
    none of it.
  - `allexposure/Project.toml`: the `AllExposure` package (HDF5, CSV, JSON, ArgParse, ERFA_jll, Distributed).
  - `alltargets/Project.toml`: the allTargets builder (as developed; HDF5 0.17.4, ArgParse, JSON, SHA).
  Instantiate each once: `julia --project=<dir> -e 'using Pkg; Pkg.instantiate()'`.
- **Python** (3.12) is **not needed** for a default (0.4.1) run. It is needed only for the 0.4.0
  bit-regression (`NUMERICS=legacy040`, the frozen `legacy_python/` scripts) and the optional astropy
  cross-check (`BARY=astropy`): `requirements.txt` (h5py, numpy, astropy, pyerfa, astropy-iers-data pinned
  to the environment that reproduced 0.4.0). With uv: `uv venv --python 3.12 && uv pip install -r requirements.txt`.
- **Environment variables on Ceph** (set by the driver, as `workup/run_workup.sh` does for the
  workup): `CEPHTWEAKS_LAZYIO=1`, `HDF5_USE_FILE_LOCKING=FALSE`.

## Regression against 0.4.0 (2026-11-24)

Every step was run from this directory on the 0.4.0 inputs and compared with the released files in
`data/work/mwm/0.4.0/public/` (bitwise for HDF5 values, plus dtype, shape, chunking, filters and every
attribute including its HDF5 string type; byte-for-byte for CSV/MD). To test each step on its own,
the visit, fiber-offset and offtarget steps were fed the released upstream tables (the inputs they
had in 0.4.0).

| product | coverage | result |
|---|---|---|
| identity.h5 (step 1) | all 26,715,540 rows | identical |
| bad_nights.csv, bad_fibers.csv, rv_railed.csv (step 2) | full files | byte-identical |
| allExposure (Julia builder, `--numerics legacy040`, pinned IERS tables) | all 52 datasets x 26,715,540 rows, all attributes, layouts | **identical** (ALL EQUAL) |
| allExposure (Julia builder, default) | all 52 datasets | `v_rel`, `v_rad`, `e_v_rad`, `v_rad_corr` differ as tabulated above; new root attributes; everything else identical |
| allExposure, default build run twice | all 52 datasets and attributes | identical (deterministic; the files differ only in HDF5 header timestamps) |
| allExposure (legacy Python, run with today's bundled IERS table) | its 48 datasets x all rows | bit-identical to the Julia build on all 48 shared datasets except 1,387 rows of `v_barycentric_correction` (MJD 61224-61230), the IERS table (below) |
| allExposure, Julia cold vs warm cache | all 52 datasets | identical; the Julia time harvest equals the legacy exposure_times.h5 value for value |
| exposure_times.h5 (legacy step 2) | full file | byte-identical |
| seqmap.jls (step 4) | full file | byte-identical |
| allVisitRV (steps 5-7, 12, 16) | all 25 datasets x 6,670,336 rows, attributes, layouts | identical |
| fiber_offsets.csv, FIBER_OFFSETS_TABLE.md (steps 8-11) | full files | byte-identical |
| offtarget_exposures.csv (step 15) | full file | byte-identical |
| restrict-ids .npy (step 17) | full file | byte-identical |
| allTargets (step 19, 0.4.1 builder, `--release 0.4.0`, 0.4.0 title, 0.4.0 `/meta`) | all 57 released datasets (44 columns + `carton_pks/` + `carton_table/`), dtype, shape, layout, values, dataset attributes | 57 / 57 identical; differs only in provenance file attributes (builder path and sha256, creation time, mapping/carton-table paths) and in the 48 extra `/meta` columns it now passes through (by design) |
| recommended cuts | full table | 25,516,821 (5 cuts) and 25,207,025 (6 cuts) rows kept, as documented |

**The IERS table, found and pinned.** The first version of this PR reported one difference:
`v_barycentric_correction` in 1,387 rows (27 exposures, MJD 61224-61230, at most 1.0e-10 km/s), from
astropy's run-time choice of IERS Earth-orientation table, which 0.4.0 did not record. Testing the
weekly `astropy-iers-data` releases around the 0.4.0 build date (2026-09-13) against all 26.7 M rows,
the tables of 0.2026.9.7.0.56.14, combined as astropy's `IERS_Auto` does (IERS-A
`finals2000A.all` with the IERS-B `eopc04.1962-now` values substituted), reproduce every bit. They are
vendored in `allexposure/iers/` and pinned (details and the tables that do not match:
`allexposure/iers/fetch_iers.sh`).

## Things 0.4.0 does that look like bugs (preserved for bit-identity; not fixed)

- allVisitRV's root attribute `provenance` says "STAGED, NOT RELEASED. Columns 1-21 are
  byte-identical copies of <scifest path>..." although the file is released
  (`fiber_offsets/05_apply_allvisitrv.jl`, kept via its 4th argument).
- allExposure `offtarget` description and allVisitRV `offtarget` description cite an internal task
  number ("Task 94") and the pass-1 almanac path.
- allExposure root attributes are a mix of fixed-length and variable-length strings (a side effect of
  the Julia copy-and-swap that added `v_rad_corr`); reproduced exactly (schema `FixedStr` / `VlenStr`).
- `fiber_offset_convention` (root attribute) was added to the released allExposure by a step with no
  surviving script; its text is in `allexposure/src/schema.jl`.
- 0.4.0's allExposure `v_rel`, `v_rad`, `e_v_rad` bits depend on the CPU: numpy evaluates
  `10 ** x` with Intel SVML on AVX-512 machines (1-ulp differences from a correctly rounded pow in
  about 7 % of rows). The 0.4.0 Python chain would not reproduce 0.4.0 on a rome node. **Fixed for
  0.4.1** (Julia `exp10`); kept behind `--numerics legacy040` for the regression.
- 0.4.0 used an unrecorded IERS table. **Fixed for 0.4.1**: identified and pinned (above).
- `join_status` 5 (fiber index out of range) cannot occur: the almanac index is computed from
  `adjfiberindx` and is always 1-300.
- `sdss_id_agrees` treats `sdss_id <= 0` as untestable (3), so the solar-system target (sdss_id 0,
  Vesta) is never cross-checked; `make_restrict_ids.jl` and allVisitRV drop `sdss_id <= 0`.
- allVisitRV and the fiber-offset fit use the 5-cut mask; the released `RECOMMENDED_CUTS.md` adds a
  sixth cut (`offtarget == 0`) that they predate.
- The fiber-offset fit, the allVisitRV mask and `v_rad_corr_flag` all hard-code the 5-cut mask
  independently (four copies).

## Not in this directory

- The analyses behind the exclusion lists (wavecal census, trace-jump census, chi2 pathology,
  RV-grid rail) are research products in dated folders; their CSVs are inputs here.
- The 0.4.0 verification scripts that compared each staged file against the previous release
  (`s7`, `s11`, `s12`, `s5_verify_staged`, `s11_final_verify`, ...): they were tied to one copy-and-swap
  transition. `validation/compare_h5.jl` replaces them for any rebuild.
- Diagnostic and plotting scripts (`s1-s3`, `s5`, `s7`, `s8` of allexposure; `s1_diag`, `s3-s5` of
  allvisitrv; `s5-s8`, `s5b` of fiber_offsets; `s6` of offtarget).
- `targetdb_carton.csv` was dumped once from targetdb without a kept script (see alltargets/README.md).
- `RECOMMENDED_CUTS.md`, `allVisitRV_README.md` and the release `README` prose are hand-written
  documents; `allvisitrv/allVisitRV_README.md` is the 0.4.0 text, kept for reference.
