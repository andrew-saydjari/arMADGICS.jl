# Release summary tables: allExposure, allVisitRV, allTargets

The scripts that build the APOGEE release summary tables and their side products
(`data/work/mwm/<release>/public/`): `allExposure_<name>.h5`, `allVisitRV_<name>.h5`,
`allTargets_<name>.h5`, `bad_nights.csv`, `bad_fibers.csv`, `rv_railed.csv`, `fiber_offsets.csv`,
`FIBER_OFFSETS_TABLE.md`, `offtarget_exposures.csv`.

**The tables are products. Never commit them.** Every script writes to a path given as an argument,
and `.gitignore` here refuses h5/fits/csv/npy/parquet/jls outputs. The only committed data are three
small fixtures/reference inputs (< 100 kB each, whitelisted in `.gitignore`).

These scripts were moved here from the dated work folders that built 0.4.0. Except for allExposure,
the logic is unchanged: only hard-coded paths and release names became arguments, scratch code was
removed, and each script got a header (purpose, inputs, outputs, place in the chain, 0.4.0 source).
**allExposure was ported to Julia** (`allexposure/`, a single-pass builder; see
[allexposure/README.md](allexposure/README.md)); the 0.4.0 Python scripts are kept in
`allexposure/legacy_python/` as the reference implementation for the bit-identity regression.

## The chain

Inputs for a release: the almanac file the workup used, the arMADGICS raw run (`full_list_info.h5`)
and workup outputs (`arMADGICS_out_*.h5`), the apred 1D outputs (mid-exposure times), the wavecal /
trace census CSVs behind the exclusion lists, and the target-metadata file written by
`almanac add metadata`.

| # | step | script | inputs | outputs | lang |
|---|---|---|---|---|---|
| 1 | row identity | `allexposure/legacy_python/01_build_identity.jl` (also done inside the builder) | full_list_info.h5 | identity.h5 | Julia |
| 2 | exclusion lists | `exclusions/build_exclusion_lists.jl` | wavecal census CSVs, trace-jump CSV, identity.h5, workup (RV_pixoff_final, RV_pixoff_disc_final, snr, RV_flag, ingestBit) | bad_nights.csv, bad_fibers.csv, rv_railed.csv | Julia |
| 3 | allExposure, pass A | `allexposure/bin/build_allexposure.jl` (no `--fiber-offsets`) | full_list_info.h5, almanac, apred, workup, exclusion lists | allExposure without the 3 v_rad_corr columns | Julia (+2 Python numerics kernels) |
| 4 | visit sequences | `allvisitrv/01_seqmap.jl` | almanac (`sequences/objects`) | seqmap.jls | Julia |
| 5 | visits | `allvisitrv/02_build_visits.jl` | allExposure (pass A), seqmap.jls | v2_cols.jls | Julia |
| 6 | write allVisitRV | `allvisitrv/03_write_allvisitrv.jl` | v2_cols.jls | allVisitRV (21 datasets) | Julia |
| 7 | allVisitRV attributes | `allvisitrv/04_attrs_allvisitrv.py` | allExposure, allVisitRV | (in place) descriptions, `categories` | Python |
| 8 | fiber offsets: working tables | `fiber_offsets/01_build_working_tables.jl` | allExposure (pass A), allVisitRV (21), seqmap.jls | build.h5 | Julia |
| 9 | fiber offsets: solve | `fiber_offsets/02_solve_offsets.jl` | build.h5 | offsets_raw.csv, solve.h5 | Julia |
| 10 | fiber offsets: table | `fiber_offsets/03_offset_table.jl` | offsets_raw.csv, bad_fibers.csv | fiber_offsets.csv (review header) | Julia |
| 11 | fiber offsets: release files | `fiber_offsets/04_release_header.jl` | fiber_offsets.csv | **fiber_offsets_bare.csv (= released fiber_offsets.csv)**, FIBER_OFFSETS_TABLE.md | Julia |
| 12 | allVisitRV corrected velocities | `fiber_offsets/05_apply_allvisitrv.jl` | fiber_offsets.csv, build.h5, allVisitRV (21) | allVisitRV (24) | Julia |
| 13 | allExposure, final | `allexposure/bin/build_allexposure.jl --fiber-offsets` | as step 3 + fiber_offsets.csv (cache reused) | **allExposure** (52 datasets) | Julia |
| 14 | offtarget: FPS scan | `offtarget/01_almanac_scan.jl` | almanac | almanac_fps_flags.h5 | Julia |
| 15 | offtarget: sidecar + visits | `offtarget/02_sidecar_and_visits.jl` | allExposure, allVisitRV (24), almanac_fps_flags.h5 | **offtarget_exposures.csv**, offtarget_visit_flag.h5 | Julia |
| 16 | allVisitRV offtarget | `offtarget/03_stage_allvisitrv.py` | allVisitRV (24), offtarget_visit_flag.h5 | **allVisitRV** (25 datasets) | Python |
| 17 | metadata ids | `metadata/make_restrict_ids.py` | allExposure | sdss_id .npy | Python |
| 18 | target metadata | `metadata/build_metadata.sh` (`almanac add metadata`, needs the SDSS database) | almanac copy, ids | decorated almanac (`/meta`, ...) | shell + almanac |
| 19 | allTargets | `alltargets/build_allTargets.jl` (0.4.1 version: passes through every `/meta` column) | allExposure, decorated almanac, `alltargets/reference/*` | **allTargets** | Julia |
| 20 | checks | `cuts/recommended_cuts.py` (or `.jl`), `validation/compare_h5.py` | products | counts / bitwise comparison | Python, Julia |

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
       META_FILE=<decorated almanac> PYTHON=<venv>/bin/python IERS_A=<pinned finals2000A.all>
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
- **Python** (3.12): `requirements.txt` (h5py, numpy, astropy, pyerfa, astropy-iers-data pinned to the
  environment that reproduced 0.4.0). With uv: `uv venv --python 3.12 && uv pip install -r requirements.txt`.
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
| allExposure (Julia builder) | all 52 datasets x 26,715,540 rows, all attributes, layouts | identical except 3 columns in 1,387 rows of the last 7 nights; see below |
| allExposure (legacy Python, run today) | its 48 datasets x all rows | same 1,387-row difference, and bit-identical to the Julia build on all 48 shared datasets |
| allExposure, Julia cold vs warm cache | all 52 datasets | identical; the Julia time harvest equals the legacy exposure_times.h5 value for value |
| exposure_times.h5 (legacy step 2) | full file | byte-identical |
| seqmap.jls (step 4) | full file | byte-identical |
| allVisitRV (steps 5-7, 12, 16) | all 25 datasets x 6,670,336 rows, attributes, layouts | identical |
| fiber_offsets.csv, FIBER_OFFSETS_TABLE.md (steps 8-11) | full files | byte-identical |
| offtarget_exposures.csv (step 15) | full file | byte-identical |
| restrict-ids .npy (step 17) | full file | byte-identical |
| allTargets (step 19, 0.4.1 builder, `--release 0.4.0`, 0.4.0 title, 0.4.0 `/meta`) | all 57 released datasets (44 columns + `carton_pks/` + `carton_table/`), dtype, shape, layout, values, dataset attributes | 57 / 57 identical; differs only in provenance file attributes (builder path and sha256, creation time, mapping/carton-table paths) and in the 48 extra `/meta` columns it now passes through (by design) |
| recommended cuts | full table | 25,516,821 (5 cuts) and 25,207,025 (6 cuts) rows kept, as documented |

**The one difference: the barycentric correction of the last week.** `v_barycentric_correction`
differs in 1,387 rows (27 exposures, MJD 61224-61230), by at most 1.0e-10 km/s (149,897 ulp), and it
carries into `v_rad` (1,178 rows) and `v_rad_corr` (1,150 rows). The 0.4.0 Python chain run today
gives exactly the same difference, so it is not the port. Cause: astropy's IERS Earth-orientation
table. The barycentric correction depends on UT1-UTC and polar motion, and astropy picks its table
at run time (`IERS_Auto`). The table astropy uses today (the bundled `astropy-iers-data`
2026.8.31) differs, for those July 2026 dates, from the one in use when 0.4.0 was built on
2026-09-13. Neither the table nor its version was recorded. A freshly downloaded IERS-A table differs
on more nights (61222-61230), so the 0.4.0 table cannot be reconstructed from today's sources. The
Julia builder now records the table it used (path and sha256, in the build report) and can pin one
(`--iers-a` / `IERS_A`). It should be pinned for 0.4.1.

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
  about 7 % of rows). The 0.4.0 Python chain would not reproduce 0.4.0 on a rome node.
- The barycentric correction depends on an unrecorded IERS table (above).
- `join_status` 5 (fiber index out of range) cannot occur: the almanac index is computed from
  `adjfiberindx` and is always 1-300.
- `sdss_id_agrees` treats `sdss_id <= 0` as untestable (3), so the solar-system target (sdss_id 0,
  Vesta) is never cross-checked; `make_restrict_ids.py` and allVisitRV drop `sdss_id <= 0`.
- allVisitRV and the fiber-offset fit use the 5-cut mask; the released `RECOMMENDED_CUTS.md` adds a
  sixth cut (`offtarget == 0`) that they predate.
- The fiber-offset fit, the allVisitRV mask and `v_rad_corr_flag` all hard-code the 5-cut mask
  independently (four copies).

## Not in this directory

- The analyses behind the exclusion lists (wavecal census, trace-jump census, chi2 pathology,
  RV-grid rail) are research products in dated folders; their CSVs are inputs here.
- The 0.4.0 verification scripts that compared each staged file against the previous release
  (`s7`, `s11`, `s12`, `s5_verify_staged`, `s11_final_verify`, ...): they were tied to one copy-and-swap
  transition. `validation/compare_h5.py` replaces them for any rebuild.
- Diagnostic and plotting scripts (`s1-s3`, `s5`, `s7`, `s8` of allexposure; `s1_diag`, `s3-s5` of
  allvisitrv; `s5-s8`, `s5b` of fiber_offsets; `s6` of offtarget).
- `targetdb_carton.csv` was dumped once from targetdb without a kept script (see alltargets/README.md).
- `RECOMMENDED_CUTS.md`, `allVisitRV_README.md` and the release `README` prose are hand-written
  documents; `allvisitrv/allVisitRV_README.md` is the 0.4.0 text, kept for reference.
