# allTargets builder for 0.4.1 (all `/meta` columns, no column list)

2026-11-24. The 0.4.1 version of `2026_11_19/carton_flags/build_allTargets.jl`, the builder of
the 0.4.0 product `data/work/mwm/0.4.0/public/allTargets_pass1_dr21.h5`. The 0.4.0 script was
copied, not edited. Developed in a dated work folder (2026_11_24/alltargets_041) and moved into
this repository with its tests and tools; the `test_*` output folders stay in the work folder (they
are products) and are only cited below.

**Chain position:** last table of the release chain (../README.md). Inputs: the release's final
allExposure (for the sdss_id set and n_exposures), the targets file written by
`almanac add metadata` (../metadata/), and the two small reference inputs in `reference/`.

The request (2026-11-24): "I want all targets to just get all of the metadata columns by default,
I don't want to have to update it everytime more metadata columns are added".

## Run it for 0.4.1

```bash
cd scripts/release_tables/alltargets
C=reference                     # carton table + semaphore mapping used for 0.4.0 (committed fixtures)
julia --project=. -t 1 build_allTargets.jl --release 0.4.1 \
   --allexposure <0.4.1 allExposure .h5> \
   --meta        <0.4.1 targets file from `almanac add metadata --output ...` (default or --extended)> \
   --carton-table      $C/targetdb_carton.csv \
   --semaphore-mapping $C/sdss5_target_4_with_groups.csv \
   --semaphore-version 0.2.9 \
   --out <output dir>/<allTargets name>.h5
```

- Every input is an argument. Nothing points at 0.4.0 except the carton table and mapping
  above, and those can be replaced (see "Needs attention").
- Optional: `--title`, `--allexposure-label` (name used in the `n_exposures` description; by
  default the allExposure file name without `.h5`), `--name-fallback false` (exact-pk flags only,
  astra's rule), `--include-solar-system false`, `--keep-meta-target-flags true`, and
  `--checkpoint DIR` instead of `--meta` (partial preview from a running almanac build; carried
  over from 0.4.0 and not tested here).
- Cost (0.4.0 inputs, 4.02 M rows): see "Tests". It is one process on one core, so the head node
  is fine for it.
- The script writes `<out>.tmp-<pid>` and renames it to `<out>` when it finishes, so a crash
  never leaves a half-written `<out>`.

## What changed from 0.4.0

| | 0.4.0 builder | this builder |
|---|---|---|
| `/meta` columns read | the 41 in a hard-coded `COLUMNS` list | **every dataset in `/meta`**, found when the file is read |
| missing 0.4.0 column in `/meta` | warning, column left out | **error** (0.4.1 must hold every 0.4.0 column) |
| output name clash (a `/meta` column named `n_cartons`, `carton_pks`, ...) | n/a | **error** |
| `sdss5_target_flags` | rebuilt from the carton pairs | the same. The `/meta` column of that name is not copied; its width goes in attribute `meta_sdss5_target_flags_nbytes`. `--keep-meta-target-flags true` also copies it, as `meta_sdss5_target_flags` |
| rows | sdss_id > 0 in allExposure | the same, plus the solar-system sdss_id (0, Vesta). That row is kept only when allExposure has it and `/meta` marks it with `lead == "solar_system"` |
| HDF5 order | not tracked (alphabetical) | creation order tracked = column order below. h5py lists datasets in this order; h5ls and HDF5.jl list them alphabetically |
| release name in text | "0.4.0" and "allExposure_pass1_dr21" written into the code | `--release`, `--allexposure-label`, `--title` |
| `units` from description | trailing `[...]` | the same, except a trailing abundance label such as `[Fe/H]` is not taken as a unit |
| new file attributes | | `release`, `column_order`, `column_order_rule`, `columns_from_040`, `meta_passthrough_columns`, `n_meta_passthrough_columns`, `columns_added_since_040`, `meta_columns_not_passed_through`, `meta_n_columns`, `meta_sdss5_target_flags_nbytes`, `meta_extended`, `meta_extended_source`, `meta_extended_column_groups`, `multidim_columns`, `solar_system_sdss_id`, `n_solar_system_rows`, `solar_system_rule` |

These are unchanged from 0.4.0:

- the derived columns `n_cartons`, `n_exposures` and `metadata_status` (same bits)
- the carton decoding and `flags_rule`
- the groups `carton_pks/` and `carton_table/`
- the chunking and compression
- the `meta_*` copies of the `/meta` scalar attributes. These include `meta_almanac_version`,
  `meta_almanac_commit`, `meta_input_path`, `meta_extended`, `meta_mode` and `meta_restrict_to`.
  `meta_source` is the path of the input file.

### Column order rule

1. The 44 columns of 0.4.0, in the order of the 0.4.0 product registry
   (`column_order` in `allTargets_pass1_dr21.toml`). The `ORDER_040` list in the script sets the
   order of these columns only. It does not choose which columns are copied.
2. Then every other `/meta` column, in its order in the `/meta` file. That order comes from the
   `/meta` group's `columns` attribute, which is checked against the datasets actually present.
   If the attribute is missing, the HDF5 order of the group is used.
3. Then `meta_sdss5_target_flags`, if requested. The groups `carton_pks/` and `carton_table/`
   come last.

### How each column type is handled

- **Numbers**: the element type is kept. Booleans from h5py (an int8 enum) become `Int8` 0/1,
  as in 0.4.0.
- **Values `/meta` already marks as missing** (NaN, -1, "", 0): copied unchanged.
- **allExposure ids that are not in `/meta`**: filled with NaN for floats, -1 for signed
  integers, 0 for unsigned integers and "" for strings, as in 0.4.0. This did not happen in any
  test.
- **Strings**: `/meta` stores them as fixed-length null-padded strings, and the width changes
  from build to build. For example, `allstar_dr17_synspec_rev1` is 58 bytes in 0.4.0 and 56 in
  the smoke file. The output stores them as variable-length UTF-8, as in 0.4.0. Each new column
  records its original type in the attribute `meta_dtype`.
- **Multi-dimensional columns**: kept as they are. Only the row axis is re-indexed. Each new
  column gets the attribute `shape_per_row`, and the file attribute `multidim_columns` lists
  them all. Today the only one is `sdss5_target_flags`. Its output version is rebuilt and is
  73 bytes wide (mapping v4). `/meta` holds 71 bytes (almanac's v3), and that width is recorded
  in `meta_sdss5_target_flags_nbytes`.
- **Dataset attributes**: new columns get every attribute of the `/meta` dataset (today only
  `description`), plus `units`, `source = "/meta/<name>"` and `meta_dtype`. The 0.4.0 columns
  keep exactly their 0.4.0 attributes.
- **The solar-system row** (sdss_id 0, Vesta): almanac writes it without querying the
  database. It has `lead = "solar_system"`, catalogid and version_id -1, and every catalog value
  NaN, -1 or "". The builder copies it like any other row:
  - `sdss5_target_flags` is all zero, because it has no carton pairs;
  - `n_exposures` counts its allExposure rows;
  - `metadata_status` = 28 (NO_GAIA_DR3 | NO_2MASS | NO_CARTON).
  To find it, use `lead == "solar_system"` or the file attribute `solar_system_sdss_id`.

## Readers: do not assume a fixed column set

The columns now follow the `/meta` input. A default-mode almanac gives 92 columns, an
`--extended` one gives 171, and future almanac versions may add more. Readers should:

- select columns by name;
- use the file attributes `column_order`, `columns_added_since_040` and
  `meta_passthrough_columns` to find out what is present;
- not index columns by position;
- not assume that a column absent from the file is all fill values.

The 44 columns of 0.4.0 are always present: the builder fails if any of them is missing.

## Tests

All test outputs are in this folder. The inputs were opened read-only.

### (a) Regression on the 0.4.0 inputs: every 0.4.0 dataset reproduced exactly

The inputs were:

- `/meta`: `2026_11_18/targets_040/allobs_57618_61230_meta.h5` (default mode, 89 columns, almanac 4e3c218)
- allExposure: the 0.4.0 `allExposure_pass1_dr21.h5`
- the carton table and semaphore mapping from `2026_11_19/carton_flags/`

The run passed `--release 0.4.0` and the 0.4.0 `--title`, so the release text is the same as
0.4.0. Results are in `test_regression_040/`: `build.log`, `compare_to_release.log` and
`verify_passthrough.log`.

- **Build**: 188 s wall, one core, peak RSS 6.3 GB. 4,019,650 rows. The carton pair counts and
  status counts are the same as in the 0.4.0 `rebuild.log`.
- **Comparison with the released file** (`tests/compare_to_release.jl`):
  - **57 of 57 released datasets are identical**: the 44 top-level columns, plus the 4
    `carton_pks/` and 9 `carton_table/` datasets.
  - The comparison covered the HDF5 datatype (H5Tequal), shape, chunking and filters, the
    values (NaN-aware) and every dataset attribute (`description`, `units`). Group attributes
    are also identical.
  - Only three file attributes differ: `builder`, `builder_sha256` and `created_utc`. No
    0.4.0 attribute is missing.
  - 16 file attributes are new (listed above).
- **New columns: 48**. The 0.4.0 `/meta` has 89 columns. 41 of them were already in 0.4.0
  (`sdss_id`, `sdss5_target_flags` and 39 others), so 89 − 41 = 48 are new:
  `allwise, catwise, catwise2020, gaia_ag_gspphot, gaia_azero_gspphot, gaia_distance_gspphot,
  gaia_dr2_source, gaia_ebpminrp_gspphot, gaia_grvs_mag, gaia_grvs_mag_error,
  gaia_grvs_mag_nb_transits, gaia_logg_gspphot, gaia_mh_gspphot, gaia_pm,
  gaia_rv_amplitude_robust, gaia_rv_atm_param_origin, gaia_rv_chisq_pvalue,
  gaia_rv_expected_sig_to_noise, gaia_rv_nb_deblended_transits, gaia_rv_nb_transits,
  gaia_rv_renormalised_gof, gaia_rv_template_fe_h, gaia_rv_template_logg, gaia_rv_template_teff,
  gaia_rv_time_duration, gaia_rv_visibility_periods_used, gaia_rvs_spec_sig_to_noise,
  gaia_source_id, gaia_teff_gspphot, gaia_vbroad, gaia_vbroad_error, gaia_vbroad_nb_transits,
  glimpse, guvcat, panstarrs1, ps1_g18, sdss_dr13_photoobj, sdss_dr17_specobj, skymapper_dr2,
  supercosmos, twomass_h_cmsig, twomass_h_snr, twomass_j_cmsig, twomass_j_snr, twomass_k_cmsig,
  twomass_k_snr, tycho2, unwise`. That makes 92 columns.
- **Pass-through check** (`tests/verify_passthrough.jl`): all 87 pass-through columns are equal
  to `/meta` after joining on sdss_id. For each one the element type, values and description
  match.
- **Size**: the output is 2.6 GB, against 1.5 GB for the 0.4.0 file. This test file can be
  deleted.

### (b) An `--extended` targets file (168 columns, with the Vesta row)

The inputs were:

- `/meta`: `2026_11_24/almanac_041_prep/smoke/targets_extended.h5`. It has 18,001 ids, one of
  them sdss_id 0, `extended = 1`, almanac b215068.
- allExposure: none exists for the smoke data. `tools/make_allexposure_like.jl` built a
  stand-in, `test_extended/allExposure_like_smoke041.h5`, from the smoke almanac file
  `smoke_041.h5`. It has one row per (object exposure × fiber of that exposure's fiber table):
  120,300 rows, the 18,000 sdss_id > 0 and 8 rows with sdss_id 0 (the 8 Vesta exposures at APO
  on MJD 55755).
  - The 0.4.0 allExposure could not be used: it has no sdss_id 0 rows, and it does not cover
    the smoke ids.
  - So `n_exposures` in this test is not what a real allExposure would give.
- the same carton table and semaphore mapping as for (a). Every smoke carton_pk was in the
  table.

Results (`test_extended/`):

- **All 168 `/meta` columns came through**: 166 copied verbatim and checked (0 problems),
  `sdss_id` as the key, and `sdss5_target_flags` rebuilt (71 bytes in `/meta`, 73 in the
  output).
- The output has **171 columns**: the 44 of 0.4.0 plus 127 new ones (the 48 above and the 79
  `--extended` columns).
- The file attributes record the build: `meta_extended = 1`, the 7
  `meta_extended_column_groups`, `meta_almanac_commit = b215068...` and `meta_source`.
- **Vesta row** (sdss_id 0):
  - the only values that are not fill values are `lead = "solar_system"`, `n_exposures = 8` and
    `metadata_status = 28` (NO_GAIA_DR3, NO_2MASS, NO_CARTON);
  - catalogid and version_id are -1, and every catalog value is NaN, -1 or "";
  - it has 0 carton pairs, so its flags are all zero;
  - `solar_system_sdss_id = [0]`, and `row_definition` mentions the row.
- **Column order**: checked with h5py. It is the 44 columns of 0.4.0, then `allwise`, ...,
  `sdss4_apogee_extra_target_flags`, then `carton_pks` and `carton_table`.

### Other checks (`test_misc/`)

- The default smoke file (89 columns) gives 92 output columns.
- `--include-solar-system false` drops sdss_id 0.
- `--keep-meta-target-flags true` adds `meta_sdss5_target_flags` (71 × N, with
  `shape_per_row = [71]`).
- A `/meta` with an added `n_cartons` dataset stops with
  `name collision: /meta column(s) ["n_cartons"] would collide ...`.

### Not tested

- The `--checkpoint` mode.
- A real 0.4.1 allExposure (none exists yet).
- A full-size extended targets file. The memory estimate for one is about 10 GB [I]. The
  builder reads one column at a time.

## Needs attention for 0.4.1 (outside this task; listed, not done)

1. **The allExposure builder** (`2026_09_13/allexposure/s6_build_allexposure.py` and the steps
   after it) is a separate script. It reads the pass-1 almanac fiber tables. It must be rebuilt
   for 0.4.1, including the new arMADGICS `exposure_state` column. It also decides whether the
   Vesta fibers (sdss_id 0, category science) appear in allExposure. allTargets gets a Vesta
   row only if they do. Its `sdss_id_agrees` code 3 treats sdss_id ≤ 0 as untestable, which
   now covers Vesta too.
2. **`targetdb_carton.csv`** was dumped once for 0.4.0. If the 0.4.1 targets have a carton_pk
   that is not in it, the builder stops with `carton_pk N not in --carton-table`. In that case,
   dump the table again from targetdb.
3. **The semaphore mapping**: v4 from sdss-semaphore 0.2.9, the newest on 2026-11-19. Check
   that it is still the newest. 22.9 M of the 63.1 M pairs had no bit at all (see
   `2026_11_19/carton_flags/CARTON_FLAGS.md` and `SEMAPHORE_ISSUE_DRAFT.md`).
4. **The product registry** `allTargets_pass1_dr21.toml` (the arMADGICS datamodel draft) hard-codes
   the 44 columns of 0.4.0 in `column_order` and in the per-column descriptions. For 0.4.1 it
   should take the column list from the file's `column_order` attribute, and fall back to the
   dataset `description` for new columns (its header already says it does this).
5. **Duplicate content**: `gaia_source_id` and `gaia_dr3_source` are both "Gaia DR3 source ID"
   in `/meta`, and both are copied.
6. **Decisions for the maintainer**: the 0.4.1 output file name; default or `--extended` for the 0.4.1
   targets build; whether to ship `meta_sdss5_target_flags`.

## Files

- `build_allTargets.jl`: the builder. `Project.toml` and `Manifest.toml` are copied from the 0.4.0
  folder (HDF5.jl 0.17.4, ArgParse, JSON, SHA).
- `tests/compare_to_release.jl`: regression against a released file. For each dataset it
  compares the datatype, shape, layout, values and attributes, and then lists the datasets and
  file attributes that differ.
- `tests/verify_passthrough.jl`: checks that each `/meta` column reached the output unchanged,
  and prints the solar-system row.
- `tools/inspect_meta.jl`: lists a `/meta` group with types, shapes and attributes.
- `tools/make_allexposure_like.jl`: test helper only. It builds an allExposure-like `sdss_id`
  table from an almanac file.
- `reference/targetdb_carton.csv`: targetdb.carton joined to mapper and version (941 cartons,
  54 kB), dumped once from targetdb on 2026-11-18; no dump script was kept. Committed as a fixture.
  To regenerate, dump `carton_pk, carton, program, mapper, plan, tag` from targetdb (query not
  recorded; reconstruct and diff against this file).
- `reference/sdss5_target_4_with_groups.csv`: the sdss-semaphore 0.2.9 carton-to-bit mapping,
  version 4, byte-identical to the file shipped upstream in sdss-semaphore 0.2.9. Committed as a
  fixture; regenerate by copying it from the sdss-semaphore release.
- `reference/allTargets_pass1_dr21.toml`: the 0.4.0 product registry entry (datamodel draft),
  source of the 0.4.0 `column_order` that `ORDER_040` in the builder reproduces.
- `test_regression_040/`, `test_extended/`, `test_misc/` (work folder only, not in git): test
  outputs and logs.
