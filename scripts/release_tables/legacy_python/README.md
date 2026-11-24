# legacy_python: FROZEN Python helpers (regression reference only)

These are the Python versions of the small release-table helpers, as they were when 0.4.0 was rebuilt
on 2026-11-24. **They are not used by `run_release_tables.sh`.** Each was ported to Julia, and the Julia
port reproduces it exactly on the 0.4.0 inputs (2026-11-24):

| Python (frozen) | Julia (used) | reproduction |
|---|---|---|
| `04_attrs_allvisitrv.py` | `../allvisitrv/04_attrs_allvisitrv.jl` | all values, attributes and HDF5 string types equal (`compare_h5`) |
| `03_stage_allvisitrv.py` | `../offtarget/03_stage_allvisitrv.jl` | equal to the Python output and to the released allVisitRV (`compare_h5`); file bytes differ only in HDF5 object-header timestamps/layout |
| `make_restrict_ids.py` | `../metadata/make_restrict_ids.jl` | `.npy` byte-identical |
| `recommended_cuts.py` | `../cuts/recommended_cuts.jl` | same counts (25,516,821 / 25,207,025), same output and `--expect` exit codes |
| `compare_h5.py` | `../validation/compare_h5.jl` | same verdicts, mismatch counts, max ulp, stdout and JSON on the test pairs |

Do not edit them. The 0.4.0 allExposure reference chain is in `../allexposure/legacy_python/`.
