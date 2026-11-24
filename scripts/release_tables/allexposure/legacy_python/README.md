# legacy_python: the 0.4.0 reference implementation of allExposure

These are the scripts that built the 0.4.0 `allExposure_pass1_dr21.h5`, kept so the bit-identity
regression of the Julia builder stays reproducible. **Use the Julia builder (`../bin/build_allexposure.jl`)
for new releases.**

| step | script | 0.4.0 source (work area) | what it does |
|---|---|---|---|
| 1 | `01_build_identity.jl` | `2026_09_13/validation_anomaly_hunt/s1_build_identity.jl` | full_list_info.h5 strings -> identity.h5 integers |
| 2 | `02_harvest_times.py` | `2026_09_13/allexposure/s4_harvest_times.py` | mid-exposure times from ar1Dunical files (process pool) |
| 3 | `03_build_allexposure.py` | `2026_09_13/allexposure/s6_build_allexposure.py` | almanac join, barycentric correction, velocities (44 datasets) |
| 4 | `04_add_bary_flag.py` | `2026_09_13/allexposure/s9_add_flags.py` | adds `bary_flag` in place |
| 5 | `05_add_galactic_qa.py` | `2026_09_13/allexposure/s10_add_galactic_qa.py` | adds `glon`, `glat`, `qa_bitmask` in place |

After step 5 the 0.4.0 chain added, by copy-and-swap, `v_rad_corr`, `v_rad_offset_applied`,
`v_rad_corr_flag` (`2026_09_16/fiber_offsets/s10_final_build.jl`) and `offtarget`
(`2026_09_16/offtarget_sidecar/s1,s2,s4`), and one root attribute, `fiber_offset_convention`, for which
no script was found (its text is in `../src/schema.jl`). Those steps are stages of the Julia builder;
the Julia copy-and-swap rewrote the root attributes as fixed-length strings, which is why the 0.4.0
file mixes fixed- and variable-length string attributes.

Changes from the dated-folder originals: hard-coded paths and two release-specific assert constants
became arguments; header comments added. No computed value changed. The regression
(`../../README.md`, "Regression") ran this chain on 2026-11-24 and compared it with the release.
