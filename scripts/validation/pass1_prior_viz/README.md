# pass-1 prior gallery

`plot_pass1_priors.jl` renders, for N random fibers at each telescope, the samples and the
top-3 eigenvectors of **every** prior family the pass-1 run consumes (starCont pass-1c,
starLines E7 per-fiber LSF, skyCont E5, skyLine faint plain + faint GSPICE), plus
cross-fiber summaries and a fleet-wide (all 600 fibers) scan of leading-mode localization.

Run (any env with CairoMakie, ColorSchemes, HDF5):

```
julia --project=<env> scripts/validation/pass1_prior_viz/plot_pass1_priors.jl
```

Env overrides: `PRIORVIZ_OUT`, `PRIORVIZ_SEED` (default 20260907), `PRIORVIZ_NFIB`
(default 5 per telescope), `PRIORVIZ_FIBERS` (explicit comma list), `PRIORVIZ_FLEET=0`
to skip the 600-fiber scan. Every run also writes `prior_viz_qa.txt` with the structural
invariants, the eigenvalue accounting, and any anomalies flagged.

2026-09-07 output (seed 20260907; fibers APO 59/95/124/153/196, LCO 323/330/338/476/585)
is written to `/mnt/ceph/users/asaydjari/working/2026_09_07/plots/pass1_priors/`
(internal path, deliberately not linked publicly).

Two findings from that run are written up in that output directory: the deployed faint
sky-line prior spends u2/u3 on chip-edge and bright-line-wing pixels on 513/600 fibers,
and `src/priors.jl` still defaults `ARM_SKY_PRIOR_DIR` to the pre-rebuild `built/`.
