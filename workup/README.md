# arMADGICS workup

Streams the per-fiber batch files (`NNN/arMADGICS_fiber_NNN_batch_SSSSSSS.h5`)
into one preallocated `arMADGICS_out_<key>.h5` per key, with row placement
derived from IDENTITY (`full_list_info.h5` + the batch filename), never from
file-list position — the W1 row contract (`RowContract.jl`).

**MPI is the production default tier** (user decision 2026-09-03, after the
cross-node Slurm smoke passed all three file-locking cases and the fiber-2
MPI-vs-serial comparison was verdicted IDENTICAL). The serial tier remains
(a) the `WORKUP_TIER=serial` fallback when the module MPI stack is
unavailable and (b) the reference implementation for W4-style regressions.

## Entrypoint contract — `run_workup.sh` (call this, nothing else)

```
workup/run_workup.sh <rawdir> <redux> <outdir> [extra passthrough args]
```

| positional | meaning | mapped to |
|---|---|---|
| `<rawdir>` | raw batch dir: contains `NNN/` fiber subdirs, `batch_info.txt`, `full_list_info.h5` | writer `--rawdir` |
| `<redux>` | reduxBase (dir containing `apred/<mjd>/ar1Dunical_*.h5`) | NOT consumed by the writers; logged and exported as `WORKUP_REDUX` so the caller runs the W3 validator (`validate_workup.jl --rawdir <rawdir> --redux <redux>`) against the same inputs |
| `<outdir>` | output dir (created if absent) for the `arMADGICS_out_<key>.h5` files | writer `--outdir` |
| `[extra]` | appended verbatim to the writer CLI | see below |

Passthrough flags valid for BOTH tiers: `--fibers F1:F2` (contiguous fiber
window), `--allow-missing`, `--batchsize N` (default 100),
`--progress-every N`. Serial tier only: `--resume`, `--ckpt-every N`.

Environment knobs:

| var | default | effect |
|---|---|---|
| `WORKUP_TIER` | `mpi` | `serial` → run `workup_serial.jl` instead (no modules needed) |
| `WORKUP_RANKS` | `auto` | `auto` = size from data + node resources (below); an explicit integer is honored as the TOTAL rank count (the auto arithmetic is still printed for audit) |
| `WORKUP_MEM_FRACTION` | `0.90` | memory headroom factor in the auto sizing |
| `WORKUP_SIZING_DRYRUN` | unset | `1` → print the sizing block and exit 0 (audit mode, launches nothing) |
| `HDF5_USE_FILE_LOCKING` | `FALSE` | forwarded to HDF5 (MPI tier) |
| `CEPHTWEAKS_LAZYIO` | `1` | ceph lazy I/O (MPI tier) |

### Rank auto-sizing (`WORKUP_RANKS=auto`, the default)

Memory does not pool across nodes, so sizing is per node first:

```
per_rank_est   = BASE (1500 MB) + INFLIGHT (3) × per_batch_payload_MB
ranks_per_node = min( floor(WORKUP_MEM_FRACTION × mem_per_node / per_rank_est),
                      cpus_per_node ),  ≥ 1
total_ranks    = ranks_per_node × nnodes, capped by the number of batch files, ≥ 1
```

- `per_batch_payload` comes from a startup size probe (`size_probe.jl`):
  the same key/shape/dtype discovery RowContract performs, on a sample
  batch, summed over all keys — fatter future keys automatically shrink
  the rank count. The probe also counts the batch files in the requested
  `--fibers` window (the work-unit cap).
- `BASE` = 1500 MB/rank is the measured runtime baseline (Julia + MPI +
  parallel HDF5; max RSS 1.45 GB/rank in the 2026-09-02 fiber-2 MPI run,
  rank 0 carrying the 26.5M-row identity index); `INFLIGHT` = 3 batch
  payloads per rank (read buffer + write-path copy + slack). Both are
  env-overridable (`WORKUP_BASE_MB`, `WORKUP_INFLIGHT`) but should not
  normally be touched.
- Resources: inside Slurm, `SLURM_MEM_PER_NODE` (or
  `SLURM_MEM_PER_CPU × cpus/node`), `SLURM_JOB_CPUS_PER_NODE` (minimum
  entry on heterogeneous allocations — memory heterogeneity is not visible
  in the environment, so homogeneous memory per node is assumed), and
  `SLURM_NNODES`. Outside Slurm: `/proc/meminfo` MemAvailable + `nproc`
  (nnodes = 1).
- The full arithmetic (per-rank estimate breakdown, per-node resources,
  the min(), the ×nnodes multiplication, and the work-unit cap) is printed
  as `[run_workup:sizing]` lines on EVERY run — also when `WORKUP_RANKS`
  is explicit, so the log always shows what auto would have chosen.
- Note: current workloads are I/O-bound and per-batch payloads are ~91 MB,
  so auto typically lands **cpu-capped** (e.g. 32 on ccalin051, where the
  memory cap would allow ~237). The RAM cap exists to protect future
  fatter keys / leaner nodes; the mocked-env unit test
  (`test/test_run_workup_sizing.sh`) covers both regimes plus multi-node.

Launch logic (encapsulated — callers do NOT load modules or pick launchers):

- `WORKUP_TIER=serial` → `julia --project=workup workup_serial.jl ...`
  (the sized total becomes `--nworkers`)
- MPI tier: loads `cephtweaks openmpi/5.0.6 hdf5/mpi-1.14.5`, exports the
  ceph env, then
  - multi-node Slurm allocation (`SLURM_JOB_NUM_NODES > 1`) →
    `srun --mpi=pmix --ntasks=<total> --ntasks-per-node=<ranks/node>
    julia --project=workup/mpi_env workup_mpi.jl ...` — the explicit
    `--ntasks-per-node` enforces the sized distribution so ranks land as
    computed instead of packing onto one node
  - otherwise (ccalin051 or a 1-node allocation) →
    `mpiexec -np <total> julia --project=workup/mpi_env workup_mpi.jl ...`
- If the modules fail to load (wrong cluster, no Lmod), it exits 3 with a
  message naming the `WORKUP_TIER=serial` fallback. It never silently
  degrades tiers on its own.

Exit status: `0` only on a fully written workup. Non-zero on ANY failure,
including: missing batches without `--allow-missing` (row-shift protection),
unexpected extra batch files, and any batch failing the per-batch integrity
check (truncated/corrupt file → the run ABORTS; nothing further is written).
Callers must treat non-zero as "no usable workup".

Outputs in `<outdir>`: one `arMADGICS_out_<key>.h5` per key discovered from
the batches (36 for current production batches). Each file: dataset `<key>`
of shape `(leading..., nsamp)`, a `hdr` dataset carrying the batch
producer's git-provenance attrs, file-level `workup_*` attrs (rawdir, fiber
window, row offset, nsamp, timestamp, writer), and — only with
`--allow-missing` — a `missing_row_mask` dataset (UInt8, 1 = row belongs to
a missing batch and holds sentinel values: NaN for floats, typemin for
ints). Raw batch files are NEVER deleted (W5: deletion is a separate,
validation-gated step).

Example (the dailies pattern):

```bash
workup/run_workup.sh \
    /path/to/outdir/arMADGICS/raw \
    /path/to/outdir \
    /path/to/outdir/arMADGICS/wu_th \
    --allow-missing
```

## Pieces

| file | role |
|---|---|
| `run_workup.sh` | THE entrypoint (above) |
| `size_probe.jl` | data-side inputs for the rank auto-sizing (per-batch payload bytes + batch count) |
| `RowContract.jl` | W1 identity/row-range/expected-set/integrity logic (shared by everything) |
| `WorkupSerial.jl` | planning + output preallocation + streaming engine |
| `workup_serial.jl` | serial-tier CLI: Distributed readers, single writer |
| `workup_mpi.jl` | MPI-tier CLI: rank-strided readers, shared-file MPIO writes |
| `validate_workup.jl` | W3 row-matching validator (run after every workup) |
| `compare_workup_outputs.jl` | streaming bitwise comparison of two workup output dirs |
| `mpi_env/` | pinned module-stack Julia env for the MPI tier (see its README) |
| `regression/` | W4 reference: the old in-memory workup, minimally adapted |
| `test/` | unit tests (`julia --project=. test/runtests.jl` from `workup/`) |
