# W2 MPI tier — multi-node smoke test (PREPARED, NOT SUBMITTED)

`submit_mpi_smoke.sbatch` in this directory is ready for you to submit
manually; nothing here has been (or will be) submitted automatically.

## Why a Slurm job at all

Everything else about the MPI tier was validated locally on ccalin051
(single node, `mpiexec` with ≤4 ranks): parallel HDF5 loads, MPIO
create/write/close on ceph works, the locking matrix was measured, and a
one-fiber `workup_mpi.jl` run was compared bitwise against the serial tier.
The one thing local testing CANNOT prove is cross-node behavior: with ranks
on two different ceph clients, independent MPIO writes to one shared file
exercise ceph's client cache-coherence, cross-node HDF5 file locking, and
collective metadata operations over InfiniBand. That is exactly what this
job tests, at deliberately small scale (2 nodes × 4 ranks, ~7 GB per case).

## What the job does

1. `mpio_smoke.jl` cross-node, three locking settings:
   `HDF5_USE_FILE_LOCKING=FALSE` (expected PASS), `TRUE` (the known ceph
   hazard — may fail; its failure is a result, not a problem), and unset.
   Each case writes an (8700 × 100000) Float64 dataset rank-strided in
   100-row batches and verifies sampled rows bitwise after serial reopen.
2. A real `workup_mpi.jl` run of fiber 2 of the 2026_05_01 partial corpus
   (550 batches, ~50 GB read / ~50 GB written) into the job's scratch dir.
   Compare it against a serial-tier run of the same fiber with
   `compare_workup_outputs.jl` (command printed at the end of the job log)
   — the verdict must be IDENTICAL.

## How to submit

```bash
cd /mnt/home/asaydjari/gitcode/arMADGICS.jl     # or the worktree
mkdir -p slurm_logs
sbatchAKS workup/mpi_env/submit_mpi_smoke.sbatch
```

Notes:
* `REPO` at the top of the sbatch script points at `~/gitcode/arMADGICS.jl`;
  edit it if you submit from a worktree.
* The job writes only under
  `/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/w2_mpi_smoke/slurm_<jobid>/`
  and deletes its smoke files; the fiber-2 workup output (~50 GB) is left
  for the comparison and can be removed afterwards.
* `srun --mpi=pmix` is the expected launch for the module OpenMPI; if PMIx
  negotiation fails on the chosen nodes, `mpirun -np 8` from within the
  allocation is the fallback.
* Interpretation: if lockFALSE passes cross-node and the fiber-2 comparison
  is IDENTICAL, the MPI tier is cleared for production consideration; if
  only single-node passes, the per-rank-files + concat fallback (or the
  serial tier, which the W2 profiling already sized) is the path.
