#!/bin/bash
# ------------------------------------------------------------------------------
# Cross-node MPIO-on-ceph smoke test for the W2 MPI workup tier.
# PREPARED, NOT SUBMITTED — see SUBMIT_MPI_TEST.md. Submit manually with
# sbatchAKS (or sbatch) after reading that file.
#
# 2 nodes x 4 ranks: the thing local (ccalin051) testing cannot prove is that
# INDEPENDENT MPIO writes to one shared HDF5 file on ceph remain correct when
# the ranks are on DIFFERENT clients (page-cache coherence across ceph
# clients, file locking across nodes, collective metadata ops over IB).
# ------------------------------------------------------------------------------
#SBATCH --partition=cca
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --constraint="[genoa|icelake|rome]"
#SBATCH --time=0-01:00
#SBATCH --job-name=w2_mpio_smoke
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------

set -euo pipefail

REPO=${REPO:-/mnt/home/asaydjari/gitcode/arMADGICS.jl}   # env-overridable; worktree default until PR #26 merges
SCRATCH=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_02/w2_mpi_smoke/slurm_${SLURM_JOB_ID}
mkdir -p "$SCRATCH"

module load cephtweaks openmpi/5.0.6 hdf5/mpi-1.14.5
export CEPHTWEAKS_LAZYIO=1

run_case () {
    local label=$1; shift
    echo "=========== case: $label ==========="
    # ~7 GB dataset so writes clearly exceed any single client's cache
    srun --mpi=pmix julia --project=$REPO/workup/mpi_env \
        $REPO/workup/mpi_env/mpio_smoke.jl \
        --out $SCRATCH/smoke_${label}.h5 \
        --npix 8700 --rows-per-batch 100 --batches-per-rank 125
    rm -f $SCRATCH/smoke_${label}.h5
}

# locking matrix, cross-node
export HDF5_USE_FILE_LOCKING=FALSE; run_case lockFALSE
export HDF5_USE_FILE_LOCKING=TRUE;  run_case lockTRUE  || echo "lockTRUE FAILED (expected hazard on ceph)"
unset HDF5_USE_FILE_LOCKING;        run_case lockUNSET || echo "lockUNSET FAILED"

# small real cross-node workup: fiber 2 of the 2026_05_01 partial corpus
export HDF5_USE_FILE_LOCKING=FALSE
srun --mpi=pmix julia --project=$REPO/workup/mpi_env \
    $REPO/workup/workup_mpi.jl \
    --rawdir /mnt/ceph/users/sdssv/work/asaydjari/2026_05_01/outdir/arMADGICS/raw_57600_61160 \
    --outdir $SCRATCH/wu_fiber2_mpi \
    --fibers 2:2

echo "Now compare against a serial-tier reference of the same fiber, e.g.:"
echo "  julia --project=$REPO/workup workup/compare_workup_outputs.jl \\"
echo "      --a <serial fiber-2 outdir> --b $SCRATCH/wu_fiber2_mpi"
echo "DONE"
