#!/bin/bash
# ------------------------------------------------------------------------------
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --constraint="[genoa|icelake|rome]"
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=0-04:00
#SBATCH --job-name=release_tables
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
# Template, never submitted by the scripts. Set the variables listed in run_release_tables.sh
# (RELEASE, NAME, OUTDIR, WORKDIR, FULL_LIST, ALMANAC, APRED, WORKUP, ...) in the environment
# before `sbatch`. THREADS / WORKERS default to SLURM_CPUS_PER_TASK / SLURM_CPUS_ON_NODE.
# The worker processes of the allExposure reads are I/O-latency bound (a few % CPU each), so
# WORKERS may exceed the core count there.
# NOTE: NUMERICS=legacy reproduces 0.4.0 bits only on AVX-512 CPUs (numpy SVML); rome has no
# AVX-512, so use "[genoa|icelake]" for a bit-identity rebuild, or NUMERICS=native.
set -euo pipefail
export THREADS=${THREADS:-$SLURM_CPUS_PER_TASK}
: "${RELEASE_TABLES_DIR:?set to the scripts/release_tables checkout (sbatch copies this file to a spool dir)}"
srun --cpus-per-task="$SLURM_CPUS_PER_TASK" "$RELEASE_TABLES_DIR/run_release_tables.sh"
