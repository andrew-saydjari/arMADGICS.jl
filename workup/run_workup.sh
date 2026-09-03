#!/bin/bash
# =============================================================================
# run_workup.sh — THE single entrypoint for the arMADGICS workup.
#
#   run_workup.sh <rawdir> <redux> <outdir> [extra passthrough args]
#
#     <rawdir>  arMADGICS raw batch dir (NNN/ fiber dirs + batch_info.txt +
#               full_list_info.h5)                       → writer --rawdir
#     <redux>   reduxBase (the outdir containing apred/<mjd>/ar1Dunical_*.h5).
#               NOT consumed by the writers; recorded in the log and exported
#               as WORKUP_REDUX so the caller / a following W3 validation step
#               (validate_workup.jl --redux) uses the same value.
#     <outdir>  output dir for arMADGICS_out_<key>.h5   → writer --outdir
#     [extra]   appended verbatim to the writer CLI. Common: --fibers F1:F2,
#               --allow-missing, --batchsize N, --progress-every N.
#               Serial-tier only: --resume, --ckpt-every N, --nworkers N.
#
# Tier selection (user decision 2026-09-03: MPI is the production default):
#   WORKUP_TIER=mpi     (default) workup_mpi.jl under the module MPI stack
#   WORKUP_TIER=serial  workup_serial.jl (no modules needed) — fallback and
#                       the W4 reference implementation
#
# MPI launch logic:
#   * inside a multi-node Slurm allocation → srun --mpi=pmix
#   * single node (ccalin051, or a 1-node allocation) → mpiexec -np N
#   * N = WORKUP_RANKS (default 8)
#
# Environment ceremony (encapsulated here; see workup/mpi_env/README.md):
#   module load cephtweaks openmpi/5.0.6 hdf5/mpi-1.14.5
#   export CEPHTWEAKS_LAZYIO=1 HDF5_USE_FILE_LOCKING=FALSE
# If the modules cannot be loaded, this script REFUSES to run the MPI tier
# and tells you to rerun with WORKUP_TIER=serial.
#
# Exit status: non-zero on ANY failure, including a batch integrity abort
# (truncated/corrupt batch) or a missing-batch reconciliation error.
# Raw batch files are never deleted (W5).
# =============================================================================
set -euo pipefail

usage() {
    echo "usage: run_workup.sh <rawdir> <redux> <outdir> [extra passthrough args]" >&2
    echo "  env: WORKUP_TIER=mpi|serial (default mpi), WORKUP_RANKS=N (default 8)" >&2
    exit 2
}

[ $# -ge 3 ] || usage
RAWDIR=$1
REDUX=$2
OUTDIR=$3
shift 3

WORKUP_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MPI_ENV=$WORKUP_DIR/mpi_env
TIER=${WORKUP_TIER:-mpi}
RANKS=${WORKUP_RANKS:-8}

[ -d "$RAWDIR" ] || { echo "run_workup.sh: rawdir not found: $RAWDIR" >&2; exit 2; }
export WORKUP_REDUX=$REDUX

echo "[run_workup] tier=$TIER rawdir=$RAWDIR outdir=$OUTDIR redux=$REDUX extra: $*"

if [ "$TIER" = "serial" ]; then
    exec julia --project="$WORKUP_DIR" "$WORKUP_DIR/workup_serial.jl" \
        --rawdir "$RAWDIR" --outdir "$OUTDIR" --nworkers "$RANKS" "$@"
elif [ "$TIER" != "mpi" ]; then
    echo "run_workup.sh: unknown WORKUP_TIER='$TIER' (expected 'mpi' or 'serial')" >&2
    exit 2
fi

# ---- MPI tier (default) -----------------------------------------------------
# make `module` available in non-interactive shells
if ! type module &>/dev/null; then
    for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh; do
        # shellcheck disable=SC1090
        [ -r "$f" ] && source "$f" && break
    done
fi
if ! type module &>/dev/null; then
    echo "run_workup.sh: FATAL: no 'module' command available — cannot set up the MPI stack." >&2
    echo "  Fallback: rerun with WORKUP_TIER=serial (no modules required)." >&2
    exit 3
fi

if ! module load cephtweaks openmpi/5.0.6 hdf5/mpi-1.14.5 2>/tmp/run_workup_module_err.$$; then
    echo "run_workup.sh: FATAL: 'module load cephtweaks openmpi/5.0.6 hdf5/mpi-1.14.5' failed:" >&2
    cat /tmp/run_workup_module_err.$$ >&2
    rm -f /tmp/run_workup_module_err.$$
    echo "  Fallback: rerun with WORKUP_TIER=serial (no modules required)." >&2
    exit 3
fi
rm -f /tmp/run_workup_module_err.$$

if ! command -v mpiexec &>/dev/null; then
    echo "run_workup.sh: FATAL: mpiexec not on PATH after module load — MPI stack broken." >&2
    echo "  Fallback: rerun with WORKUP_TIER=serial." >&2
    exit 3
fi

export CEPHTWEAKS_LAZYIO=${CEPHTWEAKS_LAZYIO:-1}
export HDF5_USE_FILE_LOCKING=${HDF5_USE_FILE_LOCKING:-FALSE}

WRITER=("$WORKUP_DIR/workup_mpi.jl" --rawdir "$RAWDIR" --outdir "$OUTDIR" "$@")

if [ -n "${SLURM_JOB_ID:-}" ] && [ "${SLURM_JOB_NUM_NODES:-1}" -gt 1 ]; then
    echo "[run_workup] multi-node Slurm allocation (${SLURM_JOB_NUM_NODES} nodes) → srun --mpi=pmix"
    exec srun --mpi=pmix julia --project="$MPI_ENV" "${WRITER[@]}"
else
    echo "[run_workup] single node → mpiexec -np $RANKS"
    exec mpiexec -np "$RANKS" julia --project="$MPI_ENV" "${WRITER[@]}"
fi
