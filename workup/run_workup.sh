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
# Rank sizing (user decision 2026-09-03: no blind static default):
#   WORKUP_RANKS=auto   (default) size from data + node resources, PER NODE
#                       (memory does not pool across nodes):
#         per_rank_est  = BASE + INFLIGHT × per-batch payload bytes
#         ranks/node    = min( floor(WORKUP_MEM_FRACTION × mem_per_node
#                                     / per_rank_est), cpus_per_node ), ≥1
#         total ranks   = ranks/node × nnodes, capped by #batches, ≥1
#   WORKUP_RANKS=<int>  explicit total (the auto arithmetic is still printed
#                       for the audit trail)
#   WORKUP_MEM_FRACTION (default 0.90) memory headroom factor
#   WORKUP_SIZING_DRYRUN=1  print the sizing block and exit 0 (audit mode)
#
# MPI launch logic:
#   * multi-node Slurm allocation → srun --mpi=pmix --ntasks=<total>
#     --ntasks-per-node=<ranks/node>  (enforces the sized distribution —
#     ranks must not pack onto one node)
#   * single node (ccalin051, or a 1-node allocation) → mpiexec -np <total>
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

# ---------------------------------------------------------------- sizing lib
# Runtime baseline per rank (Julia + MPI + parallel HDF5 + plan/identity
# share). Measured 2026-09-02/03 on ccalin051: workup_mpi.jl fiber-2 run,
# 4 ranks, /usr/bin/time max RSS 1.45 GB per rank process (rank 0, which
# also holds the 26.5M-row identity index, is the max); serial-tier readers
# ~1 GB, serial master 3.4 GB. 1500 MB is the conservative per-rank constant.
WORKUP_BASE_MB=${WORKUP_BASE_MB:-1500}
# In-flight batch payloads a rank holds at once: read buffer + the HDF5
# write path's copy + one batch of slack (serial tier additionally buffers
# up to 2×nworkers batches AT THE MASTER, which the base + headroom absorb).
WORKUP_INFLIGHT=${WORKUP_INFLIGHT:-3}

# min of a Slurm cpus-per-node spec like "72(x2),36" → 36. Heterogeneous
# allocations: we size for the SMALLEST node (memory heterogeneity is not
# visible in the environment, so homogeneous memory is assumed — see README).
_min_cpus_per_node() {
    tr ',' '\n' <<< "$1" | sed 's/(x[0-9]*)//' | sort -n | head -1
}

# sizing_report_and_choose <batch_bytes> <nbatch>
# Prints the auditable sizing block; sets AUTO_RANKS (total),
# AUTO_RANKS_PER_NODE, and SIZING_NNODES.
sizing_report_and_choose() {
    local batch_bytes=$1 nbatch=$2
    local frac=${WORKUP_MEM_FRACTION:-0.90}
    local batch_mb=$(( (batch_bytes + 1048575) / 1048576 ))
    local per_rank_mb=$(( WORKUP_BASE_MB + WORKUP_INFLIGHT * batch_mb ))

    local nnodes=1 cpus_pn mem_pn mem_src het_note=""
    if [ -n "${SLURM_JOB_ID:-}" ]; then
        nnodes=${SLURM_NNODES:-${SLURM_JOB_NUM_NODES:-1}}
        if [ -n "${SLURM_JOB_CPUS_PER_NODE:-}" ]; then
            cpus_pn=$(_min_cpus_per_node "$SLURM_JOB_CPUS_PER_NODE")
            [[ "$SLURM_JOB_CPUS_PER_NODE" == *,* ]] && \
                het_note=" (heterogeneous cpus '$SLURM_JOB_CPUS_PER_NODE' → using min)"
        else
            cpus_pn=${SLURM_CPUS_ON_NODE:-$(nproc)}
        fi
        if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
            mem_pn=$SLURM_MEM_PER_NODE; mem_src="SLURM_MEM_PER_NODE"
        elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
            mem_pn=$(( SLURM_MEM_PER_CPU * cpus_pn )); mem_src="SLURM_MEM_PER_CPU × cpus/node"
        else
            mem_pn=$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)
            mem_src="/proc/meminfo MemAvailable (no Slurm mem env)"
        fi
    else
        cpus_pn=${WORKUP_TEST_NPROC:-$(nproc)}
        mem_pn=${WORKUP_TEST_MEMAVAIL_MB:-$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)}
        mem_src="/proc/meminfo MemAvailable"
    fi

    local mem_cap_pn rpn total_raw total
    mem_cap_pn=$(awk -v a="$mem_pn" -v f="$frac" -v p="$per_rank_mb" 'BEGIN{print int(a*f/p)}')
    rpn=$mem_cap_pn
    [ "$cpus_pn" -lt "$rpn" ] && rpn=$cpus_pn
    [ "$rpn" -lt 1 ] && rpn=1
    total_raw=$(( rpn * nnodes ))
    total=$total_raw
    [ "$nbatch" -lt "$total" ] && total=$nbatch
    [ "$total" -lt 1 ] && total=1
    # if the work-unit cap bit, shrink ranks/node to match (ceil division)
    local rpn_launch=$(( (total + nnodes - 1) / nnodes ))

    echo "[run_workup:sizing] per-rank estimate: base ${WORKUP_BASE_MB} MB + ${WORKUP_INFLIGHT} in-flight × ${batch_mb} MB/batch = ${per_rank_mb} MB"
    echo "[run_workup:sizing] per-node: mem ${mem_pn} MB (${mem_src}), cpus ${cpus_pn}${het_note}"
    echo "[run_workup:sizing]   ranks/node = min( floor(${frac} × ${mem_pn} / ${per_rank_mb}) = ${mem_cap_pn}, cpus ${cpus_pn} ) = ${rpn}"
    echo "[run_workup:sizing] nodes: ${nnodes} → total = ${rpn} × ${nnodes} = ${total_raw}; work units: ${nbatch} batches → total = ${total} (${rpn_launch}/node)"

    AUTO_RANKS=$total
    AUTO_RANKS_PER_NODE=$rpn_launch
    SIZING_NNODES=$nnodes
}

# When sourced as a library (unit tests), stop here.
if [ "${WORKUP_SIZING_LIB:-0}" = "1" ]; then return 0 2>/dev/null || exit 0; fi

# ----------------------------------------------------------------------- main
usage() {
    echo "usage: run_workup.sh <rawdir> <redux> <outdir> [extra passthrough args]" >&2
    echo "  env: WORKUP_TIER=mpi|serial (default mpi), WORKUP_RANKS=auto|N (default auto)," >&2
    echo "       WORKUP_MEM_FRACTION=F (default 0.90), WORKUP_SIZING_DRYRUN=1 (audit only)" >&2
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
RANKS_SETTING=${WORKUP_RANKS:-auto}

[ -d "$RAWDIR" ] || { echo "run_workup.sh: rawdir not found: $RAWDIR" >&2; exit 2; }
export WORKUP_REDUX=$REDUX

echo "[run_workup] tier=$TIER rawdir=$RAWDIR outdir=$OUTDIR redux=$REDUX extra: $*"

# ---- rank sizing (always printed, even with explicit WORKUP_RANKS) ----------
F1=1; F2=600
_prev=""
for a in "$@"; do
    if [ "$_prev" = "--fibers" ] && [[ "$a" =~ ^([0-9]+):([0-9]+)$ ]]; then
        F1=${BASH_REMATCH[1]}; F2=${BASH_REMATCH[2]}
    fi
    _prev=$a
done

AUTO_RANKS=""; AUTO_RANKS_PER_NODE=1; SIZING_NNODES=1
if probe_out=$(julia --project="$WORKUP_DIR" "$WORKUP_DIR/size_probe.jl" "$RAWDIR" "$F1" "$F2" 2>&1); then
    BATCH_BYTES=$(sed -n 's/^BATCH_BYTES=//p' <<< "$probe_out")
    NBATCH=$(sed -n 's/^NBATCH=//p' <<< "$probe_out")
    SAMPLE=$(sed -n 's/^SAMPLE=//p' <<< "$probe_out")
    echo "[run_workup:sizing] probe: fibers ${F1}:${F2} → ${NBATCH} batches; sample ${SAMPLE} payload ${BATCH_BYTES} bytes"
    sizing_report_and_choose "$BATCH_BYTES" "$NBATCH"
else
    echo "[run_workup:sizing] WARNING: size probe failed:" >&2
    sed 's/^/[run_workup:sizing]   /' <<< "$probe_out" >&2
fi

if [ "$RANKS_SETTING" = "auto" ]; then
    if [ -z "$AUTO_RANKS" ]; then
        RANKS=8; RANKS_PER_NODE=8
        echo "[run_workup:sizing] auto requested but probe failed → falling back to ${RANKS} ranks" >&2
    else
        RANKS=$AUTO_RANKS; RANKS_PER_NODE=$AUTO_RANKS_PER_NODE
        echo "[run_workup:sizing] chosen: ${RANKS} ranks (${RANKS_PER_NODE}/node) [auto]"
    fi
else
    RANKS=$RANKS_SETTING
    RANKS_PER_NODE=$(( (RANKS + SIZING_NNODES - 1) / SIZING_NNODES ))
    if [ -n "$AUTO_RANKS" ]; then
        echo "[run_workup:sizing] chosen: ${RANKS} ranks (${RANKS_PER_NODE}/node) [explicit WORKUP_RANKS; auto would have chosen ${AUTO_RANKS} (${AUTO_RANKS_PER_NODE}/node)]"
    else
        echo "[run_workup:sizing] chosen: ${RANKS} ranks (${RANKS_PER_NODE}/node) [explicit WORKUP_RANKS; auto unavailable (probe failed)]"
    fi
fi

if [ "${WORKUP_SIZING_DRYRUN:-0}" = "1" ]; then
    echo "[run_workup] WORKUP_SIZING_DRYRUN=1 → exiting after sizing (nothing launched)"
    exit 0
fi

# ---- serial tier ------------------------------------------------------------
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
    echo "[run_workup] multi-node Slurm allocation (${SLURM_JOB_NUM_NODES} nodes) → srun --mpi=pmix --ntasks=$RANKS --ntasks-per-node=$RANKS_PER_NODE"
    exec srun --mpi=pmix --ntasks="$RANKS" --ntasks-per-node="$RANKS_PER_NODE" \
        julia --project="$MPI_ENV" "${WRITER[@]}"
else
    echo "[run_workup] single node → mpiexec -np $RANKS"
    exec mpiexec -np "$RANKS" julia --project="$MPI_ENV" "${WRITER[@]}"
fi
