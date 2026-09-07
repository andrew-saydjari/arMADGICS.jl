#!/bin/bash
# Shard entrypoint for the full telluric refit.
# Usage:  run_shard.sh [SHARD_ID]      (defaults to $SLURM_ARRAY_TASK_ID)
#
# Reads shards/shard_map.tsv for the shard's input list + output h5, runs
# fit_domeflats.py against the per-telescope support JSON
# ($PROJ/support_<tele>.json, produced by telluric_support.py) with the
# validated production settings, and maintains
# status/shard_<id>.STATUS (RUNNING / DONE / FAILED* states).  Restart-safe:
# the fitter skips rows already marked success in the shard output, so
# resubmitting a killed/failed task resumes where it left off.
set -u
# PROJ is the RUN directory: it holds the frozen .pkl inputs, the input lists,
# shards/, out/, status/ and the support JSONs.  CODE is where this script and
# fit_domeflats.py live (the repo).  They are usually different.
PROJ="${PROJ:-$(pwd)}"
CODE="${CODE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
UV="${UV:-$HOME/.local/bin/uv}"

SHARD_ID="${1:-${SLURM_ARRAY_TASK_ID:-}}"
if [ -z "$SHARD_ID" ]; then
    echo "ERROR: no shard id (arg 1 or SLURM_ARRAY_TASK_ID)" >&2
    exit 2
fi

cd "$PROJ" || exit 2   # frozen inputs (T_init.pkl etc.) are read from CWD
mkdir -p status out

row=$(awk -F'\t' -v s="$SHARD_ID" '$1 == s' shards/shard_map.tsv)
if [ -z "$row" ]; then
    echo "ERROR: shard $SHARD_ID not in shards/shard_map.tsv" >&2
    exit 2
fi
TELE=$(echo "$row" | cut -f2)
NFILES=$(echo "$row" | cut -f3)
LIST=$(echo "$row" | cut -f4)
OUT=$(echo "$row" | cut -f5)
STATUS="status/shard_${SHARD_ID}.STATUS"

write_status() {
    {
        echo "state=$1"
        echo "shard=$SHARD_ID"
        echo "tele=$TELE"
        echo "n_files=$NFILES"
        echo "host=$(hostname)"
        echo "slurm_job=${SLURM_JOB_ID:-none}"
        echo "updated=$(date -Is)"
        [ -n "${2:-}" ] && echo "detail=$2"
    } > "$STATUS".tmp && mv "$STATUS".tmp "$STATUS"
}

# Environment: jax[cuda12] wheels bundle their CUDA runtime -> no modules needed.
export JAX_PLATFORMS=""            # let jax pick GPU; preflight below verifies
if [ -n "${SLURM_JOB_ID:-}" ]; then
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
    RUNPFX=""
else
    # Local smoke on ccalin051: shared-node etiquette (<=4 cores, nice, capped GPU mem)
    export OMP_NUM_THREADS=4
    export XLA_PYTHON_CLIENT_MEM_FRACTION="${XLA_PYTHON_CLIENT_MEM_FRACTION:-0.4}"
    RUNPFX="nice -n 10 taskset -c 0-3"
fi

write_status RUNNING "preflight"
if ! $RUNPFX "$UV" run python -c '
import jax
devs = jax.devices()
print("jax", jax.__version__, "devices:", devs)
assert any(d.platform == "gpu" for d in devs), "no GPU visible to jax"
'; then
    write_status FAILED_NOGPU "jax could not see a GPU (driver/CUDA compat?)"
    exit 1
fi

write_status RUNNING "fitting"
SUPPORT="$PROJ/support_${TELE}.json"
if [ ! -f "$SUPPORT" ]; then
    write_status FAILED_NOSUPPORT "missing $SUPPORT (run telluric_support.py first)"
    exit 1
fi

$RUNPFX "$UV" run python "$CODE/fit_domeflats.py" \
    --support "$SUPPORT" \
    --input  "$LIST" \
    --output "$OUT" \
    --no-figures \
    --s2 --s2-iters 50 --lambda-t 1e4 --update-t --s-pixels 1000
rc=$?

# Completion check: every row of the shard output marked success
ndone=$($RUNPFX "$UV" run python -c "
import h5py, sys
try:
    with h5py.File('$OUT', 'r') as f:
        print(int(f['success'][:].sum()))
except Exception as e:
    print(-1)
")
if [ "$rc" -eq 0 ] && [ "$ndone" -eq "$NFILES" ]; then
    write_status DONE "n_success=$ndone/$NFILES"
    echo "shard $SHARD_ID DONE ($ndone/$NFILES)"
    exit 0
else
    write_status INCOMPLETE "rc=$rc n_success=$ndone/$NFILES (resubmit to resume)"
    echo "shard $SHARD_ID INCOMPLETE rc=$rc n_success=$ndone/$NFILES" >&2
    exit 1
fi
