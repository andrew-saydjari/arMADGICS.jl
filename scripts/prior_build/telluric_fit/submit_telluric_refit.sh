#!/bin/bash
# ------------------------------------------------------------------------------
# FULL telluric refit — sbatch, one GPU per shard on one 8-GPU node.
#
# PREREQUISITE: support_apo.json and support_lco.json must exist in $PROJ
#   (uv run python telluric_support.py --list list_<t>_full.txt \
#        --telescope <t> --out support_<t>.json --npz live_<t>.npz)
#
# Submit (from the RUN directory):   PROJ=$PWD sbatchAKS <this script>
# Resubmit only unfinished shards, e.g.:  sbatchAKS --array=2,5 submit_telluric_refit.sh
#   (restart-safe: each shard resumes from its own out/shard_<k>_<tele>.h5,
#    skipping rows already marked success)
#
# Partition: gpuxl Hopper pool (H100-94GB 4x/node, H200-141GB 8x/node).
#   - The fit runs in float64 (jax_enable_x64): Hopper FP64 >> RTX-class FP64,
#     so gpuxl is both the most available AND the fastest option for this job.
#   - NO --reservation: the rocky9 reservation EXPIRED 2026-08-31. The active
#     rocky8 reservation locks the A100 nodes (workergpu037-054) instead.
#   - gpuxl QoS per-user caps: 64 GPUs / 24 jobs / 3-day wall — 8x 1-GPU tasks
#     at 20 h fits comfortably.
#
# Sizing: largest shard = 2,359 files x 21.8 s/file (validated A6000 rate)
#   ~ 14.3 h; Hopper should be faster. -t 20:00:00 gives margin either way.
# ------------------------------------------------------------------------------
#SBATCH --job-name=telluric_refit
#SBATCH --partition=gpuxl
#SBATCH --constraint=h200
#SBATCH --nodes=1
#SBATCH --gpus=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=800G
#SBATCH --time=20:00:00
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
# NOTE (2026-09-03): the gpuxl QOS enforces MinTRES gres/gpu=4, so the original
# 8x 1-GPU job array was rejected (QOSMinGRES). Repacked as ONE 8-GPU h200-node
# job running all 8 shards concurrently, one GPU each. Rerun a subset with e.g.:
#   SHARDS="2 5" sbatchAKS submit_telluric_refit.sh
# ------------------------------------------------------------------------------

# PROJ = run directory (frozen inputs, lists, support JSONs, out/, status/).
# CODE = directory holding run_shard.sh / fit_domeflats.py (this repo).
PROJ="${PROJ:?set PROJ to the run directory}"
CODE="${CODE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export PROJ CODE
SHARDS="${SHARDS:-0 1 2 3 4 5 6 7}"
gpu=0
pids=()
for k in $SHARDS; do
    CUDA_VISIBLE_DEVICES=$gpu bash "$CODE/run_shard.sh" "$k" > "$PROJ/slurm_logs/shard_${k}.out" 2>&1 &
    pids+=($!)
    gpu=$((gpu+1))
done
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done
echo "all shards finished, overall rc=$rc"
grep -H 'state=' "$PROJ"/status/shard_*.STATUS
exit $rc
