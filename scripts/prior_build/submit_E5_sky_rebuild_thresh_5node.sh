#!/bin/bash
#SBATCH --partition=cca
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=2
#SBATCH --constraint=genoa

#SBATCH --mem=0
#SBATCH --time=4:00:00
#SBATCH --job-name=ar_E5_sky_rebuild
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
# 5-NODE genoa variant of submit_E5_sky_rebuild_thresh.sh (AKS 2026-09-07), for when
# fewer than 7 genoa nodes are free:
#
#   cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E5_run && \
#     sbatchAKS /mnt/home/asaydjari/gitcode/worktrees/arM-E5b/scripts/prior_build/submit_E5_sky_rebuild_thresh_5node.sh "combined:telemaj_union"
#
# NOTE, so this file does not get used out of habit: it is a CONVENIENCE, not a
# necessity. --nodes can be overridden directly on the sbatchAKS command line, which
# passes "$@" through to sbatch, and a CLI --nodes beats the #SBATCH directive:
#
#   sbatchAKS .../submit_E5_sky_rebuild_thresh.sh --nodes=5 "combined:telemaj_union"
#
# That is exactly how AKS launched job 6995960. The only thing this file buys is that
# the node count is recorded in the script rather than in shell history. Use whichever
# is clearer at the time -- they produce identical jobs.
#
# Like the icelake variant this is a THIN WRAPPER over the shared body: header only, no
# copied logic. The SLURM_NTASKS bug fixed on 2026-09-07 came from copying a body
# between scripts and changing the header without re-checking the body; more copies
# would invite it back. The body derives everything from SLURM_NNODES,
# SLURM_NTASKS_PER_NODE and /proc/meminfo, so it needs no edit for a node-count change.
#
# SIZING (unchanged from the 7-node parent; genoa is still REQUIRED, never rome):
#   MEASURED peak RSS 19.0 GB/worker; 48 x 19 = 912 GB of genoa's 1538 GB (59%).
#   The same 48 workers on a 1 TB rome node would be 89-93% -- a live OOM risk.
#
# WALL TIME: 5 x 48 = 240 workers vs 600 fibers => 3 waves at 83% efficiency.
#   3 waves x ~38.7 min/fiber (MEASURED) = ~1.9-2.1 h, vs ~1.3-1.6 h for 7 nodes.
#   The extra ~0.5 h is worth it whenever 5 genoa nodes are free and 7 are not, since
#   the queue wait for the 6th and 7th node dominates that difference.
# ------------------------------------------------------------------------------
exec bash "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/submit_E5_sky_rebuild_thresh.sh" "$@"
