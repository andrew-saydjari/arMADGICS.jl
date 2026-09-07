#!/bin/bash
#SBATCH --partition=cca
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2
#SBATCH --constraint=icelake

#SBATCH --mem=0
#SBATCH --time=4:00:00
#SBATCH --job-name=ar_E5_sky_rebuild
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
# ICELAKE variant of submit_E5_sky_rebuild_thresh.sh (AKS 2026-09-07), for when genoa
# is unavailable. Fire exactly like the genoa version:
#
#   cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E5_run && \
#     sbatchAKS /mnt/home/asaydjari/gitcode/worktrees/arM-E5b/scripts/prior_build/submit_E5_sky_rebuild_thresh_icelake.sh "combined:telemaj_union"
#
# This is a THIN WRAPPER: it changes ONLY the #SBATCH header and then execs the shared
# body. It deliberately does NOT copy the body. The SLURM_NTASKS bug fixed on
# 2026-09-07 was caused precisely by copying a body between scripts and changing the
# header without re-checking the body against it; a third copy would invite the same
# failure. Everything the body needs it derives from the Slurm environment
# (SLURM_NNODES, SLURM_NTASKS_PER_NODE, /proc/meminfo), so it is node-type agnostic.
#
# WHY --ntasks-per-node=32 AND WHY IT IS NOT OPTIONAL:
#   icelake is 64 cores / 1 TB. At the required --cpus-per-task=2 that is 32 tasks per
#   node, not 48. Do NOT simply swap the constraint on the genoa script: at 48 tasks x 2
#   CPUs = 96 CPUs the job could never be allocated on a 64-core node, and even if it
#   were, the body's memory guard would refuse to start (48 x 19 GB = 912 GB against a
#   memory-safe max of ~42 workers on 1 TB). That refusal is the guard working, not a
#   misconfiguration to be overridden.
#   MEASURED sizing basis (2026_09_06/gspice_resource_profile/GSPICE_RESOURCE_REPORT.md):
#   peak RSS 19.0 GB/worker, CPU-bound at 93.3% of 2 cores, iowait 0.00.
#   32 workers x 19 GB = 608 GB of icelake's ~1000 GB (61%) -- the same headroom fraction
#   as 48 x 19 = 912 GB of genoa's 1538 GB (59%). Memory-safe.
#
# WALL TIME (7 nodes x 32 = 224 workers vs 600 fibers => 3 waves, 89% efficiency):
#   MEASURED per-fiber rebuild cost on the profiled hardware: ~38.7 min.
#   3 waves x 38.7 min = ~1.9 h IF icelake matched that per-core throughput.
#   INFERRED: it does not. The profile is Zen-4-class; icelake (Xeon Ice Lake) is
#   typically 1.2-1.4x slower on this BLAS-heavy per-fiber work, so budget ~2.3-2.7 h.
#   That still fits the 4 h walltime, but the margin is real -- if it matters, raise the
#   node count rather than the walltime:
#     --nodes=10 -> 320 workers, 2 waves, ~1.3 h nominal / ~1.6-1.8 h inferred (94% eff)
#     --nodes=19 -> 608 workers, 1 wave,  ~0.65 h nominal / ~0.8-0.9 h inferred (99% eff)
#   --nodes can be overridden on the sbatchAKS command line; it beats this directive.
#
# DO NOT use a bracketed "[genoa|icelake]" constraint to hedge between the two.
#   sbatch(1) "Matching OR" is unambiguous: square brackets mean ALL allocated nodes
#   share ONE feature, so the allocation would be homogeneous (the mixed case is the
#   BARE "genoa|icelake", without brackets -- the opposite of the common assumption).
#   Homogeneity is not the problem here; --ntasks-per-node is. One value cannot be right
#   for both a 96-core and a 64-core node. In practice --ntasks-per-node=48 --cpus-per-task=2
#   demands 96 CPUs/node, which no icelake node can satisfy, so a bracketed constraint
#   would silently resolve to genoa-only and buy nothing. Pick ONE constraint explicitly.
# ------------------------------------------------------------------------------
exec bash "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/submit_E5_sky_rebuild_thresh.sh" "$@"
