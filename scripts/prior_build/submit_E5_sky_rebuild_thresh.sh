#!/bin/bash
#SBATCH --partition=cca
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=2
#SBATCH --constraint=genoa

#SBATCH --mem=0
#SBATCH --time=4:00:00
#SBATCH --job-name=ar_E5_sky_rebuild
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
# E5 finding #35 REBUILD launcher: rebuild the skyLines priors under a chosen
# bright/faint threshold policy. Ready to fire for ANY of the four options in
# prior_outputs/sky_pass1/THRESHOLD_FINDING.md — no edits needed, pass the policy
# as the first argument (or set E5_THRESH_POLICY before sbatch; the value is
# re-read here at run time, so it is safe against sbatch's env snapshot):
#
#   sbatchAKS submit_E5_sky_rebuild_thresh.sh "<policy>"   # one 7-node job
#
#   sbatchAKS submit_E5_sky_rebuild_thresh.sh "linedetect"      # RECOMMENDED (calibrated)
#   sbatchAKS submit_E5_sky_rebuild_thresh.sh "abs:35,8"        # (A) match DR17 ~8.3% bright
#   sbatchAKS submit_E5_sky_rebuild_thresh.sh "abs:150,40"      # (B) physical bright lines, ~5%
#   sbatchAKS submit_E5_sky_rebuild_thresh.sh "quantile:0.083"  # (C) unit-free per-fiber quantile
#   sbatchAKS submit_E5_sky_rebuild_thresh.sh "off"             # (D) declared no-split
#
# Option (D) needs NO rebuild if the current job finishes: "off" is numerically
# identical to what the inherited constants already produce (they flag 0 pixels).
# It is provided so the intent can be recorded in the products rather than implied.
#
# Safety properties:
#  - writes to a POLICY-TAGGED dir (built_abs_apo35p0_lco8p0, built_q0p083,
#    built_nosplit, ...), never overwriting the current built/ products;
#  - per-fiber skip-if-built + atomic .claims/ inside that dir, so it resumes and
#    can run concurrently with anything else;
#  - reuses the existing threshold-invariant skyCont priors by symlink, so only the
#    skyLines half is recomputed (verified by e5_assert_skycont_invariant.jl, which
#    this script runs first and refuses to proceed without);
#  - bright-fraction guard set to ERROR here: a policy that flags outside 1-20% of
#    pixels stops the build instead of silently producing another no-op.
#
# DISTRIBUTION: SlurmManager, one multi-node job (AKS 2026-09-06).
# e5_sky_run.jl calls addprocs(SlurmManager()) when SLURM_JOB_ID is set, following the
# apMADGICS original (src/prior_build/build_skyLines.jl:9-10) and arM's port of it
# (scripts/prior_build/build_skyCont.jl:26-31). SlurmManager spans the whole allocation,
# so the 7 nodes requested here are all used. The header MUST therefore set
# --ntasks-per-node (one Slurm task per Julia worker) and --cpus-per-task=2 (the BLAS
# threads each worker needs). The script asserts, loudly, that workers land on all 7
# nodes: addprocs(::Int) would silently confine the job to one node, which is exactly the
# bug that went undetected in the earlier E5 build.
#
# SIZING RATIONALE (measured, 2026-09-06 GSPICE resource profile at
# 2026_09_06/gspice_resource_profile/GSPICE_RESOURCE_REPORT.md -- keep this with the script):
#   * CPU-bound: 93.3% of 2 cores/worker, 99.1% user time, iowait 0.00. No shared working
#     set between fibers, so nodes scale cleanly and I/O is a non-factor.
#   * Peak RSS 19.0 GB/worker (VmHWM 18.51-19.37, only 4% spread) = 9.5 GB/core.
#   * genoa is REQUIRED, not preferred: 16 GB/core vs rome's 8 GB/core against a 9.5
#     GB/core need. 48 workers x 19 GB = 912 GB of genoa's 1538 (59%); the same 48 on a
#     1024 GB rome node would be 89-93% -- a live OOM risk. NEVER rome.
#     Safe fallback if genoa queues badly: --constraint="[genoa|icelake]" with
#     --ntasks-per-node=32 (icelake has 64 cores / 1024 GB), ~1.9 h.
#   * 48 workers/node x 2 cores = 96 cores = a full genoa node.
#   * 7 nodes x 48 = 336 workers vs 600 fibers => 2 waves at ~38.7 min/fiber (rebuild cost;
#     skyCont is symlinked and genuinely skipped) => ~1.1-1.6 h at ~89% efficiency.
#     An 8th node buys ZERO wall-clock (still 2 waves) while idling more slots, so 7 is
#     the answer under AKS's effectiveness condition. AKS: "Let's do 7 for that in future."
#
# The atomic mkdir .claims/ machinery is KEPT even though SlurmManager now handles
# distribution: it costs nothing and still gives restart/resume safety, skip-if-built, and
# protection against an accidentally overlapping run.
# ------------------------------------------------------------------------------
SLURM_NTASKS=$(($SLURM_CPUS_ON_NODE * $SLURM_NNODES))
export SLURM_NTASKS
env | grep SLURM | while read -r line; do
    echo "$line"
done

set -e
set -o pipefail

hostname
echo ${SLURM_JOB_NODELIST:-no_nodelist}
echo "running from $(pwd)"

if [ -n "${SLURM_JOB_ID:-}" ] ; then
    script_path=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}' | awk '{print $1}')
else
    script_path=$(realpath "$0")
fi
base_dir="$(dirname "$(dirname "$(dirname "$script_path")")")"
echo "base_dir: $base_dir"

julia_version="1.11.0"
juliaup add $julia_version || true

print_elapsed_time() {
    current_seconds=$SECONDS
    if [ -n "${LAST_TIME:-}" ]; then
        diff_seconds=$((current_seconds - LAST_TIME))
        diff_time=$(printf '%dd %dh:%dm:%ds\n' $(($diff_seconds/86400)) $(($diff_seconds%86400/3600)) $(($diff_seconds%3600/60)) $(($diff_seconds%60)))
        echo
        echo "Time since last step: $diff_time"
    fi
    formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($current_seconds/86400)) $(($current_seconds%86400/3600)) $(($current_seconds%3600/60)) $(($current_seconds%60)))
    echo "Elapsed time: $formatted_time"
    echo
    echo "--------- $1 ---------"
    echo
    LAST_TIME=$current_seconds
}

# --- configuration (kept IN the script; sbatch reads this file at run time) ---
export E5_OUT=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1
export E5_REDUX=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/testbed_run/
export E5_ALMANAC=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/testbed_prep/almanac_testbed_dr21.h5
export E5_STARCONT=$E5_OUT/starcont_fresh/APOGEE_starcont_svd_60_f
export E5_SAMPLE_SUBDIR=samples
export E5_UNSCREENED=0
export E5_FIBERS=1:600
export E5_BUILD_ORDER=asc
export E5_USE_SLURM=auto
# 2 BLAS threads per worker is REQUIRED, not a tuning choice: 1 thread x 96 workers is
# memory-infeasible on genoa at 19 GB/worker.
export E5_BLAS_THREADS=2

# Workers now come from the Slurm allocation (one task per worker), so instead of
# computing a count we VERIFY the allocation is memory-safe at the measured 19 GB/worker
# peak and refuse to start if it is not.
PEAK_GB_PER_WORKER=19
TPN=${SLURM_NTASKS_PER_NODE:-48}
if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "${SLURM_MEM_PER_NODE:-0}" -gt 0 ] 2>/dev/null; then
    MEM_GB=$(( SLURM_MEM_PER_NODE / 1024 ))
else
    MEM_GB=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) / 1024 / 1024 ))
fi
SAFE_TPN=$(( (MEM_GB * 80 / 100) / PEAK_GB_PER_WORKER ))
echo "node memory ${MEM_GB} GB; ${TPN} workers/node x ${PEAK_GB_PER_WORKER} GB = $(( TPN * PEAK_GB_PER_WORKER )) GB (memory-safe max ${SAFE_TPN})"
if [ "$TPN" -gt "$SAFE_TPN" ]; then
    echo "ERROR: --ntasks-per-node=$TPN exceeds the memory-safe maximum $SAFE_TPN for a ${MEM_GB} GB node"
    echo "       at the measured ${PEAK_GB_PER_WORKER} GB/worker peak. Lower --ntasks-per-node or use genoa."
    exit 3
fi
# E5_NWORKERS is only the OFF-Slurm fallback; under Slurm the allocation decides.
export E5_NWORKERS=$TPN

# threshold policy: first CLI arg wins, else pre-set E5_THRESH_POLICY, else fail loudly
export E5_THRESH_POLICY="${1:-${E5_THRESH_POLICY:-}}"
if [ -z "$E5_THRESH_POLICY" ] || [ "$E5_THRESH_POLICY" = "legacy" ]; then
    echo "ERROR: pass a threshold policy, e.g. 'abs:35,8' | 'abs:150,40' | 'quantile:0.083' | 'off'."
    echo "       'legacy' is the inherited no-op (finding #35) and is refused here."
    exit 2
fi
export E5_BRIGHT_GUARD=error
# calibrated policies land at ~8% bright (DR17: 8.35% APO / 8.15% LCO); 4-15% is loose
# enough for per-fiber scatter and tight enough to catch a no-op or a runaway mask
export E5_BRIGHT_FRAC_LO=0.04
export E5_BRIGHT_FRAC_HI=0.15
echo "E5_THRESH_POLICY=$E5_THRESH_POLICY E5_BRIGHT_GUARD=$E5_BRIGHT_GUARD E5_NWORKERS=$E5_NWORKERS"

print_elapsed_time "assert skyCont is threshold-policy invariant (gates the symlink reuse)"
julia +$julia_version --project=$base_dir $base_dir/scripts/prior_build/e5_assert_skycont_invariant.jl

print_elapsed_time "rebuild skyLines under policy $E5_THRESH_POLICY"
julia +$julia_version --project=$base_dir $base_dir/scripts/prior_build/e5_sky_run.jl --stage=build

print_elapsed_time "Job Completed"
