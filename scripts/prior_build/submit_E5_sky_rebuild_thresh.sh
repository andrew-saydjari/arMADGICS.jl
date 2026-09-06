#!/bin/bash
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --constraint="[genoa|icelake|rome]"

#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --job-name=ar_E5_sky_rebuild
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
# E5 finding #35 REBUILD launcher: rebuild the skyLines priors under a chosen
# bright/faint threshold policy. Ready to fire for ANY of the four options in
# prior_outputs/sky_pass1/THRESHOLD_FINDING.md — no edits needed, pass the policy
# as the first argument (or set E5_THRESH_POLICY before sbatch; the value is
# re-read here at run time, so it is safe against sbatch's env snapshot):
#
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
# Expected wall time: skyLines only, ~30-50 min/fiber/worker => ~5-8 h for 600
# fibers at 32-48 workers on a cca node (skyCont adds nothing, it is symlinked).
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
export E5_BLAS_THREADS=2
export E5_NWORKERS=$(( SLURM_CPUS_ON_NODE / 2 < 48 ? SLURM_CPUS_ON_NODE / 2 : 48 ))

# threshold policy: first CLI arg wins, else pre-set E5_THRESH_POLICY, else fail loudly
export E5_THRESH_POLICY="${1:-${E5_THRESH_POLICY:-}}"
if [ -z "$E5_THRESH_POLICY" ] || [ "$E5_THRESH_POLICY" = "legacy" ]; then
    echo "ERROR: pass a threshold policy, e.g. 'abs:35,8' | 'abs:150,40' | 'quantile:0.083' | 'off'."
    echo "       'legacy' is the inherited no-op (finding #35) and is refused here."
    exit 2
fi
export E5_BRIGHT_GUARD=error
echo "E5_THRESH_POLICY=$E5_THRESH_POLICY E5_BRIGHT_GUARD=$E5_BRIGHT_GUARD E5_NWORKERS=$E5_NWORKERS"

print_elapsed_time "assert skyCont is threshold-policy invariant (gates the symlink reuse)"
julia +$julia_version --project=$base_dir $base_dir/scripts/prior_build/e5_assert_skycont_invariant.jl

print_elapsed_time "rebuild skyLines under policy $E5_THRESH_POLICY"
julia +$julia_version --project=$base_dir $base_dir/scripts/prior_build/e5_sky_run.jl --stage=build

print_elapsed_time "Job Completed"
