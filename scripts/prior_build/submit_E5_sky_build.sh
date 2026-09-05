#!/bin/bash
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --constraint="[genoa|icelake|rome]"

#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --job-name=ar_E5_sky_build
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
# E5 pass-1 sky-prior BUILD (skyCont + skyLines faint/GSPICE) on the DR21
# testbed pooled samples — sbatch launcher, mirroring the apMADGICS original
# (src/prior_build/submit_priors_np.sh at Utah / this repo's
# scripts/prior_build/submit_priors.sh port). What changed vs the original:
#   - Utah original ran the whole prior sequence (sample -> build) with
#     hardcoded date paths; this script runs ONLY the E5 build stage against
#     the pre-assembled 2026_09_04 sky_pass1 samples via e5_sky_run.jl
#     --stage=build (env-configured I/O, no Pkg ops needed at run time).
#   - per-fiber resume: skip-if-built (checkpoint on the three output files)
#     PLUS atomic .claims/ dibs, so this job is safe to run CONCURRENTLY with
#     the on-node ccalin051 run (which ascends adjfib; this job descends via
#     E5_BUILD_ORDER=desc — the two meet in the middle and never double-build).
#   - --mem=0 (all node memory), as in the Utah original: GSPICE peaks
#     ~15 GB/worker at nspec~19k, so workers are sized to min(cores/2, 48)
#     (48 workers x 2 BLAS threads = 96 cores on genoa; memory-capped on rome).
# Historical timing (this corpus, w5-3435X, 2 BLAS threads): ~30-50 min/fiber;
# full 600 from scratch on genoa ~7-8 h; with the on-node run concurrent it
# only picks up the unbuilt tail. Stale claims after a crash: release with
# scripts/prior_build/e5_release_claims.sh <built_dir> <dead-hostname>.
#
# Submit (AKS) from a directory containing slurm_logs/, e.g.:
#   cd /mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/E5_run && \
#     sbatchAKS /mnt/home/asaydjari/gitcode/worktrees/arM-E5b/scripts/prior_build/submit_E5_sky_build.sh
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

# --- E5 build configuration (kept IN the script: sbatch snapshots env at
# --- submission but reads this file at run time, so edits here take effect) ---
export E5_OUT=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_04/prior_outputs/sky_pass1
export E5_REDUX=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/testbed_run/
export E5_ALMANAC=/mnt/ceph/users/sdssv/work/asaydjari/2026_09_03/testbed_prep/almanac_testbed_dr21.h5
export E5_STARCONT=$E5_OUT/starcont_fresh/APOGEE_starcont_svd_60_f
export E5_SAMPLE_SUBDIR=samples
export E5_UNSCREENED=0
export E5_FIBERS=1:600
export E5_BUILD_ORDER=desc        # on-node ccalin051 run ascends; this job descends
export E5_BLAS_THREADS=2
export E5_NWORKERS=$(( SLURM_CPUS_ON_NODE / 2 < 48 ? SLURM_CPUS_ON_NODE / 2 : 48 ))
echo "E5_NWORKERS=$E5_NWORKERS E5_BLAS_THREADS=$E5_BLAS_THREADS E5_BUILD_ORDER=$E5_BUILD_ORDER"

print_elapsed_time "build_skyCont + build_skyLines (E5 pass-1, claims+resume)"
# ~30-50 min/fiber/worker at 2 BLAS threads; 600 fibers from scratch ~7-8 h at 48 workers
julia +$julia_version --project=$base_dir $base_dir/scripts/prior_build/e5_sky_run.jl --stage=build

print_elapsed_time "Job Completed"
