#!/bin/bash
#SBATCH --account=sdss-np-fast
#SBATCH --partition=sdss-np
#SBATCH --qos=sdss-np-req
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=64

#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=96:00:00
#SBATCH --job-name=arMADGICS
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------

echo $SLURM_JOB_NODELIST

# Warm the precompile cache ONCE, serially, before pipeline.jl spawns its worker pool. This
# is the launcher with the largest exposure in the repo: SlurmManager sizes the pool from
# SLURM_NTASKS = 8 nodes x 64 tasks, and `@everywhere` then executes on the head AND all 512
# workers at once, so every package first touched inside that block is resolved by 513
# processes simultaneously against one shared depot.
# MEASURED precedent (prior-build job 6995969, only 224 workers on 7 nodes): 219 of them
# precompiled ApogeeReduction independently and the phase cost 15 min 3 s; contention gets
# worse, not better, as workers are added. Rationale and numbers: scripts/lib_julia_warm.sh.
# The version argument MUST match the `julia +1.11.0` below -- precompile caches are keyed
# on the exact Julia build, so warming the wrong one is pure cost with no benefit.
source "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/scripts/lib_julia_warm.sh"
warm_julia_precompile "1.11.0" "./" "./pipeline.jl"

julia +1.11.0 --project="./" pipeline.jl

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"