#!/bin/bash
#SBATCH --account=sdss-np
#SBATCH --partition=sdss-shared-np
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64

#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=96:00:00
#SBATCH --job-name=arMADGICS
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------

echo $SLURM_JOB_NODELIST

almanac_name=$1
reduxBase="/uufs/chpc.utah.edu/common/home/u6039752/scratch1/sandbox51/airflow-ApogeeReduction.jl/daily/outdir"
almanacFile=${reduxBase}/almanac/${almanac_name}.h5

# Warm the precompile cache ONCE, serially, before either driver spawns workers. One node
# here, so the cross-node pidfile storm measured in prior-build job 6995969 does not apply
# in full -- but `@everywhere` still fans out to 64 workers on this node plus the head, and
# a cold ApogeeReduction is ~92 s of redundant work per worker. Warming also front-loads
# workup.jl's set, and it verifies on the head that every package the workers need actually
# loads. Rationale and measured numbers: scripts/lib_julia_warm.sh.
# The version argument MUST match the `julia +1.11.0` calls below: precompile caches are
# keyed on the exact Julia build.
source "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/scripts/lib_julia_warm.sh"
warm_julia_precompile "1.11.0" "./" "./pipeline.jl" "./workup.jl"

julia +1.11.0 --project="./" pipeline.jl --redux_base $reduxBase --almanac_file $almanacFile

julia +1.11.0 --project="./" workup.jl --outdir "../outdir/arMADGICS/raw/"



# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"