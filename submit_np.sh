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

julia +1.11.6 --project="./" pipeline.jl --redux_base $reduxBase --almanac_file $almanacFile

julia +1.11.6 --project="./" workup.jl --outdir "../outdir/arMADGICS/raw/"

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"