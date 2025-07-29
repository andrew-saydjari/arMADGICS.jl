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

julia +1.11.0 --project="./" pipeline.jl

# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"