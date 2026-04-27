#!/bin/bash
#SBATCH --partition=cca
#SBATCH --nodes=1
#SBATCH --constraint="[genoa|icelake|rome]"

#SBATCH --time=4-00:00
#SBATCH --job-name=ar_priors
#SBATCH --output=slurm_logs/%x_%j.out
# ------------------------------------------------------------------------------
SLURM_NTASKS=$(($SLURM_CPUS_ON_NODE * $SLURM_NNODES))
export SLURM_NTASKS
# Print all SLURM environment variables
env | grep SLURM | while read -r line; do
    echo "$line"
done

# exit the script immediately if a command fails
set -e
set -o pipefail

if [ -n "$SLURM_JOB_NODELIST" ]; then
    hostname
    echo $SLURM_JOB_NODELIST
else
    hostname
fi
echo "running from $(pwd)"

if [ -n "${SLURM_JOB_ID:-}" ] ; then
    script_path=$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')
else
    script_path=$(realpath "$0")
fi
base_dir="$(dirname "$(dirname "$(dirname "$script_path")")")"
echo "base_dir: $base_dir"

julia_version="1.11.0" # 1.11.6
juliaup add $julia_version

# nice function to print the elapsed time at each step
print_elapsed_time() {
    current_seconds=$SECONDS
    if [ -n "$LAST_TIME" ]; then
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

# ----- sequential -----
# julia +1.10.0 sample_sky.jl # 3k core-h, 7.7h on 6 nodes, 10% CPU usage (2 corrupted skySpec_tellDiv_ files had to be manually rm-ed, switch to pout exit code model next time)
# julia +1.10.0 build_skyCont.jl # 672 core-h, 1.75h on 6 nodes, 100% cpu usage [OOM possible with Krylov]
# julia +1.10.0 build_skyLines.jl # 2.7k core-h, 7h on 6 nodes, 100% cpu usage [never use Krylov]
# julia +1.10.0 sample_Tfun.jl # ~2.3k, 6h on 6 nodes, 2-20% cpu usage (1 restart, no manual intervention, switch to pout exit code model next time)
julia +$julia_version --project=$base_dir $base_dir/src/prior_build/sample_starCont.jl # 18 core-h, 0.6 h on 1 node (30 cores)
# julia +1.10.0 build_starCont.jl # 346 core-h, 0.9 h on 6 nodes, 100% cpu usage [OOM possible with Krylov]
# ----- sequential -----
# julia +1.10.0 sample_Korg.jl # 966.4 core-h, 2.5h on 6 nodes, 34.8 core-s/spec, 100% cpu usage
# julia +1.10.0 build_starLines.jl # 40 core-h, 40 min on 1 node, 50% cpu usage
# ----- indep -----
# julia +1.10.0 build_DIB.jl # 145 core-h, 2.3 h on 1 node, 100% cpu usage
# ----- second pass -----
# julia +1.10.2 build_starLines_dd.jl # 300 core-h, 2.8 h on 4 nodes, 35% cpu usage APO/100% cpu usage LCO (core-h are low, but problem is memory usage is high)

print_elapsed_time "Job Completed"