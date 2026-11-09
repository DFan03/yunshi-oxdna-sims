#!/bin/sh
#
#SBATCH --job-name="job_name"
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=research-as-bn

unset CONDA_SHLVL
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate oxdna

# Get the correct line from cmds
CMD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" 2state_cmds)

echo "Running task $SLURM_ARRAY_TASK_ID: $CMD"

srun bash -c "$CMD"


