#!/bin/sh
#SBATCH --job-name="prod"
#SBATCH --partition=gpu-a100-small
#SBATCH --time=4:00:00

#SBATCH --array=0-9   # 25 jobs (0–23)

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus-per-task=1
#SBATCH --account=research-as-bn
unset CONDA_SHLVL
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate oxdna

# Get the correct line from cmds
CMD=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" cmds)

echo "Running task $SLURM_ARRAY_TASK_ID: $CMD"

srun bash -c "$CMD"
