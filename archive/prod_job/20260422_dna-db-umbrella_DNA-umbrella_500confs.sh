#!/bin/sh
#
#SBATCH --job-name="prod"
#SBATCH --partition=gpu-a100-small
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus-per-task=1
unset CONDA_SHLVL
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate oxdna
srun bash cmds
