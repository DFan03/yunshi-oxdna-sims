#!/bin/sh
#
#SBATCH --job-name="K1IVT_md"
#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=research-as-bn
#SBATCH --gpus-per-task=1
srun bash cmds_md
