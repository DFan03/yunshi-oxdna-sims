#!/bin/sh
#
#SBATCH --job-name="K1IVT_mc"
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=research-as-bn

srun bash cmds_mc
