#!/bin/bash
#SBATCH --job-name="default_run_custom_young"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00
#SBATCH --output=run_default_custom_young.out

srun ../../default/fdtd param_custom_young.txt
