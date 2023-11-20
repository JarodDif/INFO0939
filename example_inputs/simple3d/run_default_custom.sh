#!/bin/bash
#SBATCH --job-name="default_run_custom"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00
#SBATCH --output=run_default_custom.out

srun ../../default/fdtd param_custom.txt
