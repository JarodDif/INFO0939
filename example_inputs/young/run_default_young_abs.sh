#!/bin/bash
#SBATCH --job-name="default_run_young_abs"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00
#SBATCH --output=run_default_young_abs.out

srun ../../default/fdtd param_young_abs.txt
