#!/bin/bash
#SBATCH --job-name="simple/3d"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00
#SBATCH --output=run_3d.out

module load OpenMPI
srun ../../fdtd param_3d.txt
