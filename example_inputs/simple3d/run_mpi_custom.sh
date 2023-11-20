#!/bin/bash
#SBATCH --job-name="simple/custom"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00
#SBATCH --output=run_custom.out

module load OpenMPI
srun ../../mpi/fdtd param_custom.txt
