#!/bin/bash
#SBATCH --job-name="mpi/simple/3d"
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00
#SBATCH --output=run_mpi.out

module load OpenMPI
srun ../../fdtd param_3d.txt