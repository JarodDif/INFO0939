#!/bin/bash
#SBATCH --job-name="mpi_run_custom"
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00
#SBATCH --output=run_mpi_custom.out

module load OpenMPI
srun ../../mpi/fdtd param_custom.txt
