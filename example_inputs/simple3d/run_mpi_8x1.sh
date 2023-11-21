#!/bin/bash
#SBATCH --job-name="mpi_8x1"
#SBATCH --ntasks=8
#
#SBATCH --time=00:04:00 # hh:mm:ss
#SBATCH --output=run_mpi_8x1.out
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

module load OpenMPI

srun ../../mpi/fdtd param_3d.txt