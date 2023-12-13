#!/bin/bash
#SBATCH --job-name="mpi_point"
#SBATCH --output="run_mpi.out"
#SBATCH --ntasks=8
#SBATCH --time=00:04:00 # hh:mm:ss
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

module load OpenMPI

srun ../../mpi/fdtd param_point_3d.txt