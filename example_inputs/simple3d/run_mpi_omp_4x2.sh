#!/bin/bash
#SBATCH --job-name="mpi_4x2"
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --time=00:04:00 # hh:mm:ss
#SBATCH --output=run_mpi_4x2.out
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

module load OpenMPI
export OMP_NUM_THREADS=2

srun  ../../mpi_omp/fdtd param_3d.txt