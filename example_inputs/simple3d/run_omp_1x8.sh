#!/bin/bash
#SBATCH --job-name="omp_1x8"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:04:00 # hh:mm:ss
#SBATCH --output=run_omp_1x8.out
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

module load OpenMPI
export OMP_NUM_THREADS=8

srun  ../../omp/fdtd param_3d.txt