#!/bin/bash
#SBATCH --job-name="test"
#SBATCH --output="run_mpi_aware.out"
#SBATCH --partition=debug-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1:00
#SBATCH --account=ulghpsc

module load OpenMPI/4.1.5-Clang-16.0.6-GFortran-11.3.0-CUDA-11.7.0

export OMP_TARGET_OFFLOAD=MANDATORY

srun ./mpi_aware