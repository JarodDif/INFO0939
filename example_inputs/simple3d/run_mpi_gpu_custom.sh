#!/bin/bash
#SBATCH --job-name="mpi_gpu_custom"
#SBATCH --output="run_mpi_gpu_custom.out"
#SBATCH --partition=gpu
#SBATCH --ntasks=2
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1:00:00
#SBATCH --account=ulghpsc

module load OpenMPI/4.1.5-Clang-16.0.6-GFortran-11.3.0-CUDA-11.7.0

export OMP_TARGET_OFFLOAD=MANDATORY

srun ../../mpi_gpu_single_point/fdtd param_custom.txt
