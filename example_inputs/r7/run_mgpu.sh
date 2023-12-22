#!/bin/bash
#SBATCH --job-name="mgpu_r7"
#SBATCH --output="run_mgpu.out"
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=1:00:00
#SBATCH --account=ulghpsc

module load OpenMPI/4.1.5-Clang-16.0.6-GFortran-11.3.0-CUDA-11.7.0

export OMP_TARGET_OFFLOAD=MANDATORY

srun ../../mpi_gpu_subdomain_save/fdtd param_r7.txt
