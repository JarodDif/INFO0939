#!/bin/bash
#SBATCH --job-name="mpgu_point"
#SBATCH --output="run_mgpu.out"
#SBATCH --partition=gpu
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=15:00
#SBATCH --account=ulghpsc

module load OpenMPI/4.1.5-Clang-16.0.6-GFortran-11.3.0-CUDA-11.7.0

export OMP_TARGET_OFFLOAD=MANDATORY

srun ../../mpi_gpu_single_point/fdtd param_point_3d.txt
