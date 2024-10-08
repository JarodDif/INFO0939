#!/bin/bash
#SBATCH --job-name="gpu_point"
#SBATCH --output="run_gpu.out"
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus=1
#SBATCH --time=15:00
#SBATCH --account=ulghpsc

module load Clang/16.0.6-GCCcore-11.3.0-CUDA-11.7.0

export OMP_TARGET_OFFLOAD=MANDATORY

srun ../../gpu/fdtd param_point_3d.txt
