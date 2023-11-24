#!/bin/bash
#SBATCH --job-name="job_custom"
#SBATCH --output="run_gpu_custom.out"
#SBATCH --partition=gpu
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --gpus=1
#SBATCH --time=15:00
#SBATCH --account=ulghpsc

module load Clang/16.0.6-GCCcore-11.3.0-CUDA-11.7.0

export OMP_TARGET_OFFLOAD=MANDATORY

../../gpu/fdtd param_custom.txt
