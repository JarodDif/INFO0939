#!/bin/bash
#SBATCH --job-name="ss_11164"
#SBATCH --output="run_11164.out"
#SBATCH --ntasks=64
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

srun ../fdtd param_3d_100_100_100.txt