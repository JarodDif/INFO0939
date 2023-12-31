#!/bin/bash
#SBATCH --job-name="ws_64"
#SBATCH --output="run_1e6_64.out"
#SBATCH --ntasks=64
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

srun ../fdtd param_3d_400_400_400.txt