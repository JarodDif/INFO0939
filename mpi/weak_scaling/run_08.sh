#!/bin/bash
#SBATCH --job-name="ws_08"
#SBATCH --output="run_08.out"
#SBATCH --ntasks=8
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

srun ../fdtd param_3d_200_200_200.txt