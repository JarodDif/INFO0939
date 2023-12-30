#!/bin/bash
#SBATCH --job-name="ws_04"
#SBATCH --output="run_04.out"
#SBATCH --ntasks=4
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=batch

module load OpenMPI

srun ../fdtd param_3d_200_200_100.txt