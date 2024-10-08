#!/bin/bash
#SBATCH --job-name="ws_02"
#SBATCH --output="run_8e6_02.out"
#SBATCH --ntasks=2
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

srun ../fdtd param_3d_400_200_200.txt