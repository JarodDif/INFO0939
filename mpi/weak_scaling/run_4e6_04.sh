#!/bin/bash
#SBATCH --job-name="ws_04"
#SBATCH --output="run_4e6_04.out"
#SBATCH --ntasks=4
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

srun ../fdtd param_3d_400_200_200.txt