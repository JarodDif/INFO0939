#!/bin/bash
#SBATCH --job-name="ws_32"
#SBATCH --output="run_8e6_32.out"
#SBATCH --ntasks=32
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

srun ../fdtd param_3d_800_800_400.txt