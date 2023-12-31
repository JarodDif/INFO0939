#!/bin/bash
#SBATCH --job-name="ss_22216"
#SBATCH --output="run_22216.out"
#SBATCH --ntasks=16
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

srun ../fdtd param_3d_200_200_200.txt