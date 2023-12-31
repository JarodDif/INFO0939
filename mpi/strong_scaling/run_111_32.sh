#!/bin/bash
#SBATCH --job-name="ss_11132"
#SBATCH --output="run_11132.out"
#SBATCH --ntasks=32
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

srun ../fdtd param_3d_100_100_100.txt