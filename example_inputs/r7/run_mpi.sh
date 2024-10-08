#!/bin/bash
#SBATCH --job-name="mpi_r7"
#SBATCH --output="run_mpi.out"
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --time=01:00:00

module load OpenMPI
srun ../../mpi/fdtd param_r7.txt
