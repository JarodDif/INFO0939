#!/bin/bash
#SBATCH --job-name="mpi_run_big_r7.sh"
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --time=01:00:00
#SBATCH --output=run_mpi_custom.out

module load OpenMPI
srun ../../mpi/fdtd param_custom.txt
