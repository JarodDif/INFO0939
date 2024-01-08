#!/bin/bash
#SBATCH --job-name="omp_r7"
#SBATCH --output="run_omp.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=20G
#SBATCH --time=01:00:00

module load OpenMPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

../../omp/fdtd param_r7.txt
