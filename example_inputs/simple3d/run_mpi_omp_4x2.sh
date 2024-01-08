#!/bin/bash
#SBATCH --job-name="mpi_omp_4x2"
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --time=00:04:00 # hh:mm:ss
#SBATCH --output=run_mpi_omp_4x2.out
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

module load OpenMPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

srun  ../../mpi_omp/fdtd param_3d.txt