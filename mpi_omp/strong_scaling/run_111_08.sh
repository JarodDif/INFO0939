#!/bin/bash
#SBATCH --job-name="ss_11108"
#SBATCH --output="run_11108.out"
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

srun ../fdtd param_3d_100_100_100.txt