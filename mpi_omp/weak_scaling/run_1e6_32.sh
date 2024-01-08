#!/bin/bash
#SBATCH --job-name="ws_32"
#SBATCH --output="run_1e6_32.out"
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

srun ../fdtd param_3d_400_400_200.txt