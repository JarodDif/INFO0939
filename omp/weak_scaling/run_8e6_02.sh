#!/bin/bash
#SBATCH --job-name="ws_02"
#SBATCH --output="run_8e6_02.out"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=20480
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

srun ../fdtd param_3d_400_200_200.txt