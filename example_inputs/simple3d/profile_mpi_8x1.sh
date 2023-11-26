#!/bin/bash -l
#
#SBATCH --job-name="Scalasca profile"
#SBATCH --ntasks=8
#
#SBATCH --time=00:04:00 # hh:mm:ss
#SBATCH --output=run_mpi_8x1.out
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

# Loading Scalasca will also load OpenMPI
module load Scalasca

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

scalasca -analyze srun ./myapp_instr