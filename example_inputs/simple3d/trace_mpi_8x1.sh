#!/bin/bash -l
#
#SBATCH --job-name="trace_mpi_8x1"
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:04:00 # hh:mm:ss
#SBATCH --output=trace_mpi_8x1.out
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

module load Scalasca
module load OpenMPI

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

scalasca -analyze -q -t -f scorep.filter srun ../../mpi/fdtd_instr param_3d.txt
