#!/bin/bash -l
#
#SBATCH --job-name="trace_mpi_omp_4x2"
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --time=00:04:00 # hh:mm:ss
#SBATCH --output=trace_mpi_omp_4x2.out
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

module load Scalasca
module load OpenMPI

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

scalasca -analyze -q -t -f scorep.filter srun ../../mpi_omp/fdtd param_3d.txt
