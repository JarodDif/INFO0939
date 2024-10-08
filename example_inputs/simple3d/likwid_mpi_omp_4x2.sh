#!/bin/bash -l
#
#SBATCH --job-name="likwid_mpi_omp_4x2"
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --time=00:10:00 # hh:mm:ss
#SBATCH --output=likwid_mpi_omp_4x2.out
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

module load likwid
module load OpenMPI

unset OMP_NUM_THREADS

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

likwid-mpirun --mpi slurm -np ${SLURM_NTASKS} -t ${SLURM_CPUS_PER_TASK} -g CACHE ../../mpi_omp/fdtd param_3d.txt
