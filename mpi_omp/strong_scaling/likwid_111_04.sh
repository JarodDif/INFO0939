#!/bin/bash -l
#SBATCH --job-name="likwid_ss_11104"
#SBATCH --output="likwid_11104.out"
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00 # hh:mm:ss
#SBATCH --mem-per-cpu=10240
#SBATCH --partition=hmem
#SBATCH --exclusive

module load OpenMPI
module load likwid

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

likwid-mpirun --mpi slurm -np ${SLURM_NTASKS} -t ${SLURM_CPUS_PER_TASK} -g CACHE ../fdtd param_3d_100_100_100.txt