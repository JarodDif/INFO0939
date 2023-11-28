#!/bin/bash -l
#
#SBATCH --job-name="profile_omp_1x8"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:04:00 # hh:mm:ss
#SBATCH --output=profile_omp_1x8.out
#SBATCH --mem-per-cpu=1024 # 1GB
#SBATCH --partition=batch

module load Scalasca
module load OpenMPI

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

scalasca -analyze -q -t -f scorep.filter srun ../../omp/fdtd param_3d.txt
