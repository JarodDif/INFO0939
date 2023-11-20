#!/bin/bash
#SBATCH --job-name="audio"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00
#SBATCH --output=run_audio.out

module load OpenMPI
srun ../../fdtd param_audio.txt
