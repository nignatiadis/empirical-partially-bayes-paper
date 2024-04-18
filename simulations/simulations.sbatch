#!/bin/bash

#SBATCH --job-name=epb
#SBATCH --array=1-108
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20000
#SBATCH --account=pi-ignat
#SBATCH --partition=caslake

# Execute the Julia script, passing the SLURM_ARRAY_TASK_ID as an argument
julia simulations.jl $SLURM_ARRAY_TASK_ID
