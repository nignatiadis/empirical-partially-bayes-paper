#!/bin/bash

#SBATCH --job-name=noiserdd
#SBATCH --array=1-18
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20000


# Execute the Julia script, passing the SLURM_ARRAY_TASK_ID as an argument
julia simulations.jl $SLURM_ARRAY_TASK_ID

