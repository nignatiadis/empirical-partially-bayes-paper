#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=23:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --output=simulations.log
#SBATCH --account=stats




# execute script
julia simulations.jl
