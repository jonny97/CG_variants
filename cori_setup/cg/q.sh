#!/bin/bash
#SBATCH -C haswell
#SBATCH -p debug            # change this option for non-debug runs
#SBATCH -N 1                # you'll never need more than 1 node for the serial code
#SBATCH -t 00:20:00         # adjust the amount of time as necessary
#SBATCH --qos=debug
#SBATCH --time=5
#SBATCH --nodes=2
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell

srun check-mpi.intel.cori
