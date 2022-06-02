#!/bin/bash
#SBATCH --partition=kamiak   # Partition (like a queue in PBS)
#SBATCH --job-name=dturtle_bwp # Job Name
#SBATCH --output=R_%j.out
#SBATCH --error=R_%j.err
#SBATCH --time=2-00:00:00    # Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1            # Node count required for the job
#SBATCH --ntasks-per-node=1  # Number of tasks to be launched per Node
#SBATCH --cpus-per-task=14    # Number of threads per task (OMP threads)

module load r/4.1.0

Rscript AllTissues_dTurtle_01.05.21.R
