#!/bin/bash
#SBATCH --job-name=go_plot
#SBATCH --partition=RM-shared
#SBATCH --output=logs/plot_%j.out
#SBATCH --error=logs/plot_%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1900M

cd $SLURM_SUBMIT_DIR

python top10_GO_BP_Plot.py