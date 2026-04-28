#!/bin/bash
#SBATCH -J rgreat_task3
#SBATCH -p RM
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -o logs/out.%j.txt
#SBATCH -e logs/err.%j.txt

cd $SLURM_SUBMIT_DIR
source ~/.bashrc
conda activate base

export LD_LIBRARY_PATH=/jet/home/dommalap/anaconda3/lib:$LD_LIBRARY_PATH

Rscript run_rgreat.R
