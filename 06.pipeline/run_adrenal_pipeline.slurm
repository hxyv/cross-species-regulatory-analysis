#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1-00:00:00
#SBATCH --job-name=adrenal_pipeline
#SBATCH -A bio230007p
#SBATCH --output=06.pipeline/run_adrenal_pipeline.%j.out
#SBATCH --error=06.pipeline/run_adrenal_pipeline.%j.err
#SBATCH --export=ALL

repo_root="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${repo_root}"

bash 06.pipeline/run_adrenal_pipeline.sh "$@"
