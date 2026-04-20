#!/bin/bash
#SBATCH -N 2
#SBATCH -p RM-shared
#SBATCH -t 1-00:00:00
#SBATCH --job-name=adrenal_pipeline
#SBATCH -A bio230007p
#SBATCH --output=07.pipeline/run_adrenal_pipeline.%j.out
#SBATCH --error=07.pipeline/run_adrenal_pipeline.%j.err
#SBATCH --export=ALL

cd /ocean/projects/bio230007p/xhu15/cross-species-regulatory-analysis

bash 07.pipeline/run_adrenal_pipeline.sh "$@"
