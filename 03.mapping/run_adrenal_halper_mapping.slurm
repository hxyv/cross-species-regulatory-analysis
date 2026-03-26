#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 12:00:00
#SBATCH --job-name=adrenal_halper_mapping
#SBATCH -A bio230007p
#SBATCH --output=/ocean/projects/bio230007p/xhu15/cross-species-regulatory-analysis/03.mapping/run_adrenal_halper_mapping.%j.out
#SBATCH --error=/ocean/projects/bio230007p/xhu15/cross-species-regulatory-analysis/03.mapping/run_adrenal_halper_mapping.%j.err
#SBATCH --export=ALL

set -euo pipefail

cd /ocean/projects/bio230007p/xhu15/cross-species-regulatory-analysis

bash 03.mapping/run_adrenal_halper_mapping.sh
