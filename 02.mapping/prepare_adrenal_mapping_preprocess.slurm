#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 1-00:00:00
#SBATCH --job-name=adrenal_mapping_preprocess
#SBATCH -A bio230007p
#SBATCH --output=02.mapping/prepare_adrenal_mapping_preprocess.%j.out
#SBATCH --error=02.mapping/prepare_adrenal_mapping_preprocess.%j.err
#SBATCH --export=ALL

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

bash 02.mapping/prepare_adrenal_mapping_preprocess.sh
