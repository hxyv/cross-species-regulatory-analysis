#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 12:00:00
#SBATCH --job-name=adrenal_halper_mapping
#SBATCH -A bio230007p
#SBATCH --output=02.mapping/run_adrenal_halper_mapping.%j.out
#SBATCH --error=02.mapping/run_adrenal_halper_mapping.%j.err
#SBATCH --export=ALL

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

bash 02.mapping/run_adrenal_halper_mapping.sh
