#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 02:00:00
#SBATCH --job-name=adrenal_bedtools_intersect
#SBATCH -A bio230007p
#SBATCH --output=02.mapping/run_adrenal_bedtools_intersection.%j.out
#SBATCH --error=02.mapping/run_adrenal_bedtools_intersection.%j.err
#SBATCH --export=ALL

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

bash 02.mapping/run_adrenal_bedtools_intersection.sh
