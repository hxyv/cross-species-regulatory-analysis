#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mapping_dir="${repo_root}/03.mapping"
hal_file=${HAL_FILE:-/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal}
hal_bin="/ocean/projects/bio230007p/xhu15/tools/hal/bin"
halper_repo="/ocean/projects/bio230007p/xhu15/tools/halLiftover-postprocessing"
halper_script="${halper_repo}/halper_map_peak_orthologs.sh"

export PATH=${hal_bin}:${PATH}
export PYTHONPATH=${halper_repo}:${PYTHONPATH:-}

${hal_bin}/halStats --genomes ${hal_file}

bash ${halper_script} \
  -b ${mapping_dir}/human_adrenal_idr_optimal.mapping_preprocess.bed.gz \
  -o ${mapping_dir} \
  -s Human \
  -t Mouse \
  -n human_adrenal_idr_optimal \
  -c ${hal_file}
