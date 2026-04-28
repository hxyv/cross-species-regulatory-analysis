#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mapping_dir="${repo_root}/02.mapping"
config_file=${PIPELINE_CONFIG:-${repo_root}/pipeline.conf}

if [[ -f ${config_file} ]]; then
  source "${config_file}"
fi

hal_file=${HAL_FILE:-/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal}
hal_bin=${HAL_BIN:-}
halper_repo=${HALPER_DIR:-}
halper_script=${HALPER_SCRIPT:-${halper_repo}/halper_map_peak_orthologs.sh}

export PATH=${hal_bin}:${PATH}
export PYTHONPATH=${halper_repo}:${PYTHONPATH:-}

# Confirm the configured HAL install can see the alignment and supported genomes.
${hal_bin}/halStats --genomes ${hal_file}

bash ${halper_script} \
  -b ${mapping_dir}/human_adrenal_idr_optimal.mapping_preprocess.bed.gz \
  -o ${mapping_dir} \
  -s Human \
  -t Mouse \
  -n human_adrenal_idr_optimal \
  -c ${hal_file}
