#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mapping_dir="${repo_root}/03.mapping"

human_source=${HUMAN_PEAKS:-/ocean/projects/bio230007p/ikaplow/HumanAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz}
mouse_source=${MOUSE_PEAKS:-/ocean/projects/bio230007p/ikaplow/MouseAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz}

human_link=${mapping_dir}/human_adrenal_idr_optimal.narrowPeak.gz
mouse_link=${mapping_dir}/mouse_adrenal_idr_optimal.narrowPeak.gz

human_output=${mapping_dir}/human_adrenal_idr_optimal.mapping_preprocess.bed
mouse_output=${mapping_dir}/mouse_adrenal_idr_optimal.mapping_preprocess.bed

ln -sfn ${human_source} ${human_link}
ln -sfn ${mouse_source} ${mouse_link}

gzip -dc ${human_link} \
  | awk 'BEGIN { OFS="\t" }
    {
      peak_id = "human_adrenal_gland_idr_optimal_peak_" sprintf("%07d", NR)
      name = ($4 == "." || $4 == "") ? peak_id : $4
      print $1, $2, $3, name, $5, $6, $7, $8, $9, $10
    }' \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
  > ${human_output}

gzip -f ${human_output}

gzip -dc ${mouse_link} \
  | awk 'BEGIN { OFS="\t" }
    {
      peak_id = "mouse_adrenal_gland_idr_optimal_peak_" sprintf("%07d", NR)
      name = ($4 == "." || $4 == "") ? peak_id : $4
      print $1, $2, $3, name, $5, $6, $7, $8, $9, $10
    }' \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
  > ${mouse_output}

gzip -f ${mouse_output}
