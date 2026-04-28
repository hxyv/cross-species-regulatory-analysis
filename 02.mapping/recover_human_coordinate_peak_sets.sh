#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mapping_dir="${repo_root}/02.mapping"
config_file=${PIPELINE_CONFIG:-${repo_root}/pipeline.conf}

if [[ -f ${config_file} ]]; then
  source "${config_file}"
fi

human_peaks=${mapping_dir}/human_adrenal_idr_optimal.mapping_preprocess.bed.gz
human_specific_mapped=${mapping_dir}/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed.gz
human_shared_mapped=${mapping_dir}/human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed.gz

human_specific_names=${mapping_dir}/human_specific_peak_names.txt
human_shared_names=${mapping_dir}/human_shared_peak_names.txt

human_specific_original=${mapping_dir}/human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed
human_shared_original=${mapping_dir}/human_adrenal_idr_optimal.shared.original_human_coordinates.bed
summary=${mapping_dir}/human_coordinate_recovered_peak_sets_summary.tsv

# Recover original human-coordinate rows by matching mapped peak names back to
# the human preprocessing file. This keeps downstream human analyses in hg38.
gzip -dc ${human_specific_mapped} | cut -f4 | sort -u > ${human_specific_names}
gzip -dc ${human_shared_mapped} | cut -f4 | sort -u > ${human_shared_names}

gzip -dc ${human_peaks} \
  | awk 'BEGIN { OFS="\t" } NR==FNR { keep[$1]=1; next } ($4 in keep)' ${human_specific_names} - \
  > ${human_specific_original}

gzip -dc ${human_peaks} \
  | awk 'BEGIN { OFS="\t" } NR==FNR { keep[$1]=1; next } ($4 in keep)' ${human_shared_names} - \
  > ${human_shared_original}

gzip -f ${human_specific_original}
gzip -f ${human_shared_original}

printf "set\tcount\n" > ${summary}
printf "human_specific_mapped_names\t" >> ${summary}
wc -l < ${human_specific_names} >> ${summary}
printf "human_specific_original_human_coordinates\t" >> ${summary}
gzip -dc ${human_specific_original}.gz | wc -l >> ${summary}
printf "human_shared_mapped_names\t" >> ${summary}
wc -l < ${human_shared_names} >> ${summary}
printf "human_shared_original_human_coordinates\t" >> ${summary}
gzip -dc ${human_shared_original}.gz | wc -l >> ${summary}
