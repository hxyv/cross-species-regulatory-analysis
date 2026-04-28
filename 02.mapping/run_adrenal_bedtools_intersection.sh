#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mapping_dir="${repo_root}/02.mapping"
config_file=${PIPELINE_CONFIG:-${repo_root}/pipeline.conf}

if [[ -f ${config_file} ]]; then
  source "${config_file}"
fi

mapped_human=${mapping_dir}/human_adrenal_idr_optimal.HumanToMouse.HALPER.narrowPeak.gz
native_mouse=${mapping_dir}/mouse_adrenal_idr_optimal.mapping_preprocess.bed.gz

shared_mapped=${mapping_dir}/human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed
human_mapped_only=${mapping_dir}/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed
shared_mouse=${mapping_dir}/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed
mouse_only=${mapping_dir}/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed
summary_tsv=${mapping_dir}/adrenal_human_to_mouse_intersection_summary.tsv

module load bedtools/2.30.0 >/dev/null 2>&1

run_intersect() {
  a_file=$1
  b_file=$2
  mode=$3
  out_file=$4

  bedtools intersect \
    -a ${a_file} \
    -b ${b_file} \
    ${mode} \
    > ${out_file}
}

write_count() {
  label=$1
  input_file=$2

  printf "%s\t" "${label}" >> ${summary_tsv}
  gzip -dc ${input_file} | wc -l >> ${summary_tsv}
}

# Query mapped human peaks against native mouse peaks to define shared and
# non-overlapping mapped-human intervals.
run_intersect ${mapped_human} ${native_mouse} -u ${shared_mapped}
run_intersect ${mapped_human} ${native_mouse} -v ${human_mapped_only}

# Reverse the query direction so we can count overlap from the native mouse side.
run_intersect ${native_mouse} ${mapped_human} -u ${shared_mouse}
run_intersect ${native_mouse} ${mapped_human} -v ${mouse_only}

gzip -f ${shared_mapped}
gzip -f ${human_mapped_only}
gzip -f ${shared_mouse}
gzip -f ${mouse_only}

printf "set\tcount\n" > ${summary_tsv}
write_count mapped_human_total ${mapped_human}
write_count native_mouse_total ${native_mouse}
write_count mapped_human_shared_with_mouse ${shared_mapped}.gz
write_count mapped_human_no_mouse_overlap ${human_mapped_only}.gz
write_count mouse_shared_with_human_mapped ${shared_mouse}.gz
write_count mouse_no_human_mapped_overlap ${mouse_only}.gz
