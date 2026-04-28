#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
config_file="${PIPELINE_CONFIG:-${repo_root}/pipeline.conf}"

run_rgreat=1
run_homer=1

while [[ $# -gt 0 ]]; do
  case $1 in
    --config)
      config_file=$2
      shift 2
      ;;
    --skip-rgreat)
      run_rgreat=0
      shift
      ;;
    --skip-homer)
      run_homer=0
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -f ${config_file} ]]; then
  source "${config_file}"
fi

export PIPELINE_CONFIG=${config_file}
export HUMAN_PEAKS=${HUMAN_PEAKS:-/ocean/projects/bio230007p/ikaplow/HumanAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz}
export MOUSE_PEAKS=${MOUSE_PEAKS:-/ocean/projects/bio230007p/ikaplow/MouseAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz}
export HAL_FILE=${HAL_FILE:-/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal}
export HAL_BIN=${HAL_BIN:-}
export HALPER_DIR=${HALPER_DIR:-}
export HUMAN_TSS=${HUMAN_TSS:-/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed}
export MOUSE_TSS=${MOUSE_TSS:-/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed}
export HOMER_DIR=${HOMER_DIR:-}

require_file() {
  label=$1
  path=$2
  if [[ ! -f ${path} ]]; then
    echo "Missing ${label}: ${path}" >&2
    exit 1
  fi
}

require_dir() {
  label=$1
  path=$2
  if [[ ! -d ${path} ]]; then
    echo "Missing ${label}: ${path}" >&2
    exit 1
  fi
}

require_command() {
  tool=$1
  if ! command -v ${tool} >/dev/null 2>&1; then
    echo "Missing required command: ${tool}" >&2
    exit 1
  fi
}

require_python_module() {
  module_name=$1
  python3 -c "import ${module_name}" >/dev/null 2>&1 || {
    echo "Missing required Python module: ${module_name}" >&2
    exit 1
  }
}

decompress_to_file() {
  in_file=$1
  out_file=$2
  gzip -dc ${in_file} > ${out_file}
}

sort_bed_file() {
  in_file=$1
  out_file=$2
  awk 'BEGIN{FS="[[:space:]]+"; OFS="\t"} NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3 {print}' ${in_file} \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
    > ${out_file}
}

require_file "human peak file" ${HUMAN_PEAKS}
require_file "mouse peak file" ${MOUSE_PEAKS}
require_file "HAL alignment" ${HAL_FILE}
require_file "human TSS BED" ${HUMAN_TSS}
require_file "mouse TSS BED" ${MOUSE_TSS}
require_dir "HAL binary directory" ${HAL_BIN}
require_dir "HALPER directory" ${HALPER_DIR}

# Use the default Bridges bedtools module if bedtools is not already on PATH.
if ! command -v bedtools >/dev/null 2>&1; then
  module load bedtools/2.30.0 >/dev/null 2>&1 || true
fi

require_command bash
require_command awk
require_command sort
require_command cut
require_command wc
require_command gzip
require_command python3
require_command Rscript
require_command bedtools
require_python_module pandas
require_python_module matplotlib

if [[ ! -x ${HAL_BIN}/halStats ]]; then
  echo "Missing HAL executable: ${HAL_BIN}/halStats" >&2
  exit 1
fi

if [[ ! -x ${HAL_BIN}/halLiftover ]]; then
  echo "Missing HAL executable: ${HAL_BIN}/halLiftover" >&2
  exit 1
fi

if [[ ! -f ${HALPER_DIR}/halper_map_peak_orthologs.sh ]]; then
  echo "Missing HALPER script: ${HALPER_DIR}/halper_map_peak_orthologs.sh" >&2
  exit 1
fi

mapping_dir=${repo_root}/02.mapping
bp_dir=${repo_root}/03.biological_processes
pe_dir=${repo_root}/04.promoter_enhancer
motif_dir=${repo_root}/05.motifs
motif_input_dir=${motif_dir}/downstream

# Tasks 3 to 5 depend on the shared/species-specific peak sets created in Task 2.
echo "Task 2: map adrenal peaks and define shared/species-specific sets"
bash ${mapping_dir}/prepare_adrenal_mapping_preprocess.sh
bash ${mapping_dir}/run_adrenal_halper_mapping.sh
bash ${mapping_dir}/run_adrenal_bedtools_intersection.sh
bash ${mapping_dir}/recover_human_coordinate_peak_sets.sh

echo "Stage Task 3 inputs"
decompress_to_file ${mapping_dir}/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed.gz ${bp_dir}/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed
decompress_to_file ${mapping_dir}/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed.gz ${bp_dir}/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed

if [[ ${run_rgreat} -eq 1 ]]; then
  echo "Task 3: run rGREAT biological process analysis"
  cd ${bp_dir}
  Rscript run_rgreat.R
  python top10_GO_BP_Plot.py
fi

echo "Stage Task 4 inputs"
decompress_to_file ${mapping_dir}/human_adrenal_idr_optimal.mapping_preprocess.bed.gz ${pe_dir}/human_adrenal_idr_optimal.narrowPeak
decompress_to_file ${mapping_dir}/mouse_adrenal_idr_optimal.mapping_preprocess.bed.gz ${pe_dir}/mouse_adrenal_idr_optimal.narrowPeak
decompress_to_file ${mapping_dir}/human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed.gz ${pe_dir}/human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed
decompress_to_file ${mapping_dir}/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed.gz ${pe_dir}/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed
decompress_to_file ${mapping_dir}/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed.gz ${pe_dir}/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed
decompress_to_file ${mapping_dir}/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed.gz ${pe_dir}/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed
decompress_to_file ${mapping_dir}/human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed.gz ${pe_dir}/human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed
decompress_to_file ${mapping_dir}/human_adrenal_idr_optimal.shared.original_human_coordinates.bed.gz ${pe_dir}/human_adrenal_idr_optimal.shared.original_human_coordinates.bed

awk 'BEGIN{OFS="\t"} {s=$2-2000; if(s<0) s=0; print $1,s,$3+2000}' ${HUMAN_TSS} | LC_ALL=C sort -k1,1 -k2,2n > ${pe_dir}/human_tss_2kb.bed
awk 'BEGIN{OFS="\t"} {s=$2-2000; if(s<0) s=0; print $1,s,$3+2000}' ${MOUSE_TSS} | LC_ALL=C sort -k1,1 -k2,2n > ${pe_dir}/mouse_tss_2kb.bed

echo "Task 4: classify promoter and enhancer peaks"
cd ${pe_dir}
bash classifyingpeaks.sh

echo "Stage Task 5 inputs"
sort_bed_file ${pe_dir}/humanenhancers.bed ${motif_input_dir}/humanenhancers.sorted.bed
sort_bed_file ${pe_dir}/humanpromoters.bed ${motif_input_dir}/humanpromoters.sorted.bed
sort_bed_file ${pe_dir}/mouse_enhancers.bed ${motif_input_dir}/mouse_enhancers.sorted.bed
sort_bed_file ${pe_dir}/mouse_promoters.bed ${motif_input_dir}/mouse_promoters.sorted.bed
sort_bed_file ${pe_dir}/human_specific_enhancers_hg38.bed ${motif_input_dir}/human_specific_enhancers_hg38.sorted.bed
sort_bed_file ${pe_dir}/human_specific_promoters_hg38.bed ${motif_input_dir}/human_specific_promoters_hg38.sorted.bed
sort_bed_file ${pe_dir}/conserved_enhancers_hg38.bed ${motif_input_dir}/conserved_enhancers_hg38.sorted.bed
sort_bed_file ${pe_dir}/conserved_promoters_hg38.bed ${motif_input_dir}/conserved_promoters_hg38.sorted.bed
sort_bed_file ${pe_dir}/mousespecific_enhancers.bed ${motif_input_dir}/mousespecific_enhancers.sorted.bed
sort_bed_file ${pe_dir}/mousespecific_promoters.bed ${motif_input_dir}/mousespecific_promoters.sorted.bed
sort_bed_file ${pe_dir}/mouse_conserved_enhancers.bed ${motif_input_dir}/mouse_conserved_enhancers.sorted.bed
sort_bed_file ${pe_dir}/mouse_conserved_promoters.bed ${motif_input_dir}/mouse_conserved_promoters.sorted.bed

if [[ ${run_homer} -eq 1 ]]; then
  require_dir "HOMER directory" ${HOMER_DIR}
  ln -sfn ${HOMER_DIR} ${motif_dir}/homer

  echo "Task 5: run HOMER motif analyses"
  cd ${motif_input_dir}
  bash build_centered_beds.sh
  bash run_all.sh
  python3 summarize_known_motifs.py
fi

echo "Pipeline finished"
