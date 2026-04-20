#!/usr/bin/env bash
set -euo pipefail

# Example:
# bash 07.pipeline/run_adrenal_pipeline.sh --human /path/human.narrowPeak.gz --mouse /path/mouse.narrowPeak.gz --hal /path/10plusway-master.hal --human-tss /path/human_tss.bed --mouse-tss /path/mouse_tss.bed --homer /path/homer

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

human_peaks=/ocean/projects/bio230007p/ikaplow/HumanAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz
mouse_peaks=/ocean/projects/bio230007p/ikaplow/MouseAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz
hal_file=/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal
human_tss=/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed
mouse_tss=/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed
homer_dir=${HOMER_DIR:-}

run_rgreat=1
run_homer=1

while [[ $# -gt 0 ]]; do
  case $1 in
    --human)
      human_peaks=$2
      shift 2
      ;;
    --mouse)
      mouse_peaks=$2
      shift 2
      ;;
    --hal)
      hal_file=$2
      shift 2
      ;;
    --human-tss)
      human_tss=$2
      shift 2
      ;;
    --mouse-tss)
      mouse_tss=$2
      shift 2
      ;;
    --homer)
      homer_dir=$2
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
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

mapping_dir=${repo_root}/03.mapping
bp_dir=${repo_root}/04.biological_processes
pe_dir=${repo_root}/05.promoter_enhancer
motif_dir=${repo_root}/06.motifs
motif_input_dir=${motif_dir}/downstream

echo "Task 2: map adrenal peaks and define shared/species-specific sets"
HUMAN_PEAKS=${human_peaks} MOUSE_PEAKS=${mouse_peaks} bash ${mapping_dir}/prepare_adrenal_mapping_preprocess.sh
HAL_FILE=${hal_file} bash ${mapping_dir}/run_adrenal_halper_mapping.sh
bash ${mapping_dir}/run_adrenal_bedtools_intersection.sh
bash ${mapping_dir}/recover_human_coordinate_peak_sets.sh

echo "Stage Task 3 inputs"
gzip -dc ${mapping_dir}/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed.gz > ${bp_dir}/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed
gzip -dc ${mapping_dir}/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed.gz > ${bp_dir}/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed

if [[ ${run_rgreat} -eq 1 ]]; then
  echo "Task 3: run rGREAT biological process analysis"
  cd ${bp_dir}
  Rscript run_rgreat.R
  python top10_GO_BP_Plot.py
fi

echo "Stage Task 4 inputs"
gzip -dc ${mapping_dir}/human_adrenal_idr_optimal.mapping_preprocess.bed.gz > ${pe_dir}/human_adrenal_idr_optimal.narrowPeak
gzip -dc ${mapping_dir}/mouse_adrenal_idr_optimal.mapping_preprocess.bed.gz > ${pe_dir}/mouse_adrenal_idr_optimal.narrowPeak
gzip -dc ${mapping_dir}/human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed.gz > ${pe_dir}/human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed
gzip -dc ${mapping_dir}/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed.gz > ${pe_dir}/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed
gzip -dc ${mapping_dir}/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed.gz > ${pe_dir}/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed
gzip -dc ${mapping_dir}/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed.gz > ${pe_dir}/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed
gzip -dc ${mapping_dir}/human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed.gz > ${pe_dir}/human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed
gzip -dc ${mapping_dir}/human_adrenal_idr_optimal.shared.original_human_coordinates.bed.gz > ${pe_dir}/human_adrenal_idr_optimal.shared.original_human_coordinates.bed

awk 'BEGIN{OFS="\t"} {s=$2-2000; if(s<0) s=0; print $1,s,$3+2000}' ${human_tss} | LC_ALL=C sort -k1,1 -k2,2n > ${pe_dir}/human_tss_2kb.bed
awk 'BEGIN{OFS="\t"} {s=$2-2000; if(s<0) s=0; print $1,s,$3+2000}' ${mouse_tss} | LC_ALL=C sort -k1,1 -k2,2n > ${pe_dir}/mouse_tss_2kb.bed

echo "Task 4: classify promoter and enhancer peaks"
cd ${pe_dir}
bash classifyingpeaks.sh

echo "Stage Task 5 inputs"
sort_bed() {
  in_file=$1
  out_file=$2
  awk 'BEGIN{FS="[[:space:]]+"; OFS="\t"} NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3 {print}' ${in_file} \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
    > ${out_file}
}

sort_bed ${pe_dir}/humanenhancers.bed ${motif_input_dir}/humanenhancers.sorted.bed
sort_bed ${pe_dir}/humanpromoters.bed ${motif_input_dir}/humanpromoters.sorted.bed
sort_bed ${pe_dir}/mouse_enhancers.bed ${motif_input_dir}/mouse_enhancers.sorted.bed
sort_bed ${pe_dir}/mouse_promoters.bed ${motif_input_dir}/mouse_promoters.sorted.bed
sort_bed ${pe_dir}/human_specific_enhancers_hg38.bed ${motif_input_dir}/human_specific_enhancers_hg38.sorted.bed
sort_bed ${pe_dir}/human_specific_promoters_hg38.bed ${motif_input_dir}/human_specific_promoters_hg38.sorted.bed
sort_bed ${pe_dir}/conserved_enhancers_hg38.bed ${motif_input_dir}/conserved_enhancers_hg38.sorted.bed
sort_bed ${pe_dir}/conserved_promoters_hg38.bed ${motif_input_dir}/conserved_promoters_hg38.sorted.bed
sort_bed ${pe_dir}/mousespecific_enhancers.bed ${motif_input_dir}/mousespecific_enhancers.sorted.bed
sort_bed ${pe_dir}/mousespecific_promoters.bed ${motif_input_dir}/mousespecific_promoters.sorted.bed
sort_bed ${pe_dir}/mouse_conserved_enhancers.bed ${motif_input_dir}/mouse_conserved_enhancers.sorted.bed
sort_bed ${pe_dir}/mouse_conserved_promoters.bed ${motif_input_dir}/mouse_conserved_promoters.sorted.bed

if [[ ${run_homer} -eq 1 ]]; then
  if [[ -n ${homer_dir} ]]; then
    ln -sfn ${homer_dir} ${motif_dir}/homer
  fi

  echo "Task 5: run HOMER motif analyses"
  cd ${motif_input_dir}
  bash build_centered_beds.sh
  bash run_all.sh
  python3 summarize_known_motifs.py
fi

echo "Pipeline finished"
