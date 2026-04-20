#!/usr/bin/env bash

# BED/narrowPeak standardization script
# What it does for each input file:
# Keep only valid BED-like lines (>=3 columns, numeric start/end, start < end)
# Accept any whitespace as input separator, normalize output to tab
# Sort by chromosome + start + end (explicit tab delimiter)
# Deduplicate by first 3 columns; for duplicated coordinates keep higher score
# (score = column 5 when numeric). If score is missing/non-numeric, keep first.
# Write all outputs into task4/task5/

set -euo pipefail

# Locate this script's directory and run from there (task4)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Output directory: all sorted results are saved here
OUT_DIR="$SCRIPT_DIR/task5"
mkdir -p "$OUT_DIR"

# Input file list to process
FILES=(
  "conserved_enhancers_hg38.bed"
  "conserved_promoters_hg38.bed"
  "conservedenhancers.bed"
  "conservedpromoters.bed"
  "mouse_conserved_enhancers.bed"
  "mouse_conserved_promoters.bed"
  "mousespecific_enhancers.bed"
  "mousespecific_promoters.bed"
  "human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed"
  "human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed"
  "human_adrenal_idr_optimal.narrowPeak"
  "human_specific_enhancers_hg38.bed"
  "human_specific_promoters_hg38.bed"
  "human_tss_2kb.bed"
  "humanenhancers.bed"
  "humanpromoters.bed"
  "humanspecific_enhancers.bed"
  "humanspecific_promoters.bed"
  "mouse_adrenal_idr_optimal.narrowPeak"
  "mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed"
  "mouse_adrenal_idr_optimal.shared_with_human_mapped.bed"
  "mouse_enhancers.bed"
  "mouse_promoters.bed"
  "mouse_tss_2kb.bed"
)

echo "Working directory: $SCRIPT_DIR"
echo "Output directory:  $OUT_DIR"
echo "Start sorting and deduplicating BED-like files..."

# Process each file one by one
for in_file in "${FILES[@]}"; do
  # Skip file if not found
  if [[ ! -f "$in_file" ]]; then
    echo "Skip missing file: $in_file"
    continue
  fi

  # Build output filename:
  # remove original extension (.bed or .narrowPeak)
  # append .sorted.bed
  base_name="${in_file%.*}"
  out_file="$OUT_DIR/${base_name}.sorted.bed"

  # Main pipeline:
  # Parse robustly with any whitespace, filter invalid rows, normalize to tab
  # Sort by chr/start/end using tab as delimiter
  # For duplicated coordinates, keep row with higher numeric column-5 score
  # Fallback to first row when score is not numeric
  awk 'BEGIN{FS="[[:space:]]+"; OFS="\t"}
       NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3 {print}' "$in_file" \
    | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n \
    | awk '
        BEGIN{
          FS=OFS="\t"
          cur_key=""
          best_line=""
          best_has_score=0
          best_score=0
        }
        {
          key=$1 OFS $2 OFS $3
          score_is_num=($5 ~ /^-?[0-9]+([.][0-9]+)?$/)
          score_val=score_is_num ? ($5+0) : 0

          if (cur_key=="") {
            cur_key=key
            best_line=$0
            best_has_score=score_is_num
            best_score=score_val
            next
          }

          if (key != cur_key) {
            print best_line
            cur_key=key
            best_line=$0
            best_has_score=score_is_num
            best_score=score_val
            next
          }

          if (score_is_num && (!best_has_score || score_val > best_score)) {
            best_line=$0
            best_has_score=1
            best_score=score_val
          }
        }
        END{
          if (best_line != "") print best_line
        }' \
    > "$out_file"

  # Simple row counts for logging
  in_count=$(awk 'BEGIN{c=0} NF>0 {c++} END{print c}' "$in_file")
  out_count=$(awk 'BEGIN{c=0} NF>0 {c++} END{print c}' "$out_file")
  echo "$in_file -> $out_file (in=$in_count out=$out_count)"
done

echo "Done."

# Reference

# LLM was used for debugging

# [1] OpenAI, “GPT-5.4 Thinking System Card,” Mar. 5, 2026. Accessed: Apr. 13, 2026.