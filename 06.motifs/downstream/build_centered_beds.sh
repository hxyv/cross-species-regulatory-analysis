#!/usr/bin/env bash

# Build summit-centered 200bp BED windows for downstream motif analysis
#
# Algorithm (coordinate transformation, not a statistical method):
#
#   For each ATAC-seq peak in a narrowPeak file, column 10 stores the
#   "summit offset" — the distance (in bp) from the peak start to the
#   position of maximum read pileup (the summit).
#
#   This script re-centers every peak on its summit and trims it to a
#   fixed 200bp window:
#
#     Original peak:  chr1  [start .................... end]   (variable length)
#                                       ↑
#                                    summit = start + col10
#
#     Summit200:      chr1  [summit-100,  summit+100]          (fixed 200bp)
# For example:
# chr1    1000    1500    peak1    500    .    ...    ...    ...    230
# start = 1000 $2
# summit_offset = 230 $10
# summit = start + summit_offset = 1230 ($2 + $10)
# new_start = summit - 100 = 1130
# new_end = summit + 100 = 1330
# chr1    1130    1330    peak1    500    .
#
#   Formula applied to each line:
#     new_start = $2 + $10 - 100      ($2 = start, $10 = summit offset)
#     new_end   = $2 + $10 + 100
#     if new_start < 0, clamp to 0
#
#   - Raw peaks vary from ~200bp to several thousand bp.
#   - Motif analysis on variable-length regions introduces size bias
#     (longer regions contain more background sequence).
#   - The summit is where the actual TF binding or open-chromatin signal
#     is strongest, so motifs are most likely found near the summit.
#   - A uniform 200bp window eliminates length differences and focuses
#     the motif search on the biologically most relevant region.
#
# This script ONLY generates centered BED files. It does NOT run HOMER.
# Run build_centered_beds.sh first, then feed the outputs into
# run_homer_mm10_centered_bidir.sh / run_homer_hg38_centered_bidir.sh.
#
# All 8 input files are narrowPeak-like (10 columns) with a reliable
# summit offset in column 10.
#
# Output files (8 total):
#
#   mm10 group (mouse coordinates):
#     mousespecific_enhancers.summit200.bed
#     mousespecific_promoters.summit200.bed
#     mouse_conserved_enhancers.summit200.bed
#     mouse_conserved_promoters.summit200.bed
#
#   hg38 group (human coordinates):
#     human_specific_enhancers_hg38.summit200.bed
#     human_specific_promoters_hg38.summit200.bed
#     conserved_enhancers_hg38.summit200.bed
#     conserved_promoters_hg38.summit200.bed
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Input file validation
# Fail fast if any required sorted BED is missing or empty.

check_file() {
  local f="$1"
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing or empty file: $f"
    exit 1
  fi
}

# mm10 inputs (native mouse ATAC peaks with reliable summit in col 10)
check_file "mousespecific_enhancers.sorted.bed"
check_file "mousespecific_promoters.sorted.bed"
check_file "mouse_conserved_enhancers.sorted.bed"   # mouse-side conserved peaks (not liftOver)
check_file "mouse_conserved_promoters.sorted.bed"   # mouse-side conserved peaks (not liftOver)

# hg38 inputs (native human ATAC peaks with reliable summit in col 10)
check_file "human_specific_enhancers_hg38.sorted.bed"
check_file "human_specific_promoters_hg38.sorted.bed"
check_file "conserved_enhancers_hg38.sorted.bed"
check_file "conserved_promoters_hg38.sorted.bed"

# Summit-centering function
# summit_center <input.sorted.bed> <output.summit200.bed>
#
# For each line where NF >= 10 (ensures column 10 exists):
#   $2  = peak start coordinate
#   $10 = summit offset (bp from start to summit)
#   summit = $2 + $10           ← absolute genomic position of the summit
#   s = summit - 100            ← new start (200bp window left edge)
#   e = summit + 100            ← new end   (200bp window right edge)
#   if s < 0 then s = 0        ← clamp to avoid negative coordinates
#   Output: chr, s, e, name, score, strand (first 6 BED columns)

summit_center() {
  local in_file="$1"
  local out_file="$2"
  awk 'BEGIN{OFS="\t"} NF>=10{
    summit=$2+$10;             # absolute summit position
    s=summit-100; if(s<0) s=0; # left edge of 200bp window (clamped)
    e=summit+100;              # right edge of 200bp window
    print $1,s,e,$4,$5,$6      # chr, new_start, new_end, name, score, strand
  }' "$in_file" > "$out_file"
}

# A. mm10 group — mouse-specific peaks (native mouse ATAC data)
echo "A. Build summit200 from mouse-specific peaks (mm10) ..."

# Mouse-specific enhancers: mm10 coords, summit from mouse ATAC narrowPeak
summit_center "mousespecific_enhancers.sorted.bed"  "mousespecific_enhancers.summit200.bed"

# Mouse-specific promoters: mm10 coords, summit from mouse ATAC narrowPeak
summit_center "mousespecific_promoters.sorted.bed"  "mousespecific_promoters.summit200.bed"

# B. mm10 group — mouse-side conserved peaks (native mouse ATAC data)
echo "B. Build summit200 from mouse-conserved peaks (mm10) ..."

# Conserved enhancers (mouse side): mm10 coords, real mouse scores/summit
summit_center "mouse_conserved_enhancers.sorted.bed"  "mouse_conserved_enhancers.summit200.bed"

# Conserved promoters (mouse side): mm10 coords, real mouse scores/summit
summit_center "mouse_conserved_promoters.sorted.bed"  "mouse_conserved_promoters.summit200.bed"

# C. hg38 group — human-specific peaks (native human ATAC data)
echo "C. Build summit200 from human-specific peaks (hg38) ..."

# Human-specific enhancers: hg38 coords, summit from human ATAC narrowPeak
summit_center "human_specific_enhancers_hg38.sorted.bed"  "human_specific_enhancers_hg38.summit200.bed"

# Human-specific promoters: hg38 coords, summit from human ATAC narrowPeak
summit_center "human_specific_promoters_hg38.sorted.bed"  "human_specific_promoters_hg38.summit200.bed"

# D. hg38 group — conserved peaks in human coordinates (native human ATAC data)
echo "D. Build summit200 from conserved peaks (hg38) ..."

# Conserved enhancers (human side): hg38 coords, real human scores/summit
summit_center "conserved_enhancers_hg38.sorted.bed"  "conserved_enhancers_hg38.summit200.bed"

# Conserved promoters (human side): hg38 coords, real human scores/summit
summit_center "conserved_promoters_hg38.sorted.bed"  "conserved_promoters_hg38.summit200.bed"

# Validate all 8 outputs exist and are non-empty
echo "Validating generated outputs ..."

for f in \
  mousespecific_enhancers.summit200.bed \
  mousespecific_promoters.summit200.bed \
  mouse_conserved_enhancers.summit200.bed \
  mouse_conserved_promoters.summit200.bed \
  human_specific_enhancers_hg38.summit200.bed \
  human_specific_promoters_hg38.summit200.bed \
  conserved_enhancers_hg38.summit200.bed \
  conserved_promoters_hg38.summit200.bed
do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: output file is empty or missing: $f"
    exit 1
  fi
  echo "$(wc -l < "$f") lines -> $f"   # print line count for quick sanity check
done

echo "Done. All 8 summit-centered BED files generated successfully."

# Reference
# LLM was used for debugging
# [1] OpenAI, "GPT-5.4 Thinking System Card," Mar. 5, 2026. Accessed: Apr. 13, 2026.
