#!/usr/bin/env bash

# Bidirectional HOMER motif enrichment — mm10 (mouse) coordinates
#
# Purpose:
#   Run 4 HOMER findMotifsGenome.pl analyses comparing mouse-specific peaks
#   against conserved (shared) peaks, in BOTH directions, all on the mm10
#   reference genome.
#
# Why bidirectional?
#   - Forward  (specific vs conserved): finds motifs enriched in mouse-specific
#     regions relative to conserved regions.
#   - Reversed (conserved vs specific): finds motifs enriched in conserved
#     regions relative to mouse-specific regions.
#   Together they give a symmetric view of differential motif usage.
#
# 4 comparisons:
#   mousespecific_enhancers (target) vs mouse_conserved_enhancers (bg)
#   mousespecific_promoters (target) vs mouse_conserved_promoters (bg)
#   mouse_conserved_enhancers (target) vs mousespecific_enhancers (bg)
#   mouse_conserved_promoters (target) vs mousespecific_promoters (bg)
#
# Important — conserved file choice:
#   We use mouse_conserved_*.sorted.bed (native mouse ATAC peaks that overlap
#   conserved regions). These have:
#     - mm10 coordinates ✓
#     - mouse peak names (mouse_adrenal_gland_...) ✓
#     - real scores and signal values ✓
#     - reliable summit offset in column 10 ✓
#   We do NOT use conservedenhancers.sorted.bed / conservedpromoters.sorted.bed
#   because those were liftOver'd from human peaks (scores = -1, summit
#   offset may be unreliable after cross-species coordinate mapping).
#
# Bash strict mode:
#   -e: exit immediately on any command failure
#   -u: treat unset variables as errors
#   -o pipefail: pipeline returns the exit code of the first failing command

set -euo pipefail

# Environment setup

# Resolve script directory so all relative paths work regardless of
# where the script is invoked from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# HOMER installation lives one level up from downstream/ (i.e., task5/homer).
HOMER_HOME="$SCRIPT_DIR/../homer"
FIND_MOTIFS="$HOMER_HOME/bin/findMotifsGenome.pl"

# Add HOMER's bin/ to PATH so internal helper scripts (bed2pos.pl, homer2,
# findKnownMotifs.pl, etc.) are discoverable by findMotifsGenome.pl.
export PATH="$HOMER_HOME/bin:$PATH"

# Preflight checks


# Verify the main HOMER executable exists and is executable.
if [[ ! -x "$FIND_MOTIFS" ]]; then
  echo "ERROR: findMotifsGenome.pl not found at: $FIND_MOTIFS"
  exit 1
fi

# Verify every internal HOMER helper that findMotifsGenome.pl calls.
# Without these, HOMER will silently fail or produce empty results.
for tool in bed2pos.pl checkPeakFile.pl cleanUpPeakFile.pl mergePeaks homerTools findKnownMotifs.pl homer2 compareMotifs.pl; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: required HOMER helper not found in PATH: $tool"
    exit 1
  fi
done

# Analysis parameters


# mm10: mouse reference genome (already installed in HOMER).
GENOME="mm10"

# Number of parallel threads for HOMER's motif search.
THREADS=8

# Motif lengths to scan: 8-mer, 10-mer, and 12-mer.
LENS="8,10,12"

# Input file validation

# Fail fast if any required input file is missing or empty.
check_file() {
  local f="$1"
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing or empty file: $f"
    exit 1
  fi
}

# Mouse-specific peaks (mm10 native, narrowPeak format, reliable summit)
check_file "mousespecific_enhancers.sorted.bed"
check_file "mousespecific_promoters.sorted.bed"

# Conserved/shared peaks from the mouse side (mm10 native, real scores/summit)
check_file "mouse_conserved_enhancers.sorted.bed"
check_file "mouse_conserved_promoters.sorted.bed"

# HOMER runner function

# run_cmd <target_bed> <output_dir> <background_bed> <step_label>
#
# Runs one findMotifsGenome.pl analysis with:
#   -size given : use actual BED interval sizes (no resizing)
#   -bg         : user-defined background sequences
#   -p          : thread count
#   -len        : motif lengths to test
#
# Post-run checks verify that HOMER produced non-empty known motif results.
run_cmd() {
  local target_bed="$1"
  local out_dir="$2"
  local bg_bed="$3"
  local step_label="$4"

  echo "[$(date '+%F %T')] START [$step_label]: $target_bed vs $bg_bed -> $out_dir"

  # Remove any leftover output from a previous interrupted run.
  rm -rf "$out_dir"

  "$FIND_MOTIFS" "$target_bed" "$GENOME" "$out_dir" \
    -size given \
    -bg "$bg_bed" \
    -p "$THREADS" \
    -len "$LENS"

  # Verify HOMER produced the expected result files.
  if [[ ! -s "$out_dir/knownResults.txt" ]]; then
    echo "ERROR [$step_label]: missing or empty $out_dir/knownResults.txt"
    exit 1
  fi
  if [[ ! -s "$out_dir/knownResults.html" ]]; then
    echo "ERROR [$step_label]: missing or empty $out_dir/knownResults.html"
    exit 1
  fi

  echo "[$(date '+%F %T')] DONE  [$step_label]: $out_dir"
}

# Run 4 comparisons

# Forward direction: what motifs are enriched in mouse-specific vs conserved?
run_cmd "mousespecific_enhancers.sorted.bed" "out_mouseSpecEnh_vs_conserved_bidir"  "mouse_conserved_enhancers.sorted.bed"  "1/4"
run_cmd "mousespecific_promoters.sorted.bed" "out_mouseSpecProm_vs_conserved_bidir" "mouse_conserved_promoters.sorted.bed"  "2/4"

# Reversed direction: what motifs are enriched in conserved vs mouse-specific?
run_cmd "mouse_conserved_enhancers.sorted.bed"  "out_conservedEnh_vs_mouseSpec"  "mousespecific_enhancers.sorted.bed" "3/4"
run_cmd "mouse_conserved_promoters.sorted.bed"  "out_conservedProm_vs_mouseSpec" "mousespecific_promoters.sorted.bed" "4/4"

echo "All four mm10 bidirectional HOMER analyses completed successfully."

# Reference
# LLM was used for debugging
# [1] OpenAI, "GPT-5.4 Thinking System Card," Mar. 5, 2026. Accessed: Apr. 13, 2026.
