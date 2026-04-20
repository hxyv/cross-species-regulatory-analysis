#!/usr/bin/env bash

# Bidirectional HOMER motif enrichment — hg38 (human) coordinates
#
# Purpose:
#   Run 4 HOMER findMotifsGenome.pl analyses comparing human-specific peaks
#   against conserved (shared) peaks, in BOTH directions, all on the hg38
#   reference genome.
#
# Why bidirectional?
#   Forward  (specific vs conserved): finds motifs enriched in human-specific
#   regions relative to conserved regions.
#   Reversed (conserved vs specific): finds motifs enriched in conserved
#   regions relative to human-specific regions.
#   Together they give a symmetric view of differential motif usage.
#
# 4 comparisons:
#   human_specific_enhancers_hg38 (target) vs conserved_enhancers_hg38 (bg)
#   human_specific_promoters_hg38 (target) vs conserved_promoters_hg38 (bg)
#   conserved_enhancers_hg38 (target) vs human_specific_enhancers_hg38 (bg)
#   conserved_promoters_hg38 (target) vs human_specific_promoters_hg38 (bg)
#
# All input files are in native hg38 coordinates with full 10-column
# narrowPeak format (reliable score, signal, and summit offset).
#
# Prerequisite:
#   hg38 genome must be installed in HOMER before running this script:
#     cd /path/to/homer && perl configureHomer.pl -install hg38
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

# hg38: human reference genome (must be pre-installed in HOMER).
GENOME="hg38"

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

# Human-specific peaks (hg38 native, narrowPeak format)
check_file "human_specific_enhancers_hg38.sorted.bed"
check_file "human_specific_promoters_hg38.sorted.bed"

# Conserved/shared peaks (hg38 native, narrowPeak format)
check_file "conserved_enhancers_hg38.sorted.bed"
check_file "conserved_promoters_hg38.sorted.bed"

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

# Forward direction: what motifs are enriched in human-specific vs conserved?
run_cmd "human_specific_enhancers_hg38.sorted.bed" "out_humanSpecEnh_hg38_vs_conserved"  "conserved_enhancers_hg38.sorted.bed"  "1/4"
run_cmd "human_specific_promoters_hg38.sorted.bed" "out_humanSpecProm_hg38_vs_conserved" "conserved_promoters_hg38.sorted.bed"  "2/4"

# Reversed direction: what motifs are enriched in conserved vs human-specific?
run_cmd "conserved_enhancers_hg38.sorted.bed"  "out_conservedEnh_hg38_vs_humanSpec"  "human_specific_enhancers_hg38.sorted.bed" "3/4"
run_cmd "conserved_promoters_hg38.sorted.bed"  "out_conservedProm_hg38_vs_humanSpec" "human_specific_promoters_hg38.sorted.bed" "4/4"

echo "All four hg38 bidirectional HOMER analyses completed successfully."

# Reference
# LLM was used for debugging
# [1] OpenAI, "GPT-5.4 Thinking System Card," Mar. 5, 2026. Accessed: Apr. 13, 2026.
