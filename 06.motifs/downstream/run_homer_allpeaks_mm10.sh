#!/usr/bin/env bash

# HOMER motif analysis on ALL mouse peaks — mm10, default background
#
# Purpose:
#   Run HOMER findMotifsGenome.pl on the complete (unfiltered) set of mouse
#   enhancers and promoters, using HOMER's default automatically generated
#   background (no -bg argument).
#
# Why "all peaks"?
#   The species-specific and conserved analyses focus on subsets of peaks.
#   This script runs motif discovery on the FULL enhancer and promoter peak
#   sets to get a general motif profile of all open chromatin in mouse
#   adrenal tissue — answering the question: "what TF motifs are present
#   across all regulatory regions, regardless of conservation status?"
#
# What "default background" means:
#   When no -bg is provided, HOMER automatically generates background
#   sequences matched by GC content and repeat structure from the reference
#   genome. This is NOT a custom biological control — it is a statistical
#   baseline from the genome itself. The enrichment results tell you which
#   motifs appear more often in your peaks than expected by chance in
#   random genomic regions with similar sequence composition.
#
# What HOMER does internally (findMotifsGenome.pl):
#   1. Reads the target BED → extracts DNA sequences from the mm10 genome
#   2. Auto-generates background sequences from mm10 (GC/repeat-matched)
#   3. Known motif enrichment: scans both sets against a database of ~470
#      known vertebrate motifs, computes enrichment (binomial test)
#   4. De novo motif discovery: searches for novel over-represented sequence
#      patterns in target relative to background using the HOMER2 algorithm
#   5. Outputs: knownResults.txt/html (known), homerMotifs.all.motifs (de novo)
#
# 2 analyses:
#   1) mouse_enhancers.sorted.bed → all mouse enhancer peaks (16,544 regions)
#   2) mouse_promoters.sorted.bed → all mouse promoter peaks (13,034 regions)
#
# Bash strict mode:
#   -e: exit immediately on any command failure
#   -u: treat unset variables as errors
#   -o pipefail: pipeline returns the exit code of the first failing command
set -euo pipefail

# Environment setup

# Resolve this script's directory so all relative paths (input BEDs, HOMER)
# work the same regardless of where the user invokes the script from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# HOMER installation is one level up: task5/homer/
HOMER_HOME="$SCRIPT_DIR/../homer"
FIND_MOTIFS="$HOMER_HOME/bin/findMotifsGenome.pl"

# Add HOMER bin/ to PATH so that findMotifsGenome.pl can locate its
# internal helpers (bed2pos.pl, homer2, findKnownMotifs.pl, etc.).
export PATH="$HOMER_HOME/bin:$PATH"

# Preflight checks

# 1) Verify the main HOMER executable exists and is executable.
if [[ ! -x "$FIND_MOTIFS" ]]; then
  echo "ERROR: findMotifsGenome.pl not found at: $FIND_MOTIFS"
  exit 1
fi

# 2) Verify every internal helper that findMotifsGenome.pl calls at runtime.
#    If any is missing, HOMER will silently fail or produce empty results.
for tool in bed2pos.pl checkPeakFile.pl cleanUpPeakFile.pl mergePeaks homerTools findKnownMotifs.pl homer2 compareMotifs.pl; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: required HOMER helper not found in PATH: $tool"
    exit 1
  fi
done

# Analysis parameters

GENOME="mm10"    # Mouse reference genome (already installed in HOMER)
THREADS=8        # Parallel threads for motif scanning
LENS="8,10,12"   # Motif widths to search: 8bp, 10bp, and 12bp

# Input file validation

# check_file: abort immediately if a required input is missing or empty (0 bytes).
# Uses -s (exists AND size > 0) rather than -f (exists only).
check_file() {
  local f="$1"
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing or empty file: $f"
    exit 1
  fi
}

# Full (unfiltered) mouse enhancer and promoter peak sets in mm10 coordinates.
# These are the complete ATAC-seq peak sets — not subsetted by species-specificity
# or conservation status.
check_file "mouse_enhancers.sorted.bed"    # 16,544 peaks
check_file "mouse_promoters.sorted.bed"    # 13,034 peaks

# HOMER runner function

# run_cmd <target_bed> <output_dir> <step_label>
#
# Workflow inside this function:
#   1. Log the start time and target being analyzed
#   2. rm -rf the output directory to prevent stale results from a prior run
#   3. Call findMotifsGenome.pl WITHOUT -bg:
#        - Reads target_bed → extracts sequences from $GENOME
#        - HOMER auto-generates background sequences from the genome
#          (matched by GC content and repeat structure)
#        - First run may preparse the genome (takes a few minutes);
#          subsequent runs reuse the cached preparsed data
#        - Runs known motif enrichment (binomial test, ~470 vertebrate motifs)
#        - Runs de novo motif discovery (HOMER2 algorithm, 25 iterations/length)
#        - Writes results to output_dir/
#   4. Post-run validation: check that knownResults.txt and .html exist and
#      are non-empty; abort if HOMER silently failed
#   5. Log completion time
run_cmd() {
  local target_bed="$1"   # BED file of all peaks to analyze
  local out_dir="$2"      # Directory where HOMER writes all output files
  local step_label="$3"   # Human-readable label for logging (e.g. "1/2")

  echo "[$(date '+%F %T')] START [$step_label]: $target_bed -> $out_dir (default bg)"

  # Clean slate: remove any leftover output from a previous interrupted run
  rm -rf "$out_dir"

  # Core HOMER command (NO -bg → default background):
  #   -size given   : use each BED interval's actual coordinates
  #   -p <threads>  : number of CPUs for parallel motif scanning
  #   -len <widths> : motif lengths to search (comma-separated)
  "$FIND_MOTIFS" "$target_bed" "$GENOME" "$out_dir" \
    -size given \
    -p "$THREADS" \
    -len "$LENS"

  # Post-run validation: HOMER can exit 0 but produce no results if
  # internal helpers are missing or genome data is incomplete.
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

# Run 2 analyses

# All mouse enhancers: general motif profile of mouse adrenal enhancers
# → Answers: what TF motifs are enriched across ALL enhancers vs random genome?
run_cmd "mouse_enhancers.sorted.bed"  "out_allMouseEnh_mm10_defaultBg"  "1/2"

# All mouse promoters: general motif profile of mouse adrenal promoters
# → Answers: what TF motifs are enriched across ALL promoters vs random genome?
run_cmd "mouse_promoters.sorted.bed"  "out_allMouseProm_mm10_defaultBg" "2/2"

echo "All mouse all-peaks HOMER analyses completed successfully."

# Reference
# LLM was used for debugging
# [1] OpenAI, "GPT-5.4 Thinking System Card," Mar. 5, 2026. Accessed: Apr. 13, 2026.
