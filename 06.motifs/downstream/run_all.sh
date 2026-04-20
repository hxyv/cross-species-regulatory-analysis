#!/usr/bin/env bash

# Master runner: execute all HOMER motif analysis scripts sequentially
#
# This script runs all 6 analysis scripts in the correct order:
#
#   Step 0: build_centered_beds.sh        → generate 8 summit200 BED files
#   Step 1: run_homer_mm10_centered_bidir → 4 centered specific ↔ conserved (mm10)
#   Step 2: run_homer_hg38_centered_bidir → 4 centered specific ↔ conserved (hg38)
#   Step 3: run_homer_mm10_bidirectional  → 4 non-centered specific ↔ conserved (mm10)
#   Step 4: run_homer_hg38_bidirectional  → 4 non-centered specific ↔ conserved (hg38)
#   Step 5: run_homer_allpeaks_mm10       → 2 all peaks vs default bg (mm10)
#   Step 6: run_homer_allpeaks_hg38       → 2 all peaks vs default bg (hg38)
#
# Total: 20 HOMER experiments
#
# Prerequisites:
#   - HOMER installed at ../homer/
#   - mm10 genome installed in HOMER
#   - hg38 genome installed in HOMER (perl configureHomer.pl -install hg38)
#   - All required .sorted.bed files present in this directory
#
# Usage:
#   bash run_all.sh | tee run_all.log
#
# If any step fails, the script stops immediately (set -e).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "  Master runner started at $(date '+%F %T')"
echo "  Working directory: $SCRIPT_DIR"

# Step 0: Generate summit-centered 200bp BED files
# Must run BEFORE the centered HOMER scripts.
# Produces 8 .summit200.bed files from the 8 .sorted.bed inputs.
# echo ""
# echo ">>> Step 0/6: build_centered_beds.sh"
# bash "./build_centered_beds.sh"

# Step 1: mm10 centered bidirectional (4 experiments)
# specific ↔ conserved on summit-centered 200bp windows (mm10)
echo ""
echo ">>> Step 1/6: run_homer_mm10_centered_bidir.sh"
bash "./run_homer_mm10_centered_bidir.sh"

# Step 2: hg38 centered bidirectional (4 experiments)
# specific ↔ conserved on summit-centered 200bp windows (hg38)
echo ""
echo ">>> Step 2/6: run_homer_hg38_centered_bidir.sh"
bash "./run_homer_hg38_centered_bidir.sh"

# Step 3: mm10 non-centered bidirectional (4 experiments)
# specific ↔ conserved on original (variable-length) peak regions (mm10)
echo ""
echo ">>> Step 3/6: run_homer_mm10_bidirectional.sh"
bash "./run_homer_mm10_bidirectional.sh"

# Step 4: hg38 non-centered bidirectional (4 experiments)
# specific ↔ conserved on original (variable-length) peak regions (hg38)
echo ""
echo ">>> Step 4/6: run_homer_hg38_bidirectional.sh"
bash "./run_homer_hg38_bidirectional.sh"

# Step 5: mm10 all peaks vs default background (2 experiments)
# Full mouse enhancer/promoter sets vs HOMER auto-generated background
echo ""
echo ">>> Step 5/6: run_homer_allpeaks_mm10.sh"
bash "./run_homer_allpeaks_mm10.sh"

# Step 6: hg38 all peaks vs default background (2 experiments)
# Full human enhancer/promoter sets vs HOMER auto-generated background
echo ""
echo ">>> Step 6/6: run_homer_allpeaks_hg38.sh"
bash "./run_homer_allpeaks_hg38.sh"

# Done
echo ""
echo "  All 20 HOMER experiments completed successfully!"
echo "  Finished at $(date '+%F %T')"

# Reference
# LLM was used for debugging
# [1] OpenAI, "GPT-5.4 Thinking System Card," Mar. 5, 2026. Accessed: Apr. 13, 2026.
