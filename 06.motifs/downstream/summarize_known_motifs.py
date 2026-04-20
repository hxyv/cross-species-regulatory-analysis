#!/usr/bin/env python3
"""
Per-experiment HOMER known motif summary

Purpose:
  Parse knownResults.txt from each of the 20 HOMER experiment directories
  and produce:
    1. One Markdown overview file with Top 10 tables for every experiment
    2. Per-experiment TSV files for detailed downstream use

This script is a "flat summary" — it treats each experiment independently
and does NOT do cross-experiment comparisons (that's build_integrated_report.py).

Output files (all in results_summary/):
  - motif_summary_overview.md           : combined Markdown with all 20 experiments
  - <experiment>_top_by_qval.tsv        : Top 10 motifs ranked by q-value
  - <experiment>_top_by_delta.tsv       : Top 10 motifs ranked by delta_pct
                                          (only if any pass both filters)

Filters applied:
  - q-value (Benjamini) < 0.05          : FDR-corrected significance
  - delta_pct > 5.0                     : effect size (%Target - %Background)

Usage:
  python3 summarize_known_motifs.py

Relationship to build_integrated_report.py:
  - This script = per-experiment reference / appendix data
  - build_integrated_report.py = cross-experiment integrated analysis
  Both should be run; they serve complementary roles.
"""

import os
import csv
import sys

#
# Setup: resolve paths, create output directory
#

# Work from the script's own directory so that relative paths to out_* dirs
# resolve correctly regardless of where the user invokes the script.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

# All output files go into this subdirectory.
OUT_DIR = os.path.join(SCRIPT_DIR, "results_summary")
os.makedirs(OUT_DIR, exist_ok=True)

#
# Filter thresholds
#

QVAL_CUTOFF = 0.05      # Benjamini-Hochberg corrected p-value threshold
DELTA_PCT_CUTOFF = 5.0   # Minimum effect size: %Target - %Background
TOP_N = 10               # Number of top motifs to show per ranking

#
# Discover experiment directories
#

# Auto-detect all directories named out_* that contain knownResults.txt.
# This makes the script work even if experiments are added or removed later.
EXPERIMENT_DIRS = sorted([
    d for d in os.listdir(".")
    if d.startswith("out_") and os.path.isfile(os.path.join(d, "knownResults.txt"))
])


#
# Parsing function
#

def parse_known_results(filepath):
    """
    Parse a HOMER knownResults.txt file into a list of dictionaries.

    HOMER's knownResults.txt is a tab-separated file with columns:
      [0] Motif Name           - e.g. "Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer"
      [1] Consensus            - e.g. "NNAYTTCCTGHN"
      [2] P-value              - e.g. "1e-78"
      [3] Log P-value          - e.g. "-1.810e+02"
      [4] q-value (Benjamini)  - e.g. "0.0000"
      [5] # Target with Motif  - e.g. "6117.0"
      [6] % Target with Motif  - e.g. "54.77%"
      [7] # Background w/ Motif- e.g. "2481.5"
      [8] % Background w/ Motif- e.g. "45.88%"

    We compute:
      delta_pct = %Target - %Background
    This is the "effect size" — how much more frequent the motif is in
    target vs background. A positive delta means enriched in target.

    Returns: list of dicts, one per motif row.
    """
    rows = []
    with open(filepath, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)  # skip header row
        for line in reader:
            if len(line) < 9:
                continue
            try:
                motif_name = line[0]
                consensus = line[1]
                pval = line[2]
                log_pval = float(line[3])
                qval_str = line[4]
                # Handle missing or non-numeric q-values gracefully
                qval = float(qval_str) if qval_str not in ("", "NA") else 1.0
                # Strip the "%" character before converting to float
                target_pct = float(line[6].replace("%", ""))
                bg_pct = float(line[8].replace("%", ""))
                delta_pct = target_pct - bg_pct

                rows.append({
                    "motif": motif_name,
                    "consensus": consensus,
                    "p_value": pval,
                    "log_p_value": log_pval,
                    "q_value": qval,
                    "target_pct": target_pct,
                    "bg_pct": bg_pct,
                    "delta_pct": round(delta_pct, 2),
                })
            except (ValueError, IndexError):
                # Skip malformed rows rather than crashing
                continue
    return rows


#
# TSV writer
#

def write_tsv(filepath, rows, fieldnames):
    """Write a list of dicts to a tab-separated file with a header row."""
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# Column order for the per-experiment TSV files
FIELDS = ["motif", "consensus", "p_value", "log_p_value", "q_value",
          "target_pct", "bg_pct", "delta_pct"]

#
# Main loop: process each experiment directory
#

# Accumulate Markdown lines for the combined overview file
overview_lines = []
overview_lines.append("# HOMER Known Motif Summary — All 20 Experiments")
overview_lines.append(f"# Filter: q-value < {QVAL_CUTOFF}, delta_pct > {DELTA_PCT_CUTOFF}")
overview_lines.append(f"# Top {TOP_N} shown per experiment per ranking")
overview_lines.append("")

for exp_dir in EXPERIMENT_DIRS:
    txt_path = os.path.join(exp_dir, "knownResults.txt")
    rows = parse_known_results(txt_path)

    # Handle empty or unparseable files
    if not rows:
        overview_lines.append(f"## {exp_dir}")
        overview_lines.append("No motifs parsed.\n")
        continue

    # --- Filter: motifs passing BOTH q-value AND delta_pct thresholds ---
    sig_rows = [r for r in rows
                if r["q_value"] < QVAL_CUTOFF and r["delta_pct"] > DELTA_PCT_CUTOFF]

    # --- Ranking 1: Top N by q-value (regardless of delta filter) ---
    # This shows the statistically most significant motifs, even if their
    # effect size is small. Useful for seeing what HOMER considers most robust.
    all_by_qval = sorted(rows, key=lambda r: r["q_value"])[:TOP_N]

    # --- Ranking 2: Top N by delta_pct (only among those passing both filters) ---
    # This shows the motifs with the largest enrichment difference between
    # target and background. These are the biologically most interesting hits.
    sig_by_delta = sorted(sig_rows, key=lambda r: -r["delta_pct"])[:TOP_N]

    # --- Write per-experiment TSV files ---
    safe_name = exp_dir.replace("/", "_").replace("\\", "_")
    write_tsv(os.path.join(OUT_DIR, f"{safe_name}_top_by_qval.tsv"), all_by_qval, FIELDS)
    if sig_by_delta:
        write_tsv(os.path.join(OUT_DIR, f"{safe_name}_top_by_delta.tsv"), sig_by_delta, FIELDS)

    # --- Append to Markdown overview ---
    overview_lines.append(f"## {exp_dir}")
    overview_lines.append(f"Total motifs tested: {len(rows)}")
    overview_lines.append(f"Significant (q<{QVAL_CUTOFF} & delta>{DELTA_PCT_CUTOFF}): {len(sig_rows)}")
    overview_lines.append("")

    # Table 1: Top by q-value
    overview_lines.append(f"### Top {TOP_N} by q-value")
    overview_lines.append("| Rank | Motif | q-value | %Target | %Bg | delta |")
    overview_lines.append("|------|-------|---------|---------|-----|-------|")
    for i, r in enumerate(all_by_qval, 1):
        overview_lines.append(
            f"| {i} | {r['motif'][:60]} | {r['q_value']:.4f} | {r['target_pct']:.1f}% | {r['bg_pct']:.1f}% | {r['delta_pct']:+.1f} |"
        )
    overview_lines.append("")

    # Table 2: Top by delta (only if any pass filters)
    if sig_by_delta:
        overview_lines.append(f"### Top {TOP_N} by delta_pct (q<{QVAL_CUTOFF})")
        overview_lines.append("| Rank | Motif | q-value | %Target | %Bg | delta |")
        overview_lines.append("|------|-------|---------|---------|-----|-------|")
        for i, r in enumerate(sig_by_delta, 1):
            overview_lines.append(
                f"| {i} | {r['motif'][:60]} | {r['q_value']:.4f} | {r['target_pct']:.1f}% | {r['bg_pct']:.1f}% | {r['delta_pct']:+.1f} |"
            )
        overview_lines.append("")
    else:
        overview_lines.append(f"No motifs passed both q<{QVAL_CUTOFF} and delta>{DELTA_PCT_CUTOFF}.\n")

    overview_lines.append("---\n")

#
# Write the combined Markdown overview
#
overview_path = os.path.join(OUT_DIR, "motif_summary_overview.md")
with open(overview_path, "w", encoding="utf-8") as f:
    f.write("\n".join(overview_lines))

print(f"Done. Summary written to: {OUT_DIR}/")
print(f"  Overview: motif_summary_overview.md")
print(f"  Per-experiment TSVs: {len(EXPERIMENT_DIRS)} experiments x 1-2 files each")

# Reference
# LLM was used for debugging
# [1] OpenAI, "GPT-5.4 Thinking System Card," Mar. 5, 2026. Accessed: Apr. 13, 2026.
