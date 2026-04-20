#!/usr/bin/env python3
"""
Integrated HOMER Motif Analysis Report (v2)

Purpose:
  Cross-experiment analysis of all 20 HOMER motif enrichment results.
  This script reads knownResults.txt and de novo motif files from every
  experiment directory and produces a single integrated Markdown report
  with 10 analytical sections.

Sections produced:
  1)  Global motif landscape       — all peaks vs default bg (baseline)
  2)  Reciprocal validation        — forward vs reverse consistency
  3)  Centered vs non-centered     — summit200 vs original region consistency
  4)  De novo motif highlights     — top de novo per experiment
  5)  High-confidence motif list   — 4-criteria filtered final list
  6)  Cross-species comparison     — mouse-specific vs human-specific motifs
  7)  Enhancer vs promoter         — within-species regulatory element diff
  8)  TF family grouping           — aggregate by TF family (ETS, bZIP, etc.)
  9)  Cross-experiment matrix      — presence/absence TSV for heatmap
  10) De novo vs known validation  — mutual support check

Output files (all in results_summary/):
  - integrated_motif_report.md     : the full 10-section Markdown report
  - motif_presence_matrix.tsv      : N motifs × 20 experiments binary matrix

Relationship to summarize_known_motifs.py:
  - summarize_known_motifs.py = per-experiment flat summary + TSVs
  - This script = cross-experiment integrated analysis
  Both are complementary; run both.

Usage:
  python3 build_integrated_report.py
"""

import os
import csv
import re
from collections import defaultdict

#
# Setup: resolve paths, create output directory
#

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

# Significance and effect-size thresholds used throughout
QVAL_CUTOFF = 0.05    # Benjamini-Hochberg corrected p-value
DELTA_CUTOFF = 5.0     # Minimum %Target - %Background
TOP_N = 15             # Max rows per table in the report

OUT_DIR = os.path.join(SCRIPT_DIR, "results_summary")
os.makedirs(OUT_DIR, exist_ok=True)
HIGHCONF_AUDIT_PATH = os.path.join(OUT_DIR, "high_confidence_criteria_audit.tsv")

#
# Experiment directory configuration
#

# Bidirectional pairs for reciprocal validation (Section 2)
# Each tuple: (forward_dir, reverse_dir, human-readable label)
# "Forward" = species-specific as target, conserved as bg
# "Reverse" = conserved as target, species-specific as bg
# If a motif is enriched in BOTH directions, it's likely an artifact.
RECIPROCAL_PAIRS = [
    # mm10 non-centered
    ("out_mouseSpecEnh_vs_conserved_bidir", "out_conservedEnh_vs_mouseSpec", "Mouse Enhancers (mm10, original)"),
    ("out_mouseSpecProm_vs_conserved_bidir", "out_conservedProm_vs_mouseSpec", "Mouse Promoters (mm10, original)"),
    # mm10 centered
    ("out_mouseSpecEnh_centered_bidir", "out_conservedEnh_centered_vs_mouseSpec", "Mouse Enhancers (mm10, summit200)"),
    ("out_mouseSpecProm_centered_bidir", "out_conservedProm_centered_vs_mouseSpec", "Mouse Promoters (mm10, summit200)"),
    # hg38 non-centered
    ("out_humanSpecEnh_hg38_vs_conserved", "out_conservedEnh_hg38_vs_humanSpec", "Human Enhancers (hg38, original)"),
    ("out_humanSpecProm_hg38_vs_conserved", "out_conservedProm_hg38_vs_humanSpec", "Human Promoters (hg38, original)"),
    # hg38 centered
    ("out_humanSpecEnh_hg38_centered_bidir", "out_conservedEnh_hg38_centered_vs_humanSpec", "Human Enhancers (hg38, summit200)"),
    ("out_humanSpecProm_hg38_centered_bidir", "out_conservedProm_hg38_centered_vs_humanSpec", "Human Promoters (hg38, summit200)"),
]

# All-peaks experiments (Section 1)
# These use HOMER default background (random genome), giving a global baseline
# of what motifs are enriched in ALL enhancers/promoters.
ALL_PEAKS_DIRS = [
    ("out_allMouseEnh_mm10_defaultBg", "All Mouse Enhancers"),
    ("out_allMouseProm_mm10_defaultBg", "All Mouse Promoters"),
    ("out_allHumanEnh_hg38_defaultBg", "All Human Enhancers"),
    ("out_allHumanProm_hg38_defaultBg", "All Human Promoters"),
]

# Centered vs non-centered pairs (Section 3)
# Each pair: (non-centered_dir, centered_dir, label)
# Motifs appearing in BOTH are considered more robust (size-independent).
CONSISTENCY_PAIRS = [
    ("out_mouseSpecEnh_vs_conserved_bidir", "out_mouseSpecEnh_centered_bidir", "Mouse Spec Enh vs Conserved"),
    ("out_mouseSpecProm_vs_conserved_bidir", "out_mouseSpecProm_centered_bidir", "Mouse Spec Prom vs Conserved"),
    ("out_humanSpecEnh_hg38_vs_conserved", "out_humanSpecEnh_hg38_centered_bidir", "Human Spec Enh vs Conserved"),
    ("out_humanSpecProm_hg38_vs_conserved", "out_humanSpecProm_hg38_centered_bidir", "Human Spec Prom vs Conserved"),
]

# Cross-species comparison dirs (Section 6)
# Using centered (summit200), forward direction only.
# Compares which motifs are enriched in mouse-specific vs human-specific.
CROSS_SPECIES = {
    "mouse_enh": "out_mouseSpecEnh_centered_bidir",
    "human_enh": "out_humanSpecEnh_hg38_centered_bidir",
    "mouse_prom": "out_mouseSpecProm_centered_bidir",
    "human_prom": "out_humanSpecProm_hg38_centered_bidir",
}

# Enhancer vs Promoter comparison dirs (Section 7)
# Within the same species, compares enhancer-enriched vs promoter-enriched motifs.
ENH_VS_PROM = [
    ("out_mouseSpecEnh_centered_bidir", "out_mouseSpecProm_centered_bidir", "Mouse Specific (centered)"),
    ("out_humanSpecEnh_hg38_centered_bidir", "out_humanSpecProm_hg38_centered_bidir", "Human Specific (centered)"),
]


#
# Parsing functions
#

def parse_known_results(dirpath):
    """
    Parse HOMER knownResults.txt → dict keyed by motif name.

    HOMER's knownResults.txt is the primary output for known motif enrichment.
    It is a tab-separated file where each row represents one known motif from
    HOMER's internal database (~470 vertebrate motifs). The columns are:

      [0] Motif Name            e.g. "Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer"
      [1] Consensus             e.g. "NNAYTTCCTGHN" (IUPAC degenerate sequence)
      [2] P-value               e.g. "1e-78" (raw binomial test p-value)
      [3] Log P-value           e.g. "-1.810e+02" (natural log of p-value)
      [4] q-value (Benjamini)   e.g. "0.0000" (FDR-corrected p-value)
      [5] # Target with Motif   e.g. "6117.0" (count of target seqs containing motif)
      [6] % Target with Motif   e.g. "54.77%" (percentage of target seqs)
      [7] # Background w/ Motif e.g. "2481.5" (count of bg seqs containing motif)
      [8] % Background w/ Motif e.g. "45.88%" (percentage of bg seqs)

    We compute:
      delta_pct = %Target - %Background
    A positive delta means the motif is MORE frequent in target than background
    (i.e., enriched). A negative delta would mean depleted (rare in this report
    since HOMER sorts by enrichment).

    Returns:
      dict keyed by motif name → {q_value, target_pct, bg_pct, delta_pct, log_p, consensus}
      Returns empty dict if the file doesn't exist or can't be parsed.
    """
    filepath = os.path.join(dirpath, "knownResults.txt")
    if not os.path.isfile(filepath):
        return {}
    results = {}
    with open(filepath, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)  # skip the header row
        for line in reader:
            if len(line) < 9:  # need at least 9 columns for a valid row
                continue
            try:
                motif = line[0]              # full HOMER motif name
                # q-value: FDR-corrected significance; treat missing/NA as non-significant
                qval = float(line[4]) if line[4] not in ("", "NA") else 1.0
                # Strip "%" suffix before converting to float
                tgt_pct = float(line[6].replace("%", ""))  # % of target seqs with motif
                bg_pct = float(line[8].replace("%", ""))   # % of background seqs with motif
                results[motif] = {
                    "q_value": qval,
                    "target_pct": tgt_pct,
                    "bg_pct": bg_pct,
                    "delta_pct": round(tgt_pct - bg_pct, 2),  # effect size
                    "log_p": float(line[3]),    # log p-value for secondary ranking
                    "consensus": line[1],       # IUPAC consensus sequence
                }
            except (ValueError, IndexError):
                # Skip malformed rows (e.g. truncated lines) rather than crashing
                continue
    return results


def parse_denovo_motifs(dirpath):
    """
    Parse homerMotifs.all.motifs → list of dicts with rank, consensus, p-value, id.

    HOMER's de novo output file contains one entry per discovered motif.
    Each motif starts with a ">" header line (like FASTA format) followed by
    the position weight matrix (PWM) rows. We only parse the header lines.

    Header line format (tab-separated fields):
      [0] Consensus sequence   e.g. "GNACTTCCTG"
      [1] Motif ID             e.g. "1-GNACTTCCTG" (rank-consensus)
      [2] Log odds score       e.g. "7.678998"
      [3] Log p-value          e.g. "-225.309589"
      [4] Unknown/zero         e.g. "0"
      [5] Statistics string    e.g. "T:3903.0(34.95%),B:1392.5(25.94%),P:1e-97"
           T = target count/pct, B = background count/pct, P = p-value

    Returns:
      List of dicts: [{rank, consensus, p_value, id}, ...]
      Ordered by HOMER's internal ranking (best first).
    """
    filepath = os.path.join(dirpath, "homerMotifs.all.motifs")
    if not os.path.isfile(filepath):
        return []
    motifs = []
    rank = 0
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):  # header line = new motif entry
                rank += 1
                parts = line.strip().lstrip(">").split("\t")
                consensus = parts[0] if parts else "?"
                # Extract p-value from the statistics string in field [5]
                # Format: "T:3903.0(34.95%),B:1392.5(25.94%),P:1e-97"
                info = parts[5] if len(parts) > 5 else ""
                pval_match = re.search(r"P:([\d.e\-+]+)", info)
                pval_str = pval_match.group(1) if pval_match else "NA"
                # Field [1] is HOMER's internal motif ID (rank-consensus)
                best_match = parts[1] if len(parts) > 1 else "?"
                motifs.append({
                    "rank": rank,          # 1-based discovery order
                    "consensus": consensus, # IUPAC degenerate sequence
                    "p_value": pval_str,   # enrichment p-value as string
                    "id": best_match,      # HOMER internal motif identifier
                })
    return motifs


def parse_denovo_html_bestmatch(dirpath):
    """
    Parse homerResults.html for de novo best-match known motif names.

    After de novo discovery, HOMER compares each novel motif against its
    known motif database and reports the best match in the HTML output.
    Lines look like:
      "similar to Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer(0.953)"
    where 0.953 is the similarity score.

    We extract just the motif name portion (before the parenthesized score)
    for cross-validation: if a de novo motif's best match is ALSO significant
    in the known motif enrichment results, that's mutual confirmation.

    Returns:
      List of known motif name strings that de novo motifs matched to.
    """
    filepath = os.path.join(dirpath, "homerResults.html")
    if not os.path.isfile(filepath):
        return []
    matches = []
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    # Regex captures text after "similar to " up to the next "<" or "(" character
    for m in re.finditer(r"similar to ([^<(]+)", content):
        name = m.group(1).strip()
        if name:
            matches.append(name)
    return matches


def get_sig_motifs(results, qval_cut=QVAL_CUTOFF, delta_cut=DELTA_CUTOFF):
    """
    Return set of motif names passing both significance filters.

    A motif is considered "significant" if:
      1. q-value < qval_cut (statistically significant after FDR correction)
      2. delta_pct > delta_cut (biologically meaningful effect size)

    Both conditions must be met simultaneously. This avoids:
      - Statistically significant but tiny effects (q passes, delta doesn't)
      - Large effects that are not statistically robust (delta passes, q doesn't)

    Args:
      results: dict from parse_known_results()
      qval_cut: FDR threshold (default 0.05)
      delta_cut: minimum %Target - %Background (default 5.0)

    Returns:
      Set of motif name strings.
    """
    return {m for m, v in results.items()
            if v["q_value"] < qval_cut and v["delta_pct"] > delta_cut}


def top_motifs_table(results, n=TOP_N, sort_key="delta_pct"):
    """
    Return top N significant motifs sorted by the given key.

    First filters to only significant motifs (q < cutoff AND delta > cutoff),
    then sorts by the specified key:
      - "delta_pct": largest effect size first (default, most biologically relevant)
      - "q_value": most statistically significant first

    Args:
      results: dict from parse_known_results()
      n: max number of motifs to return
      sort_key: which field to sort by

    Returns:
      List of (motif_name, value_dict) tuples, length <= n.
    """
    items = [(m, v) for m, v in results.items()
             if v["q_value"] < QVAL_CUTOFF and v["delta_pct"] > DELTA_CUTOFF]
    items.sort(key=lambda x: -x[1][sort_key] if sort_key == "delta_pct" else x[1]["q_value"])
    return items[:n]


def extract_tf_family(motif_name):
    """
    Extract TF family from HOMER motif name.

    HOMER motif names follow a consistent pattern:
      "MotifName(Family)/ChIP-Seq-Dataset/Source"

    Examples:
      "Etv2(ETS)/ES-ER71-ChIP-Seq(GSE59402)/Homer"      → "ETS"
      "Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer"   → "bZIP"
      "CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer" → "Zf"

    The family name is in the FIRST set of parentheses after the motif name.
    We use regex to extract the content between the first "(" and ")".

    Common families in HOMER's vertebrate database:
      ETS, bZIP, Zf (zinc finger), Homeobox, NR (nuclear receptor),
      bHLH, Forkhead, CCAAT, HTH (helix-turn-helix), etc.

    Returns:
      Family string, or "Unknown" if no parenthesized family is found.
    """
    m = re.search(r"\(([^)]+)\)", motif_name)
    return m.group(1) if m else "Unknown"


# Auto-detect all experiment directories (out_*) that contain knownResults.txt.
# Sorted alphabetically for consistent ordering across runs.
# This makes the script work even if experiments are added or removed later.
all_dirs = sorted([
    d for d in os.listdir(".")
    if d.startswith("out_") and os.path.isfile(os.path.join(d, "knownResults.txt"))
])

#
# Build report — accumulate Markdown lines
#

lines = []

lines.append("# Integrated HOMER Motif Analysis Report (v2)")
lines.append("")
lines.append(f"Filters: q-value < {QVAL_CUTOFF}, delta_pct (=%Target - %Background) > {DELTA_CUTOFF}")
lines.append(f"")
lines.append(f"Total experiments analyzed: {len(all_dirs)}")
lines.append("")

# Section 1: Global Landscape
# Establishes a baseline: what motifs are common across ALL open chromatin
# peaks. Motifs found here are "globally common" and will be excluded from
# the high-confidence list in Section 5 to avoid false positives.

lines.append("---")
lines.append("## 1. Global Motif Landscape (All Peaks vs Default Background)")
lines.append("")
lines.append("Which TF motifs are generally enriched across ALL enhancers/promoters")
lines.append("relative to random genomic regions. These are 'globally common' motifs.")
lines.append("")

global_common = set()
for d, label in ALL_PEAKS_DIRS:
    res = parse_known_results(d)
    sig = get_sig_motifs(res)
    global_common.update(sig)
    top = top_motifs_table(res, n=10)
    lines.append(f"### {label}")
    lines.append(f"Significant motifs: **{len(sig)}**")
    lines.append("")
    if top:
        lines.append("| Rank | Motif | %Target | %Bg | delta | q-value |")
        lines.append("|------|-------|---------|-----|-------|---------|")
        for i, (m, v) in enumerate(top, 1):
            lines.append(f"| {i} | {m[:55]} | {v['target_pct']:.1f}% | {v['bg_pct']:.1f}% | {v['delta_pct']:+.1f} | {v['q_value']:.4f} |")
        lines.append("")
    else:
        lines.append("No motifs passed filters.\n")

lines.append(f"**Total globally common motifs (union across 4 all-peaks runs): {len(global_common)}**")
lines.append("")

# Section 2: Reciprocal Validation
# Core quality control: for each bidirectional pair, a motif enriched in the
# forward direction (specific > conserved) should NOT be enriched in the
# reverse direction (conserved > specific). Motifs in BOTH are artifacts.
# We also exclude globally common motifs to focus on truly differential ones.

lines.append("---")
lines.append("## 2. Reciprocal Validation (Specific vs Conserved, Bidirectional)")
lines.append("")
lines.append("Motifs enriched in forward (specific > conserved) should NOT be enriched")
lines.append("in reverse (conserved > specific). Motifs in BOTH directions are unreliable.")
lines.append("")

for fwd_dir, rev_dir, label in RECIPROCAL_PAIRS:
    fwd_res = parse_known_results(fwd_dir)
    rev_res = parse_known_results(rev_dir)
    fwd_sig = get_sig_motifs(fwd_res)
    rev_sig = get_sig_motifs(rev_res)
    both = fwd_sig & rev_sig           # enriched in BOTH = unreliable
    fwd_only = fwd_sig - rev_sig       # enriched only in forward = specific-enriched
    rev_only = rev_sig - fwd_sig       # enriched only in reverse = conserved-enriched
    fwd_not_global = fwd_only - global_common  # further exclude globally common
    rev_not_global = rev_only - global_common

    lines.append(f"### {label}")
    lines.append(f"- Forward: **{len(fwd_sig)}** sig | Reverse: **{len(rev_sig)}** sig | Both (unreliable): **{len(both)}**")
    lines.append(f"- **Specific-only** (excl. global): **{len(fwd_not_global)}** | **Conserved-only** (excl. global): **{len(rev_not_global)}**")
    lines.append("")

    if fwd_not_global:
        top_fwd = sorted(fwd_not_global, key=lambda m: -fwd_res[m]["delta_pct"])[:10]
        lines.append("**Top species-specific motifs:**")
        lines.append("")
        lines.append("| Motif | %Target | %Bg | delta | q-value |")
        lines.append("|-------|---------|-----|-------|---------|")
        for m in top_fwd:
            v = fwd_res[m]
            lines.append(f"| {m[:55]} | {v['target_pct']:.1f}% | {v['bg_pct']:.1f}% | {v['delta_pct']:+.1f} | {v['q_value']:.4f} |")
        lines.append("")

    if rev_not_global:
        top_rev = sorted(rev_not_global, key=lambda m: -rev_res[m]["delta_pct"])[:10]
        lines.append("**Top conserved motifs:**")
        lines.append("")
        lines.append("| Motif | %Target | %Bg | delta | q-value |")
        lines.append("|-------|---------|-----|-------|---------|")
        for m in top_rev:
            v = rev_res[m]
            lines.append(f"| {m[:55]} | {v['target_pct']:.1f}% | {v['bg_pct']:.1f}% | {v['delta_pct']:+.1f} | {v['q_value']:.4f} |")
        lines.append("")
    lines.append("---\n")

# Section 3: Consistency
# Motifs appearing in BOTH centered (summit200) and non-centered (original)
# analyses are more robust — their enrichment is independent of peak size.
# Jaccard index measures the overlap fraction.

lines.append("## 3. Centered vs Non-centered Consistency")
lines.append("")

for nc_dir, c_dir, label in CONSISTENCY_PAIRS:
    nc_res = parse_known_results(nc_dir)
    c_res = parse_known_results(c_dir)
    nc_sig = get_sig_motifs(nc_res)
    c_sig = get_sig_motifs(c_res)
    overlap = nc_sig & c_sig
    union = nc_sig | c_sig
    pct = len(overlap) / len(union) * 100 if union else 0

    lines.append(f"### {label}")
    lines.append(f"- Non-centered: **{len(nc_sig)}** | Centered: **{len(c_sig)}** | Overlap: **{len(overlap)}** | Jaccard: **{pct:.0f}%**")
    lines.append("")

    if overlap:
        top_overlap = sorted(overlap, key=lambda m: -c_res[m]["delta_pct"])[:10]
        lines.append("| Motif | delta (centered) | delta (non-centered) |")
        lines.append("|-------|-----------------|---------------------|")
        for m in top_overlap:
            lines.append(f"| {m[:55]} | {c_res[m]['delta_pct']:+.1f} | {nc_res[m]['delta_pct']:+.1f} |")
        lines.append("")
    lines.append("---\n")

# Section 4: De Novo
# HOMER discovers novel motifs (de novo) in addition to scanning known ones.
# These are ranked by HOMER internally. We show top 5 per experiment.

lines.append("## 4. De Novo Motif Highlights")
lines.append("")
lines.append("Top 5 de novo motifs per experiment.")
lines.append("")

for d in all_dirs:
    denovo = parse_denovo_motifs(d)
    if not denovo:
        continue
    lines.append(f"### {d}")
    lines.append("")
    lines.append("| Rank | Consensus | P-value | ID |")
    lines.append("|------|-----------|---------|----|")
    for m in denovo[:5]:
        lines.append(f"| {m['rank']} | `{m['consensus']}` | {m['p_value']} | {m['id'][:40]} |")
    lines.append("")

# Section 5: High Confidence Motif Summary
# The strictest filter: motifs must pass ALL four criteria simultaneously.
# This is the "final answer" list for each comparison.

lines.append("---")
lines.append("## 5. High-Confidence Motif Summary")
lines.append("")
lines.append("Motifs passing: (1) significant, (2) reciprocal validated,")
lines.append("(3) centered+non-centered consistent, (4) not globally common.")
lines.append("")

# Per-motif audit rows for Section 5.
# One row = one motif under one comparison label, with binary flags for each
# criterion. This is designed for manual QC: users can sort/filter this TSV
# to see exactly where motifs fail (e.g., reverse overlap vs global-common).
highconf_audit_rows = []
# Column definitions for the audit TSV.
# All criterion columns are encoded as 1/0 so the file is easy to filter in
# Excel/R/Python without additional parsing.
highconf_audit_fields = [
    "comparison",
    "motif",
    "criterion1_significant_any",
    "criterion1_sig_noncentered",
    "criterion1_sig_centered",
    "criterion2_not_in_reverse",
    "criterion3_consistent_nc_and_centered",
    "criterion4_not_global_common",
    "pass_all_4",
    "q_noncentered",
    "delta_noncentered",
    "q_centered",
    "delta_centered",
]

for nc_dir, c_dir, label in CONSISTENCY_PAIRS:
    # Forward (specific vs conserved) results in non-centered and centered runs.
    fwd_nc = parse_known_results(nc_dir)
    fwd_c = parse_known_results(c_dir)

    # Significant motif sets under the global thresholds (q and delta cutoffs).
    nc_sig = get_sig_motifs(fwd_nc)
    c_sig = get_sig_motifs(fwd_c)

    # Criterion 3: must be significant in both non-centered and centered runs.
    # This guards against motifs that only appear under one windowing strategy.
    consistent = nc_sig & c_sig

    # Find the reverse-direction paired runs for reciprocal validation.
    # Example: if forward is specific>conserved, reverse is conserved>specific.
    rev_nc_dir = [r for f, r, _ in RECIPROCAL_PAIRS if f == nc_dir]
    rev_c_dir = [r for f, r, _ in RECIPROCAL_PAIRS if f == c_dir]
    rev_sig = set()
    for rd in rev_nc_dir + rev_c_dir:
        rev_sig |= get_sig_motifs(parse_known_results(rd))

    # Final Section 5 set:
    #   consistent (criterion 3) minus reverse-significant motifs (criterion 2)
    #   minus globally common motifs from all-peaks (criterion 4).
    # Criterion 1 is implicit because "consistent" is built from significant sets.
    high_conf = consistent - rev_sig - global_common

    # Build an explicit per-motif audit table so manual review can see
    # exactly which criterion each motif passed/failed.
    # Candidate universe for auditing = all motifs observed in either forward
    # result table (centered or non-centered), not only significant ones.
    candidate_motifs = sorted(set(fwd_nc.keys()) | set(fwd_c.keys()))
    for motif in candidate_motifs:
        # Criterion 1 (significance): store per-run and merged flags so users can
        # distinguish "only one run significant" vs "both runs significant".
        c1_nc = motif in nc_sig
        c1_c = motif in c_sig
        c1_any = c1_nc or c1_c

        # Criterion 2: motif should NOT be significant in reverse direction.
        c2 = motif not in rev_sig

        # Criterion 3: motif should be significant in both centered/non-centered.
        c3 = motif in consistent

        # Criterion 4: motif should NOT be globally common in all-peaks baseline.
        c4 = motif not in global_common

        # Final pass flag for this motif under this comparison.
        pass_all = c1_any and c2 and c3 and c4

        # Keep raw q/delta values to support manual review with exact numbers.
        vnc = fwd_nc.get(motif, {})
        vc = fwd_c.get(motif, {})
        highconf_audit_rows.append({
            "comparison": label,
            "motif": motif,
            "criterion1_significant_any": int(c1_any),
            "criterion1_sig_noncentered": int(c1_nc),
            "criterion1_sig_centered": int(c1_c),
            "criterion2_not_in_reverse": int(c2),
            "criterion3_consistent_nc_and_centered": int(c3),
            "criterion4_not_global_common": int(c4),
            "pass_all_4": int(pass_all),
            "q_noncentered": vnc.get("q_value", "NA"),
            "delta_noncentered": vnc.get("delta_pct", "NA"),
            "q_centered": vc.get("q_value", "NA"),
            "delta_centered": vc.get("delta_pct", "NA"),
        })

    lines.append(f"### {label}")
    lines.append(f"High-confidence motifs: **{len(high_conf)}**")
    lines.append("")
    if high_conf:
        top_hc = sorted(high_conf, key=lambda m: -fwd_c[m]["delta_pct"])[:15]
        lines.append("| Motif | Consensus | delta (centered) | delta (non-centered) | q-value |")
        lines.append("|-------|-----------|-----------------|---------------------|---------|")
        for m in top_hc:
            vc = fwd_c[m]
            vnc = fwd_nc.get(m, {"delta_pct": 0})
            lines.append(f"| {m[:50]} | `{vc['consensus']}` | {vc['delta_pct']:+.1f} | {vnc['delta_pct']:+.1f} | {vc['q_value']:.4f} |")
        lines.append("")
    else:
        lines.append("No motifs passed all four criteria.\n")

    # Compact per-comparison failure counts shown directly in Markdown.
    # These are quick diagnostics to explain why high_conf can be zero.
    failed_c2 = len([m for m in candidate_motifs if m in rev_sig])
    failed_c3 = len([m for m in candidate_motifs if m not in consistent])
    failed_c4 = len([m for m in candidate_motifs if m in global_common])
    lines.append(
        f"Audit snapshot: total motifs evaluated={len(candidate_motifs)}, "
        f"failed criterion2={failed_c2}, failed criterion3={failed_c3}, failed criterion4={failed_c4}."
    )
    lines.append("")
    lines.append("---\n")

# Persist Section 5 motif-level audit table.
# Output: results_summary/high_confidence_criteria_audit.tsv
with open(HIGHCONF_AUDIT_PATH, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=highconf_audit_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(highconf_audit_rows)

lines.append("A per-motif criterion audit table has been written to:")
lines.append("`results_summary/high_confidence_criteria_audit.tsv`")
lines.append("")

# Section 6: Cross-Species
# Key biological question: do mouse-specific and human-specific regions
# use the same TFs? Shared motifs = convergent evolution or shared ancestry.
# Species-only motifs = candidate drivers of species-specific regulation.

lines.append("## 6. Cross-Species Comparison")
lines.append("")
lines.append("Do mouse-specific and human-specific regions share the same enriched TF motifs?")
lines.append("Comparison uses centered (summit200) forward-direction results.")
lines.append("")

for region_type in ["enh", "prom"]:
    # Cross-species comparison is performed on forward centered runs only.
    # This keeps region windowing consistent between species (summit200).
    mouse_dir = CROSS_SPECIES[f"mouse_{region_type}"]
    human_dir = CROSS_SPECIES[f"human_{region_type}"]
    m_res = parse_known_results(mouse_dir)
    h_res = parse_known_results(human_dir)
    m_sig = get_sig_motifs(m_res)
    h_sig = get_sig_motifs(h_res)

    # Set algebra:
    #   shared     = motifs passing thresholds in both species
    #   mouse_only = motifs passing thresholds only in mouse
    #   human_only = motifs passing thresholds only in human
    # Shared=0 therefore means "empty intersection under current thresholds",
    # not necessarily "no related TF biology at family level".
    shared = m_sig & h_sig
    mouse_only = m_sig - h_sig
    human_only = h_sig - m_sig

    label = "Enhancers" if region_type == "enh" else "Promoters"
    lines.append(f"### Species-Specific {label}")
    lines.append(f"- Mouse-specific enriched: **{len(m_sig)}** | Human-specific enriched: **{len(h_sig)}**")
    lines.append(f"- **Shared across species: {len(shared)}**")
    lines.append(f"- Mouse-only: {len(mouse_only)} | Human-only: {len(human_only)}")
    lines.append("")

    if shared:
        # Rank shared motifs by combined effect size across the two species.
        # This prioritizes motifs strong in both, not only one species.
        top_shared = sorted(shared, key=lambda m: -(m_res[m]["delta_pct"] + h_res[m]["delta_pct"]))[:10]
        lines.append("**Motifs enriched in BOTH species' specific regions:**")
        lines.append("")
        lines.append("| Motif | delta (mouse) | delta (human) |")
        lines.append("|-------|--------------|--------------|")
        for m in top_shared:
            lines.append(f"| {m[:55]} | {m_res[m]['delta_pct']:+.1f} | {h_res[m]['delta_pct']:+.1f} |")
        lines.append("")

    if mouse_only:
        # Show strongest mouse-only motifs by mouse delta effect size.
        top_mo = sorted(mouse_only, key=lambda m: -m_res[m]["delta_pct"])[:5]
        lines.append("**Mouse-only enriched (top 5):**")
        lines.append("")
        lines.append("| Motif | delta (mouse) | Family |")
        lines.append("|-------|--------------|--------|")
        for m in top_mo:
            lines.append(f"| {m[:55]} | {m_res[m]['delta_pct']:+.1f} | {extract_tf_family(m)} |")
        lines.append("")

    if human_only:
        # Show strongest human-only motifs by human delta effect size.
        top_ho = sorted(human_only, key=lambda m: -h_res[m]["delta_pct"])[:5]
        lines.append("**Human-only enriched (top 5):**")
        lines.append("")
        lines.append("| Motif | delta (human) | Family |")
        lines.append("|-------|--------------|--------|")
        for m in top_ho:
            lines.append(f"| {m[:55]} | {h_res[m]['delta_pct']:+.1f} | {extract_tf_family(m)} |")
        lines.append("")
    lines.append("---\n")

# Section 7: Enhancer vs Promoter
# Within the same species, enhancers and promoters often have different TF
# usage patterns. Enhancer-only motifs may indicate distal regulatory logic;
# promoter-only motifs may relate to core promoter elements.

lines.append("## 7. Enhancer vs Promoter Comparison")
lines.append("")
lines.append("Within the same species, do enhancers and promoters enrich different motifs?")
lines.append("")

for enh_dir, prom_dir, label in ENH_VS_PROM:
    e_res = parse_known_results(enh_dir)
    p_res = parse_known_results(prom_dir)
    e_sig = get_sig_motifs(e_res)
    p_sig = get_sig_motifs(p_res)
    shared = e_sig & p_sig
    enh_only = e_sig - p_sig
    prom_only = p_sig - e_sig

    lines.append(f"### {label}")
    lines.append(f"- Enhancer enriched: **{len(e_sig)}** | Promoter enriched: **{len(p_sig)}**")
    lines.append(f"- **Shared: {len(shared)}** | Enhancer-only: {len(enh_only)} | Promoter-only: {len(prom_only)}")
    lines.append("")

    if enh_only:
        top_eo = sorted(enh_only, key=lambda m: -e_res[m]["delta_pct"])[:5]
        lines.append("**Enhancer-only (top 5):**")
        lines.append("")
        lines.append("| Motif | delta | Family |")
        lines.append("|-------|-------|--------|")
        for m in top_eo:
            lines.append(f"| {m[:55]} | {e_res[m]['delta_pct']:+.1f} | {extract_tf_family(m)} |")
        lines.append("")

    if prom_only:
        top_po = sorted(prom_only, key=lambda m: -p_res[m]["delta_pct"])[:5]
        lines.append("**Promoter-only (top 5):**")
        lines.append("")
        lines.append("| Motif | delta | Family |")
        lines.append("|-------|-------|--------|")
        for m in top_po:
            lines.append(f"| {m[:55]} | {p_res[m]['delta_pct']:+.1f} | {extract_tf_family(m)} |")
        lines.append("")
    lines.append("---\n")

# Section 8: TF Family Summary
# Individual HOMER motifs are specific variants (e.g. Etv2, ETS1, ERG are
# all ETS family). Grouping by family avoids over-counting and reveals
# which TF families dominate across all experiments.

lines.append("## 8. TF Family Summary")
lines.append("")
lines.append("Motifs grouped by TF family. Count = number of distinct motifs from that")
lines.append("family that are significant across all 20 experiments (union).")
lines.append("")

family_counts = defaultdict(set)
for d in all_dirs:
    res = parse_known_results(d)
    sig = get_sig_motifs(res)
    for m in sig:
        fam = extract_tf_family(m)
        family_counts[fam].add(m)

sorted_families = sorted(family_counts.items(), key=lambda x: -len(x[1]))

lines.append("| TF Family | # Distinct Motifs | Example Motifs |")
lines.append("|-----------|-------------------|----------------|")
for fam, motifs in sorted_families[:25]:
    examples = ", ".join(sorted(motifs)[:3])
    if len(examples) > 70:
        examples = examples[:67] + "..."
    lines.append(f"| {fam} | {len(motifs)} | {examples} |")
lines.append("")

# Section 9: Cross-Experiment Matrix
# Binary presence/absence matrix: rows = motifs, columns = experiments.
# 1 = motif is significant in that experiment, 0 = not.
# This TSV can be imported into R or Python for heatmap visualization.
# We also report the most ubiquitous (appear in many experiments) and
# most specific (appear in ≤2 experiments) motifs.

lines.append("---")
lines.append("## 9. Cross-Experiment Motif Presence Matrix")
lines.append("")
lines.append("A TSV file has been written to `results_summary/motif_presence_matrix.tsv`.")
lines.append("Each row is a motif, each column is an experiment. Values: 1 = significant, 0 = not.")
lines.append("Use this for heatmap visualization.")
lines.append("")

# Build the matrix data
all_sig_union = set()
exp_sig_dict = {}
for d in all_dirs:
    res = parse_known_results(d)
    sig = get_sig_motifs(res)
    exp_sig_dict[d] = sig
    all_sig_union.update(sig)

sorted_motifs = sorted(all_sig_union)

# Write the matrix TSV
matrix_path = os.path.join(OUT_DIR, "motif_presence_matrix.tsv")
with open(matrix_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Motif"] + all_dirs)
    for m in sorted_motifs:
        row = [m] + [1 if m in exp_sig_dict[d] else 0 for d in all_dirs]
        writer.writerow(row)

lines.append(f"Matrix dimensions: **{len(sorted_motifs)} motifs × {len(all_dirs)} experiments**")
lines.append("")

# Most ubiquitous motifs (significant in the most experiments)
exp_counts = [(m, sum(1 for d in all_dirs if m in exp_sig_dict[d])) for m in sorted_motifs]
exp_counts.sort(key=lambda x: -x[1])
lines.append("**Most ubiquitous motifs** (significant in the most experiments):")
lines.append("")
lines.append("| Motif | # Experiments (of 20) | Family |")
lines.append("|-------|-----------------------|--------|")
for m, cnt in exp_counts[:15]:
    lines.append(f"| {m[:55]} | {cnt} | {extract_tf_family(m)} |")
lines.append("")

# Most specific motifs (significant in ≤2 experiments)
exp_counts_rare = [x for x in exp_counts if x[1] <= 2]
exp_counts_rare.sort(key=lambda x: x[1])
if exp_counts_rare:
    lines.append("**Most specific motifs** (significant in ≤2 experiments):")
    lines.append("")
    lines.append("| Motif | # Experiments | Family |")
    lines.append("|-------|---------------|--------|")
    for m, cnt in exp_counts_rare[:15]:
        which = [d for d in all_dirs if m in exp_sig_dict[d]]
        lines.append(f"| {m[:45]} | {cnt} ({', '.join(w[:25] for w in which)}) | {extract_tf_family(m)} |")
    lines.append("")

# Section 10: De Novo vs Known
# Cross-validation: if a de novo motif's best match to a known motif also
# appears in the significant known motif list for the same experiment,
# that's mutual support — the de novo discovery independently confirms
# the known motif enrichment.

lines.append("---")
lines.append("## 10. De Novo vs Known Cross-Validation")
lines.append("")
lines.append("For each experiment, check if the de novo motifs' best known-motif matches")
lines.append("also appear in the significant known motif list. Overlap = mutual support.")
lines.append("")

for d in all_dirs:
    known_res = parse_known_results(d)
    known_sig = get_sig_motifs(known_res)

    # Get de novo best matches from the HTML report
    denovo_matches = parse_denovo_html_bestmatch(d)
    if not denovo_matches:
        continue

    # Check if any de novo best match overlaps with the significant known list
    confirmed = []
    for dm in denovo_matches:
        for ks in known_sig:
            if dm.lower() in ks.lower() or ks.lower().startswith(dm.lower().split("/")[0]):
                confirmed.append((dm, ks))
                break

    lines.append(f"### {d}")
    lines.append(f"- De novo best matches parsed: {len(denovo_matches)}")
    lines.append(f"- Matches also in known significant list: **{len(confirmed)}**")
    if confirmed:
        lines.append("")
        lines.append("| De Novo Best Match | Confirmed Known Motif |")
        lines.append("|--------------------|----------------------|")
        for dm, ks in confirmed[:10]:
            lines.append(f"| {dm[:45]} | {ks[:50]} |")
    lines.append("")


#
# Write final report
#
report_path = os.path.join(OUT_DIR, "integrated_motif_report.md")
with open(report_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

print(f"Done. Integrated report: {report_path}")
print(f"  Motif presence matrix: {matrix_path}")
print(f"  High-confidence audit: {HIGHCONF_AUDIT_PATH}")

# Reference
# LLM was used for debugging
# [1] OpenAI, "GPT-5.4 Thinking System Card," Mar. 5, 2026. Accessed: Apr. 13, 2026.
