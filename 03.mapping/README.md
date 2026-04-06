# Task 2 Workflow

## Goal

Map human adrenal ATAC-seq peaks to mouse and identify shared versus non-overlapping regions in mouse coordinates.

## Input peak sets

- Human species-level peaks:
  `/ocean/projects/bio230007p/ikaplow/HumanAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz`
- Mouse species-level peaks:
  `/ocean/projects/bio230007p/ikaplow/MouseAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz`

`idr.optimal_peak` was used as the replicate-aware species-level peak set for each species.

## Step 1: preprocess peak files

Script:

- `03.mapping/prepare_adrenal_mapping_preprocess.sh`

Outputs:

- `03.mapping/human_adrenal_idr_optimal.mapping_preprocess.bed.gz`
- `03.mapping/mouse_adrenal_idr_optimal.mapping_preprocess.bed.gz`

This step creates project-local symlinks, normalizes the name column, and sorts the intervals.

## Step 2: map human peaks to mouse with HALPER

Script:

- `03.mapping/run_adrenal_halper_mapping.sh`

Alignment:

- `/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal`

Main output:

- `03.mapping/human_adrenal_idr_optimal.HumanToMouse.HALPER.narrowPeak.gz`

Supporting outputs:

- `03.mapping/human_adrenal_idr_optimal.HumanToMouse.halLiftover.sFile.bed.gz`
- `03.mapping/human_adrenal_idr_optimal.HumanToMouse.halLiftover.tFile.bed.gz`

This workflow uses one HALPER direction only:

- `Human -> Mouse`

## Step 3: intersect mapped human peaks with native mouse peaks

Script:

- `03.mapping/run_adrenal_bedtools_intersection.sh`

Tool:

- `bedtools intersect`

Outputs:

- `03.mapping/human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed.gz`
- `03.mapping/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed.gz`
- `03.mapping/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed.gz`
- `03.mapping/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed.gz`
- `03.mapping/adrenal_human_to_mouse_intersection_summary.tsv`

The current overlap rule is default `bedtools intersect`, so any overlap counts.

## Summary counts

From `03.mapping/adrenal_human_to_mouse_intersection_summary.tsv`:

- mapped human total: `99076`
- native mouse total: `48263`
- mapped human shared with mouse: `35321`
- mapped human no mouse overlap: `63755`
- mouse shared with human mapped: `27967`
- mouse no human mapped overlap: `20296`

## Coordinate system

All final comparison sets are in mouse coordinates.

## Most useful Task 2 outputs

- shared set:
  `03.mapping/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed.gz`
- mouse-specific set:
  `03.mapping/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed.gz`
- mapped human-only set:
  `03.mapping/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed.gz`
