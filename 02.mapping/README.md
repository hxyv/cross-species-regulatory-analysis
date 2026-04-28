# Task 2 Workflow

## Goal

Map human adrenal ATAC-seq peaks to mouse and identify shared versus non-overlapping regions in mouse coordinates.

## Input peak sets

These scripts can read shared path settings from the root-level `pipeline.conf` file. For a new environment, copy and edit:

```bash
cp pipeline.conf.example pipeline.conf
```

- Human species-level peaks:
  `/ocean/projects/bio230007p/ikaplow/HumanAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz`
- Mouse species-level peaks:
  `/ocean/projects/bio230007p/ikaplow/MouseAtac/AdrenalGland/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz`

`idr.optimal_peak` was used as the replicate-aware species-level peak set for each species.

## Step 1: preprocess peak files

Script:

- `02.mapping/prepare_adrenal_mapping_preprocess.sh`

Outputs:

- `02.mapping/human_adrenal_idr_optimal.mapping_preprocess.bed.gz`
- `02.mapping/mouse_adrenal_idr_optimal.mapping_preprocess.bed.gz`

This step creates project-local symlinks, normalizes the name column, and sorts the intervals.

## Step 2: map human peaks to mouse with HALPER

Script:

- `02.mapping/run_adrenal_halper_mapping.sh`

Alignment:

- `/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal`

Main output:

- `02.mapping/human_adrenal_idr_optimal.HumanToMouse.HALPER.narrowPeak.gz`

Supporting outputs:

- `02.mapping/human_adrenal_idr_optimal.HumanToMouse.halLiftover.sFile.bed.gz`
- `02.mapping/human_adrenal_idr_optimal.HumanToMouse.halLiftover.tFile.bed.gz`

This workflow uses one HALPER direction only:

- `Human -> Mouse`

Task 2 must complete successfully before downstream biological process, promoter/enhancer, or motif analyses are run, because those later steps depend on the shared and species-specific peak sets produced here.

## Step 3: intersect mapped human peaks with native mouse peaks

Script:

- `02.mapping/run_adrenal_bedtools_intersection.sh`

Tool:

- `bedtools intersect`

Outputs:

- `02.mapping/human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed.gz`
- `02.mapping/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed.gz`
- `02.mapping/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed.gz`
- `02.mapping/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed.gz`
- `02.mapping/adrenal_human_to_mouse_intersection_summary.tsv`

The current overlap rule is default `bedtools intersect`, so any overlap counts.

## Step 4: recover human-coordinate peak sets for human analyses

The intersection files above are in mouse coordinates because the workflow maps human peaks to mouse before comparing them with native mouse peaks. Downstream human motif analyses need human-coordinate intervals, so this step maps the shared and human-specific peak names back to the original preprocessed human peak file.

Script:

- `02.mapping/recover_human_coordinate_peak_sets.sh`

Inputs:

- original human peaks:
  `02.mapping/human_adrenal_idr_optimal.mapping_preprocess.bed.gz`
- mapped human peaks with mouse overlap:
  `02.mapping/human_adrenal_idr_optimal.HumanToMouse.shared_with_mouse_native.bed.gz`
- mapped human peaks without mouse overlap:
  `02.mapping/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed.gz`

Outputs:

- human-specific peaks in original human coordinates:
  `02.mapping/human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed.gz`
- shared human peaks in original human coordinates:
  `02.mapping/human_adrenal_idr_optimal.shared.original_human_coordinates.bed.gz`
- recovery summary:
  `02.mapping/human_coordinate_recovered_peak_sets_summary.tsv`

## Summary counts

From `02.mapping/adrenal_human_to_mouse_intersection_summary.tsv`:

- mapped human total: `99076`
- native mouse total: `48263`
- mapped human shared with mouse: `35321`
- mapped human no mouse overlap: `63755`
- mouse shared with human mapped: `27967`
- mouse no human mapped overlap: `20296`

From `02.mapping/human_coordinate_recovered_peak_sets_summary.tsv`:

- human-specific peaks recovered in original human coordinates: `63755`
- shared human peaks recovered in original human coordinates: `35321`

## Coordinate system

The main Task 2 comparison sets are in mouse coordinates. The two recovered human-coordinate files are only for downstream analyses that specifically require original human coordinates, such as human enhancer/promoter motif analysis.

## Most useful Task 2 outputs in mouse coordinates

- shared set:
  `02.mapping/mouse_adrenal_idr_optimal.shared_with_human_mapped.bed.gz`
- mouse-specific set:
  `02.mapping/mouse_adrenal_idr_optimal.no_human_mapped_overlap.bed.gz`
- mapped human-only set:
  `02.mapping/human_adrenal_idr_optimal.HumanToMouse.no_mouse_native_overlap.bed.gz`

## Human-coordinate handoff outputs

- human-specific peaks in original human coordinates:
  `02.mapping/human_adrenal_idr_optimal.human_specific.original_human_coordinates.bed.gz`
- shared human peaks in original human coordinates:
  `02.mapping/human_adrenal_idr_optimal.shared.original_human_coordinates.bed.gz`


---

[Go back to main.](https://github.com/BioinformaticsDataPracticum2026/cross-species-regulatory-analysis#usage-step-by-step)
