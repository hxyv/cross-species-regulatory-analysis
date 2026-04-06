# Cross-Species Regulatory Analysis

Course project repository for 03-713 Bioinformatics Data Integration Practicum (Spring 2026). This team analyzes cross-species conservation of transcriptional regulatory activity using human and mouse open chromatin data from adrenal gland and ovary.

## Project Overview

The project is designed around six required analyses from the course project description:

1. Evaluate the quality of the human and mouse datasets for both assigned tissues.
2. Map open chromatin regions between species and identify conserved versus species-specific accessibility.
3. Identify candidate biological processes linked to open chromatin regions in each species, in shared regions, and in species-specific regions.
4. Classify open chromatin regions into likely enhancers and promoters, then compare their conservation across species.
5. Discover sequence motifs enriched in relevant peak sets.
6. Build an automated command-line pipeline for Tasks 2 to 5 that runs on a Linux cluster.

Task 1 should be completed for both tissues. Downstream analyses may focus on the tissue with the best-quality human and mouse datasets after the QC stage.

## Team Responsibilities

- Chester Xiao: Task 1 dataset QC for adrenal gland and ovary; primary workspace in `02.qc/`
- Xingyu Hu: Task 2 cross-species mapping and regulatory conservation; primary workspace in `03.mapping/`
- Lekhya Dommalapati: Task 3 biological process analysis; primary workspace in `04.biological_processes/`
- Maitreyee Karne: Task 4 enhancer/promoter classification and conservation comparison; primary workspace in `05.promoter_enhancer/`
- All team members: Task 5 motif analysis in `06.motifs/` and Task 6 pipeline development in `07.pipeline/`

## Repository Structure

This repository uses numbered top-level folders so the project reads like a bioinformatics workflow.

- `01.data/`: dataset manifests, sample sheets, metadata, Bridges-2 path notes, and small helper files. Large sequencing and reference files are not committed.
- `02.qc/`: quality-control summaries, plots, scripts, and the decision on which tissue pair to use downstream.
- `03.mapping/`: cross-species peak mapping, shared/open/closed region sets, and related code or outputs.
- `04.biological_processes/`: biological process enrichment analyses for each species, shared regions, and species-specific regions.
- `05.promoter_enhancer/`: promoter/enhancer classification and conservation comparison outputs.
- `06.motifs/`: motif enrichment analyses for enhancers, promoters, and shared or species-specific peak sets.
- `07.pipeline/`: cluster-runnable workflow code, configuration, and the single command entrypoint for Tasks 2 to 5.
- `08.results/`: polished figures, tables, and report-ready integrated outputs.

## Data And Compute Environment

Primary data and references are stored on Bridges-2 and should be referenced from this repository rather than copied into git.

- Human ATAC-seq data: `/ocean/projects/bio230007p/ikaplow/HumanAtac`
- Mouse ATAC-seq data: `/ocean/projects/bio230007p/ikaplow/MouseAtac`
- Human genome sequence and annotations: `/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo`
- Mouse genome sequence and annotations: `/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo`
- Multi-species alignment: `/ocean/projects/bio230007p/ikaplow/Alignments`
- TF motif data: `/ocean/projects/bio230007p/ikaplow/CIS-BP_2.00`

All jobs for this project should be designed to run on a Linux cluster environment such as Bridges-2.

## Workflow Overview

The expected analysis flow is:

`01.data -> 02.qc -> 03.mapping -> 04.biological_processes -> 05.promoter_enhancer -> 06.motifs -> 07.pipeline -> 08.results`

Typical progression:

1. Record sample metadata and dataset locations in `01.data/`.
2. Evaluate dataset quality for human and mouse adrenal gland and ovary in `02.qc/`.
3. Select the better tissue pair for downstream analysis if needed.
4. Perform cross-species mapping and define shared versus species-specific open chromatin sets in `03.mapping/`.
5. Run biological process enrichment in `04.biological_processes/`.
6. Split peaks into likely promoters and enhancers in `05.promoter_enhancer/`.
7. Run motif analyses in `06.motifs/`.
8. Consolidate reusable workflow code in `07.pipeline/` so Tasks 2 to 5 can run in a single command.
9. Copy report-ready summaries, figures, and tables to `08.results/`.

## Expected Tools

The course materials specifically point to the following tools and resources:

- `bedtools`
- `HALPER`
- `rGREAT`
- `HOMER`

Additional scripts may be written in Python, R, or shell as needed, but stable logic should move out of notebooks and into reusable scripts or pipeline code.

## Reproducibility Notes

- Do not commit large raw data, genome files, BAM files, or other heavy intermediates.
- Keep small metadata, configs, and method notes in the repo.
- Organize files directly inside each task folder in whatever simple layout the team finds easiest.
- Use `08.results/` for polished outputs that are ready for presentation or the final report.
- The final automated workflow should run Tasks 2 to 5 with a single command on a Linux cluster where required tools are installed.

---

LLM was used for debugging and version control:

[1] OpenAI, “GPT-5.4 Model,” OpenAI API. [Online]. Available: https://developers.openai.com/api/docs/models/gpt-5.4. [Accessed: Apr. 6, 2026].

[2] OpenAI, “Codex,” OpenAI Developers. [Online]. Available: https://developers.openai.com/codex/. [Accessed: Apr. 6, 2026].