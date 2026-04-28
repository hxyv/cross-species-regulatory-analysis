# Automated Pipeline

This folder contains the single-command workflow for Tasks 2 to 5.

Before running the pipeline, copy the root config template and edit the paths for your environment:

```bash
cp pipeline.conf.example pipeline.conf
```

The pipeline reads configuration from the root-level `pipeline.conf` file. Task 2 must succeed before Tasks 3 to 5 can run, and the full pipeline stages those dependencies automatically.

Run on Bridges-2 with the configured inputs:

```bash
sbatch 06.pipeline/run_adrenal_pipeline.slurm
```

Run directly with a config file:

```bash
bash 06.pipeline/run_adrenal_pipeline.sh \
  --config pipeline.conf
```

The pipeline performs preflight checks for required files, commands, HAL/HALPER locations, and Python modules before starting the analysis. It then stages intermediate files between task folders and runs the downstream workflow in project order.

---

[Go back to main.](https://github.com/BioinformaticsDataPracticum2026/cross-species-regulatory-analysis#usage-step-by-step)
