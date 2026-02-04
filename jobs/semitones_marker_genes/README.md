# SEMITONES Marker Gene Analysis Jobs

This directory contains job scripts for running SEMITONES marker gene analysis on the MCC dataset.

## Overview

The analysis compares marker gene rankings between:
- **Full MCC dataset** (ground truth)
- **SPS subsample** (100K cells)
- **Random subsample** (100K cells)

For each rare cell type, we demonstrate that SPS preserves marker gene signal better than random sampling.

## Configuration

- **Dataset**: mcc (ref=30)
- **Sample size**: 100,000 cells
- **Methods**: sps, random
- **Environment**: facs_sampling

## Running the Analysis

```bash
cd /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/jobs/semitones_marker_genes
bash run_marker_gene_analysis.sh
```

## Output

Results are saved to:
- `analysis/semitones_marker_genes/results/mcc_marker_comparison_metrics.csv`
- `analysis/semitones_marker_genes/results/mcc_marker_comparison_figure.png`
- `analysis/semitones_marker_genes/results/mcc_{celltype}_{source}_marker_ranking.csv`

## Monitoring

```bash
squeue -u $USER
```

## Logs

Logs are saved to:
- `jobs/semitones_marker_genes/logs/output_marker_gene_*.stdlog`
- `jobs/semitones_marker_genes/logs/output_marker_gene_*.stderr`
