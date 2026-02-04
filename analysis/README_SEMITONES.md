# SEMITONES Analysis

This directory contains scripts for using SEMITONES (Similarity Enrichment for Multiple Interrogation Techniques in Omics via Neighborhoods of Exemplars) to evaluate how well sampled cells represent the full dataset.

## Setup

SEMITONES requires Python < 3.10. A dedicated conda environment has been created:

```bash
conda activate semitones_env
```

If the environment doesn't exist or SEMITONES is not installed:

```bash
conda create -n semitones_env python=3.9 -y
conda activate semitones_env
pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip scanpy pandas numpy scipy
```

## Usage

### Running on SLURM Cluster (Recommended)

For production runs, use the SLURM job scripts:

```bash
# Test with a single job first
cd jobs/semitones_analysis
bash 00_test_single_job.sh

# Run full job array (all combinations)
bash 01_parallel_jobs_semitones.sh
```

See `jobs/semitones_analysis/README.md` for details.

### Quick Example (Local Testing)

Try the simple example script first to verify everything works:

```bash
conda activate semitones_env
python example_semitones_usage.py
```

### Basic Usage (Local)

Run SEMITONES analysis on a single configuration:

```bash
conda activate semitones_env
python semitones_analysis.py \
    --dataset mcc_01 \
    --ref 30 \
    --method sps \
    --size 200000 \
    --rep 0
```

### Parameters

- `--dataset`: Dataset name (`mcc_01`, `mcc_05`, or `lcmv`)
- `--ref`: Reference size (e.g., 30 for mcc_01, 34 for lcmv)
- `--method`: Sampling method (`random`, `sps`, `hopper`, `atomic`, `scsampler`)
- `--size`: Sample size (e.g., 50000, 100000, 200000, 300000)
- `--rep`: Replication number (default: 0)

### What the Script Does

1. **Loads the full dataset** from the benchmark directory
2. **Loads the sampled cell indices** (these become the reference cells)
3. **Calculates similarity scores** between all query cells (non-sampled) and reference cells (sampled)
4. **Calculates enrichment scores** using SEMITONES
5. **Analyzes cell type coverage** to see how well each cell type is represented in the reference
6. **Saves results** to `analysis/results/semitones_analysis/`

### Output

The script generates:
- A pickle file with similarity scores, enrichment scores, and coverage metrics
- Log file in `logs/semitones_analysis_*.log`

### Results Structure

The results dictionary contains:
- `similarity_scores`: Matrix of similarity scores (n_query × n_reference)
- `enrichment_scores`: Enrichment scores from SEMITONES
- `coverage_metrics`: Per-cell-type coverage statistics
- `type_similarities`: Average similarity scores per cell type
- Summary statistics (mean, median, max, min similarity)

## Example: Comparing Sampling Methods

To compare how well different sampling methods represent the dataset:

```bash
# SPS method
python semitones_analysis.py --dataset mcc_01 --ref 30 --method sps --size 200000 --rep 0

# Random method
python semitones_analysis.py --dataset mcc_01 --ref 30 --method random --size 200000 --rep 0

# Hopper method
python semitones_analysis.py --dataset mcc_01 --ref 30 --method hopper --size 200000 --rep 0
```

## Understanding the Results

- **Higher similarity scores** indicate that query cells are well-represented by the reference cells
- **Cell type coverage metrics** show which cell types are well/poorly represented in the reference
- **Type similarities** help identify which cell types have good representation in the sampled reference

## Notes

- The script uses the sampled cells as **references** and all other cells as **queries**
- SEMITONES calculates similarity based on gene expression profiles
- The analysis helps evaluate how well the sampling methods preserve the diversity of the full dataset

