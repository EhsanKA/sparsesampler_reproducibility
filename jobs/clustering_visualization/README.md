# Clustering Visualization Jobs

This directory contains scripts to visualize and compare UMAP and Leiden clustering results on sampled data (SPS vs Random).

## Overview

The visualization script (`visualize_sampled_clustering.py`) performs the following:

1. **Loads reference data** and sampling indices for both SPS and Random methods
2. **Subsets data** to only sampled cells for each method
3. **Computes UMAP** on each sampled dataset separately
4. **Computes Leiden clustering** on each sampled dataset
5. **Compares clustering quality** using ground truth cell types (ARI, NMI, cluster purity)
6. **Generates visualization figures** showing:
   - UMAP plots colored by cell type (SPS vs Random)
   - UMAP plots colored by Leiden clusters (SPS vs Random)
   - Confusion matrices (cell type vs cluster)
   - Metrics comparison bar plots

## Files

- `01_parallel_jobs_clustering.sh` - SLURM job array script to run analysis for multiple combinations
- `visualize_sampled_clustering.py` - Main Python script (located in `analysis/` directory)

## Usage

### Running via SLURM (recommended)

```bash
cd jobs/clustering_visualization
bash 01_parallel_jobs_clustering.sh
```

This will submit a job array for:
- Datasets: mcc, lcmv
- Sizes: 100000, 200000
- Reps: 0
- Total: 4 combinations

### Running directly (for testing)

```bash
cd analysis
python visualize_sampled_clustering.py \
    --dataset mcc \
    --ref 30 \
    --size 100000 \
    --rep 0 \
    --n-neighbors 15 \
    --resolution 0.5
```

### Parameters

- `--dataset`: Dataset name (mcc or lcmv)
- `--ref`: Reference size (30 for mcc, 34 for lcmv)
- `--size`: Sample size (e.g., 100000, 200000)
- `--rep`: Replication number (default: 0)
- `--n-neighbors`: Number of neighbors for UMAP (default: 15)
- `--resolution`: Leiden clustering resolution (default: 0.5)

## Output

### Metrics CSV
Location: `analysis/results/sampling_clustering/{dataset}_size{size}_rep{rep}_metrics.csv`

Contains:
- ARI (Adjusted Rand Index)
- NMI (Normalized Mutual Information)
- Average cluster purity
- Number of clusters
- Number of cell types

### Visualization Figures
Location: `figures/sampling_clustering/{dataset}_size{size}_rep{rep}_clustering_comparison.png`

Shows:
- UMAP plots for SPS and Random (cell types and clusters)
- Confusion matrices
- Metrics comparison

## Metrics Explained

- **ARI (Adjusted Rand Index)**: Measures agreement between clustering and ground truth (range: -1 to 1, higher is better)
- **NMI (Normalized Mutual Information)**: Measures mutual information between clustering and ground truth (range: 0 to 1, higher is better)
- **Cluster Purity**: Average proportion of the most common cell type in each cluster (range: 0 to 1, higher is better)





