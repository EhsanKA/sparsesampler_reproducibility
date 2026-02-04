# Downstream Analysis Comparison for Reviewer Response

This directory contains scripts to address reviewer comments about:
1. **Cell type/state diversity preservation** at different sampling depths (instead of just rare cell counts)
2. **Downstream analysis comparison** between original and subsampled populations:
   - Differential expression (DE) gene overlap
   - Clustering/cell-type annotation results comparison

## Files

### Analysis Scripts

1. **`downstream_analysis_comparison.py`**
   - Main analysis script that:
     - Calculates diversity metrics (Shannon index, Simpson index, cell type counts)
     - Compares DE genes between original and subsampled populations
     - Compares clustering results (ARI, NMI) between original and subsampled populations
   - Outputs CSV files with results

2. **`plot_downstream_analysis.py`**
   - Generates figures showing:
     - Diversity preservation across reference sizes
     - Diversity preservation across sample sizes
     - DE gene overlap comparisons
     - Clustering similarity metrics (ARI, NMI)
   - Creates combined figure for reviewer

3. **`run_downstream_analysis.sh`**
   - Convenience script to run both analysis and plotting

## Usage

### Quick Start

```bash
# Set project root (if not already set)
export PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"

# Run the analysis
cd $PROJECT_ROOT/analysis
bash run_downstream_analysis.sh
```

### Manual Execution

```bash
# Run analysis
python downstream_analysis_comparison.py

# Generate figures
python plot_downstream_analysis.py
```

## Output

### Results (CSV files)
Located in: `analysis/results/downstream_analysis/`

- `{dataset}_diversity_preservation.csv` - Diversity metrics for each dataset
- `combined_diversity_preservation.csv` - Combined results across datasets
- `{dataset}_downstream_comparison.csv` - DE and clustering comparison results
- `combined_downstream_comparison.csv` - Combined downstream analysis results

### Figures
Located in: `figures/revision/downstream_analysis/`

- `{dataset}_diversity_preservation.jpg` - Diversity preservation metrics
- `{dataset}_diversity_by_sample_size.jpg` - Diversity vs sample size
- `{dataset}_downstream_comparison.jpg` - DE and clustering comparisons
- `reviewer_response_combined_figure.jpg` - Combined figure for reviewer

## Metrics Calculated

### Diversity Metrics
1. **Number of Cell Types Preserved**: Fraction of original cell types found in subsample
2. **Shannon Diversity Index**: Measures diversity considering both richness and evenness
3. **Simpson Diversity Index**: Measures probability that two randomly selected cells are of different types
4. **Cell Type Overlap**: Fraction of original cell types present in subsample

### Downstream Analysis Metrics
1. **DE Gene Overlap**: Fraction of top DE genes that overlap between original and subsampled
2. **Adjusted Rand Index (ARI)**: Clustering similarity metric (0-1, higher is better)
3. **Normalized Mutual Information (NMI)**: Clustering similarity metric (0-1, higher is better)
4. **Cell Type Annotation Overlap**: Fraction of cell types correctly identified in subsample

## Datasets Analyzed

- **mcc_01**: Main dataset with 0.1% frequency threshold for rare cells
- **mcc_05**: Dataset with 0.5% frequency threshold for rare cells

## Methods Compared

- `random`: Random sampling (baseline)
- `sps`: SparseSampler method
- `hopper`: Hopper method
- `atomic`: Atomic method
- `scsampler`: scSampler method

## Notes

- The analysis uses rep 0 by default (can be extended to multiple reps)
- DE analysis is performed on a subset of cell types (first 5) for efficiency
- Clustering uses Leiden algorithm with default resolution
- The downstream comparison analysis is run on a subset of configurations to save computation time

## Reviewer Response Summary

This analysis addresses the reviewer's concerns by:

1. **Diversity Preservation**: Shows how well cell type diversity is preserved at different sampling depths using multiple metrics (Shannon, Simpson, cell type counts)

2. **DE Gene Comparison**: Demonstrates that subsampled populations preserve the key differential expression patterns found in the original data

3. **Clustering Comparison**: Shows that clustering results from subsampled data are similar to those from the original data (using ARI and NMI metrics)

These analyses prove the practical value of the subsampling method by showing that downstream analyses remain valid when performed on subsampled populations.







