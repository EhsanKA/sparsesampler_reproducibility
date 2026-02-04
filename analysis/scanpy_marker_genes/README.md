# Scanpy Marker Gene Analysis

This directory contains scripts for identifying marker genes for cell types in the MCC dataset using Scanpy's `rank_genes_groups`, comparing results between **Full dataset**, **SPS 100k**, and **Random 100k** samples.

## Overview

The analysis demonstrates that **SPS sampling better preserves marker gene discovery** compared to random sampling. This is evaluated by:

1. Using Scanpy's differential expression methods to identify marker genes
2. Comparing marker genes across sampling methods (Full vs SPS vs Random)
3. Evaluating recovery of **literature-based marker genes**
4. Comparing multiple DE methods (Wilcoxon, t-test, logistic regression)

## Dataset

- **MCC Dataset**: ~3.16 million cells
- **Cell Types**:
  - mesodermal cell (1,463,100 cells)
  - lateral mesodermal cell (745,494 cells)
  - fibroblast (648,047 cells)
  - chondrocyte (274,756 cells)
  - osteoblast (30,000 cells - **rare**)

## Methods Compared

### Sampling Methods

| Method | Description |
|--------|-------------|
| **Full** | Full MCC dataset (~3.16M cells) - ground truth |
| **SPS 100k** | 100,000 cells selected by SparseSampler (rep 0) |
| **Random 100k** | 100,000 randomly selected cells (rep 0) |

### DE Methods

| Method | Description |
|--------|-------------|
| **Wilcoxon** | Wilcoxon rank-sum test (non-parametric, default) |
| **t-test** | Welch's t-test (parametric) |
| **logreg** | Logistic regression (multivariate) |

## Directory Structure

```
analysis/scanpy_marker_genes/
├── README.md                          # This file
├── scanpy_marker_pipeline.py          # Main analysis script
├── compare_marker_genes.py            # Comparison and visualization script
├── run_marker_gene_analysis.sh        # SLURM job submission script
└── results/
    ├── {method}_{de_method}_marker_genes.csv
    ├── {method}_{de_method}_marker_genes.pkl
    ├── {method}_{de_method}_literature_evaluation.csv
    ├── {method}_{de_method}_literature_evaluation.pkl
    ├── comparison_{de_method}.csv
    ├── comparison_all.csv
    ├── overlap_{de_method}.pkl
    ├── correlations_{de_method}.pkl
    └── figures/
        ├── recall_comparison_{de_method}.png
        ├── rank_comparison_{de_method}.png
        ├── sps_improvement_{de_method}.png
        ├── overlap_heatmap_{de_method}.png
        ├── de_comparison_{method}.png
        └── combined_recall_heatmap.png
```

## Prerequisites

### Conda Environment

Use the `facs_sampling` environment:

```bash
conda activate facs_sampling
```

This environment should have: scanpy, pandas, numpy, scipy, matplotlib, seaborn

### Data Paths

- MCC dataset: `/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/data/mcc/benchmark/30/`
- SPS indices: `.../sps/100000/0/results.pkl`
- Random indices: `.../random/100000/0/results.pkl`

## Usage

### Option 1: Submit SLURM Jobs (Recommended)

```bash
cd /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/scanpy_marker_genes

# Run all analyses (full, sps, random) with all DE methods
bash run_marker_gene_analysis.sh all

# After all jobs complete, run comparison
bash run_marker_gene_analysis.sh compare
```

Or run individual methods:

```bash
# Run only full dataset
bash run_marker_gene_analysis.sh full

# Run only SPS subsample
bash run_marker_gene_analysis.sh sps

# Run only Random subsample
bash run_marker_gene_analysis.sh random
```

### Option 2: Run Locally

```bash
# Run everything locally (takes longer)
bash run_marker_gene_analysis.sh local
```

### Option 3: Run Interactively

```bash
# Activate environment
conda activate facs_sampling

# Run analysis for each sampling method
python scanpy_marker_pipeline.py --method full --de-methods wilcoxon t-test logreg
python scanpy_marker_pipeline.py --method sps --de-methods wilcoxon t-test logreg
python scanpy_marker_pipeline.py --method random --de-methods wilcoxon t-test logreg

# Run comparison
python compare_marker_genes.py --de-methods wilcoxon t-test logreg
```

## Pipeline Details

### Step 1: Data Loading
- Loads MCC dataset (~11GB)
- Loads sampling indices (for SPS/Random methods)
- Subsets to sampled cells

### Step 2: Marker Gene Identification
For each DE method:
- Runs `sc.tl.rank_genes_groups` with one-vs-rest comparison
- Extracts scores, p-values, log fold changes
- Reports top 50 marker genes per cell type

### Step 3: Literature Marker Evaluation
- Compares identified markers against literature-based markers
- Computes:
  - **Recall@10**: % of literature markers in top 10 genes
  - **Recall@50**: % of literature markers in top 50 genes
  - **Mean Rank**: Average rank of literature markers

### Step 4: Comparison (Full vs SPS vs Random)
- Computes improvement metrics for SPS over Random
- Calculates marker gene overlap (Jaccard similarity)
- Computes Spearman rank correlations
- Generates visualizations

## Literature Marker Sources

Marker genes are sourced from (same as SEMITONES analysis):

1. **CellMarker 2.0** (Hu et al., 2023)
   - DOI: 10.1093/nar/gkac947

2. **PanglaoDB** (Franzen et al., 2019)
   - DOI: 10.1093/database/baz046

3. **Lateral Mesoderm** (Loh et al., 2016)
   - DOI: 10.1016/j.cell.2016.06.011

See `../semitones_marker_genes/MARKER_GENE_REFERENCES.md` for full citations.

## Expected Results

**Hypothesis**: SPS sampling should show:
- Higher recall of literature markers than random sampling
- Better (lower) mean ranks for literature markers
- Especially strong improvement for rare cell types (osteoblast)

**Key Metrics**:
- Recall@10, Recall@50 per cell type
- Mean rank of literature markers
- Jaccard similarity of top marker genes

## Comparison with SEMITONES

| Aspect | SEMITONES | Scanpy |
|--------|-----------|--------|
| Method | Enrichment scoring via similarity | Differential expression testing |
| Reference cells | Sampled cells as references | One-vs-rest comparison |
| Statistics | Linear regression on distance-expression | Wilcoxon/t-test/logreg |
| Speed | Slower (pairwise similarities) | Faster (vectorized) |
| Conda env | `semitones_env` | `facs_sampling` |

## SLURM Configuration

| Job | Memory | CPUs | Time |
|-----|--------|------|------|
| Full dataset | 128GB | 4 | 8 hours |
| SPS/Random | 64GB | 4 | 4 hours |
| Comparison | 16GB | 2 | 1 hour |

- **Partition**: normal
- **Account**: ohler

## Output Files

| File | Description |
|------|-------------|
| `*_marker_genes.pkl` | Full marker gene results (pickle) |
| `*_marker_genes.csv` | Top 50 marker genes per cell type |
| `*_literature_evaluation.pkl` | Recall and rank metrics (pickle) |
| `*_literature_evaluation.csv` | Recall and rank metrics vs literature |
| `comparison_*.csv` | SPS vs Random comparison tables |
| `figures/*.png` | Visualization plots |

## Troubleshooting

### Import errors for literature_markers
Ensure the script can find the semitones_marker_genes module:
```bash
cd /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/scanpy_marker_genes
python -c "from semitones_marker_genes.literature_markers import LITERATURE_MARKERS; print('OK')"
```

### Memory issues
- Increase SLURM memory request
- For full dataset, use at least 128GB

### Missing results for comparison
Ensure all three analyses (full, sps, random) completed successfully before running comparison.

## References

- Scanpy: Wolf F.A., et al. (2018). Genome Biology. DOI: 10.1186/s13059-017-1382-0
- CellMarker 2.0: Hu et al. (2023). Nucleic Acids Research. DOI: 10.1093/nar/gkac947
- PanglaoDB: Franzen et al. (2019). Database. DOI: 10.1093/database/baz046
