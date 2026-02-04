# SEMITONES Marker Gene Analysis

This directory contains scripts for identifying marker genes for cell types in the MCC dataset using SEMITONES, comparing results when using **SPS 100k** vs **Random 100k** reference cells.

## Overview

The analysis demonstrates that **SPS sampling better captures cells with enriched marker genes** compared to random sampling. This is evaluated by:

1. Using SEMITONES to compute enrichment scores for all genes
2. Identifying marker genes per cell type based on enrichment
3. Comparing how well each sampling method recovers **literature-based marker genes**

## Dataset

- **MCC Dataset**: ~3.16 million cells
- **Cell Types**:
  - mesodermal cell (1,463,100 cells)
  - lateral mesodermal cell (745,494 cells)
  - fibroblast (648,047 cells)
  - chondrocyte (274,756 cells)
  - osteoblast (30,000 cells - **rare**)

## Methods Compared

| Method | Description |
|--------|-------------|
| **SPS 100k** | 100,000 cells selected by SparseSampler (rep 0) |
| **Random 100k** | 100,000 randomly selected cells (rep 0) |

## Directory Structure

```
analysis/semitones_marker_genes/
├── README.md                          # This file
├── MARKER_GENE_REFERENCES.md          # Literature references for marker genes
├── literature_markers.py              # Python module with marker gene definitions
├── semitones_marker_pipeline.py       # Main SEMITONES analysis script
├── compare_marker_genes.py            # Comparison and visualization script
├── run_marker_gene_analysis.sh        # SLURM job submission script
└── results/
    ├── sps_enrichment_scores.pkl      # SEMITONES enrichment scores (SPS)
    ├── random_enrichment_scores.pkl   # SEMITONES enrichment scores (Random)
    ├── sps_marker_genes.pkl           # Marker genes per cell type (SPS)
    ├── random_marker_genes.pkl        # Marker genes per cell type (Random)
    ├── sps_marker_genes_by_celltype.csv
    ├── random_marker_genes_by_celltype.csv
    ├── sps_literature_evaluation.csv  # Comparison to literature markers
    ├── random_literature_evaluation.csv
    ├── marker_gene_comparison.csv     # SPS vs Random comparison
    └── figures/
        ├── recall_comparison.png
        ├── rank_comparison.png
        ├── marker_overlap_heatmap.png
        ├── sps_improvement_summary.png
        └── literature_marker_recovery.png
```

## Prerequisites

### Conda Environments

1. **semitones_env**: For SEMITONES analysis
   ```bash
   conda activate semitones_env
   # If not installed:
   conda create -n semitones_env python=3.9 -y
   conda activate semitones_env
   pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip scanpy pandas numpy scipy
   ```

2. **facs_sampling**: For comparison/visualization (matplotlib, seaborn)

### Data Paths

- MCC dataset: `/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/data/mcc/benchmark/30/`
- SPS indices: `.../sps/100000/0/results.pkl`
- Random indices: `.../random/100000/0/results.pkl`

## Usage

### Option 1: Submit SLURM Jobs (Recommended)

```bash
cd /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/semitones_marker_genes

# Run both SPS and Random analyses
bash run_marker_gene_analysis.sh both

# After both jobs complete, run comparison
bash run_marker_gene_analysis.sh compare
```

### Option 2: Run Interactively

```bash
# Activate environment
conda activate semitones_env

# Run SPS analysis
python semitones_marker_pipeline.py --method sps --ncpu 8 --skip-permutation

# Run Random analysis
python semitones_marker_pipeline.py --method random --ncpu 8 --skip-permutation

# Run comparison (activate facs_sampling for plotting)
conda activate facs_sampling
python compare_marker_genes.py
```

## Pipeline Details

### Step 1: Data Loading
- Loads full MCC dataset (~11GB)
- Loads sampling indices (SPS or Random, 100k cells, rep 0)

### Step 2: SEMITONES Enrichment Scoring
- Uses 100k sampled cells as **reference cells**
- Computes pairwise similarities (all 3.16M cells → 100k references)
- For each gene, computes enrichment score via linear regression
- Result: enrichment score matrix (genes × reference cells)

### Step 3: Marker Gene Identification
- For each cell type, selects reference cells of that type
- Aggregates enrichment scores across those reference cells
- Ranks genes by mean enrichment score
- Reports top 50 marker genes per cell type

### Step 4: Literature Marker Evaluation
- Compares identified markers against literature-based markers
- Computes:
  - **Recall@10**: % of literature markers in top 10 genes
  - **Recall@50**: % of literature markers in top 50 genes
  - **Mean Rank**: Average rank of literature markers

### Step 5: SPS vs Random Comparison
- Computes improvement metrics for SPS over Random
- Calculates marker gene overlap (Jaccard similarity)
- Generates visualizations

## Literature Marker Sources

Marker genes are sourced from:

1. **CellMarker 2.0** (Hu et al., 2023)
   - DOI: 10.1093/nar/gkac947

2. **PanglaoDB** (Franzén et al., 2019)
   - DOI: 10.1093/database/baz046

3. **Lateral Mesoderm** (Loh et al., 2016)
   - DOI: 10.1016/j.cell.2016.06.011

See `MARKER_GENE_REFERENCES.md` for full citations and marker gene lists.

## Expected Results

**Hypothesis**: SPS sampling should show:
- Higher recall of literature markers
- Better (lower) mean ranks for literature markers
- Especially strong improvement for rare cell types (osteoblast)

**Key Metrics**:
- Recall@10, Recall@50 per cell type
- Mean rank of literature markers
- Jaccard similarity of top marker genes

## Output Files

| File | Description |
|------|-------------|
| `*_marker_genes.pkl` | Full marker gene results (pickle) |
| `*_marker_genes_by_celltype.csv` | Top 50 marker genes per cell type |
| `*_literature_evaluation.csv` | Recall and rank metrics vs literature |
| `marker_gene_comparison.csv` | SPS vs Random comparison table |
| `figures/*.png` | Visualization plots |

## SLURM Configuration

- **Memory**: 128GB (full MCC dataset is ~11GB, plus SEMITONES computations)
- **CPUs**: 8 (for parallel enrichment scoring)
- **Time**: 24 hours
- **Partition**: normal
- **Account**: ohler

## Troubleshooting

### SEMITONES not found
```bash
conda activate semitones_env
pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip
```

### Memory issues
- Increase SLURM memory request
- Use `--skip-permutation` flag to skip permutation testing

### Missing results for comparison
Ensure both SPS and Random analyses completed successfully before running comparison.

## References

- SEMITONES: Vlot et al. (2022). Nucleic Acids Research. DOI: 10.1093/nar/gkac605
- CellMarker 2.0: Hu et al. (2023). Nucleic Acids Research. DOI: 10.1093/nar/gkac947
- PanglaoDB: Franzén et al. (2019). Database. DOI: 10.1093/database/baz046
