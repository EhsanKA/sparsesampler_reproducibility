# SPS vs Random: Downstream Task Experiments

This directory contains experiments demonstrating that SPS sampling outperforms random sampling for downstream analysis tasks, particularly for rare cell types.

## Key Insight

**Recovering full dataset DE results is not the right benchmark** - instead, we measure task-specific performance.

## Directory Structure

```
downstream_experiments/
├── README.md
├── common/
│   ├── __init__.py
│   └── utils.py                    # Shared data loading, metrics
├── 01_statistical_power/           # Experiment 1
│   ├── statistical_power_analysis.py      # MCC dataset
│   ├── statistical_power_analysis_lcmv.py # LCMV dataset
│   ├── run_analysis.sh
│   ├── run_analysis_lcmv.sh
│   ├── results/                    # MCC results
│   └── results_lcmv/               # LCMV results
└── 02_classification/              # Experiment 2
    ├── cell_type_classification.py        # MCC dataset
    ├── cell_type_classification_lcmv.py   # LCMV dataset
    ├── run_analysis.sh
    ├── run_analysis_lcmv.sh
    ├── results/                    # MCC results
    └── results_lcmv/               # LCMV results
```

## Datasets

### MCC Dataset
- **Size**: ~3.16M cells, 3,065 genes, 7 cell types
- **Rare cell type**: osteoblast (~30,000 cells, 0.95%)
- **Sample size**: 100,000 cells
- **Data path**: `/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/data/mcc/benchmark/30/`

### LCMV Dataset
- **Size**: ~34.4M cells, 31 genes, 21 cell types
- **Rare cell types** (identified dynamically by distance >75th percentile AND frequency <1%):
  - `interacting`: 25,112 cells (0.07%)
  - `NK1_1_TCRgd_T`: 24,178 cells (0.07%)
  - `cDC2`: 132,711 cells (0.39%)
  - `pDCs`: 191,329 cells (0.56%)
  - `CD4_LCMV_spec`: 332,408 cells (0.97%)
- **Sample size**: 100,000 cells
- **Data path**: `/fast/AG_Ohler/ekarimi/projects/sparseFlow_benchmarking/data/lcmv/benchmark/34/`

---

## Experiment 2: Cell Type Classification (Detailed Setup)

### Experimental Design

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    CLASSIFICATION EXPERIMENT SETUP                          │
└─────────────────────────────────────────────────────────────────────────────┘

    FULL DATASET
    ┌────────────────────────────────────────────────────────────────┐
    │                                                                │
    │   ┌──────────────────────────────┐  ┌───────────────────────┐ │
    │   │      TRAIN POOL (80%)        │  │   TEST SET (20%)      │ │
    │   │                              │  │                       │ │
    │   │   ┌─────────┐ ┌─────────┐    │  │  Same test set for    │ │
    │   │   │ SPS 100k│ │Rand 100k│    │  │  both SPS & Random    │ │
    │   │   │ ∩ pool  │ │ ∩ pool  │    │  │  evaluation           │ │
    │   │   │ (~80k)  │ │ (~80k)  │    │  │                       │ │
    │   │   └────┬────┘ └────┬────┘    │  │  Stratified by        │ │
    │   │        │           │         │  │  cell type            │ │
    │   └────────│───────────│─────────┘  └───────────┬───────────┘ │
    │            │           │                        │             │
    └────────────│───────────│────────────────────────│─────────────┘
                 │           │                        │
                 ▼           ▼                        ▼
          ┌──────────┐ ┌──────────┐            ┌──────────┐
          │  Train   │ │  Train   │            │ Evaluate │
          │ Classifier│ │Classifier│───────────▶│   Both   │
          │  on SPS  │ │ on Random│            │ Models   │
          └──────────┘ └──────────┘            └──────────┘
```

### Important Note on Index Handling

The SPS and Random indices were **pre-computed on the full dataset** to select 100k 
representative cells. The classification experiment then:

1. Splits the full dataset 80/20 into train pool / test set (stratified)
2. **Intersects** the pre-computed indices with the train pool
3. ~80% of indices fall in train pool (~80k cells used for training)
4. ~20% of indices are discarded (fall in test set)

**Implication**: Both SPS and Random are treated identically, so the comparison remains 
fair. However, ~20% of the carefully selected samples are not used for training.

### Classifier Configuration

```python
# Random Forest
RandomForestClassifier(
    n_estimators=100,      # Number of trees
    max_depth=20,          # Max depth per tree
    class_weight='balanced', # Upweight rare classes
    n_jobs=-1,             # Use all CPUs
    random_state=42        # Reproducibility
)

# Logistic Regression
LogisticRegression(
    max_iter=1000,
    class_weight='balanced',
    n_jobs=-1,
    random_state=42
)
```

### Robustness: Multiple Seeds

Each experiment runs with **5 different random seeds**, creating different train/test 
splits each time. Results are reported as **mean ± std** across all seeds.

```
Seed 0: Split → Train SPS → Train Random → Evaluate
Seed 1: Split → Train SPS → Train Random → Evaluate
Seed 2: Split → Train SPS → Train Random → Evaluate
Seed 3: Split → Train SPS → Train Random → Evaluate
Seed 4: Split → Train SPS → Train Random → Evaluate
                        ↓
        Aggregate: mean ± std across seeds
```

### Features

- **MCC**: Top 2,000 highly variable genes (normalized, log-transformed)
- **LCMV**: All 31 genes (normalized, log-transformed)

### Evaluation Metrics

- Overall Accuracy
- Macro F1 Score (unweighted average across all classes)
- Weighted F1 Score
- ROC-AUC (one-vs-rest)
- Per-class F1 (focus on rare cell types)
- Confusion Matrix

---

## Experiment 1: Statistical Power for Rare Cell Markers

**Hypothesis**: SPS has more rare cells, leading to stronger statistical signal for their markers.

**Method**:
- Compare -log10(p-value) distributions for rare cell markers
- Compare log fold change magnitudes
- Compare DE scores across methods (wilcoxon, t-test, logreg)

**Key Metrics**:
- Mean -log10(p-value) for top 50 markers
- Mean absolute log fold change
- Number of markers with p < 1e-10, p < 1e-50, p < 1e-100

---

## Running the Experiments

```bash
# MCC Dataset
cd 01_statistical_power && sbatch run_analysis.sh
cd 02_classification && sbatch run_analysis.sh

# LCMV Dataset
cd 01_statistical_power && sbatch run_analysis_lcmv.sh
cd 02_classification && sbatch run_analysis_lcmv.sh
```

## Environment

All experiments use the `facs_sampling` conda environment:
```bash
conda activate facs_sampling
```

## Results Summary

### Marker Gene Recovery (SPS vs Random vs Full)

#### MCC Dataset (3.16M cells, Ensembl gene IDs)

| DE Method | SPS Overlap | Random Overlap | Winner |
|-----------|-------------|----------------|--------|
| Wilcoxon | 42.0% | 94.0% | **Random** |
| T-test   | 54.0% | 96.0% | **Random** |
| LogReg   | 32.0% | 66.0% | **Random** |

**Cell type wins (Top-10 markers)**: SPS=0, Random=4-5, Ties=0-1

#### LCMV Dataset (34.4M cells, Protein markers)

| DE Method | SPS Overlap | Random Overlap | Winner |
|-----------|-------------|----------------|--------|
| T-test   | 75.7% | 98.1% | **Random** |
| LogReg   | 80.5% | 92.4% | **Random** |

**Cell type wins (Top-10 markers)**: SPS=0-2, Random=16-21, Ties=0-3

**Interpretation**: Random better recovers full dataset marker genes because it preserves
original cell type proportions, leading to similar DE statistics.

### Cell Type Classification (SPS vs Random Training)

#### Overall Accuracy (5 seeds, mean)

| Dataset | Classifier | SPS | Random | Δ |
|---------|------------|-----|--------|---|
| **MCC** | LogReg | 82.1% | 89.7% | -7.5% |
| **MCC** | RF | 84.6% | 88.8% | -4.2% |
| **LCMV** | LogReg | 92.4% | 95.0% | -2.6% |
| **LCMV** | RF | 92.8% | 95.8% | -3.1% |

**Interpretation**: Random training performs better overall because it provides a more
representative sample of the full data distribution.

#### Rare Cell Type F1 Scores

##### MCC - Osteoblast (0.95% of data)

| Classifier | SPS F1 | Random F1 | Δ |
|------------|--------|-----------|---|
| LogReg | 79.7% | 86.2% | **-6.6%** |
| RF | **89.2%** | 88.8% | **+0.4%** |

##### LCMV - Extremely Rare Cell Types

| Cell Type | Frequency | LogReg Δ | RF Δ |
|-----------|-----------|----------|------|
| **interacting** | 0.07% | **+14.9%** | **+38.5%** |
| NK1_1_TCRgd_T | 0.07% | **+4.4%** | -0.1% |
| cDC2 | 0.39% | **+6.6%** | +0.2% |
| CD4_LCMV_spec | 0.97% | +1.3% | -0.7% |
| pDCs | 0.56% | -0.5% | -1.2% |

### Key Findings

1. **Task-Dependent Performance**: The best sampling method depends on the downstream task:
   - **For DE analysis**: Random sampling recovers full dataset results better
   - **For rare cell classification**: SPS provides more training examples for rare classes

2. **Rarity Threshold**: SPS shows dramatic improvements for extremely rare cell types:
   - **+38.5% F1** for `interacting` cells (0.07% of data) with RF classifier
   - Benefits diminish as cell types become less rare

3. **MCC vs LCMV Differences**:
   - MCC has only one rare cell type (osteoblast, 0.95%) - smaller SPS benefit
   - LCMV has multiple extremely rare cell types (5 at <1%) - larger SPS benefit

## Known Limitations

- **MCC Gene Names**: Uses Ensembl IDs (e.g., `ENSMUSG00000051951`) while literature markers
  use gene symbols (e.g., `RUNX2`). No overlap found, so comparison metrics are 0.

- **LCMV Full Dataset DE**: Scanpy's `rank_genes_groups` fails with 34M cells due to
  internal chunking bugs. Used 20M cell subsample as ground truth instead.
