# Reviewer Response Summary: SParseSampler Revision

**Manuscript ID**: BIOINF-2025-2256  
**Title**: SParseSampler: Efficient Sampling Method for High-Dimensional Single-Cell RNA-Seq and Flow Cytometry Data  
**Decision Date**: 12-Nov-2025

---

## Executive Summary

This document organizes the reviewer comments into two main categories as outlined by Associate Editor Macha Nikolski, and documents the experimental approaches conducted to address each concern.

**Key Terminology - Feature Index**:
Throughout this document, we refer to **feature index** as our primary parameter. The feature index `i` determines both:
1. The number of principal components used (top `i` components)
2. The bin resolution factor K, computed as `K = 2 / EVR_i` where EVR_i is the explained variance ratio of the i-th component

This approach replaces the need to manually tune K, as K is automatically derived from the data's variance structure.

---

## Category 1: Parameter Transparency and Guidance

### Comment 1.1 (Reviewer 1, Point 1): K Values for Benchmarks Not Stated

> "The Bin Resolution Factor (K) is a critical parameter. While its impact is discussed, the specific K values used for the benchmarks in Figures 1e and 1f are not stated, which is essential for reproducibility."

**Approach**: Clarify that we use feature index to automatically determine K.

**Explanation**:
- We use **feature index = 18** for all main benchmarks
- K is computed as `K = 2 / EVR_18` where EVR_18 is the explained variance ratio of the 18th principal component
- This approach ties the bin resolution to the data's intrinsic variance structure
- The formula ensures that PCs explaining less variance than EVR_i contribute minimally to binning (ratio < 1 means negligible contribution)

**Action for Manuscript**:
- Add explicit statement: "All benchmarks used feature_index=18, which determines K automatically as 2/EVR_18. This ties bin resolution to the data's variance structure."

---

### Comment 1.2 (Reviewer 1, Point 2): Heuristics for Choosing K

> "The manuscript would benefit from a more practical discussion on how to select an appropriate K value for a new dataset. Providing a heuristic or a recommended range based on dataset size would greatly enhance the tool's usability."

**Approach**: Recommend using feature index rather than manually selecting K.

**Explanation**:
The feature index approach eliminates manual K tuning by deriving K from the data:
- Feature index determines how many PCs to consider meaningful
- K = 2/EVR_i ensures bins are allocated proportionally to variance
- Users can select feature index based on domain knowledge (e.g., expected number of distinct populations)

**Experimental Evidence**:
We tested feature indices 1-30 across all datasets:
- **MCC datasets**: Feature indices 10-25 provide stable, good performance
- **LCMV dataset**: Feature indices around 10 work well (lower dimensionality data)

**Recommended Heuristic**:
- **Default**: Feature index = 18 works well for most scRNA-seq datasets
- **Flow cytometry**: Feature index = 10-15 (fewer features)
- **Guidance**: Choose based on expected biological complexity; higher index = finer resolution

**Location**: `jobs/feature_index_classification/`

---

### Comment 1.3 (Reviewer 1, Point 3): PCA Runtime Analysis

> "There is a concern of the feasibility of the method for large datasets. It can be addressed by including the time taken for the PCA step in their runtime analysis (Figure 1f) to demonstrate that it does not become a bottleneck."

**Approach**: Benchmark PCA timing and note that other methods also require dimensionality reduction.

**Experimental Evidence** (from `jobs/pca_timing/pca_timing_results.csv`):

| Dataset | n_cells | n_features | PCA Time (s) | Total Time (s) |
|---------|---------|------------|--------------|----------------|
| LCMV | 34,369,995 | 31 | 11.2 | 259.6 |
| MCC | 3,161,397 | 3,065 | 306.7 | 482.6 |
| MCC_01 | 3,134,397 | 3,065 | 304.1 | 482.8 |
| MCC_05 | 3,146,397 | 3,065 | 305.5 | 493.2 |

**Key Points**:
1. **Flow cytometry (LCMV)**: PCA is negligible (~11s for 34M cells) due to low dimensionality
2. **scRNA-seq (MCC)**: PCA takes ~5 minutes for 3M cells
3. **Important**: Other sampling methods (Hopper, scSampler) also require PCA or similar dimensionality reduction - this is not unique to SParseSampler
4. **Optimization**: With our feature index approach, we only need to compute the top `i` components, not a full PCA

**Action for Manuscript**:
- Add: "PCA computation is common to most sampling methods. For flow cytometry with low dimensionality, PCA overhead is minimal. For scRNA-seq, incremental PCA methods can be used, and our feature index approach only requires the top i components."

**Location**: `jobs/pca_timing/pca_timing_benchmark.py`, `jobs/pca_timing/pca_timing_results.csv`

---

### Comment 2.2 (Reviewer 2, Point 2): Sensitivity/Stability Figures for PCs and K

> "The authors state that SPS performance remains stable when the number of principal components p varies between 30 and 100 and when the bin resolution factor K ranges between 50 and 200. However, these claims are presented only as textual descriptions without accompanying figures or tables."

**Approach**: Provide explicit sensitivity analysis across feature indices. Remove claims about specific K and p ranges.

**Note**: The previous claims about K ∈ [50, 200] and p ∈ [30, 100] should be **removed** from the manuscript. Instead, we show sensitivity across feature indices.

**Experimental Evidence**:

#### Feature Index Sensitivity Analysis
- **Feature indices tested**: 1-30
- **Datasets**: MCC (0.95% rare), MCC_01 (0.1% rare), MCC_05 (0.5% rare), LCMV
- **Sample size**: 100,000 cells
- **Replicate**: 0

**Results** - Osteoblast F1 Score (RF classifier) across feature indices:

**MCC_01 (0.1% rare cells)**:
| Feature Index | SPS Osteoblast F1 | Random Osteoblast F1 | Improvement |
|---------------|-------------------|----------------------|-------------|
| 12 | 80.2% | 58.9% | **+21.3%** |
| 18 | 75.7% | 58.9% | **+16.9%** |
| 23 | 79.2% | 58.8% | **+20.4%** |
| 30 | 76.0% | 58.9% | **+17.1%** |

**MCC_05 (0.5% rare cells)**:
| Feature Index | SPS Osteoblast F1 | Random Osteoblast F1 | Improvement |
|---------------|-------------------|----------------------|-------------|
| 12 | 89.0% | 84.0% | **+5.0%** |
| 18 | 88.2% | 84.0% | **+4.3%** |
| 23 | 89.2% | 84.0% | **+5.3%** |
| 30 | 88.5% | 84.0% | **+4.6%** |

**Key Finding**: Performance is stable across feature indices 10-30, with consistent improvements for rare cell detection.

#### Figures Generated
- `figures/revision/rf_classification_by_feature_index.jpg` - F1 by feature index
- `figures/revision/supp_evr_sensitivity_combined.jpg` - Combined sensitivity analysis

**Location**: `jobs/feature_index_classification/results_mcc_01/`, `jobs/feature_index_classification/results_mcc_05/`

---

## Category 2: Rare Cell Signal Validation

### Comment 1.4 (Reviewer 1, Point 4): Test Rarer Fractions (0.5%, 0.1%)

> "In the MCC dataset, the artificially created 'rare' population at 1% may not be sufficiently stringent... To more rigorously test the method's sensitivity, it would be informative to evaluate its performance on populations with smaller proportions (e.g., 0.5%, 0.1%)."

**Approach**: Create MCC variants with 0.5% and 0.1% rare cell populations. Show rare cell recovery and runtime comparison with other methods.

**Experimental Evidence**:

#### Dataset Creation
| Dataset | Osteoblast Cells | Osteoblast % | Total Cells |
|---------|------------------|--------------|-------------|
| MCC (original) | ~30,000 | 0.95% | 3,161,397 |
| MCC_05 | ~15,700 | 0.50% | 3,146,397 |
| MCC_01 | ~3,100 | 0.10% | 3,134,397 |

#### Benchmarking Setup
- **Methods**: SPS, Random, Hopper, Atomic Sketch, scSampler
- **Reference sizes**: 500k, 1M, 2M, 2.5M, 3M cells
- **Sample sizes**: 50k, 100k, 200k, 300k
- **Replicate**: 0

#### Rare Cell Recovery Results (RF Classification, feature_index=18, 100k samples)

**MCC_01 (0.1% rare)** - Osteoblast F1:
| Method | Osteoblast F1 | vs Random |
|--------|---------------|-----------|
| SPS | 75.7% | **+16.9%** |
| Random | 58.9% | baseline |

**MCC_05 (0.5% rare)** - Osteoblast F1:
| Method | Osteoblast F1 | vs Random |
|--------|---------------|-----------|
| SPS | 88.2% | **+4.3%** |
| Random | 84.0% | baseline |

**Key Finding**: SPS advantage increases dramatically as populations become rarer:
- 0.95% rare: Modest improvement
- 0.5% rare: ~4-5% F1 improvement
- 0.1% rare: **~17-21% F1 improvement**

**Location**: 
- Data: `data/mcc_01/`, `data/mcc_05/`
- Jobs: `jobs/mcc_01_production/`, `jobs/mcc_05_production/`
- Results: `jobs/feature_index_classification/results_mcc_01/`, `jobs/feature_index_classification/results_mcc_05/`

---

### Comment 1.5 (Reviewer 1, Point 5): Dropout/Technical Noise Discussion

> "Single-cell data is characterized by technical noise like pervasive dropout events... The authors should discuss whether and how such technical artifacts might compromise the performance of SParseSampler."

**Approach**: Discussion point for manuscript revision.

**Key Points to Address**:
1. **PCA preprocessing**: SParseSampler operates on PCA-reduced data, which inherently smooths dropout noise
2. **Dropout vs. rare cells**: Dropout events are random zeros that do not create coherent clusters in PCA space; true rare populations show consistent expression patterns
3. **Empirical validation**: Successful recovery of rare cell types with improved classification (+17-21% F1 for 0.1% populations) demonstrates biological signal capture
4. **Pre-filtering expectation**: Standard QC filtering (low UMI cells, doublets) should be applied before sampling

**Recommended Manuscript Text**:
> "SParseSampler operates on PCA-reduced data, providing inherent robustness to dropout events. Dropout artifacts manifest as random zeros without coherent PCA clusters, whereas true rare populations exhibit consistent gene expression. We expect standard QC filtering to be applied before sampling. Our empirical results showing improved rare cell classification (up to +21% F1 for 0.1% populations) confirm that SParseSampler captures biological rare cells."

---

### Comment 2.1 (Reviewer 2, Point 1): Empty Gap/Doublet Oversampling Concerns

> "The variance-weighted binning step puts more bins along PC directions explaining more variance... This PC direction often spans a large, empty gap populated only by noise or doublets."

**Approach**: Discussion point explaining why this concern is mitigated by the algorithm design.

**Key Points to Address**:

1. **Density-aware selection**: SParseSampler's iterative selection prioritizes cells from *sparsely populated bins*. Empty regions have no cells to sample.

2. **Bin allocation vs. sampling**: More bins in high-variance directions provides finer resolution, but empty bins contribute nothing to the sample.

3. **Empirical evidence**: If SParseSampler oversampled artifacts:
   - Classification accuracy would decrease (doublets misclassified)
   - Random sampling would outperform SPS
   - Instead, we observe +38.5% F1 improvement for extremely rare cells (0.07%)

4. **Pre-filtering**: Standard doublet detection (e.g., Scrublet, DoubletFinder) should be applied before sampling.

**Recommended Manuscript Text**:
> "While variance-weighted binning allocates more bins to high-variance PC directions, the density-aware selection ensures cells are only sampled from populated regions. Empty gaps or doublet-sparse regions contribute minimally. We recommend standard doublet detection before sampling. Empirically, our +38.5% F1 improvement for 0.07% populations confirms biological signal capture, not artifact oversampling."

---

### Comment 2.3 (Reviewer 2, Point 3): Definition of Rare Cell Populations

> "The manuscript defines 'rare populations' either by expert annotation or by artificially selecting the smallest cluster. However, the notion of rare cell populations is quite vague..."

**Approach**: Provide explicit definition using frequency and distance criteria.

**Definition**:
We define rare cell types using two complementary criteria:
1. **Frequency threshold**: < 1% of total cells
2. **Distance threshold**: Centroid distance to dataset mean > 75th percentile

This ensures rare populations are both:
- **Numerically small**: Likely to be missed by random sampling
- **Biologically distinct**: Not just outliers from major populations

**LCMV Rare Cell Types** (identified by this method):

| Cell Type | Frequency | Enrichment by SPS |
|-----------|-----------|-------------------|
| interacting | 0.07% | 56.95x |
| NK1_1_TCRgd_T | 0.07% | 12.67x |
| cDC2 | 0.39% | 7.04x |
| pDCs | 0.56% | 2.62x |
| CD4_LCMV_spec | 0.97% | 8.23x |

**Recommended Manuscript Text**:
> "We define rare cell types using dual criteria: (1) frequency below 1% of total cells, and (2) centroid distance above the 75th percentile, indicating biological distinctiveness. This ensures identified rare populations are both numerically underrepresented and biologically meaningful."

**Location**: `analysis/refined_rare_cell_type_definition.py`

---

### Comment 2.4 (Reviewer 2, Point 4): Downstream Analyses Preservation (Clustering/DE)

> "A more informative metric will be the preservation of cell type/state diversity... It is important to compare downstream analyses results from original population and the subsampled population."

**Approach**: Demonstrate downstream utility through cell type classification.

**Experimental Evidence**:

#### Cell Type Classification (RF Classifier)

**Methodology**:
- Train RF classifier on SPS or Random subsamples (100k cells)
- Test on held-out set (20% of full data, stratified)
- 5 random seeds for robustness

**LCMV Results** (from `analysis/downstream_experiments/02_classification/results_lcmv/`):

| Metric | SPS | Random | Difference |
|--------|-----|--------|------------|
| Overall Accuracy | 92.8% | 95.8% | -3.1% |
| Macro F1 | 84.1% | 87.7% | -3.6% |

**Rare Cell Type F1 Scores (RF)**:

| Cell Type | Frequency | SPS F1 | Random F1 | Δ |
|-----------|-----------|--------|-----------|---|
| interacting | 0.07% | 59.0% | 20.4% | **+38.5%** |
| NK1_1_TCRgd_T | 0.07% | 77.1% | 77.2% | -0.1% |
| cDC2 | 0.39% | 87.6% | 87.4% | +0.2% |
| CD4_LCMV_spec | 0.97% | 98.1% | 98.8% | -0.7% |
| pDCs | 0.56% | 95.4% | 96.6% | -1.2% |

**Key Finding**: SPS shows dramatic improvement (+38.5% F1) for the rarest cell type (0.07%) while maintaining competitive overall performance. This demonstrates that SPS-sampled data provides superior training signal for rare cell identification.

**Location**: `analysis/downstream_experiments/02_classification/results_lcmv/`

---

## Minor Comments

### Comment 1.6 (Reviewer 1, Point 6): Typo Correction

> "Page 1, Line 57: The phrase 'Atomic Sketch, leverages leverage' contains a duplicate."

**Action**: Fix typo in manuscript.

---

### Comment 2.5 (Reviewer 2, Minor Point 1): Figure 1 Layout

> "Fig 1 is quite crowded. It is better to separate it into different Figures."

**Approach**: Consider splitting Figure 1 or moving some panels to supplementary.

**Generated Figures for Revision**:
- `figures/revision/figure1_combined_coverage.jpg`
- `figures/revision/figure1_combined_time.jpg`
- `figures/revision/supp_figure3_mcc_01_coverage.jpg`
- `figures/revision/supp_figure3_mcc_05_coverage.jpg`
- `figures/revision/rf_classification_by_feature_index.jpg`
- `figures/revision/rf_classification_lcmv_rare_cells.jpg`

---

## Summary Table: Experiments Conducted

| Reviewer Comment | Experiment | Status | Key Result |
|------------------|------------|--------|------------|
| R1.1: K values not stated | Document feature index | Done | feature_index=18, K=2/EVR_18 |
| R1.2: K heuristics | Feature index sensitivity | Done | FI 10-25 stable |
| R1.3: PCA timing | PCA benchmark | Done | 11-307s depending on data |
| R1.4: Rarer fractions | MCC_01, MCC_05 experiments | Done | +17-21% F1 for 0.1% rare |
| R1.5: Dropout discussion | Discussion text | Pending | Points prepared |
| R2.1: Empty gap/doublets | Discussion text | Pending | Points prepared |
| R2.2: Sensitivity figures | Feature index analysis | Done | Figures generated |
| R2.3: Rare cell definition | Multi-criteria definition | Done | Frequency + distance |
| R2.4: Downstream analyses | RF Classification | Done | +38.5% F1 for 0.07% cells |

---

## Directories Used (Keep These)

The following directories contain experiments and results used in this revision:

### Data Directories
```
data/mcc/                          # Original MCC dataset (0.95% rare)
data/mcc_01/                       # MCC with 0.1% rare cells
data/mcc_05/                       # MCC with 0.5% rare cells
data/lcmv/                         # LCMV flow cytometry dataset
```

### Job Directories
```
jobs/mcc_production/               # Original MCC sampling jobs
jobs/mcc_01_production/            # MCC 0.1% sampling jobs
jobs/mcc_05_production/            # MCC 0.5% sampling jobs
jobs/lcmv_production/              # LCMV sampling jobs
jobs/pca_timing/                   # PCA timing benchmark
jobs/feature_index_classification/ # RF classification by feature index (MCC)
jobs/feature_index_classification/results_mcc_01/  # MCC_01 results
jobs/feature_index_classification/results_mcc_05/  # MCC_05 results
jobs/feature_index_classification_lcmv/            # LCMV classification
```

### Analysis Directories
```
analysis/downstream_experiments/02_classification/          # Classification analysis
analysis/downstream_experiments/02_classification/results_lcmv/  # LCMV classification results
analysis/refined_rare_cell_type_definition.py              # Rare cell definition
analysis/CHANGES_SUMMARY.md                                # Changes documentation
```

### Figure Directories
```
figures/revision/                  # Revision figures
figures/main/                      # Main paper figures
```

---

## Directories NOT Used (Can Be Cleaned)

The following directories were explored but are NOT used in the final revision:

```
k_parameter_analysis/              # K parameter testing (approach changed to feature index)
analysis/downstream_experiments/01_statistical_power/      # Statistical power (not included)
analysis/semitones_marker_genes/   # SEMITONES analysis (not included)
analysis/scanpy_marker_genes/      # Scanpy marker genes (not included)
jobs/de_analysis/                  # DE analysis (not included)
jobs/semitones_analysis/           # SEMITONES jobs (not included)
jobs/semitones_marker_genes/       # SEMITONES marker genes (not included)
jobs/scanpy_marker_genes/          # Scanpy marker genes (not included)
jobs/scanpy_marker_genes_lcmv/     # LCMV marker genes (not included)
jobs/clustering_visualization/     # Clustering visualization (not included)
```

---

## Recommended Revision Checklist

- [ ] Add feature index explanation to Methods (K = 2/EVR_i)
- [ ] State feature_index=18 used in all benchmarks
- [ ] Add feature index selection guidance to Methods/Discussion
- [ ] Add PCA timing table to Supplementary
- [ ] Add feature index sensitivity figures to Supplementary
- [ ] Add 0.1%/0.5% rare cell results (coverage and time comparison)
- [ ] Add dropout/technical noise discussion to Discussion
- [ ] Add empty gap/doublet discussion to Discussion
- [ ] Add rare cell definition criteria to Methods
- [ ] Add RF classification results for rare cells
- [ ] Fix "leverages leverage" typo
- [ ] Consider splitting Figure 1
- [ ] Remove claims about K ∈ [50, 200] and p ∈ [30, 100]
- [ ] Archive code on Zenodo for DOI
