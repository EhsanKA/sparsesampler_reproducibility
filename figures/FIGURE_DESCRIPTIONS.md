# SparseSampler Figure Descriptions

This document provides detailed explanations of all figures generated for the SparseSampler reproducibility study.

---

## Table of Contents

1. [Main Figures](#main-figures)
   - [Figure 1: Combined Coverage](#figure-1-combined-coverage)
   - [Figure 1: Combined Time Performance](#figure-1-combined-time-performance)
   - [UMAP Visualizations](#umap-visualizations)
2. [Supplementary Figures](#supplementary-figures)
   - [Supplementary Figure 3: MCC Coverage](#supplementary-figure-3-mcc-coverage)
   - [Supplementary Figure 4: LCMV Coverage](#supplementary-figure-4-lcmv-coverage)
   - [Supplementary Figures 5-6: Time Performance](#supplementary-figures-5-6-time-performance)
3. [Revision Figures](#revision-figures)
   - [Revision Coverage Figures](#revision-coverage-figures)
   - [Revision Time Figures](#revision-time-figures)
   - [EVR Analysis Figures](#evr-analysis-figures)

---

## Main Figures

### Figure 1: Combined Coverage

**File:** `figures/figure1_combined_coverage.jpg`

**Description:**  
This figure compares the rare cell type coverage achieved by different sampling methods across both MCC and LCMV datasets. It shows the percentage of rare cells captured in sampled subsets relative to the total rare cells in the reference dataset.

**Datasets:**
- **MCC (Mouse Cell Compendium):** References of 0.5M, 1M, 2M, 2.5M, and 3.2M cells
- **LCMV (Lymphocytic choriomeningitis virus):** References of 1M, 5M, 10M, 20M, and 34M cells

**Sampling Methods Compared:**
- **Random:** Uniform random sampling (baseline)
- **SPS (SparseSampler):** Our proposed method using explained variance ratio (EVR) index 12
- **Hopper:** Diversity-based sampling using geometric sketching
- **Atomic:** Atomic sketching method
- **scSampler:** Single-cell specific sampling method

**Rare Cell Definition:**
- Cells with distance > 75th percentile from cluster centroid
- AND cell type frequency < 1% of total population

**Key Metrics:**
- Coverage percentage = (rare cells in sample / total rare cells in reference) × 100
- Higher coverage indicates better preservation of rare cell populations

---

### Figure 1: Combined Time Performance

**File:** `figures/figure1_combined_time.jpg`

**Description:**  
This figure compares the computational runtime of different sampling methods as a function of dataset size. It demonstrates the scalability advantage of SparseSampler over other methods.

**Axes:**
- **X-axis:** Reference dataset size (number of cells)
- **Y-axis:** Runtime in seconds (log scale)

**Sample Sizes Evaluated:**
- MCC: 50k, 100k, 200k, 300k cells
- LCMV: 50k, 100k, 200k cells

**Key Observations:**
- SPS shows near-constant or sub-linear scaling with dataset size
- Other methods (Hopper, Atomic, scSampler) show super-linear scaling
- The runtime advantage of SPS becomes more pronounced at larger dataset sizes

---

### UMAP Visualizations

**File:** `figures/main/umaps/combined_umaps.jpg`

**Description:**  
This figure shows UMAP (Uniform Manifold Approximation and Projection) visualizations comparing the spatial distribution of sampled cells across methods.

**Layout:**
- Two rows: MCC dataset (top) and LCMV dataset (bottom)
- Three columns per row showing:
  1. **SPS samples:** Cells selected by SparseSampler (red) vs. unselected (gray)
  2. **Rare cells:** Ground truth rare cell types (blue) vs. common cells (gray)
  3. **Random samples:** Cells selected by random sampling (red) vs. unselected (gray)

**Purpose:**
- Visually demonstrates that SPS preferentially samples from regions containing rare cell types
- Shows that random sampling misses sparse regions of the manifold
- Illustrates the biological relevance of SPS's selection strategy

---

## Supplementary Figures

### Supplementary Figure 3: MCC Coverage

**File:** `figures/supp_figure3_mcc_coverage.jpg`

**Description:**  
Detailed breakdown of rare cell coverage for the MCC dataset across all sample sizes and reference sizes.

**Layout:** Multi-panel grid showing:
- Rows: Different sample sizes (50k, 100k, 200k, 300k)
- Columns: Different reference sizes (0.5M to 3.2M cells)

**Content per Panel:**
- Bar chart comparing coverage percentage across all sampling methods
- Error bars showing variability across replicates (if applicable)

**Rare Cell Type Focus:** Osteoblasts (frequency < 1%)

---

### Supplementary Figure 4: LCMV Coverage

**File:** `figures/supp_figure4_lcmv_coverage.jpg`

**Description:**  
Detailed breakdown of rare cell coverage for the LCMV dataset across all sample sizes and reference sizes.

**Layout:** Multi-panel grid showing:
- Rows: Different sample sizes (50k, 100k, 200k)
- Columns: Different reference sizes (1M to 34M cells)

**Rare Cell Types in LCMV:**
- Interacting cells
- NK1.1+ TCRγδ T cells
- cDC2 (conventional dendritic cells type 2)
- pDCs (plasmacytoid dendritic cells)
- CD4+ LCMV-specific T cells

---

### Supplementary Figures 5-6: Time Performance

**Files:**
- `figures/supp_figure5-6_mcc_time.jpg`
- `figures/supp_figure5-6_lcmv_time.jpg`

**Description:**  
Detailed runtime analysis for MCC and LCMV datasets separately, showing how computation time scales with both sample size and reference size.

**Content:**
- Line plots with reference size on x-axis and runtime on y-axis
- Separate lines for each sampling method
- Multiple panels for different sample sizes

---

## Revision Figures

### Revision Coverage Figures

**Files:**
- `figures/revision/figure1_combined_coverage.jpg`
- `figures/revision/supp_figure3_mcc_01_coverage.jpg`
- `figures/revision/supp_figure3_mcc_05_coverage.jpg`

**Description:**  
These figures extend the coverage analysis to different rare cell frequency thresholds, addressing reviewer concerns about the generalizability of our rare cell definition.

**Frequency Thresholds:**
- **MCC_01:** Rare cells defined as frequency < 0.1% (more stringent)
- **MCC_05:** Rare cells defined as frequency < 0.5% (intermediate)

**Purpose:**
- Demonstrates that SPS maintains superior coverage across different rarity thresholds
- Shows robustness of the method to different definitions of "rare"

---

### Revision Time Figures

**Files:**
- `figures/revision/figure1_combined_time.jpg`
- `figures/revision/supp_figure5-6_mcc_01_time.jpg`
- `figures/revision/supp_figure5-6_mcc_05_time.jpg`

**Description:**  
Time performance analysis for the MCC_01 and MCC_05 datasets, confirming that runtime characteristics are consistent across different rare cell definitions.

---

### EVR Analysis Figures

**Files:**
- `figures/revision/combined_evr_rf_rank_figure_50k.jpg`
- `figures/revision/combined_evr_rf_rank_figure_100k.jpg`
- `figures/revision/combined_evr_rf_rank_figure_200k.jpg`

**Description:**  
Comprehensive analysis of how the Explained Variance Ratio (EVR) index affects SparseSampler performance. These figures justify the selection of EVR index 12 as the default parameter.

**Layout:** 3 rows × 4 columns
- **Columns:** MCC (1%), MCC (0.5%), MCC (0.1%), LCMV
- **Rows:**
  1. Coverage % by EVR Index
  2. RF Classification ΔF1 (SPS - Random) by EVR Index
  3. SPS Rank among sampling methods by EVR Index

**Row 1: Coverage Analysis**
- Shows how rare cell coverage varies with EVR index (1-30)
- Different lines represent different reference sizes
- Red dashed vertical line marks EVR index 12 (selected default)

**Row 2: Classification Performance**
- Shows improvement in Random Forest classification F1 score when using SPS vs. random sampling
- Positive values indicate SPS produces better training data
- Blue line: Macro F1 improvement
- Red line: Rare cell type F1 improvement (e.g., Osteoblast for MCC, multiple types for LCMV)

**Row 3: Rank Analysis**
- Shows SPS's rank (1-5) compared to other methods at each EVR index
- Rank 1 = best (most rare cells captured)
- Green shading: EVR indices where SPS is #1
- Red shading: EVR indices where SPS is not #1

**Key Finding:**
EVR index 12 provides a good balance across all metrics:
- High rare cell coverage across datasets
- Consistent classification improvement
- Competitive or top ranking among methods

---

## Technical Notes

### Data Sources
- **MCC:** Mouse Cell Compendium from the Tabula Muris Senis project
- **LCMV:** Single-cell data from LCMV infection study

### Sampling Configuration
- **EVR Index:** 12 (default for SPS)
- **Replicates:** Results averaged over multiple random seeds where applicable

### Figure Generation
All figures are generated using:
- Python with matplotlib for plotting
- Scanpy for single-cell data processing
- Scripts located in `figures/main/` and `figures/revision/`

### File Formats
- `.jpg` files at 350 DPI for manuscript inclusion
- `.pdf` files available for vector graphics (revision figures)

---

## Citation

If you use these figures or the SparseSampler method, please cite:

```
[Citation information to be added]
```
