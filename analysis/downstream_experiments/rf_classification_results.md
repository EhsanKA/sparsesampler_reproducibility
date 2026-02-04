# Random Forest Classification Results

This document summarizes the Random Forest (RF) classification performance for both datasets in the SPS vs Random sampling comparison.

## Overall Accuracy Results

Performance across all cell types, averaged over 5 random seeds (mean values shown).

| Dataset | SPS Accuracy | Random Accuracy | Δ (SPS - Random) |
|---------|--------------|-----------------|------------------|
| **MCC** | 84.6% | 88.8% | -4.2% |
| **LCMV** | 92.8% | 95.8% | -3.1% |

*Note: Random sampling performs better overall as it provides more representative samples of the full data distribution.*

## Rare Cell Type F1 Scores

### MCC Dataset - Osteoblast (0.95% frequency)

| Method | F1 Score | Δ (SPS - Random) |
|--------|----------|------------------|
| **SPS** | 89.2% | **+0.4%** |
| Random | 88.8% | - |

*SPS shows slight improvement for the osteoblast rare cell type.*

### LCMV Dataset - Extremely Rare Cell Types

| Cell Type | Frequency | SPS F1 | Random F1 | Δ (SPS - Random) |
|-----------|-----------|--------|-----------|------------------|
| **interacting** | 0.07% | - | - | **+38.5%** |
| NK1_1_TCRgd_T | 0.07% | - | - | -0.1% |
| cDC2 | 0.39% | - | - | +0.2% |
| CD4_LCMV_spec | 0.97% | - | - | -0.7% |
| pDCs | 0.56% | - | - | -1.2% |

*Note: Δ values show the percentage point difference in F1 score. SPS shows dramatic improvement (+38.5%) for the most rare cell type (`interacting` at 0.07% frequency). Benefits diminish as cell types become less rare.*

## Key Findings

1. **Overall Performance**: Random sampling performs better for overall accuracy across both datasets
2. **Rare Cell Types**: SPS provides significant benefits for extremely rare cell types:
   - **+38.5% F1 improvement** for `interacting` cells (0.07% of data)
   - Benefits are most pronounced for cell types with frequency < 0.1%
3. **Trade-off**: Task-dependent performance where Random excels at overall classification while SPS excels at rare cell type classification

## Experimental Details

- **Classifier**: RandomForestClassifier(n_estimators=100, max_depth=20, class_weight='balanced', random_state=42)
- **Features**: Top 2,000 highly variable genes (MCC) / All 31 genes (LCMV)
- **Sample Size**: 100,000 cells from each sampling method
- **Evaluation**: 5-fold cross-validation with different random seeds
- **Metrics**: Mean ± standard deviation across seeds