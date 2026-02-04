# How to Interpret SEMITONES Results

## Overview

SEMITONES evaluates how well your **100k SPS-sampled cells** (reference cells) represent the full dataset. The enrichment scores tell you how similar each query cell is to the reference set.

## Understanding Enrichment Scores

### What is an Enrichment Score?

- **Enrichment Score** = Mean cosine similarity between a query cell and all reference cells
- **Range**: 0 to 1 (higher = better)
- **Interpretation**:
  - **High score (e.g., >0.8)**: Query cell is very similar to reference cells → Well represented ✅
  - **Medium score (e.g., 0.5-0.8)**: Query cell is moderately similar → Reasonably represented
  - **Low score (e.g., <0.5)**: Query cell is different from references → Poorly represented ⚠️

### What This Means for Your Sampling

- **High enrichment across most cells** → Your SPS samples are good references (captured diversity well)
- **Low enrichment for many cells** → Your SPS samples might be missing some cell types
- **Variable enrichment by cell type** → Some types are better represented than others

## Interpreting the Plots

### Plot 1: Distribution of Enrichment Scores (Top Left)

**What it shows**: Histogram of all enrichment scores across all query cells

**How to interpret**:
- **Right-shifted distribution** (most scores >0.7) → Good! Most cells are well represented
- **Left-shifted distribution** (most scores <0.5) → Problem! Many cells poorly represented
- **Bimodal distribution** → Some cell types well represented, others not
- **Mean vs Median**: 
  - If mean > median → Some very high scores (good)
  - If mean < median → Some very low scores (concerning)

**Example from your results**:
- **mcc_01**: Mean ~0.25, Median ~0.24 → Most cells have moderate-low enrichment
- **lcmv**: Mean ~0.83, Median ~0.83 → Most cells have very high enrichment! ✅

### Plot 2: Enrichment Scores by Cell Type - Boxplot (Top Right)

**What it shows**: Distribution of enrichment scores for the top 10 most abundant cell types

**How to interpret**:
- **High boxes** → That cell type is well represented in your reference
- **Low boxes** → That cell type is poorly represented (might be missing from reference)
- **Wide boxes** → High variability (some cells of this type well represented, others not)
- **Narrow boxes** → Consistent representation across all cells of this type

**What to look for**:
- Are all cell types above a certain threshold (e.g., >0.5)?
- Are there any cell types with very low enrichment (bottom of plot)?
- Is there high variability (wide boxes)?

### Plot 3: Mean Enrichment by Cell Type (Bottom Left)

**What it shows**: Bar chart of mean enrichment score for top 15 cell types (sorted by enrichment)

**How to interpret**:
- **Taller bars** → Better represented cell types
- **Shorter bars** → Poorly represented cell types
- **Order**: Sorted from highest to lowest enrichment

**What to look for**:
- Which cell types have the highest enrichment? (These are well captured)
- Which cell types have the lowest enrichment? (These might be missing or rare)
- Is there a big gap between highest and lowest? (Indicates uneven representation)

**Example from your results**:
- **mcc_01**: fibroblast (0.275) > chondrocyte (0.264) > osteoblast (0.254) → All relatively low
- **lcmv**: Monocytes (0.852) > pDC_like (0.851) > early_and_mature_Granulocytesta (0.848) → All very high! ✅

### Plot 4: Cumulative Distribution (Bottom Right)

**What it shows**: What fraction of cells have enrichment scores below a given value

**How to interpret**:
- **Steep curve** → Most cells have similar enrichment scores
- **Gradual curve** → Wide range of enrichment scores
- **5th percentile line** (red): 5% of cells have enrichment below this value
- **95th percentile line** (blue): 95% of cells have enrichment below this value

**What to look for**:
- Where does the curve start? (Low start = many cells with low enrichment)
- How steep is the curve? (Steep = consistent representation)
- Gap between 5th and 95th percentile? (Large gap = high variability)

## Interpreting the Summary CSV

The `*_summary.csv` file contains per-cell-type statistics:

### Columns Explained

1. **cell_type**: Name of the cell type
2. **mean_enrichment**: Average enrichment score for this cell type
   - Higher = better represented
   - Compare across cell types to see which are well/poorly represented
3. **median_enrichment**: Median enrichment score (less affected by outliers)
4. **n_cells**: Number of cells of this type in the query set
5. **pct_high_enrichment**: Percentage of cells with enrichment > 75th percentile
   - High percentage = most cells of this type are well represented

### How to Use This

1. **Sort by mean_enrichment** (descending) to see:
   - Which cell types are best represented
   - Which cell types are worst represented

2. **Look for patterns**:
   - Are rare cell types (low n_cells) well represented? (Good sign for SPS!)
   - Are common cell types well represented? (Expected)
   - Are there any cell types with very low enrichment? (Might be missing from reference)

3. **Compare datasets**:
   - mcc_01: Lower enrichment scores overall (~0.25)
   - lcmv: Higher enrichment scores overall (~0.83)
   - This suggests lcmv SPS samples are better references (or the dataset is more homogeneous)

## Key Questions to Answer

### 1. Are SPS samples good references?

**Good signs**:
- ✅ High mean enrichment (>0.7)
- ✅ Most cell types have reasonable enrichment (>0.5)
- ✅ Low variability (narrow distributions)
- ✅ Even rare cell types have some enrichment

**Bad signs**:
- ❌ Low mean enrichment (<0.5)
- ❌ Many cell types with very low enrichment (<0.3)
- ❌ High variability (wide distributions)
- ❌ Some cell types with zero enrichment (completely missing)

### 2. Which cell types are well represented?

Look at the **Mean Enrichment by Cell Type** plot:
- Top 15 cell types = best represented
- Bottom of the list = worst represented

### 3. Are there missing cell types?

Look for:
- Cell types with enrichment < 0.3 (very low)
- Cell types with high variability (wide boxes in boxplot)
- Cell types that appear in query but have zero enrichment

### 4. How does this compare to random sampling?

(You would need to run random sampling jobs to compare, but they failed due to pickle issues)

Expected comparison:
- **SPS should have**: Higher mean enrichment, better rare cell coverage, more consistent scores
- **Random should have**: Lower mean enrichment, poorer rare cell coverage, more variable scores

## Example Interpretation: Your Results

### mcc_01 Results

**Overall**: Mean enrichment ~0.25 (moderate-low)
- Most cell types have enrichment between 0.24-0.28
- Relatively consistent across cell types
- Suggests: SPS samples represent the dataset moderately well, but there's room for improvement

**Best represented**: fibroblast (0.275), chondrocyte (0.264)
**Worst represented**: mesodermal cell (0.240)

### lcmv Results

**Overall**: Mean enrichment ~0.83 (very high!) ✅
- Most cell types have enrichment >0.80
- Very consistent across cell types
- Suggests: SPS samples are excellent references for this dataset!

**Best represented**: Monocytes (0.852), pDC_like (0.851)
**Worst represented**: CD4_LCMV_spec (0.778) - but still quite high!

## What This Tells You About Your Sampling

### If enrichment is high (like lcmv):
✅ **SPS sampling worked well!**
- The 100k samples capture the diversity of the full dataset
- Most cell types are well represented
- You can confidently use these samples as references

### If enrichment is moderate (like mcc_01):
⚠️ **SPS sampling is okay but could be better**
- The samples represent the dataset reasonably well
- Some cell types might be underrepresented
- Consider: Are there specific rare cell types missing?

### If enrichment is low:
❌ **SPS sampling might have missed important cell types**
- Many cells are poorly represented
- Some cell types might be completely missing
- Consider: Increase sample size? Adjust sampling parameters?

## Next Steps

1. **Examine the plots** to see which cell types are well/poorly represented
2. **Check the summary CSV** for specific cell types of interest
3. **Compare with random sampling** (if you can fix the pickle issues)
4. **Identify gaps**: Which cell types have low enrichment?
5. **Validate**: Do low-enrichment cell types actually exist in your reference samples?

## Summary

- **Enrichment score** = How similar a query cell is to your reference cells
- **High scores** = Good representation
- **Low scores** = Poor representation (might be missing cell types)
- **Your lcmv results** show excellent representation (mean ~0.83)
- **Your mcc_01 results** show moderate representation (mean ~0.25)

The goal is to have high enrichment scores across most cells and cell types, indicating that your SPS samples are good references that capture the diversity of your dataset!



