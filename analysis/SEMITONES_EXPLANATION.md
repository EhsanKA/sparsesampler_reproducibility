# Understanding SEMITONES and Using SPS Samples as References

## What is SEMITONES?

**SEMITONES** stands for **Similarity Enrichment for Multiple Interrogation Techniques in Omics via Neighborhoods of Exemplars**. It's a computational tool for analyzing single-cell RNA-seq data using a **reference-based approach**.

### Core Concept

SEMITONES works by:
1. **Using reference cells** - A set of well-characterized cells that serve as "exemplars" or benchmarks
2. **Comparing query cells** - All other cells in your dataset are compared against these references
3. **Calculating enrichment scores** - For each query cell, SEMITONES calculates how similar it is to the reference cells
4. **Enabling downstream analysis** - These scores can be used for cell type identification, trajectory analysis, and more

## How SEMITONES Works (Step by Step)

### 1. Reference Cells (Your SPS Samples)
- **What they are**: A curated set of cells that represent the diversity of your dataset
- **In your case**: The cells selected by SPS (SparseSampler) sampling
- **Why they work**: SPS samples are designed to preserve rare cell types and maintain diversity, making them ideal references

### 2. Query Cells (The Rest of Your Dataset)
- **What they are**: All other cells in your full dataset that are NOT in the reference
- **Purpose**: These cells will be compared against the references to understand their similarity

### 3. Similarity Calculation
SEMITONES calculates:
- **Pairwise similarities**: How similar each query cell is to each reference cell
- Uses metrics like cosine similarity, Euclidean distance, etc.
- Creates a similarity matrix: `(n_query_cells × n_reference_cells)`

### 4. Enrichment Scoring
For each query cell, SEMITONES calculates:
- **Enrichment scores**: A score indicating how well-represented that cell is by the reference set
- Higher scores = the cell is well-represented by references
- Lower scores = the cell is poorly represented (might be a rare/novel cell type)

### 5. Downstream Analysis
The enrichment scores enable:
- **Cell type identification**: Find cells similar to known reference types
- **Rare cell detection**: Identify cells with low enrichment (novel/rare types)
- **Trajectory analysis**: Understand developmental relationships
- **Quality control**: Identify outliers or problematic cells

## What Your Boss Meant: Using SPS Samples as References

### The Idea

Your boss suggested using **SPS-sampled cells as reference cells** because:

1. **SPS preserves diversity**: SPS sampling is designed to maintain rare cell types and overall diversity
2. **Representative samples**: The SPS samples should capture the full spectrum of cell types in your dataset
3. **Validation approach**: By using SPS samples as references, you can:
   - Validate that SPS samples are indeed representative
   - Show that other cells in the dataset are well-represented by SPS samples
   - Demonstrate the quality of your sampling method

### The Workflow

```
Full Dataset (e.g., 3M cells)
    │
    ├── SPS Sampled Cells (e.g., 200K cells) → REFERENCE CELLS
    │
    └── Remaining Cells (e.g., 2.8M cells) → QUERY CELLS
            │
            └── SEMITONES compares query cells to reference cells
                    │
                    └── Enrichment scores show how well represented
                        each query cell is by the reference set
```

## Downstream Analysis as Suggested in SEMITONES Paper

Based on the SEMITONES methodology, here's what downstream analysis typically includes:

### 1. Enrichment Score Distribution
- **What**: Histogram/boxplot of enrichment scores across all query cells
- **Why**: Shows how well the reference set represents the full dataset
- **Interpretation**: 
  - High scores = good representation
  - Low scores = poor representation (might indicate missing cell types in reference)

### 2. Cell Type Coverage Analysis
- **What**: Calculate enrichment scores per cell type
- **Why**: Shows which cell types are well/poorly represented in the reference
- **Interpretation**:
  - High enrichment for a cell type = that type is well-represented in SPS samples
  - Low enrichment = that type might be missing or underrepresented

### 3. Rare Cell Type Detection
- **What**: Identify cells with very low enrichment scores
- **Why**: These might be rare cell types not captured in the reference
- **Action**: Check if these are truly rare types or if SPS sampling missed them

### 4. Similarity-Based Clustering
- **What**: Use enrichment scores to cluster cells
- **Why**: Groups cells based on their similarity to reference cells
- **Result**: Can identify cell type relationships and transitions

### 5. Trajectory Analysis
- **What**: Use enrichment scores to infer developmental trajectories
- **Why**: Cells with similar enrichment patterns might be developmentally related
- **Application**: Understand cell state transitions

### 6. Quality Metrics
- **What**: Calculate summary statistics
  - Mean/median enrichment scores
  - Fraction of cells with high enrichment
  - Coverage of cell types
- **Why**: Quantify how well SPS samples represent the full dataset

## Practical Implementation

### Step 1: Prepare Your Data
```python
# Load full dataset
adata_full = load_ref_data('mcc_01', ref=30, base_path=data_path)

# Load SPS sampled indices (these become references)
sps_indices = load_sampling_indices('mcc_01', ref=30, method='sps', 
                                    size=200000, rep=0, data_path=data_path)

# Split into reference and query
reference_cells = adata_full[sps_indices]  # SPS samples = references
query_cells = adata_full[~sps_indices]     # Rest = queries
```

### Step 2: Run SEMITONES
```python
from SEMITONES.enrichment_scoring import calculate_escores

# Calculate enrichment scores
enrichment_scores = calculate_escores(
    X=reference_cells.X,      # Reference: SPS samples
    query=query_cells.X,       # Query: all other cells
    metric='cosine'
)
```

### Step 3: Analyze Results
```python
# 1. Overall enrichment distribution
import matplotlib.pyplot as plt
plt.hist(enrichment_scores.flatten(), bins=50)
plt.xlabel('Enrichment Score')
plt.ylabel('Number of Cells')
plt.title('Distribution of Enrichment Scores')

# 2. Per-cell-type enrichment
for cell_type in unique_cell_types:
    type_mask = query_cells.obs['celltype'] == cell_type
    type_enrichment = enrichment_scores[type_mask].mean()
    print(f"{cell_type}: {type_enrichment:.3f}")

# 3. Identify poorly represented cells
low_enrichment_threshold = np.percentile(enrichment_scores, 5)
poorly_represented = enrichment_scores < low_enrichment_threshold
```

## What This Tells You About Your Sampling

### If SPS Samples Are Good References:
✅ **High enrichment scores** across most query cells
✅ **Good cell type coverage** - all types have reasonable enrichment
✅ **Low variance** in enrichment scores
✅ **Rare cells are represented** - even rare types have some enrichment

### If SPS Samples Are Poor References:
❌ **Low enrichment scores** for many cells
❌ **Missing cell types** - some types have very low enrichment
❌ **High variance** - large differences between cell types
❌ **Rare cells missing** - rare types have zero enrichment

## Comparison with Other Sampling Methods

You can compare how well different sampling methods work as references:

```python
# Compare SPS vs Random as references
sps_enrichment = calculate_escores(X=sps_samples.X, query=query_cells.X)
random_enrichment = calculate_escores(X=random_samples.X, query=query_cells.X)

# SPS should have:
# - Higher mean enrichment (better representation)
# - Better coverage of rare cell types
# - More consistent scores across cell types
```

## Key Insights from This Analysis

1. **Validation of Sampling Quality**: If SPS samples work well as references, it validates that your sampling method preserves diversity

2. **Identification of Gaps**: Low enrichment scores reveal which cell types might be missing from your samples

3. **Quantitative Metrics**: Provides numerical evidence that SPS sampling is better than random sampling

4. **Biological Insights**: Enrichment patterns can reveal relationships between cell types

## Next Steps

1. **Run SEMITONES analysis** on your SPS samples
2. **Compare with random sampling** as a baseline
3. **Analyze enrichment by cell type** to see coverage
4. **Identify any missing cell types** in your samples
5. **Create visualizations** showing enrichment distributions
6. **Report metrics** in your paper showing SPS samples are good references

## Summary

**SEMITONES** uses reference cells to evaluate how well they represent a dataset. By using **SPS samples as references**, you can:
- Validate that SPS sampling preserves diversity
- Show quantitative evidence that SPS is better than random sampling
- Identify any gaps in your sampling
- Provide downstream analysis metrics as suggested in the SEMITONES paper

This approach gives you a rigorous, quantitative way to demonstrate the quality of your sampling method!

