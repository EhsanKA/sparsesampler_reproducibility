# SEMITONES Summary: How It Works and How to Use SPS Samples

## Quick Answer

**SEMITONES** evaluates how well a set of **reference cells** represents your dataset by:
1. Using your **SPS-sampled cells as references** (the "exemplars")
2. Comparing all other cells (**query cells**) to these references
3. Calculating **enrichment scores** showing how well each query cell is represented
4. Providing **downstream analysis** metrics to validate your sampling

## The Big Picture

```
Your Full Dataset (e.g., 3M cells)
│
├── SPS Sampled Cells (200K) → REFERENCE CELLS
│   └── These are your "exemplars" - well-characterized, diverse samples
│
└── Remaining Cells (2.8M) → QUERY CELLS
    └── SEMITONES asks: "How similar is each query cell to the references?"
        └── High similarity = well represented by SPS samples ✅
        └── Low similarity = poorly represented (maybe missing cell type?) ⚠️
```

## What Your Boss Meant

**"Use SPS samples as reference cells"** means:
- Your SPS-sampled cells become the **reference set** (the "gold standard")
- All other cells are **query cells** that get compared to this reference
- This validates that SPS sampling captured the diversity of your dataset
- If query cells have high enrichment scores → SPS samples are good references ✅
- If query cells have low enrichment scores → SPS samples might be missing something ⚠️

## How SEMITONES Works (Simple Explanation)

### Step 1: Similarity Calculation
For each query cell, SEMITONES calculates:
- How similar is it to each reference cell?
- Uses gene expression profiles (cosine similarity, Euclidean distance, etc.)
- Creates a similarity matrix: `(n_query_cells × n_reference_cells)`

### Step 2: Enrichment Scoring
For each query cell, SEMITONES calculates:
- **Enrichment score** = how well that cell is represented by the reference set
- High score = cell is very similar to references (well represented)
- Low score = cell is different from references (might be rare/novel type)

### Step 3: Downstream Analysis
SEMITONES enables:
- **Cell type coverage**: Which types are well/poorly represented?
- **Rare cell detection**: Which cells have very low enrichment?
- **Quality metrics**: Overall how good are your references?

## Downstream Analysis (As in SEMITONES Paper)

### 1. Enrichment Score Distribution
- **What**: Histogram showing distribution of enrichment scores
- **Why**: Shows overall quality of reference set
- **Good sign**: Most cells have high enrichment scores

### 2. Per-Cell-Type Analysis
- **What**: Calculate mean enrichment for each cell type
- **Why**: Shows which types are well represented
- **Good sign**: All cell types have reasonable enrichment

### 3. Rare Cell Identification
- **What**: Find cells with very low enrichment scores
- **Why**: These might be rare types not in your reference
- **Action**: Check if these are truly rare or if sampling missed them

### 4. Comparison Metrics
- **What**: Compare SPS vs Random sampling as references
- **Why**: Quantitatively show SPS is better
- **Expected**: SPS should have higher mean enrichment

## What This Tells You

### If SPS Samples Are Good References:
✅ High enrichment scores across most cells
✅ Good coverage of all cell types
✅ Low variance (consistent representation)
✅ Even rare cell types have some enrichment

### If SPS Samples Are Poor References:
❌ Low enrichment for many cells
❌ Some cell types have zero enrichment (missing!)
❌ High variance (inconsistent)
❌ Many rare cells have very low enrichment

## Practical Workflow

### 1. Prepare Data
```python
# Load full dataset
adata_full = load_ref_data('mcc_01', ref=30)

# Load SPS sample indices
sps_indices = load_sampling_indices('mcc_01', ref=30, method='sps', 
                                    size=200000, rep=0)

# Split into reference (SPS) and query (rest)
reference_cells = adata_full[sps_indices]
query_cells = adata_full[~sps_indices]
```

### 2. Run SEMITONES
```python
from SEMITONES.enrichment_scoring import calculate_escores

# Calculate enrichment scores
# Note: Check SEMITONES documentation for exact API usage
enrichment_scores = calculate_escores(
    X=reference_cells.X,  # Reference data
    query=query_cells.X,  # Query data (check API)
    metric='cosine'
)
```

### 3. Analyze Results
```python
# Overall distribution
plt.hist(enrichment_scores.flatten())

# Per cell type
for cell_type in unique_types:
    type_enrichment = enrichment_scores[cell_type_mask].mean()
    print(f"{cell_type}: {type_enrichment:.3f}")

# Identify poorly represented
low_enrichment = enrichment_scores < threshold
```

## Key Insights

1. **Validation**: If SPS samples work well as references, it validates your sampling method
2. **Gap Detection**: Low enrichment reveals missing cell types
3. **Quantitative Proof**: Provides numbers showing SPS > Random
4. **Biological Insights**: Enrichment patterns reveal cell type relationships

## Files Created

1. **`SEMITONES_EXPLANATION.md`** - Detailed explanation of SEMITONES
2. **`semitones_analysis.py`** - Basic SEMITONES analysis script
3. **`semitones_downstream_analysis.py`** - Full downstream analysis following SEMITONES paper
4. **`example_semitones_usage.py`** - Simple example script

## Next Steps

1. **Read** `SEMITONES_EXPLANATION.md` for detailed understanding
2. **Try** `example_semitones_usage.py` to see it in action
3. **Run** `semitones_downstream_analysis.py` for full analysis
4. **Compare** SPS vs Random sampling as references
5. **Report** metrics showing SPS samples are good references

## Summary

**SEMITONES** = Tool to evaluate how well reference cells represent a dataset
**Your SPS samples** = The reference cells (the "exemplars")
**All other cells** = Query cells (compared to references)
**Enrichment scores** = How well each query cell is represented
**Downstream analysis** = Metrics and visualizations showing sampling quality

This approach gives you **rigorous, quantitative evidence** that your SPS sampling method preserves diversity and creates representative samples! 🎯

