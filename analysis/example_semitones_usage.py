#!/usr/bin/env python
# coding: utf-8

"""
Simple example of using SEMITONES with sampled data as references.

This is a minimal example showing how to use SEMITONES to evaluate
how well sampled cells represent the full dataset.
"""

import os
import sys
import numpy as np
import scanpy as sc

# Add analysis directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

from refined_rare_cell_type_definition import load_reference_data as load_ref_data
from semitones_analysis import load_sampling_indices, get_data_path

# Try to import SEMITONES
try:
    from SEMITONES.enrichment_scoring import calculate_escores
    from SEMITONES.support_funcs import pairwise_similarities
    SEMITONES_AVAILABLE = True
except ImportError as e:
    SEMITONES_AVAILABLE = False
    print(f"Error: SEMITONES not available: {e}")
    print("Please activate the semitones_env conda environment:")
    print("  conda activate semitones_env")
    sys.exit(1)

def example_usage():
    """Simple example of using SEMITONES."""
    
    # Configuration
    dataset = 'mcc_01'
    ref = 30
    method = 'sps'
    size = 200000
    rep = 0
    
    print(f"Example: Using {method} sampled cells as references for {dataset}")
    print("="*60)
    
    # Get data path
    data_path = get_data_path()
    
    # Load full dataset
    print(f"\n1. Loading full dataset ({dataset}, ref {ref})...")
    adata_full = load_ref_data(dataset, ref, base_path=data_path)
    print(f"   Full dataset shape: {adata_full.shape}")
    
    # Load sampling indices (these will be our references)
    print(f"\n2. Loading sampled indices (method={method}, size={size}, rep={rep})...")
    indices = load_sampling_indices(dataset, ref, method, size, rep, data_path)
    
    if indices is None:
        print(f"   Error: No indices found!")
        return
    
    print(f"   Loaded {len(indices)} reference cell indices")
    
    # Prepare reference and query data
    print(f"\n3. Preparing reference and query data...")
    if isinstance(indices[0], (str, np.str_)):
        ref_mask = adata_full.obs.index.isin(indices)
        ref_idx = np.where(ref_mask)[0]
    else:
        ref_idx = indices
    
    query_idx = np.setdiff1d(np.arange(adata_full.shape[0]), ref_idx)
    
    X_ref = adata_full[ref_idx].X
    X_query = adata_full[query_idx].X
    
    # Convert to dense if sparse
    import scipy.sparse
    if scipy.sparse.issparse(X_ref):
        X_ref = X_ref.toarray()
    if scipy.sparse.issparse(X_query):
        X_query = X_query.toarray()
    
    print(f"   Reference cells: {len(ref_idx)}")
    print(f"   Query cells: {len(query_idx)}")
    
    # Calculate enrichment scores using SEMITONES
    print(f"\n4. Calculating enrichment scores with SEMITONES...")
    enrichment_scores = calculate_escores(
        X=X_ref,  # Reference data (sampled cells)
        query=X_query,  # Query data (all other cells)
        metric='cosine',
        ncpu=1
    )
    
    print(f"   Enrichment scores shape: {enrichment_scores.shape}")
    print(f"   Mean enrichment score: {np.mean(enrichment_scores):.4f}")
    print(f"   Median enrichment score: {np.median(enrichment_scores):.4f}")
    
    # Calculate pairwise similarities
    print(f"\n5. Calculating pairwise similarities...")
    similarity_scores = pairwise_similarities(
        X_query,
        X_ref,
        metric='cosine'
    )
    
    print(f"   Similarity scores shape: {similarity_scores.shape}")
    print(f"   Mean similarity: {np.mean(similarity_scores):.4f}")
    print(f"   Median similarity: {np.median(similarity_scores):.4f}")
    
    # Analyze by cell type
    print(f"\n6. Analyzing by cell type...")
    cell_types = adata_full.obs['celltype'].values
    query_cell_types = cell_types[query_idx]
    
    unique_types = np.unique(query_cell_types)
    print(f"   Found {len(unique_types)} cell types in query set")
    
    for cell_type in unique_types[:5]:  # Show first 5 types
        type_mask = query_cell_types == cell_type
        if np.any(type_mask):
            type_similarity = np.mean(similarity_scores[type_mask, :])
            print(f"   {cell_type}: mean similarity = {type_similarity:.4f}")
    
    print("\n" + "="*60)
    print("Example completed successfully!")
    print("\nTo run full analysis, use:")
    print(f"  python semitones_analysis.py --dataset {dataset} --ref {ref} --method {method} --size {size} --rep {rep}")

if __name__ == '__main__':
    example_usage()

