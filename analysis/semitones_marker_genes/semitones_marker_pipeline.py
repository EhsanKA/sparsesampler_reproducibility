#!/usr/bin/env python
# coding: utf-8

"""
SEMITONES Marker Gene Pipeline

This script identifies marker genes for each cell type in the MCC dataset using SEMITONES,
comparing results when using SPS 100k vs Random 100k reference cells.

The pipeline:
1. Loads the full MCC dataset and sampling indices
2. Runs SEMITONES enrichment scoring with sampled cells as references
3. Identifies marker genes per cell type based on enrichment scores
4. Compares results against literature-based marker genes

Usage:
    conda activate semitones_env
    python semitones_marker_pipeline.py --method sps
    python semitones_marker_pipeline.py --method random

Requirements:
    - SEMITONES library (pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip)
    - scanpy, pandas, numpy, scipy
"""

import os
import sys
import argparse
import pickle
import logging
from datetime import datetime
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse

# Add parent directory to path for imports
script_dir = os.path.dirname(os.path.abspath(__file__))
analysis_dir = os.path.dirname(script_dir)
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

# Import literature markers
from semitones_marker_genes.literature_markers import (
    LITERATURE_MARKERS,
    get_markers_for_celltype,
    get_all_markers
)

# Try to import SEMITONES
try:
    from SEMITONES.enrichment_scoring import calculate_escores, permute, sig_interval, pvals_per_cell
    from SEMITONES.support_funcs import pairwise_similarities
    SEMITONES_AVAILABLE = True
except ImportError as e:
    SEMITONES_AVAILABLE = False
    print(f"Warning: SEMITONES not available: {e}")
    print("Please activate the semitones_env conda environment:")
    print("  conda activate semitones_env")

# ============================================================================
# Configuration
# ============================================================================

# Paths
PROJECT_ROOT = os.environ.get(
    'PROJECT_ROOT',
    '/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility'
)
DATA_PATH = os.path.join(PROJECT_ROOT, 'data', 'mcc', 'benchmark', '30')
OUTPUT_DIR = os.path.join(script_dir, 'results')

# Analysis parameters
SAMPLE_SIZE = 100000
REP = 0
LABEL_KEY = 'celltype'
N_PERMUTATIONS = 100  # For significance testing
N_SDS = 5  # Number of standard deviations for significance threshold
TOP_N_GENES = 50  # Number of top marker genes to report per cell type
# Cell type balancing for fair comparison
RARE_CELL_TYPE = 'osteoblast'  # The rare cell type to preserve
# Note: We'll select all rare cells and equal numbers from other types
# This ensures rare cells are well-represented in the analysis


# ============================================================================
# Logging Setup
# ============================================================================

def setup_logger(method):
    """Set up logging for the analysis."""
    log_dir = os.path.join(PROJECT_ROOT, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'semitones_marker_gene_{timestamp}.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


# ============================================================================
# Data Loading Functions
# ============================================================================

def load_mcc_dataset():
    """Load the full MCC dataset."""
    adata_path = os.path.join(DATA_PATH, 'adata.h5ad')
    logger.info(f"Loading MCC dataset from {adata_path}")
    adata = sc.read_h5ad(adata_path)
    logger.info(f"Dataset shape: {adata.shape}")
    logger.info(f"Cell types: {adata.obs[LABEL_KEY].unique().tolist()}")
    return adata


def load_sampling_indices(method):
    """Load sampling indices for the specified method."""
    base_path = os.path.join(DATA_PATH, method, str(SAMPLE_SIZE), str(REP))
    
    # Try loading from npy file first (more reliable across numpy versions)
    npy_path = os.path.join(base_path, 'indices.npy')
    pkl_path = os.path.join(base_path, 'results.pkl')
    
    if os.path.exists(npy_path):
        logger.info(f"Loading {method} indices from {npy_path}")
        indices = np.load(npy_path)
        logger.info(f"Loaded {len(indices)} {method} indices from npy file")
        return indices
    
    if os.path.exists(pkl_path):
        logger.info(f"Loading {method} indices from {pkl_path}")
        try:
            with open(pkl_path, 'rb') as f:
                data = pickle.load(f)
            
            # Handle different pickle formats
            if isinstance(data, tuple):
                indices = data[0]
            else:
                indices = data
            
            # Convert to numpy array if needed
            if isinstance(indices, list):
                indices = np.array(indices)
            
            logger.info(f"Loaded {len(indices)} {method} indices from pkl file")
            return indices
        except Exception as e:
            logger.error(f"Failed to load pickle file: {e}")
            raise
    
    raise FileNotFoundError(f"No indices file found at {base_path}")


# ============================================================================
# SEMITONES Analysis Functions
# ============================================================================

def subsample_cells_balanced(adata, reference_indices, rare_cell_type=None, seed=42):
    """
    Subsample cells with balanced representation across cell types.
    
    Strategy:
    - Reference cells: All rare cells + equal number from each other type
    - Query cells: All rare cells + equal number per other type (matching rare count)
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    reference_indices : array-like
        Indices of reference cells from sampling method (SPS or Random 100k)
    rare_cell_type : str, optional
        Name of rare cell type. If None, uses RARE_CELL_TYPE constant.
    seed : int
        Random seed
    
    Returns:
    --------
    tuple : (subsampled adata, new reference indices for SEMITONES)
    """
    if rare_cell_type is None:
        rare_cell_type = RARE_CELL_TYPE
    
    np.random.seed(seed)
    
    n_total = adata.shape[0]
    cell_types = adata.obs[LABEL_KEY].values
    unique_types = np.unique(cell_types)
    other_types = [ct for ct in unique_types if ct != rare_cell_type]
    
    logger.info(f"Balanced subsampling for SEMITONES...")
    logger.info(f"  Original dataset: {n_total} cells")
    logger.info(f"  Original reference cells (from sampling): {len(reference_indices)}")
    logger.info(f"  Rare cell type: {rare_cell_type}")
    logger.info(f"  Other cell types: {other_types}")
    
    # ========== REFERENCE CELLS ==========
    # Strategy: All rare cells from reference + equal number from each other type
    
    ref_cell_types = cell_types[reference_indices]
    
    # Count rare cells in reference
    rare_mask_ref = ref_cell_types == rare_cell_type
    rare_ref_indices = reference_indices[rare_mask_ref]
    n_rare_ref = len(rare_ref_indices)
    
    logger.info(f"\n  REFERENCE CELLS:")
    logger.info(f"    Rare cells ({rare_cell_type}) in reference: {n_rare_ref}")
    
    # For other types, sample equal number as rare cells
    sampled_ref_indices = list(rare_ref_indices)
    
    for ct in other_types:
        ct_mask = ref_cell_types == ct
        ct_ref_indices = reference_indices[ct_mask]
        n_available = len(ct_ref_indices)
        n_sample = min(n_rare_ref, n_available)
        
        if n_sample > 0:
            sampled = np.random.choice(ct_ref_indices, size=n_sample, replace=False)
            sampled_ref_indices.extend(sampled)
            logger.info(f"    {ct}: sampled {n_sample} from {n_available} (target: {n_rare_ref})")
        else:
            logger.info(f"    {ct}: 0 cells available in reference")
    
    reference_indices = np.array(sampled_ref_indices)
    logger.info(f"    Total reference cells: {len(reference_indices)}")
    
    # ========== QUERY CELLS ==========
    # Strategy: All rare cells + equal number per other type (matching rare count in query)
    
    ref_set = set(reference_indices)
    all_indices = np.arange(n_total)
    non_ref_indices = np.array([i for i in all_indices if i not in ref_set])
    non_ref_types = cell_types[non_ref_indices]
    
    # Count rare cells in query pool
    rare_mask_query = non_ref_types == rare_cell_type
    rare_query_indices = non_ref_indices[rare_mask_query]
    n_rare_query = len(rare_query_indices)
    
    logger.info(f"\n  QUERY CELLS:")
    logger.info(f"    Rare cells ({rare_cell_type}) in query pool: {n_rare_query}")
    
    # Select all rare cells from query
    sampled_query_indices = list(rare_query_indices)
    
    # For other types, sample equal number as rare cells in query
    for ct in other_types:
        ct_mask = non_ref_types == ct
        ct_query_indices = non_ref_indices[ct_mask]
        n_available = len(ct_query_indices)
        n_sample = min(n_rare_query, n_available)
        
        if n_sample > 0:
            sampled = np.random.choice(ct_query_indices, size=n_sample, replace=False)
            sampled_query_indices.extend(sampled)
            logger.info(f"    {ct}: sampled {n_sample} from {n_available} (target: {n_rare_query})")
        else:
            logger.info(f"    {ct}: 0 cells available in query pool")
    
    query_indices = np.array(sampled_query_indices)
    logger.info(f"    Total query cells: {len(query_indices)}")
    
    # ========== COMBINE ==========
    all_selected = np.concatenate([reference_indices, query_indices])
    all_selected = np.unique(all_selected)
    all_selected = np.sort(all_selected)
    
    logger.info(f"\n  FINAL:")
    logger.info(f"    Total cells selected: {len(all_selected)}")
    
    # Create mapping from old indices to new indices
    old_to_new = {old: new for new, old in enumerate(all_selected)}
    new_ref_indices = np.array([old_to_new[i] for i in reference_indices])
    
    # Subset adata
    adata_sub = adata[all_selected].copy()
    
    # Log final cell type distribution
    final_types = adata_sub.obs[LABEL_KEY].value_counts()
    logger.info(f"    Final cell type distribution:")
    for ct, count in final_types.items():
        logger.info(f"      {ct}: {count}")
    
    logger.info(f"    Final reference cells for SEMITONES: {len(new_ref_indices)}")
    
    return adata_sub, new_ref_indices


def run_semitones_enrichment(adata, reference_indices, metric='cosine', ncpu=4,
                             rare_cell_type=None):
    """
    Run SEMITONES enrichment scoring using sampled cells as references.
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    reference_indices : array-like
        Indices of reference cells (sampled cells from SPS or Random)
    metric : str
        Similarity metric to use
    ncpu : int
        Number of CPUs for parallel processing
    rare_cell_type : str, optional
        Name of rare cell type for balanced sampling
    
    Returns:
    --------
    tuple : (enrichment scores DataFrame, subsampled adata, new reference indices)
    """
    if not SEMITONES_AVAILABLE:
        raise ImportError("SEMITONES is not available")
    
    logger.info("Preparing data for SEMITONES enrichment scoring...")
    
    # Balanced subsampling: all rare cells + equal number from other types
    adata, reference_indices = subsample_cells_balanced(
        adata, reference_indices, 
        rare_cell_type=rare_cell_type
    )
    
    # Get expression matrix
    X = adata.X
    if scipy.sparse.issparse(X):
        logger.info("Converting sparse matrix to dense...")
        X = X.toarray()
    
    # Ensure X is a pandas DataFrame for SEMITONES
    gene_names = adata.var_names.tolist()
    cell_names = adata.obs_names.tolist()
    
    logger.info(f"Expression matrix shape: {X.shape}")
    logger.info(f"Number of reference cells: {len(reference_indices)}")
    
    # SEMITONES expects the query parameter to be the reference cell indices
    # X contains all cells, query specifies which cells are references
    logger.info("Running SEMITONES enrichment scoring...")
    logger.info(f"  Metric: {metric}")
    logger.info(f"  CPUs: {ncpu}")
    
    # Convert to DataFrame for SEMITONES
    X_df = pd.DataFrame(X, index=cell_names, columns=gene_names)
    
    # SEMITONES expects cell names (row labels), not integer indices
    # Convert reference indices to cell names
    reference_cell_names = [cell_names[i] for i in reference_indices]
    
    logger.info(f"Reference cell names (first 5): {reference_cell_names[:5]}")
    
    # Calculate enrichment scores
    # In SEMITONES, 'query' refers to the reference cells
    escores = calculate_escores(
        X=X_df,
        query=reference_cell_names,
        metric=metric,
        ncpu=ncpu,
        scale_exp=True
    )
    
    logger.info(f"Enrichment scores shape: {escores.shape}")
    
    return escores, adata, reference_indices


def run_permutation_testing(adata, reference_indices, n_permutations=100, metric='cosine', ncpu=4):
    """
    Run permutation testing for significance.
    
    Parameters:
    -----------
    adata : AnnData
        Dataset (possibly subsampled)
    reference_indices : array-like
        Indices of reference cells (in subsampled dataset)
    n_permutations : int
        Number of permutations
    metric : str
        Similarity metric
    ncpu : int
        Number of CPUs
    
    Returns:
    --------
    pd.DataFrame : Permutation enrichment scores
    """
    if not SEMITONES_AVAILABLE:
        raise ImportError("SEMITONES is not available")
    
    logger.info(f"Running permutation testing with {n_permutations} permutations...")
    
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    
    gene_names = adata.var_names.tolist()
    cell_names = adata.obs_names.tolist()
    X_df = pd.DataFrame(X, index=cell_names, columns=gene_names)
    
    # Permute the expression matrix
    X_permuted = permute(X_df, n=n_permutations, axis=0, seed=42)
    
    # Calculate enrichment scores on permuted data
    pscores = calculate_escores(
        X=X_permuted,
        query=list(reference_indices),
        metric=metric,
        ncpu=ncpu,
        scale_exp=True
    )
    
    logger.info(f"Permutation scores shape: {pscores.shape}")
    
    return pscores


def identify_marker_genes_per_celltype(escores, adata, reference_indices, 
                                        pscores=None, n_sds=5, top_n=50):
    """
    Identify marker genes for each cell type based on enrichment scores.
    
    Parameters:
    -----------
    escores : pd.DataFrame
        Enrichment scores (genes x reference cells)
    adata : AnnData
        Full dataset (for cell type labels)
    reference_indices : array-like
        Indices of reference cells
    pscores : pd.DataFrame, optional
        Permutation scores for significance testing
    n_sds : int
        Number of standard deviations for significance threshold
    top_n : int
        Number of top genes to report
    
    Returns:
    --------
    dict : Marker genes per cell type with scores and rankings
    """
    logger.info("Identifying marker genes per cell type...")
    
    # Get cell type labels for reference cells
    ref_celltypes = adata.obs[LABEL_KEY].iloc[reference_indices].values
    unique_celltypes = np.unique(ref_celltypes)
    
    logger.info(f"Cell types in reference set: {unique_celltypes.tolist()}")
    
    # Calculate significance thresholds if permutation scores provided
    sig_thresholds = None
    if pscores is not None:
        logger.info("Calculating significance thresholds from permutation scores...")
        sig_thresholds = sig_interval(pscores, n_sds=n_sds)
    
    results = {}
    
    for celltype in unique_celltypes:
        logger.info(f"Processing cell type: {celltype}")
        
        # Get reference cells of this cell type
        celltype_mask = ref_celltypes == celltype
        celltype_ref_indices = np.array(reference_indices)[celltype_mask]
        n_refs = len(celltype_ref_indices)
        
        logger.info(f"  Number of reference cells: {n_refs}")
        
        if n_refs == 0:
            logger.warning(f"  No reference cells for {celltype}, skipping")
            continue
        
        # Get enrichment scores for this cell type's reference cells
        # escores columns are reference cell indices
        celltype_escores = escores.loc[:, escores.columns.isin(celltype_ref_indices)]
        
        if celltype_escores.shape[1] == 0:
            # Try using positional indices
            col_mask = [i for i, c in enumerate(escores.columns) if c in celltype_ref_indices]
            if len(col_mask) > 0:
                celltype_escores = escores.iloc[:, col_mask]
            else:
                logger.warning(f"  Could not find enrichment scores for {celltype} reference cells")
                continue
        
        # Calculate mean enrichment score per gene across this cell type's references
        mean_escores = celltype_escores.mean(axis=1)
        std_escores = celltype_escores.std(axis=1)
        max_escores = celltype_escores.max(axis=1)
        
        # Rank genes by mean enrichment score
        gene_ranking = mean_escores.sort_values(ascending=False)
        
        # Get top N genes
        top_genes = gene_ranking.head(top_n)
        
        # Determine significance if thresholds available
        significant_genes = None
        if sig_thresholds is not None:
            # A gene is significant if its mean score is above threshold for any reference
            significant_genes = []
            for gene in gene_ranking.index:
                for ref_idx in celltype_ref_indices:
                    if ref_idx in sig_thresholds:
                        lower, upper = sig_thresholds[ref_idx]
                        gene_score = escores.loc[gene, ref_idx] if ref_idx in escores.columns else None
                        if gene_score is not None and (gene_score > upper or gene_score < lower):
                            significant_genes.append(gene)
                            break
        
        results[celltype] = {
            'n_reference_cells': n_refs,
            'top_genes': top_genes.index.tolist(),
            'top_scores': top_genes.values.tolist(),
            'mean_enrichment': mean_escores.to_dict(),
            'std_enrichment': std_escores.to_dict(),
            'gene_ranking': gene_ranking.index.tolist(),
            'significant_genes': significant_genes,
        }
        
        logger.info(f"  Top 5 marker genes: {top_genes.head().index.tolist()}")
    
    return results


def evaluate_literature_markers(marker_results, adata):
    """
    Evaluate how well the identified markers match literature-based markers.
    
    Parameters:
    -----------
    marker_results : dict
        Results from identify_marker_genes_per_celltype
    adata : AnnData
        Full dataset (to check which genes are present)
    
    Returns:
    --------
    dict : Evaluation metrics per cell type
    """
    logger.info("Evaluating against literature-based markers...")
    
    # Get genes present in the dataset
    dataset_genes = set(adata.var_names.tolist())
    
    evaluation = {}
    
    for celltype, results in marker_results.items():
        logger.info(f"Evaluating {celltype}...")
        
        # Get literature markers for this cell type
        lit_markers = get_markers_for_celltype(celltype)
        
        if not lit_markers:
            logger.warning(f"  No literature markers found for {celltype}")
            evaluation[celltype] = {
                'literature_markers': [],
                'literature_markers_in_dataset': [],
                'recall_at_10': None,
                'recall_at_50': None,
                'mean_rank': None,
            }
            continue
        
        # Check which literature markers are in the dataset
        lit_markers_in_dataset = [m for m in lit_markers if m in dataset_genes]
        
        logger.info(f"  Literature markers: {lit_markers}")
        logger.info(f"  Present in dataset: {lit_markers_in_dataset}")
        
        if not lit_markers_in_dataset:
            logger.warning(f"  No literature markers found in dataset for {celltype}")
            evaluation[celltype] = {
                'literature_markers': lit_markers,
                'literature_markers_in_dataset': [],
                'recall_at_10': 0,
                'recall_at_50': 0,
                'mean_rank': None,
            }
            continue
        
        # Get gene ranking
        gene_ranking = results['gene_ranking']
        
        # Calculate recall at different cutoffs
        top_10 = set(gene_ranking[:10])
        top_50 = set(gene_ranking[:50])
        
        recall_10 = len(set(lit_markers_in_dataset) & top_10) / len(lit_markers_in_dataset)
        recall_50 = len(set(lit_markers_in_dataset) & top_50) / len(lit_markers_in_dataset)
        
        # Calculate mean rank of literature markers
        ranks = []
        for marker in lit_markers_in_dataset:
            if marker in gene_ranking:
                rank = gene_ranking.index(marker) + 1  # 1-indexed
                ranks.append(rank)
        
        mean_rank = np.mean(ranks) if ranks else None
        
        # Individual marker ranks
        marker_ranks = {m: gene_ranking.index(m) + 1 if m in gene_ranking else None 
                       for m in lit_markers_in_dataset}
        
        evaluation[celltype] = {
            'literature_markers': lit_markers,
            'literature_markers_in_dataset': lit_markers_in_dataset,
            'n_literature_markers_in_dataset': len(lit_markers_in_dataset),
            'recall_at_10': recall_10,
            'recall_at_50': recall_50,
            'mean_rank': mean_rank,
            'marker_ranks': marker_ranks,
            'markers_in_top_10': [m for m in lit_markers_in_dataset if m in top_10],
            'markers_in_top_50': [m for m in lit_markers_in_dataset if m in top_50],
        }
        
        logger.info(f"  Recall@10: {recall_10:.2%}")
        logger.info(f"  Recall@50: {recall_50:.2%}")
        logger.info(f"  Mean rank: {mean_rank}")
    
    return evaluation


# ============================================================================
# Main Pipeline
# ============================================================================

def run_pipeline(method, skip_permutation=False, ncpu=4, rare_cell_type=None, data_path=None):
    """
    Run the full SEMITONES marker gene pipeline.
    
    Parameters:
    -----------
    method : str
        Sampling method ('sps' or 'random')
    skip_permutation : bool
        Whether to skip permutation testing (faster)
    ncpu : int
        Number of CPUs for parallel processing
    rare_cell_type : str, optional
        Rare cell type name for balanced sampling. Default uses RARE_CELL_TYPE.
    data_path : str, optional
        Path to data directory. If None, uses DATA_PATH constant.
    """
    global DATA_PATH
    if data_path is not None:
        DATA_PATH = data_path
        logger.info(f"Using custom data path: {DATA_PATH}")
    
    if rare_cell_type is None:
        rare_cell_type = RARE_CELL_TYPE
    
    logger.info("=" * 80)
    logger.info(f"SEMITONES Marker Gene Pipeline - {method.upper()}")
    logger.info("=" * 80)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load data
    logger.info("\n" + "=" * 40)
    logger.info("Step 1: Loading Data")
    logger.info("=" * 40)
    
    adata_full = load_mcc_dataset()
    reference_indices = load_sampling_indices(method)
    
    # Run SEMITONES enrichment scoring (with subsampling if needed)
    logger.info("\n" + "=" * 40)
    logger.info("Step 2: SEMITONES Enrichment Scoring")
    logger.info("=" * 40)
    
    escores, adata, reference_indices = run_semitones_enrichment(
        adata_full, 
        reference_indices, 
        metric='cosine',
        ncpu=ncpu,
        rare_cell_type=rare_cell_type
    )
    
    # Save enrichment scores
    escores_path = os.path.join(OUTPUT_DIR, f'{method}_enrichment_scores.pkl')
    with open(escores_path, 'wb') as f:
        pickle.dump(escores, f)
    logger.info(f"Saved enrichment scores to {escores_path}")
    
    # Permutation testing (optional)
    pscores = None
    if not skip_permutation:
        logger.info("\n" + "=" * 40)
        logger.info("Step 3: Permutation Testing")
        logger.info("=" * 40)
        
        pscores = run_permutation_testing(
            adata,
            reference_indices,
            n_permutations=N_PERMUTATIONS,
            metric='cosine',
            ncpu=ncpu
        )
        
        pscores_path = os.path.join(OUTPUT_DIR, f'{method}_permutation_scores.pkl')
        with open(pscores_path, 'wb') as f:
            pickle.dump(pscores, f)
        logger.info(f"Saved permutation scores to {pscores_path}")
    
    # Identify marker genes per cell type
    logger.info("\n" + "=" * 40)
    logger.info("Step 4: Identifying Marker Genes per Cell Type")
    logger.info("=" * 40)
    
    marker_results = identify_marker_genes_per_celltype(
        escores,
        adata,
        reference_indices,
        pscores=pscores,
        n_sds=N_SDS,
        top_n=TOP_N_GENES
    )
    
    # Save marker results
    markers_path = os.path.join(OUTPUT_DIR, f'{method}_marker_genes.pkl')
    with open(markers_path, 'wb') as f:
        pickle.dump(marker_results, f)
    logger.info(f"Saved marker results to {markers_path}")
    
    # Also save as CSV for easy viewing
    marker_df_rows = []
    for celltype, results in marker_results.items():
        for i, (gene, score) in enumerate(zip(results['top_genes'], results['top_scores'])):
            marker_df_rows.append({
                'cell_type': celltype,
                'rank': i + 1,
                'gene': gene,
                'mean_enrichment_score': score,
            })
    
    marker_df = pd.DataFrame(marker_df_rows)
    marker_csv_path = os.path.join(OUTPUT_DIR, f'{method}_marker_genes_by_celltype.csv')
    marker_df.to_csv(marker_csv_path, index=False)
    logger.info(f"Saved marker genes CSV to {marker_csv_path}")
    
    # Evaluate against literature markers
    logger.info("\n" + "=" * 40)
    logger.info("Step 5: Evaluating Against Literature Markers")
    logger.info("=" * 40)
    
    evaluation = evaluate_literature_markers(marker_results, adata)
    
    # Save evaluation results
    eval_path = os.path.join(OUTPUT_DIR, f'{method}_literature_evaluation.pkl')
    with open(eval_path, 'wb') as f:
        pickle.dump(evaluation, f)
    logger.info(f"Saved evaluation to {eval_path}")
    
    # Save evaluation as CSV
    eval_df_rows = []
    for celltype, metrics in evaluation.items():
        eval_df_rows.append({
            'cell_type': celltype,
            'n_literature_markers': len(metrics.get('literature_markers', [])),
            'n_in_dataset': metrics.get('n_literature_markers_in_dataset', 0),
            'recall_at_10': metrics.get('recall_at_10'),
            'recall_at_50': metrics.get('recall_at_50'),
            'mean_rank': metrics.get('mean_rank'),
            'markers_in_top_10': ', '.join(metrics.get('markers_in_top_10', [])),
            'markers_in_top_50': ', '.join(metrics.get('markers_in_top_50', [])),
        })
    
    eval_df = pd.DataFrame(eval_df_rows)
    eval_csv_path = os.path.join(OUTPUT_DIR, f'{method}_literature_evaluation.csv')
    eval_df.to_csv(eval_csv_path, index=False)
    logger.info(f"Saved evaluation CSV to {eval_csv_path}")
    
    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    
    for celltype, metrics in evaluation.items():
        logger.info(f"\n{celltype}:")
        logger.info(f"  Literature markers in dataset: {metrics.get('n_literature_markers_in_dataset', 0)}")
        logger.info(f"  Recall@10: {metrics.get('recall_at_10', 'N/A')}")
        logger.info(f"  Recall@50: {metrics.get('recall_at_50', 'N/A')}")
        logger.info(f"  Mean rank: {metrics.get('mean_rank', 'N/A')}")
    
    logger.info("\n" + "=" * 80)
    logger.info("Pipeline completed successfully!")
    logger.info("=" * 80)
    
    return {
        'enrichment_scores': escores,
        'marker_results': marker_results,
        'evaluation': evaluation,
    }


# ============================================================================
# Entry Point
# ============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='SEMITONES Marker Gene Pipeline for MCC Dataset'
    )
    parser.add_argument(
        '--method', 
        type=str, 
        required=True,
        choices=['sps', 'random'],
        help='Sampling method (sps or random)'
    )
    parser.add_argument(
        '--skip-permutation',
        action='store_true',
        help='Skip permutation testing (faster but no significance)'
    )
    parser.add_argument(
        '--ncpu',
        type=int,
        default=4,
        help='Number of CPUs for parallel processing (default: 4)'
    )
    parser.add_argument(
        '--rare-cell-type',
        type=str,
        default=None,
        help=f'Rare cell type name for balanced sampling (default: {RARE_CELL_TYPE})'
    )
    parser.add_argument(
        '--data-path',
        type=str,
        default=None,
        help='Path to data directory (default: uses DATA_PATH constant)'
    )
    
    args = parser.parse_args()
    
    # Check SEMITONES availability
    if not SEMITONES_AVAILABLE:
        print("ERROR: SEMITONES is not available.")
        print("Please activate the semitones_env conda environment:")
        print("  conda activate semitones_env")
        sys.exit(1)
    
    # Setup logging
    logger = setup_logger(args.method)
    
    # Run pipeline
    try:
        results = run_pipeline(
            method=args.method,
            skip_permutation=args.skip_permutation,
            ncpu=args.ncpu,
            rare_cell_type=args.rare_cell_type,
            data_path=args.data_path
        )
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
