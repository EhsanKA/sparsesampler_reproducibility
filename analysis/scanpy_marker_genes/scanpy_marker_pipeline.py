#!/usr/bin/env python
# coding: utf-8

"""
Scanpy Marker Gene Pipeline

This script identifies marker genes for each cell type in the MCC dataset using Scanpy's
rank_genes_groups, comparing results when using Full dataset vs SPS 100k vs Random 100k.

The pipeline:
1. Loads the full MCC dataset and sampling indices
2. Subsets to sampled cells (or uses full dataset)
3. Runs Scanpy rank_genes_groups with multiple DE methods
4. Identifies marker genes per cell type
5. Compares results against literature-based marker genes

Usage:
    conda activate facs_sampling
    python scanpy_marker_pipeline.py --method full --de-methods wilcoxon t-test logreg
    python scanpy_marker_pipeline.py --method sps --de-methods wilcoxon t-test logreg
    python scanpy_marker_pipeline.py --method random --de-methods wilcoxon t-test logreg

Requirements:
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
import warnings

warnings.filterwarnings('ignore')

# Add parent directory to path for imports
script_dir = os.path.dirname(os.path.abspath(__file__))
analysis_dir = os.path.dirname(script_dir)
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

# Import literature markers from semitones_marker_genes
from semitones_marker_genes.literature_markers import (
    LITERATURE_MARKERS,
    get_markers_for_celltype,
    get_all_markers
)

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
TOP_N_GENES = 50  # Number of top marker genes to report per cell type

# Available DE methods
AVAILABLE_DE_METHODS = ['wilcoxon', 't-test', 't-test_overestim_var', 'logreg']

# Initialize logger as None (will be set up in main)
logger = None


# ============================================================================
# Logging Setup
# ============================================================================

def setup_logger(method):
    """Set up logging for the analysis."""
    global logger
    log_dir = os.path.join(PROJECT_ROOT, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'scanpy_marker_gene_{method}_{timestamp}.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)
    return logger


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
    
    # Log cell type distribution
    cell_type_counts = adata.obs[LABEL_KEY].value_counts()
    logger.info("Cell type distribution:")
    for ct, count in cell_type_counts.items():
        logger.info(f"  {ct}: {count:,} cells ({100*count/adata.n_obs:.2f}%)")
    
    return adata


def load_sampling_indices(method):
    """Load sampling indices for the specified method."""
    if method == 'full':
        return None  # Use full dataset
    
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


def subset_to_sampled_cells(adata, indices):
    """Subset AnnData to sampled cells."""
    if indices is None:
        logger.info("Using full dataset (no subsetting)")
        return adata
    
    logger.info(f"Subsetting to {len(indices)} sampled cells...")
    adata_sub = adata[indices].copy()
    
    # Log subsetted cell type distribution
    cell_type_counts = adata_sub.obs[LABEL_KEY].value_counts()
    logger.info("Subsetted cell type distribution:")
    for ct, count in cell_type_counts.items():
        logger.info(f"  {ct}: {count:,} cells ({100*count/adata_sub.n_obs:.2f}%)")
    
    return adata_sub


# ============================================================================
# Scanpy Marker Gene Analysis Functions
# ============================================================================

def run_marker_gene_analysis(adata, de_method='wilcoxon'):
    """
    Run Scanpy rank_genes_groups for marker gene identification.
    
    Parameters:
    -----------
    adata : AnnData
        Dataset to analyze
    de_method : str
        Differential expression method ('wilcoxon', 't-test', 't-test_overestim_var', 'logreg')
    
    Returns:
    --------
    AnnData : adata with rank_genes_groups results stored
    """
    logger.info(f"Running marker gene analysis with method: {de_method}")
    
    # Check if we need to filter genes with zero variance
    # This can cause issues with some DE methods
    if scipy.sparse.issparse(adata.X):
        gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
    else:
        gene_var = np.var(adata.X, axis=0)
    
    n_zero_var = np.sum(gene_var == 0)
    if n_zero_var > 0:
        logger.warning(f"Found {n_zero_var} genes with zero variance")
    
    # Run rank_genes_groups
    key_added = f'rank_genes_{de_method}'
    
    try:
        sc.tl.rank_genes_groups(
            adata,
            groupby=LABEL_KEY,
            method=de_method,
            pts=True,  # Include fraction of cells expressing
            key_added=key_added,
            n_genes=adata.n_vars  # Rank all genes
        )
        logger.info(f"Successfully ran rank_genes_groups with {de_method}")
    except Exception as e:
        logger.error(f"Error running rank_genes_groups with {de_method}: {e}")
        raise
    
    return adata


def extract_marker_genes(adata, de_method='wilcoxon', top_n=50):
    """
    Extract top marker genes per cell type from rank_genes_groups results.
    
    Parameters:
    -----------
    adata : AnnData
        Dataset with rank_genes_groups results
    de_method : str
        DE method used (to find the correct key)
    top_n : int
        Number of top genes to extract per cell type
    
    Returns:
    --------
    dict : Marker genes per cell type with scores and rankings
    """
    key = f'rank_genes_{de_method}'
    
    if key not in adata.uns:
        raise ValueError(f"rank_genes_groups results not found for key '{key}'")
    
    logger.info(f"Extracting top {top_n} marker genes per cell type...")
    
    results = {}
    cell_types = adata.uns[key]['names'].dtype.names
    
    for celltype in cell_types:
        logger.info(f"Processing cell type: {celltype}")
        
        # Get results for this cell type
        genes = adata.uns[key]['names'][celltype]
        scores = adata.uns[key]['scores'][celltype]
        
        # Some methods (like logreg) don't have pvals, pvals_adj, or logfoldchanges
        pvals = None
        pvals_adj = None
        logfoldchanges = None
        
        if 'pvals' in adata.uns[key]:
            pvals = adata.uns[key]['pvals'][celltype]
        if 'pvals_adj' in adata.uns[key]:
            pvals_adj = adata.uns[key]['pvals_adj'][celltype]
        if 'logfoldchanges' in adata.uns[key]:
            logfoldchanges = adata.uns[key]['logfoldchanges'][celltype]
        
        # Get pts (fraction expressing) if available
        pts = None
        pts_rest = None
        if 'pts' in adata.uns[key]:
            pts = adata.uns[key]['pts'][celltype]
        if 'pts_rest' in adata.uns[key]:
            pts_rest = adata.uns[key]['pts_rest'][celltype]
        
        # Store results
        results[celltype] = {
            'top_genes': list(genes[:top_n]),
            'top_scores': list(scores[:top_n]),
            'gene_ranking': list(genes),  # Full ranking
            'all_scores': {g: s for g, s in zip(genes, scores)},
        }
        
        # Add optional fields if available
        if pvals is not None:
            results[celltype]['top_pvals'] = list(pvals[:top_n])
            results[celltype]['all_pvals'] = {g: p for g, p in zip(genes, pvals)}
        if pvals_adj is not None:
            results[celltype]['top_pvals_adj'] = list(pvals_adj[:top_n])
        if logfoldchanges is not None:
            results[celltype]['top_logfoldchanges'] = list(logfoldchanges[:top_n])
            results[celltype]['all_logfoldchanges'] = {g: lfc for g, lfc in zip(genes, logfoldchanges)}
        if pts is not None:
            results[celltype]['top_pts'] = list(pts[:top_n])
        if pts_rest is not None:
            results[celltype]['top_pts_rest'] = list(pts_rest[:top_n])
        
        logger.info(f"  Top 5 marker genes: {list(genes[:5])}")
    
    return results


def evaluate_literature_markers(marker_results, adata):
    """
    Evaluate how well the identified markers match literature-based markers.
    
    Parameters:
    -----------
    marker_results : dict
        Results from extract_marker_genes
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


def save_results(marker_results, evaluation, method, de_method, output_dir):
    """
    Save marker gene results and evaluation to files.
    
    Parameters:
    -----------
    marker_results : dict
        Results from extract_marker_genes
    evaluation : dict
        Results from evaluate_literature_markers
    method : str
        Sampling method ('full', 'sps', 'random')
    de_method : str
        DE method used
    output_dir : str
        Output directory
    """
    os.makedirs(output_dir, exist_ok=True)
    prefix = f"{method}_{de_method}"
    
    # Save marker results as pickle
    markers_path = os.path.join(output_dir, f'{prefix}_marker_genes.pkl')
    with open(markers_path, 'wb') as f:
        pickle.dump(marker_results, f)
    logger.info(f"Saved marker results to {markers_path}")
    
    # Save marker results as CSV
    marker_df_rows = []
    for celltype, results in marker_results.items():
        for i, gene in enumerate(results['top_genes']):
            row = {
                'cell_type': celltype,
                'rank': i + 1,
                'gene': gene,
                'score': results['top_scores'][i],
            }
            # Add optional fields if available
            if 'top_pvals' in results:
                row['pval'] = results['top_pvals'][i]
            if 'top_pvals_adj' in results:
                row['pval_adj'] = results['top_pvals_adj'][i]
            if 'top_logfoldchanges' in results:
                row['logfoldchange'] = results['top_logfoldchanges'][i]
            if 'top_pts' in results:
                row['pct_in_group'] = results['top_pts'][i]
            if 'top_pts_rest' in results:
                row['pct_in_rest'] = results['top_pts_rest'][i]
            marker_df_rows.append(row)
    
    marker_df = pd.DataFrame(marker_df_rows)
    marker_csv_path = os.path.join(output_dir, f'{prefix}_marker_genes.csv')
    marker_df.to_csv(marker_csv_path, index=False)
    logger.info(f"Saved marker genes CSV to {marker_csv_path}")
    
    # Save evaluation as pickle
    eval_path = os.path.join(output_dir, f'{prefix}_literature_evaluation.pkl')
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
    eval_csv_path = os.path.join(output_dir, f'{prefix}_literature_evaluation.csv')
    eval_df.to_csv(eval_csv_path, index=False)
    logger.info(f"Saved evaluation CSV to {eval_csv_path}")


# ============================================================================
# Main Pipeline
# ============================================================================

def run_pipeline(method, de_methods, output_dir=None):
    """
    Run the full Scanpy marker gene pipeline.
    
    Parameters:
    -----------
    method : str
        Sampling method ('full', 'sps', 'random')
    de_methods : list
        List of DE methods to run
    output_dir : str, optional
        Output directory. If None, uses default OUTPUT_DIR.
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR
    
    logger.info("=" * 80)
    logger.info(f"Scanpy Marker Gene Pipeline - {method.upper()}")
    logger.info(f"DE methods: {de_methods}")
    logger.info("=" * 80)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    logger.info("\n" + "=" * 40)
    logger.info("Step 1: Loading Data")
    logger.info("=" * 40)
    
    adata_full = load_mcc_dataset()
    indices = load_sampling_indices(method)
    
    # Subset to sampled cells
    logger.info("\n" + "=" * 40)
    logger.info("Step 2: Subsetting Data")
    logger.info("=" * 40)
    
    adata = subset_to_sampled_cells(adata_full, indices)
    
    # Run marker gene analysis with each DE method
    all_results = {}
    
    for de_method in de_methods:
        logger.info("\n" + "=" * 40)
        logger.info(f"Step 3: Running {de_method} analysis")
        logger.info("=" * 40)
        
        # Run analysis
        adata = run_marker_gene_analysis(adata, de_method=de_method)
        
        # Extract marker genes
        logger.info("\n" + "=" * 40)
        logger.info(f"Step 4: Extracting marker genes ({de_method})")
        logger.info("=" * 40)
        
        marker_results = extract_marker_genes(adata, de_method=de_method, top_n=TOP_N_GENES)
        
        # Evaluate against literature markers
        logger.info("\n" + "=" * 40)
        logger.info(f"Step 5: Evaluating against literature markers ({de_method})")
        logger.info("=" * 40)
        
        evaluation = evaluate_literature_markers(marker_results, adata)
        
        # Save results
        logger.info("\n" + "=" * 40)
        logger.info(f"Step 6: Saving results ({de_method})")
        logger.info("=" * 40)
        
        save_results(marker_results, evaluation, method, de_method, output_dir)
        
        all_results[de_method] = {
            'marker_results': marker_results,
            'evaluation': evaluation,
        }
        
        # Print summary for this method
        logger.info(f"\n{de_method.upper()} Summary:")
        for celltype, metrics in evaluation.items():
            logger.info(f"  {celltype}:")
            logger.info(f"    Recall@10: {metrics.get('recall_at_10', 'N/A')}")
            logger.info(f"    Recall@50: {metrics.get('recall_at_50', 'N/A')}")
            logger.info(f"    Mean rank: {metrics.get('mean_rank', 'N/A')}")
    
    # Print final summary
    logger.info("\n" + "=" * 80)
    logger.info("FINAL SUMMARY")
    logger.info("=" * 80)
    
    for de_method, results in all_results.items():
        logger.info(f"\n{de_method.upper()}:")
        for celltype, metrics in results['evaluation'].items():
            recall_10 = metrics.get('recall_at_10')
            recall_50 = metrics.get('recall_at_50')
            mean_rank = metrics.get('mean_rank')
            logger.info(f"  {celltype}: Recall@10={recall_10:.2%} if recall_10 else 'N/A', "
                       f"Recall@50={recall_50:.2%} if recall_50 else 'N/A', "
                       f"MeanRank={mean_rank:.1f} if mean_rank else 'N/A'")
    
    logger.info("\n" + "=" * 80)
    logger.info("Pipeline completed successfully!")
    logger.info(f"Results saved to: {output_dir}")
    logger.info("=" * 80)
    
    return all_results


# ============================================================================
# Entry Point
# ============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Scanpy Marker Gene Pipeline for MCC Dataset'
    )
    parser.add_argument(
        '--method', 
        type=str, 
        required=True,
        choices=['full', 'sps', 'random'],
        help='Sampling method (full, sps, or random)'
    )
    parser.add_argument(
        '--de-methods',
        type=str,
        nargs='+',
        default=['wilcoxon'],
        choices=AVAILABLE_DE_METHODS,
        help=f'DE methods to use (default: wilcoxon). Options: {AVAILABLE_DE_METHODS}'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory (default: analysis/scanpy_marker_genes/results)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logger(args.method)
    
    # Run pipeline
    try:
        results = run_pipeline(
            method=args.method,
            de_methods=args.de_methods,
            output_dir=args.output_dir
        )
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
