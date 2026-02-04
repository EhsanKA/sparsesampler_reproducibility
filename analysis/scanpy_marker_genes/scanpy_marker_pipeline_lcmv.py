#!/usr/bin/env python
# coding: utf-8

"""
Scanpy Marker Gene Pipeline for LCMV Dataset

This script identifies marker genes for each cell type in the LCMV dataset using Scanpy's
rank_genes_groups, comparing results when using Full dataset vs SPS 100k vs Random 100k.

The pipeline:
1. Loads the full LCMV dataset and sampling indices
2. Subsets to sampled cells (or uses full dataset)
3. Runs Scanpy rank_genes_groups with multiple DE methods
4. Identifies marker genes per cell type
5. Saves results for downstream analysis

Usage:
    conda activate facs_sampling
    python scanpy_marker_pipeline_lcmv.py --method full --de-methods wilcoxon t-test logreg
    python scanpy_marker_pipeline_lcmv.py --method sps --de-methods wilcoxon t-test logreg
    python scanpy_marker_pipeline_lcmv.py --method random --de-methods wilcoxon t-test logreg

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

# ============================================================================
# Configuration
# ============================================================================

# Paths
PROJECT_ROOT = os.environ.get(
    'PROJECT_ROOT',
    '/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility'
)
# LCMV data is in the sparseFlow_benchmarking directory
DATA_PATH_LCMV = '/fast/AG_Ohler/ekarimi/projects/sparseFlow_benchmarking/data/lcmv/benchmark/34'

script_dir = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(script_dir, 'results_lcmv')

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
    log_file = os.path.join(log_dir, f'scanpy_marker_gene_lcmv_{method}_{timestamp}.log')
    
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

def load_lcmv_dataset():
    """Load the full LCMV dataset."""
    adata_path = os.path.join(DATA_PATH_LCMV, 'adata.h5ad')
    logger.info(f"Loading LCMV dataset from {adata_path}")
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
    
    base_path = os.path.join(DATA_PATH_LCMV, method, str(SAMPLE_SIZE), str(REP))
    
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


def subset_to_sampled_cells(adata, indices, max_cells=None):
    """Subset AnnData to sampled cells.
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    indices : np.ndarray or None
        Indices to subset to. If None and max_cells is set, random subsample.
    max_cells : int or None
        Maximum number of cells. If set and indices is None, randomly subsample
        to this many cells (stratified by cell type).
    """
    if indices is None and max_cells is None:
        logger.info("Using full dataset (no subsetting)")
        return adata
    
    if indices is None and max_cells is not None:
        # Stratified random subsample to max_cells
        logger.info(f"Subsampling full dataset to {max_cells:,} cells (stratified)...")
        
        from sklearn.model_selection import train_test_split
        
        n_total = adata.n_obs
        if max_cells >= n_total:
            logger.info(f"max_cells ({max_cells:,}) >= total cells ({n_total:,}), using full dataset")
            return adata
        
        # Calculate sampling fraction
        frac = max_cells / n_total
        
        # Stratified subsample
        all_indices = np.arange(n_total)
        labels = adata.obs[LABEL_KEY].values
        
        # Use train_test_split to get stratified subsample
        _, subsample_idx = train_test_split(
            all_indices,
            test_size=frac,
            stratify=labels,
            random_state=42
        )
        
        logger.info(f"Selected {len(subsample_idx):,} cells via stratified sampling")
        adata_sub = adata[subsample_idx].copy()
    else:
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

def run_marker_gene_analysis(adata, de_method='wilcoxon', use_pts=True):
    """
    Run Scanpy rank_genes_groups for marker gene identification.
    
    Parameters:
    -----------
    adata : AnnData
        Dataset to analyze
    de_method : str
        Differential expression method ('wilcoxon', 't-test', 't-test_overestim_var', 'logreg')
    use_pts : bool
        Whether to compute fraction of cells expressing each gene.
        Set to False for very large datasets (>10M cells) to avoid chunking issues.
    
    Returns:
    --------
    AnnData : adata with rank_genes_groups results stored
    """
    logger.info(f"Running marker gene analysis with method: {de_method}, pts={use_pts}")
    
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
            pts=use_pts,  # Include fraction of cells expressing (disable for large datasets)
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


def save_results(marker_results, method, de_method, output_dir):
    """
    Save marker gene results to files.
    
    Parameters:
    -----------
    marker_results : dict
        Results from extract_marker_genes
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


# ============================================================================
# Main Pipeline
# ============================================================================

def run_pipeline(method, de_methods, output_dir=None, max_cells=None):
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
    max_cells : int, optional
        Maximum number of cells for full dataset. If set and method='full',
        will randomly subsample to this many cells (stratified by cell type).
        Useful for avoiding scanpy memory issues with very large datasets.
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR
    
    logger.info("=" * 80)
    logger.info(f"Scanpy Marker Gene Pipeline (LCMV) - {method.upper()}")
    logger.info(f"DE methods: {de_methods}")
    if max_cells and method == 'full':
        logger.info(f"Max cells: {max_cells:,} (subsampling full dataset)")
    logger.info("=" * 80)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    logger.info("\n" + "=" * 40)
    logger.info("Step 1: Loading Data")
    logger.info("=" * 40)
    
    adata_full = load_lcmv_dataset()
    indices = load_sampling_indices(method)
    
    # Subset to sampled cells
    logger.info("\n" + "=" * 40)
    logger.info("Step 2: Subsetting Data")
    logger.info("=" * 40)
    
    # For full dataset, apply max_cells limit if specified
    effective_max_cells = max_cells if method == 'full' else None
    adata = subset_to_sampled_cells(adata_full, indices, max_cells=effective_max_cells)
    
    # Free memory
    del adata_full
    import gc
    gc.collect()
    
    # Run marker gene analysis with each DE method
    all_results = {}
    
    for de_method in de_methods:
        logger.info("\n" + "=" * 40)
        logger.info(f"Step 3: Running {de_method} analysis")
        logger.info("=" * 40)
        
        # For very large datasets (full LCMV ~34M cells), disable pts calculation
        # to avoid scanpy's chunking bug with wilcoxon test
        use_pts = True
        if method == 'full' and adata.n_obs > 10_000_000:
            logger.info(f"Large dataset detected ({adata.n_obs:,} cells), disabling pts calculation")
            use_pts = False
        
        # Run analysis
        adata = run_marker_gene_analysis(adata, de_method=de_method, use_pts=use_pts)
        
        # Extract marker genes
        logger.info("\n" + "=" * 40)
        logger.info(f"Step 4: Extracting marker genes ({de_method})")
        logger.info("=" * 40)
        
        marker_results = extract_marker_genes(adata, de_method=de_method, top_n=TOP_N_GENES)
        
        # Save results
        logger.info("\n" + "=" * 40)
        logger.info(f"Step 5: Saving results ({de_method})")
        logger.info("=" * 40)
        
        save_results(marker_results, method, de_method, output_dir)
        
        all_results[de_method] = {
            'marker_results': marker_results,
        }
        
        # Print summary for this method
        logger.info(f"\n{de_method.upper()} Summary:")
        for celltype in marker_results.keys():
            top_genes = marker_results[celltype]['top_genes'][:5]
            logger.info(f"  {celltype}: {top_genes}")
    
    # Print final summary
    logger.info("\n" + "=" * 80)
    logger.info("FINAL SUMMARY")
    logger.info("=" * 80)
    
    logger.info(f"\nAnalyzed {len(marker_results)} cell types with {len(de_methods)} DE methods")
    
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
        description='Scanpy Marker Gene Pipeline for LCMV Dataset'
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
        help='Output directory (default: analysis/scanpy_marker_genes/results_lcmv)'
    )
    parser.add_argument(
        '--max-cells',
        type=int,
        default=None,
        help='Maximum cells for full dataset (stratified subsample). '
             'Use to avoid scanpy memory issues with very large datasets (e.g., --max-cells 20000000)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logger(args.method)
    
    # Run pipeline
    try:
        results = run_pipeline(
            method=args.method,
            de_methods=args.de_methods,
            output_dir=args.output_dir,
            max_cells=args.max_cells
        )
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
