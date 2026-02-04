#!/usr/bin/env python
# coding: utf-8

"""
SEMITONES Analysis using Sampled Data as References

This script uses SEMITONES (Similarity Enrichment for Multiple Interrogation Techniques 
in Omics via Neighborhoods of Exemplars) to evaluate how well sampled cells represent 
the full dataset. The sampled cells are used as references, and similarity scores are 
calculated for all cells in the dataset.

Usage:
    conda activate semitones_env
    python semitones_analysis.py --dataset mcc_01 --ref 30 --method sps --size 200000 --rep 0
"""

import os
import sys
import numpy as np
import pandas as pd
import pickle
import logging
from datetime import datetime
import argparse
import scanpy as sc
import scipy.sparse

# Add analysis directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

from refined_rare_cell_type_definition import (
    load_reference_data as load_ref_data
)

# Try to import SEMITONES
try:
    from SEMITONES import enrichment_scoring
    from SEMITONES import cell_selection
    from SEMITONES.support_funcs import pairwise_similarities
    SEMITONES_AVAILABLE = True
except ImportError as e:
    SEMITONES_AVAILABLE = False
    print(f"Warning: SEMITONES not available: {e}")
    print("Please activate the semitones_env conda environment:")
    print("  conda activate semitones_env")

# Set up logging
def setup_logger():
    log_dir = os.path.join(os.path.dirname(script_dir), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'semitones_analysis_{timestamp}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

# Constants
METHODS = ['random', 'sps', 'hopper', 'atomic', 'scsampler']
label_key = 'celltype'

# Dataset configurations
DATASET_CONFIGS = {
    'mcc_01': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
    },
    'mcc_05': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
    },
    'lcmv': {
        'references': [1, 5, 10, 20, 34],
        'sizes': [50000, 100000, 200000],
    }
}

REPS = [0]  # Using rep 0 for now

def get_data_path():
    """Get data path from environment variable or derive from script location."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
    return os.path.join(project_root, 'data')

def get_output_dir():
    """Get output directory for results."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, 'results', 'semitones_analysis')
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def load_sampling_indices(dataset, ref, method, size, rep, data_path):
    """Load sampling indices for a given configuration."""
    directory = f"{dataset}/benchmark"
    path = os.path.join(data_path, directory)
    
    try:
        if method == 'atomic':
            atomic_address = os.path.join(path, f"{ref}/atomic/{size}/{rep}/results.csv")
            if os.path.isfile(atomic_address):
                if dataset in ['mcc_01', 'mcc_05']:
                    indices = pd.read_csv(atomic_address)['x'].values.astype(str)
                else:
                    indices = pd.read_csv(atomic_address)['x'].values.astype(int)
                return indices
            else:
                return None
        else:
            samples_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/results.pkl")
            if os.path.isfile(samples_address):
                with open(samples_address, 'rb') as handle:
                    samples = pickle.load(handle)
                return samples[0]
            else:
                return None
    except Exception as e:
        logger.warning(f"Error loading indices for {dataset}, ref {ref}, method {method}, size {size}, rep {rep}: {str(e)}")
        return None

def prepare_data_for_semitones(adata, reference_indices):
    """
    Prepare data for SEMITONES analysis.
    
    Parameters:
    -----------
    adata : sc.AnnData
        Full AnnData object
    reference_indices : np.ndarray
        Indices of reference cells (sampled cells)
    
    Returns:
    --------
    tuple: (reference_adata, query_adata, reference_idx, query_idx)
    """
    # Handle string vs integer indices
    if isinstance(reference_indices[0], (str, np.str_)):
        # String indices - use obs.index
        ref_mask = adata.obs.index.isin(reference_indices)
        reference_idx = np.where(ref_mask)[0]
    else:
        # Integer indices - assume positional
        reference_idx = reference_indices
    
    # Get query indices (all cells not in reference)
    query_idx = np.setdiff1d(np.arange(adata.shape[0]), reference_idx)
    
    # Create reference and query AnnData objects
    reference_adata = adata[reference_idx].copy()
    query_adata = adata[query_idx].copy()
    
    logger.info(f"Reference cells: {reference_adata.shape[0]}")
    logger.info(f"Query cells: {query_adata.shape[0]}")
    
    return reference_adata, query_adata, reference_idx, query_idx

def run_semitones_analysis(adata_full, reference_indices, n_neighbors=15, metric='cosine'):
    """
    Run SEMITONES analysis using sampled cells as references.
    
    Parameters:
    -----------
    adata_full : sc.AnnData
        Full dataset
    reference_indices : np.ndarray
        Indices of reference cells (sampled cells)
    n_neighbors : int
        Number of neighbors for similarity calculation
    metric : str
        Distance metric to use
    
    Returns:
    --------
    dict: Results dictionary with similarity scores and metrics
    """
    if not SEMITONES_AVAILABLE:
        raise ImportError("SEMITONES is not available. Please activate semitones_env environment.")
    
    logger.info("Preparing data for SEMITONES...")
    reference_adata, query_adata, ref_idx, query_idx = prepare_data_for_semitones(
        adata_full, reference_indices
    )
    
    # Ensure data is in the right format (dense or sparse)
    # SEMITONES typically works with dense arrays or sparse matrices
    X_ref = reference_adata.X
    X_query = query_adata.X
    
    if scipy.sparse.issparse(X_ref):
        X_ref = X_ref.toarray()
    if scipy.sparse.issparse(X_query):
        X_query = X_query.toarray()
    
    logger.info("Calculating enrichment scores using SEMITONES...")
    # SEMITONES calculate_escores takes reference data (X) and query data
    # It calculates similarity and enrichment scores internally
    enrichment_scores = enrichment_scoring.calculate_escores(
        X=X_ref,  # Reference data (sampled cells)
        query=X_query,  # Query data (all other cells)
        metric=metric,
        ncpu=1  # Number of CPUs to use
    )
    
    logger.info(f"Enrichment scores shape: {enrichment_scores.shape}")
    
    # Calculate pairwise similarities separately if needed
    logger.info("Calculating pairwise similarities...")
    similarity_scores = pairwise_similarities(
        X_query, 
        X_ref, 
        metric=metric
    )
    
    logger.info(f"Similarity scores shape: {similarity_scores.shape}")
    
    # Calculate summary statistics
    results = {
        'similarity_scores': similarity_scores,
        'enrichment_scores': enrichment_scores,
        'reference_indices': ref_idx,
        'query_indices': query_idx,
        'n_reference': len(ref_idx),
        'n_query': len(query_idx),
        'mean_similarity': np.mean(similarity_scores),
        'median_similarity': np.median(similarity_scores),
        'max_similarity': np.max(similarity_scores),
        'min_similarity': np.min(similarity_scores),
    }
    
    return results

def analyze_cell_type_coverage(adata_full, reference_indices, similarity_results, label_key='celltype'):
    """
    Analyze how well reference cells cover different cell types.
    
    Parameters:
    -----------
    adata_full : sc.AnnData
        Full dataset
    reference_indices : np.ndarray
        Indices of reference cells
    similarity_results : dict
        Results from run_semitones_analysis
    label_key : str
        Key for cell type labels
    
    Returns:
    --------
    dict: Cell type coverage metrics
    """
    # Get cell types
    all_cell_types = adata_full.obs[label_key].values
    
    # Handle string vs integer indices
    if isinstance(reference_indices[0], (str, np.str_)):
        ref_mask = adata_full.obs.index.isin(reference_indices)
    else:
        ref_mask = np.zeros(len(adata_full), dtype=bool)
        ref_mask[reference_indices] = True
    
    ref_cell_types = all_cell_types[ref_mask]
    query_cell_types = all_cell_types[~ref_mask]
    
    # Calculate coverage per cell type
    unique_types = np.unique(all_cell_types)
    coverage_metrics = {}
    
    for cell_type in unique_types:
        total_count = np.sum(all_cell_types == cell_type)
        ref_count = np.sum(ref_cell_types == cell_type)
        query_count = np.sum(query_cell_types == cell_type)
        
        coverage_metrics[cell_type] = {
            'total': total_count,
            'in_reference': ref_count,
            'in_query': query_count,
            'reference_fraction': ref_count / total_count if total_count > 0 else 0,
            'query_fraction': query_count / total_count if total_count > 0 else 0,
        }
    
    # Calculate average similarity per cell type
    similarity_scores = similarity_results['similarity_scores']
    query_idx = similarity_results['query_indices']
    query_cell_types_array = query_cell_types
    
    type_similarities = {}
    for cell_type in unique_types:
        type_mask = query_cell_types_array == cell_type
        if np.any(type_mask):
            type_similarities[cell_type] = {
                'mean_similarity': np.mean(similarity_scores[type_mask, :]),
                'median_similarity': np.median(similarity_scores[type_mask, :]),
                'max_similarity': np.max(similarity_scores[type_mask, :]),
            }
        else:
            type_similarities[cell_type] = {
                'mean_similarity': np.nan,
                'median_similarity': np.nan,
                'max_similarity': np.nan,
            }
    
    return {
        'coverage_metrics': coverage_metrics,
        'type_similarities': type_similarities,
    }

def run_single_semitones_analysis(dataset, ref, method, size, rep, data_path, output_dir):
    """
    Run SEMITONES analysis for a single configuration.
    
    Parameters:
    -----------
    dataset : str
        Dataset name
    ref : int
        Reference size
    method : str
        Sampling method
    size : int
        Sample size
    rep : int
        Replication number
    data_path : str
        Path to data directory
    output_dir : str
        Path to output directory
    
    Returns:
    --------
    dict: Analysis results
    """
    logger.info(f"\n{'='*80}")
    logger.info(f"SEMITONES Analysis: {dataset}, ref={ref}, method={method}, size={size}, rep={rep}")
    logger.info(f"{'='*80}")
    
    # Load full dataset
    logger.info(f"Loading full dataset: {dataset}, ref {ref}")
    adata_full = load_ref_data(dataset, ref, base_path=data_path)
    logger.info(f"Full dataset shape: {adata_full.shape}")
    
    # Load sampling indices
    logger.info(f"Loading sampling indices for method={method}, size={size}, rep={rep}")
    indices = load_sampling_indices(dataset, ref, method, size, rep, data_path)
    
    if indices is None:
        logger.warning(f"No indices found for {dataset}, ref {ref}, method {method}, size {size}, rep {rep}")
        return None
    
    logger.info(f"Loaded {len(indices)} reference cell indices")
    
    # Run SEMITONES analysis
    try:
        similarity_results = run_semitones_analysis(
            adata_full, 
            indices,
            n_neighbors=15,
            metric='cosine'
        )
        
        # Analyze cell type coverage
        coverage_results = analyze_cell_type_coverage(
            adata_full,
            indices,
            similarity_results,
            label_key=label_key
        )
        
        # Combine results
        results = {
            'dataset': dataset,
            'ref': ref,
            'method': method,
            'sample_size': size,
            'rep': rep,
            'n_reference': similarity_results['n_reference'],
            'n_query': similarity_results['n_query'],
            'mean_similarity': similarity_results['mean_similarity'],
            'median_similarity': similarity_results['median_similarity'],
            'max_similarity': similarity_results['max_similarity'],
            'min_similarity': similarity_results['min_similarity'],
            'coverage_metrics': coverage_results['coverage_metrics'],
            'type_similarities': coverage_results['type_similarities'],
        }
        
        # Save results
        output_file = os.path.join(
            output_dir,
            f"{dataset}_ref{ref}_{method}_size{size}_rep{rep}_semitones_results.pkl"
        )
        with open(output_file, 'wb') as f:
            pickle.dump(results, f)
        logger.info(f"Saved results to {output_file}")
        
        return results
        
    except Exception as e:
        logger.error(f"Error running SEMITONES analysis: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return None

def main():
    parser = argparse.ArgumentParser(
        description='Run SEMITONES analysis using sampled cells as references'
    )
    parser.add_argument('--dataset', type=str, required=True,
                       choices=['mcc_01', 'mcc_05', 'lcmv'],
                       help='Dataset name')
    parser.add_argument('--ref', type=int, required=True,
                       help='Reference size')
    parser.add_argument('--method', type=str, required=True,
                       choices=METHODS,
                       help='Sampling method')
    parser.add_argument('--size', type=int, required=True,
                       help='Sample size')
    parser.add_argument('--rep', type=int, default=0,
                       help='Replication number (default: 0)')
    
    args = parser.parse_args()
    
    if not SEMITONES_AVAILABLE:
        logger.error("SEMITONES is not available. Please activate semitones_env environment:")
        logger.error("  conda activate semitones_env")
        sys.exit(1)
    
    # Get paths
    data_path = get_data_path()
    output_dir = get_output_dir()
    
    # Run analysis
    results = run_single_semitones_analysis(
        args.dataset,
        args.ref,
        args.method,
        args.size,
        args.rep,
        data_path,
        output_dir
    )
    
    if results is None:
        logger.error("Analysis failed")
        sys.exit(1)
    
    logger.info("Analysis completed successfully")

if __name__ == '__main__':
    main()

