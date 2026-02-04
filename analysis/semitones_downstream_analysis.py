#!/usr/bin/env python
# coding: utf-8

"""
SEMITONES Downstream Analysis - Following the SEMITONES Paper Methodology

This script implements the downstream analysis workflow suggested in the SEMITONES paper:
1. Use SPS samples as reference cells
2. Calculate enrichment scores for all query cells
3. Analyze enrichment by cell type
4. Identify rare/poorly represented cells
5. Create visualizations and metrics

Usage:
    conda activate semitones_env
    python semitones_downstream_analysis.py --dataset mcc_01 --ref 30 --method sps --size 200000 --rep 0
"""

import os
import sys
import numpy as np
import pandas as pd
import pickle
import logging
from datetime import datetime
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
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
    from SEMITONES.enrichment_scoring import calculate_escores
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
    log_file = os.path.join(log_dir, f'semitones_downstream_{timestamp}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

label_key = 'celltype'

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
    output_dir = os.path.join(script_dir, 'results', 'semitones_downstream')
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
            # First try to load from .npy file (more compatible)
            npy_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/indices.npy")
            if os.path.isfile(npy_address):
                try:
                    import numpy as np
                    indices = np.load(npy_address)
                    logger.info(f"Loaded indices from .npy file: {len(indices)} indices")
                    return indices
                except Exception as e:
                    logger.warning(f"Failed to load .npy file: {e}, trying pickle...")
            
            # Fall back to pickle file
            samples_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/results.pkl")
            if os.path.isfile(samples_address):
                try:
                    with open(samples_address, 'rb') as handle:
                        samples = pickle.load(handle)
                    return samples[0]
                except (ModuleNotFoundError, AttributeError, ImportError, ValueError) as e:
                    # Handle numpy compatibility issues with pickle files
                    # Try loading with different numpy version handling
                    logger.warning(f"Pickle loading error (numpy compatibility): {str(e)}")
                    logger.info("Attempting to load with encoding='latin1' for compatibility...")
                    try:
                        # Try with encoding='latin1' for Python 2/3 compatibility
                        with open(samples_address, 'rb') as handle:
                            samples = pickle.load(handle, encoding='latin1')
                        return samples[0]
                    except Exception as e2:
                        logger.info("Trying with errors='ignore'...")
                        try:
                            # Try with errors='ignore' for numpy compatibility
                            import sys
                            if sys.version_info >= (3, 7):
                                with open(samples_address, 'rb') as handle:
                                    samples = pickle.load(handle, encoding='latin1', errors='ignore')
                            else:
                                with open(samples_address, 'rb') as handle:
                                    samples = pickle.load(handle, encoding='latin1')
                            return samples[0]
                        except Exception as e3:
                            # Try using pickle5 if available
                            try:
                                import pickle5
                                with open(samples_address, 'rb') as handle:
                                    samples = pickle5.load(handle)
                                return samples[0]
                            except (ImportError, Exception) as e4:
                                # Last resort: try with numpy compatibility workaround
                                logger.info("Trying numpy compatibility workaround...")
                                try:
                                    import numpy as np
                                    # Temporarily patch numpy if needed
                                    old_getattr = None
                                    if not hasattr(np, '_core'):
                                        # Workaround for numpy version mismatch
                                        class NumpyCompat:
                                            pass
                                        np._core = NumpyCompat()
                                    
                                    with open(samples_address, 'rb') as handle:
                                        samples = pickle.load(handle)
                                    return samples[0]
                                except Exception as e5:
                                    logger.error(f"All pickle loading methods failed. Last error: {str(e5)}")
                                    logger.error("You may need to regenerate the pickle files or load them in the original environment")
                                    return None
            else:
                return None
    except Exception as e:
        logger.warning(f"Error loading indices: {str(e)}")
        return None

def prepare_reference_and_query(adata_full, reference_indices):
    """
    Split dataset into reference (SPS samples) and query (rest) cells.
    
    Returns:
    --------
    reference_adata, query_adata, ref_idx, query_idx
    """
    if isinstance(reference_indices[0], (str, np.str_)):
        ref_mask = adata_full.obs.index.isin(reference_indices)
        ref_idx = np.where(ref_mask)[0]
    else:
        ref_idx = reference_indices
    
    query_idx = np.setdiff1d(np.arange(adata_full.shape[0]), ref_idx)
    
    reference_adata = adata_full[ref_idx].copy()
    query_adata = adata_full[query_idx].copy()
    
    return reference_adata, query_adata, ref_idx, query_idx

def calculate_enrichment_scores(reference_adata, query_adata, metric='cosine', batch_size=20000):
    """
    Calculate SEMITONES enrichment scores using the proper SEMITONES calculate_escores function.
    
    Processes in batches to avoid memory issues with large datasets.
    
    Parameters:
    -----------
    batch_size : int
        Number of query cells to process at once (for batching if needed)
    """
    n_query = query_adata.shape[0]
    n_ref = reference_adata.shape[0]
    
    logger.info(f"Calculating SEMITONES enrichment scores...")
    logger.info(f"Query cells: {n_query}, Reference cells: {n_ref}")
    logger.info(f"Using metric: {metric}")
    
    # Convert reference to dense (smaller, can fit in memory)
    X_ref = reference_adata.X
    if scipy.sparse.issparse(X_ref):
        logger.info("Converting reference data from sparse to dense...")
        X_ref = X_ref.toarray()
    
    # For large datasets, use batch processing by default
    # Threshold: if query cells > 500k, use batching
    use_batching = n_query > 500000
    
    if use_batching:
        logger.info(f"Large dataset detected ({n_query} query cells), using batch processing (batch_size={batch_size})...")
        X_query_sparse = query_adata.X
        n_batches = (n_query + batch_size - 1) // batch_size
        logger.info(f"Processing {n_batches} batches...")
        
        enrichment_scores = np.zeros(n_query)
        
        for i in range(n_batches):
            start_idx = i * batch_size
            end_idx = min((i + 1) * batch_size, n_query)
            
            logger.info(f"Processing batch {i+1}/{n_batches} (cells {start_idx}-{end_idx})...")
            
            # Convert batch to dense
            if scipy.sparse.issparse(X_query_sparse):
                X_query_batch = X_query_sparse[start_idx:end_idx].toarray()
            else:
                X_query_batch = X_query_sparse[start_idx:end_idx]
            
            # Calculate pairwise similarities for this batch
            try:
                # Use sklearn's cosine_similarity directly (same as SEMITONES uses internally)
                from sklearn.metrics.pairwise import cosine_similarity
                similarity_batch = cosine_similarity(X_query_batch, X_ref)
                
                # Enrichment score is the mean similarity to all reference cells
                batch_scores = np.mean(similarity_batch, axis=1)
                enrichment_scores[start_idx:end_idx] = batch_scores
                
            except Exception as e:
                logger.error(f"Error processing batch {i+1}: {e}")
                raise
            
            # Free memory
            del X_query_batch, batch_scores
        
    else:
        # For smaller datasets, try processing all at once
        logger.info("Attempting to calculate enrichment scores for all query cells at once...")
        
        X_query = query_adata.X
        if scipy.sparse.issparse(X_query):
            logger.info("Converting query data from sparse to dense...")
            X_query = X_query.toarray()
        
        try:
            # Use sklearn's cosine_similarity directly (same as SEMITONES uses internally)
            from sklearn.metrics.pairwise import cosine_similarity
            similarity_matrix = cosine_similarity(X_query, X_ref)
            
            # Enrichment score is the mean similarity to all reference cells
            enrichment_scores = np.mean(similarity_matrix, axis=1)
            logger.info(f"Enrichment scores shape: {enrichment_scores.shape}")
            logger.info("Successfully calculated enrichment scores using cosine similarity")
            
        except MemoryError:
            logger.warning("Memory error with full dataset, falling back to batch processing...")
            # Fall back to batch processing
            X_query_sparse = query_adata.X
            n_batches = (n_query + batch_size - 1) // batch_size
            logger.info(f"Processing {n_batches} batches...")
            
            enrichment_scores = np.zeros(n_query)
            
            for i in range(n_batches):
                start_idx = i * batch_size
                end_idx = min((i + 1) * batch_size, n_query)
                
                logger.info(f"Processing batch {i+1}/{n_batches} (cells {start_idx}-{end_idx})...")
                
                # Convert batch to dense
                if scipy.sparse.issparse(X_query_sparse):
                    X_query_batch = X_query_sparse[start_idx:end_idx].toarray()
                else:
                    X_query_batch = X_query_sparse[start_idx:end_idx]
                
                # Calculate enrichment scores for this batch using cosine similarity
                from sklearn.metrics.pairwise import cosine_similarity
                similarity_batch = cosine_similarity(X_query_batch, X_ref)
                batch_scores = np.mean(similarity_batch, axis=1)
                
                enrichment_scores[start_idx:end_idx] = batch_scores
                
                # Free memory
                del X_query_batch, batch_scores
    
    # For compatibility, return None for similarity_matrix (not needed)
    similarity_matrix = None
    
    logger.info("Finished calculating enrichment scores")
    logger.info(f"Enrichment scores stats: mean={np.mean(enrichment_scores):.4f}, "
               f"median={np.median(enrichment_scores):.4f}, "
               f"min={np.min(enrichment_scores):.4f}, max={np.max(enrichment_scores):.4f}")
    
    return similarity_matrix, enrichment_scores

def analyze_enrichment_by_cell_type(query_adata, enrichment_scores, label_key='celltype'):
    """Analyze enrichment scores by cell type."""
    cell_types = query_adata.obs[label_key].values
    unique_types = np.unique(cell_types)
    
    results = {}
    for cell_type in unique_types:
        type_mask = cell_types == cell_type
        type_scores = enrichment_scores[type_mask]
        
        results[cell_type] = {
            'mean_enrichment': np.mean(type_scores),
            'median_enrichment': np.median(type_scores),
            'std_enrichment': np.std(type_scores),
            'min_enrichment': np.min(type_scores),
            'max_enrichment': np.max(type_scores),
            'n_cells': np.sum(type_mask),
            'pct_high_enrichment': np.sum(type_scores > np.percentile(enrichment_scores, 75)) / len(type_scores) * 100
        }
    
    return results

def identify_poorly_represented_cells(enrichment_scores, threshold_percentile=5):
    """Identify cells with low enrichment scores (poorly represented)."""
    threshold = np.percentile(enrichment_scores, threshold_percentile)
    poorly_represented = enrichment_scores < threshold
    return poorly_represented, threshold

def create_enrichment_visualizations(enrichment_scores, query_adata, cell_type_enrichment, 
                                    output_dir, dataset, method, size):
    """Create visualizations as suggested in SEMITONES paper."""
    
    # 1. Overall enrichment distribution
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Histogram of enrichment scores
    axes[0, 0].hist(enrichment_scores, bins=50, edgecolor='black', alpha=0.7)
    axes[0, 0].axvline(np.mean(enrichment_scores), color='red', linestyle='--', 
                      label=f'Mean: {np.mean(enrichment_scores):.3f}')
    axes[0, 0].axvline(np.median(enrichment_scores), color='blue', linestyle='--', 
                      label=f'Median: {np.median(enrichment_scores):.3f}')
    axes[0, 0].set_xlabel('Enrichment Score')
    axes[0, 0].set_ylabel('Number of Cells')
    axes[0, 0].set_title('Distribution of Enrichment Scores')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Boxplot by cell type (top 10 types)
    cell_types = query_adata.obs[label_key].values
    type_df = pd.DataFrame({
        'cell_type': cell_types,
        'enrichment': enrichment_scores
    })
    top_types = type_df.groupby('cell_type').size().nlargest(10).index
    type_df_filtered = type_df[type_df['cell_type'].isin(top_types)]
    
    sns.boxplot(data=type_df_filtered, x='cell_type', y='enrichment', ax=axes[0, 1])
    axes[0, 1].set_xlabel('Cell Type')
    axes[0, 1].set_ylabel('Enrichment Score')
    axes[0, 1].set_title('Enrichment Scores by Cell Type (Top 10)')
    axes[0, 1].tick_params(axis='x', rotation=45)
    axes[0, 1].grid(True, alpha=0.3)
    
    # Mean enrichment per cell type
    type_means = [(ct, cell_type_enrichment[ct]['mean_enrichment']) 
                  for ct in sorted(cell_type_enrichment.keys())]
    types, means = zip(*sorted(type_means, key=lambda x: x[1], reverse=True)[:15])
    
    axes[1, 0].barh(range(len(types)), means)
    axes[1, 0].set_yticks(range(len(types)))
    axes[1, 0].set_yticklabels(types)
    axes[1, 0].set_xlabel('Mean Enrichment Score')
    axes[1, 0].set_title('Mean Enrichment by Cell Type (Top 15)')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Cumulative distribution
    sorted_scores = np.sort(enrichment_scores)
    axes[1, 1].plot(sorted_scores, np.arange(len(sorted_scores)) / len(sorted_scores), 
                   linewidth=2)
    axes[1, 1].axvline(np.percentile(enrichment_scores, 5), color='red', 
                       linestyle='--', label='5th percentile')
    axes[1, 1].axvline(np.percentile(enrichment_scores, 95), color='blue', 
                       linestyle='--', label='95th percentile')
    axes[1, 1].set_xlabel('Enrichment Score')
    axes[1, 1].set_ylabel('Cumulative Fraction')
    axes[1, 1].set_title('Cumulative Distribution of Enrichment Scores')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 
                               f"{dataset}_{method}_size{size}_enrichment_analysis.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved enrichment visualization to {output_file}")

def run_downstream_analysis(dataset, ref, method, size, rep, data_path, output_dir):
    """Run complete SEMITONES downstream analysis."""
    
    logger.info(f"\n{'='*80}")
    logger.info(f"SEMITONES Downstream Analysis")
    logger.info(f"Dataset: {dataset}, Ref: {ref}, Method: {method}, Size: {size}, Rep: {rep}")
    logger.info(f"{'='*80}")
    
    # Load data
    logger.info("Loading full dataset...")
    adata_full = load_ref_data(dataset, ref, base_path=data_path)
    logger.info(f"Full dataset: {adata_full.shape[0]} cells, {adata_full.shape[1]} genes")
    
    # Load SPS sample indices (these become references)
    logger.info(f"Loading {method} sample indices...")
    reference_indices = load_sampling_indices(dataset, ref, method, size, rep, data_path)
    
    if reference_indices is None:
        logger.error("Failed to load sample indices")
        return None
    
    logger.info(f"Loaded {len(reference_indices)} reference cell indices")
    
    # Split into reference and query
    logger.info("Preparing reference and query cells...")
    reference_adata, query_adata, ref_idx, query_idx = prepare_reference_and_query(
        adata_full, reference_indices
    )
    
    logger.info(f"Reference cells: {reference_adata.shape[0]}")
    logger.info(f"Query cells: {query_adata.shape[0]}")
    
    # Calculate enrichment scores
    logger.info("Calculating enrichment scores...")
    similarity_matrix, enrichment_scores = calculate_enrichment_scores(
        reference_adata, query_adata
    )
    
    logger.info(f"Enrichment scores calculated: mean={np.mean(enrichment_scores):.4f}, "
               f"median={np.median(enrichment_scores):.4f}")
    
    # Analyze by cell type
    logger.info("Analyzing enrichment by cell type...")
    cell_type_enrichment = analyze_enrichment_by_cell_type(
        query_adata, enrichment_scores, label_key
    )
    
    # Identify poorly represented cells
    poorly_represented, threshold = identify_poorly_represented_cells(enrichment_scores)
    logger.info(f"Identified {np.sum(poorly_represented)} poorly represented cells "
               f"(threshold: {threshold:.4f})")
    
    # Create visualizations
    logger.info("Creating visualizations...")
    create_enrichment_visualizations(
        enrichment_scores, query_adata, cell_type_enrichment,
        output_dir, dataset, method, size
    )
    
    # Compile results
    # Note: enrichment_scores is saved separately to avoid memory issues
    results = {
        'dataset': dataset,
        'ref': ref,
        'method': method,
        'sample_size': size,
        'rep': rep,
        'n_reference': len(ref_idx),
        'n_query': len(query_idx),
        'enrichment_stats': {
            'mean': float(np.mean(enrichment_scores)),
            'median': float(np.median(enrichment_scores)),
            'std': float(np.std(enrichment_scores)),
            'min': float(np.min(enrichment_scores)),
            'max': float(np.max(enrichment_scores)),
            'pct_high': float(np.sum(enrichment_scores > np.percentile(enrichment_scores, 75)) / len(enrichment_scores) * 100),
            'pct_low': float(np.sum(enrichment_scores < np.percentile(enrichment_scores, 25)) / len(enrichment_scores) * 100)
        },
        'cell_type_enrichment': cell_type_enrichment,
        'n_poorly_represented': int(np.sum(poorly_represented)),
        'poorly_represented_threshold': float(threshold),
        'enrichment_scores': enrichment_scores  # Save enrichment scores (1D array, manageable size)
    }
    
    # Save results
    results_file = os.path.join(output_dir, 
                                f"{dataset}_{method}_size{size}_rep{rep}_results.pkl")
    with open(results_file, 'wb') as f:
        pickle.dump(results, f)
    logger.info(f"Saved results to {results_file}")
    
    # Save summary CSV
    summary_df = pd.DataFrame({
        'cell_type': list(cell_type_enrichment.keys()),
        'mean_enrichment': [cell_type_enrichment[ct]['mean_enrichment'] 
                           for ct in cell_type_enrichment.keys()],
        'median_enrichment': [cell_type_enrichment[ct]['median_enrichment'] 
                            for ct in cell_type_enrichment.keys()],
        'n_cells': [cell_type_enrichment[ct]['n_cells'] 
                   for ct in cell_type_enrichment.keys()],
        'pct_high_enrichment': [cell_type_enrichment[ct]['pct_high_enrichment'] 
                               for ct in cell_type_enrichment.keys()]
    })
    summary_df = summary_df.sort_values('mean_enrichment', ascending=False)
    
    csv_file = os.path.join(output_dir, 
                           f"{dataset}_{method}_size{size}_rep{rep}_summary.csv")
    summary_df.to_csv(csv_file, index=False)
    logger.info(f"Saved summary to {csv_file}")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Run SEMITONES downstream analysis using SPS samples as references'
    )
    parser.add_argument('--dataset', type=str, required=True,
                       choices=['mcc', 'mcc_01', 'mcc_05', 'lcmv'],
                       help='Dataset name')
    parser.add_argument('--ref', type=int, required=True,
                       help='Reference size')
    parser.add_argument('--method', type=str, required=True,
                       choices=['random', 'sps', 'hopper', 'atomic', 'scsampler'],
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
    
    data_path = get_data_path()
    output_dir = get_output_dir()
    
    results = run_downstream_analysis(
        args.dataset, args.ref, args.method, args.size, args.rep,
        data_path, output_dir
    )
    
    if results is None:
        logger.error("Analysis failed")
        sys.exit(1)
    
    logger.info("Downstream analysis completed successfully!")

if __name__ == '__main__':
    main()

