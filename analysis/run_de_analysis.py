#!/usr/bin/env python
# coding: utf-8

"""
DE Gene Analysis - Focused script to run DE analysis comparison
between full datasets and subsampled data.
"""

import os
import sys
import numpy as np
import pandas as pd
import pickle
import logging
import argparse
from datetime import datetime
import scanpy as sc
import scipy.sparse

# Add analysis directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

from refined_rare_cell_type_definition import (
    load_reference_data as load_ref_data
)

# Set up logging
def setup_logger():
    log_dir = os.path.join(os.path.dirname(script_dir), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'de_analysis_{timestamp}.log')
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

# Dataset configurations - focusing on full datasets
DATASET_CONFIGS = {
    'mcc': {
        'full_ref': 30,  # Full MCC dataset
        'sizes': [100000, 200000],  # Sample sizes to test
        'methods': ['sps', 'random'],  # Compare best method vs baseline
    },
    'lcmv': {
        'full_ref': 34,  # Full LCMV dataset
        'sizes': [100000, 200000],  # Sample sizes to test
        'methods': ['sps', 'random'],
    }
}

REPS = [0]  # Using rep 0

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
    output_dir = os.path.join(script_dir, 'results', 'de_analysis')
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

def get_n_top_genes_for_dataset(dataset, adata):
    """
    Get appropriate number of top genes based on dataset.
    LCMV has fewer genes, so use fewer top genes.
    """
    total_genes = adata.shape[1]
    
    if dataset == 'lcmv':
        # LCMV has < 45 genes, so use 5-10 top genes
        n_top = min(10, max(5, int(total_genes * 0.2)))  # 20% of genes, but between 5-10
        logger.info(f"LCMV dataset: {total_genes} total genes, using {n_top} top genes")
    else:
        # MCC and others: use 50 genes (or 20% if fewer genes available)
        n_top = min(50, max(10, int(total_genes * 0.2)))
        logger.info(f"{dataset.upper()} dataset: {total_genes} total genes, using {n_top} top genes")
    
    return n_top

def perform_de_analysis_single(adata, cell_type, label_key='celltype', n_top_genes=50):
    """
    Perform DE analysis for a single cell type.
    
    Returns:
    --------
    list: Top n_top_genes gene names
    """
    try:
        # Split into target cell type vs all others
        adata_target = adata[adata.obs[label_key] == cell_type].copy()
        adata_other = adata[adata.obs[label_key] != cell_type].copy()
        
        if adata_target.shape[0] < 10 or adata_other.shape[0] < 10:
            return None
        
        # Combine and label
        adata_combined = adata_target.concatenate(adata_other, batch_key='batch_temp')
        adata_combined.obs['group'] = 'target'
        adata_combined.obs.loc[adata_combined.obs[label_key] != cell_type, 'group'] = 'other'
        
        # Perform DE - use use_raw=False if raw doesn't exist
        try:
            sc.tl.rank_genes_groups(
                adata_combined, 
                groupby='group', 
                groups=['target'],
                reference='other', 
                method='wilcoxon', 
                n_genes=n_top_genes,
                use_raw=False
            )
        except:
            # Fallback to t-test if wilcoxon fails
            sc.tl.rank_genes_groups(
                adata_combined, 
                groupby='group', 
                groups=['target'],
                reference='other', 
                method='t-test', 
                n_genes=n_top_genes,
                use_raw=False
            )
        
        # Get top genes
        genes_df = pd.DataFrame(adata_combined.uns['rank_genes_groups']['names'])
        top_genes = genes_df.iloc[:n_top_genes, 0].tolist()
        
        return top_genes
        
    except Exception as e:
        logger.warning(f"Error in DE analysis for cell type {cell_type}: {str(e)}")
        return None

def compare_de_results(adata_full, adata_sampled, label_key='celltype', n_top_genes=50, max_cell_types=None):
    """
    Compare DE results between full and sampled datasets.
    
    Parameters:
    -----------
    adata_full : AnnData
        Full dataset
    adata_sampled : AnnData
        Sampled dataset
    label_key : str
        Key for cell type labels
    n_top_genes : int
        Number of top genes to compare
    max_cell_types : int, optional
        Maximum number of cell types to analyze. If None, analyzes all cell types.
        
    Returns:
    --------
    dict: Comparison results
    """
    # Get common cell types
    full_types = set(adata_full.obs[label_key].unique())
    sampled_types = set(adata_sampled.obs[label_key].unique())
    common_types = sorted(list(full_types.intersection(sampled_types)))
    
    if len(common_types) < 2:
        logger.warning("Not enough common cell types for DE analysis")
        return None
    
    # Limit to max_cell_types if specified, otherwise analyze all
    if max_cell_types is not None:
        common_types = common_types[:max_cell_types]
        logger.info(f"Analyzing {len(common_types)} cell types (limited from {len(full_types.intersection(sampled_types))} total)")
    else:
        logger.info(f"Analyzing all {len(common_types)} cell types")
    
    logger.info(f"Analyzing DE for {len(common_types)} cell types: {common_types}")
    
    results = []
    
    for cell_type in common_types:
        logger.info(f"  Processing cell type: {cell_type}")
        
        # Get DE genes from full dataset
        full_genes = perform_de_analysis_single(adata_full, cell_type, label_key, n_top_genes)
        
        if full_genes is None:
            continue
        
        # Get DE genes from sampled dataset
        sampled_genes = perform_de_analysis_single(adata_sampled, cell_type, label_key, n_top_genes)
        
        if sampled_genes is None:
            continue
        
        # Calculate overlap
        full_gene_set = set(full_genes)
        sampled_gene_set = set(sampled_genes)
        overlap = len(full_gene_set.intersection(sampled_gene_set)) / max(len(full_gene_set), 1)
        
        results.append({
            'cell_type': cell_type,
            'n_full_genes': len(full_gene_set),
            'n_sampled_genes': len(sampled_gene_set),
            'n_overlap': len(full_gene_set.intersection(sampled_gene_set)),
            'overlap_fraction': overlap,
            'full_genes': ','.join(full_genes),
            'sampled_genes': ','.join(sampled_genes),
            'overlapping_genes': ','.join(full_gene_set.intersection(sampled_gene_set))
        })
        
        logger.info(f"    Overlap: {overlap:.2%} ({len(full_gene_set.intersection(sampled_gene_set))}/{len(full_gene_set)})")
    
    if len(results) == 0:
        return None
    
    # Calculate summary statistics
    overlaps = [r['overlap_fraction'] for r in results]
    
    summary = {
        'mean_overlap': np.mean(overlaps),
        'median_overlap': np.median(overlaps),
        'std_overlap': np.std(overlaps),
        'min_overlap': np.min(overlaps),
        'max_overlap': np.max(overlaps),
        'n_cell_types_analyzed': len(results)
    }
    
    return {
        'summary': summary,
        'by_cell_type': results
    }

def analyze_de_for_dataset(dataset, data_path, output_dir):
    """Run DE analysis for a single dataset."""
    logger.info(f"\n{'='*80}")
    logger.info(f"Starting DE analysis for {dataset}")
    logger.info(f"{'='*80}")
    
    config = DATASET_CONFIGS[dataset]
    full_ref = config['full_ref']
    sizes = config['sizes']
    methods = config['methods']
    
    base_path = data_path
    
    # Load full dataset
    logger.info(f"Loading full dataset: {dataset}, ref {full_ref}")
    adata_full = load_ref_data(dataset, full_ref, base_path=base_path)
    logger.info(f"Full dataset shape: {adata_full.shape}")
    
    # Determine appropriate number of top genes for this dataset
    n_top_genes = get_n_top_genes_for_dataset(dataset, adata_full)
    
    all_results = []
    
    for size in sizes:
        for method in methods:
            for rep in REPS:
                logger.info(f"\nProcessing: {dataset}, method={method}, size={size}, rep={rep}")
                
                # Load sampling indices
                indices = load_sampling_indices(dataset, full_ref, method, size, rep, base_path)
                
                if indices is None:
                    logger.warning(f"  No indices found, skipping")
                    continue
                
                # Create sampled adata
                if isinstance(indices[0], (str, np.str_)):
                    adata_sampled = adata_full[adata_full.obs.index.intersection(indices)].copy()
                else:
                    adata_sampled = adata_full[indices].copy()
                
                logger.info(f"  Sampled dataset shape: {adata_sampled.shape}")
                
                if adata_sampled.shape[0] < 100:
                    logger.warning(f"  Too few cells, skipping")
                    continue
                
                # Compare DE results - analyze all cell types
                de_comparison = compare_de_results(
                    adata_full, 
                    adata_sampled, 
                    label_key=label_key,
                    n_top_genes=n_top_genes,
                    max_cell_types=None  # Analyze all cell types
                )
                
                if de_comparison is None:
                    logger.warning(f"  No DE comparison results")
                    continue
                
                # Add metadata
                result_entry = {
                    'dataset': dataset,
                    'method': method,
                    'sample_size': size,
                    'rep': rep,
                    'full_n_cells': adata_full.shape[0],
                    'sampled_n_cells': adata_sampled.shape[0],
                    'n_top_genes': n_top_genes,  # Store the number of top genes used
                }
                result_entry.update(de_comparison['summary'])
                
                all_results.append(result_entry)
                
                # Save detailed results for this configuration
                detailed_file = os.path.join(
                    output_dir, 
                    f'{dataset}_{method}_size{size}_rep{rep}_de_detailed.csv'
                )
                detailed_df = pd.DataFrame(de_comparison['by_cell_type'])
                detailed_df.to_csv(detailed_file, index=False)
                logger.info(f"  Saved detailed results to {detailed_file}")
                logger.info(f"  Mean overlap: {de_comparison['summary']['mean_overlap']:.2%}")
                logger.info(f"  Used {n_top_genes} top genes for comparison")
    
    return all_results

def run_single_de_analysis(dataset, method, size, rep, data_path, output_dir):
    """Run DE analysis for a single configuration."""
    logger.info(f"\n{'='*80}")
    logger.info(f"Starting DE analysis: {dataset}, method={method}, size={size}, rep={rep}")
    logger.info(f"{'='*80}")
    
    config = DATASET_CONFIGS[dataset]
    full_ref = config['full_ref']
    
    base_path = data_path
    
    # Load full dataset
    logger.info(f"Loading full dataset: {dataset}, ref {full_ref}")
    adata_full = load_ref_data(dataset, full_ref, base_path=base_path)
    logger.info(f"Full dataset shape: {adata_full.shape}")
    
    # Determine appropriate number of top genes for this dataset
    n_top_genes = get_n_top_genes_for_dataset(dataset, adata_full)
    
    # Load sampling indices
    indices = load_sampling_indices(dataset, full_ref, method, size, rep, base_path)
    
    if indices is None:
        logger.warning(f"No indices found for {dataset}, ref {full_ref}, method {method}, size {size}, rep {rep}")
        return None
    
    # Create sampled adata
    if isinstance(indices[0], (str, np.str_)):
        adata_sampled = adata_full[adata_full.obs.index.intersection(indices)].copy()
    else:
        adata_sampled = adata_full[indices].copy()
    
    logger.info(f"Sampled dataset shape: {adata_sampled.shape}")
    
    if adata_sampled.shape[0] < 100:
        logger.warning(f"Too few cells ({adata_sampled.shape[0]}), skipping")
        return None
    
    # Compare DE results
    de_comparison = compare_de_results(
        adata_full, 
        adata_sampled, 
        label_key=label_key,
        n_top_genes=n_top_genes,
        max_cell_types=None  # Analyze all cell types
    )
    
    if de_comparison is None:
        logger.warning(f"No DE comparison results")
        return None
    
    # Add metadata
    result_entry = {
        'dataset': dataset,
        'method': method,
        'sample_size': size,
        'rep': rep,
        'full_n_cells': adata_full.shape[0],
        'sampled_n_cells': adata_sampled.shape[0],
        'n_top_genes': n_top_genes,  # Store the number of top genes used
    }
    result_entry.update(de_comparison['summary'])
    
    # Save detailed results for this configuration
    detailed_file = os.path.join(
        output_dir, 
        f'{dataset}_{method}_size{size}_rep{rep}_de_detailed.csv'
    )
    detailed_df = pd.DataFrame(de_comparison['by_cell_type'])
    detailed_df.to_csv(detailed_file, index=False)
    logger.info(f"Saved detailed results to {detailed_file}")
    logger.info(f"Mean overlap: {de_comparison['summary']['mean_overlap']:.2%}")
    logger.info(f"Used {n_top_genes} top genes for comparison")
    
    return result_entry

def main():
    """Main function to run DE analysis."""
    parser = argparse.ArgumentParser(description='Run DE gene analysis comparison')
    parser.add_argument('--dataset', type=str, choices=['mcc', 'lcmv'], 
                       help='Dataset to analyze')
    parser.add_argument('--method', type=str, 
                       choices=['random', 'sps', 'hopper', 'atomic', 'scsampler'],
                       help='Sampling method')
    parser.add_argument('--size', type=int, help='Sample size')
    parser.add_argument('--rep', type=int, help='Replicate number')
    
    args = parser.parse_args()
    
    data_path = get_data_path()
    output_dir = get_output_dir()
    
    # If arguments provided, run single analysis
    if args.dataset and args.method and args.size is not None and args.rep is not None:
        logger.info("Running single DE analysis with provided arguments")
        result = run_single_de_analysis(
            args.dataset, args.method, args.size, args.rep, data_path, output_dir
        )
        
        if result:
            # Append to summary file
            summary_file = os.path.join(output_dir, 'de_analysis_summary.csv')
            if os.path.exists(summary_file):
                summary_df = pd.read_csv(summary_file)
                summary_df = pd.concat([summary_df, pd.DataFrame([result])], ignore_index=True)
            else:
                summary_df = pd.DataFrame([result])
            summary_df.to_csv(summary_file, index=False)
            logger.info(f"Updated summary file: {summary_file}")
    else:
        # Run for all configurations (original behavior)
        logger.info("Running DE analysis for all configurations")
        all_results = []
        
        # Run for each dataset
        for dataset in ['mcc', 'lcmv']:
            try:
                results = analyze_de_for_dataset(dataset, data_path, output_dir)
                all_results.extend(results)
            except Exception as e:
                logger.error(f"Error processing {dataset}: {str(e)}")
                import traceback
                traceback.print_exc()
                continue
        
        # Save summary results
        if all_results:
            summary_df = pd.DataFrame(all_results)
            summary_file = os.path.join(output_dir, 'de_analysis_summary.csv')
            summary_df.to_csv(summary_file, index=False)
            logger.info(f"\nSaved summary results to {summary_file}")
            logger.info(f"\nSummary:\n{summary_df[['dataset', 'method', 'sample_size', 'mean_overlap', 'n_cell_types_analyzed']]}")
        else:
            logger.warning("No results to save")
    
    logger.info("\nDE analysis completed!")

if __name__ == "__main__":
    main()

