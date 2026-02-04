#!/usr/bin/env python
# coding: utf-8

"""
Downstream Analysis Comparison for Reviewer Response

This script addresses reviewer comments by:
1. Calculating cell type/state diversity preservation at different sampling depths
2. Comparing DE genes between original and subsampled populations
3. Comparing clustering/cell-type annotation results between original and subsampled populations
"""

import os
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from datetime import datetime
from itertools import product
from scipy import stats
# Removed sklearn clustering metrics - using simpler cell type annotation comparison instead
import scanpy as sc
import scipy.sparse

# Add analysis directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

from refined_rare_cell_type_definition import (
    calculate_cell_type_distances,
    identify_rare_cell_types_distance_and_frequency,
    load_reference_data as load_ref_data
)

# Set up logging
def setup_logger():
    log_dir = os.path.join(os.path.dirname(script_dir), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'downstream_analysis_comparison_{timestamp}.log')
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
        'frequency_threshold_pct': 0.1,
        'x_values': [0.5, 1.0, 2.0, 2.5, 3.2]
    },
    'mcc_05': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'frequency_threshold_pct': 0.5,
        'x_values': [0.5, 1.0, 2.0, 2.5, 3.2]
    },
    'lcmv': {
        'references': [1, 5, 10, 20, 34],
        'sizes': [50000, 100000, 200000],
        'frequency_threshold_pct': 1.0,
        'x_values': [0.1, 0.5, 1.0, 2.0, 3.4]
    }
}

REPS = [0]  # Using rep 0 for now, can be extended

def get_data_path():
    """Get data path from environment variable or derive from script location."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
    return os.path.join(project_root, 'data')

# Removed get_sparseflow_data_path - all data is in the main data directory

def get_output_dir():
    """Get output directory for results."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, 'results', 'downstream_analysis')
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

def calculate_rare_cell_type_metrics(obs, rare_types, label_key='celltype'):
    """
    Calculate metrics that emphasize rare cell type preservation.
    
    Parameters:
    -----------
    obs : pd.DataFrame
        Observation dataframe with cell type labels
    rare_types : list
        List of rare cell type names
    label_key : str
        Key in obs containing cell type labels
        
    Returns:
    --------
    dict: Dictionary of metrics
    """
    cell_type_counts = obs[label_key].value_counts()
    total_cells = len(obs)
    
    # Rare cell type metrics
    rare_cell_types_present = [ct for ct in rare_types if ct in cell_type_counts.index]
    n_rare_types_present = len(rare_cell_types_present)
    n_rare_types_total = len(rare_types)
    rare_type_coverage = n_rare_types_present / max(n_rare_types_total, 1)
    
    # Count rare cells
    rare_cell_count = cell_type_counts[cell_type_counts.index.isin(rare_types)].sum()
    rare_cell_fraction = rare_cell_count / max(total_cells, 1)
    
    # Cell type richness (total number of unique types)
    n_cell_types = len(cell_type_counts)
    
    # Weighted cell type preservation (gives more weight to rare types)
    # Calculate a score where rare types contribute more
    if len(rare_types) > 0:
        rare_weight = 10.0  # Weight for rare types
        common_weight = 1.0  # Weight for common types
        
        weighted_score = 0
        for cell_type, count in cell_type_counts.items():
            if cell_type in rare_types:
                weighted_score += rare_weight * (count > 0)  # Binary: present or not
            else:
                weighted_score += common_weight * (count > 0)
        
        # Normalize by maximum possible score
        max_score = len(rare_types) * rare_weight + (n_cell_types - len(rare_types)) * common_weight
        weighted_preservation = weighted_score / max(max_score, 1)
    else:
        weighted_preservation = 0
    
    metrics = {
        'n_cell_types': n_cell_types,
        'n_rare_types_present': n_rare_types_present,
        'n_rare_types_total': n_rare_types_total,
        'rare_type_coverage': rare_type_coverage,
        'rare_cell_count': rare_cell_count,
        'rare_cell_fraction': rare_cell_fraction,
        'weighted_preservation': weighted_preservation,
        'total_cells': total_cells,
        'min_cell_type_count': cell_type_counts.min(),
        'max_cell_type_count': cell_type_counts.max(),
        'median_cell_type_count': cell_type_counts.median(),
    }
    
    return metrics, cell_type_counts

def calculate_diversity_metrics(obs, rare_types, label_key='celltype'):
    """
    Calculate diversity metrics focused on rare cell type preservation.
    Replaces Shannon/Simpson with metrics that emphasize rare cell types.
    """
    return calculate_rare_cell_type_metrics(obs, rare_types, label_key)

def analyze_diversity_preservation(dataset, data_path):
    """Analyze cell type diversity preservation across different sampling methods and depths."""
    logger.info(f"Analyzing diversity preservation for {dataset}")
    
    config = DATASET_CONFIGS[dataset]
    references = config['references']
    sizes = config['sizes']
    frequency_threshold_pct = config['frequency_threshold_pct']
    
    # All data is in the main data directory
    base_path = data_path
    
    results = []
    
    for ref in references:
        logger.info(f"Processing reference {ref} for {dataset}")
        
        # Load original data
        adata = load_ref_data(dataset, ref, base_path=base_path)
        obs_original = adata.obs.copy()
        
        # Identify rare cell types
        from refined_rare_cell_type_definition import (
            calculate_cell_type_distances,
            identify_rare_cell_types_distance_and_frequency
        )
        cell_type_df = calculate_cell_type_distances(adata, label_key=label_key)
        total_cells = obs_original.shape[0]
        rare_types, _, _ = identify_rare_cell_types_distance_and_frequency(
            cell_type_df, 
            total_cells, 
            distance_percentile=75, 
            frequency_threshold_pct=frequency_threshold_pct
        )
        
        # Calculate original diversity metrics (with rare types)
        original_metrics, original_counts = calculate_diversity_metrics(
            obs_original, rare_types, label_key
        )
        
        # Load obs.csv if available for better cell type info
        directory = f"{dataset}/benchmark"
        path = os.path.join(base_path, directory)
        obs_file = os.path.join(path, f"{ref}/obs.csv")
        if os.path.isfile(obs_file):
            obs_df = pd.read_csv(obs_file, index_col=0)
            obs_df[label_key] = obs_df[label_key].astype('category')
        else:
            obs_df = obs_original
        
        for size, method, rep in product(sizes, METHODS, REPS):
            # Load sampling indices
            indices = load_sampling_indices(dataset, ref, method, size, rep, base_path)
            
            if indices is None:
                continue
            
            # Get sampled observations
            if method == 'atomic' and dataset in ['mcc_01', 'mcc_05']:
                # For atomic method with string indices
                sampled_obs = obs_df.loc[obs_df.index.intersection(indices)]
            else:
                # For other methods with integer indices
                if isinstance(indices[0], (str, np.str_)):
                    sampled_obs = obs_df.loc[obs_df.index.intersection(indices)]
                else:
                    sampled_obs = obs_df.iloc[indices]
            
            if len(sampled_obs) == 0:
                logger.warning(f"No cells sampled for {dataset}, ref {ref}, method {method}, size {size}, rep {rep}")
                continue
            
            # Calculate sampled diversity metrics (with rare types)
            sampled_metrics, sampled_counts = calculate_diversity_metrics(
                sampled_obs, rare_types, label_key
            )
            
            # Calculate preservation ratios
            n_cell_types_preserved = sampled_metrics['n_cell_types'] / max(original_metrics['n_cell_types'], 1)
            rare_type_coverage_preserved = sampled_metrics['rare_type_coverage'] / max(original_metrics['rare_type_coverage'], 1)
            
            # Rare cell enrichment: ratio of rare cell fraction in sample vs original
            original_rare_fraction = original_metrics['rare_cell_fraction']
            sampled_rare_fraction = sampled_metrics['rare_cell_fraction']
            rare_cell_enrichment = sampled_rare_fraction / max(original_rare_fraction, 1e-6) if original_rare_fraction > 0 else 0
            
            # Calculate cell type overlap
            original_types = set(original_counts.index)
            sampled_types = set(sampled_counts.index)
            cell_type_overlap = len(original_types.intersection(sampled_types)) / max(len(original_types), 1)
            
            result = {
                'dataset': dataset,
                'ref': ref,
                'method': method,
                'sample_size': size,
                'rep': rep,
                'original_n_cell_types': original_metrics['n_cell_types'],
                'sampled_n_cell_types': sampled_metrics['n_cell_types'],
                'n_cell_types_preserved': n_cell_types_preserved,
                'original_n_rare_types': original_metrics['n_rare_types_total'],
                'sampled_n_rare_types_present': sampled_metrics['n_rare_types_present'],
                'rare_type_coverage': sampled_metrics['rare_type_coverage'],
                'rare_type_coverage_preserved': rare_type_coverage_preserved,
                'original_rare_cell_fraction': original_rare_fraction,
                'sampled_rare_cell_fraction': sampled_rare_fraction,
                'rare_cell_enrichment': rare_cell_enrichment,
                'weighted_preservation': sampled_metrics['weighted_preservation'],
                'cell_type_overlap': cell_type_overlap,
                'original_total_cells': original_metrics['total_cells'],
                'sampled_total_cells': sampled_metrics['total_cells'],
            }
            
            results.append(result)
    
    results_df = pd.DataFrame(results)
    return results_df

def perform_de_analysis(adata_original, adata_sampled, label_key='celltype', n_top_genes=50):
    """Perform differential expression analysis and compare results."""
    # Get common cell types
    original_types = set(adata_original.obs[label_key].unique())
    sampled_types = set(adata_sampled.obs[label_key].unique())
    common_types = original_types.intersection(sampled_types)
    
    if len(common_types) < 2:
        logger.warning("Not enough common cell types for DE analysis")
        return None
    
    de_results = {}
    
    # For each cell type, find DE genes in original vs sampled
    for cell_type in list(common_types)[:5]:  # Limit to first 5 for efficiency
        try:
            # Original data
            adata_orig_subset = adata_original[adata_original.obs[label_key] == cell_type].copy()
            adata_orig_other = adata_original[adata_original.obs[label_key] != cell_type].copy()
            
            if adata_orig_subset.shape[0] < 10 or adata_orig_other.shape[0] < 10:
                continue
            
            # Sampled data
            adata_samp_subset = adata_sampled[adata_sampled.obs[label_key] == cell_type].copy()
            adata_samp_other = adata_sampled[adata_sampled.obs[label_key] != cell_type].copy()
            
            if adata_samp_subset.shape[0] < 10 or adata_samp_other.shape[0] < 10:
                continue
            
            # Create temporary adata for DE - use concatenate properly
            adata_orig_de = adata_orig_subset.concatenate(adata_orig_other, batch_key='batch_temp')
            adata_orig_de.obs['group'] = 'target'
            adata_orig_de.obs.loc[adata_orig_de.obs[label_key] != cell_type, 'group'] = 'other'
            
            adata_samp_de = adata_samp_subset.concatenate(adata_samp_other, batch_key='batch_temp')
            adata_samp_de.obs['group'] = 'target'
            adata_samp_de.obs.loc[adata_samp_de.obs[label_key] != cell_type, 'group'] = 'other'
            
            # Perform DE - use use_raw=False if raw doesn't exist
            try:
                sc.tl.rank_genes_groups(adata_orig_de, groupby='group', groups=['target'], 
                                       reference='other', method='wilcoxon', n_genes=n_top_genes,
                                       use_raw=False)
                sc.tl.rank_genes_groups(adata_samp_de, groupby='group', groups=['target'], 
                                       reference='other', method='wilcoxon', n_genes=n_top_genes,
                                       use_raw=False)
            except:
                # Fallback to t-test if wilcoxon fails
                sc.tl.rank_genes_groups(adata_orig_de, groupby='group', groups=['target'], 
                                       reference='other', method='t-test', n_genes=n_top_genes,
                                       use_raw=False)
                sc.tl.rank_genes_groups(adata_samp_de, groupby='group', groups=['target'], 
                                       reference='other', method='t-test', n_genes=n_top_genes,
                                       use_raw=False)
            
            # Get top genes
            orig_genes_df = pd.DataFrame(adata_orig_de.uns['rank_genes_groups']['names'])
            samp_genes_df = pd.DataFrame(adata_samp_de.uns['rank_genes_groups']['names'])
            
            # Get top n_top_genes
            orig_gene_list = orig_genes_df.iloc[:n_top_genes, 0].tolist()
            samp_gene_list = samp_genes_df.iloc[:n_top_genes, 0].tolist()
            
            # Calculate overlap
            orig_gene_set = set(orig_gene_list)
            samp_gene_set = set(samp_gene_list)
            overlap = len(orig_gene_set.intersection(samp_gene_set)) / max(len(orig_gene_set), 1)
            
            de_results[cell_type] = {
                'original_genes': list(orig_gene_set),
                'sampled_genes': list(samp_gene_set),
                'overlap': overlap,
                'n_original_genes': len(orig_gene_set),
                'n_sampled_genes': len(samp_gene_set),
                'n_overlap': len(orig_gene_set.intersection(samp_gene_set))
            }
        except Exception as e:
            logger.warning(f"Error in DE analysis for cell type {cell_type}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    return de_results

def perform_cell_type_annotation_comparison(adata_original, adata_sampled, label_key='celltype'):
    """
    Compare cell type annotations between original and sampled data.
    This is a simpler approach that doesn't require clustering the full dataset.
    """
    try:
        # Compare cell type annotation preservation
        orig_types = set(adata_original.obs[label_key].unique())
        samp_types = set(adata_sampled.obs[label_key].unique())
        type_overlap = len(orig_types.intersection(samp_types)) / max(len(orig_types), 1)
        
        # Calculate cell type frequency preservation
        orig_counts = adata_original.obs[label_key].value_counts()
        samp_counts = adata_sampled.obs[label_key].value_counts()
        
        # For common cell types, calculate correlation of frequencies
        common_types = orig_types.intersection(samp_types)
        if len(common_types) > 1:
            orig_freqs = [orig_counts.get(ct, 0) / len(adata_original) for ct in common_types]
            samp_freqs = [samp_counts.get(ct, 0) / len(adata_sampled) for ct in common_types]
            frequency_correlation = np.corrcoef(orig_freqs, samp_freqs)[0, 1] if len(orig_freqs) > 1 else None
        else:
            frequency_correlation = None
        
        # Calculate Jaccard similarity for cell type presence
        jaccard_similarity = len(common_types) / max(len(orig_types.union(samp_types)), 1)
        
        return {
            'cell_type_overlap': type_overlap,
            'frequency_correlation': frequency_correlation,
            'jaccard_similarity': jaccard_similarity,
            'n_original_types': len(orig_types),
            'n_sampled_types': len(samp_types),
            'n_common_types': len(common_types)
        }
    except Exception as e:
        logger.warning(f"Error in cell type annotation comparison: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

def analyze_downstream_comparisons(dataset, data_path, max_configs=5):
    """Analyze downstream comparisons (DE and clustering) for a subset of configurations."""
    logger.info(f"Analyzing downstream comparisons for {dataset}")
    
    config = DATASET_CONFIGS[dataset]
    references = config['references']
    sizes = config['sizes']
    
    # All data is in the main data directory
    base_path = data_path
    
    results = []
    config_count = 0
    
    # Sample a few configurations to analyze (to save computation time)
    for ref in references[-2:]:  # Use last 2 references (largest datasets)
        if config_count >= max_configs:
            break
        
        logger.info(f"Processing reference {ref} for {dataset}")
        
        # Load original data
        adata_original = load_ref_data(dataset, ref, base_path=base_path)
        
        for size in sizes[-2:]:  # Use largest 2 sample sizes
            if config_count >= max_configs:
                break
            
            for method in ['sps', 'random']:  # Compare best method vs baseline
                if config_count >= max_configs:
                    break
                
                for rep in REPS:
                    # Load sampling indices
                    indices = load_sampling_indices(dataset, ref, method, size, rep, base_path)
                    
                    if indices is None:
                        continue
                    
                    # Create sampled adata
                    if method == 'atomic' and dataset in ['mcc_01', 'mcc_05']:
                        adata_sampled = adata_original[adata_original.obs.index.intersection(indices)].copy()
                    else:
                        if isinstance(indices[0], (str, np.str_)):
                            adata_sampled = adata_original[adata_original.obs.index.intersection(indices)].copy()
                        else:
                            adata_sampled = adata_original[indices].copy()
                    
                    if adata_sampled.shape[0] < 100:
                        continue
                    
                    logger.info(f"Analyzing {dataset}, ref {ref}, method {method}, size {size}, rep {rep}")
                    
                    # Perform cell type annotation comparison (no clustering needed)
                    annotation_results = perform_cell_type_annotation_comparison(
                        adata_original, adata_sampled, label_key=label_key
                    )
                    
                    # Perform DE analysis (simplified - just calculate overlap metric)
                    de_results = perform_de_analysis(adata_original, adata_sampled, label_key=label_key)
                    
                    result = {
                        'dataset': dataset,
                        'ref': ref,
                        'method': method,
                        'sample_size': size,
                        'rep': rep,
                    }
                    
                    if annotation_results:
                        result.update(annotation_results)
                    
                    if de_results:
                        # Average overlap across cell types
                        avg_de_overlap = np.mean([v['overlap'] for v in de_results.values()])
                        result['avg_de_gene_overlap'] = avg_de_overlap
                        result['n_cell_types_analyzed'] = len(de_results)
                    else:
                        result['avg_de_gene_overlap'] = None
                        result['n_cell_types_analyzed'] = 0
                    
                    results.append(result)
                    config_count += 1
                    
                    if config_count >= max_configs:
                        break
    
    results_df = pd.DataFrame(results)
    return results_df

def main():
    """Main function to run all analyses."""
    logger.info("Starting downstream analysis comparison")
    
    data_path = get_data_path()
    output_dir = get_output_dir()
    
    # Analyze diversity preservation for all datasets
    all_diversity_results = []
    for dataset in ['mcc_01', 'mcc_05']:  # Focus on main datasets
        logger.info(f"Processing diversity preservation for {dataset}")
        diversity_df = analyze_diversity_preservation(dataset, data_path)
        all_diversity_results.append(diversity_df)
        
        # Save individual results
        output_file = os.path.join(output_dir, f'{dataset}_diversity_preservation.csv')
        diversity_df.to_csv(output_file, index=False)
        logger.info(f"Saved diversity results to {output_file}")
    
    # Combine and save all diversity results
    if all_diversity_results:
        combined_diversity = pd.concat(all_diversity_results, ignore_index=True)
        combined_file = os.path.join(output_dir, 'combined_diversity_preservation.csv')
        combined_diversity.to_csv(combined_file, index=False)
        logger.info(f"Saved combined diversity results to {combined_file}")
    
    # Analyze downstream comparisons (DE and clustering) for a subset
    all_downstream_results = []
    for dataset in ['mcc_01']:  # Focus on one dataset for detailed analysis
        logger.info(f"Processing downstream comparisons for {dataset}")
        downstream_df = analyze_downstream_comparisons(dataset, data_path, max_configs=10)
        all_downstream_results.append(downstream_df)
        
        # Save individual results
        output_file = os.path.join(output_dir, f'{dataset}_downstream_comparison.csv')
        downstream_df.to_csv(output_file, index=False)
        logger.info(f"Saved downstream results to {output_file}")
    
    # Combine and save all downstream results
    if all_downstream_results:
        combined_downstream = pd.concat(all_downstream_results, ignore_index=True)
        combined_file = os.path.join(output_dir, 'combined_downstream_comparison.csv')
        combined_downstream.to_csv(combined_file, index=False)
        logger.info(f"Saved combined downstream results to {combined_file}")
    
    logger.info("All analyses completed successfully")

if __name__ == "__main__":
    main()

