#!/usr/bin/env python
# coding: utf-8

import os
import sys
import numpy as np
import pandas as pd
import pickle
import logging
from datetime import datetime
import scanpy as sc

# Add analysis directory to path to import rare cell type functions
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
analysis_dir = os.path.join(project_root, 'analysis')
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

from refined_rare_cell_type_definition import (
    calculate_cell_type_distances,
    identify_rare_cell_types_distance_and_frequency
)

# Set up logging
def setup_logger():
    log_dir = os.path.join(project_root, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'create_sampling_methods_table_{timestamp}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

# Constants
DATASET_CONFIGS = {
    'lcmv': {
        'references': [1, 5, 10, 20, 34],
        'sizes': [50000, 100000, 200000],
        'frequency_threshold_pct': 1.0
    },
    'mcc': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'frequency_threshold_pct': 1.0
    },
    'mcc_01': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'frequency_threshold_pct': 0.1
    },
    'mcc_05': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'frequency_threshold_pct': 0.5
    }
}

# Sampling methods to include
METHODS = ['random', 'sps', 'hopper', 'atomic', 'scsampler']
REP = 0  # Use repetition 0
label_key = 'celltype'

# File paths
def get_data_path():
    """Get data path from environment variable or derive from script location."""
    project_root_env = os.environ.get('PROJECT_ROOT')
    if project_root_env is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root_env = os.path.dirname(os.path.dirname(script_dir))
    return os.path.join(project_root_env, 'data')

def get_output_dir():
    """Get output directory for tables."""
    project_root_env = os.environ.get('PROJECT_ROOT')
    if project_root_env is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root_env = os.path.dirname(os.path.dirname(script_dir))
    return os.path.join(project_root_env, 'jobs', 'test_sampling_methods', 'tables')

file_path_env = get_data_path()
OUTPUT_DIR = get_output_dir()
os.makedirs(OUTPUT_DIR, exist_ok=True)


def get_benchmark_path(dataset):
    """Get benchmark path for a dataset."""
    return os.path.join(file_path_env, f"{dataset}/benchmark")


def load_reference_obs(dataset, ref):
    """Load reference observation data."""
    benchmark_path = get_benchmark_path(dataset)
    obs_path = os.path.join(benchmark_path, f"{ref}/obs.csv")
    if os.path.exists(obs_path):
        df = pd.read_csv(obs_path, index_col=0)
        df[label_key] = df[label_key].astype('category')
        return df
    else:
        adata_path = os.path.join(benchmark_path, f"{ref}/adata.h5ad")
        if os.path.exists(adata_path):
            adata = sc.read_h5ad(adata_path)
            return adata.obs
        else:
            logger.warning(f"Could not load obs for {dataset}, ref {ref}")
            return None


def identify_rare_cell_types(dataset, ref, obs):
    """Identify rare cell types for a dataset."""
    benchmark_path = get_benchmark_path(dataset)
    adata_path = os.path.join(benchmark_path, f"{ref}/adata.h5ad")
    
    if not os.path.exists(adata_path):
        logger.warning(f"Could not find adata for {dataset}, ref {ref}")
        return []
    
    adata = sc.read_h5ad(adata_path)
    config = DATASET_CONFIGS[dataset]
    
    cell_type_df = calculate_cell_type_distances(adata, label_key=label_key)
    
    total_cells = obs.shape[0]
    rare_types, _, _ = identify_rare_cell_types_distance_and_frequency(
        cell_type_df, 
        total_cells, 
        distance_percentile=75, 
        frequency_threshold_pct=config['frequency_threshold_pct']
    )
    
    return rare_types


def load_method_indices(dataset, ref, method, sample_size, rep):
    """Load cell indices for a given method.
    
    Args:
        dataset: Dataset name
        ref: Reference number
        method: Method name ('random', 'sps', 'hopper', 'atomic', 'scsampler')
        sample_size: Sample size
        rep: Repetition number
    
    Returns:
        For pickle methods: numpy array of integer indices (positional)
        For atomic: list of index labels (strings or integers depending on dataset)
        Returns None if file not found or error
    """
    benchmark_path = get_benchmark_path(dataset)
    
    try:
        if method == 'atomic':
            # Atomic uses CSV format
            results_path = os.path.join(benchmark_path, f"{ref}/atomic/{sample_size}/{rep}/results.csv")
            if not os.path.exists(results_path):
                logger.debug(f"Results not found: {results_path}")
                return None
            
            df = pd.read_csv(results_path)
            if dataset == 'mcc' or dataset.startswith('mcc_'):
                # For MCC datasets, indices are strings (index labels)
                indices = df['x'].values.astype(str).tolist()
            else:
                # For other datasets, indices might be integers (could be positional or labels)
                indices = df['x'].values.astype(int).tolist()
            return indices
        else:
            # Other methods use pickle format
            results_path = os.path.join(benchmark_path, f"{ref}/{method}/{sample_size}/{rep}/results.pkl")
            if not os.path.exists(results_path):
                logger.debug(f"Results not found: {results_path}")
                return None
            
            with open(results_path, 'rb') as f:
                results_data = pickle.load(f)
            
            # Results.pkl contains (indices, elapsed_time)
            # indices is the first element of the tuple, and is a numpy array of integer positions
            cell_indices = results_data[0]
            return cell_indices
    except Exception as e:
        logger.warning(f"Error loading indices for {dataset}, ref {ref}, method {method}, size {sample_size}, rep {rep}: {str(e)}")
        return None


def count_rare_cells(obs, indices, rare_types, method):
    """Count rare cells in the sampled indices.
    
    Args:
        obs: Observation DataFrame
        indices: For pickle methods: numpy array of integer positions
                 For atomic: list of index labels (strings or integers)
        rare_types: List of rare cell type names
        method: Method name to determine how to handle indices
    
    Returns:
        Number of rare cells found
    """
    if indices is None or len(indices) == 0:
        return 0
    
    try:
        if method == 'atomic':
            # Atomic provides index labels directly
            sampled_obs = obs.loc[indices]
        else:
            # Other methods provide integer positions (numpy array)
            # Use iloc with the indices directly (they are positional integer indices)
            sampled_obs = obs.iloc[indices]
        
        # Count rare cells in sample - count by cell type first, then filter to rare types
        value_counts = sampled_obs[label_key].value_counts()
        
        # Filter to only rare cell types and sum their counts
        rare_types_set = set(rare_types)
        sample_cell_types = set(value_counts.index)
        matching_rare_types = rare_types_set.intersection(sample_cell_types)
        
        if len(matching_rare_types) > 0:
            rare_cell_type_counts = value_counts[value_counts.index.isin(matching_rare_types)]
            rare_count = rare_cell_type_counts.sum()
        else:
            rare_count = 0
        
        return rare_count
    except Exception as e:
        logger.warning(f"Error counting rare cells: {str(e)}")
        return 0


def create_sampling_methods_table(dataset, sample_size, ref_data):
    """Create a table with methods as rows and references as columns.
    
    Args:
        dataset: Dataset name (e.g., 'lcmv', 'mcc', 'mcc_01', 'mcc_05')
        sample_size: Specific sample size to use (required)
        ref_data: Pre-loaded reference data dictionary with keys as ref numbers
                  and values as dicts with 'obs' and 'rare_types' keys
    
    Returns:
        DataFrame with method as rows and references as columns.
        Entries are the number of rare cells collected.
    """
    config = DATASET_CONFIGS[dataset]
    references = config['references']
    
    # Initialize table
    table_data = []
    
    for method in METHODS:
        row = {'method': method}
        
        for ref in references:
            if ref_data[ref]['obs'] is None:
                row[f'ref_{ref}M'] = 0
                continue
            
            obs = ref_data[ref]['obs']
            rare_types = ref_data[ref]['rare_types']
            
            if len(rare_types) == 0:
                row[f'ref_{ref}M'] = 0
                continue
            
            # Load indices for this method (rep 0)
            indices = load_method_indices(dataset, ref, method, sample_size, REP)
            
            if indices is None:
                row[f'ref_{ref}M'] = 0
                continue
            
            # Count rare cells
            rare_count = count_rare_cells(obs, indices, rare_types, method)
            row[f'ref_{ref}M'] = rare_count
            
            if rare_count > 0:
                logger.info(f"Dataset {dataset}, method {method}, ref {ref}, size {sample_size}, rep {REP}: Found {rare_count} rare cells")
            else:
                logger.debug(f"Dataset {dataset}, method {method}, ref {ref}, size {sample_size}, rep {REP}: No rare cells found")
        
        table_data.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(table_data)
    df.set_index('method', inplace=True)
    
    # Rename columns to be more readable (references as column names)
    df.columns = [f'{ref}M' for ref in references]
    
    return df


def main():
    logger.info("Creating sampling methods comparison tables for all datasets")
    
    # Only process mcc, mcc_01, mcc_05 (lcmv already completed)
    datasets = ['mcc', 'mcc_01', 'mcc_05']
    
    for dataset in datasets:
        logger.info(f"Processing {dataset} dataset...")
        config = DATASET_CONFIGS[dataset]
        references = config['references']
        
        # Load references once per dataset (not per size)
        logger.info(f"Loading references for {dataset} dataset...")
        ref_data = {}
        for ref in references:
            obs = load_reference_obs(dataset, ref)
            if obs is None:
                logger.warning(f"Could not load obs for {dataset}, ref {ref}")
                ref_data[ref] = {'obs': None, 'rare_types': []}
                continue
            
            rare_types = identify_rare_cell_types(dataset, ref, obs)
            ref_data[ref] = {'obs': obs, 'rare_types': rare_types}
            logger.info(f"Dataset {dataset}, ref {ref}: Found {len(rare_types)} rare cell types")
        
        # Create a table for each sample size
        for size in config['sizes']:
            logger.info(f"Creating table for {dataset}, sample size {size}...")
            df = create_sampling_methods_table(dataset, sample_size=size, ref_data=ref_data)
            
            output_file = os.path.join(OUTPUT_DIR, f'{dataset}_sampling_methods_table_size_{size}.csv')
            df.to_csv(output_file)
            logger.info(f"Saved table to {output_file}")
            
            # Print summary for first dataset and first size as example
            if dataset == datasets[0] and size == config['sizes'][0]:
                print("\n" + "="*80)
                print(f"{dataset.upper()} Sampling Methods Comparison Table (Sample Size: {size})")
                print("="*80)
                print("Rows: Sampling Methods")
                print("Columns: References")
                print("Entries: Number of Rare Cells Collected (Rep 0)")
                print("="*80)
                print(df.to_string())
                print("="*80)
    
    logger.info("All tables created successfully.")


if __name__ == "__main__":
    main()

