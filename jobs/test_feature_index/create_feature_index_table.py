#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
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
    log_file = os.path.join(log_dir, f'create_feature_index_table_{timestamp}.log')
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

FEATURE_INDICES = list(range(1, 31))  # 1 to 30
REPS = [0]
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
    return os.path.join(project_root_env, 'jobs', 'test_feature_index', 'tables')

file_path_env = get_data_path()
OUTPUT_DIR = get_output_dir()
RESULTS_DIR = os.path.join(file_path_env, 'test_feature_index')
os.makedirs(OUTPUT_DIR, exist_ok=True)


def get_benchmark_path(dataset):
    """Get benchmark path for a dataset."""
    return file_path_env


def load_reference_obs(dataset, ref):
    """Load reference observation data."""
    benchmark_path = get_benchmark_path(dataset)
    obs_path = os.path.join(benchmark_path, f"{dataset}/benchmark/{ref}/obs.csv")
    if os.path.exists(obs_path):
        df = pd.read_csv(obs_path, index_col=0)
        df[label_key] = df[label_key].astype('category')
        return df
    else:
        adata_path = os.path.join(benchmark_path, f"{dataset}/benchmark/{ref}/adata.h5ad")
        if os.path.exists(adata_path):
            adata = sc.read_h5ad(adata_path)
            return adata.obs
        else:
            logger.warning(f"Could not load obs for {dataset}, ref {ref}")
            return None


def identify_rare_cell_types(dataset, ref, obs):
    """Identify rare cell types for a dataset."""
    benchmark_path = get_benchmark_path(dataset)
    adata_path = os.path.join(benchmark_path, f"{dataset}/benchmark/{ref}/adata.h5ad")
    
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


def create_feature_index_table(dataset, sample_size, ref_data):
    """Create a table with references as columns and feature_index as rows.
    
    Args:
        dataset: Dataset name (e.g., 'lcmv', 'mcc', 'mcc_01', 'mcc_05')
        sample_size: Specific sample size to use (required)
        ref_data: Pre-loaded reference data dictionary with keys as ref numbers
                  and values as dicts with 'obs' and 'rare_types' keys
    
    Returns:
        DataFrame with feature_index as rows and references as columns.
        Entries are the number of rare cells collected.
    """
    config = DATASET_CONFIGS[dataset]
    references = config['references']
    
    # Initialize table
    table_data = []
    
    for feature_index in FEATURE_INDICES:
        row = {'feature_index': feature_index}
        
        for ref in references:
            if ref_data[ref]['obs'] is None:
                row[f'ref_{ref}M'] = 0
                continue
            
            obs = ref_data[ref]['obs']
            rare_types = ref_data[ref]['rare_types']
            
            if len(rare_types) == 0:
                row[f'ref_{ref}M'] = 0
                continue
            
            # Aggregate across reps for this specific sample size
            total_rare_cells = 0
            count = 0
            
            for rep in REPS:
                try:
                    results_path = os.path.join(RESULTS_DIR, dataset, str(ref), 
                                               str(feature_index), str(sample_size), str(rep), 'results.pkl')
                    
                    if not os.path.exists(results_path):
                        logger.debug(f"Results not found: {results_path}")
                        continue
                    
                    with open(results_path, 'rb') as f:
                        results_data = pickle.load(f)
                    
                    cell_indices = results_data[0][0]
                    # # Results are saved as tuple (indices, elapsed_time)
                    # # The indices are in the first element of the tuple
                    # # Based on notebook: results_data is a tuple, results_data[0] is also a tuple of length 2,
                    # # and we need results_data[0][0] to get the actual list of indices
                    # if isinstance(results_data, tuple):
                    #     cell_indices = results_data[0]
                    #     # cell_indices is itself a tuple, get the first element (the actual list of indices)
                    #     # This matches the notebook behavior: obs.iloc[cell_indices[0]]
                    #     cell_indices = cell_indices[0] if isinstance(cell_indices, (tuple, list)) else cell_indices
                    # elif isinstance(results_data, dict):
                    #     cell_indices = results_data.get('indices', results_data)
                    # else:
                    #     cell_indices = results_data
                    
                    # Use iloc with the indices directly (they are positional integer indices)
                    sampled_obs = obs.iloc[cell_indices]
                    
                    # Count rare cells in sample - count by cell type first, then filter to rare types
                    value_counts = sampled_obs[label_key].value_counts()
                    
                    # Filter to only rare cell types and sum their counts
                    # rare_types is a list of cell type names - ensure both are same type for comparison
                    rare_types_set = set(rare_types)
                    sample_cell_types = set(value_counts.index)
                    matching_rare_types = rare_types_set.intersection(sample_cell_types)
                    
                    if len(matching_rare_types) > 0:
                        rare_cell_type_counts = value_counts[value_counts.index.isin(matching_rare_types)]
                        rare_count = rare_cell_type_counts.sum()
                    else:
                        logger.warning(f"No matching rare types found for {dataset}, ref {ref}, FI {feature_index}, size {sample_size}, rep {rep}")
                        print(f"No matching rare types found for {dataset}, ref {ref}, FI {feature_index}, size {sample_size}, rep {rep}")
                        rare_count = 0
                    
                    total_rare_cells += rare_count
                    count += 1
                    
                    if rare_count > 0:
                        logger.info(f"Dataset {dataset}, ref {ref}, FI {feature_index}, size {sample_size}, rep {rep}: Found {rare_count} rare cells (sampled: {len(sampled_obs)}, matching rare types: {list(matching_rare_types)})")
                    elif count == 1:  # Log once per feature_index/ref combo for debugging
                        logger.info(f"Dataset {dataset}, ref {ref}, FI {feature_index}, size {sample_size}, rep {rep}: No rare cells (sampled: {len(sampled_obs)}, rare_types: {list(rare_types)[:3]}, sample_types: {list(sample_cell_types)[:5]})")
                except Exception as e:
                    logger.warning(f"Error processing {dataset}, ref {ref}, feature_index {feature_index}, size {sample_size}, rep {rep}: {str(e)}")
                    print(f"Error processing {dataset}, ref {ref}, feature_index {feature_index}, size {sample_size}, rep {rep}: {str(e)}")
                    continue
            
            # Average across reps
            if count > 0:
                row[f'ref_{ref}M'] = total_rare_cells / count
            else:
                row[f'ref_{ref}M'] = 0

        logger.info("="*80)
        logger.info("DEBUGGING-- new row created")
        logger.info(f"Dataset {dataset}, ref {ref}, FI {feature_index}, size {sample_size}: Total rare cells: {total_rare_cells}, count: {count}")
        logger.info(f"Dataset {dataset}, ref {ref}, FI {feature_index}, size {sample_size}: Row: {row}")
        table_data.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(table_data)
    df.set_index('feature_index', inplace=True)
    
    # Rename columns to be more readable (references as column names)
    df.columns = [f'{ref}M' for ref in references]
    
    return df


def process_single_table(dataset, sample_size):
    """Process a single dataset and sample size combination."""
    logger.info(f"Creating feature_index comparison table for {dataset}, size {sample_size}")
    
    if dataset not in DATASET_CONFIGS:
        logger.error(f"Unknown dataset: {dataset}. Available datasets: {list(DATASET_CONFIGS.keys())}")
        sys.exit(1)
    
    config = DATASET_CONFIGS[dataset]
    references = config['references']
    
    if sample_size not in config['sizes']:
        logger.error(f"Invalid sample size {sample_size} for dataset {dataset}. Available sizes: {config['sizes']}")
        sys.exit(1)
    
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
    
    # Create the table for this sample size
    logger.info(f"Creating table for {dataset}, sample size {sample_size}...")
    df = create_feature_index_table(dataset, sample_size=sample_size, ref_data=ref_data)
    
    output_file = os.path.join(OUTPUT_DIR, f'{dataset}_feature_index_table_size_{sample_size}.csv')
    df.to_csv(output_file)
    logger.info(f"Saved table to {output_file}")
    
    # Print summary
    print("\n" + "="*80)
    print(f"{dataset.upper()} Feature Index Comparison Table (Sample Size: {sample_size})")
    print("="*80)
    print("Rows: Feature Indices (1-30)")
    print("Columns: References")
    print("Entries: Number of Rare Cells Collected")
    print("="*80)
    print(df.to_string())
    print("="*80)
    
    logger.info(f"Successfully created table for {dataset}, size {sample_size}")


def get_all_table_combinations():
    """Generate all dataset and size combinations.
    
    Returns:
        list of tuples: [(dataset, size), ...]
        Order matches the bash script: lcmv, mcc, mcc_01, mcc_05
    """
    combinations = []
    # Use explicit order to match bash script
    datasets = ['lcmv', 'mcc', 'mcc_01', 'mcc_05']
    for dataset in datasets:
        config = DATASET_CONFIGS[dataset]
        for size in config['sizes']:
            combinations.append((dataset, size))
    return combinations


def main():
    parser = argparse.ArgumentParser(description='Create feature index comparison tables')
    parser.add_argument('--dataset', type=str, default=None,
                        choices=['lcmv', 'mcc', 'mcc_01', 'mcc_05'],
                        help='Dataset name (required if --array-index not provided)')
    parser.add_argument('--size', type=int, default=None,
                        help='Sample size (required if --array-index not provided)')
    parser.add_argument('--array-index', type=int, default=None,
                        help='SLURM array task ID (0-indexed) to map to dataset/size combination')
    
    args = parser.parse_args()
    
    # If array-index is provided, use it to get dataset and size
    if args.array_index is not None:
        combinations = get_all_table_combinations()
        if args.array_index < 0 or args.array_index >= len(combinations):
            logger.error(f"Array index {args.array_index} is out of range. Valid range: 0-{len(combinations)-1}")
            sys.exit(1)
        dataset, size = combinations[args.array_index]
        logger.info(f"Array index {args.array_index} maps to dataset={dataset}, size={size}")
    elif args.dataset is not None and args.size is not None:
        dataset = args.dataset
        size = args.size
    else:
        parser.error("Either --array-index OR both --dataset and --size must be provided")
        sys.exit(1)
    
    process_single_table(dataset, size)


if __name__ == "__main__":
    main()

