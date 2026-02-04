#!/usr/bin/env python
# coding: utf-8

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from sparsesampler.sampling import sample
import time
import scipy.sparse
import pickle
import argparse

# Get PROJECT_ROOT from environment or derive from script location
project_root = os.environ.get('PROJECT_ROOT')
if project_root is None:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(script_dir))

file_path_env = os.path.join(project_root, 'data')

# Dataset configurations
DATASET_CONFIGS = {
    'lcmv': {
        'references': [1, 5, 10, 20, 34],
        'sizes': [50000, 100000, 200000],
        'benchmark_path': 'lcmv/benchmark',
    },
    'mcc': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'benchmark_path': 'mcc/benchmark',
    },
    'mcc_01': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'benchmark_path': 'mcc_01/benchmark',
    },
    'mcc_05': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'benchmark_path': 'mcc_05/benchmark',
    }
}

FEATURE_INDICES = list(range(1, 31))  # 1 to 30 inclusive
REPS = [0]  # Only 1 replication
label_key = 'celltype'

# Output directory
OUTPUT_DIR = os.path.join(file_path_env, 'test_feature_index')


def get_data_path(dataset):
    """Get the data path for a given dataset."""
    # All datasets use the same data path
    return file_path_env


def generate_sps_with_feature_index(adata, size, feature_index, seed=1234):
    """Generate SPS samples with a specific feature_index."""
    # Convert sparse matrix to dense if needed
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    
    # Use the optimized function that works with X directly
    return generate_sps_with_feature_index_from_X(X, size, feature_index, seed)


def save_results(dataset, ref, feature_index, size, rep, indices, elapsed_time):
    """Save results to disk incrementally.
    
    The results are saved in the same format as the original parallel.py scripts:
    - Tuple: (indices, elapsed_time)
    - This matches what the plotting scripts expect
    """
    output_dir = os.path.join(OUTPUT_DIR, dataset, str(ref), str(feature_index), str(size), str(rep))
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, 'results.pkl')
    
    # Save in the same format as original scripts: (indices, elapsed_time)
    save_data = (indices, elapsed_time)
    
    with open(output_file, 'wb') as f:
        pickle.dump(save_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    print(f"Saved results to {output_file}")


def process_single_combination(dataset, ref, feature_index, size, rep, seed):
    """Process a single combination of parameters."""
    print(f"\n{'='*60}")
    print(f"Processing: dataset={dataset}, ref={ref}, feature_index={feature_index}, size={size}, rep={rep}")
    print(f"{'='*60}")
    
    # Check if already processed
    output_file = os.path.join(OUTPUT_DIR, dataset, str(ref), str(feature_index), str(size), str(rep), 'results.pkl')
    if os.path.exists(output_file):
        print(f"Results already exist, skipping: {output_file}")
        return
    
    # Get data path
    data_path = get_data_path(dataset)
    config = DATASET_CONFIGS[dataset]
    benchmark_path = os.path.join(data_path, config['benchmark_path'])
    
    # Load adata
    adata_path = os.path.join(benchmark_path, f"{ref}/adata.h5ad")
    if not os.path.exists(adata_path):
        print(f"Warning: adata file not found: {adata_path}")
        return
    
    print(f"Loading data from: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    adata.obs[label_key] = adata.obs[label_key].astype('category')
    adata.var.index = adata.var.index.astype('object')
    
    # Generate samples
    try:
        indices, elapsed_time = generate_sps_with_feature_index(adata, size, feature_index, seed)
        
        # Save results incrementally
        save_results(dataset, ref, feature_index, size, rep, indices, elapsed_time)
        
        print(f"Successfully completed: dataset={dataset}, ref={ref}, feature_index={feature_index}, size={size}, rep={rep}")
        
    except Exception as e:
        print(f"Error processing {dataset}, ref={ref}, feature_index={feature_index}, size={size}, rep={rep}: {str(e)}")
        import traceback
        traceback.print_exc()
        raise


def process_all_for_dataset_ref(dataset, ref, seed_base=10000):
    """Process all feature_index and size combinations for a given dataset and reference.
    This is more efficient as it loads the data once."""
    print(f"\n{'#'*60}")
    print(f"Processing all combinations for dataset={dataset}, ref={ref}")
    print(f"{'#'*60}")
    
    # Get data path
    data_path = get_data_path(dataset)
    config = DATASET_CONFIGS[dataset]
    benchmark_path = os.path.join(data_path, config['benchmark_path'])
    
    # Load adata once
    adata_path = os.path.join(benchmark_path, f"{ref}/adata.h5ad")
    if not os.path.exists(adata_path):
        print(f"Warning: adata file not found: {adata_path}")
        return
    
    print(f"Loading data from: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    adata.obs[label_key] = adata.obs[label_key].astype('category')
    adata.var.index = adata.var.index.astype('object')
    
    # Convert sparse matrix to dense if needed (once)
    if scipy.sparse.issparse(adata.X):
        print("Converting sparse matrix to dense...")
        adata.X = adata.X.toarray()
    
    # Extract X once to avoid repeated access
    X = adata.X
    
    # Process all combinations
    for feature_index in FEATURE_INDICES:
        for size in config['sizes']:
            for rep in REPS:
                seed = seed_base + ref * 1000 + feature_index * 10 + size // 1000 + rep
                
                # Check if already processed
                output_file = os.path.join(OUTPUT_DIR, dataset, str(ref), str(feature_index), str(size), str(rep), 'results.pkl')
                if os.path.exists(output_file):
                    print(f"Skipping (already exists): dataset={dataset}, ref={ref}, feature_index={feature_index}, size={size}, rep={rep}")
                    continue
                
                print(f"\nProcessing: feature_index={feature_index}, size={size}, rep={rep}")
                
                try:
                    # Generate samples (pass X directly to avoid repeated access)
                    indices, elapsed_time = generate_sps_with_feature_index_from_X(X, size, feature_index, seed)
                    
                    # Save results incrementally
                    save_results(dataset, ref, feature_index, size, rep, indices, elapsed_time)
                    
                except Exception as e:
                    print(f"Error processing feature_index={feature_index}, size={size}, rep={rep}: {str(e)}")
                    import traceback
                    traceback.print_exc()
                    # Continue with next combination instead of failing completely
                    continue
    
    # Explicitly free memory
    del adata
    del X
    import gc
    gc.collect()
    print(f"Completed processing for dataset={dataset}, ref={ref}. Memory freed.")


def generate_sps_with_feature_index_from_X(X, size, feature_index, seed=1234):
    """Generate SPS samples with a specific feature_index, using pre-loaded X matrix.
    This avoids repeated data loading/conversion."""
    print(f'Running SPS with feature_index={feature_index}, size={size}, seed={seed}')
    start_time = time.time()
    
    # Use auto_k=False with feature_index parameter (direct parameter in sample function)
    # k parameter is required when auto_k=False
    k = size // 100
    result = sample(X=X, size=size, seed=seed, auto_k=False, k=k, feature_index=feature_index)
    
    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")
    
    return result, elapsed_time


def main():
    parser = argparse.ArgumentParser(description="Test feature_index parameter for SPS sampling")
    parser.add_argument("--dataset", type=str, required=True, 
                       choices=['lcmv', 'mcc', 'mcc_01', 'mcc_05'],
                       help="Dataset to process")
    parser.add_argument("--ref", type=int, required=True, help="Reference to process")
    parser.add_argument("--feature_index", type=int, default=None, 
                       help="Specific feature_index to test (if None, tests all 10-30)")
    parser.add_argument("--size", type=int, default=None,
                       help="Specific size to test (if None, tests all sizes)")
    parser.add_argument("--rep", type=int, default=0, help="Replicate number")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--all", action='store_true', 
                       help="Process all feature_index and size combinations for this dataset/ref (more efficient)")
    
    args = parser.parse_args()
    
    print("###############################")
    print("************ New run **********")
    print(f"Dataset: {args.dataset}")
    print(f"Reference: {args.ref}")
    print(f"Node: {os.uname().nodename}")
    print("###############################")
    
    # Validate dataset
    if args.dataset not in DATASET_CONFIGS:
        print(f"Error: Unknown dataset {args.dataset}")
        return
    
    config = DATASET_CONFIGS[args.dataset]
    
    # Validate ref
    if args.ref not in config['references']:
        print(f"Error: Reference {args.ref} not in valid references {config['references']}")
        return
    
    # Set seed
    if args.seed is None:
        seed_base = 10000 + args.dataset.__hash__() % 1000
    else:
        seed_base = args.seed
    
    # Process based on arguments
    if args.all:
        # Process all combinations for this dataset/ref (efficient mode - loads data once)
        process_all_for_dataset_ref(args.dataset, args.ref, seed_base)
    elif args.feature_index is not None and args.size is not None:
        # Process single combination
        seed = seed_base + args.ref * 1000 + args.feature_index * 10 + args.size // 1000 + args.rep
        process_single_combination(args.dataset, args.ref, args.feature_index, args.size, args.rep, seed)
    else:
        print("Error: Use --all flag or specify both --feature_index and --size")
        return


if __name__ == "__main__":
    main()

