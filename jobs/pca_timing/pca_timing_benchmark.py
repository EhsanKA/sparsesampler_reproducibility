#!/usr/bin/env python
# coding: utf-8
"""
PCA Timing Benchmark for all datasets
Measures PCA computation time using the sparsesampler pipeline
"""

import os
import sys
import numpy as np
import scanpy as sc
import scipy.sparse
import time
from sparsesampler.sampling import sample

# Get PROJECT_ROOT from environment or derive from script location
project_root = os.environ.get('PROJECT_ROOT')
if project_root is None:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(script_dir))

data_path = os.path.join(project_root, 'data')

# Dataset configurations - using FULL datasets (largest reference)
DATASETS = {
    'lcmv': {
        'path': 'lcmv/benchmark/34/adata.h5ad',  # 34M cells (full)
        'label_key': 'celltype'
    },
    'mcc': {
        'path': 'mcc/benchmark/30/adata.h5ad',  # 3M cells (full)
        'label_key': 'celltype'
    },
    'mcc_01': {
        'path': 'mcc_01/benchmark/30/adata.h5ad',  # 3M cells (full)
        'label_key': 'celltype'
    },
    'mcc_05': {
        'path': 'mcc_05/benchmark/30/adata.h5ad',  # 3M cells (full)
        'label_key': 'celltype'
    },
}


def measure_pca_time(X, size=50000, seed=1234, n_runs=3):
    """
    Measure PCA time by running the sample function and capturing timing.
    The sparsesampler prints PCA time internally.
    We'll run it multiple times and return the average.
    """
    import io
    import re
    from contextlib import redirect_stdout
    
    pca_times = []
    total_times = []
    
    for i in range(n_runs):
        # Capture stdout to get PCA time
        f = io.StringIO()
        with redirect_stdout(f):
            start = time.time()
            result = sample(X=X, size=min(size, X.shape[0]-1), seed=seed+i, 
                          auto_k=False, k=None, feature_index=18)
            total_time = time.time() - start
        
        output = f.getvalue()
        
        # Parse PCA time from output
        pca_match = re.search(r'Elapsed time after PCA: ([\d.]+)', output)
        if pca_match:
            pca_times.append(float(pca_match.group(1)))
        
        total_times.append(total_time)
    
    avg_pca_time = np.mean(pca_times) if pca_times else None
    avg_total_time = np.mean(total_times)
    
    return avg_pca_time, avg_total_time


def main():
    print("=" * 80)
    print("PCA Timing Benchmark")
    print("=" * 80)
    print(f"Node: {os.uname().nodename}")
    print()
    
    results = []
    
    for dataset_name, config in DATASETS.items():
        adata_path = os.path.join(data_path, config['path'])
        
        if not os.path.exists(adata_path):
            print(f"WARNING: {adata_path} not found, skipping {dataset_name}")
            continue
        
        print(f"\nLoading {dataset_name}...")
        adata = sc.read_h5ad(adata_path)
        
        # Convert sparse to dense if needed
        X = adata.X
        if scipy.sparse.issparse(X):
            X = X.toarray()
        
        n_cells = X.shape[0]
        n_features = X.shape[1]
        
        print(f"  n_cells: {n_cells:,}")
        print(f"  n_features: {n_features:,}")
        print(f"  Running PCA timing (3 runs)...")
        
        pca_time, total_time = measure_pca_time(X, size=50000, n_runs=3)
        
        results.append({
            'dataset': dataset_name,
            'n_cells': n_cells,
            'n_features': n_features,
            'pca_time': pca_time,
            'total_time': total_time
        })
        
        print(f"  PCA time: {pca_time:.2f}s" if pca_time else "  PCA time: N/A")
        print(f"  Total time: {total_time:.2f}s")
        
        # Free memory
        del adata, X
        import gc
        gc.collect()
    
    # Print summary table
    print("\n" + "=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)
    print(f"{'Dataset':<10} | {'n_cells':>12} | {'n_features':>12} | {'PCA Time (s)':>14}")
    print("-" * 60)
    
    for r in results:
        pca_str = f"{r['pca_time']:.2f}" if r['pca_time'] else "N/A"
        print(f"{r['dataset']:<10} | {r['n_cells']:>12,} | {r['n_features']:>12,} | {pca_str:>14}")
    
    print("=" * 80)
    
    # Save results to CSV
    output_dir = os.path.join(project_root, 'jobs', 'pca_timing')
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'pca_timing_results.csv')
    
    import csv
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['dataset', 'n_cells', 'n_features', 'pca_time', 'total_time'])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
