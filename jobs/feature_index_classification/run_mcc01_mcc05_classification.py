#!/usr/bin/env python
# coding: utf-8

"""
Run RF classification for all feature indices for mcc_01 and mcc_05 datasets.

This script processes all feature indices sequentially.

Usage:
    python run_mcc01_mcc05_classification.py --dataset mcc_01
    python run_mcc01_mcc05_classification.py --dataset mcc_05
    python run_mcc01_mcc05_classification.py --dataset all
"""

import os
import sys
import argparse
import subprocess
from datetime import datetime

PROJECT_ROOT = "/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
SCRIPT_DIR = os.path.join(PROJECT_ROOT, "jobs/feature_index_classification")

def run_classification(dataset, feature_indices):
    """Run classification for all feature indices."""
    
    script_path = os.path.join(SCRIPT_DIR, "classify_by_feature_index_unified.py")
    
    print(f"\n{'='*70}")
    print(f"Running classification for {dataset}")
    print(f"Feature indices: {feature_indices[0]}-{feature_indices[-1]}")
    print(f"Start time: {datetime.now()}")
    print(f"{'='*70}\n")
    
    for fi in feature_indices:
        print(f"\n--- Processing feature index {fi} ---")
        
        cmd = [
            "python", script_path,
            "--dataset", dataset,
            "--feature-index", str(fi),
            "--size", "100000",
            "--rep", "0"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=False, text=True)
            if result.returncode != 0:
                print(f"Warning: Feature index {fi} failed with return code {result.returncode}")
        except Exception as e:
            print(f"Error processing feature index {fi}: {e}")
            continue
    
    print(f"\n{'='*70}")
    print(f"Completed {dataset}")
    print(f"End time: {datetime.now()}")
    print(f"{'='*70}\n")


def combine_results(dataset):
    """Combine all results for a dataset into a single CSV."""
    import pandas as pd
    import glob
    
    # Set output directory based on dataset
    if dataset == 'mcc':
        results_dir = os.path.join(SCRIPT_DIR, "results")
    else:
        results_dir = os.path.join(SCRIPT_DIR, f"results_{dataset}")
    
    # Find all summary CSV files
    csv_files = sorted(glob.glob(os.path.join(results_dir, "summary_feature_index_*.csv")))
    
    if not csv_files:
        print(f"No results found for {dataset} in {results_dir}")
        return None
    
    print(f"Found {len(csv_files)} result files for {dataset}")
    
    # Combine all results
    dfs = []
    for f in csv_files:
        df = pd.read_csv(f)
        dfs.append(df)
    
    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.sort_values('feature_index')
    
    # Save combined results
    combined_path = os.path.join(results_dir, "all_feature_indices_summary.csv")
    combined.to_csv(combined_path, index=False)
    print(f"Combined results saved to: {combined_path}")
    
    return combined


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run RF classification for mcc_01 and mcc_05'
    )
    parser.add_argument(
        '--dataset',
        type=str,
        required=True,
        choices=['mcc_01', 'mcc_05', 'all'],
        help='Dataset to process (mcc_01, mcc_05, or all)'
    )
    parser.add_argument(
        '--combine-only',
        action='store_true',
        help='Only combine existing results, do not run classification'
    )
    parser.add_argument(
        '--start-index',
        type=int,
        default=1,
        help='Starting feature index (default: 1)'
    )
    parser.add_argument(
        '--end-index',
        type=int,
        default=30,
        help='Ending feature index (default: 30)'
    )
    
    args = parser.parse_args()
    
    feature_indices = list(range(args.start_index, args.end_index + 1))
    
    datasets = ['mcc_01', 'mcc_05'] if args.dataset == 'all' else [args.dataset]
    
    for dataset in datasets:
        if not args.combine_only:
            run_classification(dataset, feature_indices)
        
        # Combine results
        combine_results(dataset)
    
    print("\nAll done!")
