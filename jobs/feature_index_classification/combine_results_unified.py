#!/usr/bin/env python
# coding: utf-8

"""
Combine and analyze results from all feature index classifications for all datasets.

Usage:
    python combine_results_unified.py
    python combine_results_unified.py --dataset mcc_01
"""

import os
import pandas as pd
import numpy as np
import glob
import argparse

# Paths
PROJECT_ROOT = "/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
SCRIPT_DIR = os.path.join(PROJECT_ROOT, "jobs/feature_index_classification")

DATASETS = {
    'mcc': {
        'results_dir': os.path.join(SCRIPT_DIR, 'results'),
        'rare_cell_type': 'osteoblast',
        'label': 'MCC (3.2M, 0.95% rare)',
    },
    'mcc_01': {
        'results_dir': os.path.join(SCRIPT_DIR, 'results_mcc_01'),
        'rare_cell_type': 'osteoblast',
        'label': 'MCC_01 (3.1M, 0.1% rare)',
    },
    'mcc_05': {
        'results_dir': os.path.join(SCRIPT_DIR, 'results_mcc_05'),
        'rare_cell_type': 'osteoblast',
        'label': 'MCC_05 (3.1M, 0.5% rare)',
    },
}


def combine_results(dataset):
    """Combine all results for a dataset into a single CSV."""
    
    config = DATASETS[dataset]
    results_dir = config['results_dir']
    
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


def print_summary(dataset, df, rare_cell_type):
    """Print summary of results."""
    
    print(f"\n" + "=" * 100)
    print(f"{dataset.upper()} FEATURE INDEX CLASSIFICATION RESULTS")
    print("=" * 100)
    
    # Check column names for rare cell type improvement
    rare_imp_col = f'{rare_cell_type}_f1_improvement'
    if rare_imp_col not in df.columns:
        # Try with different column name pattern
        for col in df.columns:
            if rare_cell_type in col and 'improvement' in col:
                rare_imp_col = col
                break
    
    print(f"\n{'FI':>4} {'EVR':>10} {'SPS Acc':>10} {'Rand Acc':>10} {'Δ Acc':>8} "
          f"{'SPS F1':>10} {'Rand F1':>10} {'Δ F1':>8} "
          f"{'SPS Rare':>10} {'Rand Rare':>10} {'Δ Rare':>8}")
    print("-" * 100)
    
    for _, row in df.iterrows():
        rare_sps_col = f'sps_{rare_cell_type}_f1_mean'
        rare_rand_col = f'random_{rare_cell_type}_f1_mean'
        
        print(f"{int(row['feature_index']):>4} {row['evr']:>10.2f} "
              f"{row['sps_accuracy_mean']:>10.4f} {row['random_accuracy_mean']:>10.4f} "
              f"{row['accuracy_improvement']:>+8.2f}% "
              f"{row['sps_macro_f1_mean']:>10.4f} {row['random_macro_f1_mean']:>10.4f} "
              f"{row['macro_f1_improvement']:>+8.2f}% "
              f"{row.get(rare_sps_col, 0):>10.4f} {row.get(rare_rand_col, 0):>10.4f} "
              f"{row.get(rare_imp_col, 0):>+8.2f}%")
    
    # Find best indices
    print("\n" + "-" * 80)
    print("BEST FEATURE INDICES:")
    print("-" * 80)
    
    # Best by rare cell type F1 improvement
    if rare_imp_col in df.columns:
        best_rare = df.loc[df[rare_imp_col].idxmax()]
        print(f"  Best {rare_cell_type} F1 improvement: FI {int(best_rare['feature_index'])} "
              f"(Δ={best_rare[rare_imp_col]:+.2f}%)")
    
    # Best by macro F1
    best_f1 = df.loc[df['sps_macro_f1_mean'].idxmax()]
    print(f"  Best SPS Macro F1: FI {int(best_f1['feature_index'])} "
          f"(F1={best_f1['sps_macro_f1_mean']:.4f})")


def main(datasets_to_process=None):
    """Process all specified datasets."""
    
    if datasets_to_process is None:
        datasets_to_process = list(DATASETS.keys())
    
    all_results = {}
    
    for dataset in datasets_to_process:
        if dataset not in DATASETS:
            print(f"Unknown dataset: {dataset}")
            continue
        
        config = DATASETS[dataset]
        df = combine_results(dataset)
        
        if df is not None:
            all_results[dataset] = df
            print_summary(dataset, df, config['rare_cell_type'])
    
    return all_results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Combine and analyze feature index classification results'
    )
    parser.add_argument(
        '--dataset',
        type=str,
        nargs='+',
        choices=list(DATASETS.keys()) + ['all'],
        default=['all'],
        help='Dataset(s) to process'
    )
    
    args = parser.parse_args()
    
    if 'all' in args.dataset:
        datasets = list(DATASETS.keys())
    else:
        datasets = args.dataset
    
    main(datasets)
