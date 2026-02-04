#!/usr/bin/env python
"""
Analyze optimal EVR index across different sample sizes and datasets.
Evaluate if EVR index 12 is consistently optimal and assess its stability.
"""

import os
import pandas as pd
import numpy as np

# Paths
TABLES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tables')

# Dataset configurations
DATASETS = {
    'lcmv': {'sizes': [50000, 100000, 200000], 'label': 'LCMV'},
    'mcc': {'sizes': [50000, 100000, 200000], 'label': 'MCC (1%)'},
    'mcc_05': {'sizes': [50000, 100000, 200000], 'label': 'MCC (0.5%)'},
    'mcc_01': {'sizes': [50000, 100000, 200000], 'label': 'MCC (0.1%)'},
}


def load_table(dataset, size):
    """Load a feature index table."""
    path = os.path.join(TABLES_DIR, f'{dataset}_feature_index_table_size_{size}.csv')
    if os.path.exists(path):
        return pd.read_csv(path)
    return None


def find_optimal_evr(df, exclude_evr2=True):
    """Find the optimal EVR index (highest mean coverage across references).
    
    Args:
        df: DataFrame with feature_index column and reference columns
        exclude_evr2: If True, exclude EVR=2 which often has anomalously high values
    """
    # Get reference columns (all except feature_index)
    ref_cols = [c for c in df.columns if c != 'feature_index']
    
    # Calculate mean across references for each EVR
    df_analysis = df.copy()
    df_analysis['mean_coverage'] = df_analysis[ref_cols].mean(axis=1)
    
    if exclude_evr2:
        # Exclude EVR=2 from consideration (often anomalous)
        df_filtered = df_analysis[df_analysis['feature_index'] != 2]
    else:
        df_filtered = df_analysis
    
    # Find the EVR with maximum mean coverage
    optimal_idx = df_filtered['mean_coverage'].idxmax()
    optimal_evr = df_filtered.loc[optimal_idx, 'feature_index']
    optimal_coverage = df_filtered.loc[optimal_idx, 'mean_coverage']
    
    return int(optimal_evr), optimal_coverage


def get_evr_performance(df, evr_index):
    """Get the mean coverage for a specific EVR index."""
    ref_cols = [c for c in df.columns if c != 'feature_index']
    row = df[df['feature_index'] == evr_index]
    if len(row) == 0:
        return None
    return row[ref_cols].values.mean()


def analyze_evr_stability():
    """Analyze EVR index stability across datasets and sample sizes."""
    
    print("=" * 80)
    print("EVR Index Optimal Analysis Across Sample Sizes")
    print("=" * 80)
    
    # Store results
    all_results = []
    optimal_evrs = []
    
    for dataset, config in DATASETS.items():
        print(f"\n{'='*60}")
        print(f"Dataset: {config['label']}")
        print("="*60)
        
        for size in config['sizes']:
            df = load_table(dataset, size)
            if df is None:
                print(f"  Size {size:,}: No data")
                continue
            
            optimal_evr, optimal_cov = find_optimal_evr(df, exclude_evr2=True)
            
            # Also get performance at EVR=12 and EVR=18 (default)
            evr12_cov = get_evr_performance(df, 12)
            evr18_cov = get_evr_performance(df, 18)
            
            # Calculate relative performance of EVR 12 vs optimal
            if evr12_cov and optimal_cov:
                evr12_rel = (evr12_cov / optimal_cov) * 100
            else:
                evr12_rel = None
                
            print(f"\n  Sample Size: {size:,}")
            print(f"    Optimal EVR: {optimal_evr} (mean coverage: {optimal_cov:.1f})")
            print(f"    EVR 12 coverage: {evr12_cov:.1f} ({evr12_rel:.1f}% of optimal)" if evr12_cov else "    EVR 12: N/A")
            print(f"    EVR 18 coverage: {evr18_cov:.1f}" if evr18_cov else "    EVR 18: N/A")
            
            all_results.append({
                'dataset': dataset,
                'label': config['label'],
                'size': size,
                'optimal_evr': optimal_evr,
                'optimal_coverage': optimal_cov,
                'evr12_coverage': evr12_cov,
                'evr12_relative': evr12_rel,
                'evr18_coverage': evr18_cov,
            })
            optimal_evrs.append(optimal_evr)
    
    # Summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY: Optimal EVR Index Distribution")
    print("=" * 80)
    
    # Count frequency of each optimal EVR
    evr_counts = pd.Series(optimal_evrs).value_counts().sort_index()
    print("\nOptimal EVR frequency across all conditions:")
    for evr, count in evr_counts.items():
        pct = count / len(optimal_evrs) * 100
        print(f"  EVR {evr}: {count} times ({pct:.1f}%)")
    
    # Overall statistics
    print(f"\nMost frequent optimal EVR: {evr_counts.idxmax()} (appeared {evr_counts.max()} times)")
    print(f"EVR range of optimal values: {min(optimal_evrs)} - {max(optimal_evrs)}")
    
    # EVR 12 performance summary
    print("\n" + "=" * 80)
    print("EVR 12 Performance Analysis")
    print("=" * 80)
    
    results_df = pd.DataFrame(all_results)
    evr12_relative = results_df['evr12_relative'].dropna()
    
    print(f"\nEVR 12 relative performance (% of optimal):")
    print(f"  Mean: {evr12_relative.mean():.1f}%")
    print(f"  Min:  {evr12_relative.min():.1f}%")
    print(f"  Max:  {evr12_relative.max():.1f}%")
    print(f"  Std:  {evr12_relative.std():.1f}%")
    
    # Check consistency by dataset
    print("\n" + "=" * 80)
    print("Optimal EVR by Dataset (across sample sizes)")
    print("=" * 80)
    
    for dataset in DATASETS.keys():
        dataset_results = results_df[results_df['dataset'] == dataset]
        evrs = dataset_results['optimal_evr'].tolist()
        label = DATASETS[dataset]['label']
        print(f"\n{label}:")
        print(f"  Optimal EVRs across sizes: {evrs}")
        print(f"  Consistent: {'Yes' if len(set(evrs)) == 1 else 'No'}")
    
    # Check consistency by sample size
    print("\n" + "=" * 80)
    print("Optimal EVR by Sample Size (across datasets)")
    print("=" * 80)
    
    for size in [50000, 100000, 200000]:
        size_results = results_df[results_df['size'] == size]
        evrs = size_results['optimal_evr'].tolist()
        print(f"\n{size:,} samples:")
        print(f"  Optimal EVRs across datasets: {evrs}")
        print(f"  Most common: {pd.Series(evrs).mode().iloc[0] if len(evrs) > 0 else 'N/A'}")
    
    # Return results for further analysis
    return results_df


def analyze_top_k_evrs(k=5):
    """Analyze the top-K performing EVR indices for each condition."""
    
    print("\n" + "=" * 80)
    print(f"Top-{k} EVR Indices Analysis")
    print("=" * 80)
    
    top_evrs_all = []
    
    for dataset, config in DATASETS.items():
        for size in config['sizes']:
            df = load_table(dataset, size)
            if df is None:
                continue
            
            ref_cols = [c for c in df.columns if c != 'feature_index']
            df_analysis = df.copy()
            df_analysis['mean_coverage'] = df_analysis[ref_cols].mean(axis=1)
            
            # Exclude EVR=2
            df_filtered = df_analysis[df_analysis['feature_index'] != 2]
            
            # Get top-k EVRs
            top_k = df_filtered.nlargest(k, 'mean_coverage')['feature_index'].tolist()
            top_evrs_all.extend(top_k)
            
            print(f"\n{config['label']} (n={size:,}):")
            print(f"  Top-{k} EVRs: {top_k}")
    
    # Overall frequency in top-k
    print("\n" + "=" * 80)
    print(f"Frequency of EVR indices appearing in top-{k}")
    print("=" * 80)
    
    evr_freq = pd.Series(top_evrs_all).value_counts().sort_values(ascending=False)
    total_conditions = len(DATASETS) * 3  # 4 datasets x 3 sample sizes
    
    print(f"\nTotal conditions analyzed: {total_conditions}")
    print(f"\nEVR index frequency in top-{k}:")
    for evr, count in evr_freq.head(15).items():
        pct = count / total_conditions * 100
        print(f"  EVR {evr:2d}: {count:2d} appearances ({pct:5.1f}%)")
    
    return evr_freq


if __name__ == "__main__":
    results_df = analyze_evr_stability()
    print("\n")
    top_k_freq = analyze_top_k_evrs(k=5)
    
    # Final recommendation
    print("\n" + "=" * 80)
    print("RECOMMENDATION")
    print("=" * 80)
    
    # Find the most robust EVR
    most_common_optimal = results_df['optimal_evr'].mode().iloc[0]
    evr12_mean_rel = results_df['evr12_relative'].mean()
    
    print(f"\n1. Most frequently optimal EVR: {most_common_optimal}")
    print(f"2. EVR 12 achieves {evr12_mean_rel:.1f}% of optimal coverage on average")
    
    if evr12_mean_rel >= 90:
        print(f"\n→ EVR 12 is a GOOD choice (within 10% of optimal across conditions)")
    elif evr12_mean_rel >= 80:
        print(f"\n→ EVR 12 is ACCEPTABLE but not optimal")
    else:
        print(f"\n→ EVR 12 may NOT be the best choice; consider EVR {most_common_optimal}")
