#!/usr/bin/env python
"""
Analyze EVR index performance to recommend optimal default value.

This script provides quantitative analysis to support the choice of 
default EVR index for the paper.
"""

import os
import pandas as pd
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.join(script_dir, 'tables')

# Total rare cells per reference (from obs.csv counts)
TOTAL_RARE_CELLS = {
    'lcmv': {'1M': 20226, '5M': 102561, '10M': 205573, '20M': 410544, '34M': 705738},
    'mcc': {'5M': 4750, '10M': 9399, '20M': 18884, '25M': 23891, '30M': 30000},
    'mcc_05': {'5M': 2396, '10M': 4652, '20M': 9490, '25M': 11984, '30M': 15000},
    'mcc_01': {'5M': 458, '10M': 929, '20M': 1912, '25M': 2377, '30M': 3000},
}

def load_and_analyze():
    """Load all tables and compute statistics."""
    datasets = ['lcmv', 'mcc', 'mcc_05', 'mcc_01']
    all_data = {}
    
    for dataset in datasets:
        table_path = os.path.join(TABLES_DIR, f'{dataset}_feature_index_table_size_100000.csv')
        if os.path.exists(table_path):
            df = pd.read_csv(table_path)
            all_data[dataset] = df
    
    return all_data

def compute_coverage_stats(all_data):
    """Compute coverage percentage statistics for each EVR index."""
    results = []
    
    for evr_idx in range(1, 31):
        row = {'evr_index': evr_idx}
        coverages = []
        
        for dataset, df in all_data.items():
            if evr_idx in df['feature_index'].values:
                evr_row = df[df['feature_index'] == evr_idx].iloc[0]
                
                for ref in TOTAL_RARE_CELLS[dataset].keys():
                    if ref in evr_row.index:
                        captured = evr_row[ref]
                        total = TOTAL_RARE_CELLS[dataset][ref]
                        coverage = (captured / total) * 100
                        coverages.append(coverage)
                        row[f'{dataset}_{ref}_coverage'] = coverage
        
        if coverages:
            row['mean_coverage'] = np.mean(coverages)
            row['std_coverage'] = np.std(coverages)
            row['min_coverage'] = np.min(coverages)
            row['max_coverage'] = np.max(coverages)
            row['cv'] = row['std_coverage'] / row['mean_coverage'] if row['mean_coverage'] > 0 else 0
        
        results.append(row)
    
    return pd.DataFrame(results)

def compute_raw_count_stats(all_data):
    """Compute raw count statistics for each EVR index."""
    results = []
    
    for evr_idx in range(1, 31):
        row = {'evr_index': evr_idx}
        all_counts = []
        
        for dataset, df in all_data.items():
            if evr_idx in df['feature_index'].values:
                evr_row = df[df['feature_index'] == evr_idx].iloc[0]
                
                for col in df.columns:
                    if col != 'feature_index' and col in evr_row.index:
                        count = evr_row[col]
                        all_counts.append(count)
                        row[f'{dataset}_{col}_count'] = count
        
        if all_counts:
            row['mean_count'] = np.mean(all_counts)
            row['total_count'] = np.sum(all_counts)
        
        results.append(row)
    
    return pd.DataFrame(results)

def find_optimal_evr():
    """Find optimal EVR index based on multiple criteria."""
    all_data = load_and_analyze()
    
    coverage_df = compute_coverage_stats(all_data)
    count_df = compute_raw_count_stats(all_data)
    
    print("=" * 80)
    print("EVR INDEX OPTIMIZATION ANALYSIS")
    print("=" * 80)
    
    # Find peaks and optimal points
    print("\n1. COVERAGE STATISTICS (% of rare cells captured)")
    print("-" * 60)
    print(f"{'EVR':<5} {'Mean%':>8} {'Std%':>8} {'Min%':>8} {'Max%':>8} {'CV':>8}")
    print("-" * 60)
    
    for _, row in coverage_df.iterrows():
        evr = int(row['evr_index'])
        mean_cov = row.get('mean_coverage', 0)
        std_cov = row.get('std_coverage', 0)
        min_cov = row.get('min_coverage', 0)
        max_cov = row.get('max_coverage', 0)
        cv = row.get('cv', 0)
        
        # Highlight notable values
        marker = ""
        if evr == 10:
            marker = " <-- Peak LCMV coverage"
        elif evr == 12:
            marker = " <-- High overall coverage"
        elif evr == 3 or evr == 4:
            marker = " <-- Peak MCC coverage"
        elif evr == 18:
            marker = " <-- Current default"
        
        print(f"{evr:<5} {mean_cov:>8.1f} {std_cov:>8.1f} {min_cov:>8.1f} {max_cov:>8.1f} {cv:>8.2f}{marker}")
    
    # Find best indices by criteria
    print("\n" + "=" * 80)
    print("2. TOP 5 EVR INDICES BY CRITERIA")
    print("=" * 80)
    
    # By mean coverage
    top_by_mean = coverage_df.nlargest(5, 'mean_coverage')[['evr_index', 'mean_coverage', 'std_coverage', 'min_coverage']]
    print("\nTop 5 by MEAN coverage:")
    print(top_by_mean.to_string(index=False))
    
    # By minimum coverage (worst-case)
    top_by_min = coverage_df.nlargest(5, 'min_coverage')[['evr_index', 'min_coverage', 'mean_coverage', 'std_coverage']]
    print("\nTop 5 by MINIMUM coverage (best worst-case):")
    print(top_by_min.to_string(index=False))
    
    # By stability (lowest CV)
    # Filter to reasonable coverage first (>5%)
    stable_df = coverage_df[coverage_df['mean_coverage'] > 5].nsmallest(5, 'cv')[['evr_index', 'cv', 'mean_coverage', 'std_coverage']]
    print("\nTop 5 by STABILITY (lowest CV, mean coverage > 5%):")
    print(stable_df.to_string(index=False))
    
    # Dataset-specific analysis
    print("\n" + "=" * 80)
    print("3. DATASET-SPECIFIC OPTIMAL EVR INDICES")
    print("=" * 80)
    
    for dataset, df in all_data.items():
        total_rare = TOTAL_RARE_CELLS[dataset]
        ref_cols = [col for col in df.columns if col != 'feature_index']
        
        # Calculate mean coverage across references for this dataset
        dataset_coverages = []
        for evr_idx in range(1, 31):
            if evr_idx in df['feature_index'].values:
                evr_row = df[df['feature_index'] == evr_idx].iloc[0]
                coverages = []
                for ref in ref_cols:
                    if ref in total_rare:
                        coverages.append((evr_row[ref] / total_rare[ref]) * 100)
                if coverages:
                    dataset_coverages.append({
                        'evr': evr_idx,
                        'mean_coverage': np.mean(coverages),
                        'min_coverage': np.min(coverages)
                    })
        
        ds_df = pd.DataFrame(dataset_coverages)
        best_mean = ds_df.loc[ds_df['mean_coverage'].idxmax()]
        best_min = ds_df.loc[ds_df['min_coverage'].idxmax()]
        
        print(f"\n{dataset.upper()}:")
        print(f"  Best by mean coverage: EVR {int(best_mean['evr'])} ({best_mean['mean_coverage']:.1f}%)")
        print(f"  Best by min coverage:  EVR {int(best_min['evr'])} ({best_min['min_coverage']:.1f}%)")
    
    # Recommendation
    print("\n" + "=" * 80)
    print("4. RECOMMENDATION")
    print("=" * 80)
    
    # Find EVR that balances mean coverage and stability
    # Score = mean_coverage - penalty*CV (penalize high variability)
    coverage_df['score'] = coverage_df['mean_coverage'] - 2 * coverage_df['cv'] * coverage_df['mean_coverage']
    
    # Focus on EVR 5-20 range (reasonable PCA components)
    reasonable_range = coverage_df[(coverage_df['evr_index'] >= 5) & (coverage_df['evr_index'] <= 20)]
    best_balanced = reasonable_range.loc[reasonable_range['score'].idxmax()]
    
    print(f"\nBased on balanced analysis (mean coverage - variability penalty):")
    print(f"  Recommended EVR index: {int(best_balanced['evr_index'])}")
    print(f"  Mean coverage: {best_balanced['mean_coverage']:.1f}%")
    print(f"  Min coverage:  {best_balanced['min_coverage']:.1f}%")
    print(f"  Std deviation: {best_balanced['std_coverage']:.1f}%")
    
    # Compare with current default (18)
    current = coverage_df[coverage_df['evr_index'] == 18].iloc[0]
    print(f"\nCurrent default (EVR 18):")
    print(f"  Mean coverage: {current['mean_coverage']:.1f}%")
    print(f"  Min coverage:  {current['min_coverage']:.1f}%")
    print(f"  Std deviation: {current['std_coverage']:.1f}%")
    
    # EVR 10-12 range analysis
    print("\n" + "=" * 80)
    print("5. DETAILED ANALYSIS OF EVR 10-12 RANGE")
    print("=" * 80)
    
    for evr in [10, 11, 12]:
        row = coverage_df[coverage_df['evr_index'] == evr].iloc[0]
        print(f"\nEVR {evr}:")
        print(f"  Mean coverage: {row['mean_coverage']:.1f}%")
        print(f"  Min coverage:  {row['min_coverage']:.1f}%")
        print(f"  Max coverage:  {row['max_coverage']:.1f}%")
        print(f"  Std deviation: {row['std_coverage']:.1f}%")
        print(f"  CV (coefficient of variation): {row['cv']:.2f}")
    
    return coverage_df, count_df

if __name__ == "__main__":
    coverage_df, count_df = find_optimal_evr()
    
    # Save results
    output_dir = os.path.join(script_dir, 'analysis_results')
    os.makedirs(output_dir, exist_ok=True)
    
    coverage_df.to_csv(os.path.join(output_dir, 'evr_coverage_analysis.csv'), index=False)
    print(f"\nResults saved to {output_dir}/evr_coverage_analysis.csv")
