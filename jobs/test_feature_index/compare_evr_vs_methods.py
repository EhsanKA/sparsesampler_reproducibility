#!/usr/bin/env python
"""
Compare EVR index performance against other sampling methods.

Find which EVR indices achieve at least 80% of the best method's performance.
"""

import os
import pandas as pd
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
EVR_TABLES_DIR = os.path.join(script_dir, 'tables')
METHODS_TABLES_DIR = os.path.join(os.path.dirname(script_dir), 'test_sampling_methods', 'tables')

# Dataset configurations
DATASETS = {
    'lcmv': {'refs': ['1M', '5M', '10M', '20M', '34M']},
    'mcc': {'refs': ['5M', '10M', '20M', '25M', '30M']},
    'mcc_01': {'refs': ['5M', '10M', '20M', '25M', '30M']},
    'mcc_05': {'refs': ['5M', '10M', '20M', '25M', '30M']},
}

def load_evr_table(dataset):
    """Load EVR index table for a dataset."""
    path = os.path.join(EVR_TABLES_DIR, f'{dataset}_feature_index_table_size_100000.csv')
    if os.path.exists(path):
        return pd.read_csv(path)
    return None

def load_methods_table(dataset):
    """Load sampling methods table for a dataset."""
    path = os.path.join(METHODS_TABLES_DIR, f'{dataset}_sampling_methods_table_size_100000.csv')
    if os.path.exists(path):
        return pd.read_csv(path)
    return None

def analyze_threshold(threshold=0.80):
    """Analyze which EVR indices pass the threshold compared to best method."""
    
    print("=" * 90)
    print(f"EVR INDEX vs BEST METHOD COMPARISON (Threshold: {threshold*100:.0f}% of best)")
    print("=" * 90)
    
    all_results = []
    
    for dataset, cfg in DATASETS.items():
        evr_df = load_evr_table(dataset)
        methods_df = load_methods_table(dataset)
        
        if evr_df is None or methods_df is None:
            print(f"\nSkipping {dataset}: missing data")
            continue
        
        print(f"\n{'='*90}")
        print(f"DATASET: {dataset.upper()}")
        print("=" * 90)
        
        refs = cfg['refs']
        
        # Find best method for each reference
        print("\n1. Best method per reference:")
        print("-" * 60)
        best_per_ref = {}
        for ref in refs:
            if ref in methods_df.columns:
                methods_df_temp = methods_df.copy()
                # Exclude methods with 0 (likely missing data)
                valid_methods = methods_df_temp[methods_df_temp[ref] > 0]
                if len(valid_methods) > 0:
                    best_idx = valid_methods[ref].idxmax()
                    best_method = valid_methods.loc[best_idx, 'method']
                    best_value = valid_methods.loc[best_idx, ref]
                    best_per_ref[ref] = {'method': best_method, 'value': best_value}
                    print(f"  {ref}: {best_method} = {best_value:.0f}")
        
        # Compare EVR indices against best
        print(f"\n2. EVR indices achieving >= {threshold*100:.0f}% of best:")
        print("-" * 60)
        
        evr_pass_count = {}  # Count how many refs each EVR passes
        
        for evr_idx in range(1, 31):
            if evr_idx not in evr_df['feature_index'].values:
                continue
            
            evr_row = evr_df[evr_df['feature_index'] == evr_idx].iloc[0]
            passes = 0
            total_refs = 0
            
            for ref in refs:
                if ref not in evr_row.index or ref not in best_per_ref:
                    continue
                
                evr_value = evr_row[ref]
                best_value = best_per_ref[ref]['value']
                threshold_value = best_value * threshold
                total_refs += 1
                
                if evr_value >= threshold_value:
                    passes += 1
                    
                all_results.append({
                    'dataset': dataset,
                    'evr_index': evr_idx,
                    'ref': ref,
                    'evr_value': evr_value,
                    'best_value': best_value,
                    'best_method': best_per_ref[ref]['method'],
                    'ratio': evr_value / best_value if best_value > 0 else 0,
                    'passes': evr_value >= threshold_value
                })
            
            if total_refs > 0:
                evr_pass_count[evr_idx] = (passes, total_refs)
        
        # Print EVR indices that pass all refs
        print(f"\n  EVR indices passing ALL references ({len(refs)} refs):")
        for evr_idx, (passes, total) in evr_pass_count.items():
            if passes == total:
                print(f"    EVR {evr_idx}: {passes}/{total} refs")
        
        # Print EVR indices that pass at least 80% of refs
        print(f"\n  EVR indices passing >= 80% of references:")
        for evr_idx, (passes, total) in evr_pass_count.items():
            pass_pct = passes / total * 100
            if pass_pct >= 80:
                print(f"    EVR {evr_idx}: {passes}/{total} refs ({pass_pct:.0f}%)")
    
    # Overall summary
    results_df = pd.DataFrame(all_results)
    
    print("\n" + "=" * 90)
    print("OVERALL SUMMARY: EVR indices passing >= 80% of best across ALL datasets")
    print("=" * 90)
    
    # Count passes per EVR index across all dataset/ref combinations
    evr_summary = results_df.groupby('evr_index').agg({
        'passes': ['sum', 'count'],
        'ratio': 'mean'
    }).reset_index()
    evr_summary.columns = ['evr_index', 'passes', 'total', 'mean_ratio']
    evr_summary['pass_rate'] = evr_summary['passes'] / evr_summary['total'] * 100
    
    print(f"\n{'EVR':<6} {'Passes':>8} {'Total':>8} {'Pass%':>10} {'Mean Ratio':>12}")
    print("-" * 50)
    
    # Sort by pass rate descending
    evr_summary_sorted = evr_summary.sort_values('pass_rate', ascending=False)
    
    for _, row in evr_summary_sorted.iterrows():
        evr = int(row['evr_index'])
        passes = int(row['passes'])
        total = int(row['total'])
        pass_rate = row['pass_rate']
        mean_ratio = row['mean_ratio']
        
        marker = ""
        if pass_rate >= 80:
            marker = " <-- PASSES"
        
        print(f"{evr:<6} {passes:>8} {total:>8} {pass_rate:>9.1f}% {mean_ratio:>11.2f}{marker}")
    
    # List EVR indices that pass
    passing_evrs = evr_summary[evr_summary['pass_rate'] >= 80]['evr_index'].tolist()
    
    print("\n" + "=" * 90)
    print(f"EVR INDICES THAT PASS (>= 80% of best in >= 80% of cases):")
    print("=" * 90)
    print(f"\n  {sorted(passing_evrs)}")
    
    # Recommendation considering stability
    if passing_evrs:
        # Among passing EVRs, find one with best average ratio
        passing_data = evr_summary[evr_summary['evr_index'].isin(passing_evrs)]
        best_evr = passing_data.loc[passing_data['mean_ratio'].idxmax(), 'evr_index']
        print(f"\n  Recommended EVR (highest mean ratio among passing): {int(best_evr)}")
    
    return results_df, evr_summary

if __name__ == "__main__":
    results_df, evr_summary = analyze_threshold(threshold=0.80)
    
    # Save results
    output_dir = os.path.join(script_dir, 'analysis_results')
    os.makedirs(output_dir, exist_ok=True)
    
    results_df.to_csv(os.path.join(output_dir, 'evr_vs_methods_detailed.csv'), index=False)
    evr_summary.to_csv(os.path.join(output_dir, 'evr_vs_methods_summary.csv'), index=False)
    print(f"\nResults saved to {output_dir}/")
