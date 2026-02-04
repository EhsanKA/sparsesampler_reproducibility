#!/usr/bin/env python
# coding: utf-8

"""
Combine all feature index classification results into a single summary CSV.

Usage:
    python combine_results.py
"""

import os
import pandas as pd
import glob

PROJECT_ROOT = "/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
RESULTS_DIR = os.path.join(PROJECT_ROOT, "jobs/feature_index_classification_lcmv/results")

def combine_results():
    """Combine all summary CSVs into one."""
    
    # Find all summary files
    pattern = os.path.join(RESULTS_DIR, "summary_feature_index_*.csv")
    files = glob.glob(pattern)
    
    if not files:
        print("No summary files found!")
        return None
    
    print(f"Found {len(files)} summary files")
    
    # Load and combine
    dfs = []
    for f in sorted(files):
        try:
            df = pd.read_csv(f)
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {f}: {e}")
    
    if not dfs:
        print("No valid data loaded!")
        return None
    
    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.sort_values('feature_index')
    
    # Save combined results
    output_file = os.path.join(RESULTS_DIR, "all_feature_indices_summary.csv")
    combined.to_csv(output_file, index=False)
    print(f"\nCombined results saved to: {output_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    print(f"\nTotal feature indices: {len(combined)}")
    print(f"Feature index range: {combined['feature_index'].min()} - {combined['feature_index'].max()}")
    
    if 'accuracy_improvement' in combined.columns:
        best_acc = combined.loc[combined['accuracy_improvement'].idxmax()]
        print(f"\nBest accuracy improvement: Feature index {int(best_acc['feature_index'])} ({best_acc['accuracy_improvement']:+.2f}%)")
    
    if 'macro_f1_improvement' in combined.columns:
        best_f1 = combined.loc[combined['macro_f1_improvement'].idxmax()]
        print(f"Best macro F1 improvement: Feature index {int(best_f1['feature_index'])} ({best_f1['macro_f1_improvement']:+.2f}%)")
    
    # Check for rare cell type improvements
    rare_ct_columns = [c for c in combined.columns if '_f1_improvement' in c and c not in ['accuracy_improvement', 'macro_f1_improvement']]
    
    for col in rare_ct_columns:
        cell_type = col.replace('_f1_improvement', '')
        best_row = combined.loc[combined[col].idxmax()]
        print(f"Best {cell_type} F1 improvement: Feature index {int(best_row['feature_index'])} ({best_row[col]:+.2f}%)")
    
    return combined


if __name__ == '__main__':
    combine_results()
