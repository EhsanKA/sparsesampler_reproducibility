#!/usr/bin/env python
# coding: utf-8

"""
Combine and analyze results from all feature index classifications.

Usage:
    python combine_results.py
"""

import os
import pandas as pd
import numpy as np
import glob

# Paths
PROJECT_ROOT = "/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
RESULTS_DIR = os.path.join(PROJECT_ROOT, "jobs/feature_index_classification/results")

def main():
    # Find all summary CSV files
    csv_files = sorted(glob.glob(os.path.join(RESULTS_DIR, "summary_feature_index_*.csv")))
    
    if not csv_files:
        print("No results found. Run the classification jobs first.")
        return
    
    print(f"Found {len(csv_files)} result files")
    
    # Combine all results
    dfs = []
    for f in csv_files:
        df = pd.read_csv(f)
        dfs.append(df)
    
    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.sort_values('feature_index')
    
    # Save combined results
    combined_path = os.path.join(RESULTS_DIR, "all_feature_indices_summary.csv")
    combined.to_csv(combined_path, index=False)
    print(f"Combined results saved to: {combined_path}")
    
    # Print summary table
    print("\n" + "=" * 120)
    print("FEATURE INDEX CLASSIFICATION RESULTS")
    print("=" * 120)
    
    print(f"\n{'Index':>6} {'EVR':>8} {'SPS Acc':>10} {'Rand Acc':>10} {'Δ Acc':>8} "
          f"{'SPS F1':>10} {'Rand F1':>10} {'Δ F1':>8} "
          f"{'SPS Osteo':>10} {'Rand Osteo':>10} {'Δ Osteo':>8}")
    print("-" * 120)
    
    for _, row in combined.iterrows():
        print(f"{row['feature_index']:>6} {row['evr']:>8.4f} "
              f"{row['sps_accuracy_mean']:>10.4f} {row['random_accuracy_mean']:>10.4f} "
              f"{row['accuracy_improvement']:>+8.2f}% "
              f"{row['sps_macro_f1_mean']:>10.4f} {row['random_macro_f1_mean']:>10.4f} "
              f"{row['macro_f1_improvement']:>+8.2f}% "
              f"{row['sps_osteoblast_f1_mean']:>10.4f} {row['random_osteoblast_f1_mean']:>10.4f} "
              f"{row['osteoblast_f1_improvement']:>+8.2f}%")
    
    # Find best indices
    print("\n" + "=" * 120)
    print("BEST FEATURE INDICES")
    print("=" * 120)
    
    # Best by macro F1 (SPS)
    best_f1 = combined.loc[combined['sps_macro_f1_mean'].idxmax()]
    print(f"\nBest SPS Macro F1: Index {int(best_f1['feature_index'])} "
          f"(F1={best_f1['sps_macro_f1_mean']:.4f})")
    
    # Best by osteoblast F1 (SPS)
    best_osteo = combined.loc[combined['sps_osteoblast_f1_mean'].idxmax()]
    print(f"Best SPS Osteoblast F1: Index {int(best_osteo['feature_index'])} "
          f"(F1={best_osteo['sps_osteoblast_f1_mean']:.4f})")
    
    # Best by accuracy (SPS)
    best_acc = combined.loc[combined['sps_accuracy_mean'].idxmax()]
    print(f"Best SPS Accuracy: Index {int(best_acc['feature_index'])} "
          f"(Acc={best_acc['sps_accuracy_mean']:.4f})")
    
    # Best improvement in osteoblast F1
    best_osteo_imp = combined.loc[combined['osteoblast_f1_improvement'].idxmax()]
    print(f"Best Osteoblast F1 Improvement: Index {int(best_osteo_imp['feature_index'])} "
          f"(Δ={best_osteo_imp['osteoblast_f1_improvement']:+.2f}%)")
    
    # Best overall (smallest gap to random while maintaining good osteoblast)
    # Score = macro_f1 + 0.5 * osteoblast_improvement (favoring osteoblast)
    combined['composite_score'] = combined['sps_macro_f1_mean'] + 0.005 * combined['osteoblast_f1_improvement']
    best_composite = combined.loc[combined['composite_score'].idxmax()]
    print(f"\nBest Composite (F1 + Osteoblast improvement): Index {int(best_composite['feature_index'])} "
          f"(Macro F1={best_composite['sps_macro_f1_mean']:.4f}, "
          f"Osteoblast Δ={best_composite['osteoblast_f1_improvement']:+.2f}%)")
    
    # Top 5 by composite score
    print("\n" + "-" * 80)
    print("TOP 5 RECOMMENDED FEATURE INDICES (by composite score):")
    print("-" * 80)
    top5 = combined.nlargest(5, 'composite_score')
    for _, row in top5.iterrows():
        print(f"  Index {int(row['feature_index']):>2}: "
              f"Macro F1={row['sps_macro_f1_mean']:.4f}, "
              f"Osteoblast F1={row['sps_osteoblast_f1_mean']:.4f} "
              f"(Δ={row['osteoblast_f1_improvement']:+.2f}%), "
              f"Accuracy={row['sps_accuracy_mean']:.4f}")
    
    # Worst 5
    print("\n" + "-" * 80)
    print("WORST 5 FEATURE INDICES (avoid these):")
    print("-" * 80)
    worst5 = combined.nsmallest(5, 'composite_score')
    for _, row in worst5.iterrows():
        print(f"  Index {int(row['feature_index']):>2}: "
              f"Macro F1={row['sps_macro_f1_mean']:.4f}, "
              f"Osteoblast F1={row['sps_osteoblast_f1_mean']:.4f} "
              f"(Δ={row['osteoblast_f1_improvement']:+.2f}%), "
              f"Accuracy={row['sps_accuracy_mean']:.4f}")
    
    print("\n" + "=" * 120)
    
    return combined


if __name__ == '__main__':
    main()
