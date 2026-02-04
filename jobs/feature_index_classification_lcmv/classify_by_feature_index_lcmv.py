#!/usr/bin/env python
# coding: utf-8

"""
Feature Index Classification Analysis (LCMV Dataset)

This script compares RF classification performance when training on SPS samples
generated with different feature indices vs Random sampling for LCMV dataset.

Goal: Find which feature_index gives:
1. Good overall F1 score
2. Better or similar performance on rare cell types

Usage:
    python classify_by_feature_index_lcmv.py --feature-index 18 --size 100000
"""

import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
import scipy.sparse
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, f1_score, precision_recall_fscore_support,
    confusion_matrix
)
from sklearn.preprocessing import LabelEncoder, StandardScaler
import warnings
import time

warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================

# Paths
PROJECT_ROOT = "/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
LCMV_DATA_PATH = "/fast/AG_Ohler/ekarimi/projects/sparseFlow_benchmarking/data/lcmv/benchmark/34"
TEST_FEATURE_INDEX_PATH = os.path.join(PROJECT_ROOT, "data/test_feature_index/lcmv/34")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "jobs/feature_index_classification_lcmv/results")

# Classification parameters
LABEL_KEY = 'celltype'
TEST_SIZE = 0.2
N_SEEDS = 5
N_TOP_GENES = 31  # LCMV has only 31 genes
REF = 34  # Reference size for LCMV feature index testing (all 34 genes)

# Rare cell types in LCMV (identified from prior analysis)
# These are cell types with <1% frequency
RARE_CELL_TYPES = [
    'interacting',      # ~0.07%
    'NK1_1_TCRgd_T',    # ~0.07%  
    'cDC2',             # ~0.39%
    'pDCs',             # ~0.56%
    'CD4_LCMV_spec',    # ~0.97%
]


# ============================================================================
# Data Loading Functions
# ============================================================================

def load_lcmv_dataset():
    """Load LCMV dataset from h5ad file."""
    print("Loading LCMV dataset...")
    adata_path = os.path.join(LCMV_DATA_PATH, "adata.h5ad")
    adata = sc.read_h5ad(adata_path)
    adata.obs[LABEL_KEY] = adata.obs[LABEL_KEY].astype('category')
    
    # Convert sparse to dense if needed
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    
    print(f"  Dataset shape: {adata.shape}")
    print(f"  Cell types: {adata.obs[LABEL_KEY].value_counts().to_dict()}")
    return adata


def load_sps_indices(feature_index, size=100000, rep=0):
    """Load SPS sampling indices for a given feature index."""
    pkl_path = os.path.join(TEST_FEATURE_INDEX_PATH, f"{feature_index}/{size}/{rep}/results.pkl")
    
    with open(pkl_path, 'rb') as f:
        res = pickle.load(f)
    
    # Structure is ((indices, evr), time)
    indices = res[0][0]
    evr = res[0][1]
    elapsed_time = res[1]
    
    print(f"  Loaded SPS indices for feature_index={feature_index}: {len(indices)} cells")
    print(f"  EVR: {evr:.4f}, Time: {elapsed_time:.2f}s")
    
    return np.array(indices), evr, elapsed_time


def load_random_indices(size=100000, rep=0):
    """Load random sampling indices."""
    pkl_path = os.path.join(LCMV_DATA_PATH, f"random/{size}/{rep}/results.pkl")
    
    with open(pkl_path, 'rb') as f:
        res = pickle.load(f)
    
    # Structure is (indices, time) or just indices
    if isinstance(res, tuple):
        indices = res[0]
    else:
        indices = res
    
    print(f"  Loaded Random indices: {len(indices)} cells")
    
    return np.array(indices)


# ============================================================================
# Feature Preparation
# ============================================================================

def prepare_features(adata):
    """Prepare feature matrix using all genes (LCMV has only 31 genes)."""
    print(f"Preparing features using all {adata.n_vars} genes...")
    
    # Work with a copy
    adata_temp = adata.copy()
    
    # Normalize if not already done
    if 'log1p' not in adata_temp.uns:
        print("  Normalizing data...")
        sc.pp.normalize_total(adata_temp, target_sum=1e4)
        sc.pp.log1p(adata_temp)
    
    # Use all genes for LCMV (only 31 genes)
    X = adata_temp.X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    
    gene_names = adata_temp.var_names.tolist()
    
    print(f"  Feature matrix shape: {X.shape}")
    
    return X, gene_names


# ============================================================================
# Classification Functions
# ============================================================================

def run_classification(X, labels, train_indices, test_indices, label_encoder):
    """Train RF classifier and evaluate."""
    
    # Get training and test data
    X_train = X[train_indices]
    y_train = labels[train_indices]
    X_test = X[test_indices]
    y_test = labels[test_indices]
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train classifier
    clf = RandomForestClassifier(
        n_estimators=100,
        max_depth=20,
        class_weight='balanced',
        n_jobs=-1,
        random_state=42
    )
    
    start_time = time.time()
    clf.fit(X_train_scaled, y_train)
    train_time = time.time() - start_time
    
    # Predictions
    y_pred = clf.predict(X_test_scaled)
    
    # Metrics
    accuracy = accuracy_score(y_test, y_pred)
    macro_f1 = f1_score(y_test, y_pred, average='macro')
    weighted_f1 = f1_score(y_test, y_pred, average='weighted')
    
    # Per-class metrics
    precision, recall, f1, support = precision_recall_fscore_support(
        y_test, y_pred, average=None
    )
    
    class_names = label_encoder.classes_
    per_class = {}
    for i, class_name in enumerate(class_names):
        per_class[class_name] = {
            'precision': precision[i],
            'recall': recall[i],
            'f1': f1[i],
            'support': int(support[i])
        }
    
    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    
    return {
        'accuracy': accuracy,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'per_class': per_class,
        'confusion_matrix': cm,
        'class_names': class_names.tolist(),
        'train_time': train_time,
        'n_train': len(train_indices)
    }


def get_sample_indices_in_pool(sample_indices, pool_indices):
    """Get indices that are in both sample and pool."""
    pool_set = set(pool_indices)
    return np.array([i for i in sample_indices if i in pool_set])


# ============================================================================
# Main Analysis
# ============================================================================

def run_analysis(feature_index, size=100000, rep=0):
    """Run the complete classification comparison for a given feature index."""
    
    print("=" * 80)
    print(f"Feature Index Classification Analysis (LCMV)")
    print(f"Feature Index: {feature_index}, Size: {size}, Rep: {rep}")
    print("=" * 80)
    
    # Load data
    print("\n" + "-" * 40)
    print("Loading data...")
    print("-" * 40)
    
    adata = load_lcmv_dataset()
    sps_indices, evr, sps_time = load_sps_indices(feature_index, size, rep)
    random_indices = load_random_indices(size, rep)
    
    # Prepare features
    print("\n" + "-" * 40)
    print("Preparing features...")
    print("-" * 40)
    
    X, gene_names = prepare_features(adata)
    
    # Encode labels
    label_encoder = LabelEncoder()
    labels = label_encoder.fit_transform(adata.obs[LABEL_KEY].values)
    print(f"Classes: {label_encoder.classes_.tolist()}")
    
    # Identify which rare cell types exist in this dataset
    available_rare_types = [ct for ct in RARE_CELL_TYPES if ct in label_encoder.classes_]
    print(f"Rare cell types to analyze: {available_rare_types}")
    
    # Cell type distribution in samples
    print("\n" + "-" * 40)
    print("Cell type distribution in samples:")
    print("-" * 40)
    
    for name, indices in [("SPS", sps_indices), ("Random", random_indices)]:
        sample_labels = adata.obs[LABEL_KEY].iloc[indices]
        print(f"\n{name}:")
        for ct in available_rare_types:
            cnt = (sample_labels == ct).sum()
            pct = 100 * cnt / len(indices)
            print(f"  {ct}: {cnt} ({pct:.3f}%)")
    
    # Run classification for multiple seeds
    all_results = {'sps': [], 'random': []}
    
    for seed in range(N_SEEDS):
        print(f"\n" + "-" * 40)
        print(f"Seed {seed}")
        print("-" * 40)
        
        # Create train/test split
        all_indices = np.arange(len(labels))
        train_pool, test_indices = train_test_split(
            all_indices,
            test_size=TEST_SIZE,
            stratify=labels,
            random_state=seed
        )
        
        print(f"Train pool: {len(train_pool)}, Test set: {len(test_indices)}")
        
        # Get sample indices in train pool
        sps_train = get_sample_indices_in_pool(sps_indices, train_pool)
        random_train = get_sample_indices_in_pool(random_indices, train_pool)
        
        print(f"SPS in train pool: {len(sps_train)}")
        print(f"Random in train pool: {len(random_train)}")
        
        # SPS classification
        print("\nSPS Training...")
        sps_result = run_classification(X, labels, sps_train, test_indices, label_encoder)
        sps_result['seed'] = seed
        all_results['sps'].append(sps_result)
        
        print(f"  Accuracy: {sps_result['accuracy']:.4f}")
        print(f"  Macro F1: {sps_result['macro_f1']:.4f}")
        for rare_ct in available_rare_types:
            if rare_ct in sps_result['per_class']:
                print(f"  {rare_ct} F1: {sps_result['per_class'][rare_ct]['f1']:.4f}")
        
        # Random classification
        print("\nRandom Training...")
        random_result = run_classification(X, labels, random_train, test_indices, label_encoder)
        random_result['seed'] = seed
        all_results['random'].append(random_result)
        
        print(f"  Accuracy: {random_result['accuracy']:.4f}")
        print(f"  Macro F1: {random_result['macro_f1']:.4f}")
        for rare_ct in available_rare_types:
            if rare_ct in random_result['per_class']:
                print(f"  {rare_ct} F1: {random_result['per_class'][rare_ct]['f1']:.4f}")
    
    # Aggregate results
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    summary = {
        'feature_index': feature_index,
        'size': size,
        'evr': evr,
        'sps_time': sps_time,
    }
    
    for method in ['sps', 'random']:
        results = all_results[method]
        
        summary[f'{method}_accuracy_mean'] = np.mean([r['accuracy'] for r in results])
        summary[f'{method}_accuracy_std'] = np.std([r['accuracy'] for r in results])
        summary[f'{method}_macro_f1_mean'] = np.mean([r['macro_f1'] for r in results])
        summary[f'{method}_macro_f1_std'] = np.std([r['macro_f1'] for r in results])
        
        # Per-class F1 means for rare cell types
        for ct in available_rare_types:
            f1_values = [r['per_class'].get(ct, {}).get('f1', 0) for r in results]
            summary[f'{method}_{ct}_f1_mean'] = np.mean(f1_values)
            summary[f'{method}_{ct}_f1_std'] = np.std(f1_values)
        
        # Per-class F1 means for all classes
        for ct in label_encoder.classes_:
            f1_values = [r['per_class'][ct]['f1'] for r in results]
            safe_ct = ct.replace(' ', '_').replace('/', '_')
            summary[f'{method}_{safe_ct}_f1_mean'] = np.mean(f1_values)
    
    # Improvements
    summary['accuracy_improvement'] = (summary['sps_accuracy_mean'] - summary['random_accuracy_mean']) * 100
    summary['macro_f1_improvement'] = (summary['sps_macro_f1_mean'] - summary['random_macro_f1_mean']) * 100
    
    for rare_ct in available_rare_types:
        sps_key = f'sps_{rare_ct}_f1_mean'
        random_key = f'random_{rare_ct}_f1_mean'
        if sps_key in summary and random_key in summary:
            summary[f'{rare_ct}_f1_improvement'] = (summary[sps_key] - summary[random_key]) * 100
    
    print(f"\nFeature Index: {feature_index}")
    print(f"EVR: {evr:.4f}")
    print(f"\nAccuracy:")
    print(f"  SPS:    {summary['sps_accuracy_mean']:.4f} ± {summary['sps_accuracy_std']:.4f}")
    print(f"  Random: {summary['random_accuracy_mean']:.4f} ± {summary['random_accuracy_std']:.4f}")
    print(f"  Improvement: {summary['accuracy_improvement']:+.2f}%")
    print(f"\nMacro F1:")
    print(f"  SPS:    {summary['sps_macro_f1_mean']:.4f} ± {summary['sps_macro_f1_std']:.4f}")
    print(f"  Random: {summary['random_macro_f1_mean']:.4f} ± {summary['random_macro_f1_std']:.4f}")
    print(f"  Improvement: {summary['macro_f1_improvement']:+.2f}%")
    
    for rare_ct in available_rare_types:
        sps_key = f'sps_{rare_ct}_f1_mean'
        random_key = f'random_{rare_ct}_f1_mean'
        sps_std_key = f'sps_{rare_ct}_f1_std'
        random_std_key = f'random_{rare_ct}_f1_std'
        imp_key = f'{rare_ct}_f1_improvement'
        
        if sps_key in summary:
            print(f"\n{rare_ct} F1:")
            print(f"  SPS:    {summary[sps_key]:.4f} ± {summary.get(sps_std_key, 0):.4f}")
            print(f"  Random: {summary[random_key]:.4f} ± {summary.get(random_std_key, 0):.4f}")
            print(f"  Improvement: {summary.get(imp_key, 0):+.2f}%")
    
    # Save results
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Save detailed results
    output_pkl = os.path.join(OUTPUT_DIR, f"classification_feature_index_{feature_index}.pkl")
    with open(output_pkl, 'wb') as f:
        pickle.dump({
            'summary': summary,
            'all_results': all_results,
            'gene_names': gene_names,
            'label_encoder_classes': label_encoder.classes_.tolist(),
            'rare_cell_types': available_rare_types,
        }, f)
    print(f"\nResults saved to: {output_pkl}")
    
    # Save summary CSV
    summary_csv = os.path.join(OUTPUT_DIR, f"summary_feature_index_{feature_index}.csv")
    pd.DataFrame([summary]).to_csv(summary_csv, index=False)
    print(f"Summary saved to: {summary_csv}")
    
    return summary


# ============================================================================
# Entry Point
# ============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Feature Index Classification Analysis (LCMV)'
    )
    parser.add_argument(
        '--feature-index',
        type=int,
        required=True,
        help='Feature index to analyze (1-25)'
    )
    parser.add_argument(
        '--size',
        type=int,
        default=100000,
        help='Sample size (default: 100000)'
    )
    parser.add_argument(
        '--rep',
        type=int,
        default=0,
        help='Replicate number (default: 0)'
    )
    
    args = parser.parse_args()
    
    try:
        summary = run_analysis(
            feature_index=args.feature_index,
            size=args.size,
            rep=args.rep
        )
        print("\nAnalysis completed successfully!")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
