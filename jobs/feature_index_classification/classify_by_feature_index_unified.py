#!/usr/bin/env python
# coding: utf-8

"""
Unified Feature Index Classification Analysis

This script compares RF classification performance for SPS vs Random sampling
across different datasets (MCC, MCC_01, MCC_05) and feature indices.

Usage:
    python classify_by_feature_index_unified.py --dataset mcc_01 --feature-index 18
    python classify_by_feature_index_unified.py --dataset mcc_05 --feature-index 18
"""

import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
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

PROJECT_ROOT = "/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
DATA_PATH = os.path.join(PROJECT_ROOT, "data")

# Dataset configurations
DATASET_CONFIGS = {
    'mcc': {
        'data_path': os.path.join(DATA_PATH, "mcc/benchmark/30"),
        'feature_index_path': os.path.join(DATA_PATH, "test_feature_index/mcc/30"),
        'output_dir': os.path.join(PROJECT_ROOT, "jobs/feature_index_classification/results"),
        'rare_cell_type': 'osteoblast',
    },
    'mcc_01': {
        'data_path': os.path.join(DATA_PATH, "mcc_01/benchmark/30"),
        'feature_index_path': os.path.join(DATA_PATH, "test_feature_index/mcc_01/30"),
        'output_dir': os.path.join(PROJECT_ROOT, "jobs/feature_index_classification/results_mcc_01"),
        'rare_cell_type': 'osteoblast',
    },
    'mcc_05': {
        'data_path': os.path.join(DATA_PATH, "mcc_05/benchmark/30"),
        'feature_index_path': os.path.join(DATA_PATH, "test_feature_index/mcc_05/30"),
        'output_dir': os.path.join(PROJECT_ROOT, "jobs/feature_index_classification/results_mcc_05"),
        'rare_cell_type': 'osteoblast',
    },
}

# Classification parameters
LABEL_KEY = 'celltype'
TEST_SIZE = 0.2
N_SEEDS = 5
N_TOP_GENES = 2000


# ============================================================================
# Data Loading Functions
# ============================================================================

def load_dataset(config):
    """Load dataset from h5ad file."""
    print(f"Loading dataset from {config['data_path']}...")
    adata_path = os.path.join(config['data_path'], "adata.h5ad")
    adata = sc.read_h5ad(adata_path)
    adata.obs[LABEL_KEY] = adata.obs[LABEL_KEY].astype('category')
    
    # Convert sparse to dense if needed
    if scipy.sparse.issparse(adata.X):
        print("  Converting sparse matrix to dense...")
        adata.X = adata.X.toarray()
    
    print(f"  Dataset shape: {adata.shape}")
    print(f"  Cell types: {adata.obs[LABEL_KEY].value_counts().to_dict()}")
    return adata


def load_sps_indices(config, feature_index, size=100000, rep=0):
    """Load SPS sampling indices for a given feature index."""
    pkl_path = os.path.join(config['feature_index_path'], f"{feature_index}/{size}/{rep}/results.pkl")
    
    if not os.path.exists(pkl_path):
        raise FileNotFoundError(f"SPS indices not found: {pkl_path}")
    
    with open(pkl_path, 'rb') as f:
        res = pickle.load(f)
    
    # Structure is ((indices, evr), time)
    indices = res[0][0]
    evr = res[0][1]
    elapsed_time = res[1]
    
    print(f"  Loaded SPS indices for feature_index={feature_index}: {len(indices)} cells")
    print(f"  EVR: {evr:.4f}, Time: {elapsed_time:.2f}s")
    
    return np.array(indices), evr, elapsed_time


def load_random_indices(config, size=100000, rep=0):
    """Load random sampling indices."""
    pkl_path = os.path.join(config['data_path'], f"random/{size}/{rep}/results.pkl")
    
    if not os.path.exists(pkl_path):
        raise FileNotFoundError(f"Random indices not found: {pkl_path}")
    
    with open(pkl_path, 'rb') as f:
        res = pickle.load(f)
    
    # Structure is (indices, time)
    indices = res[0] if isinstance(res, tuple) else res
    
    print(f"  Loaded Random indices: {len(indices)} cells")
    
    return np.array(indices)


# ============================================================================
# Feature Preparation
# ============================================================================

def prepare_features(adata, n_top_genes=2000):
    """Prepare feature matrix using highly variable genes."""
    print(f"Preparing features using top {n_top_genes} variable genes...")
    
    # Work with a copy
    adata_temp = adata.copy()
    
    # Normalize if not already done
    if 'log1p' not in adata_temp.uns:
        print("  Normalizing data...")
        sc.pp.normalize_total(adata_temp, target_sum=1e4)
        sc.pp.log1p(adata_temp)
    
    # Find highly variable genes
    print("  Finding highly variable genes...")
    sc.pp.highly_variable_genes(
        adata_temp, 
        n_top_genes=min(n_top_genes, adata_temp.n_vars),
        subset=False
    )
    
    # Get HVG indices
    hvg_mask = adata_temp.var['highly_variable'].values
    gene_names = adata_temp.var_names[hvg_mask].tolist()
    
    # Extract feature matrix
    X = adata_temp[:, hvg_mask].X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    
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

def run_analysis(dataset, feature_index, size=100000, rep=0):
    """Run the complete classification comparison for a given dataset and feature index."""
    
    if dataset not in DATASET_CONFIGS:
        raise ValueError(f"Unknown dataset: {dataset}. Available: {list(DATASET_CONFIGS.keys())}")
    
    config = DATASET_CONFIGS[dataset]
    rare_cell_type = config['rare_cell_type']
    
    print("=" * 80)
    print(f"Feature Index Classification Analysis")
    print(f"Dataset: {dataset}, Feature Index: {feature_index}, Size: {size}")
    print("=" * 80)
    
    # Create output directory
    os.makedirs(config['output_dir'], exist_ok=True)
    
    # Check if results already exist
    output_csv = os.path.join(config['output_dir'], f"summary_feature_index_{feature_index}.csv")
    if os.path.exists(output_csv):
        print(f"Results already exist: {output_csv}")
        return None
    
    # Load data
    print("\n" + "-" * 40)
    print("Loading data...")
    print("-" * 40)
    
    adata = load_dataset(config)
    sps_indices, evr, sps_time = load_sps_indices(config, feature_index, size, rep)
    random_indices = load_random_indices(config, size, rep)
    
    # Prepare features
    print("\n" + "-" * 40)
    print("Preparing features...")
    print("-" * 40)
    
    X, gene_names = prepare_features(adata, N_TOP_GENES)
    
    # Encode labels
    label_encoder = LabelEncoder()
    labels = label_encoder.fit_transform(adata.obs[LABEL_KEY].values)
    print(f"Classes: {label_encoder.classes_.tolist()}")
    
    # Cell type distribution in samples
    print("\n" + "-" * 40)
    print("Cell type distribution in samples:")
    print("-" * 40)
    
    for name, indices in [("SPS", sps_indices), ("Random", random_indices)]:
        sample_labels = adata.obs[LABEL_KEY].iloc[indices]
        print(f"\n{name}:")
        for ct in label_encoder.classes_:
            cnt = (sample_labels == ct).sum()
            pct = 100 * cnt / len(indices)
            print(f"  {ct}: {cnt} ({pct:.2f}%)")
    
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
        print(f"  {rare_cell_type} F1: {sps_result['per_class'][rare_cell_type]['f1']:.4f}")
        
        # Random classification
        print("\nRandom Training...")
        random_result = run_classification(X, labels, random_train, test_indices, label_encoder)
        random_result['seed'] = seed
        all_results['random'].append(random_result)
        
        print(f"  Accuracy: {random_result['accuracy']:.4f}")
        print(f"  Macro F1: {random_result['macro_f1']:.4f}")
        print(f"  {rare_cell_type} F1: {random_result['per_class'][rare_cell_type]['f1']:.4f}")
    
    # Aggregate results
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    summary = {
        'dataset': dataset,
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
        summary[f'{method}_{rare_cell_type}_f1_mean'] = np.mean([r['per_class'][rare_cell_type]['f1'] for r in results])
        summary[f'{method}_{rare_cell_type}_f1_std'] = np.std([r['per_class'][rare_cell_type]['f1'] for r in results])
        
        # Per-class F1 means
        for ct in label_encoder.classes_:
            summary[f'{method}_{ct}_f1_mean'] = np.mean([r['per_class'][ct]['f1'] for r in results])
    
    # Improvements
    summary['accuracy_improvement'] = (summary['sps_accuracy_mean'] - summary['random_accuracy_mean']) * 100
    summary['macro_f1_improvement'] = (summary['sps_macro_f1_mean'] - summary['random_macro_f1_mean']) * 100
    summary[f'{rare_cell_type}_f1_improvement'] = (summary[f'sps_{rare_cell_type}_f1_mean'] - summary[f'random_{rare_cell_type}_f1_mean']) * 100
    
    print(f"\nDataset: {dataset}")
    print(f"Feature Index: {feature_index}")
    print(f"EVR: {evr:.4f}")
    print(f"\nAccuracy:")
    print(f"  SPS:    {summary['sps_accuracy_mean']:.4f} ± {summary['sps_accuracy_std']:.4f}")
    print(f"  Random: {summary['random_accuracy_mean']:.4f} ± {summary['random_accuracy_std']:.4f}")
    print(f"  Improvement: {summary['accuracy_improvement']:+.2f}%")
    print(f"\nMacro F1:")
    print(f"  SPS:    {summary['sps_macro_f1_mean']:.4f} ± {summary['sps_macro_f1_std']:.4f}")
    print(f"  Random: {summary['random_macro_f1_mean']:.4f} ± {summary['random_macro_f1_std']:.4f}")
    print(f"  Improvement: {summary['macro_f1_improvement']:+.2f}%")
    print(f"\n{rare_cell_type} F1:")
    print(f"  SPS:    {summary[f'sps_{rare_cell_type}_f1_mean']:.4f} ± {summary[f'sps_{rare_cell_type}_f1_std']:.4f}")
    print(f"  Random: {summary[f'random_{rare_cell_type}_f1_mean']:.4f} ± {summary[f'random_{rare_cell_type}_f1_std']:.4f}")
    print(f"  Improvement: {summary[f'{rare_cell_type}_f1_improvement']:+.2f}%")
    
    # Save results
    # Save detailed results
    output_pkl = os.path.join(config['output_dir'], f"classification_feature_index_{feature_index}.pkl")
    with open(output_pkl, 'wb') as f:
        pickle.dump({
            'summary': summary,
            'all_results': all_results,
            'gene_names': gene_names,
            'label_encoder_classes': label_encoder.classes_.tolist(),
        }, f)
    print(f"\nResults saved to: {output_pkl}")
    
    # Save summary CSV
    pd.DataFrame([summary]).to_csv(output_csv, index=False)
    print(f"Summary saved to: {output_csv}")
    
    return summary


# ============================================================================
# Entry Point
# ============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Unified Feature Index Classification Analysis'
    )
    parser.add_argument(
        '--dataset',
        type=str,
        required=True,
        choices=list(DATASET_CONFIGS.keys()),
        help='Dataset to analyze (mcc, mcc_01, mcc_05)'
    )
    parser.add_argument(
        '--feature-index',
        type=int,
        required=True,
        help='Feature index to analyze (1-30)'
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
            dataset=args.dataset,
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
