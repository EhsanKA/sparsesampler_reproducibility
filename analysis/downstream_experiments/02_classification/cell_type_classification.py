#!/usr/bin/env python
# coding: utf-8

"""
Experiment 2: Cell Type Classification Comparison

This script compares classification performance when training on SPS vs Random
sampled data, with evaluation on a common held-out test set.

Hypothesis: Classifier trained on SPS sample generalizes better, especially 
for rare cell types (osteoblast).

Method:
1. Split full dataset: 80% train pool / 20% test (stratified)
2. From train pool, identify cells in SPS 100k OR Random 100k
3. Train classifier (Logistic Regression / Random Forest)
4. Evaluate on same held-out test set

Usage:
    conda activate facs_sampling
    python cell_type_classification.py
    
    # With specific options
    python cell_type_classification.py --n-seeds 5 --classifiers logreg rf
"""

import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, f1_score, classification_report, confusion_matrix,
    roc_auc_score, precision_recall_fscore_support
)
from sklearn.preprocessing import LabelEncoder, StandardScaler
import scipy.sparse
import warnings

warnings.filterwarnings('ignore')

# Add parent directories to path
script_dir = os.path.dirname(os.path.abspath(__file__))
experiments_dir = os.path.dirname(script_dir)
analysis_dir = os.path.dirname(experiments_dir)
project_root = os.path.dirname(analysis_dir)

sys.path.insert(0, experiments_dir)
sys.path.insert(0, analysis_dir)

from common.utils import (
    PROJECT_ROOT,
    DATA_PATH,
    RARE_CELL_TYPE,
    CELL_TYPES,
    LABEL_KEY,
    setup_logger,
    load_mcc_dataset,
    load_sampling_indices,
    get_cell_type_distribution,
)

# ============================================================================
# Configuration
# ============================================================================

OUTPUT_DIR = os.path.join(script_dir, 'results')
FIGURES_DIR = os.path.join(OUTPUT_DIR, 'figures')

# Classification parameters
TEST_SIZE = 0.2  # Fraction of data for test set
N_SEEDS = 5      # Number of random seeds for robustness
N_TOP_GENES = 2000  # Number of highly variable genes to use as features

# Available classifiers
CLASSIFIERS = {
    'logreg': lambda: LogisticRegression(
        max_iter=1000, 
        class_weight='balanced',
        n_jobs=-1,
        random_state=42
    ),
    'rf': lambda: RandomForestClassifier(
        n_estimators=100,
        max_depth=20,
        class_weight='balanced',
        n_jobs=-1,
        random_state=42
    ),
}


# ============================================================================
# Data Preparation Functions
# ============================================================================

def prepare_features(adata, n_top_genes=2000, logger=None):
    """
    Prepare feature matrix from AnnData, using highly variable genes.
    
    Parameters:
    -----------
    adata : AnnData
        Dataset with expression data
    n_top_genes : int
        Number of top variable genes to use
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    tuple : (X, gene_names) where X is the feature matrix
    """
    import scanpy as sc
    
    if logger:
        logger.info(f"Preparing features using top {n_top_genes} variable genes...")
    
    # Work with a copy
    adata_temp = adata.copy()
    
    # Normalize if not already done
    if 'log1p' not in adata_temp.uns:
        if logger:
            logger.info("  Normalizing data...")
        sc.pp.normalize_total(adata_temp, target_sum=1e4)
        sc.pp.log1p(adata_temp)
    
    # Find highly variable genes
    if logger:
        logger.info("  Finding highly variable genes...")
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
    
    if logger:
        logger.info(f"  Feature matrix shape: {X.shape}")
    
    return X, gene_names


def create_train_test_split(adata, test_size=0.2, random_state=42, logger=None):
    """
    Create stratified train/test split.
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    test_size : float
        Fraction for test set
    random_state : int
        Random seed
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    tuple : (train_indices, test_indices)
    """
    if logger:
        logger.info(f"Creating stratified train/test split (test_size={test_size})...")
    
    labels = adata.obs[LABEL_KEY].values
    indices = np.arange(len(labels))
    
    train_idx, test_idx = train_test_split(
        indices,
        test_size=test_size,
        stratify=labels,
        random_state=random_state
    )
    
    if logger:
        logger.info(f"  Train pool size: {len(train_idx)}")
        logger.info(f"  Test set size: {len(test_idx)}")
    
    return train_idx, test_idx


def get_sample_indices_in_pool(sample_indices, pool_indices, logger=None):
    """
    Get the subset of sample indices that are in the training pool.
    
    Parameters:
    -----------
    sample_indices : np.ndarray
        Indices from sampling method (SPS or Random)
    pool_indices : np.ndarray
        Indices in the training pool
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    np.ndarray : Indices that are in both sample and pool
    """
    pool_set = set(pool_indices)
    in_pool = np.array([i for i in sample_indices if i in pool_set])
    
    if logger:
        logger.info(f"  Sample indices in train pool: {len(in_pool)} / {len(sample_indices)}")
    
    return in_pool


# ============================================================================
# Classification Functions
# ============================================================================

def train_classifier(X_train, y_train, classifier_type='logreg', logger=None):
    """
    Train a classifier.
    
    Parameters:
    -----------
    X_train : np.ndarray
        Training features
    y_train : np.ndarray
        Training labels (encoded)
    classifier_type : str
        Type of classifier ('logreg' or 'rf')
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    sklearn classifier : Trained classifier
    """
    if logger:
        logger.info(f"  Training {classifier_type} classifier...")
    
    clf = CLASSIFIERS[classifier_type]()
    clf.fit(X_train, y_train)
    
    return clf


def evaluate_classifier(clf, X_test, y_test, label_encoder, logger=None):
    """
    Evaluate classifier performance.
    
    Parameters:
    -----------
    clf : sklearn classifier
        Trained classifier
    X_test : np.ndarray
        Test features
    y_test : np.ndarray
        Test labels (encoded)
    label_encoder : LabelEncoder
        Label encoder for decoding
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    dict : Evaluation metrics
    """
    # Predictions
    y_pred = clf.predict(X_test)
    y_prob = clf.predict_proba(X_test) if hasattr(clf, 'predict_proba') else None
    
    # Class names
    class_names = label_encoder.classes_
    
    # Overall metrics
    accuracy = accuracy_score(y_test, y_pred)
    macro_f1 = f1_score(y_test, y_pred, average='macro')
    weighted_f1 = f1_score(y_test, y_pred, average='weighted')
    
    # Per-class metrics
    precision, recall, f1, support = precision_recall_fscore_support(
        y_test, y_pred, average=None
    )
    
    per_class_metrics = {}
    for i, class_name in enumerate(class_names):
        per_class_metrics[class_name] = {
            'precision': precision[i],
            'recall': recall[i],
            'f1': f1[i],
            'support': support[i]
        }
    
    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    
    # ROC-AUC (if probabilities available)
    roc_auc = None
    per_class_auc = {}
    if y_prob is not None:
        try:
            roc_auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')
            # Per-class AUC
            for i, class_name in enumerate(class_names):
                y_true_binary = (y_test == i).astype(int)
                per_class_auc[class_name] = roc_auc_score(y_true_binary, y_prob[:, i])
        except Exception as e:
            if logger:
                logger.warning(f"Could not compute ROC-AUC: {e}")
    
    metrics = {
        'accuracy': accuracy,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'roc_auc': roc_auc,
        'per_class': per_class_metrics,
        'per_class_auc': per_class_auc,
        'confusion_matrix': cm,
        'class_names': class_names.tolist(),
        'y_pred': y_pred,
        'y_prob': y_prob,
    }
    
    if logger:
        logger.info(f"  Accuracy: {accuracy:.4f}")
        logger.info(f"  Macro F1: {macro_f1:.4f}")
        logger.info(f"  Weighted F1: {weighted_f1:.4f}")
        if roc_auc:
            logger.info(f"  ROC-AUC: {roc_auc:.4f}")
        logger.info(f"  Osteoblast F1: {per_class_metrics.get('osteoblast', {}).get('f1', 'N/A'):.4f}")
    
    return metrics


def run_single_comparison(adata, X, sps_indices, random_indices, 
                          train_pool, test_indices, 
                          label_encoder, classifier_type='logreg',
                          seed=42, logger=None):
    """
    Run a single comparison between SPS and Random training.
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    X : np.ndarray
        Feature matrix for all cells
    sps_indices : np.ndarray
        SPS sampling indices
    random_indices : np.ndarray
        Random sampling indices  
    train_pool : np.ndarray
        Training pool indices
    test_indices : np.ndarray
        Test set indices
    label_encoder : LabelEncoder
        Label encoder
    classifier_type : str
        Type of classifier
    seed : int
        Random seed
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    dict : Comparison results
    """
    if logger:
        logger.info(f"\n--- Seed {seed} ---")
    
    # Get labels
    labels = label_encoder.transform(adata.obs[LABEL_KEY].values)
    
    # Get sample indices that are in the training pool
    sps_train = get_sample_indices_in_pool(sps_indices, train_pool, logger)
    random_train = get_sample_indices_in_pool(random_indices, train_pool, logger)
    
    if logger:
        # Log cell type distribution
        sps_labels = labels[sps_train]
        random_labels = labels[random_train]
        
        for method, train_labels in [('SPS', sps_labels), ('Random', random_labels)]:
            unique, counts = np.unique(train_labels, return_counts=True)
            logger.info(f"  {method} training set distribution:")
            for u, c in zip(unique, counts):
                class_name = label_encoder.inverse_transform([u])[0]
                logger.info(f"    {class_name}: {c} ({100*c/len(train_labels):.2f}%)")
    
    # Prepare data
    X_test = X[test_indices]
    y_test = labels[test_indices]
    
    results = {}
    
    # Train and evaluate on SPS
    if logger:
        logger.info("\nSPS Training:")
    X_train_sps = X[sps_train]
    y_train_sps = labels[sps_train]
    
    # Scale features
    scaler_sps = StandardScaler()
    X_train_sps_scaled = scaler_sps.fit_transform(X_train_sps)
    X_test_sps_scaled = scaler_sps.transform(X_test)
    
    clf_sps = train_classifier(X_train_sps_scaled, y_train_sps, classifier_type, logger)
    results['sps'] = evaluate_classifier(clf_sps, X_test_sps_scaled, y_test, label_encoder, logger)
    results['sps']['n_train'] = len(sps_train)
    results['sps']['train_indices'] = sps_train
    
    # Train and evaluate on Random
    if logger:
        logger.info("\nRandom Training:")
    X_train_random = X[random_train]
    y_train_random = labels[random_train]
    
    scaler_random = StandardScaler()
    X_train_random_scaled = scaler_random.fit_transform(X_train_random)
    X_test_random_scaled = scaler_random.transform(X_test)
    
    clf_random = train_classifier(X_train_random_scaled, y_train_random, classifier_type, logger)
    results['random'] = evaluate_classifier(clf_random, X_test_random_scaled, y_test, label_encoder, logger)
    results['random']['n_train'] = len(random_train)
    results['random']['train_indices'] = random_train
    
    # Comparison summary
    if logger:
        logger.info("\nComparison Summary:")
        logger.info(f"  Accuracy: SPS={results['sps']['accuracy']:.4f}, "
                   f"Random={results['random']['accuracy']:.4f}, "
                   f"Diff={(results['sps']['accuracy']-results['random']['accuracy'])*100:+.2f}%")
        logger.info(f"  Macro F1: SPS={results['sps']['macro_f1']:.4f}, "
                   f"Random={results['random']['macro_f1']:.4f}, "
                   f"Diff={(results['sps']['macro_f1']-results['random']['macro_f1'])*100:+.2f}%")
        
        # Osteoblast F1
        sps_osteo_f1 = results['sps']['per_class'].get('osteoblast', {}).get('f1', 0)
        random_osteo_f1 = results['random']['per_class'].get('osteoblast', {}).get('f1', 0)
        logger.info(f"  Osteoblast F1: SPS={sps_osteo_f1:.4f}, "
                   f"Random={random_osteo_f1:.4f}, "
                   f"Diff={(sps_osteo_f1-random_osteo_f1)*100:+.2f}%")
    
    return results


# ============================================================================
# Visualization Functions
# ============================================================================

def plot_f1_comparison(all_results, classifier_type, output_dir):
    """Create F1 score comparison plot across classes."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Get class names and aggregate results
    class_names = all_results[0]['sps']['class_names']
    
    sps_f1s = {c: [] for c in class_names}
    random_f1s = {c: [] for c in class_names}
    
    for result in all_results:
        for c in class_names:
            sps_f1s[c].append(result['sps']['per_class'][c]['f1'])
            random_f1s[c].append(result['random']['per_class'][c]['f1'])
    
    # 1. Bar chart with error bars
    ax = axes[0]
    x = np.arange(len(class_names))
    width = 0.35
    
    sps_means = [np.mean(sps_f1s[c]) for c in class_names]
    sps_stds = [np.std(sps_f1s[c]) for c in class_names]
    random_means = [np.mean(random_f1s[c]) for c in class_names]
    random_stds = [np.std(random_f1s[c]) for c in class_names]
    
    bars1 = ax.bar(x - width/2, sps_means, width, yerr=sps_stds, 
                   label='SPS', color='#2ecc71', alpha=0.7, capsize=3)
    bars2 = ax.bar(x + width/2, random_means, width, yerr=random_stds,
                   label='Random', color='#e74c3c', alpha=0.7, capsize=3)
    
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('F1 Score')
    ax.set_title(f'Per-Class F1 Score ({classifier_type})\nMean ± Std over {len(all_results)} seeds')
    ax.set_xticks(x)
    ax.set_xticklabels(class_names, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.05)
    
    # 2. Box plot for osteoblast (rare cell)
    ax = axes[1]
    data = [sps_f1s['osteoblast'], random_f1s['osteoblast']]
    bp = ax.boxplot(data, labels=['SPS', 'Random'], patch_artist=True)
    bp['boxes'][0].set_facecolor('#2ecc71')
    bp['boxes'][1].set_facecolor('#e74c3c')
    for patch in bp['boxes']:
        patch.set_alpha(0.7)
    
    ax.set_ylabel('F1 Score')
    ax.set_title(f'Osteoblast (Rare Cell) F1 Score\n({len(all_results)} seeds)')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'f1_comparison_{classifier_type}.png'), 
                dpi=150, bbox_inches='tight')
    plt.close()


def plot_confusion_matrices(result, classifier_type, output_dir, seed=0):
    """Create side-by-side confusion matrices."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    class_names = result['sps']['class_names']
    
    for ax, (method, title, cmap) in zip(axes, 
                                          [('sps', 'SPS', 'Greens'),
                                           ('random', 'Random', 'Reds')]):
        cm = result[method]['confusion_matrix']
        
        # Normalize
        cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        
        sns.heatmap(cm_norm, annot=True, fmt='.2f', cmap=cmap,
                   xticklabels=class_names, yticklabels=class_names,
                   ax=ax, vmin=0, vmax=1)
        ax.set_xlabel('Predicted')
        ax.set_ylabel('True')
        ax.set_title(f'{title} Training\n(Acc={result[method]["accuracy"]:.3f}, '
                    f'F1={result[method]["macro_f1"]:.3f})')
    
    plt.suptitle(f'Confusion Matrices ({classifier_type}, seed={seed})', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'confusion_matrix_{classifier_type}_seed{seed}.png'),
                dpi=150, bbox_inches='tight')
    plt.close()


def plot_overall_comparison(all_results, classifier_type, output_dir):
    """Create overall metrics comparison plot."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    metrics = ['accuracy', 'macro_f1', 'weighted_f1']
    titles = ['Accuracy', 'Macro F1', 'Weighted F1']
    
    for ax, metric, title in zip(axes, metrics, titles):
        sps_vals = [r['sps'][metric] for r in all_results]
        random_vals = [r['random'][metric] for r in all_results]
        
        positions = [1, 2]
        bp = ax.boxplot([sps_vals, random_vals], positions=positions, 
                       labels=['SPS', 'Random'], patch_artist=True)
        bp['boxes'][0].set_facecolor('#2ecc71')
        bp['boxes'][1].set_facecolor('#e74c3c')
        for patch in bp['boxes']:
            patch.set_alpha(0.7)
        
        # Add individual points
        for i, vals in enumerate([sps_vals, random_vals]):
            x = np.random.normal(positions[i], 0.04, len(vals))
            ax.scatter(x, vals, alpha=0.5, s=20, color='black')
        
        ax.set_ylabel(title)
        ax.set_title(title)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Show mean difference
        mean_diff = np.mean(sps_vals) - np.mean(random_vals)
        ax.annotate(f'Δ = {mean_diff:+.4f}', xy=(1.5, min(sps_vals + random_vals)),
                   ha='center', fontsize=10)
    
    plt.suptitle(f'Overall Classification Metrics ({classifier_type})\n{len(all_results)} seeds',
                fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'overall_comparison_{classifier_type}.png'),
                dpi=150, bbox_inches='tight')
    plt.close()


def plot_auc_comparison(all_results, classifier_type, output_dir):
    """Create ROC-AUC comparison plot."""
    if all_results[0]['sps']['per_class_auc'] is None:
        return
    
    class_names = all_results[0]['sps']['class_names']
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(class_names))
    width = 0.35
    
    sps_aucs = {c: [] for c in class_names}
    random_aucs = {c: [] for c in class_names}
    
    for result in all_results:
        for c in class_names:
            if c in result['sps']['per_class_auc']:
                sps_aucs[c].append(result['sps']['per_class_auc'][c])
                random_aucs[c].append(result['random']['per_class_auc'][c])
    
    sps_means = [np.mean(sps_aucs[c]) if sps_aucs[c] else 0 for c in class_names]
    random_means = [np.mean(random_aucs[c]) if random_aucs[c] else 0 for c in class_names]
    
    ax.bar(x - width/2, sps_means, width, label='SPS', color='#2ecc71', alpha=0.7)
    ax.bar(x + width/2, random_means, width, label='Random', color='#e74c3c', alpha=0.7)
    
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('ROC-AUC')
    ax.set_title(f'Per-Class ROC-AUC ({classifier_type})')
    ax.set_xticks(x)
    ax.set_xticklabels(class_names, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0.5, 1.05)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='Random chance')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'auc_comparison_{classifier_type}.png'),
                dpi=150, bbox_inches='tight')
    plt.close()


def create_summary_table(all_results, classifier_types, output_dir):
    """Create summary table of results."""
    rows = []
    
    for clf_type, results in all_results.items():
        # Aggregate metrics
        sps_acc = [r['sps']['accuracy'] for r in results]
        random_acc = [r['random']['accuracy'] for r in results]
        sps_f1 = [r['sps']['macro_f1'] for r in results]
        random_f1 = [r['random']['macro_f1'] for r in results]
        sps_osteo_f1 = [r['sps']['per_class'].get('osteoblast', {}).get('f1', 0) for r in results]
        random_osteo_f1 = [r['random']['per_class'].get('osteoblast', {}).get('f1', 0) for r in results]
        
        row = {
            'classifier': clf_type,
            'n_seeds': len(results),
            'sps_accuracy_mean': np.mean(sps_acc),
            'sps_accuracy_std': np.std(sps_acc),
            'random_accuracy_mean': np.mean(random_acc),
            'random_accuracy_std': np.std(random_acc),
            'accuracy_improvement': (np.mean(sps_acc) - np.mean(random_acc)) * 100,
            'sps_macro_f1_mean': np.mean(sps_f1),
            'sps_macro_f1_std': np.std(sps_f1),
            'random_macro_f1_mean': np.mean(random_f1),
            'random_macro_f1_std': np.std(random_f1),
            'macro_f1_improvement': (np.mean(sps_f1) - np.mean(random_f1)) * 100,
            'sps_osteoblast_f1_mean': np.mean(sps_osteo_f1),
            'sps_osteoblast_f1_std': np.std(sps_osteo_f1),
            'random_osteoblast_f1_mean': np.mean(random_osteo_f1),
            'random_osteoblast_f1_std': np.std(random_osteo_f1),
            'osteoblast_f1_improvement': (np.mean(sps_osteo_f1) - np.mean(random_osteo_f1)) * 100,
        }
        rows.append(row)
    
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(output_dir, 'classification_summary.csv'), index=False)
    
    # Create readable summary
    summary_text = """
Cell Type Classification Summary
================================

"""
    for _, row in df.iterrows():
        summary_text += f"""
{row['classifier'].upper()} Classifier ({row['n_seeds']} seeds):
  Accuracy:
    SPS:    {row['sps_accuracy_mean']:.4f} ± {row['sps_accuracy_std']:.4f}
    Random: {row['random_accuracy_mean']:.4f} ± {row['random_accuracy_std']:.4f}
    Improvement: {row['accuracy_improvement']:+.2f}%
    
  Macro F1:
    SPS:    {row['sps_macro_f1_mean']:.4f} ± {row['sps_macro_f1_std']:.4f}
    Random: {row['random_macro_f1_mean']:.4f} ± {row['random_macro_f1_std']:.4f}
    Improvement: {row['macro_f1_improvement']:+.2f}%
    
  Osteoblast F1 (rare cell):
    SPS:    {row['sps_osteoblast_f1_mean']:.4f} ± {row['sps_osteoblast_f1_std']:.4f}
    Random: {row['random_osteoblast_f1_mean']:.4f} ± {row['random_osteoblast_f1_std']:.4f}
    Improvement: {row['osteoblast_f1_improvement']:+.2f}%
"""
    
    with open(os.path.join(output_dir, 'classification_summary.txt'), 'w') as f:
        f.write(summary_text)
    
    return df, summary_text


# ============================================================================
# Main Analysis
# ============================================================================

def run_analysis(classifier_types=None, n_seeds=5, n_top_genes=2000, logger=None):
    """
    Run the complete classification comparison analysis.
    
    Parameters:
    -----------
    classifier_types : list
        Classifiers to use
    n_seeds : int
        Number of random seeds
    n_top_genes : int
        Number of top variable genes for features
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    dict : All results
    """
    if classifier_types is None:
        classifier_types = ['logreg', 'rf']
    
    if logger is None:
        logger = setup_logger('cell_type_classification')
    
    # Create output directories
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)
    
    logger.info("=" * 80)
    logger.info("Experiment 2: Cell Type Classification Comparison")
    logger.info("=" * 80)
    
    # Load data
    logger.info("\n" + "=" * 40)
    logger.info("Step 1: Loading data")
    logger.info("=" * 40)
    
    adata = load_mcc_dataset(logger)
    sps_indices = load_sampling_indices('sps', logger)
    random_indices = load_sampling_indices('random', logger)
    
    # Prepare features
    logger.info("\n" + "=" * 40)
    logger.info("Step 2: Preparing features")
    logger.info("=" * 40)
    
    X, gene_names = prepare_features(adata, n_top_genes, logger)
    
    # Encode labels
    label_encoder = LabelEncoder()
    label_encoder.fit(adata.obs[LABEL_KEY].values)
    logger.info(f"Classes: {label_encoder.classes_.tolist()}")
    
    # Run classification for each classifier and seed
    all_results = {clf: [] for clf in classifier_types}
    
    for seed in range(n_seeds):
        logger.info("\n" + "=" * 40)
        logger.info(f"Step 3: Classification (Seed {seed})")
        logger.info("=" * 40)
        
        # Create train/test split
        train_pool, test_indices = create_train_test_split(
            adata, TEST_SIZE, random_state=seed, logger=logger
        )
        
        for clf_type in classifier_types:
            logger.info(f"\n--- Classifier: {clf_type} ---")
            
            result = run_single_comparison(
                adata, X, sps_indices, random_indices,
                train_pool, test_indices, label_encoder,
                classifier_type=clf_type, seed=seed, logger=logger
            )
            result['seed'] = seed
            result['test_indices'] = test_indices
            all_results[clf_type].append(result)
    
    # Create visualizations
    logger.info("\n" + "=" * 40)
    logger.info("Step 4: Creating visualizations")
    logger.info("=" * 40)
    
    for clf_type, results in all_results.items():
        logger.info(f"\nCreating plots for {clf_type}...")
        
        plot_f1_comparison(results, clf_type, FIGURES_DIR)
        plot_overall_comparison(results, clf_type, FIGURES_DIR)
        plot_auc_comparison(results, clf_type, FIGURES_DIR)
        
        # Confusion matrix for first seed
        plot_confusion_matrices(results[0], clf_type, FIGURES_DIR, seed=0)
    
    # Create summary
    logger.info("\n" + "=" * 40)
    logger.info("Step 5: Creating summary")
    logger.info("=" * 40)
    
    summary_df, summary_text = create_summary_table(all_results, classifier_types, OUTPUT_DIR)
    logger.info(summary_text)
    
    # Save all results
    with open(os.path.join(OUTPUT_DIR, 'classification_results.pkl'), 'wb') as f:
        pickle.dump({
            'results': all_results,
            'gene_names': gene_names,
            'label_encoder_classes': label_encoder.classes_.tolist(),
        }, f)
    
    logger.info("\n" + "=" * 80)
    logger.info("Analysis complete!")
    logger.info(f"Results saved to: {OUTPUT_DIR}")
    logger.info(f"Figures saved to: {FIGURES_DIR}")
    logger.info("=" * 80)
    
    return all_results


# ============================================================================
# Entry Point
# ============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Cell Type Classification Comparison'
    )
    parser.add_argument(
        '--classifiers',
        type=str,
        nargs='+',
        default=['logreg', 'rf'],
        choices=list(CLASSIFIERS.keys()),
        help='Classifiers to use'
    )
    parser.add_argument(
        '--n-seeds',
        type=int,
        default=5,
        help='Number of random seeds for robustness (default: 5)'
    )
    parser.add_argument(
        '--n-top-genes',
        type=int,
        default=2000,
        help='Number of highly variable genes for features (default: 2000)'
    )
    
    args = parser.parse_args()
    
    logger = setup_logger('cell_type_classification')
    
    try:
        results = run_analysis(
            classifier_types=args.classifiers,
            n_seeds=args.n_seeds,
            n_top_genes=args.n_top_genes,
            logger=logger
        )
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
