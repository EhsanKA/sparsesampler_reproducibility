#!/usr/bin/env python
# coding: utf-8

"""
Experiment 1: Statistical Power Analysis for Rare Cell Markers

This script compares the statistical power of SPS vs Random sampling for
detecting marker genes of rare cell types (particularly osteoblast).

Hypothesis: SPS sampling captures more osteoblast cells, leading to:
- Lower p-values for osteoblast markers (stronger statistical signal)
- Higher log fold changes 
- More markers passing stringent significance thresholds

Usage:
    conda activate facs_sampling
    python statistical_power_analysis.py
    
    # Or with custom parameters
    python statistical_power_analysis.py --de-method wilcoxon --top-n 50
"""

import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
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
    setup_logger,
    load_mcc_dataset,
    load_sampling_indices,
    load_marker_gene_results,
    get_cell_type_distribution,
    compute_osteoblast_statistics,
)

# ============================================================================
# Configuration
# ============================================================================

OUTPUT_DIR = os.path.join(script_dir, 'results')
FIGURES_DIR = os.path.join(OUTPUT_DIR, 'figures')

# DE methods available
DE_METHODS = ['wilcoxon', 't-test', 'logreg']

# Significance thresholds for analysis
P_THRESHOLDS = [1e-10, 1e-50, 1e-100]


# ============================================================================
# Analysis Functions
# ============================================================================

def load_osteoblast_markers(de_method='wilcoxon'):
    """
    Load osteoblast marker gene results for SPS and Random sampling.
    
    Parameters:
    -----------
    de_method : str
        DE method used
        
    Returns:
    --------
    tuple : (sps_markers, random_markers, full_markers) DataFrames
    """
    sps_df = load_marker_gene_results('sps', de_method)
    random_df = load_marker_gene_results('random', de_method)
    full_df = load_marker_gene_results('full', de_method)
    
    # Filter for osteoblast
    sps_osteo = sps_df[sps_df['cell_type'] == RARE_CELL_TYPE].copy()
    random_osteo = random_df[random_df['cell_type'] == RARE_CELL_TYPE].copy()
    full_osteo = full_df[full_df['cell_type'] == RARE_CELL_TYPE].copy()
    
    return sps_osteo, random_osteo, full_osteo


def compute_statistical_power_metrics(sps_markers, random_markers, full_markers, 
                                       top_n=50, logger=None):
    """
    Compute statistical power metrics comparing SPS vs Random.
    
    Parameters:
    -----------
    sps_markers : pd.DataFrame
        SPS osteoblast markers
    random_markers : pd.DataFrame
        Random osteoblast markers
    full_markers : pd.DataFrame
        Full dataset osteoblast markers
    top_n : int
        Number of top markers to analyze
    logger : logging.Logger
        Logger for output
        
    Returns:
    --------
    dict : Statistical power metrics
    """
    metrics = {}
    
    # Focus on top N markers
    sps_top = sps_markers.head(top_n)
    random_top = random_markers.head(top_n)
    full_top = full_markers.head(top_n)
    
    # 1. P-value analysis (if available)
    if 'pval' in sps_top.columns:
        # Handle zero p-values (replace with minimum float)
        sps_pvals = np.clip(sps_top['pval'].values, np.finfo(float).tiny, 1)
        random_pvals = np.clip(random_top['pval'].values, np.finfo(float).tiny, 1)
        full_pvals = np.clip(full_top['pval'].values, np.finfo(float).tiny, 1)
        
        # -log10(p-value)
        sps_neglogp = -np.log10(sps_pvals)
        random_neglogp = -np.log10(random_pvals)
        full_neglogp = -np.log10(full_pvals)
        
        metrics['pvalue'] = {
            'sps_mean_neglogp': np.mean(sps_neglogp),
            'random_mean_neglogp': np.mean(random_neglogp),
            'full_mean_neglogp': np.mean(full_neglogp),
            'sps_median_neglogp': np.median(sps_neglogp),
            'random_median_neglogp': np.median(random_neglogp),
            'full_median_neglogp': np.median(full_neglogp),
            'sps_neglogp_values': sps_neglogp,
            'random_neglogp_values': random_neglogp,
            'full_neglogp_values': full_neglogp,
        }
        
        # Count markers at different significance thresholds
        for thresh in P_THRESHOLDS:
            sps_count = np.sum(sps_pvals < thresh)
            random_count = np.sum(random_pvals < thresh)
            full_count = np.sum(full_pvals < thresh)
            metrics['pvalue'][f'sps_n_below_{thresh}'] = sps_count
            metrics['pvalue'][f'random_n_below_{thresh}'] = random_count
            metrics['pvalue'][f'full_n_below_{thresh}'] = full_count
        
        if logger:
            logger.info(f"\nP-value analysis (top {top_n} markers):")
            logger.info(f"  SPS mean -log10(p):    {metrics['pvalue']['sps_mean_neglogp']:.2f}")
            logger.info(f"  Random mean -log10(p): {metrics['pvalue']['random_mean_neglogp']:.2f}")
            logger.info(f"  Full mean -log10(p):   {metrics['pvalue']['full_mean_neglogp']:.2f}")
            for thresh in P_THRESHOLDS:
                logger.info(f"  Markers with p < {thresh}:")
                logger.info(f"    SPS: {metrics['pvalue'][f'sps_n_below_{thresh}']}, "
                          f"Random: {metrics['pvalue'][f'random_n_below_{thresh}']}, "
                          f"Full: {metrics['pvalue'][f'full_n_below_{thresh}']}")
    
    # 2. Log fold change analysis (if available)
    if 'logfoldchange' in sps_top.columns:
        sps_lfc = np.abs(sps_top['logfoldchange'].values)
        random_lfc = np.abs(random_top['logfoldchange'].values)
        full_lfc = np.abs(full_top['logfoldchange'].values)
        
        metrics['logfoldchange'] = {
            'sps_mean_abs_lfc': np.mean(sps_lfc),
            'random_mean_abs_lfc': np.mean(random_lfc),
            'full_mean_abs_lfc': np.mean(full_lfc),
            'sps_median_abs_lfc': np.median(sps_lfc),
            'random_median_abs_lfc': np.median(random_lfc),
            'full_median_abs_lfc': np.median(full_lfc),
            'sps_lfc_values': sps_lfc,
            'random_lfc_values': random_lfc,
            'full_lfc_values': full_lfc,
        }
        
        if logger:
            logger.info(f"\nLog Fold Change analysis (top {top_n} markers):")
            logger.info(f"  SPS mean |LFC|:    {metrics['logfoldchange']['sps_mean_abs_lfc']:.3f}")
            logger.info(f"  Random mean |LFC|: {metrics['logfoldchange']['random_mean_abs_lfc']:.3f}")
            logger.info(f"  Full mean |LFC|:   {metrics['logfoldchange']['full_mean_abs_lfc']:.3f}")
    
    # 3. Score analysis
    sps_scores = sps_top['score'].values
    random_scores = random_top['score'].values
    full_scores = full_top['score'].values
    
    metrics['score'] = {
        'sps_mean_score': np.mean(sps_scores),
        'random_mean_score': np.mean(random_scores),
        'full_mean_score': np.mean(full_scores),
        'sps_median_score': np.median(sps_scores),
        'random_median_score': np.median(random_scores),
        'full_median_score': np.median(full_scores),
        'sps_score_values': sps_scores,
        'random_score_values': random_scores,
        'full_score_values': full_scores,
    }
    
    if logger:
        logger.info(f"\nScore analysis (top {top_n} markers):")
        logger.info(f"  SPS mean score:    {metrics['score']['sps_mean_score']:.2f}")
        logger.info(f"  Random mean score: {metrics['score']['random_mean_score']:.2f}")
        logger.info(f"  Full mean score:   {metrics['score']['full_mean_score']:.2f}")
    
    # 4. Percent expression analysis (if available)
    if 'pct_in_group' in sps_top.columns:
        sps_pct = sps_top['pct_in_group'].values * 100  # Convert to percentage
        random_pct = random_top['pct_in_group'].values * 100
        full_pct = full_top['pct_in_group'].values * 100
        
        metrics['expression'] = {
            'sps_mean_pct': np.mean(sps_pct),
            'random_mean_pct': np.mean(random_pct),
            'full_mean_pct': np.mean(full_pct),
            'sps_pct_values': sps_pct,
            'random_pct_values': random_pct,
            'full_pct_values': full_pct,
        }
        
        if logger:
            logger.info(f"\nExpression % analysis (top {top_n} markers):")
            logger.info(f"  SPS mean %:    {metrics['expression']['sps_mean_pct']:.1f}%")
            logger.info(f"  Random mean %: {metrics['expression']['random_mean_pct']:.1f}%")
            logger.info(f"  Full mean %:   {metrics['expression']['full_mean_pct']:.1f}%")
    
    return metrics


def compute_marker_overlap(sps_markers, random_markers, full_markers, top_n=50, logger=None):
    """
    Compute overlap of top marker genes between methods.
    
    Parameters:
    -----------
    sps_markers, random_markers, full_markers : pd.DataFrame
        Marker gene results
    top_n : int
        Number of top markers to compare
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    dict : Overlap statistics
    """
    sps_genes = set(sps_markers.head(top_n)['gene'].values)
    random_genes = set(random_markers.head(top_n)['gene'].values)
    full_genes = set(full_markers.head(top_n)['gene'].values)
    
    overlap = {
        'sps_random_overlap': len(sps_genes & random_genes),
        'sps_full_overlap': len(sps_genes & full_genes),
        'random_full_overlap': len(random_genes & full_genes),
        'all_three_overlap': len(sps_genes & random_genes & full_genes),
        'sps_unique': len(sps_genes - random_genes - full_genes),
        'random_unique': len(random_genes - sps_genes - full_genes),
        'sps_genes': sps_genes,
        'random_genes': random_genes,
        'full_genes': full_genes,
    }
    
    if logger:
        logger.info(f"\nMarker gene overlap (top {top_n}):")
        logger.info(f"  SPS ∩ Random: {overlap['sps_random_overlap']} ({100*overlap['sps_random_overlap']/top_n:.1f}%)")
        logger.info(f"  SPS ∩ Full:   {overlap['sps_full_overlap']} ({100*overlap['sps_full_overlap']/top_n:.1f}%)")
        logger.info(f"  Random ∩ Full: {overlap['random_full_overlap']} ({100*overlap['random_full_overlap']/top_n:.1f}%)")
        logger.info(f"  All three:    {overlap['all_three_overlap']}")
    
    return overlap


def statistical_tests(metrics, logger=None):
    """
    Perform statistical tests comparing SPS vs Random.
    
    Parameters:
    -----------
    metrics : dict
        Metrics from compute_statistical_power_metrics
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    dict : Test results
    """
    tests = {}
    
    # Mann-Whitney U test for -log10(p-values)
    if 'pvalue' in metrics:
        sps_vals = metrics['pvalue']['sps_neglogp_values']
        random_vals = metrics['pvalue']['random_neglogp_values']
        
        stat, pval = stats.mannwhitneyu(sps_vals, random_vals, alternative='greater')
        tests['neglogp_mannwhitney'] = {
            'statistic': stat,
            'pvalue': pval,
            'interpretation': 'SPS > Random' if pval < 0.05 else 'No significant difference'
        }
        
        if logger:
            logger.info(f"\nStatistical tests:")
            logger.info(f"  -log10(p) Mann-Whitney U: stat={stat:.2f}, p={pval:.2e}")
            logger.info(f"    Interpretation: {tests['neglogp_mannwhitney']['interpretation']}")
    
    # Mann-Whitney U test for scores
    if 'score' in metrics:
        sps_scores = metrics['score']['sps_score_values']
        random_scores = metrics['score']['random_score_values']
        
        stat, pval = stats.mannwhitneyu(sps_scores, random_scores, alternative='greater')
        tests['score_mannwhitney'] = {
            'statistic': stat,
            'pvalue': pval,
            'interpretation': 'SPS > Random' if pval < 0.05 else 'No significant difference'
        }
        
        if logger:
            logger.info(f"  Score Mann-Whitney U: stat={stat:.2f}, p={pval:.2e}")
            logger.info(f"    Interpretation: {tests['score_mannwhitney']['interpretation']}")
    
    # Mann-Whitney U test for |LFC|
    if 'logfoldchange' in metrics:
        sps_lfc = metrics['logfoldchange']['sps_lfc_values']
        random_lfc = metrics['logfoldchange']['random_lfc_values']
        
        stat, pval = stats.mannwhitneyu(sps_lfc, random_lfc, alternative='greater')
        tests['lfc_mannwhitney'] = {
            'statistic': stat,
            'pvalue': pval,
            'interpretation': 'SPS > Random' if pval < 0.05 else 'No significant difference'
        }
        
        if logger:
            logger.info(f"  |LFC| Mann-Whitney U: stat={stat:.2f}, p={pval:.2e}")
            logger.info(f"    Interpretation: {tests['lfc_mannwhitney']['interpretation']}")
    
    return tests


# ============================================================================
# Visualization Functions
# ============================================================================

def plot_pvalue_comparison(metrics, de_method, output_dir):
    """Create p-value comparison plots."""
    if 'pvalue' not in metrics:
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # 1. Box plot of -log10(p-values)
    ax = axes[0]
    data = {
        'SPS': metrics['pvalue']['sps_neglogp_values'],
        'Random': metrics['pvalue']['random_neglogp_values'],
        'Full': metrics['pvalue']['full_neglogp_values'],
    }
    positions = [1, 2, 3]
    colors = ['#2ecc71', '#e74c3c', '#3498db']
    
    bp = ax.boxplot([data['SPS'], data['Random'], data['Full']], 
                    positions=positions, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax.set_xticklabels(['SPS', 'Random', 'Full'])
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'Statistical Significance ({de_method})\nOsteoblast Markers')
    ax.grid(True, alpha=0.3)
    
    # 2. Bar chart of significant markers at different thresholds
    ax = axes[1]
    x = np.arange(len(P_THRESHOLDS))
    width = 0.25
    
    sps_counts = [metrics['pvalue'][f'sps_n_below_{t}'] for t in P_THRESHOLDS]
    random_counts = [metrics['pvalue'][f'random_n_below_{t}'] for t in P_THRESHOLDS]
    full_counts = [metrics['pvalue'][f'full_n_below_{t}'] for t in P_THRESHOLDS]
    
    ax.bar(x - width, sps_counts, width, label='SPS', color='#2ecc71', alpha=0.7)
    ax.bar(x, random_counts, width, label='Random', color='#e74c3c', alpha=0.7)
    ax.bar(x + width, full_counts, width, label='Full', color='#3498db', alpha=0.7)
    
    ax.set_xlabel('P-value threshold')
    ax.set_ylabel('Number of markers')
    ax.set_title('Markers Below Significance Threshold')
    ax.set_xticks(x)
    ax.set_xticklabels([f'<{t}' for t in P_THRESHOLDS])
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # 3. Cumulative distribution
    ax = axes[2]
    for name, vals, color in [('SPS', metrics['pvalue']['sps_neglogp_values'], '#2ecc71'),
                               ('Random', metrics['pvalue']['random_neglogp_values'], '#e74c3c'),
                               ('Full', metrics['pvalue']['full_neglogp_values'], '#3498db')]:
        sorted_vals = np.sort(vals)
        cumulative = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        ax.plot(sorted_vals, cumulative, label=name, color=color, linewidth=2)
    
    ax.set_xlabel('-log10(p-value)')
    ax.set_ylabel('Cumulative Proportion')
    ax.set_title('Cumulative Distribution of P-values')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'pvalue_comparison_{de_method}.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_volcano_comparison(sps_markers, random_markers, de_method, output_dir, top_n=50):
    """Create side-by-side volcano plots."""
    if 'logfoldchange' not in sps_markers.columns or 'pval' not in sps_markers.columns:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for ax, (df, title, color) in zip(axes, [(sps_markers, 'SPS', '#2ecc71'),
                                              (random_markers, 'Random', '#e74c3c')]):
        lfc = df['logfoldchange'].values
        pvals = np.clip(df['pval'].values, np.finfo(float).tiny, 1)
        neglogp = -np.log10(pvals)
        
        # Color by significance
        colors = np.where(neglogp > 50, color, '#95a5a6')
        sizes = np.where(neglogp > 50, 50, 20)
        
        ax.scatter(lfc, neglogp, c=colors, s=sizes, alpha=0.6, edgecolors='none')
        
        # Highlight top genes
        top_genes = df.head(10)['gene'].values
        for i, gene in enumerate(df['gene'].head(10).values):
            ax.annotate(gene, (lfc[i], neglogp[i]), fontsize=8, alpha=0.8)
        
        ax.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='p=1e-50')
        ax.axhline(y=100, color='darkred', linestyle='--', alpha=0.5, label='p=1e-100')
        ax.axvline(x=0, color='gray', linestyle='-', alpha=0.3)
        
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-log10(p-value)')
        ax.set_title(f'{title} Sample: Osteoblast Markers ({de_method})')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'volcano_comparison_{de_method}.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_score_comparison(metrics, de_method, output_dir):
    """Create score comparison plots."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # 1. Box plot of scores
    ax = axes[0]
    data = {
        'SPS': metrics['score']['sps_score_values'],
        'Random': metrics['score']['random_score_values'],
        'Full': metrics['score']['full_score_values'],
    }
    colors = ['#2ecc71', '#e74c3c', '#3498db']
    
    bp = ax.boxplot([data['SPS'], data['Random'], data['Full']], 
                    positions=[1, 2, 3], patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax.set_xticklabels(['SPS', 'Random', 'Full'])
    ax.set_ylabel('Marker Score')
    ax.set_title(f'Marker Gene Scores ({de_method})\nOsteoblast')
    ax.grid(True, alpha=0.3)
    
    # 2. Rank comparison
    ax = axes[1]
    ranks = np.arange(1, len(data['SPS']) + 1)
    ax.plot(ranks, np.sort(data['SPS'])[::-1], label='SPS', color='#2ecc71', linewidth=2)
    ax.plot(ranks, np.sort(data['Random'])[::-1], label='Random', color='#e74c3c', linewidth=2)
    ax.plot(ranks, np.sort(data['Full'])[::-1], label='Full', color='#3498db', linewidth=2, linestyle='--')
    
    ax.set_xlabel('Rank')
    ax.set_ylabel('Score')
    ax.set_title('Score by Rank')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'score_comparison_{de_method}.png'), dpi=150, bbox_inches='tight')
    plt.close()


def plot_summary_comparison(all_metrics, output_dir):
    """Create summary comparison across DE methods."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    de_methods = list(all_metrics.keys())
    x = np.arange(len(de_methods))
    width = 0.35
    
    # 1. Mean -log10(p-value)
    ax = axes[0]
    if all(('pvalue' in m) for m in all_metrics.values()):
        sps_vals = [all_metrics[de]['pvalue']['sps_mean_neglogp'] for de in de_methods]
        random_vals = [all_metrics[de]['pvalue']['random_mean_neglogp'] for de in de_methods]
        
        ax.bar(x - width/2, sps_vals, width, label='SPS', color='#2ecc71', alpha=0.7)
        ax.bar(x + width/2, random_vals, width, label='Random', color='#e74c3c', alpha=0.7)
        ax.set_ylabel('Mean -log10(p-value)')
        ax.set_title('Statistical Significance')
        ax.set_xticks(x)
        ax.set_xticklabels(de_methods)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    # 2. Mean score
    ax = axes[1]
    sps_vals = [all_metrics[de]['score']['sps_mean_score'] for de in de_methods]
    random_vals = [all_metrics[de]['score']['random_mean_score'] for de in de_methods]
    
    ax.bar(x - width/2, sps_vals, width, label='SPS', color='#2ecc71', alpha=0.7)
    ax.bar(x + width/2, random_vals, width, label='Random', color='#e74c3c', alpha=0.7)
    ax.set_ylabel('Mean Score')
    ax.set_title('Marker Gene Scores')
    ax.set_xticks(x)
    ax.set_xticklabels(de_methods)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # 3. Mean |LFC|
    ax = axes[2]
    if all(('logfoldchange' in m) for m in all_metrics.values()):
        sps_vals = [all_metrics[de]['logfoldchange']['sps_mean_abs_lfc'] for de in de_methods]
        random_vals = [all_metrics[de]['logfoldchange']['random_mean_abs_lfc'] for de in de_methods]
        
        ax.bar(x - width/2, sps_vals, width, label='SPS', color='#2ecc71', alpha=0.7)
        ax.bar(x + width/2, random_vals, width, label='Random', color='#e74c3c', alpha=0.7)
        ax.set_ylabel('Mean |Log Fold Change|')
        ax.set_title('Effect Size')
        ax.set_xticks(x)
        ax.set_xticklabels(de_methods)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle('Statistical Power Comparison: SPS vs Random (Osteoblast Markers)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'summary_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()


def create_summary_table(all_metrics, all_tests, osteoblast_stats, output_dir):
    """Create summary table of results."""
    rows = []
    
    for de_method, metrics in all_metrics.items():
        row = {'de_method': de_method}
        
        # Score metrics
        row['sps_mean_score'] = metrics['score']['sps_mean_score']
        row['random_mean_score'] = metrics['score']['random_mean_score']
        row['score_improvement'] = (row['sps_mean_score'] - row['random_mean_score']) / row['random_mean_score'] * 100
        
        # P-value metrics (if available)
        if 'pvalue' in metrics:
            row['sps_mean_neglogp'] = metrics['pvalue']['sps_mean_neglogp']
            row['random_mean_neglogp'] = metrics['pvalue']['random_mean_neglogp']
            row['neglogp_improvement'] = (row['sps_mean_neglogp'] - row['random_mean_neglogp']) / row['random_mean_neglogp'] * 100
            
            for thresh in P_THRESHOLDS:
                row[f'sps_n_below_{thresh}'] = metrics['pvalue'][f'sps_n_below_{thresh}']
                row[f'random_n_below_{thresh}'] = metrics['pvalue'][f'random_n_below_{thresh}']
        
        # LFC metrics (if available)
        if 'logfoldchange' in metrics:
            row['sps_mean_lfc'] = metrics['logfoldchange']['sps_mean_abs_lfc']
            row['random_mean_lfc'] = metrics['logfoldchange']['random_mean_abs_lfc']
        
        # Test results
        if de_method in all_tests:
            tests = all_tests[de_method]
            if 'score_mannwhitney' in tests:
                row['score_test_pval'] = tests['score_mannwhitney']['pvalue']
            if 'neglogp_mannwhitney' in tests:
                row['neglogp_test_pval'] = tests['neglogp_mannwhitney']['pvalue']
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Add osteoblast counts
    summary_text = f"""
Statistical Power Analysis Summary
==================================

Osteoblast Cell Counts:
  Full dataset: {osteoblast_stats['full']['osteoblast_count']:,} cells ({osteoblast_stats['full']['osteoblast_pct']:.2f}%)
  SPS sample:   {osteoblast_stats['sps']['osteoblast_count']:,} cells ({osteoblast_stats['sps']['osteoblast_pct']:.2f}%, {osteoblast_stats['sps']['enrichment_vs_full']:.2f}x enrichment)
  Random sample: {osteoblast_stats['random']['osteoblast_count']:,} cells ({osteoblast_stats['random']['osteoblast_pct']:.2f}%, {osteoblast_stats['random']['enrichment_vs_full']:.2f}x enrichment)

Key Findings:
"""
    
    for de_method, metrics in all_metrics.items():
        summary_text += f"\n{de_method}:\n"
        summary_text += f"  Score improvement: {(metrics['score']['sps_mean_score'] - metrics['score']['random_mean_score']) / metrics['score']['random_mean_score'] * 100:.1f}%\n"
        if 'pvalue' in metrics:
            summary_text += f"  -log10(p) improvement: {(metrics['pvalue']['sps_mean_neglogp'] - metrics['pvalue']['random_mean_neglogp']) / metrics['pvalue']['random_mean_neglogp'] * 100:.1f}%\n"
    
    # Save
    df.to_csv(os.path.join(output_dir, 'statistical_power_summary.csv'), index=False)
    
    with open(os.path.join(output_dir, 'statistical_power_summary.txt'), 'w') as f:
        f.write(summary_text)
        f.write("\n\nDetailed metrics:\n")
        f.write(df.to_string())
    
    return df, summary_text


# ============================================================================
# Main Analysis
# ============================================================================

def run_analysis(de_methods=None, top_n=50, logger=None):
    """
    Run the complete statistical power analysis.
    
    Parameters:
    -----------
    de_methods : list
        DE methods to analyze
    top_n : int
        Number of top markers to analyze
    logger : logging.Logger
        Logger
        
    Returns:
    --------
    dict : All results
    """
    if de_methods is None:
        de_methods = DE_METHODS
    
    if logger is None:
        logger = setup_logger('statistical_power_analysis')
    
    # Create output directories
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)
    
    logger.info("=" * 80)
    logger.info("Experiment 1: Statistical Power Analysis for Rare Cell Markers")
    logger.info("=" * 80)
    
    # Load cell count statistics
    logger.info("\n" + "=" * 40)
    logger.info("Step 1: Computing osteoblast statistics")
    logger.info("=" * 40)
    
    try:
        adata = load_mcc_dataset(logger)
        sps_indices = load_sampling_indices('sps', logger)
        random_indices = load_sampling_indices('random', logger)
        
        osteoblast_stats = compute_osteoblast_statistics(
            adata, sps_indices, random_indices, logger
        )
    except Exception as e:
        logger.warning(f"Could not compute osteoblast statistics: {e}")
        osteoblast_stats = None
    
    # Analyze each DE method
    all_metrics = {}
    all_tests = {}
    all_overlaps = {}
    
    for de_method in de_methods:
        logger.info("\n" + "=" * 40)
        logger.info(f"Step 2: Analyzing {de_method} results")
        logger.info("=" * 40)
        
        try:
            # Load markers
            sps_markers, random_markers, full_markers = load_osteoblast_markers(de_method)
            logger.info(f"Loaded {len(sps_markers)} SPS markers, {len(random_markers)} Random markers")
            
            # Compute metrics
            metrics = compute_statistical_power_metrics(
                sps_markers, random_markers, full_markers, 
                top_n=top_n, logger=logger
            )
            all_metrics[de_method] = metrics
            
            # Compute overlap
            overlap = compute_marker_overlap(
                sps_markers, random_markers, full_markers,
                top_n=top_n, logger=logger
            )
            all_overlaps[de_method] = overlap
            
            # Statistical tests
            tests = statistical_tests(metrics, logger)
            all_tests[de_method] = tests
            
            # Create visualizations
            logger.info(f"\nCreating visualizations for {de_method}...")
            plot_pvalue_comparison(metrics, de_method, FIGURES_DIR)
            plot_volcano_comparison(sps_markers, random_markers, de_method, FIGURES_DIR, top_n)
            plot_score_comparison(metrics, de_method, FIGURES_DIR)
            
        except Exception as e:
            logger.error(f"Failed to analyze {de_method}: {e}")
            import traceback
            logger.error(traceback.format_exc())
    
    # Create summary visualizations and table
    logger.info("\n" + "=" * 40)
    logger.info("Step 3: Creating summary")
    logger.info("=" * 40)
    
    if all_metrics:
        plot_summary_comparison(all_metrics, FIGURES_DIR)
        
        if osteoblast_stats:
            summary_df, summary_text = create_summary_table(
                all_metrics, all_tests, osteoblast_stats, OUTPUT_DIR
            )
            logger.info("\n" + summary_text)
    
    # Save all results
    results = {
        'osteoblast_stats': osteoblast_stats,
        'metrics': all_metrics,
        'tests': all_tests,
        'overlaps': all_overlaps,
    }
    
    with open(os.path.join(OUTPUT_DIR, 'statistical_power_results.pkl'), 'wb') as f:
        pickle.dump(results, f)
    
    logger.info("\n" + "=" * 80)
    logger.info("Analysis complete!")
    logger.info(f"Results saved to: {OUTPUT_DIR}")
    logger.info(f"Figures saved to: {FIGURES_DIR}")
    logger.info("=" * 80)
    
    return results


# ============================================================================
# Entry Point
# ============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Statistical Power Analysis for Rare Cell Markers'
    )
    parser.add_argument(
        '--de-methods',
        type=str,
        nargs='+',
        default=DE_METHODS,
        help=f'DE methods to analyze. Options: {DE_METHODS}'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=50,
        help='Number of top markers to analyze (default: 50)'
    )
    
    args = parser.parse_args()
    
    logger = setup_logger('statistical_power_analysis')
    
    try:
        results = run_analysis(
            de_methods=args.de_methods,
            top_n=args.top_n,
            logger=logger
        )
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
