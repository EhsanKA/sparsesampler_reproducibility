#!/usr/bin/env python
# coding: utf-8

"""
Compare Marker Genes: Full vs SPS vs Random, across DE methods

This script compares marker gene results from different sampling methods and DE methods,
generating visualizations demonstrating that SPS better captures literature-based
marker genes compared to random sampling.

Usage:
    python compare_marker_genes.py
    python compare_marker_genes.py --de-methods wilcoxon t-test logreg

Requirements:
    - Results from scanpy_marker_pipeline.py for full, sps, and random methods
    - matplotlib, seaborn, pandas, numpy
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

# Add parent directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
analysis_dir = os.path.dirname(script_dir)
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

from semitones_marker_genes.literature_markers import LITERATURE_MARKERS, get_all_markers

# ============================================================================
# Configuration
# ============================================================================

RESULTS_DIR = os.path.join(script_dir, 'results')
FIGURES_DIR = os.path.join(RESULTS_DIR, 'figures')

# Sampling methods to compare
SAMPLING_METHODS = ['full', 'sps', 'random']

# Default DE methods
DEFAULT_DE_METHODS = ['wilcoxon', 't-test', 'logreg']

# Ensure directories exist
os.makedirs(FIGURES_DIR, exist_ok=True)

# Plot settings
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
sns.set_style('whitegrid')

# Colors for sampling methods
METHOD_COLORS = {
    'full': '#3498db',   # Blue
    'sps': '#2ecc71',    # Green
    'random': '#e74c3c', # Red
}

# Colors for DE methods
DE_COLORS = {
    'wilcoxon': '#9b59b6',      # Purple
    't-test': '#e67e22',        # Orange
    't-test_overestim_var': '#f39c12',  # Yellow
    'logreg': '#1abc9c',        # Teal
}


# ============================================================================
# Data Loading Functions
# ============================================================================

def load_results(sampling_method, de_method):
    """Load results for a given sampling method and DE method."""
    results = {}
    prefix = f"{sampling_method}_{de_method}"
    
    # Load marker genes
    markers_path = os.path.join(RESULTS_DIR, f'{prefix}_marker_genes.pkl')
    if os.path.exists(markers_path):
        with open(markers_path, 'rb') as f:
            results['marker_genes'] = pickle.load(f)
    
    # Load evaluation
    eval_path = os.path.join(RESULTS_DIR, f'{prefix}_literature_evaluation.pkl')
    if os.path.exists(eval_path):
        with open(eval_path, 'rb') as f:
            results['evaluation'] = pickle.load(f)
    
    return results if results else None


def load_evaluation_csv(sampling_method, de_method):
    """Load evaluation results from CSV."""
    prefix = f"{sampling_method}_{de_method}"
    csv_path = os.path.join(RESULTS_DIR, f'{prefix}_literature_evaluation.csv')
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return None


def load_marker_genes_csv(sampling_method, de_method):
    """Load marker genes from CSV."""
    prefix = f"{sampling_method}_{de_method}"
    csv_path = os.path.join(RESULTS_DIR, f'{prefix}_marker_genes.csv')
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return None


def load_all_results(de_methods):
    """Load all results for all sampling methods and DE methods."""
    all_results = {}
    
    for sampling_method in SAMPLING_METHODS:
        all_results[sampling_method] = {}
        for de_method in de_methods:
            results = load_results(sampling_method, de_method)
            if results:
                all_results[sampling_method][de_method] = results
                print(f"Loaded: {sampling_method}/{de_method}")
            else:
                print(f"Missing: {sampling_method}/{de_method}")
    
    return all_results


# ============================================================================
# Comparison Functions
# ============================================================================

def compare_sampling_methods(all_results, de_method):
    """
    Compare recall metrics between sampling methods for a given DE method.
    
    Returns:
    --------
    pd.DataFrame : Comparison table
    """
    rows = []
    
    # Get cell types from any available result
    cell_types = None
    for sampling_method in SAMPLING_METHODS:
        if de_method in all_results.get(sampling_method, {}):
            eval_results = all_results[sampling_method][de_method].get('evaluation', {})
            cell_types = list(eval_results.keys())
            break
    
    if cell_types is None:
        return pd.DataFrame()
    
    for celltype in cell_types:
        row = {'cell_type': celltype, 'de_method': de_method}
        
        for sampling_method in SAMPLING_METHODS:
            if de_method in all_results.get(sampling_method, {}):
                eval_results = all_results[sampling_method][de_method].get('evaluation', {})
                if celltype in eval_results:
                    metrics = eval_results[celltype]
                    row[f'{sampling_method}_recall_at_10'] = metrics.get('recall_at_10')
                    row[f'{sampling_method}_recall_at_50'] = metrics.get('recall_at_50')
                    row[f'{sampling_method}_mean_rank'] = metrics.get('mean_rank')
                    row[f'{sampling_method}_n_markers'] = metrics.get('n_literature_markers_in_dataset', 0)
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    
    # Calculate improvements (SPS over Random)
    if 'sps_recall_at_10' in df.columns and 'random_recall_at_10' in df.columns:
        df['recall_10_improvement'] = df['sps_recall_at_10'] - df['random_recall_at_10']
    if 'sps_recall_at_50' in df.columns and 'random_recall_at_50' in df.columns:
        df['recall_50_improvement'] = df['sps_recall_at_50'] - df['random_recall_at_50']
    if 'sps_mean_rank' in df.columns and 'random_mean_rank' in df.columns:
        df['rank_improvement'] = df['random_mean_rank'] - df['sps_mean_rank']  # Lower is better
    
    return df


def compare_de_methods(all_results, sampling_method):
    """
    Compare recall metrics between DE methods for a given sampling method.
    
    Returns:
    --------
    pd.DataFrame : Comparison table
    """
    rows = []
    
    de_methods = list(all_results.get(sampling_method, {}).keys())
    if not de_methods:
        return pd.DataFrame()
    
    # Get cell types from first DE method
    first_de = de_methods[0]
    eval_results = all_results[sampling_method][first_de].get('evaluation', {})
    cell_types = list(eval_results.keys())
    
    for celltype in cell_types:
        row = {'cell_type': celltype, 'sampling_method': sampling_method}
        
        for de_method in de_methods:
            if de_method in all_results[sampling_method]:
                eval_results = all_results[sampling_method][de_method].get('evaluation', {})
                if celltype in eval_results:
                    metrics = eval_results[celltype]
                    row[f'{de_method}_recall_at_10'] = metrics.get('recall_at_10')
                    row[f'{de_method}_recall_at_50'] = metrics.get('recall_at_50')
                    row[f'{de_method}_mean_rank'] = metrics.get('mean_rank')
        
        rows.append(row)
    
    return pd.DataFrame(rows)


def compare_marker_overlap(all_results, de_method, top_n=50):
    """
    Compare overlap of top marker genes between sampling methods.
    
    Returns:
    --------
    dict : Overlap metrics per cell type
    """
    overlap_results = {}
    
    # Get marker genes for each sampling method
    markers_by_method = {}
    for sampling_method in SAMPLING_METHODS:
        if de_method in all_results.get(sampling_method, {}):
            markers_by_method[sampling_method] = all_results[sampling_method][de_method].get('marker_genes', {})
    
    if len(markers_by_method) < 2:
        return {}
    
    # Get cell types
    cell_types = list(next(iter(markers_by_method.values())).keys())
    
    for celltype in cell_types:
        overlap_results[celltype] = {}
        
        for method1 in markers_by_method:
            for method2 in markers_by_method:
                if method1 >= method2:
                    continue
                
                top1 = set(markers_by_method[method1].get(celltype, {}).get('top_genes', [])[:top_n])
                top2 = set(markers_by_method[method2].get(celltype, {}).get('top_genes', [])[:top_n])
                
                if not top1 or not top2:
                    continue
                
                intersection = top1 & top2
                union = top1 | top2
                jaccard = len(intersection) / len(union) if union else 0
                
                key = f'{method1}_vs_{method2}'
                overlap_results[celltype][key] = {
                    'n_overlap': len(intersection),
                    'jaccard': jaccard,
                    'n_method1_unique': len(top1 - top2),
                    'n_method2_unique': len(top2 - top1),
                }
    
    return overlap_results


def compute_rank_correlation(all_results, de_method):
    """
    Compute Spearman rank correlation of gene rankings between sampling methods.
    
    Returns:
    --------
    dict : Correlation metrics per cell type
    """
    correlations = {}
    
    # Get marker genes for each sampling method
    markers_by_method = {}
    for sampling_method in SAMPLING_METHODS:
        if de_method in all_results.get(sampling_method, {}):
            markers_by_method[sampling_method] = all_results[sampling_method][de_method].get('marker_genes', {})
    
    if len(markers_by_method) < 2:
        return {}
    
    cell_types = list(next(iter(markers_by_method.values())).keys())
    
    for celltype in cell_types:
        correlations[celltype] = {}
        
        for method1 in markers_by_method:
            for method2 in markers_by_method:
                if method1 >= method2:
                    continue
                
                ranking1 = markers_by_method[method1].get(celltype, {}).get('gene_ranking', [])
                ranking2 = markers_by_method[method2].get(celltype, {}).get('gene_ranking', [])
                
                if not ranking1 or not ranking2:
                    continue
                
                # Get common genes
                common_genes = list(set(ranking1) & set(ranking2))
                
                if len(common_genes) < 10:
                    continue
                
                # Get ranks for common genes
                ranks1 = [ranking1.index(g) + 1 for g in common_genes]
                ranks2 = [ranking2.index(g) + 1 for g in common_genes]
                
                corr, pval = stats.spearmanr(ranks1, ranks2)
                
                key = f'{method1}_vs_{method2}'
                correlations[celltype][key] = {
                    'spearman_r': corr,
                    'p_value': pval,
                    'n_common_genes': len(common_genes),
                }
    
    return correlations


# ============================================================================
# Visualization Functions
# ============================================================================

def plot_recall_comparison_by_sampling(comparison_df, de_method, save_path=None):
    """
    Plot comparison of recall metrics between sampling methods for one DE method.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    cell_types = comparison_df['cell_type'].tolist()
    x = np.arange(len(cell_types))
    width = 0.25
    
    # Recall@10
    ax1 = axes[0]
    for i, method in enumerate(SAMPLING_METHODS):
        col = f'{method}_recall_at_10'
        if col in comparison_df.columns:
            values = comparison_df[col].fillna(0).values
            bars = ax1.bar(x + (i - 1) * width, values, width, 
                          label=method.upper(), color=METHOD_COLORS.get(method, 'gray'), alpha=0.8)
    
    ax1.set_ylabel('Recall@10')
    ax1.set_title(f'Literature Marker Recovery (Top 10) - {de_method}')
    ax1.set_xticks(x)
    ax1.set_xticklabels(cell_types, rotation=45, ha='right')
    ax1.legend()
    ax1.set_ylim(0, 1.1)
    
    # Recall@50
    ax2 = axes[1]
    for i, method in enumerate(SAMPLING_METHODS):
        col = f'{method}_recall_at_50'
        if col in comparison_df.columns:
            values = comparison_df[col].fillna(0).values
            bars = ax2.bar(x + (i - 1) * width, values, width,
                          label=method.upper(), color=METHOD_COLORS.get(method, 'gray'), alpha=0.8)
    
    ax2.set_ylabel('Recall@50')
    ax2.set_title(f'Literature Marker Recovery (Top 50) - {de_method}')
    ax2.set_xticks(x)
    ax2.set_xticklabels(cell_types, rotation=45, ha='right')
    ax2.legend()
    ax2.set_ylim(0, 1.1)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


def plot_rank_comparison_by_sampling(comparison_df, de_method, save_path=None):
    """
    Plot comparison of mean ranks for literature markers.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    cell_types = comparison_df['cell_type'].tolist()
    x = np.arange(len(cell_types))
    width = 0.25
    
    for i, method in enumerate(SAMPLING_METHODS):
        col = f'{method}_mean_rank'
        if col in comparison_df.columns:
            values = comparison_df[col].fillna(0).values
            ax.bar(x + (i - 1) * width, values, width,
                   label=method.upper(), color=METHOD_COLORS.get(method, 'gray'), alpha=0.8)
    
    ax.set_ylabel('Mean Rank of Literature Markers (lower is better)')
    ax.set_title(f'Literature Marker Ranking - {de_method}')
    ax.set_xticks(x)
    ax.set_xticklabels(cell_types, rotation=45, ha='right')
    ax.legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


def plot_improvement_summary(comparison_df, de_method, save_path=None):
    """
    Plot summary of SPS improvement over Random.
    """
    # Filter out rows without improvement data
    df = comparison_df.dropna(subset=['recall_10_improvement', 'recall_50_improvement'])
    
    if len(df) == 0:
        print("No improvement data available")
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Recall@10 improvement
    ax1 = axes[0]
    colors = ['#2ecc71' if x >= 0 else '#e74c3c' for x in df['recall_10_improvement']]
    ax1.bar(df['cell_type'], df['recall_10_improvement'] * 100, color=colors, alpha=0.8)
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax1.set_ylabel('Improvement (%)')
    ax1.set_title('SPS vs Random: Recall@10 Improvement')
    ax1.tick_params(axis='x', rotation=45)
    
    # Recall@50 improvement
    ax2 = axes[1]
    colors = ['#2ecc71' if x >= 0 else '#e74c3c' for x in df['recall_50_improvement']]
    ax2.bar(df['cell_type'], df['recall_50_improvement'] * 100, color=colors, alpha=0.8)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax2.set_ylabel('Improvement (%)')
    ax2.set_title('SPS vs Random: Recall@50 Improvement')
    ax2.tick_params(axis='x', rotation=45)
    
    # Rank improvement
    if 'rank_improvement' in df.columns:
        df_rank = df.dropna(subset=['rank_improvement'])
        if len(df_rank) > 0:
            ax3 = axes[2]
            colors = ['#2ecc71' if x >= 0 else '#e74c3c' for x in df_rank['rank_improvement']]
            ax3.bar(df_rank['cell_type'], df_rank['rank_improvement'], color=colors, alpha=0.8)
            ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
            ax3.set_ylabel('Rank Improvement (higher = better)')
            ax3.set_title('SPS vs Random: Mean Rank Improvement')
            ax3.tick_params(axis='x', rotation=45)
    
    plt.suptitle(f'SPS Improvement over Random Sampling - {de_method}', fontsize=14, y=1.02)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


def plot_de_method_comparison(all_results, sampling_method, save_path=None):
    """
    Plot comparison of DE methods for a given sampling method.
    """
    comparison_df = compare_de_methods(all_results, sampling_method)
    
    if comparison_df.empty:
        print(f"No data for {sampling_method}")
        return
    
    de_methods = [col.replace('_recall_at_10', '') for col in comparison_df.columns 
                  if '_recall_at_10' in col]
    
    if not de_methods:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    cell_types = comparison_df['cell_type'].tolist()
    x = np.arange(len(cell_types))
    width = 0.8 / len(de_methods)
    
    # Recall@10
    ax1 = axes[0]
    for i, de_method in enumerate(de_methods):
        col = f'{de_method}_recall_at_10'
        if col in comparison_df.columns:
            values = comparison_df[col].fillna(0).values
            offset = (i - len(de_methods)/2 + 0.5) * width
            ax1.bar(x + offset, values, width,
                   label=de_method, color=DE_COLORS.get(de_method, 'gray'), alpha=0.8)
    
    ax1.set_ylabel('Recall@10')
    ax1.set_title(f'DE Method Comparison (Top 10) - {sampling_method.upper()}')
    ax1.set_xticks(x)
    ax1.set_xticklabels(cell_types, rotation=45, ha='right')
    ax1.legend()
    ax1.set_ylim(0, 1.1)
    
    # Recall@50
    ax2 = axes[1]
    for i, de_method in enumerate(de_methods):
        col = f'{de_method}_recall_at_50'
        if col in comparison_df.columns:
            values = comparison_df[col].fillna(0).values
            offset = (i - len(de_methods)/2 + 0.5) * width
            ax2.bar(x + offset, values, width,
                   label=de_method, color=DE_COLORS.get(de_method, 'gray'), alpha=0.8)
    
    ax2.set_ylabel('Recall@50')
    ax2.set_title(f'DE Method Comparison (Top 50) - {sampling_method.upper()}')
    ax2.set_xticks(x)
    ax2.set_xticklabels(cell_types, rotation=45, ha='right')
    ax2.legend()
    ax2.set_ylim(0, 1.1)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


def plot_combined_heatmap(all_results, de_methods, save_path=None):
    """
    Plot heatmap of recall@50 across all sampling methods and DE methods.
    """
    # Collect data
    data = []
    
    for sampling_method in SAMPLING_METHODS:
        for de_method in de_methods:
            if de_method in all_results.get(sampling_method, {}):
                eval_results = all_results[sampling_method][de_method].get('evaluation', {})
                for celltype, metrics in eval_results.items():
                    recall = metrics.get('recall_at_50')
                    if recall is not None:
                        data.append({
                            'sampling_method': sampling_method,
                            'de_method': de_method,
                            'cell_type': celltype,
                            'recall_at_50': recall
                        })
    
    if not data:
        print("No data for heatmap")
        return
    
    df = pd.DataFrame(data)
    
    # Create pivot table
    df['method_combo'] = df['sampling_method'] + '_' + df['de_method']
    pivot = df.pivot(index='cell_type', columns='method_combo', values='recall_at_50')
    
    # Plot heatmap
    fig, ax = plt.subplots(figsize=(14, 8))
    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='RdYlGn', ax=ax,
                vmin=0, vmax=1, cbar_kws={'label': 'Recall@50'})
    ax.set_title('Literature Marker Recovery (Recall@50) - All Methods')
    ax.set_xlabel('Sampling Method + DE Method')
    ax.set_ylabel('Cell Type')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


def plot_overlap_heatmap(overlap_results, de_method, save_path=None):
    """
    Plot heatmap of marker gene overlap (Jaccard) between sampling methods.
    """
    # Collect Jaccard values
    cell_types = list(overlap_results.keys())
    comparison_pairs = []
    
    for celltype in cell_types:
        for pair in overlap_results[celltype].keys():
            if pair not in comparison_pairs:
                comparison_pairs.append(pair)
    
    if not comparison_pairs:
        print("No overlap data")
        return
    
    data = np.zeros((len(cell_types), len(comparison_pairs)))
    
    for i, celltype in enumerate(cell_types):
        for j, pair in enumerate(comparison_pairs):
            if pair in overlap_results[celltype]:
                data[i, j] = overlap_results[celltype][pair]['jaccard']
    
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(data, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax,
                xticklabels=comparison_pairs, yticklabels=cell_types,
                vmin=0, vmax=1)
    ax.set_title(f'Marker Gene Overlap (Jaccard) - {de_method}')
    ax.set_xlabel('Comparison Pair')
    ax.set_ylabel('Cell Type')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


# ============================================================================
# Main Comparison Function
# ============================================================================

def run_comparison(de_methods=None):
    """
    Run full comparison between sampling methods and DE methods.
    """
    if de_methods is None:
        de_methods = DEFAULT_DE_METHODS
    
    print("=" * 80)
    print("Comparing Marker Gene Results")
    print(f"Sampling methods: {SAMPLING_METHODS}")
    print(f"DE methods: {de_methods}")
    print("=" * 80)
    
    # Load all results
    print("\nLoading results...")
    all_results = load_all_results(de_methods)
    
    # Check if we have any results
    n_results = sum(len(v) for v in all_results.values())
    if n_results == 0:
        print("ERROR: No results found. Please run scanpy_marker_pipeline.py first.")
        print(f"  Expected files in: {RESULTS_DIR}")
        return
    
    print(f"\nLoaded {n_results} result sets")
    
    # Create comparison tables for each DE method
    all_comparisons = {}
    
    for de_method in de_methods:
        print(f"\n{'='*40}")
        print(f"Comparing sampling methods for {de_method}")
        print('='*40)
        
        comparison_df = compare_sampling_methods(all_results, de_method)
        
        if comparison_df.empty:
            print(f"  No data for {de_method}")
            continue
        
        all_comparisons[de_method] = comparison_df
        
        # Save comparison table
        comparison_path = os.path.join(RESULTS_DIR, f'comparison_{de_method}.csv')
        comparison_df.to_csv(comparison_path, index=False)
        print(f"Saved comparison table to {comparison_path}")
        
        # Compute marker overlap
        overlap_results = compare_marker_overlap(all_results, de_method, top_n=50)
        
        if overlap_results:
            overlap_path = os.path.join(RESULTS_DIR, f'overlap_{de_method}.pkl')
            with open(overlap_path, 'wb') as f:
                pickle.dump(overlap_results, f)
            print(f"Saved overlap results to {overlap_path}")
        
        # Compute rank correlations
        correlations = compute_rank_correlation(all_results, de_method)
        
        if correlations:
            corr_path = os.path.join(RESULTS_DIR, f'correlations_{de_method}.pkl')
            with open(corr_path, 'wb') as f:
                pickle.dump(correlations, f)
            print(f"Saved correlations to {corr_path}")
        
        # Generate visualizations
        print(f"\nGenerating visualizations for {de_method}...")
        
        plot_recall_comparison_by_sampling(
            comparison_df, de_method,
            save_path=os.path.join(FIGURES_DIR, f'recall_comparison_{de_method}.png')
        )
        
        plot_rank_comparison_by_sampling(
            comparison_df, de_method,
            save_path=os.path.join(FIGURES_DIR, f'rank_comparison_{de_method}.png')
        )
        
        plot_improvement_summary(
            comparison_df, de_method,
            save_path=os.path.join(FIGURES_DIR, f'sps_improvement_{de_method}.png')
        )
        
        if overlap_results:
            plot_overlap_heatmap(
                overlap_results, de_method,
                save_path=os.path.join(FIGURES_DIR, f'overlap_heatmap_{de_method}.png')
            )
    
    # Cross-DE method comparisons
    print(f"\n{'='*40}")
    print("Cross-DE method comparisons")
    print('='*40)
    
    for sampling_method in SAMPLING_METHODS:
        if all_results.get(sampling_method):
            plot_de_method_comparison(
                all_results, sampling_method,
                save_path=os.path.join(FIGURES_DIR, f'de_comparison_{sampling_method}.png')
            )
    
    # Combined heatmap
    plot_combined_heatmap(
        all_results, de_methods,
        save_path=os.path.join(FIGURES_DIR, 'combined_recall_heatmap.png')
    )
    
    # Save combined comparison table
    if all_comparisons:
        combined_df = pd.concat(all_comparisons.values(), ignore_index=True)
        combined_path = os.path.join(RESULTS_DIR, 'comparison_all.csv')
        combined_df.to_csv(combined_path, index=False)
        print(f"\nSaved combined comparison to {combined_path}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)
    
    for de_method, comparison_df in all_comparisons.items():
        print(f"\n{de_method.upper()}:")
        print("-" * 40)
        
        # Print recall comparison
        cols_to_show = ['cell_type']
        for method in SAMPLING_METHODS:
            if f'{method}_recall_at_50' in comparison_df.columns:
                cols_to_show.append(f'{method}_recall_at_50')
        
        if len(cols_to_show) > 1:
            print(comparison_df[cols_to_show].to_string(index=False))
        
        # Print mean improvements
        if 'recall_50_improvement' in comparison_df.columns:
            mean_improvement = comparison_df['recall_50_improvement'].mean()
            if not np.isnan(mean_improvement):
                print(f"\nMean SPS improvement (Recall@50): {mean_improvement:.2%}")
    
    print("\n" + "=" * 80)
    print("Comparison completed!")
    print(f"Results saved to: {RESULTS_DIR}")
    print(f"Figures saved to: {FIGURES_DIR}")
    print("=" * 80)


# ============================================================================
# Entry Point
# ============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare Scanpy marker gene results across methods'
    )
    parser.add_argument(
        '--de-methods',
        type=str,
        nargs='+',
        default=DEFAULT_DE_METHODS,
        help=f'DE methods to compare (default: {DEFAULT_DE_METHODS})'
    )
    
    args = parser.parse_args()
    
    run_comparison(de_methods=args.de_methods)
