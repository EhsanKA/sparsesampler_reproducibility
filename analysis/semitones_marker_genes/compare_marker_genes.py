#!/usr/bin/env python
# coding: utf-8

"""
Compare Marker Genes: SPS vs Random

This script compares marker gene results from SPS and Random sampling methods
and generates visualizations demonstrating that SPS better captures literature-based
marker genes.

Usage:
    python compare_marker_genes.py

Requirements:
    - Results from semitones_marker_pipeline.py for both sps and random methods
    - matplotlib, seaborn, pandas, numpy
"""

import os
import sys
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

# Ensure directories exist
os.makedirs(FIGURES_DIR, exist_ok=True)

# Plot settings
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
sns.set_style('whitegrid')


# ============================================================================
# Data Loading Functions
# ============================================================================

def load_results(method):
    """Load results for a given method."""
    results = {}
    
    # Load marker genes
    markers_path = os.path.join(RESULTS_DIR, f'{method}_marker_genes.pkl')
    if os.path.exists(markers_path):
        with open(markers_path, 'rb') as f:
            results['marker_genes'] = pickle.load(f)
    
    # Load evaluation
    eval_path = os.path.join(RESULTS_DIR, f'{method}_literature_evaluation.pkl')
    if os.path.exists(eval_path):
        with open(eval_path, 'rb') as f:
            results['evaluation'] = pickle.load(f)
    
    # Load enrichment scores (optional, may be large)
    escores_path = os.path.join(RESULTS_DIR, f'{method}_enrichment_scores.pkl')
    if os.path.exists(escores_path):
        results['enrichment_scores_path'] = escores_path
    
    return results


def load_evaluation_csv(method):
    """Load evaluation results from CSV."""
    csv_path = os.path.join(RESULTS_DIR, f'{method}_literature_evaluation.csv')
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return None


def load_marker_genes_csv(method):
    """Load marker genes from CSV."""
    csv_path = os.path.join(RESULTS_DIR, f'{method}_marker_genes_by_celltype.csv')
    if os.path.exists(csv_path):
        return pd.read_csv(csv_path)
    return None


# ============================================================================
# Comparison Functions
# ============================================================================

def compare_recall_metrics(sps_eval, random_eval):
    """
    Compare recall metrics between SPS and Random methods.
    
    Returns:
    --------
    pd.DataFrame : Comparison table
    """
    rows = []
    
    for celltype in sps_eval.keys():
        if celltype not in random_eval:
            continue
        
        sps_metrics = sps_eval[celltype]
        random_metrics = random_eval[celltype]
        
        rows.append({
            'cell_type': celltype,
            'sps_recall_at_10': sps_metrics.get('recall_at_10'),
            'random_recall_at_10': random_metrics.get('recall_at_10'),
            'sps_recall_at_50': sps_metrics.get('recall_at_50'),
            'random_recall_at_50': random_metrics.get('recall_at_50'),
            'sps_mean_rank': sps_metrics.get('mean_rank'),
            'random_mean_rank': random_metrics.get('mean_rank'),
            'n_literature_markers': sps_metrics.get('n_literature_markers_in_dataset', 0),
        })
    
    df = pd.DataFrame(rows)
    
    # Calculate improvement
    df['recall_10_improvement'] = df['sps_recall_at_10'] - df['random_recall_at_10']
    df['recall_50_improvement'] = df['sps_recall_at_50'] - df['random_recall_at_50']
    df['rank_improvement'] = df['random_mean_rank'] - df['sps_mean_rank']  # Lower is better
    
    return df


def compare_marker_overlap(sps_markers, random_markers, top_n=50):
    """
    Compare overlap of top marker genes between SPS and Random.
    
    Returns:
    --------
    dict : Overlap metrics per cell type
    """
    overlap_results = {}
    
    for celltype in sps_markers.keys():
        if celltype not in random_markers:
            continue
        
        sps_top = set(sps_markers[celltype]['top_genes'][:top_n])
        random_top = set(random_markers[celltype]['top_genes'][:top_n])
        
        intersection = sps_top & random_top
        union = sps_top | random_top
        
        jaccard = len(intersection) / len(union) if union else 0
        
        # Genes unique to each method
        sps_unique = sps_top - random_top
        random_unique = random_top - sps_top
        
        overlap_results[celltype] = {
            'n_sps': len(sps_top),
            'n_random': len(random_top),
            'n_overlap': len(intersection),
            'jaccard_similarity': jaccard,
            'overlapping_genes': list(intersection),
            'sps_unique_genes': list(sps_unique),
            'random_unique_genes': list(random_unique),
        }
    
    return overlap_results


def compute_rank_correlation(sps_markers, random_markers):
    """
    Compute Spearman rank correlation of gene rankings between methods.
    
    Returns:
    --------
    dict : Correlation metrics per cell type
    """
    correlations = {}
    
    for celltype in sps_markers.keys():
        if celltype not in random_markers:
            continue
        
        sps_ranking = sps_markers[celltype]['gene_ranking']
        random_ranking = random_markers[celltype]['gene_ranking']
        
        # Get common genes
        common_genes = set(sps_ranking) & set(random_ranking)
        
        if len(common_genes) < 10:
            correlations[celltype] = {'spearman_r': None, 'p_value': None}
            continue
        
        # Get ranks for common genes
        sps_ranks = [sps_ranking.index(g) + 1 for g in common_genes]
        random_ranks = [random_ranking.index(g) + 1 for g in common_genes]
        
        # Compute Spearman correlation
        corr, pval = stats.spearmanr(sps_ranks, random_ranks)
        
        correlations[celltype] = {
            'spearman_r': corr,
            'p_value': pval,
            'n_common_genes': len(common_genes),
        }
    
    return correlations


# ============================================================================
# Visualization Functions
# ============================================================================

def plot_recall_comparison(comparison_df, save_path=None):
    """
    Plot comparison of recall metrics between SPS and Random.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Recall@10
    ax1 = axes[0]
    x = np.arange(len(comparison_df))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, comparison_df['sps_recall_at_10'], width, 
                    label='SPS', color='#2ecc71', alpha=0.8)
    bars2 = ax1.bar(x + width/2, comparison_df['random_recall_at_10'], width,
                    label='Random', color='#e74c3c', alpha=0.8)
    
    ax1.set_ylabel('Recall@10')
    ax1.set_title('Recovery of Literature Markers (Top 10 Genes)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(comparison_df['cell_type'], rotation=45, ha='right')
    ax1.legend()
    ax1.set_ylim(0, 1.1)
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        if not np.isnan(height):
            ax1.annotate(f'{height:.2f}',
                        xy=(bar.get_x() + bar.get_width()/2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    
    for bar in bars2:
        height = bar.get_height()
        if not np.isnan(height):
            ax1.annotate(f'{height:.2f}',
                        xy=(bar.get_x() + bar.get_width()/2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    
    # Recall@50
    ax2 = axes[1]
    bars3 = ax2.bar(x - width/2, comparison_df['sps_recall_at_50'], width,
                    label='SPS', color='#2ecc71', alpha=0.8)
    bars4 = ax2.bar(x + width/2, comparison_df['random_recall_at_50'], width,
                    label='Random', color='#e74c3c', alpha=0.8)
    
    ax2.set_ylabel('Recall@50')
    ax2.set_title('Recovery of Literature Markers (Top 50 Genes)')
    ax2.set_xticks(x)
    ax2.set_xticklabels(comparison_df['cell_type'], rotation=45, ha='right')
    ax2.legend()
    ax2.set_ylim(0, 1.1)
    
    for bar in bars3:
        height = bar.get_height()
        if not np.isnan(height):
            ax2.annotate(f'{height:.2f}',
                        xy=(bar.get_x() + bar.get_width()/2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    
    for bar in bars4:
        height = bar.get_height()
        if not np.isnan(height):
            ax2.annotate(f'{height:.2f}',
                        xy=(bar.get_x() + bar.get_width()/2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


def plot_rank_comparison(comparison_df, save_path=None):
    """
    Plot comparison of mean ranks for literature markers.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(len(comparison_df))
    width = 0.35
    
    # Lower rank is better, so we want SPS bars to be shorter
    bars1 = ax.bar(x - width/2, comparison_df['sps_mean_rank'], width,
                   label='SPS', color='#2ecc71', alpha=0.8)
    bars2 = ax.bar(x + width/2, comparison_df['random_mean_rank'], width,
                   label='Random', color='#e74c3c', alpha=0.8)
    
    ax.set_ylabel('Mean Rank of Literature Markers (lower is better)')
    ax.set_title('Ranking of Literature Markers in SEMITONES Results')
    ax.set_xticks(x)
    ax.set_xticklabels(comparison_df['cell_type'], rotation=45, ha='right')
    ax.legend()
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        if height and not np.isnan(height):
            ax.annotate(f'{height:.0f}',
                        xy=(bar.get_x() + bar.get_width()/2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    
    for bar in bars2:
        height = bar.get_height()
        if height and not np.isnan(height):
            ax.annotate(f'{height:.0f}',
                        xy=(bar.get_x() + bar.get_width()/2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


def plot_overlap_heatmap(overlap_results, save_path=None):
    """
    Plot heatmap of marker gene overlap between methods.
    """
    # Create data for heatmap
    celltypes = list(overlap_results.keys())
    metrics = ['n_overlap', 'jaccard_similarity']
    
    data = []
    for ct in celltypes:
        data.append([
            overlap_results[ct]['n_overlap'],
            overlap_results[ct]['jaccard_similarity'] * 100,  # Convert to percentage
        ])
    
    df = pd.DataFrame(data, index=celltypes, columns=['Overlapping Genes', 'Jaccard Similarity (%)'])
    
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(df, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax)
    ax.set_title('Marker Gene Overlap: SPS vs Random (Top 50 Genes)')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


def plot_improvement_summary(comparison_df, save_path=None):
    """
    Plot summary of SPS improvement over Random.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Filter out NaN values
    df = comparison_df.dropna(subset=['recall_10_improvement', 'recall_50_improvement', 'rank_improvement'])
    
    if len(df) == 0:
        print("No data available for improvement summary")
        return
    
    # Recall@10 improvement
    ax1 = axes[0]
    colors = ['#2ecc71' if x > 0 else '#e74c3c' for x in df['recall_10_improvement']]
    ax1.bar(df['cell_type'], df['recall_10_improvement'] * 100, color=colors, alpha=0.8)
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax1.set_ylabel('Improvement (%)')
    ax1.set_title('SPS Improvement in Recall@10')
    ax1.tick_params(axis='x', rotation=45)
    
    # Recall@50 improvement
    ax2 = axes[1]
    colors = ['#2ecc71' if x > 0 else '#e74c3c' for x in df['recall_50_improvement']]
    ax2.bar(df['cell_type'], df['recall_50_improvement'] * 100, color=colors, alpha=0.8)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax2.set_ylabel('Improvement (%)')
    ax2.set_title('SPS Improvement in Recall@50')
    ax2.tick_params(axis='x', rotation=45)
    
    # Rank improvement (positive = SPS better)
    ax3 = axes[2]
    colors = ['#2ecc71' if x > 0 else '#e74c3c' for x in df['rank_improvement']]
    ax3.bar(df['cell_type'], df['rank_improvement'], color=colors, alpha=0.8)
    ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax3.set_ylabel('Rank Improvement (higher = better)')
    ax3.set_title('SPS Improvement in Mean Rank')
    ax3.tick_params(axis='x', rotation=45)
    
    plt.suptitle('SPS vs Random: Literature Marker Recovery Improvement', fontsize=14, y=1.02)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


def plot_literature_marker_recovery(sps_eval, random_eval, save_path=None):
    """
    Plot detailed recovery of individual literature markers.
    """
    celltypes = [ct for ct in sps_eval.keys() if ct in random_eval]
    
    fig, axes = plt.subplots(len(celltypes), 1, figsize=(12, 4 * len(celltypes)))
    
    if len(celltypes) == 1:
        axes = [axes]
    
    for i, celltype in enumerate(celltypes):
        ax = axes[i]
        
        sps_ranks = sps_eval[celltype].get('marker_ranks', {})
        random_ranks = random_eval[celltype].get('marker_ranks', {})
        
        if not sps_ranks or not random_ranks:
            ax.text(0.5, 0.5, f'No data for {celltype}', ha='center', va='center')
            continue
        
        # Get markers present in both
        markers = list(sps_ranks.keys())
        sps_values = [sps_ranks.get(m) for m in markers]
        random_values = [random_ranks.get(m) for m in markers]
        
        x = np.arange(len(markers))
        width = 0.35
        
        ax.bar(x - width/2, sps_values, width, label='SPS', color='#2ecc71', alpha=0.8)
        ax.bar(x + width/2, random_values, width, label='Random', color='#e74c3c', alpha=0.8)
        
        ax.set_ylabel('Rank (lower is better)')
        ax.set_title(f'{celltype}: Ranking of Literature Markers')
        ax.set_xticks(x)
        ax.set_xticklabels(markers, rotation=45, ha='right')
        ax.legend()
        
        # Invert y-axis so lower ranks appear higher
        ax.invert_yaxis()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    plt.close()


# ============================================================================
# Main Comparison Function
# ============================================================================

def run_comparison():
    """
    Run full comparison between SPS and Random methods.
    """
    print("=" * 80)
    print("Comparing SPS vs Random Marker Gene Results")
    print("=" * 80)
    
    # Load results
    print("\nLoading results...")
    sps_results = load_results('sps')
    random_results = load_results('random')
    
    if not sps_results or not random_results:
        print("ERROR: Results not found. Please run semitones_marker_pipeline.py first.")
        print(f"  Expected files in: {RESULTS_DIR}")
        return
    
    if 'evaluation' not in sps_results or 'evaluation' not in random_results:
        print("ERROR: Evaluation results not found.")
        return
    
    sps_eval = sps_results['evaluation']
    random_eval = random_results['evaluation']
    sps_markers = sps_results.get('marker_genes', {})
    random_markers = random_results.get('marker_genes', {})
    
    # Compare recall metrics
    print("\nComparing recall metrics...")
    comparison_df = compare_recall_metrics(sps_eval, random_eval)
    
    # Save comparison table
    comparison_path = os.path.join(RESULTS_DIR, 'marker_gene_comparison.csv')
    comparison_df.to_csv(comparison_path, index=False)
    print(f"Saved comparison table to {comparison_path}")
    
    # Compare marker overlap
    print("\nComparing marker gene overlap...")
    overlap_results = compare_marker_overlap(sps_markers, random_markers, top_n=50)
    
    overlap_path = os.path.join(RESULTS_DIR, 'marker_gene_overlap.pkl')
    with open(overlap_path, 'wb') as f:
        pickle.dump(overlap_results, f)
    print(f"Saved overlap results to {overlap_path}")
    
    # Compute rank correlation
    print("\nComputing rank correlations...")
    correlations = compute_rank_correlation(sps_markers, random_markers)
    
    corr_path = os.path.join(RESULTS_DIR, 'rank_correlations.pkl')
    with open(corr_path, 'wb') as f:
        pickle.dump(correlations, f)
    print(f"Saved correlations to {corr_path}")
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    
    plot_recall_comparison(
        comparison_df,
        save_path=os.path.join(FIGURES_DIR, 'recall_comparison.png')
    )
    
    plot_rank_comparison(
        comparison_df,
        save_path=os.path.join(FIGURES_DIR, 'rank_comparison.png')
    )
    
    plot_overlap_heatmap(
        overlap_results,
        save_path=os.path.join(FIGURES_DIR, 'marker_overlap_heatmap.png')
    )
    
    plot_improvement_summary(
        comparison_df,
        save_path=os.path.join(FIGURES_DIR, 'sps_improvement_summary.png')
    )
    
    plot_literature_marker_recovery(
        sps_eval,
        random_eval,
        save_path=os.path.join(FIGURES_DIR, 'literature_marker_recovery.png')
    )
    
    # Print summary
    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)
    
    print("\nRecall Metrics (SPS vs Random):")
    print(comparison_df[['cell_type', 'sps_recall_at_10', 'random_recall_at_10', 
                         'sps_recall_at_50', 'random_recall_at_50']].to_string(index=False))
    
    print("\nMean Improvement with SPS:")
    print(f"  Recall@10: {comparison_df['recall_10_improvement'].mean():.2%}")
    print(f"  Recall@50: {comparison_df['recall_50_improvement'].mean():.2%}")
    print(f"  Mean Rank: {comparison_df['rank_improvement'].mean():.1f} positions better")
    
    print("\nMarker Gene Overlap (Top 50):")
    for celltype, overlap in overlap_results.items():
        print(f"  {celltype}: {overlap['n_overlap']} overlapping, "
              f"Jaccard={overlap['jaccard_similarity']:.2f}")
    
    print("\nRank Correlations:")
    for celltype, corr in correlations.items():
        if corr['spearman_r'] is not None:
            print(f"  {celltype}: Spearman r={corr['spearman_r']:.3f}, p={corr['p_value']:.2e}")
    
    print("\n" + "=" * 80)
    print("Comparison completed!")
    print(f"Results saved to: {RESULTS_DIR}")
    print(f"Figures saved to: {FIGURES_DIR}")
    print("=" * 80)


# ============================================================================
# Entry Point
# ============================================================================

if __name__ == '__main__':
    run_comparison()
