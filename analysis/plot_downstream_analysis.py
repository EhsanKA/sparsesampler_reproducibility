#!/usr/bin/env python
# coding: utf-8

"""
Generate figures for downstream analysis comparison (Reviewer Response)

This script creates figures showing:
1. Cell type diversity preservation at different sampling depths
2. DE gene overlap between original and subsampled populations
3. Clustering/cell-type annotation comparison results
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from datetime import datetime

# Add analysis directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

# Set up logging
def setup_logger():
    log_dir = os.path.join(os.path.dirname(script_dir), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'plot_downstream_analysis_{timestamp}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

# Constants
METHODS = ['random', 'sps', 'hopper', 'atomic', 'scsampler']
METHOD_COLORS = {
    'random': 'grey',
    'sps': 'purple',
    'hopper': 'blue',
    'atomic': 'orange',
    'scsampler': 'green'
}

DATASET_CONFIGS = {
    'mcc_01': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'x_values': [0.5, 1.0, 2.0, 2.5, 3.2],
        'label': 'MCC_01'
    },
    'mcc_05': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'x_values': [0.5, 1.0, 2.0, 2.5, 3.2],
        'label': 'MCC_05'
    }
}

def get_output_dir():
    """Get output directory for figures."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
    return os.path.join(project_root, 'figures', 'revision', 'downstream_analysis')

def get_results_dir():
    """Get results directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, 'results', 'downstream_analysis')

def plot_diversity_preservation(diversity_df, dataset, output_dir):
    """Plot cell type diversity preservation metrics."""
    logger.info(f"Plotting diversity preservation for {dataset}")
    
    config = DATASET_CONFIGS[dataset]
    x_values = config['x_values']
    references = config['references']
    
    # Filter data for this dataset
    df = diversity_df[diversity_df['dataset'] == dataset].copy()
    
    if len(df) == 0:
        logger.warning(f"No data found for {dataset}")
        return
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'Cell Type Diversity Preservation: {config["label"]}', fontsize=16, fontweight='bold')
    
    # Plot 1: Number of cell types preserved
    ax1 = axes[0, 0]
    for method in METHODS:
        method_data = df[df['method'] == method]
        if len(method_data) == 0:
            continue
        
        # Average across reps and sample sizes for each ref
        y_values = []
        for ref in references:
            ref_data = method_data[method_data['ref'] == ref]
            if len(ref_data) > 0:
                y_values.append(ref_data['n_cell_types_preserved'].mean())
            else:
                y_values.append(np.nan)
        
        ax1.plot(x_values, y_values, '-o', color=METHOD_COLORS[method], 
                label=method, linewidth=2, markersize=6)
    
    ax1.set_xlabel('Reference Size (M)', fontsize=12)
    ax1.set_ylabel('Fraction of Cell Types Preserved', fontsize=12)
    ax1.set_title('Cell Type Count Preservation', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_ylim([0, 1.1])
    
    # Plot 2: Rare cell type coverage
    ax2 = axes[0, 1]
    for method in METHODS:
        method_data = df[df['method'] == method]
        if len(method_data) == 0:
            continue
        
        y_values = []
        for ref in references:
            ref_data = method_data[method_data['ref'] == ref]
            if len(ref_data) > 0:
                y_values.append(ref_data['rare_type_coverage'].mean())
            else:
                y_values.append(np.nan)
        
        ax2.plot(x_values, y_values, '-o', color=METHOD_COLORS[method], 
                label=method, linewidth=2, markersize=6)
    
    ax2.set_xlabel('Reference Size (M)', fontsize=12)
    ax2.set_ylabel('Rare Cell Type Coverage', fontsize=12)
    ax2.set_title('Rare Cell Type Coverage', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_ylim([0, 1.1])
    
    # Plot 3: Rare cell enrichment
    ax3 = axes[1, 0]
    for method in METHODS:
        method_data = df[df['method'] == method]
        if len(method_data) == 0:
            continue
        
        y_values = []
        for ref in references:
            ref_data = method_data[method_data['ref'] == ref]
            if len(ref_data) > 0:
                y_values.append(ref_data['rare_cell_enrichment'].mean())
            else:
                y_values.append(np.nan)
        
        ax3.plot(x_values, y_values, '-o', color=METHOD_COLORS[method], 
                label=method, linewidth=2, markersize=6)
    
    ax3.set_xlabel('Reference Size (M)', fontsize=12)
    ax3.set_ylabel('Rare Cell Enrichment Ratio', fontsize=12)
    ax3.set_title('Rare Cell Enrichment', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=10)
    ax3.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='No enrichment')
    ax3.legend(fontsize=10)
    
    # Plot 4: Cell type overlap
    ax4 = axes[1, 1]
    for method in METHODS:
        method_data = df[df['method'] == method]
        if len(method_data) == 0:
            continue
        
        y_values = []
        for ref in references:
            ref_data = method_data[method_data['ref'] == ref]
            if len(ref_data) > 0:
                y_values.append(ref_data['cell_type_overlap'].mean())
            else:
                y_values.append(np.nan)
        
        ax4.plot(x_values, y_values, '-o', color=METHOD_COLORS[method], 
                label=method, linewidth=2, markersize=6)
    
    ax4.set_xlabel('Reference Size (M)', fontsize=12)
    ax4.set_ylabel('Cell Type Overlap Fraction', fontsize=12)
    ax4.set_title('Cell Type Overlap', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=10)
    ax4.set_ylim([0, 1.1])
    
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, f'{dataset}_diversity_preservation.jpg')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved figure: {output_file}")

def plot_diversity_by_sample_size(diversity_df, dataset, output_dir):
    """Plot diversity preservation across different sample sizes."""
    logger.info(f"Plotting diversity by sample size for {dataset}")
    
    config = DATASET_CONFIGS[dataset]
    sizes = config['sizes']
    
    # Filter data for this dataset
    df = diversity_df[diversity_df['dataset'] == dataset].copy()
    
    if len(df) == 0:
        logger.warning(f"No data found for {dataset}")
        return
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'Diversity Preservation by Sample Size: {config["label"]}', 
                 fontsize=16, fontweight='bold')
    
    metrics = [
        ('n_cell_types_preserved', 'Fraction of Cell Types Preserved', axes[0, 0]),
        ('rare_type_coverage', 'Rare Cell Type Coverage', axes[0, 1]),
        ('rare_cell_enrichment', 'Rare Cell Enrichment Ratio', axes[1, 0]),
        ('weighted_preservation', 'Weighted Preservation Score', axes[1, 1])
    ]
    
    for metric, ylabel, ax in metrics:
        # Group by method and sample size
        for method in METHODS:
            method_data = df[df['method'] == method]
            if len(method_data) == 0:
                continue
            
            y_values = []
            for size in sizes:
                size_data = method_data[method_data['sample_size'] == size]
                if len(size_data) > 0:
                    y_values.append(size_data[metric].mean())
                else:
                    y_values.append(np.nan)
            
            ax.plot(sizes, y_values, '-o', color=METHOD_COLORS[method], 
                   label=method, linewidth=2, markersize=6)
        
        ax.set_xlabel('Sample Size', fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(ylabel.replace(' Preserved', '').replace(' Fraction', ''), 
                    fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10)
        ax.set_ylim([0, 1.1])
    
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, f'{dataset}_diversity_by_sample_size.jpg')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved figure: {output_file}")

def plot_downstream_comparison(downstream_df, dataset, output_dir):
    """Plot downstream analysis comparison results (DE and clustering)."""
    logger.info(f"Plotting downstream comparison for {dataset}")
    
    # Filter data for this dataset
    df = downstream_df[downstream_df['dataset'] == dataset].copy()
    
    if len(df) == 0:
        logger.warning(f"No data found for {dataset}")
        return
    
    config = DATASET_CONFIGS[dataset]
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'Downstream Analysis Comparison: {config["label"]}', 
                 fontsize=16, fontweight='bold')
    
    # Plot 1: DE gene overlap
    ax1 = axes[0, 0]
    if 'avg_de_gene_overlap' in df.columns and df['avg_de_gene_overlap'].notna().any():
        de_data = df[df['avg_de_gene_overlap'].notna()]
        if len(de_data) > 0:
            sns.boxplot(data=de_data, x='method', y='avg_de_gene_overlap', ax=ax1, 
                       palette=METHOD_COLORS)
            ax1.set_xlabel('Method', fontsize=12)
            ax1.set_ylabel('Average DE Gene Overlap', fontsize=12)
            ax1.set_title('Differential Expression Gene Overlap', fontsize=13, fontweight='bold')
            ax1.grid(True, alpha=0.3, axis='y')
            ax1.tick_params(axis='x', rotation=45)
    
    # Plot 2: Cell type frequency correlation
    ax2 = axes[0, 1]
    if 'frequency_correlation' in df.columns and df['frequency_correlation'].notna().any():
        freq_data = df[df['frequency_correlation'].notna()]
        if len(freq_data) > 0:
            sns.boxplot(data=freq_data, x='method', y='frequency_correlation', ax=ax2, 
                       palette=METHOD_COLORS)
            ax2.set_xlabel('Method', fontsize=12)
            ax2.set_ylabel('Frequency Correlation', fontsize=12)
            ax2.set_title('Cell Type Frequency Correlation', fontsize=13, fontweight='bold')
            ax2.grid(True, alpha=0.3, axis='y')
            ax2.tick_params(axis='x', rotation=45)
    
    # Plot 3: Jaccard similarity
    ax3 = axes[1, 0]
    if 'jaccard_similarity' in df.columns and df['jaccard_similarity'].notna().any():
        jaccard_data = df[df['jaccard_similarity'].notna()]
        if len(jaccard_data) > 0:
            sns.boxplot(data=jaccard_data, x='method', y='jaccard_similarity', ax=ax3, 
                       palette=METHOD_COLORS)
            ax3.set_xlabel('Method', fontsize=12)
            ax3.set_ylabel('Jaccard Similarity', fontsize=12)
            ax3.set_title('Cell Type Jaccard Similarity', fontsize=13, fontweight='bold')
            ax3.grid(True, alpha=0.3, axis='y')
            ax3.tick_params(axis='x', rotation=45)
    
    # Plot 4: Cell type overlap
    ax4 = axes[1, 1]
    if 'cell_type_overlap' in df.columns and df['cell_type_overlap'].notna().any():
        overlap_data = df[df['cell_type_overlap'].notna()]
        if len(overlap_data) > 0:
            sns.boxplot(data=overlap_data, x='method', y='cell_type_overlap', ax=ax4, 
                       palette=METHOD_COLORS)
            ax4.set_xlabel('Method', fontsize=12)
            ax4.set_ylabel('Cell Type Overlap Fraction', fontsize=12)
            ax4.set_title('Cell Type Annotation Overlap', fontsize=13, fontweight='bold')
            ax4.grid(True, alpha=0.3, axis='y')
            ax4.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, f'{dataset}_downstream_comparison.jpg')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved figure: {output_file}")

def create_combined_figure(diversity_df, downstream_df, output_dir):
    """Create a combined figure showing key results for the reviewer."""
    logger.info("Creating combined figure for reviewer")
    
    # Filter for mcc_01 dataset
    diversity_mcc01 = diversity_df[diversity_df['dataset'] == 'mcc_01'].copy()
    downstream_mcc01 = downstream_df[downstream_df['dataset'] == 'mcc_01'].copy()
    
    config = DATASET_CONFIGS['mcc_01']
    x_values = config['x_values']
    references = config['references']
    
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Title
    fig.suptitle('Downstream Analysis Comparison: Cell Type Diversity and Analysis Preservation', 
                 fontsize=18, fontweight='bold', y=0.98)
    
    # Plot 1: Cell type count preservation (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    for method in ['sps', 'random']:  # Focus on key methods
        method_data = diversity_mcc01[diversity_mcc01['method'] == method]
        if len(method_data) > 0:
            y_values = []
            for ref in references:
                ref_data = method_data[method_data['ref'] == ref]
                if len(ref_data) > 0:
                    y_values.append(ref_data['n_cell_types_preserved'].mean())
                else:
                    y_values.append(np.nan)
            ax1.plot(x_values, y_values, '-o', color=METHOD_COLORS[method], 
                    label=method, linewidth=2.5, markersize=8)
    ax1.set_xlabel('Reference Size (M)', fontsize=11)
    ax1.set_ylabel('Fraction Preserved', fontsize=11)
    ax1.set_title('Cell Type Count Preservation', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_ylim([0, 1.1])
    
    # Plot 2: Shannon diversity preservation (top middle)
    ax2 = fig.add_subplot(gs[0, 1])
    for method in ['sps', 'random']:
        method_data = diversity_mcc01[diversity_mcc01['method'] == method]
        if len(method_data) > 0:
            y_values = []
            for ref in references:
                ref_data = method_data[method_data['ref'] == ref]
                if len(ref_data) > 0:
                    y_values.append(ref_data['shannon_preserved'].mean())
                else:
                    y_values.append(np.nan)
            ax2.plot(x_values, y_values, '-o', color=METHOD_COLORS[method], 
                    label=method, linewidth=2.5, markersize=8)
    ax2.set_xlabel('Reference Size (M)', fontsize=11)
    ax2.set_ylabel('Shannon Diversity Preserved', fontsize=11)
    ax2.set_title('Shannon Diversity Preservation', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_ylim([0, 1.1])
    
    # Plot 3: DE gene overlap (top right)
    ax3 = fig.add_subplot(gs[0, 2])
    if len(downstream_mcc01) > 0 and 'avg_de_gene_overlap' in downstream_mcc01.columns:
        de_data = downstream_mcc01[downstream_mcc01['avg_de_gene_overlap'].notna()]
        if len(de_data) > 0:
            methods_to_plot = ['sps', 'random']
            plot_data = de_data[de_data['method'].isin(methods_to_plot)]
            if len(plot_data) > 0:
                sns.boxplot(data=plot_data, x='method', y='avg_de_gene_overlap', ax=ax3, 
                           palette=METHOD_COLORS)
                ax3.set_xlabel('Method', fontsize=11)
                ax3.set_ylabel('DE Gene Overlap', fontsize=11)
                ax3.set_title('Differential Expression Gene Overlap', fontsize=12, fontweight='bold')
                ax3.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Cell type frequency correlation (middle left)
    ax4 = fig.add_subplot(gs[1, 0])
    if len(downstream_mcc01) > 0 and 'frequency_correlation' in downstream_mcc01.columns:
        freq_data = downstream_mcc01[downstream_mcc01['frequency_correlation'].notna()]
        if len(freq_data) > 0:
            methods_to_plot = ['sps', 'random']
            plot_data = freq_data[freq_data['method'].isin(methods_to_plot)]
            if len(plot_data) > 0:
                sns.boxplot(data=plot_data, x='method', y='frequency_correlation', ax=ax4, 
                           palette=METHOD_COLORS)
                ax4.set_xlabel('Method', fontsize=11)
                ax4.set_ylabel('Frequency Correlation', fontsize=11)
                ax4.set_title('Cell Type Frequency Correlation', fontsize=12, fontweight='bold')
                ax4.grid(True, alpha=0.3, axis='y')
    
    # Plot 5: Jaccard similarity (middle middle)
    ax5 = fig.add_subplot(gs[1, 1])
    if len(downstream_mcc01) > 0 and 'jaccard_similarity' in downstream_mcc01.columns:
        jaccard_data = downstream_mcc01[downstream_mcc01['jaccard_similarity'].notna()]
        if len(jaccard_data) > 0:
            methods_to_plot = ['sps', 'random']
            plot_data = jaccard_data[jaccard_data['method'].isin(methods_to_plot)]
            if len(plot_data) > 0:
                sns.boxplot(data=plot_data, x='method', y='jaccard_similarity', ax=ax5, 
                           palette=METHOD_COLORS)
                ax5.set_xlabel('Method', fontsize=11)
                ax5.set_ylabel('Jaccard Similarity', fontsize=11)
                ax5.set_title('Cell Type Jaccard Similarity', fontsize=12, fontweight='bold')
                ax5.grid(True, alpha=0.3, axis='y')
    
    # Plot 6: Rare cell type coverage by sample size (middle right)
    ax6 = fig.add_subplot(gs[1, 2])
    sizes = config['sizes']
    for method in ['sps', 'random']:
        method_data = diversity_mcc01[diversity_mcc01['method'] == method]
        if len(method_data) > 0:
            y_values = []
            for size in sizes:
                size_data = method_data[method_data['sample_size'] == size]
                if len(size_data) > 0:
                    y_values.append(size_data['rare_type_coverage'].mean())
                else:
                    y_values.append(np.nan)
            ax6.plot(sizes, y_values, '-o', color=METHOD_COLORS[method], 
                    label=method, linewidth=2.5, markersize=8)
    ax6.set_xlabel('Sample Size', fontsize=11)
    ax6.set_ylabel('Rare Cell Type Coverage', fontsize=11)
    ax6.set_title('Rare Cell Coverage by Sample Size', fontsize=12, fontweight='bold')
    ax6.grid(True, alpha=0.3)
    ax6.legend(fontsize=10)
    ax6.set_ylim([0, 1.1])
    
    # Plot 7-9: Summary statistics (bottom row)
    # Could add summary tables or additional metrics here
    
    plt.savefig(os.path.join(output_dir, 'reviewer_response_combined_figure.jpg'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved combined figure: {os.path.join(output_dir, 'reviewer_response_combined_figure.jpg')}")

def main():
    """Main function to generate all figures."""
    logger.info("Starting figure generation for downstream analysis")
    
    output_dir = get_output_dir()
    os.makedirs(output_dir, exist_ok=True)
    
    results_dir = get_results_dir()
    
    # Load diversity results
    diversity_file = os.path.join(results_dir, 'combined_diversity_preservation.csv')
    if os.path.exists(diversity_file):
        diversity_df = pd.read_csv(diversity_file)
        logger.info(f"Loaded diversity results: {len(diversity_df)} rows")
        
        # Plot for each dataset
        for dataset in ['mcc_01', 'mcc_05']:
            plot_diversity_preservation(diversity_df, dataset, output_dir)
            plot_diversity_by_sample_size(diversity_df, dataset, output_dir)
    else:
        logger.warning(f"Diversity results file not found: {diversity_file}")
        diversity_df = pd.DataFrame()
    
    # Load downstream comparison results
    downstream_file = os.path.join(results_dir, 'combined_downstream_comparison.csv')
    if os.path.exists(downstream_file):
        downstream_df = pd.read_csv(downstream_file)
        logger.info(f"Loaded downstream results: {len(downstream_df)} rows")
        
        # Plot for each dataset
        for dataset in ['mcc_01']:
            plot_downstream_comparison(downstream_df, dataset, output_dir)
    else:
        logger.warning(f"Downstream results file not found: {downstream_file}")
        downstream_df = pd.DataFrame()
    
    # Create combined figure
    if len(diversity_df) > 0 and len(downstream_df) > 0:
        create_combined_figure(diversity_df, downstream_df, output_dir)
    
    logger.info("All figures generated successfully")

if __name__ == "__main__":
    main()

