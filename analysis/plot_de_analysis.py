#!/usr/bin/env python
# coding: utf-8

"""
Visualize DE Gene Analysis Results
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from datetime import datetime

# Set up logging
def setup_logger():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_dir = os.path.join(os.path.dirname(script_dir), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'plot_de_analysis_{timestamp}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

# Constants
METHOD_COLORS = {
    'random': 'grey',
    'sps': 'purple',
    'hopper': 'blue',
    'atomic': 'orange',
    'scsampler': 'green'
}

def get_results_dir():
    """Get results directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, 'results', 'de_analysis')

def get_figures_dir():
    """Get figures directory."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
    return os.path.join(project_root, 'figures', 'revision', 'de_analysis')

def plot_de_overlap_summary(summary_df, output_dir):
    """Plot summary of DE gene overlap."""
    logger.info("Creating DE overlap summary plots")
    
    # Set style
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'legend.fontsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.figsize': [10, 6],
        'savefig.dpi': 300
    })
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: Boxplot by method
    ax1 = axes[0]
    methods_to_plot = ['sps', 'random']
    plot_data = summary_df[summary_df['method'].isin(methods_to_plot)]
    
    if len(plot_data) > 0:
        sns.boxplot(data=plot_data, x='method', y='mean_overlap', ax=ax1, 
                   palette=METHOD_COLORS, order=methods_to_plot)
        ax1.set_xlabel('Method', fontsize=12)
        ax1.set_ylabel('Mean DE Gene Overlap', fontsize=12)
        ax1.set_title('DE Gene Overlap by Method', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')
        ax1.set_ylim([0, 1.1])
    
    # Plot 2: Overlap by sample size
    ax2 = axes[1]
    for method in methods_to_plot:
        method_data = summary_df[summary_df['method'] == method]
        if len(method_data) > 0:
            sizes = sorted(method_data['sample_size'].unique())
            overlaps = [method_data[method_data['sample_size'] == s]['mean_overlap'].mean() 
                       for s in sizes]
            ax2.plot(sizes, overlaps, '-o', color=METHOD_COLORS[method], 
                    label=method, linewidth=2, markersize=8)
    
    ax2.set_xlabel('Sample Size', fontsize=12)
    ax2.set_ylabel('Mean DE Gene Overlap', fontsize=12)
    ax2.set_title('DE Gene Overlap by Sample Size', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_ylim([0, 1.1])
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, 'de_overlap_summary.jpg')
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    logger.info(f"Saved summary plot: {output_file}")

def plot_de_overlap_by_cell_type(detailed_df, output_dir, dataset, method, size):
    """Plot DE overlap for individual cell types."""
    logger.info(f"Creating cell-type-specific plot for {dataset}, {method}, size {size}")
    
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 10,
        'axes.titlesize': 12,
        'figure.figsize': [12, 8],
        'savefig.dpi': 300
    })
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Sort by overlap
    plot_df = detailed_df.sort_values('overlap_fraction', ascending=True)
    
    # Create bar plot
    bars = ax.barh(range(len(plot_df)), plot_df['overlap_fraction'], 
                   color=METHOD_COLORS.get(method, 'blue'), alpha=0.7)
    
    # Add value labels
    for i, (idx, row) in enumerate(plot_df.iterrows()):
        ax.text(row['overlap_fraction'] + 0.01, i, 
               f"{row['overlap_fraction']:.2%} ({row['n_overlap']}/{row['n_full_genes']})",
               va='center', fontsize=9)
    
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df['cell_type'], fontsize=9)
    ax.set_xlabel('DE Gene Overlap Fraction', fontsize=11)
    ax.set_title(f'DE Gene Overlap by Cell Type\n{dataset.upper()}, {method.upper()}, Size={size:,}', 
                fontsize=13, fontweight='bold')
    ax.set_xlim([0, 1.1])
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, f'{dataset}_{method}_size{size}_by_celltype.jpg')
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    logger.info(f"Saved cell-type plot: {output_file}")

def plot_de_comparison_combined(summary_df, output_dir):
    """Create combined comparison figure."""
    logger.info("Creating combined DE comparison figure")
    
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 11,
        'axes.titlesize': 13,
        'legend.fontsize': 10,
        'figure.figsize': [18, 10],
        'savefig.dpi': 300
    })
    
    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    fig.suptitle('Differential Expression Gene Analysis: Full vs Subsampled Data', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # Plot 1: Overlap by method (all datasets)
    ax1 = axes[0, 0]
    methods_to_plot = ['sps', 'random']
    plot_data = summary_df[summary_df['method'].isin(methods_to_plot)]
    if len(plot_data) > 0:
        sns.boxplot(data=plot_data, x='method', y='mean_overlap', ax=ax1, 
                   palette=METHOD_COLORS, order=methods_to_plot)
        ax1.set_xlabel('Method', fontsize=11)
        ax1.set_ylabel('Mean DE Gene Overlap', fontsize=11)
        ax1.set_title('DE Gene Overlap by Method', fontsize=13, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')
        ax1.set_ylim([0, 1.1])
    
    # Plot 2: Overlap by dataset
    ax2 = axes[0, 1]
    if 'dataset' in summary_df.columns:
        datasets = summary_df['dataset'].unique()
        plot_data = summary_df[summary_df['method'].isin(methods_to_plot)]
        if len(plot_data) > 0:
            sns.boxplot(data=plot_data, x='dataset', y='mean_overlap', hue='method', 
                       ax=ax2, palette=METHOD_COLORS)
            ax2.set_xlabel('Dataset', fontsize=11)
            ax2.set_ylabel('Mean DE Gene Overlap', fontsize=11)
            ax2.set_title('DE Gene Overlap by Dataset', fontsize=13, fontweight='bold')
            ax2.grid(True, alpha=0.3, axis='y')
            ax2.set_ylim([0, 1.1])
            ax2.legend(title='Method', fontsize=9)
    
    # Plot 3: Overlap by sample size
    ax3 = axes[1, 0]
    for method in methods_to_plot:
        method_data = summary_df[summary_df['method'] == method]
        if len(method_data) > 0:
            sizes = sorted(method_data['sample_size'].unique())
            overlaps = [method_data[method_data['sample_size'] == s]['mean_overlap'].mean() 
                        for s in sizes]
            ax3.plot(sizes, overlaps, '-o', color=METHOD_COLORS[method], 
                    label=method, linewidth=2.5, markersize=10)
    ax3.set_xlabel('Sample Size', fontsize=11)
    ax3.set_ylabel('Mean DE Gene Overlap', fontsize=11)
    ax3.set_title('DE Gene Overlap vs Sample Size', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=10)
    ax3.set_ylim([0, 1.1])
    
    # Plot 4: Number of cell types analyzed
    ax4 = axes[1, 1]
    plot_data = summary_df[summary_df['method'].isin(methods_to_plot)]
    if len(plot_data) > 0:
        sns.boxplot(data=plot_data, x='method', y='n_cell_types_analyzed', ax=ax4, 
                   palette=METHOD_COLORS, order=methods_to_plot)
        ax4.set_xlabel('Method', fontsize=11)
        ax4.set_ylabel('Number of Cell Types Analyzed', fontsize=11)
        ax4.set_title('Cell Types Analyzed by Method', fontsize=13, fontweight='bold')
        ax4.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    output_file = os.path.join(output_dir, 'de_analysis_combined.jpg')
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    logger.info(f"Saved combined plot: {output_file}")

def main():
    """Main function to generate all DE analysis plots."""
    logger.info("Starting DE analysis visualization")
    
    results_dir = get_results_dir()
    figures_dir = get_figures_dir()
    os.makedirs(figures_dir, exist_ok=True)
    
    # Load summary results
    summary_file = os.path.join(results_dir, 'de_analysis_summary.csv')
    if not os.path.exists(summary_file):
        logger.error(f"Summary file not found: {summary_file}")
        logger.error("Please run run_de_analysis.py first to generate results")
        return
    
    summary_df = pd.read_csv(summary_file)
    logger.info(f"Loaded summary results: {len(summary_df)} rows")
    
    # Create summary plots
    plot_de_overlap_summary(summary_df, figures_dir)
    plot_de_comparison_combined(summary_df, figures_dir)
    
    # Create cell-type-specific plots for each configuration
    for _, row in summary_df.iterrows():
        dataset = row['dataset']
        method = row['method']
        size = row['sample_size']
        rep = row['rep']
        
        detailed_file = os.path.join(
            results_dir, 
            f'{dataset}_{method}_size{size}_rep{rep}_de_detailed.csv'
        )
        
        if os.path.exists(detailed_file):
            detailed_df = pd.read_csv(detailed_file)
            plot_de_overlap_by_cell_type(detailed_df, figures_dir, dataset, method, size)
        else:
            logger.warning(f"Detailed file not found: {detailed_file}")
    
    logger.info("All DE analysis plots generated successfully!")
    logger.info(f"Figures saved to: {figures_dir}")

if __name__ == "__main__":
    main()






