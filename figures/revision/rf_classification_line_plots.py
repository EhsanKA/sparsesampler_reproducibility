#!/usr/bin/env python
# coding: utf-8

"""
RF Classification Line Plots by Feature Index

Creates clean line plots showing:
- X-axis: Feature Index (1-30)
- Y-axis: Percentage (Delta F1)
- Lines: Delta Macro F1, Delta Rare Cell F1

For datasets: MCC, MCC_01, MCC_05, LCMV
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
output_dir = script_dir

# Data paths
DATA_PATHS = {
    'mcc': os.path.join(project_root, "jobs/feature_index_classification/results/all_feature_indices_summary.csv"),
    'mcc_01': os.path.join(project_root, "jobs/feature_index_classification/results_mcc_01/all_feature_indices_summary.csv"),
    'mcc_05': os.path.join(project_root, "jobs/feature_index_classification/results_mcc_05/all_feature_indices_summary.csv"),
    'lcmv': os.path.join(project_root, "jobs/feature_index_classification_lcmv/results/all_feature_indices_summary.csv"),
}

# Dataset configurations
DATASET_CONFIGS = {
    'mcc': {
        'label': 'MCC (3.2M cells, 0.95% rare)',
        'rare_col': 'osteoblast_f1_improvement',
        'rare_label': 'Δ Osteoblast F1',
    },
    'mcc_01': {
        'label': 'MCC_01 (3.1M cells, 0.1% rare)',
        'rare_col': 'osteoblast_f1_improvement',
        'rare_label': 'Δ Osteoblast F1',
    },
    'mcc_05': {
        'label': 'MCC_05 (3.1M cells, 0.5% rare)',
        'rare_col': 'osteoblast_f1_improvement',
        'rare_label': 'Δ Osteoblast F1',
    },
    'lcmv': {
        'label': 'LCMV (34M cells)',
        'rare_col': 'interacting_f1_improvement',
        'rare_label': 'Δ Interacting F1 (0.07%)',
    },
}

# Plot styling
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'legend.fontsize': 11,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'lines.linewidth': 2.5,
    'lines.markersize': 7,
    'figure.dpi': 100,
    'savefig.dpi': 350,
})


def load_all_data():
    """Load all dataset results."""
    data = {}
    for dataset, path in DATA_PATHS.items():
        if os.path.exists(path):
            df = pd.read_csv(path)
            df = df.sort_values('feature_index')
            data[dataset] = df
            print(f"Loaded {dataset}: {len(df)} feature indices")
        else:
            print(f"Warning: {dataset} not found at {path}")
    return data


def create_line_plots(data):
    """Create a 2x2 grid of line plots for all datasets."""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    # Colors
    color_macro = '#3498db'  # Blue
    color_rare = '#e74c3c'   # Red
    
    datasets_order = ['mcc', 'mcc_01', 'mcc_05', 'lcmv']
    
    for idx, dataset in enumerate(datasets_order):
        ax = axes[idx]
        
        if dataset not in data:
            ax.text(0.5, 0.5, f'{dataset} data not available', 
                   ha='center', va='center', transform=ax.transAxes)
            continue
        
        df = data[dataset]
        config = DATASET_CONFIGS[dataset]
        
        x = df['feature_index'].values
        
        # Plot Delta Macro F1
        ax.plot(x, df['macro_f1_improvement'].values, '-o', 
                color=color_macro, label='Δ Macro F1', alpha=0.9)
        
        # Plot Delta Rare Cell F1
        rare_col = config['rare_col']
        if rare_col in df.columns:
            ax.plot(x, df[rare_col].values, '-s', 
                    color=color_rare, label=config['rare_label'], alpha=0.9)
            
            # Shade positive region for rare cell type
            ax.fill_between(x, 0, df[rare_col].values,
                           where=df[rare_col].values > 0,
                           alpha=0.15, color=color_rare)
        
        # Zero line
        ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.6)
        
        # Labels and title
        ax.set_xlabel('Feature Index')
        ax.set_ylabel('ΔF1 (SPS - Random) %')
        ax.set_title(config['label'], fontweight='bold')
        
        # Legend
        ax.legend(loc='best')
        
        # Grid
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # X-axis limits
        ax.set_xlim(0, max(x) + 1)
        
        # Set x-ticks at intervals
        ax.set_xticks(np.arange(0, max(x) + 1, 5))
        
        # Annotate best rare cell improvement
        if rare_col in df.columns:
            best_idx = df[rare_col].idxmax()
            best_fi = df.loc[best_idx, 'feature_index']
            best_val = df.loc[best_idx, rare_col]
            if best_val > 0:
                ax.annotate(f'+{best_val:.1f}%', 
                           xy=(best_fi, best_val),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=10, color='darkred', fontweight='bold')
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(output_dir, 'rf_classification_line_plots.jpg')
    plt.savefig(output_path, bbox_inches='tight', dpi=350)
    print(f"\nSaved: {output_path}")
    
    pdf_path = os.path.join(output_dir, 'rf_classification_line_plots.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {pdf_path}")
    
    plt.close()


def create_combined_line_plot(data):
    """Create a single plot with all datasets showing rare cell type improvement."""
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Colors for each dataset
    colors = {
        'mcc': '#2ecc71',      # Green
        'mcc_01': '#e74c3c',   # Red
        'mcc_05': '#3498db',   # Blue
        'lcmv': '#9b59b6',     # Purple
    }
    
    markers = {
        'mcc': 'o',
        'mcc_01': 's',
        'mcc_05': '^',
        'lcmv': 'D',
    }
    
    labels = {
        'mcc': 'MCC (0.95% rare)',
        'mcc_01': 'MCC_01 (0.1% rare)',
        'mcc_05': 'MCC_05 (0.5% rare)',
        'lcmv': 'LCMV (0.07% rare)',
    }
    
    for dataset in ['mcc', 'mcc_01', 'mcc_05', 'lcmv']:
        if dataset not in data:
            continue
        
        df = data[dataset]
        config = DATASET_CONFIGS[dataset]
        rare_col = config['rare_col']
        
        if rare_col not in df.columns:
            continue
        
        x = df['feature_index'].values
        y = df[rare_col].values
        
        ax.plot(x, y, f'-{markers[dataset]}', 
                color=colors[dataset], 
                label=labels[dataset],
                alpha=0.85, linewidth=2.5, markersize=7)
    
    # Zero line
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.6)
    
    # Labels
    ax.set_xlabel('Feature Index', fontsize=14)
    ax.set_ylabel('Δ Rare Cell Type F1 (SPS - Random) %', fontsize=14)
    ax.set_title('SPS Improvement for Rare Cell Types by Feature Index', fontsize=16, fontweight='bold')
    
    # Legend
    ax.legend(loc='best', fontsize=11)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # X-axis
    ax.set_xlim(0, 31)
    ax.set_xticks(np.arange(0, 31, 5))
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(output_dir, 'rf_classification_rare_cell_comparison.jpg')
    plt.savefig(output_path, bbox_inches='tight', dpi=350)
    print(f"Saved: {output_path}")
    
    pdf_path = os.path.join(output_dir, 'rf_classification_rare_cell_comparison.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {pdf_path}")
    
    plt.close()


def create_macro_f1_comparison(data):
    """Create a plot comparing macro F1 improvement across datasets."""
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    colors = {
        'mcc': '#2ecc71',
        'mcc_01': '#e74c3c',
        'mcc_05': '#3498db',
        'lcmv': '#9b59b6',
    }
    
    markers = {
        'mcc': 'o',
        'mcc_01': 's',
        'mcc_05': '^',
        'lcmv': 'D',
    }
    
    labels = {
        'mcc': 'MCC (3.2M)',
        'mcc_01': 'MCC_01 (3.1M)',
        'mcc_05': 'MCC_05 (3.1M)',
        'lcmv': 'LCMV (34M)',
    }
    
    for dataset in ['mcc', 'mcc_01', 'mcc_05', 'lcmv']:
        if dataset not in data:
            continue
        
        df = data[dataset]
        x = df['feature_index'].values
        y = df['macro_f1_improvement'].values
        
        ax.plot(x, y, f'-{markers[dataset]}', 
                color=colors[dataset], 
                label=labels[dataset],
                alpha=0.85, linewidth=2.5, markersize=7)
    
    # Zero line
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.6)
    
    # Labels
    ax.set_xlabel('Feature Index', fontsize=14)
    ax.set_ylabel('Δ Macro F1 (SPS - Random) %', fontsize=14)
    ax.set_title('Macro F1 Improvement by Feature Index', fontsize=16, fontweight='bold')
    
    # Legend
    ax.legend(loc='best', fontsize=11)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # X-axis
    ax.set_xlim(0, 31)
    ax.set_xticks(np.arange(0, 31, 5))
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(output_dir, 'rf_classification_macro_f1_comparison.jpg')
    plt.savefig(output_path, bbox_inches='tight', dpi=350)
    print(f"Saved: {output_path}")
    
    plt.close()


if __name__ == "__main__":
    # Load data
    data = load_all_data()
    
    # Create plots
    print("\nGenerating line plots...")
    create_line_plots(data)
    create_combined_line_plot(data)
    create_macro_f1_comparison(data)
    
    print("\nDone!")
