#!/usr/bin/env python
"""
Plot EVR Index Sensitivity from pre-computed tables with coverage percentage.

This script creates a 2-row figure:
- Top row: Coverage percentage (rare cells captured / total rare cells in reference)
- Bottom row: Raw number of rare cells captured (for reference)

The coverage percentage is defined as:
  coverage = (rare cells captured at given EVR) / (total rare cells in reference size)
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Paths
script_dir = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.join(script_dir, 'tables')
project_root = os.path.dirname(os.path.dirname(script_dir))
FIGURE_DIR = os.path.join(project_root, 'figures', 'revision')
DATA_PATH = os.path.join(project_root, 'data')
os.makedirs(FIGURE_DIR, exist_ok=True)

# Add analysis directory to path for rare cell type identification
analysis_dir = os.path.join(project_root, 'analysis')
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

# Dataset configurations
DATASETS = {
    'mcc': {
        'label': 'MCC (1%)',
        'refs': ['0.5M', '1M', '2M', '2.5M', '3.2M'],
        'ref_nums': [5, 10, 20, 25, 30],
    },
    'mcc_05': {
        'label': 'MCC (0.5%)',
        'refs': ['0.5M', '1M', '2M', '2.5M', '3.2M'],
        'ref_nums': [5, 10, 20, 25, 30],
    },
    'mcc_01': {
        'label': 'MCC (0.1%)',
        'refs': ['0.5M', '1M', '2M', '2.5M', '3.2M'],
        'ref_nums': [5, 10, 20, 25, 30],
    },
    'lcmv': {
        'label': 'LCMV',
        'refs': ['1M', '5M', '10M', '20M', '34M'],
        'ref_nums': [1, 5, 10, 20, 34],
    },
}

# EVR range 1-30
EVR_RANGE = list(range(1, 31))

# Define rare cell types for each dataset
# Based on analysis - osteoblast is the rarest cell type in MCC (~1% frequency)
# LCMV rare cell types based on refined_rare_cell_type_definition.py
RARE_CELL_TYPES = {
    'lcmv': ['NK1_1_TCRgd_T', 'interacting', 'cDC2', 'pDCs', 'CD4_LCMV_spec'],
    'mcc': ['osteoblast'],
    'mcc_01': ['osteoblast'],
    'mcc_05': ['osteoblast'],
}


def count_rare_cells_in_reference(obs_path, rare_types, label_col='celltype'):
    """Count rare cells in the reference dataset from obs.csv."""
    if not os.path.exists(obs_path):
        return None
    
    obs_df = pd.read_csv(obs_path)
    
    # Check for celltype column (may be unnamed index column)
    if label_col not in obs_df.columns:
        # Try common alternatives
        for col in ['cell_type', 'CellType', 'cellType', 'Celltype']:
            if col in obs_df.columns:
                label_col = col
                break
        else:
            # Check if first unnamed column looks like cell types
            if obs_df.columns[0].startswith('Unnamed'):
                label_col = obs_df.columns[0]
            else:
                print(f"Warning: Could not find celltype column in {obs_path}")
                print(f"Available columns: {obs_df.columns.tolist()}")
                return None
    
    total_rare = 0
    for rare_type in rare_types:
        count = (obs_df[label_col] == rare_type).sum()
        total_rare += count
    
    return total_rare


def get_total_rare_cells_per_reference(dataset, cfg):
    """Get total rare cell counts for each reference size."""
    rare_types = RARE_CELL_TYPES.get(dataset, [])
    total_rare_per_ref = {}
    
    for ref_name, ref_num in zip(cfg['refs'], cfg['ref_nums']):
        obs_path = os.path.join(DATA_PATH, dataset, 'benchmark', str(ref_num), 'obs.csv')
        total_rare = count_rare_cells_in_reference(obs_path, rare_types)
        total_rare_per_ref[ref_name] = total_rare
        
        if total_rare is not None:
            print(f"  {dataset} ref {ref_name}: {total_rare} rare cells")
    
    return total_rare_per_ref


def plot_evr_coverage_for_sample_size(sample_size, total_rare_cells):
    """Generate EVR sensitivity coverage plot for a specific sample size.
    
    Args:
        sample_size: Sample size (e.g., 50000, 100000, 200000)
        total_rare_cells: Pre-computed dict of total rare cells per dataset/reference
    """
    # Plot settings
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'legend.fontsize': 9,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.figsize': [16, 12],
        'lines.markersize': 6,
        'lines.linewidth': 1.8,
        'savefig.dpi': 300,
    })
    
    # Create figure with 2 rows: top = coverage %, bottom = raw counts
    fig, axes = plt.subplots(2, 4, figsize=(18, 10))
    
    # Plot each dataset
    for idx, (dataset, cfg) in enumerate(DATASETS.items()):
        ax_coverage = axes[0, idx]  # Top row: coverage percentage
        ax_raw = axes[1, idx]       # Bottom row: raw counts
        
        # Load the table for specified sample size
        table_path = os.path.join(TABLES_DIR, f'{dataset}_feature_index_table_size_{sample_size}.csv')
        if not os.path.exists(table_path):
            ax_coverage.text(0.5, 0.5, 'No data', ha='center', va='center', 
                           transform=ax_coverage.transAxes)
            ax_coverage.set_title(cfg['label'])
            ax_raw.text(0.5, 0.5, 'No data', ha='center', va='center', 
                       transform=ax_raw.transAxes)
            continue
        
        df = pd.read_csv(table_path)
        df = df[df['feature_index'].isin(EVR_RANGE)]
        
        colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(cfg['refs'])))
        
        for j, ref in enumerate(cfg['refs']):
            if ref not in df.columns:
                continue
            
            raw_values = df[ref].values
            x_values = df['feature_index'].values
            
            # Calculate coverage percentage
            total_rare = total_rare_cells[dataset].get(ref)
            if total_rare is not None and total_rare > 0:
                coverage_pct = (raw_values / total_rare) * 100
            else:
                coverage_pct = np.zeros_like(raw_values)
            
            # Plot coverage percentage (top row)
            ax_coverage.plot(x_values, coverage_pct, '-o', color=colors[j],
                           label=f'{ref}', linewidth=1.8, markersize=5, alpha=0.85)
            
            # Plot raw counts (bottom row)
            ax_raw.plot(x_values, raw_values, '-o', color=colors[j],
                       label=f'{ref}', linewidth=1.8, markersize=5, alpha=0.85)
        
        # Mark default EVR index (18)
        ax_coverage.axvline(x=18, color='red', linestyle='--', alpha=0.6, linewidth=1.5)
        ax_raw.axvline(x=18, color='red', linestyle='--', alpha=0.6, linewidth=1.5)
        
        # Format coverage axis (top row)
        ax_coverage.set_xlabel('EVR Index')
        ax_coverage.set_ylabel('Coverage (%)')
        ax_coverage.set_title(f'{cfg["label"]}')
        ax_coverage.set_xticks([1, 5, 10, 15, 20, 25, 30])
        ax_coverage.grid(True, alpha=0.3, linestyle='--')
        ax_coverage.legend(loc='upper right', fontsize=8, title='Ref Size', 
                          ncol=1, framealpha=0.9)
        ax_coverage.set_xlim(0.5, 30.5)
        
        # Format raw counts axis (bottom row)
        ax_raw.set_xlabel('EVR Index')
        ax_raw.set_ylabel('Rare Cells Captured')
        ax_raw.set_title(f'{cfg["label"]} - Raw Counts')
        ax_raw.set_xticks([1, 5, 10, 15, 20, 25, 30])
        ax_raw.grid(True, alpha=0.3, linestyle='--')
        ax_raw.legend(loc='upper right', fontsize=8, title='Ref Size',
                     ncol=1, framealpha=0.9)
        ax_raw.set_xlim(0.5, 30.5)
    
    # Main title with sample size
    sample_size_k = sample_size // 1000
    fig.suptitle(f'EVR Index Sensitivity Analysis (EVR 1-30, n={sample_size_k}k samples)\n'
                 'Top: Coverage % (rare cells captured / total rare cells in reference)\n'
                 'Red dashed line = default EVR index (18)',
                 fontsize=13, y=1.02)
    
    plt.tight_layout()
    
    # Save figure with sample size suffix
    if sample_size == 100000:
        # Default name for backward compatibility
        out_path = os.path.join(FIGURE_DIR, 'supp_evr_sensitivity_coverage.jpg')
    else:
        out_path = os.path.join(FIGURE_DIR, f'supp_evr_sensitivity_coverage_{sample_size_k}k.jpg')
    
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    print(f"Saved: {out_path}")
    
    plt.close()
    return out_path


def main():
    # Get total rare cells for each dataset/reference (only need to compute once)
    print("Counting total rare cells in each reference...")
    total_rare_cells = {}
    for dataset, cfg in DATASETS.items():
        total_rare_cells[dataset] = get_total_rare_cells_per_reference(dataset, cfg)
    
    # Generate plots for multiple sample sizes
    sample_sizes = [50000, 100000, 200000]
    
    for sample_size in sample_sizes:
        print(f"\n{'='*60}")
        print(f"Generating plot for sample size: {sample_size:,}")
        print('='*60)
        plot_evr_coverage_for_sample_size(sample_size, total_rare_cells)
    
    print(f"\nAll plots generated successfully!")


if __name__ == "__main__":
    main()
