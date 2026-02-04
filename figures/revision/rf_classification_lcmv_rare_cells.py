#!/usr/bin/env python
# coding: utf-8

"""
LCMV RF Classification Line Plot - All 5 Rare Cell Types

Creates a line plot showing delta F1 (SPS - Random) for each of the 5 rare cell types
in the LCMV dataset across feature indices.

Rare cell types:
- interacting (0.07%)
- NK1_1_TCRgd_T (0.07%)
- cDC2 (0.39%)
- pDCs (0.56%)
- CD4_LCMV_spec (0.97%)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
output_dir = script_dir

# Data path
LCMV_DATA_PATH = os.path.join(project_root, "jobs/feature_index_classification_lcmv/results/all_feature_indices_summary.csv")

# Rare cell types configuration
RARE_CELL_TYPES = {
    'interacting': {
        'col': 'interacting_f1_improvement',
        'label': 'interacting (0.07%)',
        'color': '#e74c3c',  # Red
        'marker': 'o',
    },
    'NK1_1_TCRgd_T': {
        'col': 'NK1_1_TCRgd_T_f1_improvement',
        'label': 'NK1_1_TCRgd_T (0.07%)',
        'color': '#9b59b6',  # Purple
        'marker': 's',
    },
    'cDC2': {
        'col': 'cDC2_f1_improvement',
        'label': 'cDC2 (0.39%)',
        'color': '#2ecc71',  # Green
        'marker': '^',
    },
    'pDCs': {
        'col': 'pDCs_f1_improvement',
        'label': 'pDCs (0.56%)',
        'color': '#f39c12',  # Orange
        'marker': 'D',
    },
    'CD4_LCMV_spec': {
        'col': 'CD4_LCMV_spec_f1_improvement',
        'label': 'CD4_LCMV_spec (0.97%)',
        'color': '#3498db',  # Blue
        'marker': 'v',
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
    'lines.markersize': 8,
    'figure.dpi': 100,
    'savefig.dpi': 350,
})


def load_lcmv_data():
    """Load LCMV classification results."""
    df = pd.read_csv(LCMV_DATA_PATH)
    df = df.sort_values('feature_index')
    print(f"Loaded LCMV data: {len(df)} feature indices")
    return df


def create_lcmv_rare_cell_plot(df):
    """Create line plot for all 5 rare cell types in LCMV."""
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    x = df['feature_index'].values
    
    # Plot each rare cell type
    for cell_type, config in RARE_CELL_TYPES.items():
        col = config['col']
        if col in df.columns:
            ax.plot(x, df[col].values, 
                    f'-{config["marker"]}', 
                    color=config['color'],
                    label=config['label'],
                    alpha=0.85,
                    linewidth=2.5,
                    markersize=8)
    
    # Also add macro F1 for reference
    ax.plot(x, df['macro_f1_improvement'].values, 
            '--', color='gray', label='Δ Macro F1', 
            alpha=0.6, linewidth=2)
    
    # Zero line
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.6)
    
    # Labels
    ax.set_xlabel('Feature Index', fontsize=14)
    ax.set_ylabel('ΔF1 (SPS - Random) %', fontsize=14)
    ax.set_title('LCMV Dataset (34M cells): Rare Cell Type F1 Improvement by Feature Index', 
                 fontsize=16, fontweight='bold')
    
    # Legend
    ax.legend(loc='upper left', fontsize=11, ncol=2)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # X-axis
    ax.set_xlim(0, max(x) + 1)
    ax.set_xticks(np.arange(0, max(x) + 1, 2))
    
    # Annotate max improvement for interacting (most dramatic)
    best_idx = df['interacting_f1_improvement'].idxmax()
    best_fi = df.loc[best_idx, 'feature_index']
    best_val = df.loc[best_idx, 'interacting_f1_improvement']
    ax.annotate(f'Max: +{best_val:.1f}%\n(FI={int(best_fi)})', 
               xy=(best_fi, best_val),
               xytext=(best_fi + 2, best_val - 5),
               arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5),
               fontsize=11, color='darkred', fontweight='bold')
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(output_dir, 'rf_classification_lcmv_rare_cells.jpg')
    plt.savefig(output_path, bbox_inches='tight', dpi=350)
    print(f"\nSaved: {output_path}")
    
    pdf_path = os.path.join(output_dir, 'rf_classification_lcmv_rare_cells.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {pdf_path}")
    
    plt.close()


def print_summary(df):
    """Print summary of rare cell type improvements."""
    
    print("\n" + "="*70)
    print("LCMV Rare Cell Type F1 Improvements Summary")
    print("="*70)
    
    for cell_type, config in RARE_CELL_TYPES.items():
        col = config['col']
        if col in df.columns:
            best_idx = df[col].idxmax()
            best_fi = df.loc[best_idx, 'feature_index']
            best_val = df.loc[best_idx, col]
            
            # Values at common feature indices
            fi_18 = df[df['feature_index'] == 18][col].values[0] if 18 in df['feature_index'].values else None
            fi_25 = df[df['feature_index'] == 25][col].values[0] if 25 in df['feature_index'].values else None
            
            print(f"\n{config['label']}:")
            print(f"  Best FI: {int(best_fi)} (Δ={best_val:+.2f}%)")
            if fi_18 is not None:
                print(f"  At FI=18: Δ={fi_18:+.2f}%")
            if fi_25 is not None:
                print(f"  At FI=25: Δ={fi_25:+.2f}%")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    # Load data
    df = load_lcmv_data()
    
    # Print summary
    print_summary(df)
    
    # Create plot
    print("\nGenerating LCMV rare cell types plot...")
    create_lcmv_rare_cell_plot(df)
    
    print("\nDone!")
