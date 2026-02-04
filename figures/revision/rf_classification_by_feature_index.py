#!/usr/bin/env python
# coding: utf-8

"""
RF Classification by Feature Index Line Plot

Creates line plots showing:
1. Delta Macro F1 (SPS - Random) by feature index
2. Delta F1 for rare cell types by feature index

Uses data from feature_index_classification results for all datasets:
- MCC (3.2M cells, 0.95% rare)
- MCC_01 (3.1M cells, 0.1% rare)
- MCC_05 (3.1M cells, 0.5% rare)
- LCMV (34M cells)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
output_dir = script_dir

# Data paths for all datasets
DATA_PATHS = {
    'mcc': os.path.join(project_root, "jobs/feature_index_classification/results/all_feature_indices_summary.csv"),
    'mcc_01': os.path.join(project_root, "jobs/feature_index_classification/results_mcc_01/all_feature_indices_summary.csv"),
    'mcc_05': os.path.join(project_root, "jobs/feature_index_classification/results_mcc_05/all_feature_indices_summary.csv"),
    'lcmv': os.path.join(project_root, "jobs/feature_index_classification_lcmv/results/all_feature_indices_summary.csv"),
}

# Dataset configurations
DATASET_CONFIGS = {
    'mcc': {
        'label': 'MCC (3.2M, 0.95% rare)',
        'rare_cell_type': 'osteoblast',
        'rare_col': 'osteoblast_f1_improvement',
    },
    'mcc_01': {
        'label': 'MCC_01 (3.1M, 0.1% rare)',
        'rare_cell_type': 'osteoblast',
        'rare_col': 'osteoblast_f1_improvement',
    },
    'mcc_05': {
        'label': 'MCC_05 (3.1M, 0.5% rare)',
        'rare_cell_type': 'osteoblast',
        'rare_col': 'osteoblast_f1_improvement',
    },
    'lcmv': {
        'label': 'LCMV (34M)',
        'rare_cell_type': 'interacting',
        'rare_col': 'interacting_f1_improvement',
    },
}

# Legacy paths for backward compatibility
MCC_DATA_PATH = DATA_PATHS['mcc']
LCMV_DATA_PATH = DATA_PATHS['lcmv']

# Configure plotting style
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'legend.fontsize': 10,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'figure.figsize': [16, 10],
    'savefig.dpi': 350,
    'figure.dpi': 100,
    'savefig.format': 'jpg',
    'lines.linewidth': 2,
    'lines.markersize': 6
})


def load_data():
    """Load the feature index classification results for MCC and LCMV (legacy)."""
    mcc_df = pd.read_csv(MCC_DATA_PATH)
    lcmv_df = pd.read_csv(LCMV_DATA_PATH)
    
    # Sort by feature index
    mcc_df = mcc_df.sort_values('feature_index')
    lcmv_df = lcmv_df.sort_values('feature_index')
    
    print(f"Loaded MCC data: {len(mcc_df)} feature indices")
    print(f"Loaded LCMV data: {len(lcmv_df)} feature indices")
    
    return mcc_df, lcmv_df


def load_all_data():
    """Load the feature index classification results for all datasets."""
    data = {}
    
    for dataset, path in DATA_PATHS.items():
        if os.path.exists(path):
            df = pd.read_csv(path)
            df = df.sort_values('feature_index')
            data[dataset] = df
            print(f"Loaded {dataset} data: {len(df)} feature indices")
        else:
            print(f"Warning: {dataset} data not found at {path}")
    
    return data


def create_combined_figure(mcc_df, lcmv_df):
    """Create a combined figure with both datasets."""
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Colors for different metrics
    colors = {
        'macro_f1': '#3498db',       # Blue
        'osteoblast': '#e74c3c',     # Red
        'interacting': '#e74c3c',    # Red
        'NK1_1_TCRgd_T': '#9b59b6',  # Purple  
        'cDC2': '#2ecc71',           # Green
        'pDCs': '#f39c12',           # Orange
        'CD4_LCMV_spec': '#1abc9c',  # Teal
    }
    
    # =========================================================================
    # Panel A: MCC - Delta Macro F1 and Delta Osteoblast F1
    # =========================================================================
    ax1 = axes[0, 0]
    
    x = mcc_df['feature_index'].values
    
    # Delta Macro F1
    ax1.plot(x, mcc_df['macro_f1_improvement'].values, '-o', 
             color=colors['macro_f1'], label='Δ Macro F1', alpha=0.8)
    
    # Delta Osteoblast F1
    ax1.plot(x, mcc_df['osteoblast_f1_improvement'].values, '-s', 
             color=colors['osteoblast'], label='Δ Osteoblast F1 (rare)', alpha=0.8)
    
    ax1.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax1.set_xlabel('Feature Index')
    ax1.set_ylabel('ΔF1 (SPS - Random) %')
    ax1.set_title('A. MCC Dataset (3.2M cells)\nRare cell type: Osteoblast (0.95%)', fontweight='bold')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlim(0, 31)
    
    # Highlight the region where osteoblast F1 improves
    positive_region = mcc_df[mcc_df['osteoblast_f1_improvement'] > 0]['feature_index'].values
    if len(positive_region) > 0:
        ax1.axvspan(positive_region.min() - 0.5, positive_region.max() + 0.5, 
                    alpha=0.1, color='green', label='SPS better for rare cells')
    
    # =========================================================================
    # Panel B: MCC - EVR vs Performance
    # =========================================================================
    ax2 = axes[0, 1]
    
    # Create twin axis for EVR
    ax2_twin = ax2.twinx()
    
    # Plot Delta F1 metrics
    ax2.plot(x, mcc_df['macro_f1_improvement'].values, '-o', 
             color=colors['macro_f1'], label='Δ Macro F1', alpha=0.8)
    ax2.plot(x, mcc_df['osteoblast_f1_improvement'].values, '-s', 
             color=colors['osteoblast'], label='Δ Osteoblast F1', alpha=0.8)
    
    # Plot EVR on twin axis
    ax2_twin.fill_between(x, 0, mcc_df['evr'].values, alpha=0.2, color='gray', label='EVR')
    ax2_twin.plot(x, mcc_df['evr'].values, '--', color='gray', alpha=0.6, linewidth=1)
    
    ax2.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax2.set_xlabel('Feature Index')
    ax2.set_ylabel('ΔF1 (SPS - Random) %')
    ax2_twin.set_ylabel('Explained Variance Ratio (EVR)', color='gray')
    ax2_twin.tick_params(axis='y', labelcolor='gray')
    ax2.set_title('B. MCC: ΔF1 vs EVR', fontweight='bold')
    ax2.legend(loc='upper left')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim(0, 31)
    
    # =========================================================================
    # Panel C: LCMV - Delta Macro F1 and all rare cell types
    # =========================================================================
    ax3 = axes[1, 0]
    
    x_lcmv = lcmv_df['feature_index'].values
    
    # Delta Macro F1
    ax3.plot(x_lcmv, lcmv_df['macro_f1_improvement'].values, '-o', 
             color=colors['macro_f1'], label='Δ Macro F1', alpha=0.8)
    
    # Delta F1 for rare cell types
    ax3.plot(x_lcmv, lcmv_df['interacting_f1_improvement'].values, '-s', 
             color=colors['interacting'], label='Δ Interacting F1 (0.07%)', alpha=0.8)
    ax3.plot(x_lcmv, lcmv_df['NK1_1_TCRgd_T_f1_improvement'].values, '-^', 
             color=colors['NK1_1_TCRgd_T'], label='Δ NK1_1_TCRgd_T F1 (0.07%)', alpha=0.7)
    ax3.plot(x_lcmv, lcmv_df['cDC2_f1_improvement'].values, '-v', 
             color=colors['cDC2'], label='Δ cDC2 F1 (0.39%)', alpha=0.7)
    
    ax3.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax3.set_xlabel('Feature Index')
    ax3.set_ylabel('ΔF1 (SPS - Random) %')
    ax3.set_title('C. LCMV Dataset (34M cells)\nRare cell types: <1% frequency', fontweight='bold')
    ax3.legend(loc='best', fontsize=9)
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.set_xlim(0, 26)
    
    # =========================================================================
    # Panel D: LCMV - Focus on most rare cell type (interacting)
    # =========================================================================
    ax4 = axes[1, 1]
    
    # Focus on the most dramatic improvements
    ax4.plot(x_lcmv, lcmv_df['interacting_f1_improvement'].values, '-o', 
             color=colors['interacting'], label='Δ Interacting F1 (0.07%)', 
             linewidth=2.5, markersize=8, alpha=0.9)
    ax4.plot(x_lcmv, lcmv_df['macro_f1_improvement'].values, '-s', 
             color=colors['macro_f1'], label='Δ Macro F1', alpha=0.7)
    
    # Add EVR on twin axis
    ax4_twin = ax4.twinx()
    ax4_twin.fill_between(x_lcmv, 0, lcmv_df['evr'].values, alpha=0.15, color='gray')
    ax4_twin.plot(x_lcmv, lcmv_df['evr'].values, '--', color='gray', alpha=0.5, linewidth=1)
    ax4_twin.set_ylabel('EVR', color='gray')
    ax4_twin.tick_params(axis='y', labelcolor='gray')
    
    ax4.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax4.set_xlabel('Feature Index')
    ax4.set_ylabel('ΔF1 (SPS - Random) %')
    ax4.set_title('D. LCMV: Focus on "Interacting" (most rare)', fontweight='bold')
    ax4.legend(loc='upper left')
    ax4.grid(True, alpha=0.3, linestyle='--')
    ax4.set_xlim(0, 26)
    
    # Annotate max improvement
    max_idx = lcmv_df['interacting_f1_improvement'].idxmax()
    max_fi = lcmv_df.loc[max_idx, 'feature_index']
    max_val = lcmv_df.loc[max_idx, 'interacting_f1_improvement']
    ax4.annotate(f'Max: +{max_val:.1f}%\n(FI={int(max_fi)})', 
                 xy=(max_fi, max_val), xytext=(max_fi + 3, max_val - 5),
                 arrowprops=dict(arrowstyle='->', color='darkred'),
                 fontsize=10, color='darkred', fontweight='bold')
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'rf_classification_by_feature_index.jpg')
    plt.savefig(output_path, bbox_inches='tight', dpi=350)
    print(f"Saved figure to: {output_path}")
    
    # Also save as PDF
    pdf_path = os.path.join(output_dir, 'rf_classification_by_feature_index.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved figure to: {pdf_path}")
    
    plt.close()
    
    return output_path


def create_summary_figure(mcc_df, lcmv_df):
    """Create a simplified summary figure comparing both datasets."""
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    colors = {
        'macro_f1': '#3498db',
        'rare': '#e74c3c',
    }
    
    # =========================================================================
    # Left: MCC
    # =========================================================================
    ax1 = axes[0]
    x = mcc_df['feature_index'].values
    
    ax1.plot(x, mcc_df['macro_f1_improvement'].values, '-o', 
             color=colors['macro_f1'], label='Δ Macro F1', alpha=0.8, linewidth=2)
    ax1.plot(x, mcc_df['osteoblast_f1_improvement'].values, '-s', 
             color=colors['rare'], label='Δ Osteoblast F1 (0.95%)', alpha=0.8, linewidth=2)
    
    ax1.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax1.fill_between(x, 0, mcc_df['osteoblast_f1_improvement'].values, 
                     where=mcc_df['osteoblast_f1_improvement'] > 0, 
                     alpha=0.2, color=colors['rare'])
    
    ax1.set_xlabel('Feature Index', fontsize=14)
    ax1.set_ylabel('ΔF1 (SPS - Random) %', fontsize=14)
    ax1.set_title('MCC Dataset (3.2M cells)', fontweight='bold', fontsize=16)
    ax1.legend(loc='lower right', fontsize=11)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlim(0, 31)
    ax1.set_ylim(-50, 10)
    
    # Add annotation for optimal region
    best_osteo_idx = mcc_df['osteoblast_f1_improvement'].idxmax()
    best_fi = mcc_df.loc[best_osteo_idx, 'feature_index']
    best_val = mcc_df.loc[best_osteo_idx, 'osteoblast_f1_improvement']
    ax1.annotate(f'Best: +{best_val:.1f}%\n(FI={int(best_fi)})', 
                 xy=(best_fi, best_val), xytext=(best_fi - 8, best_val + 3),
                 arrowprops=dict(arrowstyle='->', color='darkred'),
                 fontsize=10, color='darkred', fontweight='bold')
    
    # =========================================================================
    # Right: LCMV  
    # =========================================================================
    ax2 = axes[1]
    x_lcmv = lcmv_df['feature_index'].values
    
    ax2.plot(x_lcmv, lcmv_df['macro_f1_improvement'].values, '-o', 
             color=colors['macro_f1'], label='Δ Macro F1', alpha=0.8, linewidth=2)
    ax2.plot(x_lcmv, lcmv_df['interacting_f1_improvement'].values, '-s', 
             color=colors['rare'], label='Δ Interacting F1 (0.07%)', alpha=0.8, linewidth=2)
    
    ax2.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax2.fill_between(x_lcmv, 0, lcmv_df['interacting_f1_improvement'].values, 
                     where=lcmv_df['interacting_f1_improvement'].values > 0, 
                     alpha=0.2, color=colors['rare'])
    
    ax2.set_xlabel('Feature Index', fontsize=14)
    ax2.set_ylabel('ΔF1 (SPS - Random) %', fontsize=14)
    ax2.set_title('LCMV Dataset (34M cells)', fontweight='bold', fontsize=16)
    ax2.legend(loc='lower right', fontsize=11)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim(0, 26)
    
    # Add annotation for max improvement
    max_idx = lcmv_df['interacting_f1_improvement'].idxmax()
    max_fi = lcmv_df.loc[max_idx, 'feature_index']
    max_val = lcmv_df.loc[max_idx, 'interacting_f1_improvement']
    ax2.annotate(f'Max: +{max_val:.1f}%\n(FI={int(max_fi)})', 
                 xy=(max_fi, max_val), xytext=(max_fi - 6, max_val - 8),
                 arrowprops=dict(arrowstyle='->', color='darkred'),
                 fontsize=10, color='darkred', fontweight='bold')
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'rf_classification_by_feature_index_summary.jpg')
    plt.savefig(output_path, bbox_inches='tight', dpi=350)
    print(f"Saved summary figure to: {output_path}")
    
    plt.close()
    
    return output_path


def create_all_datasets_figure(data):
    """Create a combined figure for all four datasets (MCC, MCC_01, MCC_05, LCMV)."""
    
    # Filter to available datasets
    available = {k: v for k, v in data.items() if v is not None}
    n_datasets = len(available)
    
    if n_datasets == 0:
        print("No data available!")
        return None
    
    # Determine grid layout
    if n_datasets <= 2:
        nrows, ncols = 1, n_datasets
        figsize = (7 * n_datasets, 6)
    else:
        nrows, ncols = 2, 2
        figsize = (14, 12)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if n_datasets == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Colors
    colors = {
        'macro_f1': '#3498db',  # Blue
        'rare': '#e74c3c',      # Red
    }
    
    # Plot each dataset
    for idx, (dataset, df) in enumerate(available.items()):
        ax = axes[idx]
        config = DATASET_CONFIGS[dataset]
        
        x = df['feature_index'].values
        
        # Delta Macro F1
        ax.plot(x, df['macro_f1_improvement'].values, '-o', 
                color=colors['macro_f1'], label='Δ Macro F1', alpha=0.8, linewidth=2)
        
        # Delta rare cell type F1
        rare_col = config['rare_col']
        if rare_col in df.columns:
            ax.plot(x, df[rare_col].values, '-s', 
                    color=colors['rare'], 
                    label=f'Δ {config["rare_cell_type"]} F1', 
                    alpha=0.8, linewidth=2)
            
            # Fill positive region
            ax.fill_between(x, 0, df[rare_col].values, 
                           where=df[rare_col].values > 0, 
                           alpha=0.2, color=colors['rare'])
        
        ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
        ax.set_xlabel('Feature Index', fontsize=12)
        ax.set_ylabel('ΔF1 (SPS - Random) %', fontsize=12)
        ax.set_title(config['label'], fontweight='bold', fontsize=14)
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_xlim(0, max(x) + 1)
        
        # Add annotation for best rare cell type improvement
        if rare_col in df.columns:
            best_idx = df[rare_col].idxmax()
            best_fi = df.loc[best_idx, 'feature_index']
            best_val = df.loc[best_idx, rare_col]
            if best_val > 0:
                ax.annotate(f'Best: +{best_val:.1f}%\n(FI={int(best_fi)})', 
                           xy=(best_fi, best_val), 
                           xytext=(best_fi + 3, best_val * 0.8),
                           arrowprops=dict(arrowstyle='->', color='darkred', lw=1),
                           fontsize=9, color='darkred', fontweight='bold')
    
    # Hide extra subplots if any
    for idx in range(len(available), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'rf_classification_all_datasets.jpg')
    plt.savefig(output_path, bbox_inches='tight', dpi=350)
    print(f"Saved all datasets figure to: {output_path}")
    
    # Also save as PDF
    pdf_path = os.path.join(output_dir, 'rf_classification_all_datasets.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved figure to: {pdf_path}")
    
    plt.close()
    
    return output_path


def print_summary(mcc_df, lcmv_df):
    """Print a summary of key findings."""
    
    print("\n" + "="*70)
    print("Feature Index Classification Results Summary")
    print("="*70)
    
    # MCC
    print("\n--- MCC Dataset (3.2M cells) ---")
    best_osteo_idx = mcc_df['osteoblast_f1_improvement'].idxmax()
    best_osteo = mcc_df.loc[best_osteo_idx]
    print(f"  Best Osteoblast F1 improvement: Feature Index {int(best_osteo['feature_index'])}")
    print(f"    Δ Osteoblast F1: {best_osteo['osteoblast_f1_improvement']:+.2f}%")
    print(f"    Δ Macro F1: {best_osteo['macro_f1_improvement']:+.2f}%")
    print(f"    EVR: {best_osteo['evr']:.2f}")
    
    # Feature indices where osteoblast improves
    positive_osteo = mcc_df[mcc_df['osteoblast_f1_improvement'] > 0]['feature_index'].values
    if len(positive_osteo) > 0:
        print(f"\n  Feature indices with SPS improvement for Osteoblast:")
        print(f"    {sorted(positive_osteo)}")
    
    # LCMV
    print("\n--- LCMV Dataset (34M cells) ---")
    best_inter_idx = lcmv_df['interacting_f1_improvement'].idxmax()
    best_inter = lcmv_df.loc[best_inter_idx]
    print(f"  Best Interacting F1 improvement: Feature Index {int(best_inter['feature_index'])}")
    print(f"    Δ Interacting F1 (0.07%): {best_inter['interacting_f1_improvement']:+.2f}%")
    print(f"    Δ Macro F1: {best_inter['macro_f1_improvement']:+.2f}%")
    print(f"    EVR: {best_inter['evr']:.2f}")
    
    # All rare cell types at best feature index
    print(f"\n  All rare cell type improvements at FI={int(best_inter['feature_index'])}:")
    rare_cols = ['interacting', 'NK1_1_TCRgd_T', 'cDC2', 'pDCs', 'CD4_LCMV_spec']
    for col in rare_cols:
        imp_col = f'{col}_f1_improvement'
        if imp_col in lcmv_df.columns:
            print(f"    {col}: {best_inter[imp_col]:+.2f}%")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    # Load all available data
    all_data = load_all_data()
    
    # Create figure for all datasets (if multiple available)
    if len(all_data) > 0:
        create_all_datasets_figure(all_data)
    
    # Also create legacy figures for MCC and LCMV if available
    if 'mcc' in all_data and 'lcmv' in all_data:
        mcc_df = all_data['mcc']
        lcmv_df = all_data['lcmv']
        
        # Print summary
        print_summary(mcc_df, lcmv_df)
        
        # Create figures
        create_combined_figure(mcc_df, lcmv_df)
        create_summary_figure(mcc_df, lcmv_df)
    
    print("\nDone!")
