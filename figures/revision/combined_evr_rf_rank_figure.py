#!/usr/bin/env python
# coding: utf-8

"""
Combined EVR Analysis Figure (3x4 Layout)

Creates a comprehensive figure with:
- Row 1: Coverage % (rare cells captured / total rare cells) by EVR Index
- Row 2: RF Classification Delta F1 (SPS - Random) by Feature Index
- Row 3: SPS Rank among sampling methods by EVR Index

Columns: MCC (1%), MCC (0.5%), MCC (0.1%), LCMV
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
output_dir = script_dir

# Add analysis directory to path for rare cell type identification
analysis_dir = os.path.join(project_root, 'analysis')
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

# ============================================================================
# Configuration
# ============================================================================

# Dataset order for columns (left to right)
DATASET_ORDER = ['mcc', 'mcc_05', 'mcc_01', 'lcmv']

# Dataset configurations
DATASETS = {
    'mcc': {
        'label': 'MCC (1%)',
        'refs': ['0.5M', '1M', '2M', '2.5M', '3.2M'],  # Display labels
        'ref_cols': ['5M', '10M', '20M', '25M', '30M'],  # CSV column names
        'ref_nums': [5, 10, 20, 25, 30],
        'benchmark_path': os.path.join(project_root, 'data', 'mcc', 'benchmark'),
        'rf_path': os.path.join(project_root, 'jobs/feature_index_classification/results/all_feature_indices_summary.csv'),
        'table_path': os.path.join(project_root, 'jobs/test_feature_index/tables/mcc_feature_index_table_size_100000.csv'),
        'rare_cols_rf': ['osteoblast_f1_improvement'],
        'rare_col_labels': ['Osteoblast'],
        'rare_cell_types': ['osteoblast'],
        'frequency_threshold_pct': 1.0,
    },
    'mcc_05': {
        'label': 'MCC (0.5%)',
        'refs': ['0.5M', '1M', '2M', '2.5M', '3.2M'],  # Display labels
        'ref_cols': ['5M', '10M', '20M', '25M', '30M'],  # CSV column names
        'ref_nums': [5, 10, 20, 25, 30],
        'benchmark_path': os.path.join(project_root, 'data', 'mcc_05', 'benchmark'),
        'rf_path': os.path.join(project_root, 'jobs/feature_index_classification/results_mcc_05/all_feature_indices_summary.csv'),
        'table_path': os.path.join(project_root, 'jobs/test_feature_index/tables/mcc_05_feature_index_table_size_100000.csv'),
        'rare_cols_rf': ['osteoblast_f1_improvement'],
        'rare_col_labels': ['Osteoblast'],
        'rare_cell_types': ['osteoblast'],
        'frequency_threshold_pct': 0.5,
    },
    'mcc_01': {
        'label': 'MCC (0.1%)',
        'refs': ['0.5M', '1M', '2M', '2.5M', '3.2M'],  # Display labels
        'ref_cols': ['5M', '10M', '20M', '25M', '30M'],  # CSV column names
        'ref_nums': [5, 10, 20, 25, 30],
        'benchmark_path': os.path.join(project_root, 'data', 'mcc_01', 'benchmark'),
        'rf_path': os.path.join(project_root, 'jobs/feature_index_classification/results_mcc_01/all_feature_indices_summary.csv'),
        'table_path': os.path.join(project_root, 'jobs/test_feature_index/tables/mcc_01_feature_index_table_size_100000.csv'),
        'rare_cols_rf': ['osteoblast_f1_improvement'],
        'rare_col_labels': ['Osteoblast'],
        'rare_cell_types': ['osteoblast'],
        'frequency_threshold_pct': 0.1,
    },
    'lcmv': {
        'label': 'LCMV',
        'refs': ['1M', '5M', '10M', '20M', '34M'],  # Display labels (same as CSV column names)
        'ref_cols': ['1M', '5M', '10M', '20M', '34M'],  # CSV column names
        'ref_nums': [1, 5, 10, 20, 34],
        'benchmark_path': os.path.join(project_root, 'data', 'lcmv', 'benchmark'),
        'rf_path': os.path.join(project_root, 'jobs/feature_index_classification_lcmv/results/all_feature_indices_summary.csv'),
        'table_path': os.path.join(project_root, 'jobs/test_feature_index/tables/lcmv_feature_index_table_size_100000.csv'),
        'rare_cols_rf': ['interacting_f1_improvement', 'NK1_1_TCRgd_T_f1_improvement', 
                         'cDC2_f1_improvement', 'pDCs_f1_improvement', 'CD4_LCMV_spec_f1_improvement'],
        'rare_col_labels': ['Interacting', 'NK1_1_TCRgd_T', 'cDC2', 'pDCs', 'CD4_LCMV_spec'],
        'rare_cell_types': ['interacting', 'NK1_1_TCRgd_T', 'cDC2', 'pDCs', 'CD4_LCMV_spec'],
        'frequency_threshold_pct': 1.0,
    },
}

# Methods for ranking comparison
METHODS = ['random', 'sps', 'hopper', 'atomic', 'scsampler']
SAMPLE_SIZE = 100000
EVR_RANGE = list(range(1, 31))
label_key = 'celltype'

# Pre-computed total rare cells per reference (from obs.csv counts)
# This avoids having to load large obs.csv files
TOTAL_RARE_CELLS = {
    'lcmv': {'1M': 20226, '5M': 102561, '10M': 205573, '20M': 410544, '34M': 705738},
    'mcc': {'5M': 4750, '10M': 9399, '20M': 18884, '25M': 23891, '30M': 30000},
    'mcc_05': {'5M': 2396, '10M': 4652, '20M': 9490, '25M': 11984, '30M': 15000},
    'mcc_01': {'5M': 458, '10M': 929, '20M': 1912, '25M': 2377, '30M': 3000},
}


# ============================================================================
# Data Loading Functions
# ============================================================================

def load_coverage_data(dataset_key, sample_size=SAMPLE_SIZE):
    """
    Load coverage data from feature_index tables and compute coverage percentage.
    
    Returns:
        dict: {ref_name: {'raw': array, 'pct': array, 'total_rare': int}}
    """
    cfg = DATASETS[dataset_key]
    
    # Build table path with sample size
    table_path = os.path.join(project_root, 'jobs/test_feature_index/tables',
                              f'{dataset_key}_feature_index_table_size_{sample_size}.csv')
    
    if not os.path.exists(table_path):
        print(f"Warning: Table not found for {dataset_key}: {table_path}")
        return None
    
    # Load the feature index table (SPS rare cell counts)
    df = pd.read_csv(table_path, index_col=0)
    
    # Get total rare cells in each reference
    total_rare_per_ref = get_total_rare_cells(dataset_key)
    
    coverage_data = {}
    ref_cols = cfg.get('ref_cols', cfg['refs'])  # Use ref_cols if available, else use refs
    for ref_label, ref_col, ref_num in zip(cfg['refs'], ref_cols, cfg['ref_nums']):
        if ref_col not in df.columns:
            continue
        
        raw_values = df[ref_col].values
        total_rare = total_rare_per_ref.get(ref_num, None)
        
        if total_rare is not None and total_rare > 0:
            pct_values = (raw_values / total_rare) * 100
        else:
            pct_values = np.zeros_like(raw_values)
        
        coverage_data[ref_label] = {  # Use ref_label for the key (display label)
            'raw': raw_values,
            'pct': pct_values,
            'total_rare': total_rare,
            'feature_indices': df.index.values
        }
    
    return coverage_data


def get_total_rare_cells(dataset_key):
    """
    Get total rare cell counts for each reference size from pre-computed dictionary.
    """
    cfg = DATASETS[dataset_key]
    total_rare_per_ref = {}
    
    if dataset_key not in TOTAL_RARE_CELLS:
        print(f"Warning: No pre-computed rare cell counts for {dataset_key}")
        return total_rare_per_ref
    
    rare_counts = TOTAL_RARE_CELLS[dataset_key]
    
    for ref_num in cfg['ref_nums']:
        ref_key = f'{ref_num}M'
        if ref_key in rare_counts:
            total_rare_per_ref[ref_num] = rare_counts[ref_key]
        else:
            print(f"Warning: No rare cell count for {dataset_key} ref {ref_key}")
            total_rare_per_ref[ref_num] = None
    
    return total_rare_per_ref


def load_rf_data(dataset_key):
    """
    Load RF classification results for a dataset.
    
    Returns:
        DataFrame with feature_index, macro_f1_improvement, and rare cell F1 improvements
    """
    cfg = DATASETS[dataset_key]
    rf_path = cfg['rf_path']
    
    if not os.path.exists(rf_path):
        print(f"Warning: RF data not found for {dataset_key}: {rf_path}")
        return None
    
    df = pd.read_csv(rf_path)
    df = df.sort_values('feature_index')
    
    return df


def compute_sps_rank_for_dataset(dataset_key, ref_num, sample_size=SAMPLE_SIZE):
    """
    Compute SPS rank among methods for each feature index.
    
    Uses pre-computed tables:
    - Feature index table: SPS coverage for each feature index
    - Sampling methods table: Coverage for other methods (random, hopper, atomic, scsampler)
    
    Returns:
        dict: {feature_index: rank} where rank 1 is best (most rare cells captured)
    """
    cfg = DATASETS[dataset_key]
    
    # Load other methods' coverage from pre-computed sampling methods table
    methods_table_path = os.path.join(project_root, 'jobs', 'test_sampling_methods', 'tables',
                                       f'{dataset_key}_sampling_methods_table_size_{sample_size}.csv')
    
    if not os.path.exists(methods_table_path):
        print(f"Warning: Sampling methods table not found: {methods_table_path}")
        return None
    
    methods_df = pd.read_csv(methods_table_path, index_col=0)
    ref_col = f'{ref_num}M'
    
    if ref_col not in methods_df.columns:
        print(f"Warning: Reference column {ref_col} not found in methods table")
        return None
    
    # Get other methods' coverage for this reference size
    other_methods_coverage = {}
    for method in ['random', 'hopper', 'atomic', 'scsampler']:
        if method in methods_df.index:
            other_methods_coverage[method] = methods_df.loc[method, ref_col]
        else:
            other_methods_coverage[method] = 0
    
    # Load SPS coverage for each feature index from the feature index table
    table_path = os.path.join(project_root, 'jobs/test_feature_index/tables',
                              f'{dataset_key}_feature_index_table_size_{sample_size}.csv')
    if not os.path.exists(table_path):
        print(f"Warning: Feature index table not found: {table_path}")
        return None
    
    fi_df = pd.read_csv(table_path, index_col=0)
    
    if ref_col not in fi_df.columns:
        print(f"Warning: Reference column {ref_col} not found in feature index table")
        return None
    
    sps_coverage = dict(zip(fi_df.index, fi_df[ref_col].values))
    
    # Compute rank for each feature index
    ranks = {}
    for fi in EVR_RANGE:
        if fi not in sps_coverage:
            ranks[fi] = np.nan
            continue
        
        sps_cov = sps_coverage[fi]
        
        # Collect all method coverages
        all_coverages = [
            ('random', other_methods_coverage.get('random', 0)),
            ('sps', sps_cov),
            ('hopper', other_methods_coverage.get('hopper', 0)),
            ('atomic', other_methods_coverage.get('atomic', 0)),
            ('scsampler', other_methods_coverage.get('scsampler', 0)),
        ]
        
        # Sort by coverage (descending) and find SPS rank
        sorted_methods = sorted(all_coverages, key=lambda x: x[1], reverse=True)
        
        for rank, (method, _) in enumerate(sorted_methods, start=1):
            if method == 'sps':
                ranks[fi] = rank
                break
    
    return ranks


def compute_sps_ranks_averaged(dataset_key, sample_size=SAMPLE_SIZE):
    """
    Compute SPS ranks using the largest reference size (full dataset).
    
    For MCC datasets: 30M (3.2M cells)
    For LCMV: 34M cells
    
    Returns:
        dict: {feature_index: rank}
    """
    cfg = DATASETS[dataset_key]
    
    # Use only the largest reference (full dataset)
    largest_ref = cfg['ref_nums'][-1]
    
    ranks = compute_sps_rank_for_dataset(dataset_key, largest_ref, sample_size)
    
    if ranks is None:
        return {fi: np.nan for fi in EVR_RANGE}
    
    return ranks


# ============================================================================
# Plotting Functions
# ============================================================================

def plot_combined_figure(sample_size=SAMPLE_SIZE):
    """
    Create the combined 3x4 figure for a specific sample size.
    
    Args:
        sample_size: Sample size (e.g., 50000, 100000, 200000)
    """
    # Plot settings
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'legend.fontsize': 7,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'lines.markersize': 4,
        'lines.linewidth': 1.5,
        'figure.dpi': 100,
        'savefig.dpi': 350,
    })
    
    # Create figure: 3 rows x 4 columns
    fig, axes = plt.subplots(3, 4, figsize=(16, 11))
    
    # Color schemes
    coverage_colors = plt.cm.viridis(np.linspace(0.2, 0.9, 5))
    rf_color_macro = '#3498db'  # Blue
    rf_color_rare = '#e74c3c'   # Red
    rf_rare_colors = ['#e74c3c', '#9b59b6', '#2ecc71', '#f39c12', '#1abc9c']  # For LCMV multiple rare types
    rank_color = '#2c3e50'
    
    # Plot each dataset (column)
    for col_idx, dataset_key in enumerate(DATASET_ORDER):
        cfg = DATASETS[dataset_key]
        print(f"  Processing {cfg['label']}...")
        
        # ====================================================================
        # Row 1: Coverage %
        # ====================================================================
        print(f"    Loading coverage data...")
        ax_cov = axes[0, col_idx]
        coverage_data = load_coverage_data(dataset_key, sample_size)
        
        if coverage_data is not None:
            for j, ref_name in enumerate(cfg['refs']):
                if ref_name in coverage_data:
                    data = coverage_data[ref_name]
                    x = data['feature_indices']
                    y = data['pct']
                    ax_cov.plot(x, y, '-o', color=coverage_colors[j % len(coverage_colors)],
                               label=ref_name, markersize=3, alpha=0.85)
        
        # Mark default EVR index (12)
        ax_cov.axvline(x=12, color='red', linestyle='--', alpha=0.6, linewidth=1.2)
        
        ax_cov.set_xlabel('EVR Index')
        ax_cov.set_ylabel('Coverage (%)')
        ax_cov.set_title(f'{cfg["label"]}', fontweight='bold')
        ax_cov.set_xticks([1, 5, 10, 15, 20, 25, 30])
        ax_cov.set_xlim(0.5, 30.5)
        ax_cov.grid(True, alpha=0.3, linestyle='--')
        ax_cov.legend(loc='upper right', fontsize=6, title='Ref Size', ncol=1)
        
        # ====================================================================
        # Row 2: RF Delta F1
        # ====================================================================
        print(f"    Loading RF classification data...")
        ax_rf = axes[1, col_idx]
        rf_data = load_rf_data(dataset_key)
        
        if rf_data is not None:
            x = rf_data['feature_index'].values
            
            # Plot Delta Macro F1
            ax_rf.plot(x, rf_data['macro_f1_improvement'].values, '-o',
                      color=rf_color_macro, label='Δ Macro F1', alpha=0.9, markersize=3)
            
            # Plot Delta F1 for rare cell types
            for i, (rare_col, rare_label) in enumerate(zip(cfg['rare_cols_rf'], cfg['rare_col_labels'])):
                if rare_col in rf_data.columns:
                    color = rf_rare_colors[i % len(rf_rare_colors)]
                    ax_rf.plot(x, rf_data[rare_col].values, '-s',
                              color=color, label=f'Δ {rare_label} F1', 
                              alpha=0.8, markersize=3)
        
        # Zero line
        ax_rf.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
        
        ax_rf.set_xlabel('EVR Index')
        ax_rf.set_ylabel('ΔF1 (SPS - Random) %')
        ax_rf.set_title(f'{cfg["label"]}', fontweight='bold')
        ax_rf.set_xticks([1, 5, 10, 15, 20, 25, 30])
        ax_rf.set_xlim(0.5, 30.5)
        ax_rf.grid(True, alpha=0.3, linestyle='--')
        
        # Legend placement based on number of lines
        if dataset_key == 'lcmv':
            ax_rf.legend(loc='lower right', fontsize=5, ncol=2)
        else:
            ax_rf.legend(loc='best', fontsize=7)
        
        # ====================================================================
        # Row 3: SPS Rank
        # ====================================================================
        print(f"    Computing SPS rank...")
        ax_rank = axes[2, col_idx]
        avg_ranks = compute_sps_ranks_averaged(dataset_key, sample_size)
        
        if avg_ranks:
            x = list(avg_ranks.keys())
            y = list(avg_ranks.values())
            ax_rank.plot(x, y, '-o', color=rank_color, linewidth=2, markersize=4, alpha=0.9)
            
            # Fill area above rank 1 to show where SPS is not #1
            ax_rank.fill_between(x, 1, y, where=np.array(y) > 1, 
                                alpha=0.15, color='red', label='Not #1')
            ax_rank.fill_between(x, y, 1, where=np.array(y) <= 1, 
                                alpha=0.15, color='green', label='#1')
        
        # Mark default EVR index (12)
        ax_rank.axvline(x=12, color='red', linestyle='--', alpha=0.6, linewidth=1.2)
        
        # Horizontal line at rank 1
        ax_rank.axhline(y=1, color='green', linestyle=':', alpha=0.5, linewidth=1)
        
        ax_rank.set_xlabel('EVR Index')
        ax_rank.set_ylabel('SPS Rank (1=Best)')
        ax_rank.set_title(f'{cfg["label"]}', fontweight='bold')
        ax_rank.set_xticks([1, 5, 10, 15, 20, 25, 30])
        ax_rank.set_xlim(0.5, 30.5)
        ax_rank.set_ylim(0.5, 5.5)
        ax_rank.set_yticks([1, 2, 3, 4, 5])
        ax_rank.invert_yaxis()  # Rank 1 at top
        ax_rank.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3, wspace=0.25)
    
    # Save figure with sample size in filename
    sample_size_k = sample_size // 1000
    output_path = os.path.join(output_dir, f'combined_evr_rf_rank_figure_{sample_size_k}k.jpg')
    plt.savefig(output_path, bbox_inches='tight', dpi=350)
    print(f"Saved: {output_path}")
    
    pdf_path = os.path.join(output_dir, f'combined_evr_rf_rank_figure_{sample_size_k}k.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {pdf_path}")
    
    plt.close()
    
    return output_path


# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    print("Creating Combined EVR Analysis Figures (3x4 Layout)")
    print("=" * 60)
    
    # Sample sizes to generate figures for
    SAMPLE_SIZES = [50000, 100000, 200000]
    
    # Load and check data availability
    print("\nChecking data availability...")
    for dataset_key in DATASET_ORDER:
        cfg = DATASETS[dataset_key]
        print(f"\n{cfg['label']}:")
        print(f"  RF Data: {'Found' if os.path.exists(cfg['rf_path']) else 'NOT FOUND'}")
        print(f"  Benchmark: {'Found' if os.path.exists(cfg['benchmark_path']) else 'NOT FOUND'}")
        for sample_size in SAMPLE_SIZES:
            table_path = os.path.join(project_root, 'jobs/test_feature_index/tables',
                                      f'{dataset_key}_feature_index_table_size_{sample_size}.csv')
            print(f"  Table ({sample_size//1000}k): {'Found' if os.path.exists(table_path) else 'NOT FOUND'}")
    
    # Create figures for each sample size
    for sample_size in SAMPLE_SIZES:
        print("\n" + "=" * 60)
        print(f"Generating figure for sample size: {sample_size:,}")
        print("=" * 60)
        output_path = plot_combined_figure(sample_size)
    
    print("\n" + "=" * 60)
    print("All figures generated!")
    print("Done!")
