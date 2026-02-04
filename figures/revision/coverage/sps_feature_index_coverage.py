#!/usr/bin/env python
# coding: utf-8
"""
Coverage analysis for SPS with feature_index=25 on MCC dataset.
Compares SPS (feature_index=25) against other sampling methods across different reference sizes.
Uses indices from /data/test_feature_index/mcc
"""

import os
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import logging
from datetime import datetime
from itertools import product
import matplotlib.patches as mpatches
import scanpy as sc
import scipy.sparse

# Add analysis directory to path to import rare cell type functions
script_dir = os.path.dirname(os.path.abspath(__file__))
analysis_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(script_dir))), 'analysis')
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

from refined_rare_cell_type_definition import (
    calculate_cell_type_distances,
    identify_rare_cell_types_distance_and_frequency
)

# Set up logging
def setup_logger():
    log_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'sps_feature_index_coverage_{timestamp}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

# Constants
MCC_REFERENCES = [5, 10, 20, 25, 30]
# Methods to compare: SPS with feature_index=25, plus other methods from benchmark
METHODS = ['random', 'sps_fi25', 'hopper', 'atomic', 'scsampler']
METHOD_LABELS = ['random', 'sps (fi=25)', 'hopper', 'atomic', 'scsampler']
MCC_SIZES = [50000, 100000, 200000, 300000]
REPS = [0]  # Only rep 0 available in test_feature_index data
FEATURE_INDEX = 25  # Using feature_index=25 as requested
label_key = 'celltype'

# File paths
def get_data_path():
    """Get data path from environment variable or derive from script location."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))
    return os.path.join(project_root, 'data')

def get_figures_dir():
    """Get figures directory path from environment variable or derive from script location."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))
    return os.path.join(project_root, 'figures', 'revision')

file_path_env = get_data_path()
MCC_BENCHMARK_PATH = os.path.join(file_path_env, 'mcc/benchmark')
MCC_FEATURE_INDEX_PATH = os.path.join(file_path_env, 'test_feature_index/mcc')
FIGURE_DIR = get_figures_dir()

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'legend.fontsize': 20,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'figure.figsize': [6, 6],
    "legend.markerscale": 1.5,
    'savefig.dpi': 350,
    'figure.dpi': 350,
    'savefig.format': 'jpg'
})

# --- Data Loading and Processing ---
def load_reference_data(references, path):
    """Load reference observation data from benchmark path."""
    reference_obs = {}
    for ref in references:
        address = os.path.join(path, f"{ref}/obs.csv")
        df = pd.read_csv(address, index_col=0)
        df[label_key] = df[label_key].astype('category')
        reference_obs[ref] = df
    return reference_obs

def calculate_rare_cell_coverage(reference_obs, references, sizes, benchmark_path, feature_index_path, feature_index):
    """
    Calculate rare cell coverage for SPS with specific feature_index and other methods.
    
    Uses benchmark data for random, hopper, atomic, scsampler.
    Uses test_feature_index data for SPS with the specified feature_index.
    """
    res = {}
    rare_type_counts = {}
    keys = [label_key, 'count', 'method', 'ref_size', 'sample_size', 'rep', 'real_count']
    
    # Using 1% frequency threshold (same as original MCC)
    frequency_threshold_pct = 1.0
    
    for ref in references:
        obs = reference_obs[ref]
        real_counts = obs[label_key].value_counts()
        
        # Load adata to calculate distances for proper rare cell definition
        address = os.path.join(benchmark_path, f"{ref}/adata.h5ad")
        adata = sc.read_h5ad(address)
        
        # Calculate cell type distances
        cell_type_df = calculate_cell_type_distances(adata, label_key=label_key)
        
        # Identify rare cell types using distance AND frequency criteria
        total_cells = obs.shape[0]
        rare_types, _, _ = identify_rare_cell_types_distance_and_frequency(
            cell_type_df, 
            total_cells, 
            distance_percentile=75, 
            frequency_threshold_pct=frequency_threshold_pct
        )
        
        logger.info(f"Ref: {ref}, Rare types identified: {rare_types}")
        
        # Store the sum count of rare types for this reference size
        rare_type_counts[ref] = real_counts[rare_types].sum() if len(rare_types) > 0 else 0
        obs['rare'] = obs[label_key].isin(rare_types)
        obs['real_count'] = obs[label_key].map(real_counts)
        
        group_dict = {key: [] for key in keys}
        
        for size, method, rep in product(sizes, METHODS, REPS):
            try:
                if method == 'sps_fi25':
                    # Use test_feature_index data for SPS with feature_index
                    sps_address = os.path.join(feature_index_path, f"{ref}/{feature_index}/{size}/{rep}/results.pkl")
                    if os.path.isfile(sps_address):
                        with open(sps_address, 'rb') as handle:
                            samples = pickle.load(handle)
                        # test_feature_index has nested structure: samples[0][0] contains indices
                        indices = obs.iloc[samples[0][0]].index.tolist()
                    else:
                        logger.warning(f"SPS feature_index file not found: {sps_address}")
                        continue
                elif method == 'atomic':
                    atomic_address = os.path.join(benchmark_path, f"{ref}/atomic/{size}/{rep}/results.csv")
                    if os.path.isfile(atomic_address):
                        indices = pd.read_csv(atomic_address)['x'].values.astype(str)
                    else:
                        continue
                else:
                    # Use benchmark data for other methods (random, hopper, scsampler)
                    samples_address = os.path.join(benchmark_path, f"{ref}/{method}/{size}/{rep}/results.pkl")
                    with open(samples_address, 'rb') as handle:
                        samples = pickle.load(handle)
                    indices = obs.iloc[samples[0]].index.tolist()
                
                value_counts = obs.loc[indices, label_key].value_counts()
                counts, labels = value_counts.tolist(), value_counts.index.tolist()
                method_rep, size_rep, ref_rep, rep_rep = ([x] * len(counts) for x in [method, size, ref, rep])
                real_count = real_counts[labels].values
                values = [labels, counts, method_rep, ref_rep, size_rep, rep_rep, real_count]
                for key, value in zip(keys, values):
                    group_dict[key].extend(value)
            except Exception as e:
                logger.warning(f"Error processing {method} for ref {ref}, size {size}, rep {rep}: {str(e)}")
                continue
        
        gdf = pd.DataFrame.from_dict(group_dict)
        res[ref] = gdf[gdf[label_key].isin(rare_types)]
    
    return res, rare_type_counts

# --- Plotting ---
def plot_single_rare_cell_coverage(ax, res, references, sample_size, dataset_label, colors, x_values=None):
    """Plot single rare cell coverage comparison."""
    output = {}
    for ref in references:
        df = res[ref]
        output[ref] = np.array([
            df[(df['method'] == method) & (df['sample_size'] == sample_size)]['count'].sum() / len(REPS)
            for method in METHODS
        ])
    if x_values is None:
        x_values = np.array(references)
    else:
        x_values = np.array(x_values)
    
    for j, method in enumerate(METHODS):
        y_values = np.array([output[ref][j] for ref in references])
        ax.plot(x_values, y_values, '-o', color=colors[j % len(colors)], label=METHOD_LABELS[j])
    
    ax.set_xlabel('Reference Size (M)')
    ax.set_ylabel('Number of Covered Rare Cells')
    ax.set_title(dataset_label)
    ax.grid(True, which="both", linestyle="--")


def plot_rare_cell_coverage_subplots(res, references, sizes, dataset_label, filename, colors, n_rows, n_cols, rare_type_counts, x_values=None):
    """Plot rare cell coverage subplots for different sample sizes."""
    output2 = {}
    for size in sizes:
        b = []
        for ref in references:
            df = res[ref]
            a = []
            for method in METHODS:
                sample_count = df[(df['ref_size'] == ref) & (df['method'] == method) & (df['sample_size'] == size)]['count'].sum()
                sample_count /= len(REPS)
                a.append(sample_count)
            b.append(a)
        b = np.array(b).T
        output2[size] = b
    
    plt.rcParams.update({
        'font.size': 40,
        'axes.labelsize': 40,
        'axes.titlesize': 44,
        'legend.fontsize': 36,
        'xtick.labelsize': 36,
        'ytick.labelsize': 36,
        'figure.figsize': [36, 30],
        'legend.markerscale': 2.0,
        'lines.markersize': 10,
        'lines.linewidth': 2.5
    })
    
    if x_values is None:
        x_values = np.array(references)
    else:
        x_values = np.array(x_values)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(36, 30))
    axes = axes.flatten()
    
    for i, size in enumerate(output2.keys()):
        y_streams = output2[size]
        ax = axes[i]
        ax2 = ax.twinx()
        
        # Plot lines on primary y-axis
        for j, y_stream in enumerate(y_streams):
            ax.plot(x_values, y_stream, '-o', color=colors[j % len(colors)], label=METHOD_LABELS[j])
        
        # Set bar width
        bar_width = 0.15
        
        # Plot bars on secondary y-axis
        ax2.bar(x_values, [rare_type_counts[ref] for ref in references], 
                width=bar_width, alpha=0.3, color='gray', label='Total Rare Cells', align='center')
        
        # Set labels and titles
        ax.set_ylabel('Number of Rare Cells')
        ax2.set_ylabel('Total Rare Cells', color='gray')
        ax2.tick_params(axis='y', labelcolor='gray')
        
        if len(sizes) == 4:
            if i > 1:
                ax.set_xlabel('Reference Size (M)')
        if len(sizes) == 3:
            if i > 0:
                ax.set_xlabel('Reference Size (M)')
        
        ax.set_title(f'Sample Size: {size}')
        ax.grid(True, which="both", ls="--", alpha=0.7)
        
        # Set x-ticks
        ax.set_xticks(x_values)
        ax.set_xticklabels([str(x) for x in x_values])

    for i in range(len(output2), n_rows * n_cols):
        fig.delaxes(axes[i])

    # Get handles/labels from the first axis for the methods
    handles, labels = axes[1].get_legend_handles_labels()
    # Add a custom patch for the bar
    bar_patch = mpatches.Patch(color='gray', alpha=0.3, label='Total Rare Cells')
    handles.append(bar_patch)
    labels.append('Total Rare Cells')
    fig.legend(
        handles, labels,
        loc='center right',
        bbox_to_anchor=(1.08, 0.575),
        borderaxespad=0.0,
        frameon=True
    )
    plt.subplots_adjust(right=0.85, bottom=0.27, wspace=0.6)

    output_path = os.path.join(FIGURE_DIR, filename)
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved subplot figure: {output_path}")

# --- Main ---
def main():
    logger.info(f"Starting SPS feature_index={FEATURE_INDEX} coverage analysis for MCC dataset")
    colors = ['purple', 'blue', 'orange', 'green', 'grey']
    
    # Load reference data from benchmark path
    logger.info("Loading MCC reference data...")
    mcc_reference_obs = load_reference_data(MCC_REFERENCES, MCC_BENCHMARK_PATH)
    
    # Calculate rare cell coverage using SPS with feature_index and other methods
    logger.info("Calculating rare cell coverage...")
    mcc_res, mcc_rare_type_counts = calculate_rare_cell_coverage(
        mcc_reference_obs, 
        MCC_REFERENCES, 
        MCC_SIZES, 
        MCC_BENCHMARK_PATH,
        MCC_FEATURE_INDEX_PATH,
        FEATURE_INDEX
    )
    
    # X values for MCC (in millions of cells)
    MCC_X_VALUES = [0.5, 1.0, 2.0, 2.5, 3.2]
    
    # --- Plot single figure ---
    logger.info("Creating single plot figure...")
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    plot_single_rare_cell_coverage(ax, mcc_res, MCC_REFERENCES, 100000, 
                                    f'MCC - SPS (feature_index={FEATURE_INDEX})', colors, x_values=MCC_X_VALUES)
    ax.legend(loc='best')
    plt.tight_layout()
    output_path = os.path.join(FIGURE_DIR, f'sps_fi{FEATURE_INDEX}_mcc_coverage_single.jpg')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved single figure: {output_path}")
    
    # --- Plot subplots ---
    logger.info("Creating MCC subplots...")
    plot_rare_cell_coverage_subplots(
        mcc_res, 
        MCC_REFERENCES, 
        MCC_SIZES, 
        'MCC', 
        f'sps_fi{FEATURE_INDEX}_mcc_coverage.jpg', 
        colors, 
        2, 2, 
        mcc_rare_type_counts, 
        x_values=MCC_X_VALUES
    )
    
    logger.info("All figures created successfully.")

if __name__ == "__main__":
    main()
