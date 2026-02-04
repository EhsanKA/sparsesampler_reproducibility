#!/usr/bin/env python
# coding: utf-8

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
    log_file = os.path.join(log_dir, f'rare_cell_coverage_combined_revision_{timestamp}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

# Constants
MCC_01_REFERENCES = [5, 10, 20, 25, 30]
MCC_05_REFERENCES = [5, 10, 20, 25, 30]
METHODS = ['random', 'sps', 'hopper', 'atomic', 'scsampler']
MCC_01_SIZES = [50000, 100000, 200000, 300000]
MCC_05_SIZES = [50000, 100000, 200000, 300000]
REPS = [0]  # Using only 1 rep (rep 0) instead of 5 reps
label_key = 'celltype'

# File paths
def get_data_path():
    """Get data path from environment variable or derive from script location."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        # Derive project root from script location
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))
    return os.path.join(project_root, 'data')

def get_figures_dir():
    """Get figures directory path from environment variable or derive from script location."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        # Derive project root from script location
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))
    return os.path.join(project_root, 'figures', 'revision')

file_path_env = get_data_path()
MCC_01_PATH = os.path.join(file_path_env, 'mcc_01/benchmark')
MCC_05_PATH = os.path.join(file_path_env, 'mcc_05/benchmark')
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
    reference_obs = {}
    for ref in references:
        address = os.path.join(path, f"{ref}/obs.csv")
        df = pd.read_csv(address, index_col=0)
        df[label_key] = df[label_key].astype('category')
        reference_obs[ref] = df
    return reference_obs

def calculate_rare_cell_coverage(reference_obs, references, sizes, path, dataset_type):
    res = {}
    rare_type_counts = {}  # Dictionary to store sum counts of rare types for each reference size
    keys = [label_key, 'count', 'method', 'ref_size', 'sample_size', 'rep', 'real_count']
    
    # Define frequency threshold percentage based on dataset
    if dataset_type == 'mcc_01':
        frequency_threshold_pct = 0.1  # 1/1000 = 0.1%
    elif dataset_type == 'mcc_05':
        frequency_threshold_pct = 0.5  # 5/1000 = 0.5%
    else:
        frequency_threshold_pct = 1.0  # Default 1%
    
    for ref in references:
        obs = reference_obs[ref]
        real_counts = obs[label_key].value_counts()
        
        # Load adata to calculate distances
        address = os.path.join(path, f"{ref}/adata.h5ad")
        adata = sc.read_h5ad(address)
        
        # Calculate cell type distances
        cell_type_df = calculate_cell_type_distances(adata, label_key=label_key)
        
        # Identify rare cell types using distance and frequency
        total_cells = obs.shape[0]
        rare_types, _, _ = identify_rare_cell_types_distance_and_frequency(
            cell_type_df, 
            total_cells, 
            distance_percentile=75, 
            frequency_threshold_pct=frequency_threshold_pct
        )
        
        # Store the sum count of rare types for this reference size
        rare_type_counts[ref] = real_counts[rare_types].sum() if len(rare_types) > 0 else 0
        obs['rare'] = obs[label_key].isin(rare_types)
        obs['real_count'] = obs[label_key].map(real_counts)
        label_order = obs[label_key].value_counts().index.values
        group_dict = {key: [] for key in keys}
        for size, method, rep in product(sizes, METHODS, REPS):
            try:
                if method == 'atomic':
                    atomic_address = os.path.join(path, f"{ref}/atomic/{size}/{rep}/results.csv")
                    if os.path.isfile(atomic_address):
                        if dataset_type in ['mcc_01', 'mcc_05']:
                            indices = pd.read_csv(atomic_address)['x'].values.astype(str)
                        else:
                            indices = pd.read_csv(atomic_address)['x'].values.astype(int)
                    else:
                        continue
                else:
                    samples_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/results.pkl")
                    with open(samples_address, 'rb') as handle:
                        samples = pickle.load(handle)
                    # Handle different pickle structures: SPS has nested tuple (inner_tuple, time)
                    # where inner_tuple is (indices_list, ...), other methods have (indices, time)
                    sample_indices = samples[0]
                    if isinstance(sample_indices, tuple):
                        sample_indices = sample_indices[0]
                    indices = obs.iloc[sample_indices].index.tolist()
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
        ax.plot(x_values, y_values, '-o', color=colors[j % len(colors)], label=method)
    ax.set_xlabel('Reference Size (M)')
    ax.set_ylabel('Number of Covered Rare Cells')
    ax.set_title(dataset_label)
    ax.grid(True, which="both", linestyle="--")


def plot_rare_cell_coverage_subplots(res, references, sizes, dataset_label, filename, colors, n_rows, n_cols, rare_type_counts, x_values=None):
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
        ax2 = ax.twinx()  # Create a second y-axis
        
        # Plot lines on primary y-axis
        for j, y_stream in enumerate(y_streams):
            ax.plot(x_values, y_stream, '-o', color=colors[j % len(colors)], label=METHODS[j])
        
        # Set bar width based on dataset_label
        if dataset_label == 'LCMV':
            bar_width = 2.0
        else:
            bar_width = 0.15
        # Plot bars on secondary y-axis at the actual x_values
        ax2.bar(x_values, [rare_type_counts[ref] for ref in references], 
                width=bar_width, alpha=0.3, color='gray', label='Total Rare Cells', align='center')
        
        # Set labels and titles
        # if i % 2 == 0:
        ax.set_ylabel('Number of Rare Cells')
        ax2.set_ylabel('Total Rare Cells', color='gray')
        ax2.tick_params(axis='y', labelcolor='gray')
        
        if len(sizes) == 4:
            if i>1:
                ax.set_xlabel('Reference Size (M)')
        if len(sizes) == 3:
            if i>0:
                ax.set_xlabel('Reference Size (M)')
        
        ax.set_title(f'Sample Size: {size}')
        ax.grid(True, which="both", ls="--", alpha=0.7)
        
        # Set x-ticks to match the x_values and format them nicely
        ax.set_xticks(x_values)
        ax.set_xticklabels([str(x) for x in x_values])

    for i in range(len(output2), n_rows * n_cols):
        fig.delaxes(axes[i])

    # Only get handles/labels from the first axis for the methods
    handles, labels = axes[1].get_legend_handles_labels()
    # Add a custom patch for the bar
    bar_patch = mpatches.Patch(color='gray', alpha=0.3, label='Total Rare Cells')
    handles.append(bar_patch)
    labels.append('Total Rare Cells')
    fig.legend(
        handles, labels,
        loc='center right',
        bbox_to_anchor=(1.08, 0.575),  # Move legend slightly closer
        borderaxespad=0.0,
        frameon=True
    )
    plt.subplots_adjust(right=0.85, bottom=0.27, wspace=0.6)  # Adjust right margin

    output_path = os.path.join(FIGURE_DIR, filename)
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved subplot figure: {output_path}")

# --- Main ---
def main():
    logger.info("Starting combined rare cell coverage analysis for revision (MCC_01 and MCC_05)")
    colors = ['purple', 'blue', 'orange', 'green', 'grey']
    # MCC_01
    logger.info("Processing MCC_01 data...")
    mcc_01_reference_obs = load_reference_data(MCC_01_REFERENCES, MCC_01_PATH)
    mcc_01_res, mcc_01_rare_type_counts = calculate_rare_cell_coverage(mcc_01_reference_obs, MCC_01_REFERENCES, MCC_01_SIZES, MCC_01_PATH, 'mcc_01')
    # MCC_05
    logger.info("Processing MCC_05 data...")
    mcc_05_reference_obs = load_reference_data(MCC_05_REFERENCES, MCC_05_PATH)
    mcc_05_res, mcc_05_rare_type_counts = calculate_rare_cell_coverage(mcc_05_reference_obs, MCC_05_REFERENCES, MCC_05_SIZES, MCC_05_PATH, 'mcc_05')
    # --- Plot both single plots side by side ---
    logger.info("Creating combined single plot figure...")
    MCC_X_VALUES = [0.5, 1.0, 2.0, 2.5, 3.2]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
    plot_single_rare_cell_coverage(ax1, mcc_01_res, MCC_01_REFERENCES, 100000, 'MCC_01', colors, x_values=MCC_X_VALUES)
    plot_single_rare_cell_coverage(ax2, mcc_05_res, MCC_05_REFERENCES, 100000, 'MCC_05', colors, x_values=MCC_X_VALUES)
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=5)
    plt.tight_layout(rect=[0, 0.13, 1, 1])
    output_path = os.path.join(FIGURE_DIR, 'figure1_combined_coverage.jpg')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved combined single plot figure: {output_path}")
    # --- MCC_01 subplots ---
    logger.info("Creating MCC_01 subplots...")
    plot_rare_cell_coverage_subplots(mcc_01_res, MCC_01_REFERENCES, MCC_01_SIZES, 'MCC_01', 'supp_figure3_mcc_01_coverage.jpg', colors, 2, 2, mcc_01_rare_type_counts, x_values=MCC_X_VALUES)
    # --- MCC_05 subplots ---
    logger.info("Creating MCC_05 subplots...")
    plot_rare_cell_coverage_subplots(mcc_05_res, MCC_05_REFERENCES, MCC_05_SIZES, 'MCC_05', 'supp_figure3_mcc_05_coverage.jpg', colors, 2, 2, mcc_05_rare_type_counts, x_values=MCC_X_VALUES)
    logger.info("All figures created successfully.")

if __name__ == "__main__":
    main()

