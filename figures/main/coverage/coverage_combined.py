#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import logging
from datetime import datetime
from itertools import product
import matplotlib.patches as mpatches

# Set up logging

def setup_logger():
    log_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'rare_cell_coverage_combined_{timestamp}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

# Constants
LCMV_REFERENCES = [1, 5, 10, 20, 34]
MCC_REFERENCES = [5, 10, 20, 25, 30]
METHODS = ['random', 'sps', 'hopper', 'atomic', 'scsampler']
LCMV_SIZES = [50000, 100000, 200000]
MCC_SIZES = [50000, 100000, 200000, 300000]
REPS = [i for i in range(5)]
label_key = 'celltype'

# File paths
file_path_env = 'projects/sparsesampler_reproducibility/data'
LCMV_PATH = os.path.join(file_path_env, 'lcmv/benchmark')
MCC_PATH = os.path.join(file_path_env, 'mcc/benchmark')
FIGURE_DIR = 'projects/sparsesampler_reproducibility/figures'

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
    for ref in references:
        obs = reference_obs[ref]
        real_counts = obs[label_key].value_counts()
        rare_types = real_counts[real_counts < obs.shape[0] / 100].index.tolist()
        # Store the sum count of rare types for this reference size
        rare_type_counts[ref] = real_counts[rare_types].sum()
        obs['rare'] = obs[label_key].isin(rare_types)
        obs['real_count'] = obs[label_key].map(real_counts)
        label_order = obs[label_key].value_counts().index.values
        group_dict = {key: [] for key in keys}
        for size, method, rep in product(sizes, METHODS, REPS):
            try:
                if method == 'atomic':
                    atomic_address = os.path.join(path, f"{ref}/atomic/{size}/{rep}/results.csv")
                    if os.path.isfile(atomic_address):
                        if dataset_type == 'mcc':
                            indices = pd.read_csv(atomic_address)['x'].values.astype(str)
                        else:
                            indices = pd.read_csv(atomic_address)['x'].values.astype(int)
                    else:
                        continue
                else:
                    samples_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/results.pkl")
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
    logger.info("Starting combined rare cell coverage analysis")
    colors = ['purple', 'blue', 'orange', 'green', 'grey']
    # LCMV
    logger.info("Processing LCMV data...")
    lcmv_reference_obs = load_reference_data(LCMV_REFERENCES, LCMV_PATH)
    lcmv_res, lcmv_rare_type_counts = calculate_rare_cell_coverage(lcmv_reference_obs, LCMV_REFERENCES, LCMV_SIZES, LCMV_PATH, 'lcmv')
    # MCC
    logger.info("Processing MCC data...")
    mcc_reference_obs = load_reference_data(MCC_REFERENCES, MCC_PATH)
    mcc_res, mcc_rare_type_counts = calculate_rare_cell_coverage(mcc_reference_obs, MCC_REFERENCES, MCC_SIZES, MCC_PATH, 'mcc')
    # --- Plot both single plots side by side ---
    logger.info("Creating combined single plot figure...")
    MCC_X_VALUES = [0.5, 1.0, 2.0, 2.5, 3.2]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
    plot_single_rare_cell_coverage(ax1, mcc_res, MCC_REFERENCES, 100000, 'MCC', colors, x_values=MCC_X_VALUES)
    plot_single_rare_cell_coverage(ax2, lcmv_res, LCMV_REFERENCES, 100000, 'LCMV', colors)
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=5)
    plt.tight_layout(rect=[0, 0.13, 1, 1])
    output_path = os.path.join(FIGURE_DIR, 'figure1_combined_coverage.jpg')
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved combined single plot figure: {output_path}")
    # --- LCMV subplots ---
    logger.info("Creating LCMV subplots...")
    plot_rare_cell_coverage_subplots(lcmv_res, LCMV_REFERENCES, LCMV_SIZES, 'LCMV', 'supp_figure4_lcmv_coverage.jpg', colors, 2, 2, lcmv_rare_type_counts)
    # --- MCC subplots ---
    logger.info("Creating MCC subplots...")
    plot_rare_cell_coverage_subplots(mcc_res, MCC_REFERENCES, MCC_SIZES, 'MCC', 'supp_figure3_mcc_coverage.jpg', colors, 2, 2, mcc_rare_type_counts, x_values=MCC_X_VALUES)
    logger.info("All figures created successfully.")

if __name__ == "__main__":
    main() 