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
import scanpy as sc
import matplotlib.patches as mpatches

# Add analysis directory to path to import rare cell type functions
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
analysis_dir = os.path.join(project_root, 'analysis')
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

from refined_rare_cell_type_definition import (
    calculate_cell_type_distances,
    identify_rare_cell_types_distance_and_frequency
)

# Set up logging
def setup_logger():
    log_dir = os.path.join(project_root, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'plot_feature_index_{timestamp}.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    return logging.getLogger(__name__)

logger = setup_logger()

# Constants
DATASET_CONFIGS = {
    'lcmv': {
        'references': [1, 5, 10, 20, 34],
        'sizes': [50000, 100000, 200000],
        'frequency_threshold_pct': 1.0,
        'x_values': [0.1, 0.5, 1.0, 2.0, 3.4]  # For better visualization
    },
    'mcc': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'frequency_threshold_pct': 1.0,
        'x_values': [0.5, 1.0, 2.0, 2.5, 3.0]
    },
    'mcc_01': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'frequency_threshold_pct': 0.1,
        'x_values': [0.5, 1.0, 2.0, 2.5, 3.0]
    },
    'mcc_05': {
        'references': [5, 10, 20, 25, 30],
        'sizes': [50000, 100000, 200000, 300000],
        'frequency_threshold_pct': 0.5,
        'x_values': [0.5, 1.0, 2.0, 2.5, 3.0]
    }
}

FEATURE_INDICES = list(range(10, 31))  # 10 to 30
REPS = [0]
label_key = 'celltype'

# File paths
def get_data_path():
    """Get data path from environment variable or derive from script location."""
    project_root_env = os.environ.get('PROJECT_ROOT')
    if project_root_env is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root_env = os.path.dirname(os.path.dirname(script_dir))
    return os.path.join(project_root_env, 'data')

def get_figures_dir():
    """Get figures directory path."""
    project_root_env = os.environ.get('PROJECT_ROOT')
    if project_root_env is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root_env = os.path.dirname(os.path.dirname(script_dir))
    return os.path.join(project_root_env, 'figures', 'revision')

file_path_env = get_data_path()
FIGURE_DIR = get_figures_dir()
RESULTS_DIR = os.path.join(file_path_env, 'test_feature_index')

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


def get_benchmark_path(dataset):
    """Get benchmark path for a dataset."""
    # All datasets use the same data path
    return file_path_env


def load_reference_obs(dataset, ref):
    """Load reference observation data."""
    benchmark_path = get_benchmark_path(dataset)
    obs_path = os.path.join(benchmark_path, f"{dataset}/benchmark/{ref}/obs.csv")
    if os.path.exists(obs_path):
        df = pd.read_csv(obs_path, index_col=0)
        df[label_key] = df[label_key].astype('category')
        return df
    else:
        # If obs.csv doesn't exist, load from adata
        adata_path = os.path.join(benchmark_path, f"{dataset}/benchmark/{ref}/adata.h5ad")
        if os.path.exists(adata_path):
            adata = sc.read_h5ad(adata_path)
            return adata.obs
        else:
            logger.warning(f"Could not load obs for {dataset}, ref {ref}")
            return None


def identify_rare_cell_types(dataset, ref, obs):
    """Identify rare cell types for a dataset."""
    benchmark_path = get_benchmark_path(dataset)
    adata_path = os.path.join(benchmark_path, f"{dataset}/benchmark/{ref}/adata.h5ad")
    
    if not os.path.exists(adata_path):
        logger.warning(f"Could not find adata for {dataset}, ref {ref}")
        return []
    
    adata = sc.read_h5ad(adata_path)
    config = DATASET_CONFIGS[dataset]
    
    # Calculate cell type distances
    cell_type_df = calculate_cell_type_distances(adata, label_key=label_key)
    
    # Identify rare cell types using distance and frequency
    total_cells = obs.shape[0]
    rare_types, _, _ = identify_rare_cell_types_distance_and_frequency(
        cell_type_df, 
        total_cells, 
        distance_percentile=75, 
        frequency_threshold_pct=config['frequency_threshold_pct']
    )
    
    return rare_types


def calculate_rare_cell_coverage(dataset):
    """Calculate rare cell coverage for all feature_index values."""
    config = DATASET_CONFIGS[dataset]
    res = {}
    rare_type_counts = {}
    
    keys = [label_key, 'count', 'feature_index', 'ref_size', 'sample_size', 'rep', 'real_count']
    
    for ref in config['references']:
        obs = load_reference_obs(dataset, ref)
        if obs is None:
            continue
        
        real_counts = obs[label_key].value_counts()
        
        # Identify rare cell types
        rare_types = identify_rare_cell_types(dataset, ref, obs)
        rare_type_counts[ref] = real_counts[rare_types].sum() if len(rare_types) > 0 else 0
        
        obs['rare'] = obs[label_key].isin(rare_types)
        obs['real_count'] = obs[label_key].map(real_counts)
        
        group_dict = {key: [] for key in keys}
        
        for feature_index in FEATURE_INDICES:
            for size in config['sizes']:
                for rep in REPS:
                    try:
                        results_path = os.path.join(RESULTS_DIR, dataset, str(ref), 
                                                   str(feature_index), str(size), str(rep), 'results.pkl')
                        
                        if not os.path.exists(results_path):
                            logger.warning(f"Results not found: {results_path}")
                            continue
                        
                        with open(results_path, 'rb') as f:
                            results_data = pickle.load(f)
                        
                        # Results are saved as tuple: (indices, elapsed_time)
                        if isinstance(results_data, tuple):
                            indices = results_data[0]
                        elif isinstance(results_data, dict):
                            indices = results_data['indices']
                        else:
                            indices = results_data
                        
                        # Get cell indices (handle both list and array)
                        if isinstance(indices, tuple):
                            cell_indices = indices[0]
                        else:
                            cell_indices = indices
                        
                        # Map to obs index
                        if isinstance(cell_indices, (list, np.ndarray)):
                            if len(cell_indices) > 0:
                                if isinstance(cell_indices[0], (int, np.integer)):
                                    # Numeric indices
                                    sampled_obs = obs.iloc[cell_indices]
                                else:
                                    # String indices (cell names)
                                    sampled_obs = obs.loc[cell_indices]
                            else:
                                continue
                        else:
                            continue
                        
                        value_counts = sampled_obs[label_key].value_counts()
                        counts, labels = value_counts.tolist(), value_counts.index.tolist()
                        
                        feature_index_list = [feature_index] * len(counts)
                        size_list = [size] * len(counts)
                        ref_list = [ref] * len(counts)
                        rep_list = [rep] * len(counts)
                        real_count = real_counts[labels].values
                        
                        values = [labels, counts, feature_index_list, ref_list, size_list, rep_list, real_count]
                        for key, value in zip(keys, values):
                            group_dict[key].extend(value)
                            
                    except Exception as e:
                        logger.warning(f"Error processing {dataset}, ref {ref}, feature_index {feature_index}, size {size}, rep {rep}: {str(e)}")
                        continue
        
        gdf = pd.DataFrame.from_dict(group_dict)
        if len(gdf) > 0:
            res[ref] = gdf[gdf[label_key].isin(rare_types)]
        else:
            res[ref] = pd.DataFrame()
    
    return res, rare_type_counts


def plot_feature_index_comparison_subplots(res, references, sizes, dataset_label, filename, rare_type_counts, x_values, feature_indices_to_plot):
    """Plot comparison of different feature_index values across reference sizes.
    Similar structure to coverage_combined.py"""
    output2 = {}
    
    # Organize data: for each size, for each feature_index, get values for each reference
    for size in sizes:
        b = []
        for feature_index in feature_indices_to_plot:
            a = []
            for ref in references:
                if ref not in res or len(res[ref]) == 0:
                    a.append(0)
                    continue
                
                df = res[ref]
                sample_count = df[(df['ref_size'] == ref) & 
                                 (df['feature_index'] == feature_index) & 
                                 (df['sample_size'] == size)]['count'].sum()
                sample_count /= len(REPS)
                a.append(sample_count)
            b.append(a)
        # b is now: rows = feature_indices, columns = references
        b = np.array(b)
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
    
    n_sizes = len(sizes)
    n_cols = 2
    n_rows = (n_sizes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(36, 30))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    # Use a colormap for feature_index values
    colors = plt.cm.viridis(np.linspace(0, 1, len(feature_indices_to_plot)))
    
    for i, size in enumerate(sizes):
        if size not in output2:
            continue
        
        y_streams = output2[size]  # Shape: (num_feature_indices, num_references)
        ax = axes[i]
        ax2 = ax.twinx()  # Create a second y-axis
        
        # Plot lines on primary y-axis - one line per feature_index
        for j, feature_index in enumerate(feature_indices_to_plot):
            if j < y_streams.shape[0]:
                y_stream = y_streams[j, :]  # Get values for all references for this feature_index
                ax.plot(x_values, y_stream, '-o', color=colors[j], 
                       label=f'FI={feature_index}', linewidth=2.5, markersize=10)
        
        # Plot bars on secondary y-axis
        bar_width = 0.15 if dataset_label not in ['LCMV'] else 2.0
        ax2.bar(x_values, [rare_type_counts.get(ref, 0) for ref in references], 
                width=bar_width, alpha=0.3, color='gray', label='Total Rare Cells', align='center')
        
        ax.set_ylabel('Number of Rare Cells', fontsize=40)
        ax2.set_ylabel('Total Rare Cells', color='gray', fontsize=40)
        ax2.tick_params(axis='y', labelcolor='gray')
        
        if len(sizes) == 4:
            if i > 1:
                ax.set_xlabel('Reference Size (M)', fontsize=40)
        if len(sizes) == 3:
            if i > 0:
                ax.set_xlabel('Reference Size (M)', fontsize=40)
        
        ax.set_title(f'Sample Size: {size:,}', fontsize=44)
        ax.grid(True, which="both", ls="--", alpha=0.7)
        
        # Set x-ticks to match the x_values
        ax.set_xticks(x_values)
        ax.set_xticklabels([str(x) for x in references])
    
    # Hide unused subplots
    for i in range(len(sizes), len(axes)):
        fig.delaxes(axes[i])
    
    # Create legend
    handles, labels = axes[0].get_legend_handles_labels()
    bar_patch = mpatches.Patch(color='gray', alpha=0.3, label='Total Rare Cells')
    handles.append(bar_patch)
    labels.append('Total Rare Cells')
    fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.08, 0.575), 
              borderaxespad=0.0, frameon=True, fontsize=36)
    
    plt.subplots_adjust(right=0.85, bottom=0.27, wspace=0.6)
    
    output_path = os.path.join(FIGURE_DIR, filename)
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved subplot figure: {output_path}")


def main():
    logger.info("Starting feature_index comparison plotting")
    
    datasets = ['lcmv', 'mcc', 'mcc_01', 'mcc_05']
    
    # Use all feature indices for comparison (10-30)
    feature_indices_to_plot = FEATURE_INDICES
    
    for dataset in datasets:
        logger.info(f"Processing {dataset} data...")
        res, rare_type_counts = calculate_rare_cell_coverage(dataset)
        
        config = DATASET_CONFIGS[dataset]
        x_values = config.get('x_values', config['references'])
        
        # Create plot comparing all feature_index values for all sample sizes
        logger.info(f"Creating comparison plot for {dataset} with all feature_index values (comparing rare cells for each reference)")
        
        filename = f'coverage_{dataset}_feature_index_all.jpg'
        plot_feature_index_comparison_subplots(res, config['references'], config['sizes'], 
                                              dataset.upper(), filename, rare_type_counts, 
                                              x_values, feature_indices_to_plot)
        
        logger.info(f"Completed {dataset}")
    
    logger.info("All figures created successfully.")


if __name__ == "__main__":
    main()
