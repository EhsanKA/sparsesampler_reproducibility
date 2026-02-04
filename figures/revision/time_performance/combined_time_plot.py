#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import time
import pickle
import matplotlib.pyplot as plt
import matplotlib
import logging
from datetime import datetime

# Set up logging
def setup_logger():
    """Set up and configure the logger."""
    # Create logs directory if it doesn't exist
    log_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 'logs')
    os.makedirs(log_dir, exist_ok=True)
    
    # Create a log file with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'combined_time_plot_revision_{timestamp}.log')
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

# Initialize logger
logger = setup_logger()

# Constants
MCC_01_REFERENCES = [5, 10, 20, 25, 30]
MCC_05_REFERENCES = [5, 10, 20, 25, 30]
METHODS = ['sps', 'hopper', 'atomic', 'scsampler']
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
MCC_01_DIR = "mcc_01/benchmark"
MCC_05_DIR = "mcc_05/benchmark"
MCC_01_PATH = os.path.join(file_path_env, MCC_01_DIR)
MCC_05_PATH = os.path.join(file_path_env, MCC_05_DIR)

# Output directory for figures
FIGURE_DIR = get_figures_dir()

# Plotting configuration
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

def prepare_times_sample_size_fixed(dataset_type='mcc_01'):
    """Calculate average runtime for each method across different reference sizes."""
    logger.info(f"Calculating runtime statistics for {dataset_type}...")
    times = {}
    references = MCC_01_REFERENCES if dataset_type == 'mcc_01' else MCC_05_REFERENCES
    sizes = MCC_01_SIZES if dataset_type == 'mcc_01' else MCC_05_SIZES
    path = MCC_01_PATH if dataset_type == 'mcc_01' else MCC_05_PATH
    
    for size in sizes:
        logger.info(f"Processing sample size: {size}")
        time_list_per_ref_per_method = []
        for ref in references:
            logger.debug(f"Processing reference size: {ref}")
            time_list_per_method = []
            for method in METHODS:
                time = 0
                valid_reps = 0
                for rep in REPS:
                    try:
                        if method == 'atomic':
                            atomic_address = os.path.join(path, f"{ref}/atomic/{size}/{rep}/runtimes.csv")
                            if os.path.isfile(atomic_address):
                                runtime = pd.read_csv(atomic_address)['runtime'][0]
                                time += runtime
                                valid_reps += 1
                        else:
                            samples_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/results.pkl")
                            with open(samples_address, 'rb') as handle:
                                samples = pickle.load(handle)
                            time += samples[1]
                            valid_reps += 1
                    except Exception as e:
                        logger.error(f"Error processing {method} for ref {ref}, size {size}, rep {rep}: {e}")
                        continue
                
                if valid_reps > 0:
                    time /= valid_reps
                time_list_per_method.append(time)
            
            time_list_per_ref_per_method.append(time_list_per_method)
        
        time_list_per_ref_per_method = np.array(time_list_per_ref_per_method).T
        times[size] = time_list_per_ref_per_method
    logger.info(f"Runtime statistics calculation completed for {dataset_type}")
    return times

def plot_single_runtime_comparison(times, sample_size=100000, dataset_type='mcc_01'):
    """Create a single plot comparing runtime across reference sizes for a specific sample size."""
    logger.info(f"Creating single runtime comparison plot for {dataset_type} with sample size {sample_size}")
    references = MCC_01_REFERENCES if dataset_type == 'mcc_01' else MCC_05_REFERENCES
    x_values = np.array(references).astype(float)
    if dataset_type in ['mcc_01', 'mcc_05']:
        x_values = x_values / 10
        x_values[0] = 0.5
        x_values[-1] = 3.2
    
    colors = ['blue', 'orange', 'green', 'grey']
    y_streams = times[sample_size]
    
    fig, ax = plt.subplots(figsize=(9, 6))
    
    for j, y_stream in enumerate(y_streams):
        ax.plot(x_values, y_stream, '-o', color=colors[j % len(colors)], label=METHODS[j])
    
    ax.set_xlabel('Reference Size (M)')
    ax.set_ylabel('Time (s)')
    ax.set_title(f'Runtime Comparison for {dataset_type.upper()} (Sample Size: {sample_size})')
    ax.set_yscale('log')
    ax.grid(True, which="both", linestyle="--")
    
    # Add legend
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=4)
    
    # Save figure with new name
    output_path = os.path.join(FIGURE_DIR, f'figure1_{dataset_type}_time.jpg')
    logger.info(f"Saving single plot to {output_path}")
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def plot_combined_runtime_comparison(mcc_01_times, mcc_05_times, sample_size=100000):
    """Create a combined plot with both MCC_01 and MCC_05 runtime comparisons."""
    logger.info("Creating combined runtime comparison plot")
    
    # Set up the figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
    
    # Plot MCC_01 first
    mcc_01_x_values = np.array(MCC_01_REFERENCES).astype(float) / 10
    mcc_01_x_values[0] = 0.5
    mcc_01_x_values[-1] = 3.2
    mcc_01_y_streams = mcc_01_times[sample_size]
    colors = ['blue', 'orange', 'green', 'grey']
    
    for j, y_stream in enumerate(mcc_01_y_streams):
        ax1.plot(mcc_01_x_values, y_stream, '-o', color=colors[j % len(colors)])
    
    ax1.set_xlabel('Reference Size (M)')
    ax1.set_ylabel('Time (s)')
    ax1.set_title('MCC_01')
    ax1.set_yscale('log')
    ax1.grid(True, which="both", linestyle="--")
    
    # Plot MCC_05 second
    mcc_05_x_values = np.array(MCC_05_REFERENCES).astype(float) / 10
    mcc_05_x_values[0] = 0.5
    mcc_05_x_values[-1] = 3.2
    mcc_05_y_streams = mcc_05_times[sample_size]
    
    for j, y_stream in enumerate(mcc_05_y_streams):
        ax2.plot(mcc_05_x_values, y_stream, '-o', color=colors[j % len(colors)])
    
    ax2.set_xlabel('Reference Size (M)')
    ax2.set_ylabel('Time (s)')
    ax2.set_title('MCC_05')
    ax2.set_yscale('log')
    ax2.grid(True, which="both", linestyle="--")
    
    # Add legend
    handles, labels = ax1.get_legend_handles_labels()
    if not handles:
        # Create handles from the first plot
        handles = [plt.Line2D([0], [0], marker='o', color=colors[j], linestyle='-') for j in range(len(METHODS))]
        labels = METHODS
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.05), ncol=4)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.1, 1, 1])
    
    # Save figure
    output_path = os.path.join(FIGURE_DIR, 'figure1_combined_time.jpg')
    logger.info(f"Saving combined plot to {output_path}")
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

def plot_runtime_subplots(times, dataset_type='mcc_01'):
    """Create subplots comparing runtime across different sample sizes."""
    logger.info(f"Creating runtime comparison subplots for {dataset_type}")
    references = MCC_01_REFERENCES if dataset_type == 'mcc_01' else MCC_05_REFERENCES
    sizes = MCC_01_SIZES if dataset_type == 'mcc_01' else MCC_05_SIZES
    x_values = np.array(references).astype(float)
    if dataset_type in ['mcc_01', 'mcc_05']:
        x_values = x_values / 10
        x_values[0] = 0.5
        x_values[-1] = 3.2
    
    n_rows = len(sizes) // 2 + len(sizes) % 2
    n_cols = 2
    colors = ['blue', 'orange', 'green', 'grey']
    
    # Set larger font sizes specifically for subplots
    plt.rcParams.update({
        'font.size': 24,
        'axes.labelsize': 24,
        'axes.titlesize': 26,
        'legend.fontsize': 22,
        'xtick.labelsize': 22,
        'ytick.labelsize': 22,
        'figure.figsize': [20, 25],  # Increased figure size
        'legend.markerscale': 2.0,   # Increased marker size
        'lines.markersize': 10,      # Increased line marker size
        'lines.linewidth': 2.5       # Increased line width
    })
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 25))
    axes = axes.flatten()
    
    for i, size in enumerate(times.keys()):
        logger.debug(f"Creating subplot for sample size {size}")
        y_streams = times[size]
        ax = axes[i]
        
        for j, y_stream in enumerate(y_streams):
            ax.plot(x_values, y_stream, '-o', color=colors[j % len(colors)], label=METHODS[j])
        
        if i % 2 == 0:
            ax.set_ylabel('Time (s)')
        if len(sizes) == 4:
            if i>1:
                ax.set_xlabel('Reference Size (M)')
        if len(sizes) == 3:
            if i>0:
                ax.set_xlabel('Reference Size (M)')
        
        ax.set_title(f'Sample size: {size}')
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="--", alpha=0.7)  # Made grid lines slightly more visible
    
    # Remove unused axes
    for i in range(len(times), n_rows * n_cols):
        fig.delaxes(axes[i])
    
    # Use handles and labels from the second subplot for the legend
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc='center right',
        bbox_to_anchor=(1.0, 0.575),
        borderaxespad=0.0,
        frameon=True
    )
    plt.subplots_adjust(right=0.85, bottom=0.27)  # Make space for legend
    
    # Save figure with new name
    output_path = os.path.join(FIGURE_DIR, f'supp_figure5-6_{dataset_type}_time.jpg')
    logger.info(f"Saving subplots to {output_path}")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')  # Remove bbox_inches='tight'
    plt.close()
    
    # Reset font sizes back to default
    plt.rcParams.update({
        'font.size': 20,
        'axes.labelsize': 20,
        'axes.titlesize': 20,
        'legend.fontsize': 20,
        'xtick.labelsize': 20,
        'ytick.labelsize': 20,
        'figure.figsize': [6, 6],
        'legend.markerscale': 1.5,
        'lines.markersize': 6,
        'lines.linewidth': 1.5
    })

def main():
    logger.info("Starting combined time plot generation for revision (MCC_01 and MCC_05)")
    try:
        # Calculate times for both datasets
        mcc_01_times = prepare_times_sample_size_fixed('mcc_01')
        mcc_05_times = prepare_times_sample_size_fixed('mcc_05')
        
        # Generate individual plots
        plot_single_runtime_comparison(mcc_01_times, dataset_type='mcc_01')
        plot_single_runtime_comparison(mcc_05_times, dataset_type='mcc_05')
        
        # Generate combined plot
        plot_combined_runtime_comparison(mcc_01_times, mcc_05_times)
        
        # Generate subplots
        plot_runtime_subplots(mcc_01_times, dataset_type='mcc_01')
        plot_runtime_subplots(mcc_05_times, dataset_type='mcc_05')
        
        logger.info("Combined time plot generation completed successfully")
    except Exception as e:
        logger.error(f"An error occurred during execution: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    main()

