#!/usr/bin/env python
# coding: utf-8

"""
Calculate the distance between the mean of each cell type and the mean of the whole dataset
for mcc_01, mcc_05, mcc, and lcmv datasets.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.spatial.distance import euclidean

# Get data path from environment variable or derive from script location
def get_data_path():
    """Get data path from environment variable or derive from script location."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
    return os.path.join(project_root, 'data')

def get_sparseflow_data_path():
    """Get data path for sparseFlow_benchmarking directory."""
    return '/fast/AG_Ohler/ekarimi/projects/sparseFlow_benchmarking/data'

def load_reference_data(dataset, ref, base_path=None):
    """Load reference data for a given dataset and reference size.
    
    Parameters:
    -----------
    dataset : str
        Dataset name (e.g., 'mcc_01', 'mcc_05', 'mcc', 'lcmv')
    ref : int
        Reference number
    base_path : str, optional
        Base path to data directory. If None, will use appropriate path based on dataset.
    """
    if base_path is None:
        # Determine which base path to use
        if dataset in ['mcc_01', 'mcc_05']:
            data_path = get_data_path()
        elif dataset in ['mcc', 'lcmv']:
            data_path = get_sparseflow_data_path()
        else:
            raise ValueError(f"Unknown dataset: {dataset}")
    else:
        data_path = base_path
    
    directory = f"{dataset}/benchmark"
    path = os.path.join(data_path, directory)
    address = os.path.join(path, f"{ref}/adata.h5ad")
    adata = sc.read_h5ad(address)
    return adata

def calculate_cell_type_distances(adata, label_key='celltype'):
    """
    Calculate the distance between the mean of each cell type and the mean of the whole dataset.
    
    Parameters:
    -----------
    adata : sc.AnnData
        Annotated data object
    label_key : str
        Key in adata.obs that contains cell type labels
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: cell_type, distance, n_cells
    """
    # Get expression matrix
    X = adata.X
    
    # Convert sparse matrix to dense if needed
    if scipy.sparse.issparse(X):
        X = X.toarray()
    
    # Calculate mean of the whole dataset
    dataset_mean = np.mean(X, axis=0)
    
    # Get cell type labels
    cell_types = adata.obs[label_key].values
    
    # Calculate mean for each cell type and distance to dataset mean
    results = []
    unique_cell_types = np.unique(cell_types)
    
    for cell_type in unique_cell_types:
        # Get indices for this cell type
        cell_type_mask = cell_types == cell_type
        cell_type_data = X[cell_type_mask, :]
        
        # Calculate mean of this cell type
        cell_type_mean = np.mean(cell_type_data, axis=0)
        
        # Calculate Euclidean distance between cell type mean and dataset mean
        distance = euclidean(cell_type_mean, dataset_mean)
        
        # Count number of cells
        n_cells = np.sum(cell_type_mask)
        
        results.append({
            'cell_type': cell_type,
            'distance': distance,
            'n_cells': n_cells
        })
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Sort by distance (descending - most distant first)
    df = df.sort_values('distance', ascending=False).reset_index(drop=True)
    
    return df

def main():
    """Main function to analyze mcc_01, mcc_05, mcc, and lcmv datasets."""
    # Define datasets with their reference numbers
    # Reference numbers correspond to the full dataset
    dataset_configs = [
        {'name': 'mcc_01', 'ref': 30},
        {'name': 'mcc_05', 'ref': 30},
        {'name': 'mcc', 'ref': 30},
        {'name': 'lcmv', 'ref': 34},
    ]
    
    label_key = 'celltype'
    
    # Get output directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, 'results')
    os.makedirs(output_dir, exist_ok=True)
    
    all_results = []
    
    for config in dataset_configs:
        dataset = config['name']
        ref = config['ref']
        
        print(f"\n{'='*60}")
        print(f"Processing {dataset} (reference {ref})")
        print(f"{'='*60}")
        
        try:
            # Load data
            print(f"Loading data from {dataset}/benchmark/{ref}/adata.h5ad...")
            adata = load_reference_data(dataset, ref)
            print(f"Loaded {adata.shape[0]} cells and {adata.shape[1]} genes")
            
            # Calculate distances
            print("Calculating cell type distances...")
            results_df = calculate_cell_type_distances(adata, label_key=label_key)
            
            # Add dataset column
            results_df['dataset'] = dataset
            
            # Save individual results
            output_file = os.path.join(output_dir, f'{dataset}_cell_type_distances.csv')
            results_df.to_csv(output_file, index=False)
            print(f"Results saved to: {output_file}")
            
            # Print summary
            print(f"\nSummary for {dataset}:")
            print(f"Number of cell types: {len(results_df)}")
            n_top = min(5, len(results_df))
            if n_top > 0:
                print(f"\nTop {n_top} most distant cell types:")
                print(results_df.head(n_top).to_string(index=False))
                print(f"\nTop {n_top} closest cell types:")
                print(results_df.tail(n_top).to_string(index=False))
            
            all_results.append(results_df)
        except Exception as e:
            print(f"Error processing {dataset}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    # Combine results from all datasets
    if len(all_results) > 0:
        combined_df = pd.concat(all_results, ignore_index=True)
        
        # Save combined results
        combined_output_file = os.path.join(output_dir, 'combined_cell_type_distances.csv')
        combined_df.to_csv(combined_output_file, index=False)
        print(f"\n{'='*60}")
        print(f"Combined results saved to: {combined_output_file}")
        print(f"{'='*60}")
        
        # Print combined summary
        print(f"\nCombined Summary:")
        print(f"Total datasets processed: {len(all_results)}")
        print(f"Total unique cell types across all datasets: {combined_df['cell_type'].nunique()}")
        print(f"\nOverall statistics:")
        print(f"  Mean distance: {combined_df['distance'].mean():.4f}")
        print(f"  Median distance: {combined_df['distance'].median():.4f}")
        print(f"  Min distance: {combined_df['distance'].min():.4f}")
        print(f"  Max distance: {combined_df['distance'].max():.4f}")
    else:
        print("\nNo datasets were successfully processed.")

if __name__ == "__main__":
    main()

