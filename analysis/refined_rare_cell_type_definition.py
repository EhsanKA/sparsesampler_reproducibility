#!/usr/bin/env python
# coding: utf-8

"""
Refined definitions for rare cell types using multiple criteria.

This script implements several approaches to define rare cell types:
1. Percentile-based (bottom 5th, 10th, 25th percentile by frequency)
2. Statistical (using z-scores or standard deviations)
3. Multi-criteria (combining low frequency AND high distance from mean)
4. Relative frequency (using median/mean frequency as threshold)
5. Distance-based (using distance percentiles)
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.spatial.distance import euclidean
from scipy import stats

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
    """Load reference data for a given dataset and reference size."""
    if base_path is None:
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
    """Calculate the distance between the mean of each cell type and the mean of the whole dataset."""
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    
    dataset_mean = np.mean(X, axis=0)
    cell_types = adata.obs[label_key].values
    results = []
    unique_cell_types = np.unique(cell_types)
    
    for cell_type in unique_cell_types:
        cell_type_mask = cell_types == cell_type
        cell_type_data = X[cell_type_mask, :]
        cell_type_mean = np.mean(cell_type_data, axis=0)
        distance = euclidean(cell_type_mean, dataset_mean)
        n_cells = np.sum(cell_type_mask)
        
        results.append({
            'cell_type': cell_type,
            'distance': distance,
            'n_cells': n_cells
        })
    
    df = pd.DataFrame(results)
    return df

def identify_rare_cell_types_percentile(cell_type_df, percentile=10):
    """
    Identify rare cell types using percentile-based approach.
    
    Parameters:
    -----------
    cell_type_df : pd.DataFrame
        DataFrame with cell_type, n_cells columns
    percentile : int
        Bottom percentile threshold (e.g., 10 = bottom 10%)
        
    Returns:
    --------
    list: Rare cell type names
    """
    threshold = np.percentile(cell_type_df['n_cells'], percentile)
    rare_types = cell_type_df[cell_type_df['n_cells'] < threshold]['cell_type'].tolist()
    return rare_types, threshold

def identify_rare_cell_types_zscore(cell_type_df, z_threshold=-1.5):
    """
    Identify rare cell types using z-score approach.
    
    Parameters:
    -----------
    cell_type_df : pd.DataFrame
        DataFrame with cell_type, n_cells columns
    z_threshold : float
        Z-score threshold (e.g., -1.5 means 1.5 standard deviations below mean)
        
    Returns:
    --------
    list: Rare cell type names
    """
    z_scores = stats.zscore(np.log1p(cell_type_df['n_cells']))  # log transform for better normality
    rare_mask = z_scores < z_threshold
    rare_types = cell_type_df[rare_mask]['cell_type'].tolist()
    return rare_types, z_threshold

def identify_rare_cell_types_relative_frequency(cell_type_df, relative_threshold=0.1):
    """
    Identify rare cell types using relative frequency (median-based).
    
    Parameters:
    -----------
    cell_type_df : pd.DataFrame
        DataFrame with cell_type, n_cells columns
    relative_threshold : float
        Fraction of median frequency (e.g., 0.1 = 10% of median)
        
    Returns:
    --------
    list: Rare cell type names
    """
    median_freq = cell_type_df['n_cells'].median()
    threshold = median_freq * relative_threshold
    rare_types = cell_type_df[cell_type_df['n_cells'] < threshold]['cell_type'].tolist()
    return rare_types, threshold

def identify_rare_cell_types_multi_criteria(cell_type_df, 
                                            freq_percentile=25, 
                                            distance_percentile=75):
    """
    Identify rare cell types using multi-criteria approach.
    Rare = low frequency AND high distance from mean.
    
    Parameters:
    -----------
    cell_type_df : pd.DataFrame
        DataFrame with cell_type, n_cells, distance columns
    freq_percentile : int
        Bottom percentile for frequency (e.g., 25 = bottom 25%)
    distance_percentile : int
        Top percentile for distance (e.g., 75 = top 25%)
        
    Returns:
    --------
    list: Rare cell type names
    """
    freq_threshold = np.percentile(cell_type_df['n_cells'], freq_percentile)
    distance_threshold = np.percentile(cell_type_df['distance'], distance_percentile)
    
    rare_mask = (cell_type_df['n_cells'] < freq_threshold) & \
                (cell_type_df['distance'] > distance_threshold)
    rare_types = cell_type_df[rare_mask]['cell_type'].tolist()
    return rare_types, freq_threshold, distance_threshold

def identify_rare_cell_types_distance_based(cell_type_df, distance_percentile=75):
    """
    Identify rare cell types using distance from mean.
    Assumes rare cell types are those far from the dataset mean.
    
    Parameters:
    -----------
    cell_type_df : pd.DataFrame
        DataFrame with cell_type, distance columns
    distance_percentile : int
        Top percentile for distance (e.g., 75 = top 25%)
        
    Returns:
    --------
    list: Rare cell type names
    """
    threshold = np.percentile(cell_type_df['distance'], distance_percentile)
    rare_types = cell_type_df[cell_type_df['distance'] > threshold]['cell_type'].tolist()
    return rare_types, threshold

def identify_rare_cell_types_distance_and_frequency(cell_type_df, total_cells, 
                                                     distance_percentile=75, 
                                                     frequency_threshold_pct=1.0):
    """
    Identify rare cell types using combined criteria:
    - Top percentile of distances (e.g., top 25%)
    - AND frequency less than specified percentage of total dataset (e.g., < 1%)
    
    This combines distance-based and absolute frequency thresholds.
    
    Parameters:
    -----------
    cell_type_df : pd.DataFrame
        DataFrame with cell_type, n_cells, distance columns
    total_cells : int
        Total number of cells in the dataset
    distance_percentile : int
        Top percentile for distance (e.g., 75 = top 25%)
    frequency_threshold_pct : float
        Maximum frequency as percentage of total (e.g., 1.0 = 1%)
        
    Returns:
    --------
    list: Rare cell type names
    """
    distance_threshold = np.percentile(cell_type_df['distance'], distance_percentile)
    frequency_threshold = total_cells * (frequency_threshold_pct / 100.0)
    
    rare_mask = (cell_type_df['distance'] > distance_threshold) & \
                (cell_type_df['n_cells'] < frequency_threshold)
    rare_types = cell_type_df[rare_mask]['cell_type'].tolist()
    return rare_types, distance_threshold, frequency_threshold

def analyze_rare_cell_types(dataset, ref, label_key='celltype'):
    """Analyze rare cell types using multiple definitions."""
    print(f"\n{'='*80}")
    print(f"Analyzing rare cell types for {dataset} (reference {ref})")
    print(f"{'='*80}")
    
    # Load data
    adata = load_reference_data(dataset, ref)
    total_cells = adata.shape[0]
    print(f"Total cells: {total_cells:,}")
    
    # Calculate distances
    cell_type_df = calculate_cell_type_distances(adata, label_key=label_key)
    cell_type_df['frequency_pct'] = (cell_type_df['n_cells'] / total_cells) * 100
    cell_type_df = cell_type_df.sort_values('n_cells', ascending=True)
    
    print(f"\nCell type frequencies:")
    print(cell_type_df[['cell_type', 'n_cells', 'frequency_pct', 'distance']].to_string(index=False))
    
    results = {}
    
    # 1. Percentile-based (bottom 10th percentile)
    rare_10pct, threshold_10pct = identify_rare_cell_types_percentile(cell_type_df, percentile=10)
    results['percentile_10'] = {
        'rare_types': rare_10pct,
        'threshold': threshold_10pct,
        'threshold_pct': (threshold_10pct / total_cells) * 100,
        'n_rare': len(rare_10pct)
    }
    
    # 2. Percentile-based (bottom 5th percentile)
    rare_5pct, threshold_5pct = identify_rare_cell_types_percentile(cell_type_df, percentile=5)
    results['percentile_5'] = {
        'rare_types': rare_5pct,
        'threshold': threshold_5pct,
        'threshold_pct': (threshold_5pct / total_cells) * 100,
        'n_rare': len(rare_5pct)
    }
    
    # 3. Z-score based
    rare_zscore, z_threshold = identify_rare_cell_types_zscore(cell_type_df, z_threshold=-1.5)
    results['zscore'] = {
        'rare_types': rare_zscore,
        'z_threshold': z_threshold,
        'n_rare': len(rare_zscore)
    }
    
    # 4. Relative frequency (10% of median)
    rare_rel, rel_threshold = identify_rare_cell_types_relative_frequency(
        cell_type_df, relative_threshold=0.1)
    results['relative_freq'] = {
        'rare_types': rare_rel,
        'threshold': rel_threshold,
        'threshold_pct': (rel_threshold / total_cells) * 100,
        'n_rare': len(rare_rel)
    }
    
    # 5. Multi-criteria (low freq AND high distance)
    rare_multi, freq_thresh, dist_thresh = identify_rare_cell_types_multi_criteria(
        cell_type_df, freq_percentile=25, distance_percentile=75)
    results['multi_criteria'] = {
        'rare_types': rare_multi,
        'freq_threshold': freq_thresh,
        'freq_threshold_pct': (freq_thresh / total_cells) * 100,
        'distance_threshold': dist_thresh,
        'n_rare': len(rare_multi)
    }
    
    # 6. Distance-based (top 25% by distance)
    rare_dist, dist_threshold = identify_rare_cell_types_distance_based(
        cell_type_df, distance_percentile=75)
    results['distance_based'] = {
        'rare_types': rare_dist,
        'distance_threshold': dist_threshold,
        'n_rare': len(rare_dist)
    }
    
    # 7. Distance AND frequency (top 25% distance AND < 1% frequency)
    rare_dist_freq, dist_thresh, freq_thresh = identify_rare_cell_types_distance_and_frequency(
        cell_type_df, total_cells, distance_percentile=75, frequency_threshold_pct=1.0)
    results['distance_and_frequency'] = {
        'rare_types': rare_dist_freq,
        'distance_threshold': dist_thresh,
        'frequency_threshold': freq_thresh,
        'frequency_threshold_pct': 1.0,
        'n_rare': len(rare_dist_freq)
    }
    
    # Print summary
    print(f"\n{'='*80}")
    print("RARE CELL TYPE DEFINITIONS SUMMARY")
    print(f"{'='*80}")
    
    for method, result in results.items():
        print(f"\n{method.upper().replace('_', ' ')}:")
        print(f"  Number of rare cell types: {result['n_rare']}")
        if result['n_rare'] > 0:
            print(f"  Rare cell types: {', '.join(result['rare_types'])}")
            if 'threshold_pct' in result:
                print(f"  Threshold: {result['threshold']:.0f} cells ({result['threshold_pct']:.3f}%)")
            elif 'freq_threshold_pct' in result and 'frequency_threshold_pct' not in result:
                print(f"  Frequency threshold: {result['freq_threshold']:.0f} cells ({result['freq_threshold_pct']:.3f}%)")
                print(f"  Distance threshold: {result['distance_threshold']:.2f}")
            elif 'frequency_threshold_pct' in result:
                print(f"  Distance threshold: {result['distance_threshold']:.2f}")
                print(f"  Frequency threshold: {result['frequency_threshold']:.0f} cells (< {result['frequency_threshold_pct']:.1f}%)")
            elif 'distance_threshold' in result and 'frequency_threshold' not in result:
                print(f"  Distance threshold: {result['distance_threshold']:.2f}")
            elif 'z_threshold' in result:
                print(f"  Z-score threshold: {result['z_threshold']}")
    
    return results, cell_type_df

def main():
    """Main function to analyze all datasets."""
    dataset_configs = [
        {'name': 'mcc_01', 'ref': 30},
        {'name': 'mcc_05', 'ref': 30},
        {'name': 'mcc', 'ref': 30},
        {'name': 'lcmv', 'ref': 34},
    ]
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, 'results')
    os.makedirs(output_dir, exist_ok=True)
    
    all_summaries = []
    
    for config in dataset_configs:
        dataset = config['name']
        ref = config['ref']
        
        try:
            results, cell_type_df = analyze_rare_cell_types(dataset, ref)
            
            # Save detailed results
            summary_df = pd.DataFrame([
                {
                    'dataset': dataset,
                    'method': method,
                    'n_rare_types': result['n_rare'],
                    'rare_types': ', '.join(result['rare_types']) if result['n_rare'] > 0 else '',
                    **{k: v for k, v in result.items() if k != 'rare_types'}
                }
                for method, result in results.items()
            ])
            all_summaries.append(summary_df)
            
            # Save cell type details with rare flags
            cell_type_df['rare_percentile_10'] = cell_type_df['cell_type'].isin(results['percentile_10']['rare_types'])
            cell_type_df['rare_percentile_5'] = cell_type_df['cell_type'].isin(results['percentile_5']['rare_types'])
            cell_type_df['rare_zscore'] = cell_type_df['cell_type'].isin(results['zscore']['rare_types'])
            cell_type_df['rare_relative_freq'] = cell_type_df['cell_type'].isin(results['relative_freq']['rare_types'])
            cell_type_df['rare_multi_criteria'] = cell_type_df['cell_type'].isin(results['multi_criteria']['rare_types'])
            cell_type_df['rare_distance_based'] = cell_type_df['cell_type'].isin(results['distance_based']['rare_types'])
            cell_type_df['rare_distance_and_frequency'] = cell_type_df['cell_type'].isin(results['distance_and_frequency']['rare_types'])
            cell_type_df['dataset'] = dataset
            
            output_file = os.path.join(output_dir, f'{dataset}_rare_cell_types_analysis.csv')
            cell_type_df.to_csv(output_file, index=False)
            print(f"\nDetailed results saved to: {output_file}")
            
        except Exception as e:
            print(f"Error processing {dataset}: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    # Save combined summary
    if len(all_summaries) > 0:
        combined_summary = pd.concat(all_summaries, ignore_index=True)
        summary_file = os.path.join(output_dir, 'rare_cell_type_definitions_summary.csv')
        combined_summary.to_csv(summary_file, index=False)
        print(f"\n{'='*80}")
        print(f"Combined summary saved to: {summary_file}")
        print(f"{'='*80}")

if __name__ == "__main__":
    main()

