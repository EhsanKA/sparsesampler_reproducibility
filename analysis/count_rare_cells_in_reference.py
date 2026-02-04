#!/usr/bin/env python
# coding: utf-8

"""
Count rare cells in the reference sets (sampled cells) for SPS vs Random.
"""

import pandas as pd
import numpy as np
import scanpy as sc
import pickle
import os
import sys

# Add analysis directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

from refined_rare_cell_type_definition import (
    load_reference_data,
    identify_rare_cell_types_distance_and_frequency,
    calculate_cell_type_distances
)

# Define rare cell types for lcmv
rare_cell_types = [
    'NK1_1_TCRgd_T',
    'interacting',
    'cDC2',
    'pDCs',
    'CD4_LCMV_spec'
]

def load_sampling_indices(dataset, ref, method, size, rep, data_path):
    """Load sampling indices."""
    directory = f"{dataset}/benchmark"
    path = os.path.join(data_path, directory)
    
    # Try .npy file first
    npy_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/indices.npy")
    if os.path.isfile(npy_address):
        indices = np.load(npy_address)
        return indices
    
    # Fall back to pickle
    samples_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/results.pkl")
    if os.path.isfile(samples_address):
        with open(samples_address, 'rb') as handle:
            samples = pickle.load(handle)
        return samples[0]
    
    return None

def count_rare_cells_in_reference(adata_full, reference_indices, rare_cell_types, label_key='celltype'):
    """Count rare cells in the reference set."""
    if isinstance(reference_indices[0], (str, np.str_)):
        ref_mask = adata_full.obs.index.isin(reference_indices)
        ref_idx = np.where(ref_mask)[0]
    else:
        ref_idx = reference_indices
    
    reference_adata = adata_full[ref_idx].copy()
    
    # Count rare cells by type
    cell_types = reference_adata.obs[label_key].values
    rare_counts = {}
    total_rare = 0
    
    for rare_type in rare_cell_types:
        count = np.sum(cell_types == rare_type)
        rare_counts[rare_type] = count
        total_rare += count
    
    return rare_counts, total_rare, len(ref_idx)

# Configuration
dataset = 'lcmv'
ref = 34
size = 100000
rep = 0
data_path = '/fast/AG_Ohler/ekarimi/projects/sparseFlow_benchmarking/data'

print("Loading full dataset...")
adata_full = load_reference_data(dataset, ref, base_path=data_path)
print(f"Full dataset: {adata_full.shape[0]:,} cells")

# Identify rare cell types (to verify)
cell_type_df = calculate_cell_type_distances(adata_full, label_key='celltype')
total_cells = adata_full.shape[0]
rare_types, _, _ = identify_rare_cell_types_distance_and_frequency(
    cell_type_df, total_cells, distance_percentile=75, frequency_threshold_pct=1.0
)
print(f"\nIdentified {len(rare_types)} rare cell types: {rare_types}")

# Count rare cells in full dataset
full_counts = {}
for rare_type in rare_cell_types:
    count = np.sum(adata_full.obs['celltype'] == rare_type)
    full_counts[rare_type] = count

print(f"\nRare cell counts in FULL dataset:")
for rare_type in rare_cell_types:
    print(f"  {rare_type}: {full_counts[rare_type]:,}")

# Load sampling indices and count rare cells in reference sets
print("\n" + "="*80)
print("Counting rare cells in REFERENCE sets (sampled cells):")
print("="*80)

methods = ['random', 'sps']
results = []

for method in methods:
    print(f"\n{method.upper()} method:")
    reference_indices = load_sampling_indices(dataset, ref, method, size, rep, data_path)
    
    if reference_indices is None:
        print(f"  Error: Could not load indices for {method}")
        continue
    
    rare_counts, total_rare, n_ref = count_rare_cells_in_reference(
        adata_full, reference_indices, rare_cell_types
    )
    
    print(f"  Total reference cells: {n_ref:,}")
    print(f"  Total rare cells in reference: {total_rare:,}")
    
    result_row = {'method': method, 'total_reference': n_ref, 'total_rare': total_rare}
    for rare_type in rare_cell_types:
        count = rare_counts.get(rare_type, 0)
        pct = (count / full_counts[rare_type] * 100) if full_counts[rare_type] > 0 else 0
        print(f"    {rare_type}: {count:,} ({pct:.2f}% of full dataset)")
        result_row[f'{rare_type}_count'] = count
        result_row[f'{rare_type}_pct'] = pct
    
    results.append(result_row)

# Create comparison DataFrame
comparison_df = pd.DataFrame(results)

# Save to CSV
output_file = '/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/results/semitones_downstream/lcmv_rare_cells_in_reference_comparison.csv'
comparison_df.to_csv(output_file, index=False)
print(f"\n{'='*80}")
print(f"Results saved to: {output_file}")
print(f"{'='*80}")

# Summary comparison
print("\nSUMMARY COMPARISON:")
print("-" * 80)
for rare_type in rare_cell_types:
    random_count = comparison_df.loc[comparison_df['method'] == 'random', f'{rare_type}_count'].values[0]
    sps_count = comparison_df.loc[comparison_df['method'] == 'sps', f'{rare_type}_count'].values[0]
    diff = sps_count - random_count
    pct_diff = (diff / random_count * 100) if random_count > 0 else 0
    print(f"{rare_type:20s} | Random: {random_count:6,} | SPS: {sps_count:6,} | Diff: {diff:+6,} ({pct_diff:+.1f}%)")

total_random = comparison_df.loc[comparison_df['method'] == 'random', 'total_rare'].values[0]
total_sps = comparison_df.loc[comparison_df['method'] == 'sps', 'total_rare'].values[0]
total_diff = total_sps - total_random
total_pct_diff = (total_diff / total_random * 100) if total_random > 0 else 0
print(f"{'TOTAL':20s} | Random: {total_random:6,} | SPS: {total_sps:6,} | Diff: {total_diff:+6,} ({total_pct_diff:+.1f}%)")

