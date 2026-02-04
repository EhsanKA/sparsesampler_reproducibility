#!/usr/bin/env python
# coding: utf-8

"""
Compare SEMITONES enrichment scores for rare cell types between random and SPS methods.
"""

import pandas as pd
import os

# Define rare cell types for lcmv
rare_cell_types = [
    'NK1_1_TCRgd_T',
    'interacting',
    'cDC2',
    'pDCs',
    'CD4_LCMV_spec'
]

# File paths
random_file = '/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/results/semitones_downstream/lcmv_random_size100000_rep0_summary.csv'
sps_file = '/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/results/semitones_downstream/lcmv_sps_size100000_rep0_summary.csv'
output_file = '/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/results/semitones_downstream/lcmv_rare_cells_comparison.csv'

# Load the data
print("Loading data...")
random_df = pd.read_csv(random_file)
sps_df = pd.read_csv(sps_file)

# Filter for rare cell types
random_rare = random_df[random_df['cell_type'].isin(rare_cell_types)].copy()
sps_rare = sps_df[sps_df['cell_type'].isin(rare_cell_types)].copy()

# Merge on cell_type
comparison = pd.merge(
    random_rare[['cell_type', 'mean_enrichment', 'median_enrichment', 'n_cells', 'pct_high_enrichment']],
    sps_rare[['cell_type', 'mean_enrichment', 'median_enrichment', 'n_cells', 'pct_high_enrichment']],
    on='cell_type',
    suffixes=('_random', '_sps')
)

# Calculate differences (SPS - Random)
comparison['mean_enrichment_diff'] = comparison['mean_enrichment_sps'] - comparison['mean_enrichment_random']
comparison['median_enrichment_diff'] = comparison['median_enrichment_sps'] - comparison['median_enrichment_random']
comparison['pct_high_enrichment_diff'] = comparison['pct_high_enrichment_sps'] - comparison['pct_high_enrichment_random']

# Calculate percent change
comparison['mean_enrichment_pct_change'] = (comparison['mean_enrichment_diff'] / comparison['mean_enrichment_random']) * 100
comparison['median_enrichment_pct_change'] = (comparison['median_enrichment_diff'] / comparison['median_enrichment_random']) * 100
comparison['pct_high_enrichment_pct_change'] = (comparison['pct_high_enrichment_diff'] / comparison['pct_high_enrichment_random']) * 100

# Reorder columns for better readability
column_order = [
    'cell_type',
    'n_cells_random',
    'n_cells_sps',
    'mean_enrichment_random',
    'mean_enrichment_sps',
    'mean_enrichment_diff',
    'mean_enrichment_pct_change',
    'median_enrichment_random',
    'median_enrichment_sps',
    'median_enrichment_diff',
    'median_enrichment_pct_change',
    'pct_high_enrichment_random',
    'pct_high_enrichment_sps',
    'pct_high_enrichment_diff',
    'pct_high_enrichment_pct_change'
]

comparison = comparison[column_order]

# Sort by mean_enrichment_diff (SPS - Random) to see which ones improved most
comparison = comparison.sort_values('mean_enrichment_diff', ascending=False)

# Save to CSV
comparison.to_csv(output_file, index=False)

print(f"\nComparison saved to: {output_file}")
print(f"\nSummary:")
print(f"Number of rare cell types compared: {len(comparison)}")
print(f"\nMean enrichment differences (SPS - Random):")
print(comparison[['cell_type', 'mean_enrichment_diff', 'mean_enrichment_pct_change']].to_string(index=False))
print(f"\nFull comparison:")
print(comparison.to_string(index=False))

