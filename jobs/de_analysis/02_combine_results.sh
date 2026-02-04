#!/bin/bash

# Script to combine DE analysis results after all jobs complete
# Run this after all SLURM jobs finish

# Load cluster configuration
if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found. Please copy from template and configure."
    exit 1
fi

# Set PROJECT_ROOT environment variable
export PROJECT_ROOT=${PROJECT_ROOT}

# Change directory to the analysis folder
cd ${PROJECT_ROOT}/analysis

echo "=========================================="
echo "Combining DE Analysis Results"
echo "=========================================="

# Run Python script to combine results
python <<EOF
import os
import sys
import pandas as pd
import glob

results_dir = os.path.join(os.environ['PROJECT_ROOT'], 'analysis', 'results', 'de_analysis')
os.makedirs(results_dir, exist_ok=True)

summary_file = os.path.join(results_dir, 'de_analysis_summary.csv')

# Check if summary file exists and read it
if os.path.exists(summary_file):
    print(f"Reading existing summary file: {summary_file}")
    summary_df = pd.read_csv(summary_file)
    print(f"Found {len(summary_df)} existing entries")
    
    # Remove duplicates (keep last)
    summary_df = summary_df.drop_duplicates(
        subset=['dataset', 'method', 'sample_size', 'rep'],
        keep='last'
    )
    print(f"After removing duplicates: {len(summary_df)} entries")
else:
    print("No existing summary file found, creating new one")
    summary_df = pd.DataFrame()

# Find all detailed result files
detailed_files = glob.glob(os.path.join(results_dir, '*_de_detailed.csv'))
print(f"\nFound {len(detailed_files)} detailed result files")

# Process each detailed file
new_entries = []

for detailed_file in detailed_files:
    # Extract metadata from filename
    basename = os.path.basename(detailed_file)
    # Format: {dataset}_{method}_size{size}_rep{rep}_de_detailed.csv
    parts = basename.replace('_de_detailed.csv', '').split('_')
    
    if len(parts) < 4:
        print(f"Warning: Cannot parse filename {basename}, skipping")
        continue
    
    dataset = parts[0]
    method = parts[1]
    size = int(parts[2].replace('size', ''))
    rep = int(parts[3].replace('rep', ''))
    
    # Check if this entry already exists in summary
    if len(summary_df) > 0:
        existing = summary_df[
            (summary_df['dataset'] == dataset) & 
            (summary_df['method'] == method) & 
            (summary_df['sample_size'] == size) & 
            (summary_df['rep'] == rep)
        ]
        if len(existing) > 0:
            print(f"Skipping {basename} (already in summary)")
            continue
    
    # Read detailed file
    try:
        detailed_df = pd.read_csv(detailed_file)
    except Exception as e:
        print(f"Error reading {detailed_file}: {e}")
        continue
    
    if len(detailed_df) == 0:
        print(f"Warning: {detailed_file} is empty, skipping")
        continue
    
    # Calculate summary statistics
    overlaps = detailed_df['overlap_fraction'].values
    mean_overlap = overlaps.mean()
    median_overlap = pd.Series(overlaps).median()
    std_overlap = overlaps.std()
    min_overlap = overlaps.min()
    max_overlap = overlaps.max()
    
    result_entry = {
        'dataset': dataset,
        'method': method,
        'sample_size': size,
        'rep': rep,
        'mean_overlap': mean_overlap,
        'median_overlap': median_overlap,
        'std_overlap': std_overlap,
        'min_overlap': min_overlap,
        'max_overlap': max_overlap,
        'n_cell_types_analyzed': len(detailed_df),
        'full_n_cells': None,  # Will be filled from individual job outputs
        'sampled_n_cells': None,
    }
    
    new_entries.append(result_entry)
    print(f"Processed {basename}: {len(detailed_df)} cell types, mean overlap={mean_overlap:.2%}")

# Combine with existing summary
if new_entries:
    new_df = pd.DataFrame(new_entries)
    if len(summary_df) > 0:
        summary_df = pd.concat([summary_df, new_df], ignore_index=True)
    else:
        summary_df = new_df
    
    # Remove duplicates again (keep last)
    summary_df = summary_df.drop_duplicates(
        subset=['dataset', 'method', 'sample_size', 'rep'],
        keep='last'
    )
    
    # Sort by dataset, method, sample_size, rep
    summary_df = summary_df.sort_values(['dataset', 'method', 'sample_size', 'rep']).reset_index(drop=True)
    
    summary_df.to_csv(summary_file, index=False)
    print(f"\nSaved combined summary to: {summary_file}")
    print(f"Total entries: {len(summary_df)}")
    print(f"\nSummary:")
    print(summary_df[['dataset', 'method', 'sample_size', 'mean_overlap', 'n_cell_types_analyzed']].to_string(index=False))
else:
    if len(summary_df) > 0:
        print("\nNo new entries to add, using existing summary")
        print(f"\nSummary:")
        print(summary_df[['dataset', 'method', 'sample_size', 'mean_overlap', 'n_cell_types_analyzed']].to_string(index=False))
    else:
        print("No results to combine!")
        sys.exit(1)
EOF

echo ""
echo "Results combined successfully!"
echo "Summary file: ${PROJECT_ROOT}/analysis/results/de_analysis/de_analysis_summary.csv"

