#!/usr/bin/env python
"""Plot EVR Index Sensitivity Analysis for SParseSampler revision."""

import os
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import scanpy as sc

script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
analysis_dir = os.path.join(project_root, 'analysis')
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

from refined_rare_cell_type_definition import (
    calculate_cell_type_distances,
    identify_rare_cell_types_distance_and_frequency
)

DATASET_CONFIGS = {
    'lcmv': {'references': [1, 5, 10, 20, 34], 'sizes': [50000, 100000, 200000], 'freq_pct': 1.0, 'label': 'LCMV'},
    'mcc': {'references': [5, 10, 20, 25, 30], 'sizes': [50000, 100000, 200000, 300000], 'freq_pct': 1.0, 'label': 'MCC (1%)'},
    'mcc_01': {'references': [5, 10, 20, 25, 30], 'sizes': [50000, 100000, 200000, 300000], 'freq_pct': 0.1, 'label': 'MCC (0.1%)'},
    'mcc_05': {'references': [5, 10, 20, 25, 30], 'sizes': [50000, 100000, 200000, 300000], 'freq_pct': 0.5, 'label': 'MCC (0.5%)'},
}

EVR_INDICES = list(range(15, 21))
file_path_env = os.path.join(project_root, 'data')
RESULTS_DIR = os.path.join(file_path_env, 'test_feature_index')
FIGURE_DIR = os.path.join(project_root, 'figures', 'revision')
os.makedirs(FIGURE_DIR, exist_ok=True)

def load_obs(dataset, ref):
    adata_path = os.path.join(file_path_env, f"{dataset}/benchmark/{ref}/adata.h5ad")
    if os.path.exists(adata_path):
        return sc.read_h5ad(adata_path).obs
    return None

def get_rare_types(dataset, ref, obs):
    adata_path = os.path.join(file_path_env, f"{dataset}/benchmark/{ref}/adata.h5ad")
    if not os.path.exists(adata_path):
        return []
    adata = sc.read_h5ad(adata_path)
    df = calculate_cell_type_distances(adata, label_key='celltype')
    rare, _, _ = identify_rare_cell_types_distance_and_frequency(df, obs.shape[0], 75, DATASET_CONFIGS[dataset]['freq_pct'])
    return rare

def calc_coverage(dataset):
    cfg = DATASET_CONFIGS[dataset]
    all_data = []
    for ref in cfg['references']:
        obs = load_obs(dataset, ref)
        if obs is None:
            continue
        rare_types = get_rare_types(dataset, ref, obs)
        for evr in EVR_INDICES:
            for size in cfg['sizes']:
                path = os.path.join(RESULTS_DIR, dataset, str(ref), str(evr), str(size), '0', 'results.pkl')
                if not os.path.exists(path):
                    continue
                try:
                    with open(path, 'rb') as f:
                        data = pickle.load(f)
                    idx = data[0] if isinstance(data, tuple) else data
                    if isinstance(idx, tuple):
                        idx = idx[0]
                    sampled = obs.iloc[idx]
                    rare_count = sampled[sampled['celltype'].isin(rare_types)].shape[0]
                    all_data.append({'ref': ref, 'evr': evr, 'size': size, 'rare_count': rare_count})
                except:
                    continue
    return pd.DataFrame(all_data)

def main():
    plt.rcParams.update({'font.size': 11, 'figure.figsize': [14, 10], 'savefig.dpi': 300})
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for idx, dataset in enumerate(['mcc', 'mcc_05', 'mcc_01', 'lcmv']):
        print(f"Processing {dataset}...")
        df = calc_coverage(dataset)
        cfg = DATASET_CONFIGS[dataset]
        ax = axes[idx]
        
        if len(df) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(cfg['label'])
            continue
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(cfg['references'])))
        target_size = 100000
        
        for j, ref in enumerate(cfg['references']):
            ref_df = df[(df['ref'] == ref) & (df['size'] == target_size)]
            if len(ref_df) == 0:
                continue
            avg = ref_df.groupby('evr')['rare_count'].mean()
            ax.plot(avg.index, avg.values, '-o', color=colors[j], label=f'{ref}M', lw=1.5, ms=6)
        
        ax.axvline(x=18, color='red', linestyle='--', alpha=0.5, label='Default (18)')
        ax.set_xlabel('EVR Index')
        ax.set_ylabel('Rare Cells Captured')
        ax.set_title(f'{cfg["label"]} (n=100,000)')
        ax.set_xticks(EVR_INDICES)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='best', fontsize=9, title='Ref Size')
    
    plt.suptitle('EVR Index Sensitivity Analysis (15-20)', fontsize=14, y=1.01)
    plt.tight_layout()
    
    out_path = os.path.join(FIGURE_DIR, 'supp_evr_sensitivity_combined.jpg')
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    print(f"Saved: {out_path}")
    
    # Also copy to paper directory
    paper_path = '/fast/AG_Ohler/ekarimi/papers/sps_paper/Fig/supp/supp_evr_sensitivity_combined.jpg'
    import shutil
    shutil.copy(out_path, paper_path)
    print(f"Copied to: {paper_path}")

if __name__ == "__main__":
    main()
