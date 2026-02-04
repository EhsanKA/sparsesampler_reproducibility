#!/usr/bin/env python
"""Plot EVR Index Sensitivity from pre-computed tables."""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.join(script_dir, 'tables')
project_root = os.path.dirname(os.path.dirname(script_dir))
FIGURE_DIR = os.path.join(project_root, 'figures', 'revision')
PAPER_DIR = '/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/figures/revision'
os.makedirs(FIGURE_DIR, exist_ok=True)

DATASETS = {
    'mcc': {'label': 'MCC (1%)', 'refs': ['5M', '10M', '20M', '25M', '30M']},
    'mcc_05': {'label': 'MCC (0.5%)', 'refs': ['5M', '10M', '20M', '25M', '30M']},
    'mcc_01': {'label': 'MCC (0.1%)', 'refs': ['5M', '10M', '20M', '25M', '30M']},
    'lcmv': {'label': 'LCMV', 'refs': ['1M', '5M', '10M', '20M', '34M']},
}

EVR_RANGE = list(range(15, 21))  # Focus on EVR 15-20

def main():
    plt.rcParams.update({
        'font.size': 12, 'axes.labelsize': 12, 'axes.titlesize': 14,
        'legend.fontsize': 10, 'xtick.labelsize': 11, 'ytick.labelsize': 11,
        'figure.figsize': [14, 10], 'lines.markersize': 8, 'lines.linewidth': 2,
        'savefig.dpi': 300
    })
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for idx, (dataset, cfg) in enumerate(DATASETS.items()):
        ax = axes[idx]
        
        # Load the 100k sample size table
        table_path = os.path.join(TABLES_DIR, f'{dataset}_feature_index_table_size_100000.csv')
        if not os.path.exists(table_path):
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(cfg['label'])
            continue
        
        df = pd.read_csv(table_path)
        df = df[df['feature_index'].isin(EVR_RANGE)]
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(cfg['refs'])))
        
        for j, ref in enumerate(cfg['refs']):
            if ref in df.columns:
                ax.plot(df['feature_index'], df[ref], '-o', color=colors[j], 
                       label=ref, linewidth=2, markersize=8)
        
        ax.axvline(x=18, color='red', linestyle='--', alpha=0.6, linewidth=2)
        ax.set_xlabel('EVR Index')
        ax.set_ylabel('Rare Cells Captured')
        ax.set_title(f'{cfg["label"]} (n=100,000)')
        ax.set_xticks(EVR_RANGE)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='best', fontsize=9, title='Ref Size')
    
    plt.suptitle('EVR Index Sensitivity Analysis (15-20)\nRed dashed line = default EVR index (18)', 
                 fontsize=14, y=1.02)
    plt.tight_layout()
    
    # Save to revision figures
    out_path = os.path.join(FIGURE_DIR, 'supp_evr_sensitivity_combined.jpg')
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    print(f"Saved: {out_path}")
    
    # Copy to paper directory
    import shutil
    paper_path = os.path.join(PAPER_DIR, 'evr_sensitivity_combined.jpg')
    shutil.copy(out_path, paper_path)
    print(f"Copied to: {paper_path}")

if __name__ == "__main__":
    main()
