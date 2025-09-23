#!/usr/bin/env python
# coding: utf-8

"""
Generate combined UMAP visualizations for both MCC and LCMV datasets.
This script creates a figure with two rows of UMAP visualizations showing the distribution of sampled points and rare cell types.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import pickle
from typing import List, Dict, Tuple, Any
import logging
import argparse

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CombinedUMAPVisualizer:
    # Data-specific constants
    REFERENCES = [5, 10, 20, 25, 30]
    METHODS = ['random', 'sps', 'hopper', 'atomic']
    SIZES = [50000, 100000, 200000, 300000]
    REPS = [i for i in range(5)]
    
    # Plotting and figure constants
    FIGSIZE = (18, 12)  # Increased height for two rows
    PANEL_TITLES = ['sps', 'rare', 'random']
    COLORS_TO_PLOT = ['sps', 'rare', 'random']
    AXIS_LABEL_FONT_SIZE = 20
    PANEL_TITLE_FONT_SIZE = 20
    LEGEND_FONT_SIZE = 20
    LEGEND_HANDLE_TEXTPAD = 0.2
    LEGEND_LABEL_SPACING = 0.2
    POINT_SIZE = 15

    def __init__(self, data_path: str, test_mode: bool = False, figures_dir: str = "figures"):
        """Initialize the visualizer with data paths and configuration."""
        self.data_path = data_path
        self.test_mode = test_mode
        self.figures_dir = figures_dir
        self._setup_plotting_params()
        self.palettes = {
            'sps': {False: '#d3d3d3', True: '#e74c3c'},
            'random': {False: '#d3d3d3', True: '#e74c3c'},
            'rare': {False: '#d3d3d3', True: '#3498db'}
        }

    def _setup_plotting_params(self):
        plt.rcParams.update({
            "figure.figsize": (6, 6),
            "figure.dpi": 350,
            "savefig.dpi": 350,
            "font.size": 18,
            "axes.titlesize": 18,
            "axes.labelsize": 18,
            "legend.fontsize": 18,
            "legend.markerscale": 1.5,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
        })

    def load_reference_data(self, dataset: str, ref: int) -> sc.AnnData:
        """Load and preprocess reference data for either MCC or LCMV dataset."""
        directory = f"{dataset}/benchmark"
        path = os.path.join(self.data_path, directory)
        address = os.path.join(path, f"{ref}/adata.h5ad")
        adata_ref = sc.read_h5ad(address)
        
        if self.test_mode:
            logger.info("Creating test subset of data...")
            label_key = 'celltype'
            cell_types = adata_ref.obs[label_key].unique()
            test_indices = []
            test_size = 1000
            for ct in cell_types:
                ct_indices = np.where(adata_ref.obs[label_key] == ct)[0]
                n_samples = min(len(ct_indices), test_size // len(cell_types))
                test_indices.extend(np.random.choice(ct_indices, n_samples, replace=False))
            self.subset_obs_names = adata_ref.obs_names[test_indices].tolist()
            adata_ref = adata_ref[test_indices].copy()
            logger.info(f"Test subset created with {adata_ref.shape[0]} cells")
        else:
            self.subset_obs_names = None
        adata_ref.obs['celltype'] = adata_ref.obs['celltype'].astype('category')
        adata_ref.var.index = adata_ref.var.index.astype('object')
        return adata_ref

    def identify_rare_cell_types(self, adata: sc.AnnData) -> sc.AnnData:
        """Identify and mark rare cell types in the dataset."""
        real_counts = adata.obs['celltype'].value_counts()
        rare_dict = real_counts < adata.shape[0]/100
        rare_types = [x for x in rare_dict.keys() if rare_dict[x]]
        adata.obs['rare'] = adata.obs.apply(
            lambda x: True if x['celltype'] in rare_types else False, 
            axis=1
        )
        return adata

    def load_sampling_results(self, dataset: str, ref: int, method: str, size: int, rep: int, n_cells: int = None) -> np.ndarray:
        """Load sampling results from pickle file or CSV for atomic method."""
        directory = f"{dataset}/benchmark"
        path = os.path.join(self.data_path, directory)
        
        if self.test_mode:
            logger.info(f"Creating test sampling indices for {method}")
            # Use n_cells if provided, otherwise default to 1000
            if n_cells is None:
                n_cells = 1000
            test_sample_size = min(100, n_cells)
            sampled_indices = np.random.choice(n_cells, size=test_sample_size, replace=False)
            # Return obs_names, not integer indices
            return np.array(self.subset_obs_names)[sampled_indices]

        if method == 'atomic':
            samples_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/results.csv")
            indices = pd.read_csv(samples_address)['x'].values.astype(str)
        else:
            samples_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/results.pkl")
            with open(samples_address, 'rb') as handle:
                samples = pickle.load(handle)
            indices = samples[0]
        return indices

    def mark_sampled_points(self, adata: sc.AnnData, sample_ids, method: str) -> sc.AnnData:
        """Mark sampled points in the dataset."""
        adata.obs[method] = 'False'
        if isinstance(sample_ids[0], (str, np.str_)):
            row_indices = adata.obs.index.get_indexer(sample_ids)
            row_indices = row_indices[row_indices >= 0]
        else:
            row_indices = sample_ids
        adata.obs.iloc[row_indices, adata.obs.columns.get_loc(method)] = 'True'
        logger.info(f"Sampled points for {method}: {adata.obs[method].value_counts()}")
        return adata

    def compute_umap(self, adata: sc.AnnData) -> sc.AnnData:
        """Compute PCA and UMAP embeddings."""
        logger.info("Computing PCA...")
        sc.tl.pca(adata, svd_solver='arpack')
        logger.info("Computing neighborhood graph...")
        sc.pp.neighbors(adata, n_neighbors=15)
        logger.info("Computing UMAP...")
        sc.tl.umap(adata)
        return adata

    def set_color_schemes(self, adata: sc.AnnData) -> sc.AnnData:
        """Set color schemes in the AnnData object."""
        for method in self.COLORS_TO_PLOT:
            adata.obs[method] = adata.obs[method].astype('category')
            categories = adata.obs[method].cat.categories
            def to_bool(val):
                if val in [True, False]:
                    return val
                if str(val).lower() == 'true':
                    return True
                if str(val).lower() == 'false':
                    return False
                return val
            palette = {cat: self.palettes[method][to_bool(cat)] for cat in categories}
            adata.uns[f'{method}_colors'] = [palette[cat] for cat in categories]
        return adata

    def generate_combined_umaps(self, mcc_adata: sc.AnnData, lcmv_adata: sc.AnnData):
        """Generate combined UMAP plots for both datasets in a 2x3 grid."""
        filename = 'test_combined_umaps.jpg' if self.test_mode else 'combined_umaps.jpg'
        output_file = os.path.join(self.figures_dir, filename)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # Create figure with 2 rows and 3 columns
        fig, axs = plt.subplots(2, 3, figsize=self.FIGSIZE)
        
        # Plot MCC UMAPs (top row)
        for i, color_key in enumerate(self.COLORS_TO_PLOT):
            ax = axs[0, i]
            sc.pl.umap(
                mcc_adata,
                color=color_key,
                ax=ax,
                show=False,
                marker='s',
                size=self.POINT_SIZE,
                legend_loc='none',
            )
            ax.set_title(f"MCC - {self.PANEL_TITLES[i]}", fontsize=self.PANEL_TITLE_FONT_SIZE)
            ax.set_xlabel('UMAP1', fontsize=self.AXIS_LABEL_FONT_SIZE)
            if i == 0:
                ax.set_ylabel('UMAP2', fontsize=self.AXIS_LABEL_FONT_SIZE)
            else:
                ax.set_ylabel('')
                ax.tick_params(axis='y', which='both', left=False, labelleft=False)
            ax.tick_params(axis='both', which='major', labelsize=10)
            ax.set_aspect('equal')

        # Plot LCMV UMAPs (bottom row)
        for i, color_key in enumerate(self.COLORS_TO_PLOT):
            ax = axs[1, i]
            sc.pl.umap(
                lcmv_adata,
                color=color_key,
                ax=ax,
                show=False,
                marker='s',
                size=self.POINT_SIZE,
                legend_loc='none',
            )
            ax.set_title(f"LCMV - {self.PANEL_TITLES[i]}", fontsize=self.PANEL_TITLE_FONT_SIZE)
            ax.set_xlabel('UMAP1', fontsize=self.AXIS_LABEL_FONT_SIZE)
            if i == 0:
                ax.set_ylabel('UMAP2', fontsize=self.AXIS_LABEL_FONT_SIZE)
            else:
                ax.set_ylabel('')
                ax.tick_params(axis='y', which='both', left=False, labelleft=False)
            ax.tick_params(axis='both', which='major', labelsize=10)
            ax.set_aspect('equal')
            
            # Add legend only for bottom row
            categories = lcmv_adata.obs[color_key].cat.categories
            palette = lcmv_adata.uns[f'{color_key}_colors']
            handles = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=palette[j], label=str(cat)) 
                      for j, cat in enumerate(categories)]
            ax.legend(
                handles=handles,
                labels=[str(cat) for cat in categories],
                fontsize=self.LEGEND_FONT_SIZE,
                loc='upper center',
                bbox_to_anchor=(0.5, -0.28),
                ncol=len(categories),
                frameon=False,
                handletextpad=self.LEGEND_HANDLE_TEXTPAD,
                labelspacing=self.LEGEND_LABEL_SPACING,
            )

        plt.tight_layout(rect=[0, 0.28, 1, 1])  # Increased bottom margin for more space between rows
        plt.subplots_adjust(hspace=0.45)  # Add more vertical space between rows
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

def main():
    parser = argparse.ArgumentParser(description='Generate combined UMAP visualizations for MCC and LCMV data')
    parser.add_argument('--test', action='store_true', help='Run in test mode with smaller dataset')
    args = parser.parse_args()
    
    mcc_ref = 5  # Reference size for MCC
    lcmv_ref = 1  # Reference size for LCMV
    method = 'sps'
    method2 = 'random'
    size = 100000
    rep = 0
    data_path = 'projects/sparsesampler_reproducibility/data'
    
    try:
        visualizer = CombinedUMAPVisualizer(data_path, test_mode=args.test)
        
        # Process MCC data
        mcc_adata = visualizer.load_reference_data('mcc', mcc_ref)
        mcc_adata = visualizer.identify_rare_cell_types(mcc_adata)
        mcc_n_cells = mcc_adata.shape[0] if args.test else None
        mcc_indices = visualizer.load_sampling_results('mcc', mcc_ref, method, size, rep, n_cells=mcc_n_cells)
        mcc_indices2 = visualizer.load_sampling_results('mcc', mcc_ref, method2, size, rep, n_cells=mcc_n_cells)
        mcc_adata = visualizer.mark_sampled_points(mcc_adata, mcc_indices, 'sps')
        mcc_adata = visualizer.mark_sampled_points(mcc_adata, mcc_indices2, 'random')
        mcc_adata = visualizer.compute_umap(mcc_adata)
        mcc_adata = visualizer.set_color_schemes(mcc_adata)
        
        # Process LCMV data
        lcmv_adata = visualizer.load_reference_data('lcmv', lcmv_ref)
        lcmv_adata = visualizer.identify_rare_cell_types(lcmv_adata)
        lcmv_n_cells = lcmv_adata.shape[0] if args.test else None
        lcmv_indices = visualizer.load_sampling_results('lcmv', lcmv_ref, method, size, rep, n_cells=lcmv_n_cells)
        lcmv_indices2 = visualizer.load_sampling_results('lcmv', lcmv_ref, method2, size, rep, n_cells=lcmv_n_cells)
        lcmv_adata = visualizer.mark_sampled_points(lcmv_adata, lcmv_indices, 'sps')
        lcmv_adata = visualizer.mark_sampled_points(lcmv_adata, lcmv_indices2, 'random')
        lcmv_adata = visualizer.compute_umap(lcmv_adata)
        lcmv_adata = visualizer.set_color_schemes(lcmv_adata)
        
        # Generate combined visualization
        visualizer.generate_combined_umaps(mcc_adata, lcmv_adata)
        
        logger.info("Combined visualization generation completed successfully")
    except Exception as e:
        logger.error(f"Error in main execution: {e}")
        raise

if __name__ == "__main__":
    main() 