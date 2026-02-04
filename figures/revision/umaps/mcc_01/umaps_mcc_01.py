#!/usr/bin/env python
# coding: utf-8

"""
Generate UMAP visualizations comparing different sampling methods in MCC_01 reference dataset.
This script creates visualizations showing the distribution of sampled points and rare cell types.
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import pickle
from typing import List, Dict, Tuple, Any
import logging
import argparse

# Add analysis directory to path to import rare cell type functions
script_dir = os.path.dirname(os.path.abspath(__file__))
analysis_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))), 'analysis')
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

from refined_rare_cell_type_definition import (
    calculate_cell_type_distances,
    identify_rare_cell_types_distance_and_frequency
)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class UMAPVisualizer:
    # Data-specific constants
    REFERENCES = [5, 10, 20, 25, 30]
    METHODS = ['random', 'sps', 'hopper', 'atomic']
    SIZES = [50000, 100000, 200000, 300000]
    REPS = [i for i in range(5)]
    LABEL_KEY = 'celltype'

    # Plotting and figure constants
    FIGSIZE_SINGLE = (3, 6)
    FIGSIZE_PANELS = (18, 6)
    WSPACE_PANELS = 0.5
    PANEL_TITLES = ['sps', 'rare', 'random']
    COLORS_TO_PLOT = ['sps', 'rare', 'random']
    AXIS_LABEL_FONT_SIZE = 20
    PANEL_TITLE_FONT_SIZE = 20
    LEGEND_FONT_SIZE = 20
    LEGEND_HANDLE_TEXTPAD = 0.2
    LEGEND_LABEL_SPACING = 0.2
    POINT_SIZE = 8

    # UMAP and clustering hyperparameters
    N_NEIGHBORS = 40
    LEIDEN_RESOLUTION = 0.6

    def __init__(self, data_path: str, directory: str = "mcc_01/benchmark", test_mode: bool = False, figures_dir: str = None):
        """Initialize the UMAP visualizer with data paths and configuration."""
        self.file_path_env = data_path
        self.directory = directory
        self.path = os.path.join(self.file_path_env, self.directory)
        self.test_mode = test_mode
        if figures_dir is None:
            # Get figures directory from config
            project_root = os.environ.get('PROJECT_ROOT')
            if project_root is None:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(script_dir))))
            self.figures_dir = os.path.join(project_root, 'figures', 'revision')
        else:
            self.figures_dir = figures_dir
        self.OBS_FEATURES = [self.LABEL_KEY]
        if self.test_mode:
            self.test_size = 1000
            logger.info(f"Running in test mode with {self.test_size} cells")
        self._setup_plotting_params()
        self.palettes = {
            'sps': {False: '#d3d3d3', True: '#e74c3c'},
            'random': {False: '#d3d3d3', True: '#e74c3c'},
            'rare': {False: '#d3d3d3', True: '#3498db'}
        }

    def _setup_plotting_params(self):
        plt.rcParams.update({
            "figure.figsize": (6, 6),     # default for single-panel
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

    def load_reference_data(self, ref: int) -> sc.AnnData:
        """Load and preprocess reference data."""
        address = os.path.join(self.path, f"{ref}/adata.h5ad")
        adata_ref = sc.read_h5ad(address)
        if self.test_mode:
            logger.info("Creating test subset of data...")
            cell_types = adata_ref.obs[self.LABEL_KEY].unique()
            test_indices = []
            for ct in cell_types:
                ct_indices = np.where(adata_ref.obs[self.LABEL_KEY] == ct)[0]
                n_samples = min(len(ct_indices), self.test_size // len(cell_types))
                test_indices.extend(np.random.choice(ct_indices, n_samples, replace=False))
            self.subset_obs_names = adata_ref.obs_names[test_indices].tolist()
            adata_ref = adata_ref[test_indices].copy()
            self.n_cells = adata_ref.shape[0]
            logger.info(f"Test subset created with {self.n_cells} cells")
        else:
            self.n_cells = None
            self.subset_obs_names = None
        adata_ref.obs[self.LABEL_KEY] = adata_ref.obs[self.LABEL_KEY].astype('category')
        adata_ref.var.index = adata_ref.var.index.astype('object')
        return adata_ref

    def identify_rare_cell_types(self, adata: sc.AnnData) -> sc.AnnData:
        """Identify and mark rare cell types using distance and frequency criteria."""
        # Calculate cell type distances
        cell_type_df = calculate_cell_type_distances(adata, label_key=self.LABEL_KEY)
        
        # Identify rare cell types using distance and frequency (1/1000 = 0.1% for mcc_01)
        total_cells = adata.shape[0]
        rare_types, _, _ = identify_rare_cell_types_distance_and_frequency(
            cell_type_df, 
            total_cells, 
            distance_percentile=75, 
            frequency_threshold_pct=0.1  # 1/1000 = 0.1%
        )
        
        adata.obs['rare'] = adata.obs.apply(
            lambda x: True if x[self.LABEL_KEY] in rare_types else False, 
            axis=1
        )
        return adata

    def load_sampling_results(self, ref: int, method: str, size: int, rep: int) -> np.ndarray:
        """Load sampling results from pickle file or CSV for atomic method."""
        if self.test_mode:
            logger.info(f"Creating test sampling indices for {method}")
            n_cells = self.n_cells if self.n_cells is not None else self.test_size
            test_sample_size = min(100, n_cells)
            sampled_indices = np.random.choice(n_cells, size=test_sample_size, replace=False)
            return np.array(self.subset_obs_names)[sampled_indices]

        if method == 'atomic':
            samples_address = os.path.join(self.path, f"{ref}/{method}/{size}/{rep}/results.csv")
            indices = pd.read_csv(samples_address)['x'].values.astype(str)
        else:
            samples_address = os.path.join(self.path, f"{ref}/{method}/{size}/{rep}/results.pkl")
            with open(samples_address, 'rb') as handle:
                samples = pickle.load(handle)
            indices = samples[0]
        return indices

    def mark_sampled_points(self, adata: sc.AnnData, sample_ids, method: str) -> sc.AnnData:
        """Mark sampled points in the dataset. Accepts obs_names or row indices."""
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
        sc.pp.neighbors(adata, n_neighbors=self.N_NEIGHBORS)
        logger.info("Computing UMAP...")
        sc.tl.umap(adata)
        
        # Try to compute Leiden clustering if igraph is available
        try:
            logger.info("Computing Leiden clustering...")
            sc.tl.leiden(adata, resolution=self.LEIDEN_RESOLUTION)
        except ImportError:
            logger.warning("igraph package not found. Skipping Leiden clustering. Install with: conda install -c conda-forge python-igraph")
        except Exception as e:
            logger.warning(f"Could not compute Leiden clustering: {str(e)}")
            
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

    def generate_individual_umap(self, adata: sc.AnnData, method: str):
        """Generate UMAP plot for Figure 1 (MCC_01) with three square panels (sps, rare, random)."""
        filename = 'test_figure1_mcc_01_sampling_comparison.jpg' if self.test_mode else 'figure1_mcc_01_sampling_comparison.jpg'
        output_file = os.path.join(self.figures_dir, filename)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        fig, axs = plt.subplots(1, 3, figsize=(12, 6))
        for i, color_key in enumerate(self.COLORS_TO_PLOT):
            ax = axs[i]
            sc.pl.umap(
                adata,
                color=color_key,
                ax=ax,
                show=False,
                marker='s',
                size=self.POINT_SIZE,
                legend_loc='none',
            )
            ax.set_title(self.PANEL_TITLES[i], fontsize=self.PANEL_TITLE_FONT_SIZE)
            ax.set_xlabel('UMAP1', fontsize=self.AXIS_LABEL_FONT_SIZE)
            if i == 0:
                ax.set_ylabel('UMAP2', fontsize=self.AXIS_LABEL_FONT_SIZE)
            else:
                ax.set_ylabel('')
                ax.tick_params(axis='y', which='both', left=False, labelleft=False)
            ax.tick_params(axis='both', which='major', labelsize=10)
            ax.set_aspect('equal')
            
            categories = adata.obs[color_key].cat.categories
            palette = adata.uns[f'{color_key}_colors']
            handles = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=palette[j], label=str(cat)) for j, cat in enumerate(categories)]
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
        plt.tight_layout(rect=[0, 0.22, 1, 1])
        plt.savefig(output_file, dpi=350, bbox_inches='tight')
        plt.close()

    def generate_combined_umaps(self, adata: sc.AnnData, method: str):
        """Generate combined UMAP plots for Supplementary Figure 2."""
        filename = 'test_supp_figure2_mcc_01_sampling_analysis.jpg' if self.test_mode else 'supp_figure2_mcc_01_sampling_analysis.jpg'
        output_file = os.path.join(self.figures_dir, filename)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        colors = ['sps', 'random', 'rare', self.LABEL_KEY]
        fig, axes = plt.subplots(4, 1, figsize=(10, 18))
        axes = axes.flatten()
        for i, color in enumerate(colors):
            sc.pl.umap(
                adata,
                color=color,
                ax=axes[i],
                show=False,
                marker='s',
            )
        plt.tight_layout()
        plt.savefig(output_file, dpi=350)
        plt.close()

    def generate_label_umap(self, adata: sc.AnnData):
        """Generate UMAP plot for Supplementary Figure 2 (single plot)."""
        filename = 'test_supp_figure1_mcc_01_cell_types.jpg' if self.test_mode else 'supp_figure1_mcc_01_cell_types.jpg'
        output_file = os.path.join(self.figures_dir, filename)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        fig, ax = plt.subplots(figsize=(10, 10))

        sc.pl.umap(
            adata,
            color=self.LABEL_KEY,
            ax=ax,
            show=False,
            marker='s',
            size=self.POINT_SIZE,
            title=None,
        )

        ax.legend(
            fontsize=20,
            ncol=4,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.12),
            handletextpad=0.3,
            labelspacing=0.3,
        )

        ax.set_title(" ", fontsize=20, weight="bold")
        ax.set_aspect("equal")
        ax.tick_params(labelsize=16)

        fig.subplots_adjust(bottom=0.18)
        fig.savefig(output_file, bbox_inches="tight", pad_inches=0.1)
        plt.close(fig)

def get_data_path():
    """Get data path from environment variable or derive from script location."""
    project_root = os.environ.get('PROJECT_ROOT')
    if project_root is None:
        # Derive project root from script location
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(script_dir))))
    return os.path.join(project_root, 'data')

def main():
    parser = argparse.ArgumentParser(description='Generate UMAP visualizations for MCC_01 data')
    parser.add_argument('--test', action='store_true', help='Run in test mode with smaller dataset')
    args = parser.parse_args()
    
    ref = 5  # Default reference size
    method = 'sps'
    method2 = 'random'
    size = 100000
    rep = 0
    data_path = get_data_path()
    
    try:
        visualizer = UMAPVisualizer(data_path, test_mode=args.test)
        adata_ref = visualizer.load_reference_data(ref)
        adata_ref = visualizer.identify_rare_cell_types(adata_ref)
        
        indices = visualizer.load_sampling_results(ref, method, size, rep)
        indices2 = visualizer.load_sampling_results(ref, method2, size, rep)
        
        adata_ref = visualizer.mark_sampled_points(adata_ref, indices, 'sps')
        adata_ref = visualizer.mark_sampled_points(adata_ref, indices2, 'random')
        
        adata_ref = visualizer.compute_umap(adata_ref)
        adata_ref = visualizer.set_color_schemes(adata_ref)
        
        visualizer.generate_individual_umap(adata_ref, method)
        # visualizer.generate_combined_umaps(adata_ref, method)
        visualizer.generate_label_umap(adata_ref)
        
        logger.info("Visualization generation completed successfully")
    except Exception as e:
        logger.error(f"Error in main execution: {e}")
        raise

if __name__ == "__main__":
    main()

