#!/usr/bin/env python
# coding: utf-8

"""
Visualize and compare UMAP and Leiden clustering results on sampled data (sps vs random).
This script:
1. Loads reference data and sampling indices for sps and random methods
2. Subsets data to only sampled cells for each method
3. Computes UMAP on each sampled dataset
4. Computes Leiden clustering on each sampled dataset
5. Compares clustering quality using ground truth cell types
6. Generates visualization figures
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from typing import Tuple, Dict, Any
import logging
import argparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from scipy.stats import chi2_contingency

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=350, facecolor='white')


def get_project_root():
    """Get project root directory."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    return project_root


def get_data_path():
    """Get data path from environment or default."""
    project_root = get_project_root()
    data_path = os.path.join(project_root, 'data')
    return data_path


def get_output_dir():
    """Get output directory for results."""
    project_root = get_project_root()
    output_dir = os.path.join(project_root, 'analysis', 'results', 'sampling_clustering')
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def get_figures_dir():
    """Get figures directory."""
    project_root = get_project_root()
    figures_dir = os.path.join(project_root, 'figures', 'sampling_clustering')
    os.makedirs(figures_dir, exist_ok=True)
    return figures_dir


def load_reference_data(dataset: str, ref: int, data_path: str) -> sc.AnnData:
    """Load and preprocess reference data."""
    directory = f"{dataset}/benchmark"
    path = os.path.join(data_path, directory)
    address = os.path.join(path, f"{ref}/adata.h5ad")
    
    logger.info(f"Loading reference data from {address}")
    adata = sc.read_h5ad(address)
    
    # Ensure celltype is categorical
    if 'celltype' in adata.obs.columns:
        adata.obs['celltype'] = adata.obs['celltype'].astype('category')
    
    logger.info(f"Loaded data: {adata.shape[0]} cells, {adata.shape[1]} genes")
    return adata


def load_sampling_indices(dataset: str, ref: int, method: str, size: int, rep: int, data_path: str) -> np.ndarray:
    """Load sampling indices for a given configuration."""
    directory = f"{dataset}/benchmark"
    path = os.path.join(data_path, directory)
    
    try:
        if method == 'atomic':
            atomic_address = os.path.join(path, f"{ref}/atomic/{size}/{rep}/results.csv")
            if os.path.isfile(atomic_address):
                indices = pd.read_csv(atomic_address)['x'].values
                # Convert to int if possible, otherwise keep as string
                try:
                    indices = indices.astype(int)
                except:
                    pass
                return indices
            else:
                return None
        else:
            samples_address = os.path.join(path, f"{ref}/{method}/{size}/{rep}/results.pkl")
            if os.path.isfile(samples_address):
                with open(samples_address, 'rb') as handle:
                    samples = pickle.load(handle)
                # results.pkl contains (indices, elapsed_time)
                indices = samples[0]
                return indices
            else:
                logger.warning(f"File not found: {samples_address}")
                return None
    except Exception as e:
        logger.error(f"Error loading indices for {dataset}, ref {ref}, method {method}, size {size}, rep {rep}: {str(e)}")
        return None


def subset_sampled_data(adata: sc.AnnData, indices: np.ndarray) -> sc.AnnData:
    """Subset AnnData to only sampled cells."""
    # Handle both integer indices and obs_names (strings)
    if isinstance(indices[0], (str, np.str_)):
        # If indices are obs_names, use them directly
        adata_sampled = adata[adata.obs_names.isin(indices)].copy()
    else:
        # If indices are integer positions
        if max(indices) >= adata.shape[0]:
            logger.warning(f"Some indices ({max(indices)}) exceed data size ({adata.shape[0]}), filtering...")
            indices = indices[indices < adata.shape[0]]
        adata_sampled = adata[indices].copy()
    
    logger.info(f"Subsetted to {adata_sampled.shape[0]} sampled cells")
    return adata_sampled


def compute_umap_and_clustering(adata: sc.AnnData, n_neighbors: int = 15, resolution: float = 0.5) -> sc.AnnData:
    """Compute PCA, UMAP, and Leiden clustering on the data."""
    logger.info("Computing PCA...")
    sc.tl.pca(adata, svd_solver='arpack', n_comps=min(50, adata.shape[0]-1, adata.shape[1]-1))
    
    logger.info("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=min(n_neighbors, adata.shape[0]-1))
    
    logger.info("Computing UMAP...")
    sc.tl.umap(adata)
    
    logger.info("Computing Leiden clustering...")
    try:
        sc.tl.leiden(adata, resolution=resolution, key_added='leiden')
        logger.info(f"Found {len(adata.obs['leiden'].unique())} clusters")
    except Exception as e:
        logger.error(f"Error computing Leiden clustering: {str(e)}")
        raise
    
    return adata


def calculate_clustering_metrics(adata: sc.AnnData, label_key: str = 'celltype') -> Dict[str, float]:
    """Calculate clustering quality metrics using ground truth labels."""
    if label_key not in adata.obs.columns:
        logger.warning(f"Label key '{label_key}' not found in data")
        return {}
    
    # Get labels and clusters
    labels = adata.obs[label_key].values
    clusters = adata.obs['leiden'].values
    
    # Remove any NaN values
    valid_mask = ~(pd.isna(labels) | pd.isna(clusters))
    labels_clean = labels[valid_mask]
    clusters_clean = clusters[valid_mask]
    
    if len(labels_clean) == 0:
        logger.warning("No valid labels found for metric calculation")
        return {}
    
    # Calculate metrics
    metrics = {}
    
    # Adjusted Rand Index (ARI)
    try:
        metrics['ari'] = adjusted_rand_score(labels_clean, clusters_clean)
    except Exception as e:
        logger.warning(f"Error calculating ARI: {str(e)}")
        metrics['ari'] = np.nan
    
    # Normalized Mutual Information (NMI)
    try:
        metrics['nmi'] = normalized_mutual_info_score(labels_clean, clusters_clean)
    except Exception as e:
        logger.warning(f"Error calculating NMI: {str(e)}")
        metrics['nmi'] = np.nan
    
    # Number of clusters
    metrics['n_clusters'] = len(np.unique(clusters_clean))
    metrics['n_cell_types'] = len(np.unique(labels_clean))
    
    # Cluster purity (average max class per cluster)
    cluster_purities = []
    for cluster_id in np.unique(clusters_clean):
        cluster_mask = clusters_clean == cluster_id
        cluster_labels = labels_clean[cluster_mask]
        if len(cluster_labels) > 0:
            most_common = pd.Series(cluster_labels).value_counts().iloc[0]
            purity = most_common / len(cluster_labels)
            cluster_purities.append(purity)
    
    metrics['avg_cluster_purity'] = np.mean(cluster_purities) if cluster_purities else np.nan
    
    return metrics


def visualize_clustering_comparison(adata_sps: sc.AnnData, adata_random: sc.AnnData, 
                                   metrics_sps: Dict, metrics_random: Dict,
                                   dataset: str, size: int, rep: int, output_dir: str):
    """Create visualization comparing sps and random sampling clustering results."""
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
    
    # UMAP plots colored by cell type
    ax1 = fig.add_subplot(gs[0, 0])
    sc.pl.umap(adata_sps, color='celltype', ax=ax1, show=False, title='SPS - Cell Types', 
               legend_loc='none', frameon=False, size=20)
    ax1.set_xlabel('UMAP1', fontsize=12)
    ax1.set_ylabel('UMAP2', fontsize=12)
    
    ax2 = fig.add_subplot(gs[0, 1])
    sc.pl.umap(adata_random, color='celltype', ax=ax2, show=False, title='Random - Cell Types',
               legend_loc='none', frameon=False, size=20)
    ax2.set_xlabel('UMAP1', fontsize=12)
    ax2.set_ylabel('UMAP2', fontsize=12)
    
    # UMAP plots colored by Leiden clusters
    ax3 = fig.add_subplot(gs[0, 2])
    sc.pl.umap(adata_sps, color='leiden', ax=ax3, show=False, title='SPS - Leiden Clusters',
               legend_loc='none', frameon=False, size=20)
    ax3.set_xlabel('UMAP1', fontsize=12)
    ax3.set_ylabel('UMAP2', fontsize=12)
    
    ax4 = fig.add_subplot(gs[0, 3])
    sc.pl.umap(adata_random, color='leiden', ax=ax4, show=False, title='Random - Leiden Clusters',
               legend_loc='none', frameon=False, size=20)
    ax4.set_xlabel('UMAP1', fontsize=12)
    ax4.set_ylabel('UMAP2', fontsize=12)
    
    # Confusion matrices (cell type vs cluster)
    ax5 = fig.add_subplot(gs[1, 0:2])
    confusion_sps = pd.crosstab(adata_sps.obs['celltype'], adata_sps.obs['leiden'])
    sns.heatmap(confusion_sps, ax=ax5, cmap='Blues', cbar_kws={'label': 'Count'}, 
                fmt='d', annot=True, annot_kws={'size': 8})
    ax5.set_title('SPS: Cell Type vs Cluster', fontsize=12, fontweight='bold')
    ax5.set_xlabel('Leiden Cluster', fontsize=10)
    ax5.set_ylabel('Cell Type', fontsize=10)
    
    ax6 = fig.add_subplot(gs[1, 2:4])
    confusion_random = pd.crosstab(adata_random.obs['celltype'], adata_random.obs['leiden'])
    sns.heatmap(confusion_random, ax=ax6, cmap='Blues', cbar_kws={'label': 'Count'},
                fmt='d', annot=True, annot_kws={'size': 8})
    ax6.set_title('Random: Cell Type vs Cluster', fontsize=12, fontweight='bold')
    ax6.set_xlabel('Leiden Cluster', fontsize=10)
    ax6.set_ylabel('Cell Type', fontsize=10)
    
    # Metrics comparison bar plot
    ax7 = fig.add_subplot(gs[2, :])
    metrics_df = pd.DataFrame({
        'Method': ['SPS', 'Random'],
        'ARI': [metrics_sps.get('ari', np.nan), metrics_random.get('ari', np.nan)],
        'NMI': [metrics_sps.get('nmi', np.nan), metrics_random.get('nmi', np.nan)],
        'Avg Cluster Purity': [metrics_sps.get('avg_cluster_purity', np.nan), 
                               metrics_random.get('avg_cluster_purity', np.nan)]
    })
    
    x = np.arange(len(metrics_df))
    width = 0.25
    metrics_to_plot = ['ARI', 'NMI', 'Avg Cluster Purity']
    colors = ['#e74c3c', '#3498db', '#2ecc71']
    
    for i, (metric, color) in enumerate(zip(metrics_to_plot, colors)):
        offset = (i - 1) * width
        ax7.bar(x + offset, metrics_df[metric], width, label=metric, color=color, alpha=0.8)
    
    ax7.set_xlabel('Method', fontsize=12)
    ax7.set_ylabel('Score', fontsize=12)
    ax7.set_title('Clustering Quality Metrics Comparison', fontsize=14, fontweight='bold')
    ax7.set_xticks(x)
    ax7.set_xticklabels(metrics_df['Method'])
    ax7.legend()
    ax7.grid(axis='y', alpha=0.3)
    ax7.set_ylim([0, 1.1])
    
    # Add metric values as text
    for i, method in enumerate(metrics_df['Method']):
        for j, metric in enumerate(metrics_to_plot):
            offset = (j - 1) * width
            value = metrics_df.loc[i, metric]
            if not np.isnan(value):
                ax7.text(i + offset, value + 0.02, f'{value:.3f}', 
                        ha='center', va='bottom', fontsize=9)
    
    # Add text box with additional info
    info_text = f"Dataset: {dataset.upper()}\nSize: {size:,}\nRep: {rep}\n"
    info_text += f"SPS Clusters: {metrics_sps.get('n_clusters', 'N/A')}\n"
    info_text += f"Random Clusters: {metrics_random.get('n_clusters', 'N/A')}"
    fig.text(0.02, 0.98, info_text, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Save figure
    filename = f"{dataset}_size{size}_rep{rep}_clustering_comparison.png"
    output_file = os.path.join(output_dir, filename)
    plt.savefig(output_file, dpi=350, bbox_inches='tight')
    logger.info(f"Saved figure to {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Visualize and compare clustering on sampled data')
    parser.add_argument('--dataset', type=str, required=True, choices=['mcc', 'lcmv'],
                       help='Dataset name (mcc or lcmv)')
    parser.add_argument('--ref', type=int, required=True,
                       help='Reference size')
    parser.add_argument('--size', type=int, required=True,
                       help='Sample size')
    parser.add_argument('--rep', type=int, default=0,
                       help='Replication number (default: 0)')
    parser.add_argument('--n-neighbors', type=int, default=15,
                       help='Number of neighbors for UMAP (default: 15)')
    parser.add_argument('--resolution', type=float, default=0.5,
                       help='Leiden clustering resolution (default: 0.5)')
    
    args = parser.parse_args()
    
    # Get paths
    data_path = get_data_path()
    output_dir = get_output_dir()
    figures_dir = get_figures_dir()
    
    logger.info("="*80)
    logger.info(f"Clustering Comparison Analysis")
    logger.info(f"Dataset: {args.dataset}, Ref: {args.ref}, Size: {args.size}, Rep: {args.rep}")
    logger.info("="*80)
    
    # Load reference data
    adata_ref = load_reference_data(args.dataset, args.ref, data_path)
    
    # Load sampling indices
    logger.info("Loading sampling indices...")
    indices_sps = load_sampling_indices(args.dataset, args.ref, 'sps', args.size, args.rep, data_path)
    indices_random = load_sampling_indices(args.dataset, args.ref, 'random', args.size, args.rep, data_path)
    
    if indices_sps is None:
        logger.error("Failed to load SPS sampling indices")
        sys.exit(1)
    if indices_random is None:
        logger.error("Failed to load random sampling indices")
        sys.exit(1)
    
    logger.info(f"SPS sampled {len(indices_sps)} cells")
    logger.info(f"Random sampled {len(indices_random)} cells")
    
    # Subset to sampled cells
    logger.info("Subsetting to sampled cells...")
    adata_sps = subset_sampled_data(adata_ref, indices_sps)
    adata_random = subset_sampled_data(adata_ref, indices_random)
    
    # Compute UMAP and clustering for SPS
    logger.info("Computing UMAP and clustering for SPS data...")
    adata_sps = compute_umap_and_clustering(adata_sps, n_neighbors=args.n_neighbors, 
                                           resolution=args.resolution)
    
    # Compute UMAP and clustering for Random
    logger.info("Computing UMAP and clustering for Random data...")
    adata_random = compute_umap_and_clustering(adata_random, n_neighbors=args.n_neighbors,
                                              resolution=args.resolution)
    
    # Calculate metrics
    logger.info("Calculating clustering metrics...")
    metrics_sps = calculate_clustering_metrics(adata_sps)
    metrics_random = calculate_clustering_metrics(adata_random)
    
    logger.info("SPS Metrics:")
    for key, value in metrics_sps.items():
        logger.info(f"  {key}: {value}")
    
    logger.info("Random Metrics:")
    for key, value in metrics_random.items():
        logger.info(f"  {key}: {value}")
    
    # Save metrics to CSV
    metrics_df = pd.DataFrame({
        'method': ['sps', 'random'],
        'dataset': [args.dataset, args.dataset],
        'ref': [args.ref, args.ref],
        'size': [args.size, args.size],
        'rep': [args.rep, args.rep],
        'ari': [metrics_sps.get('ari', np.nan), metrics_random.get('ari', np.nan)],
        'nmi': [metrics_sps.get('nmi', np.nan), metrics_random.get('nmi', np.nan)],
        'avg_cluster_purity': [metrics_sps.get('avg_cluster_purity', np.nan),
                               metrics_random.get('avg_cluster_purity', np.nan)],
        'n_clusters': [metrics_sps.get('n_clusters', np.nan),
                      metrics_random.get('n_clusters', np.nan)],
        'n_cell_types': [metrics_sps.get('n_cell_types', np.nan),
                        metrics_random.get('n_cell_types', np.nan)]
    })
    
    metrics_file = os.path.join(output_dir, 
                               f"{args.dataset}_size{args.size}_rep{args.rep}_metrics.csv")
    metrics_df.to_csv(metrics_file, index=False)
    logger.info(f"Saved metrics to {metrics_file}")
    
    # Create visualization
    logger.info("Creating visualization...")
    visualize_clustering_comparison(adata_sps, adata_random, metrics_sps, metrics_random,
                                   args.dataset, args.size, args.rep, figures_dir)
    
    logger.info("Analysis completed successfully!")


if __name__ == "__main__":
    main()






