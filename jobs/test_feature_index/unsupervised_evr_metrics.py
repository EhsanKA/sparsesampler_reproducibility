#!/usr/bin/env python
"""
Develop unsupervised metrics to evaluate EVR index quality for rare cell capture.

Hypothesis: A good EVR index should produce samples that:
1. Cover more of the feature space (higher diversity)
2. Include cells far from the dataset centroid (outliers = rare cells)
3. Have higher entropy/spread in the sample

Proposed Metrics:
1. Centroid Distance Score (CDS): Mean distance of sampled cells from dataset centroid
2. Diversity Score (DS): Mean pairwise distance or k-NN distance
3. Outlier Inclusion Ratio (OIR): Fraction of dataset outliers included in sample
4. Feature Space Coverage (FSC): Variance ratio of sample vs full dataset

We'll evaluate how well these correlate with actual rare cell capture.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist
from scipy.stats import pearsonr, spearmanr
import pickle

script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
data_path = os.path.join(project_root, 'data')

# Add analysis directory
analysis_dir = os.path.join(project_root, 'analysis')
if analysis_dir not in sys.path:
    sys.path.insert(0, analysis_dir)

# Total rare cells per reference
TOTAL_RARE_CELLS = {
    'lcmv': {'1M': 20226, '5M': 102561, '10M': 205573, '20M': 410544, '34M': 705738},
    'mcc': {'5M': 4750, '10M': 9399, '20M': 18884, '25M': 23891, '30M': 30000},
    'mcc_05': {'5M': 2396, '10M': 4652, '20M': 9490, '25M': 11984, '30M': 15000},
    'mcc_01': {'5M': 458, '10M': 929, '20M': 1912, '25M': 2377, '30M': 3000},
}

DATASETS = {
    'lcmv': {'refs': [34], 'label': 'LCMV'},  # Use largest reference for analysis
    'mcc': {'refs': [30], 'label': 'MCC (1%)'},
}


def load_pca_embeddings(dataset, ref):
    """Load or compute PCA embeddings for a dataset."""
    import scanpy as sc
    
    benchmark_path = os.path.join(data_path, f"{dataset}/benchmark")
    adata_path = os.path.join(benchmark_path, f"{ref}/adata.h5ad")
    
    if not os.path.exists(adata_path):
        print(f"adata not found: {adata_path}")
        return None, None
    
    print(f"Loading {dataset} ref {ref}...")
    adata = sc.read_h5ad(adata_path)
    
    # Check if PCA exists, otherwise compute
    if 'X_pca' not in adata.obsm:
        print("  Computing PCA...")
        sc.pp.pca(adata, n_comps=50)
    
    return adata.obsm['X_pca'], adata.obs


def load_sps_indices(dataset, ref, size=100000, rep=0):
    """Load SPS sampling indices."""
    benchmark_path = os.path.join(data_path, f"{dataset}/benchmark")
    results_path = os.path.join(benchmark_path, f"{ref}/sps/{size}/{rep}/results.pkl")
    
    if not os.path.exists(results_path):
        print(f"SPS results not found: {results_path}")
        return None
    
    with open(results_path, 'rb') as f:
        results = pickle.load(f)
    
    return results[0]  # indices


def compute_centroid_distance_score(X_pca, sample_indices, n_components=None):
    """
    Centroid Distance Score (CDS): Mean distance of sampled cells from dataset centroid.
    
    Higher CDS = sample includes more cells far from the average (outliers/rare cells).
    
    Parameters:
    -----------
    X_pca : array-like
        PCA embeddings of all cells (n_cells, n_components)
    sample_indices : array-like
        Indices of sampled cells
    n_components : int, optional
        Number of PCA components to use. If None, use all.
    
    Returns:
    --------
    float: Mean distance from centroid
    """
    if n_components is not None:
        X = X_pca[:, :n_components]
    else:
        X = X_pca
    
    # Compute dataset centroid
    centroid = np.mean(X, axis=0)
    
    # Get sampled cells
    X_sample = X[sample_indices]
    
    # Compute distances from centroid
    distances = np.linalg.norm(X_sample - centroid, axis=1)
    
    return np.mean(distances)


def compute_diversity_score(X_pca, sample_indices, n_components=None, subsample=1000):
    """
    Diversity Score (DS): Mean pairwise distance between sampled cells.
    
    Higher DS = more diverse sample covering more of the feature space.
    
    Parameters:
    -----------
    X_pca : array-like
        PCA embeddings
    sample_indices : array-like
        Indices of sampled cells
    n_components : int, optional
        Number of PCA components to use
    subsample : int
        Number of cells to subsample for efficiency
    
    Returns:
    --------
    float: Mean pairwise distance
    """
    if n_components is not None:
        X = X_pca[:, :n_components]
    else:
        X = X_pca
    
    X_sample = X[sample_indices]
    
    # Subsample for computational efficiency
    if len(X_sample) > subsample:
        idx = np.random.choice(len(X_sample), subsample, replace=False)
        X_sample = X_sample[idx]
    
    # Compute mean pairwise distance
    distances = pdist(X_sample)
    return np.mean(distances)


def compute_outlier_inclusion_ratio(X_pca, sample_indices, n_components=None, outlier_percentile=95):
    """
    Outlier Inclusion Ratio (OIR): Fraction of dataset outliers included in sample.
    
    Outliers defined as cells with distance from centroid > percentile threshold.
    
    Parameters:
    -----------
    X_pca : array-like
        PCA embeddings
    sample_indices : array-like
        Indices of sampled cells
    n_components : int, optional
        Number of PCA components to use
    outlier_percentile : float
        Percentile threshold for defining outliers (e.g., 95 = top 5%)
    
    Returns:
    --------
    float: Fraction of outliers in sample vs expected by random sampling
    """
    if n_components is not None:
        X = X_pca[:, :n_components]
    else:
        X = X_pca
    
    # Compute centroid and distances
    centroid = np.mean(X, axis=0)
    all_distances = np.linalg.norm(X - centroid, axis=1)
    
    # Define outliers
    threshold = np.percentile(all_distances, outlier_percentile)
    outlier_mask = all_distances > threshold
    n_outliers_total = np.sum(outlier_mask)
    
    # Count outliers in sample
    sample_set = set(sample_indices)
    outlier_indices = np.where(outlier_mask)[0]
    n_outliers_in_sample = len(set(outlier_indices) & sample_set)
    
    # Compute enrichment ratio
    expected_fraction = len(sample_indices) / len(X)
    expected_outliers = expected_fraction * n_outliers_total
    
    if expected_outliers > 0:
        enrichment = n_outliers_in_sample / expected_outliers
    else:
        enrichment = 0
    
    return enrichment


def compute_variance_coverage(X_pca, sample_indices, n_components=None):
    """
    Feature Space Coverage (FSC): Ratio of sample variance to full dataset variance.
    
    Higher FSC = sample captures more of the dataset's variability.
    
    Parameters:
    -----------
    X_pca : array-like
        PCA embeddings
    sample_indices : array-like
        Indices of sampled cells
    n_components : int, optional
        Number of PCA components to use
    
    Returns:
    --------
    float: Variance ratio
    """
    if n_components is not None:
        X = X_pca[:, :n_components]
    else:
        X = X_pca
    
    # Full dataset variance (sum across all components)
    full_variance = np.sum(np.var(X, axis=0))
    
    # Sample variance
    X_sample = X[sample_indices]
    sample_variance = np.sum(np.var(X_sample, axis=0))
    
    return sample_variance / full_variance if full_variance > 0 else 0


def evaluate_evr_indices():
    """Evaluate all proposed metrics across EVR indices."""
    
    print("=" * 80)
    print("UNSUPERVISED EVR INDEX METRICS EVALUATION")
    print("=" * 80)
    
    # Load EVR tables for comparison
    tables_dir = os.path.join(script_dir, 'tables')
    
    results = []
    
    for dataset, cfg in DATASETS.items():
        for ref in cfg['refs']:
            print(f"\n{'='*60}")
            print(f"Dataset: {cfg['label']}, Reference: {ref}M")
            print("=" * 60)
            
            # Load PCA embeddings
            X_pca, obs = load_pca_embeddings(dataset, ref)
            if X_pca is None:
                continue
            
            print(f"Loaded {X_pca.shape[0]} cells with {X_pca.shape[1]} PCs")
            
            # Load SPS indices
            sps_indices = load_sps_indices(dataset, ref)
            if sps_indices is None:
                continue
            
            print(f"SPS sample size: {len(sps_indices)}")
            
            # Load EVR table
            evr_table_path = os.path.join(tables_dir, f'{dataset}_feature_index_table_size_100000.csv')
            if not os.path.exists(evr_table_path):
                print(f"EVR table not found: {evr_table_path}")
                continue
            
            evr_df = pd.read_csv(evr_table_path)
            
            # For each EVR index, compute metrics
            print(f"\n{'EVR':<6} {'CDS':>10} {'DS':>10} {'OIR':>10} {'FSC':>10} {'RareCapt':>12}")
            print("-" * 65)
            
            for evr_idx in range(1, 31):
                if evr_idx not in evr_df['feature_index'].values:
                    continue
                
                # Compute metrics using this many PCA components
                n_comp = evr_idx
                
                cds = compute_centroid_distance_score(X_pca, sps_indices, n_components=n_comp)
                ds = compute_diversity_score(X_pca, sps_indices, n_components=n_comp)
                oir = compute_outlier_inclusion_ratio(X_pca, sps_indices, n_components=n_comp)
                fsc = compute_variance_coverage(X_pca, sps_indices, n_components=n_comp)
                
                # Get actual rare cell capture from EVR table
                ref_col = f'{ref}M'
                if ref_col in evr_df.columns:
                    rare_captured = evr_df[evr_df['feature_index'] == evr_idx][ref_col].values[0]
                    total_rare = TOTAL_RARE_CELLS[dataset].get(ref_col, 1)
                    rare_pct = (rare_captured / total_rare) * 100
                else:
                    rare_pct = 0
                
                print(f"{evr_idx:<6} {cds:>10.2f} {ds:>10.2f} {oir:>10.2f} {fsc:>10.3f} {rare_pct:>11.1f}%")
                
                results.append({
                    'dataset': dataset,
                    'ref': ref,
                    'evr_index': evr_idx,
                    'centroid_distance_score': cds,
                    'diversity_score': ds,
                    'outlier_inclusion_ratio': oir,
                    'variance_coverage': fsc,
                    'rare_cell_pct': rare_pct
                })
    
    results_df = pd.DataFrame(results)
    
    # Compute correlations
    print("\n" + "=" * 80)
    print("CORRELATION WITH RARE CELL CAPTURE")
    print("=" * 80)
    
    metrics = ['centroid_distance_score', 'diversity_score', 'outlier_inclusion_ratio', 'variance_coverage']
    
    print(f"\n{'Metric':<30} {'Pearson r':>12} {'Spearman ρ':>12} {'p-value':>12}")
    print("-" * 70)
    
    correlations = {}
    for metric in metrics:
        if len(results_df) > 2:
            pearson_r, pearson_p = pearsonr(results_df[metric], results_df['rare_cell_pct'])
            spearman_r, spearman_p = spearmanr(results_df[metric], results_df['rare_cell_pct'])
            correlations[metric] = {'pearson': pearson_r, 'spearman': spearman_r}
            print(f"{metric:<30} {pearson_r:>12.3f} {spearman_r:>12.3f} {spearman_p:>12.4f}")
    
    # Recommend best metric
    print("\n" + "=" * 80)
    print("RECOMMENDATION")
    print("=" * 80)
    
    if correlations:
        best_metric = max(correlations.keys(), key=lambda x: abs(correlations[x]['spearman']))
        print(f"\nBest unsupervised metric: {best_metric}")
        print(f"  Spearman correlation with rare cell capture: {correlations[best_metric]['spearman']:.3f}")
        print(f"  Pearson correlation with rare cell capture: {correlations[best_metric]['pearson']:.3f}")
    
    return results_df, correlations


if __name__ == "__main__":
    results_df, correlations = evaluate_evr_indices()
    
    # Save results
    output_dir = os.path.join(script_dir, 'analysis_results')
    os.makedirs(output_dir, exist_ok=True)
    
    results_df.to_csv(os.path.join(output_dir, 'unsupervised_evr_metrics.csv'), index=False)
    print(f"\nResults saved to {output_dir}/unsupervised_evr_metrics.csv")
