#!/usr/bin/env python
# coding: utf-8

"""
Common utilities for downstream experiments.

Provides shared data loading, metrics computation, and configuration
for comparing SPS vs Random sampling methods.

Supports both MCC and LCMV datasets.
"""

import os
import sys
import logging
import pickle
from datetime import datetime
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.spatial.distance import euclidean
import warnings

warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================

PROJECT_ROOT = os.environ.get(
    'PROJECT_ROOT',
    '/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility'
)

# MCC Dataset Configuration
DATA_PATH = os.path.join(PROJECT_ROOT, 'data', 'mcc', 'benchmark', '30')
DATA_PATH_MCC = DATA_PATH  # Alias for clarity

# LCMV Dataset Configuration
# LCMV data is in the sparseFlow_benchmarking directory
DATA_PATH_LCMV = '/fast/AG_Ohler/ekarimi/projects/sparseFlow_benchmarking/data/lcmv/benchmark/34'

# Analysis parameters
SAMPLE_SIZE = 100000
REP = 0
LABEL_KEY = 'celltype'

# MCC rare cell type (hardcoded for backward compatibility)
RARE_CELL_TYPE = 'osteoblast'

# Cell types in the MCC dataset
CELL_TYPES = [
    'osteoblast',
    'chondrocyte', 
    'fibroblast',
    'mesodermal cell',
    'lateral mesodermal cell'
]

# LCMV cell types (will be populated when dataset is loaded)
LCMV_CELL_TYPES = None  # Set dynamically


# ============================================================================
# Logging Setup
# ============================================================================

def setup_logger(name, log_dir=None):
    """
    Set up logging for the analysis.
    
    Parameters:
    -----------
    name : str
        Logger name (used in log filename)
    log_dir : str, optional
        Directory for log files. Defaults to PROJECT_ROOT/logs
        
    Returns:
    --------
    logging.Logger
    """
    if log_dir is None:
        log_dir = os.path.join(PROJECT_ROOT, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'{name}_{timestamp}.log')
    
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    # Clear existing handlers
    logger.handlers = []
    
    # File handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    
    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger


# ============================================================================
# Data Loading Functions
# ============================================================================

def load_mcc_dataset(logger=None):
    """
    Load the full MCC dataset.
    
    Parameters:
    -----------
    logger : logging.Logger, optional
        Logger for output messages
        
    Returns:
    --------
    AnnData : Full MCC dataset
    """
    adata_path = os.path.join(DATA_PATH, 'adata.h5ad')
    
    if logger:
        logger.info(f"Loading MCC dataset from {adata_path}")
    
    adata = sc.read_h5ad(adata_path)
    
    if logger:
        logger.info(f"Dataset shape: {adata.shape}")
        logger.info(f"Cell types: {adata.obs[LABEL_KEY].unique().tolist()}")
    
    return adata


def load_lcmv_dataset(logger=None):
    """
    Load the full LCMV dataset.
    
    Parameters:
    -----------
    logger : logging.Logger, optional
        Logger for output messages
        
    Returns:
    --------
    AnnData : Full LCMV dataset
    """
    adata_path = os.path.join(DATA_PATH_LCMV, 'adata.h5ad')
    
    if logger:
        logger.info(f"Loading LCMV dataset from {adata_path}")
    
    adata = sc.read_h5ad(adata_path)
    
    if logger:
        logger.info(f"Dataset shape: {adata.shape}")
        logger.info(f"Cell types: {adata.obs[LABEL_KEY].unique().tolist()}")
    
    return adata


def load_sampling_indices(method, logger=None, dataset='mcc'):
    """
    Load sampling indices for the specified method.
    
    Parameters:
    -----------
    method : str
        Sampling method ('sps' or 'random')
    logger : logging.Logger, optional
        Logger for output messages
    dataset : str
        Dataset name ('mcc' or 'lcmv')
        
    Returns:
    --------
    np.ndarray : Array of cell indices
    """
    if method == 'full':
        return None
    
    # Select data path based on dataset
    if dataset == 'lcmv':
        data_path = DATA_PATH_LCMV
    else:
        data_path = DATA_PATH
    
    base_path = os.path.join(data_path, method, str(SAMPLE_SIZE), str(REP))
    
    # Try loading from npy file first
    npy_path = os.path.join(base_path, 'indices.npy')
    pkl_path = os.path.join(base_path, 'results.pkl')
    
    if os.path.exists(npy_path):
        if logger:
            logger.info(f"Loading {method} indices from {npy_path}")
        indices = np.load(npy_path)
        if logger:
            logger.info(f"Loaded {len(indices)} {method} indices from npy file")
        return indices
    
    if os.path.exists(pkl_path):
        if logger:
            logger.info(f"Loading {method} indices from {pkl_path}")
        try:
            with open(pkl_path, 'rb') as f:
                data = pickle.load(f)
            
            # Handle different pickle formats
            if isinstance(data, tuple):
                indices = data[0]
            else:
                indices = data
            
            # Convert to numpy array if needed
            if isinstance(indices, list):
                indices = np.array(indices)
            
            if logger:
                logger.info(f"Loaded {len(indices)} {method} indices from pkl file")
            return indices
        except Exception as e:
            if logger:
                logger.error(f"Failed to load pickle file: {e}")
            raise
    
    raise FileNotFoundError(f"No indices file found at {base_path}")


def load_sampling_indices_lcmv(method, logger=None):
    """
    Load sampling indices for the specified method for LCMV dataset.
    
    Parameters:
    -----------
    method : str
        Sampling method ('sps' or 'random')
    logger : logging.Logger, optional
        Logger for output messages
        
    Returns:
    --------
    np.ndarray : Array of cell indices
    """
    return load_sampling_indices(method, logger, dataset='lcmv')


def get_cell_type_distribution(adata, indices=None, logger=None):
    """
    Get cell type distribution for the dataset or a subset.
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    indices : np.ndarray, optional
        Cell indices for subset. If None, uses full dataset.
    logger : logging.Logger, optional
        Logger for output messages
        
    Returns:
    --------
    pd.Series : Cell type counts
    """
    if indices is not None:
        cell_types = adata.obs[LABEL_KEY].iloc[indices]
    else:
        cell_types = adata.obs[LABEL_KEY]
    
    counts = cell_types.value_counts()
    
    if logger:
        logger.info("Cell type distribution:")
        for ct, count in counts.items():
            pct = 100 * count / len(cell_types)
            logger.info(f"  {ct}: {count:,} cells ({pct:.2f}%)")
    
    return counts


def compute_osteoblast_statistics(adata, sps_indices, random_indices, logger=None):
    """
    Compute osteoblast cell statistics for SPS vs Random samples.
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    sps_indices : np.ndarray
        SPS sampling indices
    random_indices : np.ndarray
        Random sampling indices
    logger : logging.Logger, optional
        Logger for output messages
        
    Returns:
    --------
    dict : Statistics for both methods
    """
    # Get cell types
    full_labels = adata.obs[LABEL_KEY].values
    sps_labels = full_labels[sps_indices]
    random_labels = full_labels[random_indices]
    
    # Count osteoblasts
    full_osteoblast = np.sum(full_labels == RARE_CELL_TYPE)
    sps_osteoblast = np.sum(sps_labels == RARE_CELL_TYPE)
    random_osteoblast = np.sum(random_labels == RARE_CELL_TYPE)
    
    stats = {
        'full': {
            'total': len(full_labels),
            'osteoblast_count': full_osteoblast,
            'osteoblast_pct': 100 * full_osteoblast / len(full_labels)
        },
        'sps': {
            'total': len(sps_indices),
            'osteoblast_count': sps_osteoblast,
            'osteoblast_pct': 100 * sps_osteoblast / len(sps_indices),
            'enrichment_vs_full': (sps_osteoblast / len(sps_indices)) / (full_osteoblast / len(full_labels))
        },
        'random': {
            'total': len(random_indices),
            'osteoblast_count': random_osteoblast,
            'osteoblast_pct': 100 * random_osteoblast / len(random_indices),
            'enrichment_vs_full': (random_osteoblast / len(random_indices)) / (full_osteoblast / len(full_labels))
        }
    }
    
    if logger:
        logger.info("\nOsteoblast Statistics:")
        logger.info(f"  Full dataset: {full_osteoblast:,} cells ({stats['full']['osteoblast_pct']:.2f}%)")
        logger.info(f"  SPS sample:   {sps_osteoblast:,} cells ({stats['sps']['osteoblast_pct']:.2f}%) "
                   f"[{stats['sps']['enrichment_vs_full']:.2f}x enrichment]")
        logger.info(f"  Random sample: {random_osteoblast:,} cells ({stats['random']['osteoblast_pct']:.2f}%) "
                   f"[{stats['random']['enrichment_vs_full']:.2f}x enrichment]")
    
    return stats


# ============================================================================
# Rare Cell Type Identification Functions
# ============================================================================

def calculate_cell_type_distances(adata, label_key='celltype'):
    """
    Calculate the distance between the mean of each cell type and the mean of the whole dataset.
    
    Parameters:
    -----------
    adata : AnnData
        Dataset
    label_key : str
        Column name for cell type labels
        
    Returns:
    --------
    pd.DataFrame : DataFrame with cell_type, distance, n_cells columns
    """
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


def identify_rare_cell_types(adata, label_key='celltype', 
                             distance_percentile=75, 
                             frequency_threshold_pct=1.0,
                             logger=None):
    """
    Identify rare cell types using distance + frequency criteria.
    
    Rare cell types are defined as:
    - Distance from dataset mean > specified percentile (default 75th)
    - AND frequency < specified percentage of total cells (default 1%)
    
    Parameters:
    -----------
    adata : AnnData
        Dataset
    label_key : str
        Column name for cell type labels
    distance_percentile : int
        Top percentile for distance (e.g., 75 = top 25%)
    frequency_threshold_pct : float
        Maximum frequency as percentage of total (e.g., 1.0 = 1%)
    logger : logging.Logger, optional
        Logger for output messages
        
    Returns:
    --------
    tuple : (rare_cell_types_list, cell_type_df)
        - rare_cell_types_list: List of rare cell type names
        - cell_type_df: DataFrame with all cell type statistics
    """
    if logger:
        logger.info(f"Identifying rare cell types (distance > {distance_percentile}th pct, freq < {frequency_threshold_pct}%)...")
    
    # Calculate cell type distances and frequencies
    cell_type_df = calculate_cell_type_distances(adata, label_key=label_key)
    total_cells = adata.n_obs
    cell_type_df['frequency_pct'] = (cell_type_df['n_cells'] / total_cells) * 100
    
    # Calculate thresholds
    distance_threshold = np.percentile(cell_type_df['distance'], distance_percentile)
    frequency_threshold = total_cells * (frequency_threshold_pct / 100.0)
    
    # Identify rare cell types
    rare_mask = (cell_type_df['distance'] > distance_threshold) & \
                (cell_type_df['n_cells'] < frequency_threshold)
    rare_cell_types = cell_type_df[rare_mask]['cell_type'].tolist()
    
    if logger:
        logger.info(f"  Distance threshold (>{distance_percentile}th pct): {distance_threshold:.2f}")
        logger.info(f"  Frequency threshold (<{frequency_threshold_pct}%): {frequency_threshold:.0f} cells")
        logger.info(f"  Found {len(rare_cell_types)} rare cell types: {rare_cell_types}")
        
        logger.info("\n  Cell type statistics:")
        for _, row in cell_type_df.sort_values('frequency_pct').iterrows():
            is_rare = 'RARE' if row['cell_type'] in rare_cell_types else ''
            logger.info(f"    {row['cell_type']}: {row['n_cells']:,} cells ({row['frequency_pct']:.2f}%), "
                       f"distance={row['distance']:.2f} {is_rare}")
    
    return rare_cell_types, cell_type_df


def compute_rare_cell_type_statistics(adata, sps_indices, random_indices, 
                                       rare_cell_types, logger=None):
    """
    Compute statistics for rare cell types for SPS vs Random samples.
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    sps_indices : np.ndarray
        SPS sampling indices
    random_indices : np.ndarray
        Random sampling indices
    rare_cell_types : list
        List of rare cell type names
    logger : logging.Logger, optional
        Logger for output messages
        
    Returns:
    --------
    dict : Statistics for each rare cell type
    """
    full_labels = adata.obs[LABEL_KEY].values
    sps_labels = full_labels[sps_indices]
    random_labels = full_labels[random_indices]
    
    stats = {}
    
    for rare_ct in rare_cell_types:
        full_count = np.sum(full_labels == rare_ct)
        sps_count = np.sum(sps_labels == rare_ct)
        random_count = np.sum(random_labels == rare_ct)
        
        full_pct = 100 * full_count / len(full_labels)
        sps_pct = 100 * sps_count / len(sps_indices)
        random_pct = 100 * random_count / len(random_indices)
        
        # Avoid division by zero
        if full_pct > 0:
            sps_enrichment = sps_pct / full_pct
            random_enrichment = random_pct / full_pct
        else:
            sps_enrichment = 0
            random_enrichment = 0
        
        stats[rare_ct] = {
            'full': {
                'count': full_count,
                'pct': full_pct
            },
            'sps': {
                'count': sps_count,
                'pct': sps_pct,
                'enrichment_vs_full': sps_enrichment
            },
            'random': {
                'count': random_count,
                'pct': random_pct,
                'enrichment_vs_full': random_enrichment
            }
        }
        
        if logger:
            logger.info(f"\n{rare_ct} Statistics:")
            logger.info(f"  Full dataset: {full_count:,} cells ({full_pct:.3f}%)")
            logger.info(f"  SPS sample:   {sps_count:,} cells ({sps_pct:.3f}%) [{sps_enrichment:.2f}x]")
            logger.info(f"  Random sample: {random_count:,} cells ({random_pct:.3f}%) [{random_enrichment:.2f}x]")
    
    return stats


def subset_to_indices(adata, indices, copy=True):
    """
    Subset AnnData to specified indices.
    
    Parameters:
    -----------
    adata : AnnData
        Full dataset
    indices : np.ndarray
        Cell indices to keep
    copy : bool
        Whether to return a copy
        
    Returns:
    --------
    AnnData : Subsetted dataset
    """
    if indices is None:
        return adata.copy() if copy else adata
    return adata[indices].copy() if copy else adata[indices]


def load_marker_gene_results(method, de_method='wilcoxon', dataset='mcc'):
    """
    Load marker gene results from the scanpy_marker_genes analysis.
    
    Parameters:
    -----------
    method : str
        Sampling method ('full', 'sps', or 'random')
    de_method : str
        DE method used ('wilcoxon', 't-test', or 'logreg')
    dataset : str
        Dataset name ('mcc' or 'lcmv')
        
    Returns:
    --------
    pd.DataFrame : Marker gene results
    """
    if dataset == 'lcmv':
        results_dir = os.path.join(
            PROJECT_ROOT, 'analysis', 'scanpy_marker_genes', 'results_lcmv'
        )
    else:
        results_dir = os.path.join(
            PROJECT_ROOT, 'analysis', 'scanpy_marker_genes', 'results'
        )
    
    csv_path = os.path.join(results_dir, f'{method}_{de_method}_marker_genes.csv')
    
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Marker gene results not found: {csv_path}")
    
    return pd.read_csv(csv_path)


def load_marker_gene_results_lcmv(method, de_method='wilcoxon'):
    """
    Load marker gene results for LCMV dataset.
    
    Parameters:
    -----------
    method : str
        Sampling method ('full', 'sps', or 'random')
    de_method : str
        DE method used ('wilcoxon', 't-test', or 'logreg')
        
    Returns:
    --------
    pd.DataFrame : Marker gene results
    """
    return load_marker_gene_results(method, de_method, dataset='lcmv')


def get_rare_cell_marker_genes(method, de_method='wilcoxon', cell_type='osteoblast', 
                                top_n=50, dataset='mcc'):
    """
    Get marker genes for a specific rare cell type.
    
    Parameters:
    -----------
    method : str
        Sampling method ('full', 'sps', or 'random')
    de_method : str
        DE method used
    cell_type : str
        Cell type to get markers for
    top_n : int
        Number of top markers to return
    dataset : str
        Dataset name ('mcc' or 'lcmv')
        
    Returns:
    --------
    pd.DataFrame : Top marker genes for the cell type
    """
    df = load_marker_gene_results(method, de_method, dataset=dataset)
    markers = df[df['cell_type'] == cell_type].head(top_n)
    return markers


def get_lcmv_cell_type_distribution(adata=None, indices=None, logger=None):
    """
    Get cell type distribution for the LCMV dataset or a subset.
    
    Parameters:
    -----------
    adata : AnnData, optional
        Full dataset (will be loaded if not provided)
    indices : np.ndarray, optional
        Cell indices for subset. If None, uses full dataset.
    logger : logging.Logger, optional
        Logger for output messages
        
    Returns:
    --------
    pd.Series : Cell type counts
    """
    if adata is None:
        adata = load_lcmv_dataset(logger)
    
    if indices is not None:
        cell_types = adata.obs[LABEL_KEY].iloc[indices]
    else:
        cell_types = adata.obs[LABEL_KEY]
    
    counts = cell_types.value_counts()
    
    if logger:
        logger.info("LCMV Cell type distribution:")
        for ct, count in counts.items():
            pct = 100 * count / len(cell_types)
            logger.info(f"  {ct}: {count:,} cells ({pct:.2f}%)")
    
    return counts


if __name__ == '__main__':
    # Test the utilities
    logger = setup_logger('test_utils')
    
    logger.info("Testing common utilities...")
    
    # Test MCC data loading
    try:
        adata = load_mcc_dataset(logger)
        logger.info(f"Successfully loaded MCC dataset with shape {adata.shape}")
    except Exception as e:
        logger.error(f"Failed to load MCC dataset: {e}")
    
    # Test MCC index loading
    for method in ['sps', 'random']:
        try:
            indices = load_sampling_indices(method, logger, dataset='mcc')
            logger.info(f"Successfully loaded MCC {method} indices: {len(indices)}")
        except Exception as e:
            logger.warning(f"Could not load MCC {method} indices: {e}")
    
    # Test LCMV data loading
    try:
        adata_lcmv = load_lcmv_dataset(logger)
        logger.info(f"Successfully loaded LCMV dataset with shape {adata_lcmv.shape}")
        
        # Test rare cell type identification
        rare_types, cell_type_df = identify_rare_cell_types(adata_lcmv, logger=logger)
        logger.info(f"Identified {len(rare_types)} rare cell types in LCMV: {rare_types}")
        
    except Exception as e:
        logger.error(f"Failed to load LCMV dataset: {e}")
    
    # Test LCMV index loading
    for method in ['sps', 'random']:
        try:
            indices = load_sampling_indices(method, logger, dataset='lcmv')
            logger.info(f"Successfully loaded LCMV {method} indices: {len(indices)}")
        except Exception as e:
            logger.warning(f"Could not load LCMV {method} indices: {e}")
    
    logger.info("Utility test completed.")
