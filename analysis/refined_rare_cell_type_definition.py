import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.spatial.distance import euclidean

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
LABEL = "celltype"
REFS = {"lcmv": 34, "mcc": 30, "mcc_01": 30, "mcc_05": 30}


def calculate_cell_type_distances(adata, label_key=LABEL):
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    dataset_mean = np.mean(X, axis=0)
    types = adata.obs[label_key].values
    rows = []
    for ct in np.unique(types):
        mask = types == ct
        rows.append({
            "cell_type": ct,
            "distance": euclidean(np.mean(X[mask], axis=0), dataset_mean),
            "n_cells": int(mask.sum()),
        })
    return pd.DataFrame(rows)


def identify_rare_cell_types_distance_and_frequency(
    cell_type_df, total_cells, distance_percentile=75, frequency_threshold_pct=1.0,
):
    dist_thresh = np.percentile(cell_type_df["distance"], distance_percentile)
    freq_thresh = total_cells * (frequency_threshold_pct / 100.0)
    mask = (cell_type_df["distance"] > dist_thresh) & (cell_type_df["n_cells"] < freq_thresh)
    return cell_type_df[mask]["cell_type"].tolist(), dist_thresh, freq_thresh


def main():
    os.makedirs(RES, exist_ok=True)
    summaries = []

    for dataset, ref in [("mcc_01", 30), ("mcc_05", 30), ("mcc", 30), ("lcmv", 34)]:
        adata = sc.read_h5ad(os.path.join(DATA, f"{dataset}/benchmark/{ref}/adata.h5ad"))
        total = adata.shape[0]
        df = calculate_cell_type_distances(adata)
        df["frequency_pct"] = df["n_cells"] / total * 100
        df = df.sort_values("n_cells")

        rare, dist_thresh, freq_thresh = identify_rare_cell_types_distance_and_frequency(df, total)
        df["rare_distance_and_frequency"] = df["cell_type"].isin(rare)
        df["dataset"] = dataset
        df.to_csv(os.path.join(RES, f"{dataset}_rare_cell_types_analysis.csv"), index=False)

        summaries.append({
            "dataset": dataset,
            "n_rare_types": len(rare),
            "rare_types": ", ".join(rare),
            "distance_threshold": dist_thresh,
            "frequency_threshold": freq_thresh,
            "frequency_threshold_pct": freq_thresh / total * 100,
        })

    pd.DataFrame(summaries).to_csv(os.path.join(RES, "rare_cell_type_definitions_summary.csv"), index=False)


if __name__ == "__main__":
    main()
