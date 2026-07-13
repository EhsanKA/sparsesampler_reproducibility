import os
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
LABEL = "celltype"
FI = 12
REFS = {"lcmv": 34, "mcc": 30}


def reproduce_binning_with_coordinates(X, feature_index=FI):
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X)
    n_components = min(X_std.shape[1], feature_index + 1)
    pca = PCA(n_components=n_components)
    pca.fit(X_std)
    X_pca = pca.transform(X_std)
    del X_std

    explained = pca.explained_variance_ratio_
    k = (1.0 / explained[feature_index]) * 2
    out = np.ceil(explained * k).astype(int)
    threshold = 2
    filtered = out[out > threshold]
    while filtered.size == 0 and threshold >= out.min():
        threshold -= 1
        filtered = out[out > threshold]

    num_pcs = min(filtered.shape[0], X_pca.shape[1])
    bin_coords = np.zeros((X_pca.shape[0], num_pcs), dtype=np.int64)
    for i in range(num_pcs):
        col = X_pca[:, i]
        edges = np.linspace(col.min(), col.max(), filtered[i] + 1)
        bin_coords[:, i] = np.digitize(col, edges)
    return bin_coords, filtered[:num_pcs], num_pcs


def coords_to_id(coords, bins_per_pc):
    cell_id = coords[:, 0].copy().astype(np.int64)
    mult = int(bins_per_pc[0]) + 2
    for i in range(1, coords.shape[1]):
        cell_id += coords[:, i].astype(np.int64) * mult
        mult *= (int(bins_per_pc[i]) + 2)
    return cell_id


def compute_adjacent_stats(bin_coords, bins_per_pc, celltypes):
    n_cells, num_pcs = bin_coords.shape
    grid_ids = coords_to_id(bin_coords, bins_per_pc)
    unique_ids, inverse, counts = np.unique(grid_ids, return_inverse=True, return_counts=True)
    cell_bin_size = counts[inverse]
    id_to_count = dict(zip(unique_ids, counts))

    one_cell_indices = np.where(cell_bin_size == 1)[0]
    if len(one_cell_indices) == 0:
        return pd.DataFrame()

    id_to_cells = {}
    for idx in range(n_cells):
        gid = grid_ids[idx]
        id_to_cells.setdefault(gid, []).append(idx)

    offsets = np.zeros((2 * num_pcs, num_pcs), dtype=np.int64)
    for d in range(num_pcs):
        offsets[2 * d, d] = -1
        offsets[2 * d + 1, d] = 1

    rows = []
    for cell_idx in one_cell_indices:
        coord = bin_coords[cell_idx]
        ct = celltypes[cell_idx]
        n_adjacent_cells = 0
        n_occupied_neighbors = 0
        neighbor_celltypes = []

        for offset in offsets:
            neighbor_coord = coord + offset
            if np.any(neighbor_coord < 0):
                continue
            neighbor_id = coords_to_id(neighbor_coord.reshape(1, -1), bins_per_pc)[0]
            if neighbor_id in id_to_count:
                nc = id_to_count[neighbor_id]
                n_adjacent_cells += nc
                n_occupied_neighbors += 1
                for cidx in id_to_cells[neighbor_id]:
                    neighbor_celltypes.append(celltypes[cidx])

        same_type_neighbors = sum(1 for nct in neighbor_celltypes if nct == ct)
        rows.append({
            "cell_index": cell_idx,
            "celltype": ct,
            "n_adjacent_cells": n_adjacent_cells,
            "n_occupied_neighbor_bins": n_occupied_neighbors,
            "n_same_type_neighbors": same_type_neighbors,
            "frac_same_type": same_type_neighbors / n_adjacent_cells if n_adjacent_cells else 0.0,
            "is_isolated": n_adjacent_cells == 0,
        })
    return pd.DataFrame(rows)


def plot_distribution(df, dataset_name, fig_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    vals = df["n_adjacent_cells"].values
    x_cap = max(int(np.percentile(vals, 95)), 1)
    bins = np.linspace(-0.5, x_cap + 0.5, min(50, x_cap + 1) + 1)
    ax.hist(vals[vals <= x_cap], bins=bins, color="#4C72B0", edgecolor="white", linewidth=0.5)
    ax.set_xlim(-0.5, x_cap + 0.5)
    ax.set_xlabel("Total cells in adjacent bins")
    ax.set_ylabel("Number of 1-cell bins")
    ax.set_title(f"{dataset_name.upper()}: Cells in neighborhood")
    n_isolated = (vals == 0).sum()
    n_truncated = (vals > x_cap).sum()
    note = f"Isolated (0 neighbors): {n_isolated} ({n_isolated/len(vals):.1%})"
    if n_truncated > 0:
        note += f"\n>{x_cap} neighbors: {n_truncated} ({n_truncated/len(vals):.1%}) omitted"
    ax.annotate(
        note, xy=(0.95, 0.95), xycoords="axes fraction", ha="right", va="top",
        fontsize=8, bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7),
    )

    ax = axes[1]
    vals2 = df["n_occupied_neighbor_bins"].values
    ax.hist(vals2, bins=np.arange(0, vals2.max() + 2) - 0.5,
            color="#55A868", edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Occupied neighbor bins")
    ax.set_ylabel("Number of 1-cell bins")
    ax.set_title(f"{dataset_name.upper()}: Occupied neighbor bins")

    ax = axes[2]
    non_isolated = df[df["n_adjacent_cells"] > 0]
    if len(non_isolated):
        ax.hist(non_isolated["frac_same_type"].values, bins=20,
                color="#C44E52", edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Fraction same cell type in neighbors")
    ax.set_ylabel("Number of 1-cell bins")
    ax.set_title(f"{dataset_name.upper()}: Neighborhood type coherence")

    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(fig_dir, f"adjacent_bin_distribution_{dataset_name}.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close(fig)


def run_dataset(name, ref):
    adata = sc.read_h5ad(os.path.join(DATA, f"{name}/benchmark/{ref}/adata.h5ad"))
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    celltypes = adata.obs[LABEL].values
    bin_coords, bins_per_pc, _ = reproduce_binning_with_coordinates(X, FI)
    stats_df = compute_adjacent_stats(bin_coords, bins_per_pc, celltypes)
    if stats_df.empty:
        return stats_df
    stats_df["dataset"] = name
    plot_distribution(stats_df, name, FIG)
    return stats_df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", default="all", choices=["lcmv", "mcc", "all"])
    args = p.parse_args()
    os.makedirs(RES, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)

    datasets = REFS if args.dataset == "all" else {args.dataset: REFS[args.dataset]}
    all_stats = [run_dataset(name, ref) for name, ref in datasets.items()]
    all_stats = [s for s in all_stats if not s.empty]
    if not all_stats:
        return

    combined = pd.concat(all_stats, ignore_index=True)
    combined.to_csv(os.path.join(RES, "adjacent_bin_stats.csv"), index=False)

    summary_rows = []
    for dataset_name in combined["dataset"].unique():
        sub = combined[combined["dataset"] == dataset_name]
        non_iso = sub[~sub["is_isolated"]]
        summary_rows.append({
            "dataset": dataset_name,
            "n_one_cell_bins": len(sub),
            "n_isolated": sub["is_isolated"].sum(),
            "frac_isolated": sub["is_isolated"].mean(),
            "median_adjacent_cells": sub["n_adjacent_cells"].median(),
            "mean_adjacent_cells": sub["n_adjacent_cells"].mean(),
            "p25_adjacent_cells": sub["n_adjacent_cells"].quantile(0.25),
            "p75_adjacent_cells": sub["n_adjacent_cells"].quantile(0.75),
            "median_occupied_neighbors": sub["n_occupied_neighbor_bins"].median(),
            "median_frac_same_type": non_iso["frac_same_type"].median() if len(non_iso) else 0,
            "mean_frac_same_type": non_iso["frac_same_type"].mean() if len(non_iso) else 0,
        })
    pd.DataFrame(summary_rows).to_csv(os.path.join(RES, "adjacent_bin_summary.csv"), index=False)


if __name__ == "__main__":
    main()
