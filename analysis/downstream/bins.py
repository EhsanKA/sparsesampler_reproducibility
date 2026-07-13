import os
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.stats import ks_2samp
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
LABEL = "celltype"
FI = 12
REFS = {"lcmv": 34, "mcc": 30}
SIZES = (50000, 100000, 200000)


def reproduce_binning(X, feature_index=FI, fast=False):
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
    if fast:
        bin_arrays, mult = [], 1
        for i in range(num_pcs):
            col = X_pca[:, i]
            n_bins = filtered[i]
            edges = np.linspace(col.min(), col.max(), n_bins + 1)
            digitized = np.digitize(col, edges).astype(np.int64)
            bin_arrays.append(digitized * mult)
            mult *= (n_bins + 2)
        grid_id = bin_arrays[0]
        for arr in bin_arrays[1:]:
            grid_id = grid_id + arr
        _, inverse, counts = np.unique(grid_id, return_inverse=True, return_counts=True)
        return counts[inverse], counts

    df = pd.DataFrame(X_pca[:, :num_pcs], columns=[f"PC{i+1}" for i in range(num_pcs)])

    def digitize_col(data, n_bins):
        edges = np.linspace(data.min(), data.max(), n_bins + 1)
        return np.digitize(data, edges)

    bins = [digitize_col(df.iloc[:, i], filtered[i]) for i in range(num_pcs)]
    df["grid_cell"] = list(zip(*bins))
    return df


def quality_metrics(X):
    if scipy.sparse.issparse(X):
        total_counts = np.array(X.sum(axis=1)).flatten()
        n_features = np.array((X > 0).sum(axis=1)).flatten()
    else:
        total_counts = X.sum(axis=1)
        n_features = (X > 0).sum(axis=1)
    return {
        "total_counts": total_counts,
        "n_detected_features": n_features,
        "mean_expression": total_counts / max(X.shape[1], 1),
    }


def run_quality(name, ref):
    adata = sc.read_h5ad(os.path.join(DATA, f"{name}/benchmark/{ref}/adata.h5ad"))
    X_orig = adata.X
    X_dense = X_orig.toarray() if scipy.sparse.issparse(X_orig) else X_orig
    df = reproduce_binning(X_dense, feature_index=FI, fast=False)

    one_cell_bins = df["grid_cell"].value_counts()
    one_cell_bins = one_cell_bins[one_cell_bins == 1].index
    mask_1cell = df["grid_cell"].isin(one_cell_bins).values
    mask_rest = ~mask_1cell
    metrics = quality_metrics(X_orig)

    results, percentile_rows = [], []
    for metric_name, values in metrics.items():
        vals_1cell, vals_rest = values[mask_1cell], values[mask_rest]
        ks_stat, ks_pval = ks_2samp(vals_1cell, vals_rest)
        results.append({
            "dataset": name, "metric": metric_name,
            "group_1cell_mean": np.mean(vals_1cell), "group_1cell_median": np.median(vals_1cell),
            "group_1cell_std": np.std(vals_1cell),
            "group_rest_mean": np.mean(vals_rest), "group_rest_median": np.median(vals_rest),
            "group_rest_std": np.std(vals_rest),
            "ks_statistic": ks_stat, "ks_pvalue": ks_pval,
            "n_1cell": mask_1cell.sum(), "n_rest": mask_rest.sum(),
        })
        for pct in [5, 10, 25, 50, 75, 90, 95]:
            percentile_rows.append({
                "dataset": name, "metric": metric_name, "percentile": pct,
                "value_1cell": np.percentile(values[mask_1cell], pct),
                "value_rest": np.percentile(values[mask_rest], pct),
            })
    return pd.DataFrame(results), pd.DataFrame(percentile_rows)


def analyze_bin_sizes(cell_bin_size, celltypes, dataset_name):
    grouped = (
        pd.DataFrame({"bin_size": cell_bin_size, "celltype": celltypes})
        .groupby(["bin_size", "celltype"]).size().reset_index(name="n_cells")
    )
    grouped["dataset"] = dataset_name
    n_bins_map = {size: int((cell_bin_size == size).sum() // size) for size in grouped["bin_size"].unique()}
    grouped["n_bins"] = grouped["bin_size"].map(n_bins_map)
    return grouped


def analyze_one_cell_bins(cell_bin_size, celltypes, dataset_name):
    mask_1cell = cell_bin_size == 1
    total_one_cell, total_cells = mask_1cell.sum(), len(celltypes)
    ct_all, ct_1cell = pd.Series(celltypes), pd.Series(celltypes)[mask_1cell]
    full_counts, one_cell_counts = ct_all.value_counts(), ct_1cell.value_counts()
    rows = []
    for ct, n_total in full_counts.items():
        n_in_1cell = int(one_cell_counts.get(ct, 0))
        freq_in_dataset = n_total / total_cells
        freq_in_1cell_pool = n_in_1cell / total_one_cell if total_one_cell else 0
        rows.append({
            "dataset": dataset_name, "celltype": ct,
            "n_total_in_dataset": n_total, "freq_in_dataset": freq_in_dataset,
            "n_in_1cell_bins": n_in_1cell, "freq_in_1cell_pool": freq_in_1cell_pool,
            "frac_of_type_in_1cell_bins": n_in_1cell / n_total if n_total else 0,
            "enrichment_ratio": freq_in_1cell_pool / freq_in_dataset if freq_in_dataset else 0,
        })
    return pd.DataFrame(rows).sort_values("enrichment_ratio", ascending=False)


def simulate_sample_composition(cell_bin_size, dataset_name, sample_sizes=SIZES):
    unique_sizes, size_counts = np.unique(cell_bin_size, return_counts=True)
    order = np.argsort(unique_sizes)
    unique_sizes, size_counts = unique_sizes[order], size_counts[order]
    rows = []
    for target_size in sample_sizes:
        if target_size > len(cell_bin_size):
            continue
        cumulative = 0
        for bs, n_cells in zip(unique_sizes, size_counts):
            if cumulative >= target_size:
                break
            take = min(n_cells, target_size - cumulative)
            rows.append({
                "dataset": dataset_name, "sample_size": target_size,
                "bin_size_tier": int(bs), "n_cells_from_tier": int(take),
                "frac_of_sample": take / target_size,
            })
            cumulative += take
    return pd.DataFrame(rows)


def run_characterize(name, ref):
    adata = sc.read_h5ad(os.path.join(DATA, f"{name}/benchmark/{ref}/adata.h5ad"))
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    celltypes = adata.obs[LABEL].values
    cell_bin_size, _ = reproduce_binning(X, feature_index=FI, fast=True)
    return (
        analyze_bin_sizes(cell_bin_size, celltypes, name),
        analyze_one_cell_bins(cell_bin_size, celltypes, name),
        simulate_sample_composition(cell_bin_size, name),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--task", required=True, choices=["quality", "characterize"])
    p.add_argument("--dataset", default="all", choices=["lcmv", "mcc", "all"])
    args = p.parse_args()
    os.makedirs(RES, exist_ok=True)
    datasets = REFS if args.dataset == "all" else {args.dataset: REFS[args.dataset]}

    if args.task == "quality":
        results, pctiles = [], []
        for name, ref in datasets.items():
            r, p_df = run_quality(name, ref)
            results.append(r)
            pctiles.append(p_df)
        pd.concat(results, ignore_index=True).to_csv(
            os.path.join(RES, "quality_validation_summary.csv"), index=False)
        pd.concat(pctiles, ignore_index=True).to_csv(
            os.path.join(RES, "quality_validation_percentiles.csv"), index=False)
    else:
        bin_detail, one_cell, composition = [], [], []
        for name, ref in datasets.items():
            b, o, c = run_characterize(name, ref)
            bin_detail.append(b)
            one_cell.append(o)
            composition.append(c)
        pd.concat(bin_detail, ignore_index=True).to_csv(
            os.path.join(RES, "bin_size_celltype_detail.csv"), index=False)
        pd.concat(one_cell, ignore_index=True).to_csv(
            os.path.join(RES, "one_cell_bin_characterization.csv"), index=False)
        pd.concat(composition, ignore_index=True).to_csv(
            os.path.join(RES, "sample_composition_by_tier.csv"), index=False)


if __name__ == "__main__":
    main()
