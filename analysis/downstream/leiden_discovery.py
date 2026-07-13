import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sparsesampler.sampling import sample as sps_sample

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
LABEL = "celltype"
FI = 12
MIN_CELLS, MIN_PURITY = 10, 0.30
RESOLUTIONS = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
REFS = {"lcmv": 34, "mcc": 30}
THRESHOLDS = [10, 20, 50, 100, 200, 500]


def sample_indices(X, n_cells, method, size, seed):
    if method == "sps":
        idx, _ = sps_sample(X=X, size=size, seed=seed, feature_index=FI)
        return np.asarray(idx)
    return np.random.RandomState(seed).choice(n_cells, size=size, replace=False)


def evaluate(true_labels, cluster_labels, resolution):
    ari = adjusted_rand_score(true_labels, cluster_labels)
    nmi = normalized_mutual_info_score(true_labels, cluster_labels)
    types = np.unique(true_labels)
    discoverable, missing = [], []
    for ct in types:
        found = False
        for cl in np.unique(cluster_labels):
            mask = cluster_labels == cl
            n_ct = np.sum((true_labels == ct) & mask)
            if n_ct >= MIN_CELLS and n_ct / mask.sum() >= MIN_PURITY:
                found = True
                break
        discoverable.append(found)
        if not found:
            missing.append(ct)
    n_disc = sum(discoverable)
    return {
        "resolution": resolution,
        "ari": ari,
        "nmi": nmi,
        "n_clusters": len(np.unique(cluster_labels)),
        "n_cell_types": len(types),
        "n_discoverable": n_disc,
        "discovery_rate": n_disc / len(types),
        "discoverable_types": ";".join(t for t, d in zip(types, discoverable) if d),
        "missing_types": ";".join(missing),
    }


def sufficiency(labels):
    counts = pd.Series(labels).value_counts()
    rows = []
    for thresh in THRESHOLDS:
        ok = counts[counts >= thresh]
        rows.append({
            "min_cells_threshold": thresh,
            "n_types_sufficient": len(ok),
            "sufficient_types": ";".join(ok.index),
        })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True, choices=["lcmv", "mcc"])
    p.add_argument("--method", required=True, choices=["sps", "random"])
    p.add_argument("--size", type=int, required=True)
    p.add_argument("--rep", type=int, required=True)
    p.add_argument("--resolutions", type=float, nargs="+", default=None)
    args = p.parse_args()

    os.makedirs(RES, exist_ok=True)
    seed = 1234 + args.rep * 1000
    path = os.path.join(DATA, f"{args.dataset}/benchmark/{REFS[args.dataset]}/adata.h5ad")
    adata = sc.read_h5ad(path)
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    labels = adata.obs[LABEL].values
    if args.size >= X.shape[0]:
        return

    idx = sample_indices(X, X.shape[0], args.method, args.size, seed)
    sub = sc.AnnData(X[idx])
    sub.obs[LABEL] = labels[idx]
    sc.pp.pca(sub, n_comps=min(50, sub.shape[1] - 1, sub.shape[0] - 1))
    sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_pca")

    resolutions = args.resolutions or RESOLUTIONS
    results = []
    for res in resolutions:
        sc.tl.leiden(sub, resolution=res, key_added=f"leiden_{res}")
        r = evaluate(labels[idx], sub.obs[f"leiden_{res}"].values, res)
        r.update(dataset=args.dataset, method=args.method, size=args.size, rep=args.rep, seed=seed)
        results.append(r)

    tag = f"{args.dataset}_{args.method}_{args.size}_rep{args.rep}"
    pd.DataFrame(results).to_csv(os.path.join(RES, f"discovery_{tag}.csv"), index=False)
    suff = [{**s, "dataset": args.dataset, "method": args.method, "size": args.size, "rep": args.rep}
            for s in sufficiency(labels[idx])]
    pd.DataFrame(suff).to_csv(os.path.join(RES, f"sufficiency_{tag}.csv"), index=False)


if __name__ == "__main__":
    main()
