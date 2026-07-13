import os
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans
from sparsesampler.sampling import sample as sps_sample

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
LABEL = "celltype"
FI = 12
REFS = {"lcmv": 34, "mcc": 30}
N_PCS = 50
SUBSTATE_CELLS_PER_CLUSTER = 500
SUBSTATE_K_MAX = 20
SUBSTATE_MIN_POP = 200
SUBSTATE_MIN_CELLS = 5
MAJOR_FRAC = 0.01
MAX_QUERY = 20000
MAX_REF = 100000
NN_WORKERS = -1
_CAP_SEED = 0


def sample_indices(X, n_cells, method, size, seed):
    if method == "sps":
        idx, _ = sps_sample(X=X, size=size, seed=seed, feature_index=FI)
        return np.asarray(idx)
    return np.random.RandomState(seed).choice(n_cells, size=size, replace=False)


def precompute(dataset):
    adata = sc.read_h5ad(os.path.join(DATA, f"{dataset}/benchmark/{REFS[dataset]}/adata.h5ad"))
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    X = np.asarray(X, dtype=np.float32)
    celltypes = adata.obs[LABEL].astype(str).values
    n_cells, n_genes = X.shape

    n_comps = min(N_PCS, n_genes - 1, n_cells - 1)
    X_pca = PCA(n_components=n_comps, random_state=0).fit_transform(X).astype(np.float32)
    del X

    substate = np.full(n_cells, -1, dtype=np.int32)
    for ct in np.unique(celltypes):
        mask = celltypes == ct
        n = int(mask.sum())
        if n < SUBSTATE_MIN_POP:
            substate[mask] = 0
            continue
        k = max(2, min(SUBSTATE_K_MAX, n // SUBSTATE_CELLS_PER_CLUSTER))
        km = MiniBatchKMeans(n_clusters=k, random_state=0, n_init=3, batch_size=1024)
        substate[mask] = km.fit_predict(X_pca[mask]).astype(np.int32)

    np.savez_compressed(
        os.path.join(RES, f"within_pop_ref_{dataset}.npz"),
        X_pca=X_pca, celltypes=celltypes, substate=substate,
    )


def empty_metrics():
    return {
        "hausdorff": np.nan, "hausdorff_norm": np.nan,
        "mean_nn": np.nan, "mean_nn_norm": np.nan,
        "p90_nn": np.nan, "p90_nn_norm": np.nan,
        "spread_ratio": 0.0, "centroid_shift_norm": np.nan,
        "n_substates_retained": 0, "substate_retention": 0.0,
    }


def evaluate(dataset, method, size, rep):
    ref_path = os.path.join(RES, f"within_pop_ref_{dataset}.npz")
    if not os.path.exists(ref_path):
        raise FileNotFoundError(f"Missing {ref_path}; run with --precompute first")

    seed = 1234 + rep * 1000
    ref = np.load(ref_path, allow_pickle=True)
    X_pca = ref["X_pca"]
    celltypes = ref["celltypes"].astype(str)
    substate = ref["substate"]
    n_cells = X_pca.shape[0]

    adata = sc.read_h5ad(os.path.join(DATA, f"{dataset}/benchmark/{REFS[dataset]}/adata.h5ad"))
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    X = np.asarray(X, dtype=np.float32)
    if X.shape[0] != n_cells:
        raise ValueError(f"Cell count mismatch: adata={X.shape[0]} vs ref={n_cells}")
    if size >= n_cells:
        return

    indices = sample_indices(X, n_cells, method, size, seed)
    del X
    sampled_mask = np.zeros(n_cells, dtype=bool)
    sampled_mask[indices] = True
    cap_rng = np.random.RandomState(_CAP_SEED)

    rows = []
    for ct in np.unique(celltypes):
        full_idx = np.where(celltypes == ct)[0]
        n_full = full_idx.size
        frac_full = n_full / n_cells
        A = X_pca[full_idx]
        mu_full = A.mean(axis=0)
        var_full = float(A.var(axis=0).sum())
        radius = float(np.sqrt(var_full)) if var_full > 0 else np.nan

        samp_idx = full_idx[sampled_mask[full_idx]]
        n_sampled = samp_idx.size
        sub_full = substate[full_idx]
        n_substates = int(np.unique(sub_full).size)

        row = {
            "dataset": dataset, "method": method, "size": size,
            "rep": rep, "seed": seed, "celltype": ct,
            "n_full": n_full, "frac_full": frac_full, "is_major": frac_full >= MAJOR_FRAC,
            "n_sampled": n_sampled, "n_substates": n_substates,
        }

        if n_sampled == 0:
            row.update(empty_metrics())
            rows.append(row)
            continue

        if n_sampled > MAX_REF:
            B = X_pca[samp_idx[cap_rng.choice(n_sampled, size=MAX_REF, replace=False)]]
        else:
            B = X_pca[samp_idx]
        A_q = A if n_full <= MAX_QUERY else A[cap_rng.choice(n_full, size=MAX_QUERY, replace=False)]

        nn_dist, _ = cKDTree(B).query(A_q, k=1, workers=NN_WORKERS)
        hausdorff = float(nn_dist.max())
        mean_nn = float(nn_dist.mean())
        p90_nn = float(np.percentile(nn_dist, 90))
        spread_ratio = float(B.var(axis=0).sum()) / var_full if var_full > 0 else np.nan
        centroid_shift = float(np.linalg.norm(B.mean(axis=0) - mu_full))
        centroid_shift_norm = centroid_shift / radius if radius and radius > 0 else np.nan

        sub_samp = substate[samp_idx]
        retained = sum(
            1 for s in np.unique(sub_full)
            if int(np.sum(sub_samp == s)) >= SUBSTATE_MIN_CELLS
        )
        substate_retention = retained / n_substates if n_substates > 0 else np.nan
        norm = lambda v: v / radius if radius and radius > 0 else np.nan
        row.update({
            "hausdorff": hausdorff, "hausdorff_norm": norm(hausdorff),
            "mean_nn": mean_nn, "mean_nn_norm": norm(mean_nn),
            "p90_nn": p90_nn, "p90_nn_norm": norm(p90_nn),
            "spread_ratio": spread_ratio, "centroid_shift_norm": centroid_shift_norm,
            "n_substates_retained": retained, "substate_retention": substate_retention,
        })
        rows.append(row)

    tag = f"within_pop_{dataset}_{method}_{size}_rep{rep}"
    pd.DataFrame(rows).to_csv(os.path.join(RES, f"{tag}.csv"), index=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--precompute", action="store_true")
    p.add_argument("--dataset", required=True, choices=["lcmv", "mcc"])
    p.add_argument("--method", choices=["sps", "random"])
    p.add_argument("--size", type=int)
    p.add_argument("--rep", type=int)
    args = p.parse_args()
    os.makedirs(RES, exist_ok=True)

    if args.precompute:
        precompute(args.dataset)
        return

    if args.method is None or args.size is None or args.rep is None:
        p.error("--method, --size, and --rep required unless --precompute")
    evaluate(args.dataset, args.method, args.size, args.rep)


if __name__ == "__main__":
    main()
