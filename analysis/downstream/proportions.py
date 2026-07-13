import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.spatial.distance import jensenshannon
from sparsesampler.sampling import sample as sps_sample

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
LABEL = "celltype"
FI = 12
SIZES = [50000, 100000, 200000]
REFS = {"lcmv": 34, "mcc": 30}


def proportions(labels, all_types):
    counts = pd.Series(labels).value_counts()
    props = np.array([counts.get(ct, 0) for ct in all_types], dtype=float)
    return props / props.sum()


def run_dataset(name, ref):
    adata = sc.read_h5ad(os.path.join(DATA, f"{name}/benchmark/{ref}/adata.h5ad"))
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    labels = adata.obs[LABEL].values
    types = sorted(np.unique(labels))
    full = proportions(labels, types)

    jsd_rows, prop_rows = [], []
    for size in SIZES:
        if size >= X.shape[0]:
            continue
        for rep in range(3):
            seed = 1234 + rep * 1000
            for method in ("sps", "random"):
                if method == "sps":
                    idx, _ = sps_sample(X=X, size=size, seed=seed, feature_index=FI)
                    idx = np.asarray(idx)
                else:
                    idx = np.random.RandomState(seed).choice(X.shape[0], size=size, replace=False)
                sub = proportions(labels[idx], types)
                jsd_rows.append({"dataset": name, "method": method, "size": size, "rep": rep,
                                 "jsd": jensenshannon(full, sub)})
                for ct, fp, sp in zip(types, full, sub):
                    prop_rows.append({
                        "dataset": name, "method": method, "size": size, "rep": rep,
                        "celltype": ct, "full_proportion": fp, "sample_proportion": sp,
                        "absolute_diff": abs(sp - fp),
                        "relative_diff": (sp - fp) / fp if fp > 0 else 0,
                    })
    return pd.DataFrame(jsd_rows), pd.DataFrame(prop_rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", default="all", choices=["lcmv", "mcc", "all"])
    args = p.parse_args()

    os.makedirs(RES, exist_ok=True)
    datasets = {args.dataset: REFS[args.dataset]} if args.dataset != "all" else REFS

    jsd_parts, prop_parts = [], []
    for name, ref in datasets.items():
        jsd, prop = run_dataset(name, ref)
        jsd_parts.append(jsd)
        prop_parts.append(prop)

    jsd_df = pd.concat(jsd_parts, ignore_index=True)
    prop_df = pd.concat(prop_parts, ignore_index=True)
    jsd_df.to_csv(os.path.join(RES, "proportion_jsd.csv"), index=False)
    prop_df.to_csv(os.path.join(RES, "proportion_per_celltype.csv"), index=False)
    summary = jsd_df.groupby(["dataset", "method", "size"])["jsd"].agg(["mean", "std"]).reset_index()
    summary.to_csv(os.path.join(RES, "proportion_jsd_summary.csv"), index=False)


if __name__ == "__main__":
    main()
