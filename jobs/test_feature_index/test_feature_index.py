import os
import argparse
import pickle
import time

import numpy as np
import scanpy as sc
import scipy.sparse
from sparsesampler.sampling import sample

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
OUT = os.path.join(DATA, "test_feature_index")
LABEL = "celltype"
FI_RANGE = list(range(1, 31))
REPS = [0]

CONFIGS = {
    "lcmv": {"refs": [1, 5, 10, 20, 34], "sizes": [50000, 100000, 200000]},
    "mcc": {"refs": [5, 10, 20, 25, 30], "sizes": [50000, 100000, 200000, 300000]},
    "mcc_01": {"refs": [5, 10, 20, 25, 30], "sizes": [50000, 100000, 200000, 300000]},
    "mcc_05": {"refs": [5, 10, 20, 25, 30], "sizes": [50000, 100000, 200000, 300000]},
}


def result_path(dataset, ref, fi, size, rep):
    return os.path.join(OUT, dataset, str(ref), str(fi), str(size), str(rep), "results.pkl")


def run_sps(X, size, feature_index, seed):
    k = size // 100
    t0 = time.time()
    result = sample(X=X, size=size, seed=seed, auto_k=False, k=k, feature_index=feature_index)
    return result, time.time() - t0


def save_result(path, result, elapsed):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump((result, elapsed), f, protocol=pickle.HIGHEST_PROTOCOL)


def process_all(dataset, ref, seed_base=10000):
    cfg = CONFIGS[dataset]
    adata_path = os.path.join(DATA, f"{dataset}/benchmark/{ref}/adata.h5ad")
    if not os.path.exists(adata_path):
        return

    adata = sc.read_h5ad(adata_path)
    adata.obs[LABEL] = adata.obs[LABEL].astype("category")
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X

    for fi in FI_RANGE:
        for size in cfg["sizes"]:
            for rep in REPS:
                out = result_path(dataset, ref, fi, size, rep)
                if os.path.exists(out):
                    continue
                seed = seed_base + ref * 1000 + fi * 10 + size // 1000 + rep
                result, elapsed = run_sps(X, size, fi, seed)
                save_result(out, result, elapsed)


def process_one(dataset, ref, fi, size, rep, seed):
    out = result_path(dataset, ref, fi, size, rep)
    if os.path.exists(out):
        return
    adata = sc.read_h5ad(os.path.join(DATA, f"{dataset}/benchmark/{ref}/adata.h5ad"))
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    result, elapsed = run_sps(X, size, fi, seed)
    save_result(out, result, elapsed)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True, choices=list(CONFIGS))
    p.add_argument("--ref", type=int, required=True)
    p.add_argument("--feature_index", type=int)
    p.add_argument("--size", type=int)
    p.add_argument("--rep", type=int, default=0)
    p.add_argument("--seed", type=int)
    p.add_argument("--all", action="store_true")
    args = p.parse_args()

    cfg = CONFIGS[args.dataset]
    if args.ref not in cfg["refs"]:
        return

    seed_base = args.seed if args.seed is not None else 10000 + hash(args.dataset) % 1000
    if args.all:
        process_all(args.dataset, args.ref, seed_base)
    elif args.feature_index is not None and args.size is not None:
        seed = seed_base + args.ref * 1000 + args.feature_index * 10 + args.size // 1000 + args.rep
        process_one(args.dataset, args.ref, args.feature_index, args.size, args.rep, seed)


if __name__ == "__main__":
    main()
