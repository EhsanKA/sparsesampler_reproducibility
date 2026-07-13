import os
import io
import re
import csv
import time
from contextlib import redirect_stdout

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from sparsesampler.sampling import sample

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)))
DATASETS = {
    "lcmv": "lcmv/benchmark/34/adata.h5ad",
    "mcc": "mcc/benchmark/30/adata.h5ad",
    "mcc_01": "mcc_01/benchmark/30/adata.h5ad",
    "mcc_05": "mcc_05/benchmark/30/adata.h5ad",
}


def measure_pca_time(X, size=50000, seed=1234, n_runs=3):
    pca_times, total_times = [], []
    for i in range(n_runs):
        buf = io.StringIO()
        with redirect_stdout(buf):
            t0 = time.time()
            sample(X=X, size=min(size, X.shape[0] - 1), seed=seed + i,
                   auto_k=False, k=None, feature_index=18)
            total_times.append(time.time() - t0)
        m = re.search(r"Elapsed time after PCA: ([\d.]+)", buf.getvalue())
        if m:
            pca_times.append(float(m.group(1)))
    return (np.mean(pca_times) if pca_times else None), np.mean(total_times)


def main():
    os.makedirs(RES, exist_ok=True)
    rows = []
    for name, rel in DATASETS.items():
        path = os.path.join(DATA, rel)
        if not os.path.exists(path):
            continue
        adata = sc.read_h5ad(path)
        X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
        pca_time, total_time = measure_pca_time(X)
        rows.append({
            "dataset": name,
            "n_cells": X.shape[0],
            "n_features": X.shape[1],
            "pca_time": pca_time,
            "total_time": total_time,
        })
        del adata, X

    with open(os.path.join(RES, "pca_timing_results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dataset", "n_cells", "n_features", "pca_time", "total_time"])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
