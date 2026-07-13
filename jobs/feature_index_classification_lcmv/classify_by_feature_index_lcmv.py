import os
import argparse
import glob
import pickle
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score, precision_recall_fscore_support
from sklearn.preprocessing import LabelEncoder, StandardScaler

warnings.filterwarnings("ignore")

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
OUT = os.path.join(PROJECT_ROOT, "jobs/feature_index_classification_lcmv/results")
LABEL = "celltype"
TEST_SIZE = 0.2
N_SEEDS = 5
RARE_TYPES = ["interacting", "NK1_1_TCRgd_T", "cDC2", "pDCs", "CD4_LCMV_spec"]
BENCH = os.path.join(DATA, "lcmv/benchmark/34")
FI_PATH = os.path.join(DATA, "test_feature_index/lcmv/34")


def load_indices(path, sps=False):
    with open(path, "rb") as f:
        res = pickle.load(f)
    if sps:
        return np.asarray(res[0][0]), res[0][1], res[1]
    indices = res[0] if isinstance(res, tuple) else res
    return np.asarray(indices)


def prepare_features(adata):
    tmp = adata.copy()
    if "log1p" not in tmp.uns:
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
    X = tmp.X.toarray() if scipy.sparse.issparse(tmp.X) else tmp.X
    return X, tmp.var_names.tolist()


def classify(X, labels, train_idx, test_idx, le):
    scaler = StandardScaler()
    Xtr = scaler.fit_transform(X[train_idx])
    Xte = scaler.transform(X[test_idx])
    clf = RandomForestClassifier(
        n_estimators=100, max_depth=20, class_weight="balanced", n_jobs=-1, random_state=42,
    )
    clf.fit(Xtr, labels[train_idx])
    y_pred = clf.predict(Xte)
    y_test = labels[test_idx]
    prec, rec, f1, _ = precision_recall_fscore_support(y_test, y_pred, average=None)
    return {
        "accuracy": accuracy_score(y_test, y_pred),
        "macro_f1": f1_score(y_test, y_pred, average="macro"),
        "per_class": {le.classes_[i]: {"f1": f1[i]} for i in range(len(le.classes_))},
    }


def in_pool(sample_idx, pool):
    pool_set = set(pool)
    return np.array([i for i in sample_idx if i in pool_set])


def run_analysis(feature_index, size=100000, rep=0):
    os.makedirs(OUT, exist_ok=True)
    adata = sc.read_h5ad(os.path.join(BENCH, "adata.h5ad"))
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    sps_idx, evr, sps_time = load_indices(
        os.path.join(FI_PATH, f"{feature_index}/{size}/{rep}/results.pkl"), sps=True,
    )
    rand_idx = load_indices(os.path.join(BENCH, f"random/{size}/{rep}/results.pkl"))
    X, genes = prepare_features(adata)
    le = LabelEncoder()
    labels = le.fit_transform(adata.obs[LABEL].values)
    rare_present = [ct for ct in RARE_TYPES if ct in le.classes_]

    all_results = {"sps": [], "random": []}
    for seed in range(N_SEEDS):
        pool, test_idx = train_test_split(
            np.arange(len(labels)), test_size=TEST_SIZE, stratify=labels, random_state=seed,
        )
        for method, sample in [("sps", in_pool(sps_idx, pool)), ("random", in_pool(rand_idx, pool))]:
            r = classify(X, labels, sample, test_idx, le)
            r["seed"] = seed
            all_results[method].append(r)

    summary = {"feature_index": feature_index, "size": size, "evr": evr, "sps_time": sps_time}
    for method in ("sps", "random"):
        summary[f"{method}_accuracy_mean"] = np.mean([r["accuracy"] for r in all_results[method]])
        summary[f"{method}_accuracy_std"] = np.std([r["accuracy"] for r in all_results[method]])
        summary[f"{method}_macro_f1_mean"] = np.mean([r["macro_f1"] for r in all_results[method]])
        summary[f"{method}_macro_f1_std"] = np.std([r["macro_f1"] for r in all_results[method]])
        for ct in rare_present:
            summary[f"{method}_{ct}_f1_mean"] = np.mean([r["per_class"].get(ct, {}).get("f1", 0) for r in all_results[method]])
            summary[f"{method}_{ct}_f1_std"] = np.std([r["per_class"].get(ct, {}).get("f1", 0) for r in all_results[method]])
        for ct in le.classes_:
            safe = ct.replace(" ", "_").replace("/", "_")
            summary[f"{method}_{safe}_f1_mean"] = np.mean([r["per_class"][ct]["f1"] for r in all_results[method]])

    summary["accuracy_improvement"] = (summary["sps_accuracy_mean"] - summary["random_accuracy_mean"]) * 100
    summary["macro_f1_improvement"] = (summary["sps_macro_f1_mean"] - summary["random_macro_f1_mean"]) * 100
    for ct in rare_present:
        summary[f"{ct}_f1_improvement"] = (
            summary[f"sps_{ct}_f1_mean"] - summary[f"random_{ct}_f1_mean"]
        ) * 100

    with open(os.path.join(OUT, f"classification_feature_index_{feature_index}.pkl"), "wb") as f:
        pickle.dump({
            "summary": summary, "all_results": all_results, "gene_names": genes,
            "label_encoder_classes": le.classes_.tolist(), "rare_cell_types": rare_present,
        }, f)
    pd.DataFrame([summary]).to_csv(os.path.join(OUT, f"summary_feature_index_{feature_index}.csv"), index=False)
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--feature-index", type=int, required=True)
    p.add_argument("--size", type=int, default=100000)
    p.add_argument("--rep", type=int, default=0)
    args = p.parse_args()
    run_analysis(args.feature_index, args.size, args.rep)


if __name__ == "__main__":
    main()
