import os
import argparse
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
LABEL = "celltype"
TEST_SIZE = 0.2
N_SEEDS = 5
N_TOP_GENES = 2000

CONFIGS = {
    "mcc": {
        "bench": os.path.join(DATA, "mcc/benchmark/30"),
        "fi_path": os.path.join(DATA, "test_feature_index/mcc/30"),
        "out": os.path.join(PROJECT_ROOT, "jobs/feature_index_classification/results"),
        "rare": "osteoblast",
    },
    "mcc_01": {
        "bench": os.path.join(DATA, "mcc_01/benchmark/30"),
        "fi_path": os.path.join(DATA, "test_feature_index/mcc_01/30"),
        "out": os.path.join(PROJECT_ROOT, "jobs/feature_index_classification/results_mcc_01"),
        "rare": "osteoblast",
    },
    "mcc_05": {
        "bench": os.path.join(DATA, "mcc_05/benchmark/30"),
        "fi_path": os.path.join(DATA, "test_feature_index/mcc_05/30"),
        "out": os.path.join(PROJECT_ROOT, "jobs/feature_index_classification/results_mcc_05"),
        "rare": "osteoblast",
    },
}


def load_indices(path, sps=False):
    with open(path, "rb") as f:
        res = pickle.load(f)
    if sps:
        return np.asarray(res[0][0]), res[0][1], res[1]
    return np.asarray(res[0] if isinstance(res, tuple) else res)


def prepare_features(adata, n_top=N_TOP_GENES):
    tmp = adata.copy()
    if "log1p" not in tmp.uns:
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(tmp, n_top_genes=min(n_top, tmp.n_vars), subset=False)
    mask = tmp.var["highly_variable"].values
    X = tmp[:, mask].X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    return X, tmp.var_names[mask].tolist()


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
    per_class = {
        le.classes_[i]: {"precision": prec[i], "recall": rec[i], "f1": f1[i]}
        for i in range(len(le.classes_))
    }
    return {
        "accuracy": accuracy_score(y_test, y_pred),
        "macro_f1": f1_score(y_test, y_pred, average="macro"),
        "per_class": per_class,
    }


def in_pool(sample_idx, pool):
    pool_set = set(pool)
    return np.array([i for i in sample_idx if i in pool_set])


def run_analysis(dataset, feature_index, size=100000, rep=0):
    cfg = CONFIGS[dataset]
    rare = cfg["rare"]
    os.makedirs(cfg["out"], exist_ok=True)
    out_csv = os.path.join(cfg["out"], f"summary_feature_index_{feature_index}.csv")
    if os.path.exists(out_csv):
        return None

    adata = sc.read_h5ad(os.path.join(cfg["bench"], "adata.h5ad"))
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    sps_idx, evr, sps_time = load_indices(
        os.path.join(cfg["fi_path"], f"{feature_index}/{size}/{rep}/results.pkl"), sps=True,
    )
    rand_idx = load_indices(os.path.join(cfg["bench"], f"random/{size}/{rep}/results.pkl"))
    X, genes = prepare_features(adata)
    le = LabelEncoder()
    labels = le.fit_transform(adata.obs[LABEL].values)

    all_results = {"sps": [], "random": []}
    for seed in range(N_SEEDS):
        pool, test_idx = train_test_split(
            np.arange(len(labels)), test_size=TEST_SIZE, stratify=labels, random_state=seed,
        )
        for method, sample in [("sps", in_pool(sps_idx, pool)), ("random", in_pool(rand_idx, pool))]:
            r = classify(X, labels, sample, test_idx, le)
            r["seed"] = seed
            all_results[method].append(r)

    summary = {"dataset": dataset, "feature_index": feature_index, "size": size, "evr": evr, "sps_time": sps_time}
    for method in ("sps", "random"):
        summary[f"{method}_accuracy_mean"] = np.mean([r["accuracy"] for r in all_results[method]])
        summary[f"{method}_accuracy_std"] = np.std([r["accuracy"] for r in all_results[method]])
        summary[f"{method}_macro_f1_mean"] = np.mean([r["macro_f1"] for r in all_results[method]])
        summary[f"{method}_macro_f1_std"] = np.std([r["macro_f1"] for r in all_results[method]])
        summary[f"{method}_{rare}_f1_mean"] = np.mean([r["per_class"][rare]["f1"] for r in all_results[method]])
        summary[f"{method}_{rare}_f1_std"] = np.std([r["per_class"][rare]["f1"] for r in all_results[method]])
        for ct in le.classes_:
            summary[f"{method}_{ct}_f1_mean"] = np.mean([r["per_class"][ct]["f1"] for r in all_results[method]])

    summary["accuracy_improvement"] = (summary["sps_accuracy_mean"] - summary["random_accuracy_mean"]) * 100
    summary["macro_f1_improvement"] = (summary["sps_macro_f1_mean"] - summary["random_macro_f1_mean"]) * 100
    summary[f"{rare}_f1_improvement"] = (summary[f"sps_{rare}_f1_mean"] - summary[f"random_{rare}_f1_mean"]) * 100

    with open(os.path.join(cfg["out"], f"classification_feature_index_{feature_index}.pkl"), "wb") as f:
        pickle.dump({"summary": summary, "all_results": all_results, "gene_names": genes,
                     "label_encoder_classes": le.classes_.tolist()}, f)
    pd.DataFrame([summary]).to_csv(out_csv, index=False)
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True, choices=list(CONFIGS))
    p.add_argument("--feature-index", type=int, required=True)
    p.add_argument("--size", type=int, default=100000)
    p.add_argument("--rep", type=int, default=0)
    args = p.parse_args()
    run_analysis(args.dataset, args.feature_index, args.size, args.rep)


if __name__ == "__main__":
    main()
