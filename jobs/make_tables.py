import os
import sys
import argparse
import pickle

import pandas as pd
import scanpy as sc

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
DATA = os.path.join(PROJECT_ROOT, "data")
LABEL = "celltype"
FI_RANGE = list(range(1, 31))
REPS = [0]
METHODS = ["random", "sps", "hopper", "atomic", "scsampler"]
REP = 0

sys.path.insert(0, os.path.join(PROJECT_ROOT, "analysis"))
from refined_rare_cell_type_definition import (
    calculate_cell_type_distances,
    identify_rare_cell_types_distance_and_frequency,
)

CONFIGS = {
    "lcmv": {"refs": [1, 5, 10, 20, 34], "sizes": [50000, 100000, 200000], "freq_pct": 1.0},
    "mcc": {"refs": [5, 10, 20, 25, 30], "sizes": [50000, 100000, 200000, 300000], "freq_pct": 1.0},
    "mcc_01": {"refs": [5, 10, 20, 25, 30], "sizes": [50000, 100000, 200000, 300000], "freq_pct": 0.1},
    "mcc_05": {"refs": [5, 10, 20, 25, 30], "sizes": [50000, 100000, 200000, 300000], "freq_pct": 0.5},
}

OUT = {
    "feature_index": os.path.join(PROJECT_ROOT, "jobs", "test_feature_index", "tables"),
    "sampling_methods": os.path.join(PROJECT_ROOT, "jobs", "test_sampling_methods", "tables"),
}
FI_RESULTS = os.path.join(DATA, "test_feature_index")


def load_obs(dataset, ref):
    obs_path = os.path.join(DATA, f"{dataset}/benchmark/{ref}/obs.csv")
    if os.path.exists(obs_path):
        df = pd.read_csv(obs_path, index_col=0)
        df[LABEL] = df[LABEL].astype("category")
        return df
    adata_path = os.path.join(DATA, f"{dataset}/benchmark/{ref}/adata.h5ad")
    return sc.read_h5ad(adata_path).obs if os.path.exists(adata_path) else None


def rare_types(dataset, ref, obs):
    adata = sc.read_h5ad(os.path.join(DATA, f"{dataset}/benchmark/{ref}/adata.h5ad"))
    df = calculate_cell_type_distances(adata, label_key=LABEL)
    rare, _, _ = identify_rare_cell_types_distance_and_frequency(
        df, obs.shape[0], distance_percentile=75,
        frequency_threshold_pct=CONFIGS[dataset]["freq_pct"],
    )
    return rare


def load_ref_data(dataset):
    ref_data = {}
    for ref in CONFIGS[dataset]["refs"]:
        obs = load_obs(dataset, ref)
        ref_data[ref] = {
            "obs": obs,
            "rare_types": rare_types(dataset, ref, obs) if obs is not None else [],
        }
    return ref_data


def feature_index_table(dataset, size, ref_data):
    refs = CONFIGS[dataset]["refs"]
    rows = []
    for fi in FI_RANGE:
        row = {"feature_index": fi}
        for ref in refs:
            col = f"ref_{ref}M"
            if ref_data[ref]["obs"] is None or not ref_data[ref]["rare_types"]:
                row[col] = 0
                continue
            obs = ref_data[ref]["obs"]
            rare = set(ref_data[ref]["rare_types"])
            total, count = 0, 0
            for rep in REPS:
                pkl = os.path.join(FI_RESULTS, dataset, str(ref), str(fi), str(size), str(rep), "results.pkl")
                if not os.path.exists(pkl):
                    continue
                with open(pkl, "rb") as f:
                    data = pickle.load(f)
                idx = data[0][0]
                vc = obs.iloc[idx][LABEL].value_counts()
                total += int(vc[vc.index.isin(rare)].sum())
                count += 1
            row[col] = total / count if count else 0
        rows.append(row)
    df = pd.DataFrame(rows).set_index("feature_index")
    df.columns = [f"{r}M" for r in refs]
    return df


def load_method_indices(dataset, ref, method, size):
    bench = os.path.join(DATA, f"{dataset}/benchmark", str(ref))
    if method == "atomic":
        path = os.path.join(bench, "atomic", str(size), str(REP), "results.csv")
        if not os.path.exists(path):
            return None
        df = pd.read_csv(path)
        if dataset == "mcc" or dataset.startswith("mcc_"):
            return df["x"].values.astype(str).tolist()
        return df["x"].values.astype(int).tolist()
    path = os.path.join(bench, method, str(size), str(REP), "results.pkl")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return pickle.load(f)[0]


def count_rare(obs, indices, rare_types, method):
    if indices is None or len(indices) == 0:
        return 0
    sampled = obs.loc[indices] if method == "atomic" else obs.iloc[indices]
    vc = sampled[LABEL].value_counts()
    return int(vc[vc.index.isin(set(rare_types))].sum())


def sampling_methods_table(dataset, size, ref_data):
    refs = CONFIGS[dataset]["refs"]
    rows = []
    for method in METHODS:
        row = {"method": method}
        for ref in refs:
            col = f"ref_{ref}M"
            if ref_data[ref]["obs"] is None or not ref_data[ref]["rare_types"]:
                row[col] = 0
                continue
            idx = load_method_indices(dataset, ref, method, size)
            row[col] = count_rare(ref_data[ref]["obs"], idx, ref_data[ref]["rare_types"], method)
        rows.append(row)
    df = pd.DataFrame(rows).set_index("method")
    df.columns = [f"{r}M" for r in refs]
    return df


def fi_combinations():
    order = ["lcmv", "mcc", "mcc_01", "mcc_05"]
    return [(d, s) for d in order for s in CONFIGS[d]["sizes"]]


def run_feature_index(dataset=None, size=None, array_index=None):
    os.makedirs(OUT["feature_index"], exist_ok=True)
    if array_index is not None:
        combos = fi_combinations()
        if array_index < 0 or array_index >= len(combos):
            sys.exit(1)
        dataset, size = combos[array_index]
    if dataset is None or size is None:
        raise ValueError("provide --array-index or both --dataset and --size")
    ref_data = load_ref_data(dataset)
    df = feature_index_table(dataset, size, ref_data)
    df.to_csv(os.path.join(OUT["feature_index"], f"{dataset}_feature_index_table_size_{size}.csv"))


def run_sampling_methods():
    os.makedirs(OUT["sampling_methods"], exist_ok=True)
    for dataset in ["mcc", "mcc_01", "mcc_05"]:
        ref_data = load_ref_data(dataset)
        for size in CONFIGS[dataset]["sizes"]:
            df = sampling_methods_table(dataset, size, ref_data)
            df.to_csv(os.path.join(OUT["sampling_methods"], f"{dataset}_sampling_methods_table_size_{size}.csv"))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--type", required=True, choices=["feature_index", "sampling_methods"])
    p.add_argument("--dataset", choices=list(CONFIGS))
    p.add_argument("--size", type=int)
    p.add_argument("--array-index", type=int)
    args = p.parse_args()

    if args.type == "feature_index":
        run_feature_index(args.dataset, args.size, args.array_index)
    else:
        run_sampling_methods()


if __name__ == "__main__":
    main()
