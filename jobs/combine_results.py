import os
import glob
import argparse

import pandas as pd

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
JOBS = os.path.join(PROJECT_ROOT, "jobs")

DATASETS = {
    "mcc": os.path.join(JOBS, "feature_index_classification", "results"),
    "mcc_01": os.path.join(JOBS, "feature_index_classification", "results_mcc_01"),
    "mcc_05": os.path.join(JOBS, "feature_index_classification", "results_mcc_05"),
    "lcmv": os.path.join(JOBS, "feature_index_classification_lcmv", "results"),
}


def combine(dataset):
    res_dir = DATASETS[dataset]
    files = sorted(glob.glob(os.path.join(res_dir, "summary_feature_index_*.csv")))
    if not files:
        return None
    combined = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    combined = combined.sort_values("feature_index")
    combined.to_csv(os.path.join(res_dir, "all_feature_indices_summary.csv"), index=False)
    return combined


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", nargs="+", choices=list(DATASETS) + ["all"], default=["all"])
    args = p.parse_args()
    datasets = list(DATASETS) if "all" in args.dataset else args.dataset
    for d in datasets:
        combine(d)


if __name__ == "__main__":
    main()
