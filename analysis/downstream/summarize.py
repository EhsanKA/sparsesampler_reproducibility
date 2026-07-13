import glob
import os

import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results")
RARE_DIR = os.path.join(os.path.dirname(ROOT), "results")


def _read_parts(pattern, skip=()):
    files = sorted(
        f for f in glob.glob(os.path.join(RES, pattern))
        if os.path.basename(f) not in skip
    )
    if not files:
        return None
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)


def merge_discovery():
    skip = {"discovery_rate_results.csv", "discovery_rate_summary.csv"}
    df = _read_parts("discovery_*.csv", skip)
    if df is None:
        return
    df.to_csv(os.path.join(RES, "discovery_rate_results.csv"), index=False)
    summary = df.groupby(["dataset", "method", "size", "resolution"]).agg(
        mean_ari=("ari", "mean"),
        std_ari=("ari", "std"),
        mean_nmi=("nmi", "mean"),
        mean_discovery=("n_discoverable", "mean"),
        std_discovery=("n_discoverable", "std"),
        mean_rate=("discovery_rate", "mean"),
    ).reset_index()
    summary.to_csv(os.path.join(RES, "discovery_rate_summary.csv"), index=False)


def merge_sufficiency():
    df = _read_parts("sufficiency_*.csv")
    if df is None:
        return
    df.to_csv(os.path.join(RES, "sufficiency_results.csv"), index=False)
    summary = df.groupby(["dataset", "method", "size", "min_cells_threshold"]).agg(
        mean_types=("n_types_sufficient", "mean"),
        std_types=("n_types_sufficient", "std"),
    ).reset_index()
    summary.to_csv(os.path.join(RES, "sufficiency_summary.csv"), index=False)


def recall_by_tier(resolution=1.0):
    rare_by_ds = {}
    for ds in ("lcmv", "mcc"):
        path = os.path.join(RARE_DIR, f"{ds}_rare_cell_types_analysis.csv")
        rare_df = pd.read_csv(path)
        rare = set(rare_df.loc[rare_df["rare_distance_and_frequency"], "cell_type"])
        rare_by_ds[ds] = (rare, set(rare_df["cell_type"]) - rare)

    rows = []
    for path in sorted(glob.glob(os.path.join(RES, "discovery_*.csv"))):
        base = os.path.basename(path)
        if base in ("discovery_rate_results.csv", "discovery_rate_summary.csv"):
            continue
        parts = base.replace(".csv", "").split("_")
        if len(parts) < 5:
            continue
        dataset, method, size, rep = parts[1], parts[2], int(parts[3]), int(parts[4].replace("rep", ""))
        rare, majority = rare_by_ds[dataset]
        row = pd.read_csv(path)
        row = row[row["resolution"] == resolution]
        if row.empty:
            continue
        discovered = set(str(row.iloc[0]["discoverable_types"]).split(";")) if row.iloc[0]["discoverable_types"] else set()
        discovered.discard("")
        for ct in sorted(rare | majority):
            rows.append({
                "dataset": dataset,
                "method": method,
                "size": size,
                "rep": rep,
                "resolution": resolution,
                "celltype": ct,
                "tier": "rare" if ct in rare else "majority",
                "recall": int(ct in discovered),
            })

    detail = pd.DataFrame(rows)
    detail.to_csv(os.path.join(RES, "per_celltype_recall_detail.csv"), index=False)
    summary = detail.groupby(["dataset", "method", "size", "tier"], as_index=False).agg(
        mean_recall=("recall", "mean"),
        n_types=("recall", "count"),
    )
    summary.to_csv(os.path.join(RES, "per_celltype_recall_summary.csv"), index=False)


def merge_within_pop():
    df = _read_parts("within_pop_*_rep*.csv")
    if df is None:
        return
    df.to_csv(os.path.join(RES, "within_pop_full.csv"), index=False)
    metrics = [
        "hausdorff_norm", "mean_nn_norm", "p90_nn_norm",
        "spread_ratio", "centroid_shift_norm", "substate_retention",
    ]
    summary = df.groupby(["dataset", "method", "size", "celltype"]).agg(
        {m: "mean" for m in metrics} | {"n_full": "first", "frac_full": "first", "is_major": "first"}
    ).reset_index()
    summary.to_csv(os.path.join(RES, "within_pop_summary.csv"), index=False)

    major = summary[summary["is_major"]]
    pivot = major.pivot_table(
        index=["dataset", "size", "celltype", "n_full", "frac_full"],
        columns="method",
        values=metrics,
    )
    pivot.columns = [f"{metric}_{method}" for metric, method in pivot.columns]
    pivot.reset_index().to_csv(os.path.join(RES, "within_pop_major_compare.csv"), index=False)


if __name__ == "__main__":
    os.makedirs(RES, exist_ok=True)
    merge_discovery()
    merge_sufficiency()
    recall_by_tier()
    merge_within_pop()
