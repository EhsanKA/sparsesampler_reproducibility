import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")


def leiden_resolution(dataset="lcmv", sizes=(50000, 100000, 200000)):
    df = pd.read_csv(os.path.join(RES, "discovery_rate_summary.csv"))
    sub = df[df["dataset"] == dataset]
    resolutions = sorted(sub["resolution"].unique())
    n_types = 21 if dataset == "lcmv" else 5
    colors = {"sps": "#c44e52", "random": "#4c72b0"}

    fig, axes = plt.subplots(1, len(sizes), figsize=(4 * len(sizes), 3.5), sharey=True)
    if len(sizes) == 1:
        axes = [axes]

    for ax, size in zip(axes, sizes):
        chunk = sub[sub["size"] == size]
        for method in ("sps", "random"):
            m = chunk[chunk["method"] == method].sort_values("resolution")
            ax.plot(m["resolution"], m["mean_discovery"], "o-", label=method.upper(),
                    color=colors[method], linewidth=2, markersize=8)
        ax.axhline(n_types, color="gray", linestyle="--", linewidth=1, alpha=0.6)
        ax.set_title(f"{int(size):,} cells")
        ax.set_xlabel("Leiden resolution")
        ax.set_xticks(resolutions)
        ax.set_ylim(0, n_types + 1)
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel(f"Discoverable types (of {n_types})")
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(f"{dataset.upper()}: discovery vs Leiden resolution", y=1.08)
    fig.tight_layout()

    out = os.path.join(FIG, f"leiden_resolution_{dataset}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="lcmv", choices=["lcmv", "mcc", "both"])
    args = parser.parse_args()
    os.makedirs(FIG, exist_ok=True)
    datasets = ["lcmv", "mcc"] if args.dataset == "both" else [args.dataset]
    for ds in datasets:
        leiden_resolution(ds)
