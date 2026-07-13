import os
import argparse
import subprocess

PROJECT_ROOT = os.environ.get(
    "PROJECT_ROOT",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)
SCRIPT = os.path.join(PROJECT_ROOT, "jobs/feature_index_classification/classify_by_feature_index_unified.py")


def run_classification(dataset, start, end):
    for fi in range(start, end + 1):
        subprocess.run(
            ["python", SCRIPT, "--dataset", dataset, "--feature-index", str(fi), "--size", "100000", "--rep", "0"],
            check=False,
        )


def combine(dataset):
    import subprocess
    subprocess.run(
        ["python", os.path.join(PROJECT_ROOT, "jobs/combine_results.py"), "--dataset", dataset],
        check=False,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True, choices=["mcc_01", "mcc_05", "all"])
    p.add_argument("--combine-only", action="store_true")
    p.add_argument("--start-index", type=int, default=1)
    p.add_argument("--end-index", type=int, default=30)
    args = p.parse_args()
    datasets = ["mcc_01", "mcc_05"] if args.dataset == "all" else [args.dataset]
    for d in datasets:
        if not args.combine_only:
            run_classification(d, args.start_index, args.end_index)
        combine(d)


if __name__ == "__main__":
    main()
