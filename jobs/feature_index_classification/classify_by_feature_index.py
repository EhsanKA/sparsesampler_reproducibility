import argparse

from classify_by_feature_index_unified import run_analysis


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--feature-index", type=int, required=True)
    p.add_argument("--size", type=int, default=100000)
    p.add_argument("--rep", type=int, default=0)
    args = p.parse_args()
    run_analysis("mcc", args.feature_index, args.size, args.rep)


if __name__ == "__main__":
    main()
