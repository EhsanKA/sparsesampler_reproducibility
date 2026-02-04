# Feature Index Parameter Testing

This directory contains scripts to test the effects of the `feature_index` parameter (also known as `target_threshold_index` in sparsesampler) on the number of recovered rare cells.

## Overview

The test evaluates `feature_index` values from 10 to 30 (inclusive) across:
- **Datasets**: lcmv, mcc, mcc_01, mcc_05
- **References**: All available references for each dataset
  - lcmv: [1, 5, 10, 20, 34]
  - mcc, mcc_01, mcc_05: [5, 10, 20, 25, 30]
- **Sample sizes**: All sizes for each dataset
  - lcmv: [50000, 100000, 200000]
  - mcc, mcc_01, mcc_05: [50000, 100000, 200000, 300000]
- **Replications**: 1 (rep=0)
- **Method**: Only SPS (sparse sampling)

## Directory Structure

```
jobs/test_feature_index/
├── test_feature_index.py          # Main Python script for running tests
├── run_test_feature_index.sh      # SLURM script to submit jobs
├── plot_feature_index_results.py  # Script to generate plots
└── README.md                       # This file

data/test_feature_index/
└── {dataset}/{ref}/{feature_index}/{size}/{rep}/
    └── results.pkl                 # Saved results
```

## Usage

### 1. Running the Tests

Submit all jobs using the SLURM script:

```bash
cd jobs/test_feature_index
./run_test_feature_index.sh
```

This will submit separate job arrays for each dataset (lcmv, mcc, mcc_01, mcc_05), where each job processes all feature_index and size combinations for a single reference.

### 2. Running Individual Tests

You can also run individual tests manually:

```bash
# Process all feature_index and size combinations for a specific dataset/ref (efficient)
python test_feature_index.py --dataset mcc_01 --ref 20 --all

# Process a single combination
python test_feature_index.py --dataset mcc_01 --ref 20 --feature_index 15 --size 100000 --rep 0
```

### 3. Generating Plots

After the jobs complete, generate plots similar to supp_figure3_mcc:

```bash
cd jobs/test_feature_index
python plot_feature_index_results.py
```

The plots will be saved to `figures/revision/` with names like:
- `supp_figure3_lcmv_feature_index.jpg`
- `supp_figure3_mcc_feature_index.jpg`
- `supp_figure3_mcc_01_feature_index.jpg`
- `supp_figure3_mcc_05_feature_index.jpg`

## Features

1. **Efficient Processing**: The `--all` flag loads each dataset/reference once and generates all feature_index and size combinations, reducing I/O overhead.

2. **Incremental Saving**: Results are saved immediately after each combination is processed, so if a job fails, completed runs are preserved.

3. **Automatic Skipping**: The script checks if results already exist and skips them, making it safe to re-run.

4. **Error Handling**: If one combination fails, the script continues with the next, preventing complete job failure.

## Results Format

Results are saved as pickle files in the same format as the original benchmark scripts:
- `results.pkl`: Tuple containing `(indices, elapsed_time)`
  - `indices`: Array/list of sampled cell indices
  - `elapsed_time`: Time taken for sampling in seconds

## Notes

- The `feature_index` parameter corresponds to `target_threshold_index` in the sparsesampler `sample()` function.
- Results are stored in `data/test_feature_index/` directory.
- Logs are written to `jobs/test_feature_index/logs/`.
- The plotting script uses the same rare cell type identification logic as the revision figures.

