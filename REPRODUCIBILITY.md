# Reproducibility Guide

This document describes how to reproduce all results and figures from the SparseSampler benchmarking study.

## Prerequisites

1. **Cluster Configuration**: Copy and configure the cluster settings:
   ```bash
   cp config/cluster_config.sh.template config/cluster_config.sh
   # Edit config/cluster_config.sh with your cluster settings
   ```

2. **Conda Environment**: Activate the appropriate conda environment:
   ```bash
   conda activate facs_sampling  # or 'sps' depending on your setup
   ```

3. **Data**: Ensure the reference datasets are available in `data/` directory:
   - `data/mcc/` - Mouse Cell Compendium (1% rarity)
   - `data/mcc_01/` - MCC with 0.1% rarity
   - `data/mcc_05/` - MCC with 0.5% rarity
   - `data/lcmv/` - LCMV dataset

---

## Pipeline Overview

```
Step 1: Benchmark Data Generation
         ↓
Step 2: Feature Index Testing
         ↓
Step 3: Classification Analysis
         ↓
Step 4: Figure Generation
```

---

## Step 1: Benchmark Data Generation (Production Jobs)

Generate benchmark results for all sampling methods across datasets.

### MCC Dataset (1% rarity)
```bash
cd jobs/mcc_production/
bash 01_parallel_jobs_random.sh      # Random sampling
bash 02_parallel_jobs_sps.sh         # SparseSampler (SPS)
bash 03_parallel_jobs_hopper.sh      # Hopper
bash 04_parallel_jobs_atom.sh        # ATOM
bash 06_parallel_jobs_scsampler.sh   # scSampler
```

### MCC_01 Dataset (0.1% rarity)
```bash
cd jobs/mcc_01_production/
bash 01_parallel_jobs_random.sh
bash 02_parallel_jobs_sps.sh
bash 03_parallel_jobs_hopper.sh
bash 04_parallel_jobs_atom.sh
bash 06_parallel_jobs_scsampler.sh
```

### MCC_05 Dataset (0.5% rarity)
```bash
cd jobs/mcc_05_production/
bash 01_parallel_jobs_random.sh
bash 02_parallel_jobs_sps.sh
bash 03_parallel_jobs_hopper.sh
bash 04_parallel_jobs_atom.sh
bash 06_parallel_jobs_scsampler.sh
```

### LCMV Dataset
```bash
cd jobs/lcmv_production/
bash 01_parallel_jobs_random.sh
bash 02_parallel_jobs_sps.sh
bash 03_parallel_jobs_hopper.sh
bash 04_parallel_jobs_atom.sh
bash 06_parallel_jobs_scsampler.sh
```

**Output**: Results are saved to `data/{dataset}/benchmark/{ref}/{method}/{size}/{rep}/results.pkl`

---

## Step 2: Feature Index Testing

Test different feature indices (EVR thresholds) for SPS.

### 2a. Run Feature Index Tests
```bash
cd jobs/test_feature_index/
bash run_test_feature_index.sh
```

**Output**: Results saved to `data/test_feature_index/{dataset}/{ref}/`

### 2b. Create Feature Index Tables
```bash
cd jobs/test_feature_index/
bash run_create_feature_index_table.sh
```

**Output**: Tables saved to `jobs/test_feature_index/tables/`:
- `mcc_feature_index_table_size_*.csv`
- `mcc_01_feature_index_table_size_*.csv`
- `mcc_05_feature_index_table_size_*.csv`
- `lcmv_feature_index_table_size_*.csv`

---

## Step 3: Classification Analysis (RF Classification)

Run Random Forest classification to compare SPS vs Random sampling.

### 3a. MCC Dataset
```bash
cd jobs/feature_index_classification/
bash run_all_feature_indices.sh

# After jobs complete, combine results:
python combine_results.py
```

**Output**: `jobs/feature_index_classification/results/all_feature_indices_summary.csv`

### 3b. MCC_01 and MCC_05 Datasets
```bash
cd jobs/feature_index_classification/

# For MCC_01:
bash run_mcc01_mcc05_slurm.sh mcc_01

# For MCC_05:
bash run_mcc01_mcc05_slurm.sh mcc_05

# After jobs complete, combine results:
python combine_results_unified.py --dataset mcc_01
python combine_results_unified.py --dataset mcc_05
```

**Output**: 
- `jobs/feature_index_classification/results_mcc_01/all_feature_indices_summary.csv`
- `jobs/feature_index_classification/results_mcc_05/all_feature_indices_summary.csv`

### 3c. LCMV Dataset
```bash
cd jobs/feature_index_classification_lcmv/
bash run_all_feature_indices_lcmv.sh

# After jobs complete:
python combine_results.py
```

**Output**: `jobs/feature_index_classification_lcmv/results/all_feature_indices_summary.csv`

---

## Step 4: Sampling Methods Table

Create summary tables for different sampling methods.

```bash
cd jobs/test_sampling_methods/
bash run_create_sampling_methods_table.sh
```

**Output**: Tables saved to `jobs/test_sampling_methods/tables/`

---

## Step 5: Figure Generation

Generate all figures using SLURM job array (recommended):

```bash
cd figures/
sbatch generate_figures_array.sh
```

This submits 5 parallel tasks:
- **Task 1**: Main coverage figures (MCC + LCMV)
- **Task 2**: Main time performance figures
- **Task 3**: Revision coverage figures (MCC_01 + MCC_05)
- **Task 4**: Revision time performance figures
- **Task 5**: UMAP figures

### Or run individual figure scripts manually:

#### Main Figures
```bash
# Coverage figures
python figures/main/coverage/coverage_combined.py

# Time performance figures
python figures/main/time_performance/combined_time_plot.py

# UMAP figures
python figures/main/umaps/combined_umaps.py
```

#### Revision Figures
```bash
# Coverage figures (MCC_01, MCC_05)
python figures/revision/coverage/coverage_combined.py

# Time performance figures (MCC_01, MCC_05)
python figures/revision/time_performance/combined_time_plot.py

# Combined EVR/RF/Rank figure
python figures/revision/combined_evr_rf_rank_figure.py
```

### Figure Outputs

| Figure | Script | Output |
|--------|--------|--------|
| Main coverage | `figures/main/coverage/coverage_combined.py` | `figures/figure1_combined_coverage.jpg` |
| Main time | `figures/main/time_performance/combined_time_plot.py` | `figures/figure1_combined_time.jpg` |
| Main UMAPs | `figures/main/umaps/combined_umaps.py` | `figures/main/umaps/combined_umaps.jpg` |
| Supp MCC coverage | `figures/main/coverage/coverage_combined.py` | `figures/supp_figure3_mcc_coverage.jpg` |
| Supp LCMV coverage | `figures/main/coverage/coverage_combined.py` | `figures/supp_figure4_lcmv_coverage.jpg` |
| Supp MCC time | `figures/main/time_performance/combined_time_plot.py` | `figures/supp_figure5-6_mcc_time.jpg` |
| Revision MCC_01 coverage | `figures/revision/coverage/coverage_combined.py` | `figures/revision/supp_figure3_mcc_01_coverage.jpg` |
| Revision MCC_05 coverage | `figures/revision/coverage/coverage_combined.py` | `figures/revision/supp_figure3_mcc_05_coverage.jpg` |
| EVR/RF/Rank figure | `figures/revision/combined_evr_rf_rank_figure.py` | `figures/revision/combined_evr_rf_rank_figure.jpg` |

---

## Optional: PCA Timing Benchmark

```bash
cd jobs/pca_timing/
bash run_pca_timing.sh
```

**Output**: `jobs/pca_timing/pca_timing_results.csv`

---

## Directory Structure Summary

```
sparsesampler_reproducibility/
├── config/                     # Cluster configuration
├── data/                       # Input data and benchmark results
│   ├── mcc/benchmark/          # MCC benchmark results
│   ├── mcc_01/benchmark/       # MCC_01 benchmark results
│   ├── mcc_05/benchmark/       # MCC_05 benchmark results
│   ├── lcmv/benchmark/         # LCMV benchmark results
│   └── test_feature_index/     # Feature index test results
├── jobs/                       # SLURM job scripts
│   ├── *_production/           # Benchmark data generation
│   ├── test_feature_index/     # Feature index testing
│   ├── feature_index_classification/  # RF classification (MCC)
│   ├── feature_index_classification_lcmv/  # RF classification (LCMV)
│   └── test_sampling_methods/  # Sampling methods comparison
├── figures/                    # Figure generation scripts and outputs
│   ├── main/                   # Main paper figures
│   └── revision/               # Revision figures
├── notebooks/                  # Data preprocessing notebooks
└── analysis/                   # Analysis utilities
```

---

## Troubleshooting

1. **Missing tables**: Run Step 2b before Step 5 (EVR figure requires tables)
2. **Missing classification results**: Run Step 3 before generating EVR figure
3. **Memory errors**: Increase `--mem` in SLURM scripts
4. **Module not found**: Ensure correct conda environment is activated
