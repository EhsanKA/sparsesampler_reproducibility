# SparseSampler Reproducibility

Benchmarking and reproducibility code for the SparseSampler (SPS) method on single-cell RNA-seq and flow cytometry data. Compares SPS against Random, Hopper, Atomic Sketch, and scSampler across MCC and LCMV datasets.

## Prerequisites

- SLURM cluster access
- Conda/Mamba with environment `sps`
- Python 3.7+ (scanpy, pandas, numpy, scipy, matplotlib, scikit-learn, umap-learn, scsampler)
- R with Seurat v5 (for Atomic Sketch)

## Setup

```bash
# 1. Clone and configure
cp config/cluster_config.sh.template config/cluster_config.sh
# Edit with your cluster settings

# 2. Download data
cd data && bash download.sh

# 3. Clone Hopper into each notebook directory
cd notebooks/lcmv && git clone https://github.com/bendemeo/hopper.git
cd ../mcc && git clone https://github.com/bendemeo/hopper.git
cd ../mcc_01 && git clone https://github.com/bendemeo/hopper.git
cd ../mcc_05 && git clone https://github.com/bendemeo/hopper.git
```

## Reproducing Results

### Step 1: Data Preprocessing

```bash
# LCMV
cd notebooks/lcmv && python 0_check_the_data.py && python 01_atomic_data_converter.py

# MCC
cd notebooks/mcc && python 0_check_the_data.py && python 01_atomic_data_converter.py

# MCC_01
cd notebooks/mcc_01 && python 0_check_the_data.py && python 01_atomic_data_converter.py

# MCC_05
cd notebooks/mcc_05 && python 0_check_the_data.py && python 01_atomic_data_converter.py
```

### Step 2: Run Sampling Benchmarks

Submit jobs for each dataset (`lcmv`, `mcc`, `mcc_01`, `mcc_05`):

```bash
cd jobs/{dataset}_production/
bash 01_parallel_jobs_random.sh
bash 02_parallel_jobs_sps.sh
bash 03_parallel_jobs_hopper.sh
bash 04_parallel_jobs_atom.sh
bash 06_parallel_jobs_scsampler.sh
```

Output: `data/{dataset}/benchmark/{ref}/{method}/{size}/{rep}/results.pkl`

### Step 3: Feature Index Testing

```bash
cd jobs/test_feature_index/
bash run_test_feature_index.sh           # Test feature indices 1-30
bash run_create_feature_index_table.sh   # Create summary tables
```

### Step 4: RF Classification Analysis

```bash
# MCC
cd jobs/feature_index_classification/
bash run_all_feature_indices.sh && python combine_results.py

# MCC_01 / MCC_05
bash run_mcc01_mcc05_slurm.sh mcc_01 && python combine_results_unified.py --dataset mcc_01
bash run_mcc01_mcc05_slurm.sh mcc_05 && python combine_results_unified.py --dataset mcc_05

# LCMV
cd jobs/feature_index_classification_lcmv/
bash run_all_feature_indices_lcmv.sh && python combine_results.py
```

### Step 5: Sampling Methods Tables

```bash
cd jobs/test_sampling_methods/
bash run_create_sampling_methods_table.sh
```

### Step 6: Generate Figures

Coverage, time, and UMAP figures via SLURM array:

```bash
sbatch figures/generate_figures_array.sh
```

EVR/RF/Rank analysis figure (run separately, requires Steps 3-5):

```bash
python figures/revision/combined_evr_rf_rank_figure.py
```

Or run all figure scripts individually:

```bash
# Main figures
python figures/main/coverage/coverage_combined.py
python figures/main/time_performance/combined_time_plot.py
python figures/main/umaps/combined_umaps.py

# Revision figures
python figures/revision/coverage/coverage_combined.py
python figures/revision/time_performance/combined_time_plot.py
python figures/revision/combined_evr_rf_rank_figure.py
```

### Optional: PCA Timing Benchmark

```bash
cd jobs/pca_timing/ && bash run_pca_timing.sh
```

## Project Structure

```
sparsesampler_reproducibility/
├── analysis/
│   └── refined_rare_cell_type_definition.py   # Rare cell type identification
├── config/                                     # Cluster configuration
├── data/
│   ├── lcmv/benchmark/                        # LCMV results (refs: 1M-34M)
│   ├── mcc/benchmark/                         # MCC results, 1% rarity (refs: 0.5M-3.2M)
│   ├── mcc_01/benchmark/                      # MCC results, 0.1% rarity
│   ├── mcc_05/benchmark/                      # MCC results, 0.5% rarity
│   └── test_feature_index/                    # Feature index sweep results
├── notebooks/
│   ├── lcmv/parallel.py                       # LCMV sampling driver
│   ├── mcc/parallel.py                        # MCC sampling driver
│   ├── mcc_01/parallel.py
│   └── mcc_05/parallel.py
├── jobs/
│   ├── {dataset}_production/                  # SLURM scripts per dataset
│   ├── test_feature_index/                    # Feature index sweep + tables
│   ├── test_sampling_methods/                 # Sampling methods comparison tables
│   ├── feature_index_classification/          # RF classification (MCC variants)
│   ├── feature_index_classification_lcmv/     # RF classification (LCMV)
│   ├── combined_evr_rf_rank_figure/           # EVR analysis figure job
│   └── pca_timing/                            # PCA timing benchmark
├── figures/
│   ├── generate_figures_array.sh              # SLURM array: all figures
│   ├── main/
│   │   ├── coverage/coverage_combined.py      # -> figure1_combined_coverage.jpg
│   │   ├── time_performance/combined_time_plot.py  # -> figure1_combined_time.jpg
│   │   └── umaps/combined_umaps.py            # -> combined_umaps.jpg
│   └── revision/
│       ├── coverage/coverage_combined.py      # -> revision coverage figures
│       ├── time_performance/combined_time_plot.py  # -> revision time figures
│       └── combined_evr_rf_rank_figure.py     # -> EVR/RF/Rank analysis figures
└── logs/                                      # Job logs
```

## Output Figures

| Figure | File | Description |
|--------|------|-------------|
| Main Coverage | `figure1_combined_coverage.jpg` | Rare cell coverage: SPS vs others (MCC + LCMV) |
| Main Time | `figure1_combined_time.jpg` | Runtime scalability comparison |
| UMAPs | `main/umaps/combined_umaps.jpg` | SPS vs Random spatial distribution |
| Supp MCC Coverage | `supp_figure3_mcc_coverage.jpg` | MCC coverage by sample/ref size |
| Supp LCMV Coverage | `supp_figure4_lcmv_coverage.jpg` | LCMV coverage by sample/ref size |
| Supp MCC Time | `supp_figure5-6_mcc_time.jpg` | MCC detailed runtime |
| Supp LCMV Time | `supp_figure5-6_lcmv_time.jpg` | LCMV detailed runtime |
| Rev. Coverage | `revision/figure1_combined_coverage.jpg` | Coverage at 0.1% and 0.5% rarity |
| Rev. MCC_01 Cov. | `revision/supp_figure3_mcc_01_coverage.jpg` | MCC 0.1% rarity coverage |
| Rev. MCC_05 Cov. | `revision/supp_figure3_mcc_05_coverage.jpg` | MCC 0.5% rarity coverage |
| Rev. Time | `revision/figure1_combined_time.jpg` | Runtime at 0.1% and 0.5% rarity |
| Rev. MCC_01 Time | `revision/supp_figure5-6_mcc_01_time.jpg` | MCC 0.1% rarity detailed runtime |
| Rev. MCC_05 Time | `revision/supp_figure5-6_mcc_05_time.jpg` | MCC 0.5% rarity detailed runtime |
| EVR Analysis | `revision/combined_evr_rf_rank_figure_{50k,100k,200k}.jpg` | EVR index sensitivity (coverage, RF F1, rank) |

## Key Parameters

- **EVR Feature Index**: 12 (default for all benchmarks)
- **Rare cell definition**: distance > 75th percentile AND frequency < threshold
  - MCC: < 1%, MCC_05: < 0.5%, MCC_01: < 0.1%, LCMV: < 1%

## Sampling Methods

| Method | Description |
|--------|-------------|
| `random` | Uniform random sampling (baseline) |
| `sps` | SparseSampler — this project |
| `hopper` | TreeHopper geometric sketching |
| `atomic` | Atomic Sketch (R/Seurat v5) |
| `scsampler` | scSampler algorithm |
