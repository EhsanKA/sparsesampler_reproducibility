# Sparse Sampler Reproducibility

Reproducibility workflow for sparse sampling experiments on single-cell datasets. Benchmarks random, SPS, Hopper, Atomic, and scSampler methods on LCMV and MCC datasets.

## Prerequisites

- SGE/Grid Engine cluster access
- Python 3.7+ with scanpy, pandas, numpy
- R environment
- Conda/Mamba

## Usage

### 1. Configuration

Set up cluster configuration:
```bash
cp config/cluster_config.sh.template config/cluster_config.sh
# Edit cluster_config.sh with your email and paths
```

### 2. Data Download

```bash
cd data && ./download.sh
```

### 3. Data Processing

Process LCMV data:
```bash
cd notebooks/lcmv
python 0_check_the_data.py
python 01_atomic_data_converter.py
```

Process MCC data:
```bash
cd notebooks/mcc
python 0_check_the_data.py
python 01_atomic_data_converter.py
```

### 4. Run Jobs

Submit LCMV jobs:
```bash
cd jobs/lcmv_production
./01_parallel_jobs_random.sh
./02_parallel_jobs_sps.sh
./03_parallel_jobs_hopper.sh
./04_parallel_jobs_atom.sh
./06_parallel_jobs_scsampler.sh
```

Submit MCC jobs:
```bash
cd jobs/mcc_production
./01_parallel_jobs_random.sh
./02_parallel_jobs_sps.sh
./03_parallel_jobs_hopper.sh
./04_parallel_jobs_atom.sh
./06_parallel_jobs_scsampler.sh
```

### 5. Generate Figures

```bash
cd figures/main
python coverage/coverage_combined.py
python time_performance/combined_time_plot.py
python umaps/combined_umaps.py
python umaps/lcmv/umaps_lcmv.py
python umaps/mcc/umaps_mcc.py
```

## Notes

- The directories `notebooks/lcmv/hopper` and `notebooks/mcc/hopper` are clones of the Hopper repository
- Configuration must be set before running jobs
- Jobs should be run sequentially as listed
- Check `logs/` directories for job status
