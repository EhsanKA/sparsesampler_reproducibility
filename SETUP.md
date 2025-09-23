# Sparse Sampler Reproducibility

This repository contains scripts and notebooks for reproducing sparse sampling experiments.

## Setup Instructions

### 1. Environment Setup

1. Clone this repository
2. Copy the cluster configuration template:
   ```bash
   cp config/cluster_config.sh.template config/cluster_config.sh
   ```
3. Edit `config/cluster_config.sh` with your cluster-specific settings:
   - Update `CLUSTER_EMAIL` with your email address
   - Adjust paths if needed for your system
   - Set your conda environment name

### 2. Cluster Configuration

The job scripts are designed to work with SGE/Grid Engine. You'll need to:

1. Ensure you have access to a cluster with SGE
2. Have conda/mamba installed with the required environment
3. Customize the resource requirements in job scripts if needed

### 3. Running Jobs

Each dataset (LCMV, MCC) has separate job directories:
- `jobs/lcmv_production/`
- `jobs/mcc_production/`

To submit jobs:
```bash
cd jobs/lcmv_production
./01_parallel_jobs_random.sh
```

### 4. Data Requirements

- Place your data files in the appropriate data directories
- Ensure data files follow the expected format (see notebooks for examples)



## Directory Structure

```
├── config/                 # Configuration templates
├── data/                   # Data files (gitignored)
├── figures/               # Generated figures
├── jobs/                  # Job submission scripts
├── notebooks/             # Analysis notebooks
└── README.md
```