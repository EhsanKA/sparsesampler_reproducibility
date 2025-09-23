# Paper Figures

This directory contains the code and generated figures for the paper.

## Directory Structure

```
figures/
├── main/                    # Main figures
│   ├── umaps/              # UMAP visualizations
│   │   ├── mcc/           # MCC dataset UMAPs
│   │   │   └── 5_2_umap_samples_in_reference.ipynb  # UMAPs for different sampling strategies
│   │   └── lcmv/          # LCMV dataset UMAPs
│   │       └── 5_2_umap_samples_in_reference.ipynb  # UMAPs for different sampling strategies
│   ├── coverage/          # Cell coverage plots
│   │   └── 3_rare_cell_coverage_all.ipynb  # Coverage plots for both MCC and LCMV
│   └── time_performance/  # Time performance plots
│       └── 4_time_plot.ipynb  # Time performance plots for both MCC and LCMV
└── supplementary/         # Supplementary figures
    ├── umaps/            # Cell type UMAPs
    │   ├── mcc/         # MCC cell type UMAPs
    │   └── lcmv/        # LCMV cell type UMAPs
    ├── coverage/         # Reference size dependency coverage plots
    └── time_performance/ # Reference size dependency time plots
```

## Figure Descriptions

### Main Figures

#### 1.1 UMAPs – MCC Dataset
- Location: `main/umaps/mcc/5_2_umap_samples_in_reference.ipynb`
- Content: UMAP visualizations colored by different sampling strategies (SPS, Rare, Random)
- Dataset: MCC
- Notebook: Generates 3 UMAP plots showing the distribution of cells based on different sampling strategies

#### 1.2 UMAPs – LCMV Dataset
- Location: `main/umaps/lcmv/5_2_umap_samples_in_reference.ipynb`
- Content: UMAP visualizations colored by different sampling strategies (SPS, Rare, Random)
- Dataset: LCMV
- Notebook: Generates 3 UMAP plots showing the distribution of cells based on different sampling strategies

#### 1.3 Connected Scatter Plot – Cell Coverage
- Location: `main/coverage/3_rare_cell_coverage_all.ipynb`
- Content: Rare cell coverage plots for MCC (10k samples) and LCMV (10k samples)
- Metrics: Coverage vs Reference size
- Notebook: Generates coverage plots for both datasets, focusing on rare cell populations

#### 1.4 Connected Scatter Plot – Time Performance
- Location: `main/time_performance/4_time_plot.ipynb`
- Content: Computation time plots for MCC and LCMV (10k samples)
- Metrics: Time vs Reference size
- Notebook: Generates time performance plots for both datasets

### Supplementary Figures

#### 2.1 UMAP Cell Types
- Location: `supplementary/umaps/`
- Content: UMAP visualizations of annotated cell types for both MCC and LCMV datasets
- To be implemented: New notebooks for cell type visualization

#### 2.2 Rare Cell Coverage – Reference Size Dependency
- Location: `supplementary/coverage/`
- Content: Coverage plots showing dependency on reference size for both datasets
- To be implemented: Modified version of coverage notebook focusing on reference size dependency

#### 2.3 Time Performance – Reference Size Dependency
- Location: `supplementary/time_performance/`
- Content: Time performance plots showing dependency on reference size for both datasets
- To be implemented: Modified version of time performance notebook focusing on reference size dependency

## Usage

Each subdirectory contains:
- Python scripts for generating the figures
- Generated figure files
- Any necessary data files or intermediate results

### Running the Notebooks

1. Main Figures:
   - UMAPs: Run the respective notebooks in `main/umaps/{mcc,lcmv}/`
   - Coverage: Run `main/coverage/3_rare_cell_coverage_all.ipynb`
   - Time Performance: Run `main/time_performance/4_time_plot.ipynb`

2. Supplementary Figures:
   - To be implemented: New notebooks will be added for cell type UMAPs
   - To be implemented: Modified versions of coverage and time performance notebooks

### Dependencies
- Python environment with required packages (to be listed)
- Access to the MCC and LCMV datasets
- Required R packages for atomic sketching (to be listed)

### Data Requirements
- MCC dataset
- LCMV dataset
- Cell type annotations
- Pre-computed sampling results

Note: Additional setup instructions and dependencies will be added as the notebooks are finalized. 