# SparseSampler Reproducibility - Project Structure Diagram

This document shows the relationships between Python scripts, shell scripts, and their outputs (figures, tables, pickle files).

```
================================================================================
                    SPARSESAMPLER REPRODUCIBILITY PROJECT
                        Scripts --> Data --> Figures Pipeline
================================================================================

+-----------------------------------------------------------------------------+
|                              1. DATA GENERATION                              |
|                    (notebooks/ + jobs/*_production/*.sh)                     |
+-----------------------------------------------------------------------------+

  SLURM Job Scripts                    Python Sampling Scripts
  ----------------                     ----------------------
  jobs/lcmv_production/                notebooks/lcmv/parallel.py
    |-- 01_parallel_jobs_random.sh  -------------+
    |-- 02_parallel_jobs_sps.sh     -------------|
    |-- 03_parallel_jobs_hopper.sh  -------------+--> data/lcmv/benchmark/{ref}/{method}/{size}/{rep}/
    |-- 04_parallel_jobs_atom.sh    -------------|         +-- results.pkl  (indices, runtime)
    +-- 06_parallel_jobs_scsampler.sh -----------+

  jobs/mcc_production/                 notebooks/mcc/parallel.py
    |-- 01_parallel_jobs_random.sh  -------------+
    |-- 02_parallel_jobs_sps.sh     -------------|
    |-- 03_parallel_jobs_hopper.sh  -------------+--> data/mcc/benchmark/{ref}/{method}/{size}/{rep}/
    |-- 04_parallel_jobs_atom.sh    -------------|         +-- results.pkl
    +-- 06_parallel_jobs_scsampler.sh -----------+

  jobs/mcc_01_production/              notebooks/mcc_01/parallel.py
    +-- (same structure)           -----------------> data/mcc_01/benchmark/.../results.pkl

  jobs/mcc_05_production/              notebooks/mcc_05/parallel.py
    +-- (same structure)           -----------------> data/mcc_05/benchmark/.../results.pkl


+-----------------------------------------------------------------------------+
|                          2. MAIN FIGURE GENERATION                           |
|                              (figures/main/)                                 |
+-----------------------------------------------------------------------------+

  figures/generate_all_figures.sh -----> Runs coverage + time performance scripts
  figures/generate_figures_array.sh --> SLURM array job running:
                                         - Task 1: Main coverage (MCC + LCMV)
                                         - Task 2: Main time performance (MCC + LCMV)
                                         - Task 3: Revision coverage (MCC_01 + MCC_05)
                                         - Task 4: Revision time (MCC_01 + MCC_05)
                                         - Task 5: Main UMAPs (MCC + LCMV)


  A. COVERAGE FIGURES
  -------------------
  figures/main/coverage/
    +-- coverage_combined.py
            |
            |--- Reads: data/{lcmv,mcc}/benchmark/{ref}/adata.h5ad
            |--- Reads: data/{lcmv,mcc}/benchmark/{ref}/{method}/{size}/{rep}/results.pkl
            |--- Uses:  analysis/refined_rare_cell_type_definition.py
            |
            +--> Outputs:
                 |-- figures/figure1_combined_coverage.jpg      (Main figure)
                 |-- figures/supp_figure3_mcc_coverage.jpg      (Supplementary)
                 +-- figures/supp_figure4_lcmv_coverage.jpg     (Supplementary)


  B. TIME PERFORMANCE FIGURES
  ---------------------------
  figures/main/time_performance/
    +-- combined_time_plot.py
            |
            |--- Reads: data/{lcmv,mcc}/benchmark/{ref}/{method}/{size}/{rep}/results.pkl
            |--- Reads: data/{lcmv,mcc}/benchmark/{ref}/atomic/{size}/{rep}/runtimes.csv
            |
            +--> Outputs:
                 |-- figures/figure1_lcmv_time.jpg
                 |-- figures/figure1_mcc_time.jpg
                 |-- figures/figure1_combined_time.jpg          (Main figure)
                 |-- figures/supp_figure5-6_lcmv_time.jpg       (Supplementary)
                 +-- figures/supp_figure5-6_mcc_time.jpg        (Supplementary)


  C. UMAP FIGURES (Dataset Visualization)
  ---------------------------------------
  figures/main/umaps/
    |-- lcmv/umaps_lcmv.py  -------> figures/main/umaps/lcmv/*.jpg
    |-- mcc/umaps_mcc.py    -------> figures/main/umaps/mcc/*.jpg
    +-- combined_umaps.py   -------> figures/main/umaps/combined_umaps.jpg
            |
            |--- Reads: data/{lcmv,mcc}/benchmark/{ref}/adata.h5ad
            |--- Reads: data/{lcmv,mcc}/benchmark/{ref}/{method}/{size}/{rep}/results.pkl
            |
            +--> Run via: figures/generate_figures_array.sh (Task 5)


+-----------------------------------------------------------------------------+
|                         3. REVISION FIGURE GENERATION                        |
|                             (figures/revision/)                              |
+-----------------------------------------------------------------------------+

  A. COMBINED EVR/RF/RANK FIGURE
  ------------------------------
  jobs/combined_evr_rf_rank_figure/
    +-- run_combined_figure.sh ---> figures/revision/combined_evr_rf_rank_figure.py
                                          |
                                          |--- Reads: jobs/test_feature_index/tables/*_feature_index_table_*.csv
                                          |--- Reads: jobs/feature_index_classification/results/all_feature_indices_summary.csv
                                          |
                                          +--> Outputs:
                                               |-- figures/revision/combined_evr_rf_rank_figure.jpg
                                               |-- figures/revision/combined_evr_rf_rank_figure.pdf
                                               |-- figures/revision/combined_evr_rf_rank_figure_50k.jpg
                                               |-- figures/revision/combined_evr_rf_rank_figure_100k.jpg
                                               +-- figures/revision/combined_evr_rf_rank_figure_200k.jpg


  B. COVERAGE (Revision)
  ----------------------
  figures/revision/coverage/
    +-- coverage_combined.py   ---> figures/revision/figure1_combined_coverage.jpg
                                    figures/revision/supp_figure3_mcc_01_coverage.jpg
                                    figures/revision/supp_figure3_mcc_05_coverage.jpg


  C. TIME PERFORMANCE (Revision)
  ------------------------------
  figures/revision/time_performance/
    +-- combined_time_plot.py ---> figures/revision/figure1_mcc_01_time.jpg
                                   figures/revision/figure1_mcc_05_time.jpg
                                   figures/revision/supp_figure5-6_mcc_01_time.jpg
                                   figures/revision/supp_figure5-6_mcc_05_time.jpg


+-----------------------------------------------------------------------------+
|                            4. FEATURE INDEX TESTING                          |
|                           (jobs/test_feature_index/)                         |
+-----------------------------------------------------------------------------+

  SLURM Scripts                        Python Scripts                      Outputs
  -------------                        --------------                      -------
  run_test_feature_index.sh  --------> test_feature_index.py
                                             |
                                             +--> data/test_feature_index/{dataset}/{ref}/.../results.pkl

  run_create_feature_index_table.sh -> create_feature_index_table.py
                                             |
                                             +--> jobs/test_feature_index/tables/
                                                    |-- mcc_feature_index_table_size_*.csv
                                                    |-- mcc_01_feature_index_table_size_*.csv
                                                    |-- mcc_05_feature_index_table_size_*.csv
                                                    +-- lcmv_feature_index_table_size_*.csv

  run_plot_feature_index.sh ---------> plot_feature_index_results.py
                                             +--> figures/revision/supp_figure3_*_feature_index.jpg

  run_evr_stability.sh --------------> analyze_evr_stability.py
                                             +--> figures/revision/supp_evr_stability*.jpg

  run_evr_coverage_plot.sh ----------> plot_evr_from_tables_with_coverage.py
                                             +--> figures/revision/supp_evr_sensitivity_coverage*.jpg


+-----------------------------------------------------------------------------+
|                         5. CLASSIFICATION BY FEATURE INDEX                   |
|                        (jobs/feature_index_classification/)                  |
+-----------------------------------------------------------------------------+

  SLURM Scripts                        Python Scripts                      Outputs
  -------------                        --------------                      -------
  run_all_feature_indices.sh --------> classify_by_feature_index.py
                                             |
                                             +--> jobs/feature_index_classification/results/
                                                    |-- feature_index_{fi}_summary.csv
                                                    +-- feature_index_{fi}_per_class.csv

                                       combine_results.py
                                             |
                                             +--> jobs/feature_index_classification/results/
                                                    +-- all_feature_indices_summary.csv

  run_mcc01_mcc05_slurm.sh ----------> run_mcc01_mcc05_classification.py
                                             +--> results_mcc_01/, results_mcc_05/

  (LCMV version)
  jobs/feature_index_classification_lcmv/
    +-- run_all_feature_indices_lcmv.sh -> classify_by_feature_index_lcmv.py
                                                 +--> results/all_feature_indices_summary.csv


+-----------------------------------------------------------------------------+
|                          6. PCA TIMING BENCHMARK                             |
|                           (jobs/pca_timing/)                                 |
+-----------------------------------------------------------------------------+

  SLURM Scripts                        Python Scripts                      Outputs
  -------------                        --------------                      -------
  run_pca_timing.sh ----------------> pca_timing_benchmark.py
                                             |
                                             +--> jobs/pca_timing/pca_timing_results.csv


+-----------------------------------------------------------------------------+
|                       7. SAMPLING METHODS COMPARISON                         |
|                       (jobs/test_sampling_methods/)                          |
+-----------------------------------------------------------------------------+

  SLURM Scripts                        Python Scripts                      Outputs
  -------------                        --------------                      -------
  run_create_sampling_methods_table.sh -> create_sampling_methods_table.py
                                             |
                                             +--> jobs/test_sampling_methods/tables/
                                                    |-- {dataset}_sampling_methods_table_size_*.csv
                                                    +-- (Used by combined_evr_rf_rank_figure.py)


+-----------------------------------------------------------------------------+
|                          8. CORE ANALYSIS MODULE                             |
|                               (analysis/)                                    |
+-----------------------------------------------------------------------------+

  analysis/
    +-- refined_rare_cell_type_definition.py   # REQUIRED: Used by coverage scripts


================================================================================
                              DATA FLOW SUMMARY
================================================================================

    +----------------+
    | Raw Data       |  (data/{dataset}/benchmark/{ref}/adata.h5ad)
    | (AnnData)      |
    +-------+--------+
            |
            v
    +----------------+     notebooks/{dataset}/parallel.py
    | Sampling       | <-- (called by jobs/*_production/*.sh)
    | Methods        |
    | SPS/Random/    |
    | Hopper/Atomic/ |
    | scSampler      |
    +-------+--------+
            |
            v
    +----------------+     data/{dataset}/benchmark/{ref}/{method}/{size}/{rep}/
    | Results        |     +-- results.pkl (indices, runtime)
    | (Pickle)       |
    +-------+--------+
            |
    +-------+-------+----------------+----------------+
    |               |                |                |
    v               v                v                v
+-----------+  +-----------+  +-----------+  +---------------+
| Coverage  |  |   Time    |  |   UMAP    |  | EVR/RF/Rank   |
| Analysis  |  | Analysis  |  | Analysis  |  |   Analysis    |
|           |  |           |  |           |  |               |
|figures/   |  |figures/   |  |figures/   |  |figures/       |
|main/      |  |main/      |  |main/      |  |revision/      |
|coverage/  |  |time_      |  |umaps/     |  |combined_evr_  |
+-----+-----+  |performance|  +-----+-----+  +-------+-------+
      |        +-----+-----+        |                |
      |              |              |                |
      v              v              v                v
+-------------------------------------------------------------+
|                    OUTPUT FIGURES                           |
+-------------------------------------------------------------+
| Main Figures (figures/):                                    |
|   - figure1_combined_coverage.jpg                           |
|   - figure1_combined_time.jpg                               |
|   - figure1_lcmv_time.jpg                                   |
|   - figure1_mcc_time.jpg                                    |
|   - supp_figure3_mcc_coverage.jpg                           |
|   - supp_figure4_lcmv_coverage.jpg                          |
|   - supp_figure5-6_{mcc,lcmv}_time.jpg                      |
|   - main/umaps/combined_umaps.jpg                           |
|                                                             |
| Revision Figures (figures/revision/):                       |
|   - combined_evr_rf_rank_figure_{50k,100k,200k}.{jpg,pdf}   |
|   - supp_figure3_mcc_01_coverage.jpg                        |
|   - supp_figure3_mcc_05_coverage.jpg                        |
|   - supp_figure5-6_mcc_01_time.jpg                          |
|   - supp_figure5-6_mcc_05_time.jpg                          |
|                                                             |
| Other Outputs:                                              |
|   - jobs/pca_timing/pca_timing_results.csv                  |
+-------------------------------------------------------------+
```

## Key File Relationships

| Category | Shell Scripts | Python Scripts | Outputs |
|----------|--------------|----------------|---------|
| Data Gen | `jobs/*_production/*.sh` | `notebooks/*/parallel.py` | `data/*/benchmark/**/results.pkl` |
| Coverage | `figures/generate_figures_array.sh` | `figures/main/coverage/coverage_combined.py` | `figure1_combined_coverage.jpg`, `supp_figure3/4_*_coverage.jpg` |
| Time | `figures/generate_figures_array.sh` | `figures/main/time_performance/combined_time_plot.py` | `figure1_*_time.jpg`, `supp_figure5-6_*_time.jpg` |
| UMAPs | `figures/generate_figures_array.sh` | `figures/main/umaps/combined_umaps.py` | `figures/main/umaps/combined_umaps.jpg` |
| EVR Analysis | `jobs/test_feature_index/*.sh` | `jobs/test_feature_index/*.py` | `tables/*_feature_index_table_*.csv` |
| Classification | `jobs/feature_index_classification/*.sh` | `jobs/feature_index_classification/*.py` | `results/all_feature_indices_summary.csv` |
| Sampling Methods | `jobs/test_sampling_methods/*.sh` | `jobs/test_sampling_methods/*.py` | `tables/*_sampling_methods_table_*.csv` |
| PCA Timing | `jobs/pca_timing/run_pca_timing.sh` | `jobs/pca_timing/pca_timing_benchmark.py` | `pca_timing_results.csv` |
| Combined Figures | `jobs/combined_evr_rf_rank_figure/run_combined_figure.sh` | `figures/revision/combined_evr_rf_rank_figure.py` | `combined_evr_rf_rank_figure_*.{jpg,pdf}` |

## Datasets

| Dataset | References | Sample Sizes | Rare Cell Threshold |
|---------|------------|--------------|---------------------|
| LCMV | 1M, 5M, 10M, 20M, 34M | 50k, 100k, 200k | 1% frequency |
| MCC | 0.5M, 1M, 2M, 2.5M, 3.2M | 50k, 100k, 200k, 300k | 1% frequency |
| MCC_01 | 0.5M, 1M, 2M, 2.5M, 3.2M | 50k, 100k, 200k, 300k | 0.1% frequency |
| MCC_05 | 0.5M, 1M, 2M, 2.5M, 3.2M | 50k, 100k, 200k, 300k | 0.5% frequency |

## Sampling Methods

| Method | Description |
|--------|-------------|
| `random` | Uniform random sampling |
| `sps` | Sparse Sampler (this project) |
| `hopper` | TreeHopper algorithm |
| `atomic` | Atomic Sketching (R-based) |
| `scsampler` | scSampler algorithm |

## Cleanup Candidates (Not in Main Pipeline)

The following files/directories are NOT required to regenerate the main figures and can be removed:

### Jobs Directory
- `jobs/clustering_visualization/` - Experimental visualization
- `jobs/de_analysis/` - Differential expression analysis (unused)
- `jobs/scanpy_marker_genes/` - Marker gene analysis (unused)
- `jobs/scanpy_marker_genes_lcmv/` - Marker gene analysis LCMV (unused)
- `jobs/semitones_analysis/` - SEMITONES analysis (unused)
- `jobs/semitones_marker_genes/` - SEMITONES marker genes (unused, 32GB)

### Analysis Directory
- `analysis/downstream_experiments/` - Statistical power & classification experiments (24GB)
- `analysis/scanpy_marker_genes/` - Marker gene pipeline (unused)
- `analysis/semitones_marker_genes/` - SEMITONES marker genes (unused)
- All standalone `.py` scripts except `refined_rare_cell_type_definition.py`
- All `.md` documentation files in analysis/
- All `.txt` and `.png` output files in analysis/

### Figures/Revision - Extra Files
- `figures/revision/rf_classification_*.py` - Extra RF classification scripts
- `figures/revision/rf_classification_*.jpg/.pdf` - Extra RF outputs
- `figures/revision/coverage/sps_feature_index_coverage.py` - Extra coverage script
- `figures/revision/umaps/` - Revision UMAP scripts (mcc_01, mcc_05)
- `figures/revision/sps_fi25_*.jpg` - Extra coverage outputs
- `figures/revision/supp_evr_*.jpg` - Extra EVR outputs
