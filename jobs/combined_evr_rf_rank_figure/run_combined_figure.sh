#!/bin/bash
#SBATCH --job-name=combined_evr_rf_rank
#SBATCH --output=combined_evr_rf_rank_%j.out
#SBATCH --error=combined_evr_rf_rank_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=medium

source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh
conda activate sps
export PROJECT_ROOT=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility
cd $PROJECT_ROOT
python figures/revision/combined_evr_rf_rank_figure.py
