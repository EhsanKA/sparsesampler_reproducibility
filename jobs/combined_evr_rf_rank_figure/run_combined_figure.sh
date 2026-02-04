#!/bin/bash
#SBATCH --job-name=combined_evr_rf_rank
#SBATCH --output=combined_evr_rf_rank_%j.out
#SBATCH --error=combined_evr_rf_rank_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=medium

# Combined EVR, RF Classification, and Rank Figure Generation
# This script generates the 3x4 combined figure showing:
# - Row 1: Coverage % by EVR index
# - Row 2: RF Delta F1 by feature index
# - Row 3: SPS rank among methods by EVR index

# Activate conda environment
source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh
conda activate sps

# Set project root
export PROJECT_ROOT=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility

# Navigate to project root
cd $PROJECT_ROOT

# Run the figure generation script
echo "Starting combined EVR/RF/Rank figure generation..."
echo "Time: $(date)"
echo "========================================"

python figures/revision/combined_evr_rf_rank_figure.py

echo "========================================"
echo "Completed at: $(date)"
