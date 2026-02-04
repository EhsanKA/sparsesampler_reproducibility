#!/bin/bash
#SBATCH --job-name=sps_fi_coverage
#SBATCH --output=../../../logs/sps_fi_coverage_%j.out
#SBATCH --error=../../../logs/sps_fi_coverage_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

# Generate SPS feature_index=25 coverage figures for MCC dataset
# Compares SPS (feature_index=25) against other methods

set -e

# Set project paths
PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
export PROJECT_ROOT

echo "============================================"
echo "Starting SPS feature_index coverage analysis"
echo "Start time: $(date)"
echo "============================================"

# Create logs directory
mkdir -p "$PROJECT_ROOT/logs"

# Activate conda environment
source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh
conda activate facs_sampling
echo "Using Python: $(which python)"
echo ""

# Change to the script directory
cd "$PROJECT_ROOT/figures/revision/coverage"

echo "Running sps_feature_index_coverage.py..."
python sps_feature_index_coverage.py

echo ""
echo "============================================"
echo "Generated figures:"
echo "  - figures/revision/sps_fi25_mcc_coverage_single.jpg"
echo "  - figures/revision/sps_fi25_mcc_coverage.jpg"
echo ""
echo "Completed at $(date)"
echo "============================================"
