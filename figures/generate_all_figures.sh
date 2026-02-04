#!/bin/bash
#SBATCH --job-name=gen_figures
#SBATCH --output=logs/gen_figures_%j.out
#SBATCH --error=logs/gen_figures_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

# Generate all supplementary figures for the paper
# This script generates:
# - supp_figure3_mcc_coverage.jpg (MCC, 1% rarity with distance > 75th pct)
# - supp_figure4_lcmv_coverage.jpg (LCMV, 1% rarity with distance > 75th pct)
# - supp_figure5-6_mcc_time.jpg (MCC timing)
# - supp_figure5-6_lcmv_time.jpg (LCMV timing)
# - supp_figure3_mcc_01_coverage.jpg (MCC_01, 0.1% rarity)
# - supp_figure3_mcc_05_coverage.jpg (MCC_05, 0.5% rarity)

set -e

# Set project paths explicitly (for SLURM compatibility)
PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
SCRIPT_DIR="$PROJECT_ROOT/figures"

echo "Project root: $PROJECT_ROOT"
export PROJECT_ROOT

# Create logs directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/logs"

# Activate conda environment
source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh
conda activate sps
echo "Using Python: $(which python)"
echo "Python version: $(python --version)"

echo "============================================"
echo "Starting figure generation at $(date)"
echo "============================================"

# Generate main figures (MCC and LCMV with 1% rarity + 75th percentile distance)
echo ""
echo "--- Generating MCC and LCMV coverage figures (main) ---"
echo "Using rare cell definition: distance > 75th percentile AND frequency < 1%"
cd "$SCRIPT_DIR/main/coverage"
python coverage_combined.py

echo ""
echo "--- Generating MCC and LCMV time performance figures (main) ---"
cd "$SCRIPT_DIR/main/time_performance"
python combined_time_plot.py

# Generate revision figures (MCC_01 with 0.1% and MCC_05 with 0.5%)
echo ""
echo "--- Generating MCC_01 and MCC_05 coverage figures (revision) ---"
echo "Using rare cell definition: distance > 75th percentile AND frequency < threshold"
echo "  - MCC_01: frequency < 0.1%"
echo "  - MCC_05: frequency < 0.5%"
cd "$SCRIPT_DIR/revision/coverage"
python coverage_combined.py

echo ""
echo "============================================"
echo "Figure generation completed at $(date)"
echo "============================================"
echo ""
echo "Generated figures:"
echo "  Main (in figures/):"
echo "    - figure1_combined_coverage.jpg"
echo "    - supp_figure3_mcc_coverage.jpg"
echo "    - supp_figure4_lcmv_coverage.jpg"
echo "    - supp_figure5-6_mcc_time.jpg"
echo "    - supp_figure5-6_lcmv_time.jpg"
echo ""
echo "  Revision (in figures/revision/):"
echo "    - figure1_combined_coverage.jpg"
echo "    - supp_figure3_mcc_01_coverage.jpg"
echo "    - supp_figure3_mcc_05_coverage.jpg"
