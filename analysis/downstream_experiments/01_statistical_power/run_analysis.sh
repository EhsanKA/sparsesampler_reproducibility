#!/bin/bash
#SBATCH --job-name=stat_power
#SBATCH --output=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/statistical_power_%j.out
#SBATCH --error=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/statistical_power_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Experiment 1: Statistical Power Analysis for Rare Cell Markers
# 
# This script compares the statistical power of SPS vs Random sampling
# for detecting marker genes of rare cell types (particularly osteoblast).
#
# Usage:
#   sbatch run_analysis.sh
#   OR for interactive testing:
#   conda activate facs_sampling && python statistical_power_analysis.py

# Exit on error
set -e

# Load conda properly for SLURM
source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh

# Activate environment
conda activate facs_sampling

# Set script directory
SCRIPT_DIR="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/downstream_experiments/01_statistical_power"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "Statistical Power Analysis"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "=============================================="

# Run analysis with all DE methods
python "$SCRIPT_DIR/statistical_power_analysis.py" \
    --de-methods wilcoxon t-test logreg \
    --top-n 50

echo "=============================================="
echo "Analysis complete"
echo "End time: $(date)"
echo "=============================================="
