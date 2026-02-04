#!/bin/bash
#SBATCH --job-name=cell_classify_lcmv
#SBATCH --output=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/cell_classification_lcmv_%j.out
#SBATCH --error=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/cell_classification_lcmv_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --gres=localtmp:100G

# Experiment 2: Cell Type Classification Comparison (LCMV Dataset)
#
# This script compares classification performance when training on SPS vs Random
# sampled data, with evaluation on a common held-out test set.
#
# LCMV dataset has ~34M cells, so requires significant memory for:
# - Loading the full dataset
# - Creating train/test splits
# - Feature extraction (HVGs)
#
# Usage:
#   sbatch run_analysis_lcmv.sh
#   OR for interactive testing:
#   conda activate facs_sampling && python cell_type_classification_lcmv.py

# Exit on error
set -e

# Use local temp directory to avoid memory issues with large arrays
export TMPDIR="${LOCAL_TMPDIR:-/tmp}"
echo "Using TMPDIR: $TMPDIR"

# Load conda properly for SLURM
source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh

# Activate environment
conda activate facs_sampling

# Set script directory
SCRIPT_DIR="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/downstream_experiments/02_classification"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "Cell Type Classification Comparison (LCMV)"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "Node: $(hostname)"
echo "Memory requested: 256G"
echo "Local temp: 100G"
echo "=============================================="

# Run analysis with both classifiers, 5 random seeds
# Using 2000 highly variable genes (HVGs) as features
# Note: With 34M cells, this may take several hours
python "$SCRIPT_DIR/cell_type_classification_lcmv.py" \
    --classifiers logreg rf \
    --n-seeds 5 \
    --n-top-genes 2000

echo "=============================================="
echo "Analysis complete"
echo "End time: $(date)"
echo "=============================================="
