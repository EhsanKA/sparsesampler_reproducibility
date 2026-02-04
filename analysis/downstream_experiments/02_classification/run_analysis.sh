#!/bin/bash
#SBATCH --job-name=cell_classify
#SBATCH --output=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/cell_classification_%j.out
#SBATCH --error=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/cell_classification_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --gres=localtmp:100G

# Experiment 2: Cell Type Classification Comparison
#
# This script compares classification performance when training on SPS vs Random
# sampled data, with evaluation on a common held-out test set.
#
# Usage:
#   sbatch run_analysis.sh
#   OR for interactive testing:
#   conda activate facs_sampling && python cell_type_classification.py

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
echo "Cell Type Classification Comparison"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "Node: $(hostname)"
echo "Memory requested: 128G"
echo "Local temp: 100G"
echo "=============================================="

# Run analysis with both classifiers, 5 random seeds
# Using 2000 highly variable genes (HVGs) as features
# HVGs are selected by scanpy's highly_variable_genes() based on variance
python "$SCRIPT_DIR/cell_type_classification.py" \
    --classifiers logreg rf \
    --n-seeds 5 \
    --n-top-genes 2000

echo "=============================================="
echo "Analysis complete"
echo "End time: $(date)"
echo "=============================================="
