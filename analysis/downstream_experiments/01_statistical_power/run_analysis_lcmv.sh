#!/bin/bash
#SBATCH --job-name=stat_power_lcmv
#SBATCH --output=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/statistical_power_lcmv_%j.out
#SBATCH --error=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/statistical_power_lcmv_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal
#SBATCH --gres=localtmp:100G

# Experiment 1: Statistical Power Analysis for Rare Cell Markers (LCMV Dataset)
# 
# This script compares the statistical power of SPS vs Random sampling
# for detecting marker genes of rare cell types in the LCMV dataset.
#
# Pre-requisites:
#   - Run scanpy marker gene analysis first:
#     bash ../../../scanpy_marker_genes/run_marker_gene_analysis_lcmv.sh subsamples
#
# Usage:
#   sbatch run_analysis_lcmv.sh
#   OR for interactive testing:
#   conda activate facs_sampling && python statistical_power_analysis_lcmv.py

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
SCRIPT_DIR="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/downstream_experiments/01_statistical_power"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "Statistical Power Analysis (LCMV Dataset)"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "Node: $(hostname)"
echo "Memory requested: 128G"
echo "=============================================="

# Check if marker gene results exist
RESULTS_DIR="$SCRIPT_DIR/../../scanpy_marker_genes/results_lcmv"
if [ ! -d "$RESULTS_DIR" ]; then
    echo "ERROR: Marker gene results not found at $RESULTS_DIR"
    echo "Please run the scanpy marker gene analysis first:"
    echo "  bash ../../../scanpy_marker_genes/run_marker_gene_analysis_lcmv.sh subsamples"
    exit 1
fi

# Check for required files
for method in sps random; do
    for de_method in wilcoxon t-test logreg; do
        file="$RESULTS_DIR/${method}_${de_method}_marker_genes.csv"
        if [ ! -f "$file" ]; then
            echo "WARNING: Missing file: $file"
        fi
    done
done

# Run analysis with all DE methods
python "$SCRIPT_DIR/statistical_power_analysis_lcmv.py" \
    --de-methods wilcoxon t-test logreg \
    --top-n 50

echo "=============================================="
echo "Analysis complete"
echo "End time: $(date)"
echo "=============================================="
