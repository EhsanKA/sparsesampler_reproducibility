#!/bin/bash
#SBATCH --job-name=scanpy_marker_compare
#SBATCH --partition=normal
#SBATCH --account=ohler
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/scanpy_marker_compare_%j.out
#SBATCH --error=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/scanpy_marker_compare_%j.err

# Load conda
source ~/.bashrc
conda activate facs_sampling

# Navigate to script directory
cd /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/scanpy_marker_genes

# Run comparison
echo "Starting marker gene comparison..."
echo "Start time: $(date)"

python compare_marker_genes.py \
    --de-methods wilcoxon t-test logreg

echo "End time: $(date)"
echo "Done!"
