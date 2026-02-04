#!/bin/bash
#SBATCH --job-name=scanpy_marker_full
#SBATCH --partition=normal
#SBATCH --account=ohler
#SBATCH --time=8:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --output=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/scanpy_marker_full_%j.out
#SBATCH --error=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/scanpy_marker_full_%j.err

# Load conda
source ~/.bashrc
conda activate facs_sampling

# Navigate to script directory
cd /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/scanpy_marker_genes

# Run the pipeline
echo "Starting Scanpy marker gene analysis for full..."
echo "DE methods: wilcoxon t-test logreg"
echo "Start time: $(date)"

python scanpy_marker_pipeline.py \
    --method full \
    --de-methods wilcoxon t-test logreg

echo "End time: $(date)"
echo "Done!"
