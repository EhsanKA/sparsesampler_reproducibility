#!/bin/bash
#SBATCH --job-name=scanpy_marker_lcmv_random
#SBATCH --partition=normal
#SBATCH --account=ohler
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=localtmp:100G
#SBATCH --output=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/scanpy_marker_lcmv_random_%j.out
#SBATCH --error=/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/scanpy_marker_lcmv_random_%j.err

# Use local temp directory to avoid memory issues
export TMPDIR="${LOCAL_TMPDIR:-/tmp}"
echo "Using TMPDIR: $TMPDIR"

# Load conda
source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh
conda activate facs_sampling

# Navigate to script directory
cd /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/scanpy_marker_genes

# Run the pipeline
echo "Starting Scanpy marker gene analysis (LCMV) for random..."
echo "DE methods: wilcoxon t-test logreg"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo "Memory: 64G"

python scanpy_marker_pipeline_lcmv.py \
    --method random \
    --de-methods wilcoxon t-test logreg

echo "End time: $(date)"
echo "Done!"
