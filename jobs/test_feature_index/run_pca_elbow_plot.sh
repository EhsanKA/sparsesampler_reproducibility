#!/bin/bash
#SBATCH --job-name=pca_elbow
#SBATCH --output=logs/pca_elbow_%j.out
#SBATCH --error=logs/pca_elbow_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal

# Load cluster config
source "$(dirname "$0")/../cluster_config.sh"

# Create logs directory
mkdir -p logs

# Activate conda
source ~/.bashrc
conda activate facs_sampling

echo "Starting PCA elbow plot..."
echo "Date: $(date)"
echo "Node: $(hostname)"

# Run the script
python plot_pca_elbow.py

echo "Completed: $(date)"
