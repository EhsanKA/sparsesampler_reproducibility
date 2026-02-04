#!/bin/bash

#SBATCH -J generate_mcc_01_figure
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ekarimi@mdc-berlin.de
#SBATCH -o /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/generate_mcc_01_figure_%j.out
#SBATCH -e /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs/generate_mcc_01_figure_%j.err
#SBATCH -t 12:00:00
#SBATCH --mem=200G
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load cluster configuration
if [ -f "../../../config/cluster_config.sh" ]; then
    source ../../../config/cluster_config.sh
else
    echo "Warning: cluster_config.sh not found. Using defaults."
fi

# Create logs directory
mkdir -p /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/logs

# Load environment
source ~/.bashrc
conda activate ${CONDA_ENV:-facs_sampling}

# Set absolute path to project root
export PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"

# Change to project directory
cd "$PROJECT_ROOT"

# Run the script to generate mcc_01 figure
echo "=========================================="
echo "Generating supp_figure3_mcc_01_coverage.jpg"
echo "=========================================="
echo "Project root: ${PROJECT_ROOT}"
echo "Conda environment: ${CONDA_ENV:-facs_sampling}"
echo "Date: $(date)"
echo "=========================================="

python figures/revision/coverage/coverage_combined.py

echo "=========================================="
echo "Figure generation completed at $(date)"
echo "=========================================="

