#!/bin/bash

# =============================================================================
# SEMITONES Marker Gene Analysis - SLURM Job Script
# =============================================================================
#
# This script submits SLURM jobs to run the SEMITONES marker gene analysis
# for both SPS and Random sampling methods.
#
# Usage:
#   cd /fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility/analysis/semitones_marker_genes
#   bash run_marker_gene_analysis.sh [sps|random|both|compare]
#
# Arguments:
#   sps     - Run analysis for SPS method only
#   random  - Run analysis for Random method only
#   both    - Run analysis for both methods (default)
#   compare - Run comparison script (after both analyses complete)
#
# =============================================================================

# Load cluster configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
CONFIG_FILE="${PROJECT_ROOT}/config/cluster_config.sh"

if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    echo "Warning: cluster_config.sh not found. Using defaults."
    export CLUSTER_EMAIL="ehsan.karimiara@mdc-berlin.de"
    export PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
    export LOG_PATH="${PROJECT_ROOT}/logs"
fi

# Create log directory
LOG_DIR="${LOG_PATH}/semitones_marker_genes"
mkdir -p "$LOG_DIR"

# Parse arguments
METHOD="${1:-both}"

echo "========================================"
echo "SEMITONES Marker Gene Analysis"
echo "========================================"
echo "Project Root: ${PROJECT_ROOT}"
echo "Log Directory: ${LOG_DIR}"
echo "Method: ${METHOD}"
echo "========================================"

# Function to submit analysis job
submit_analysis_job() {
    local method=$1
    
    echo ""
    echo "Submitting ${method} analysis job..."
    
    sbatch <<EOF
#!/bin/bash

#SBATCH -J semitones_marker_${method}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_DIR}/marker_${method}_%j.out
#SBATCH -e ${LOG_DIR}/marker_${method}_%j.err
#SBATCH -t 24:00:00
#SBATCH --mem=200G
#SBATCH -A ohler
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=localtmp:100G

# =============================================================================
# SEMITONES Marker Gene Analysis - ${method^^}
# =============================================================================

echo "=========================================="
echo "SEMITONES Marker Gene Analysis - ${method^^}"
echo "=========================================="
echo "Node: \$(hostname)"
echo "Date: \$(date)"
echo "Job ID: \${SLURM_JOB_ID}"
echo "CPUs: \${SLURM_CPUS_PER_TASK}"
echo "Memory: 128GB"
echo "=========================================="

# Load conda
source ~/.bashrc

# Activate semitones environment
conda activate semitones_env

# Verify SEMITONES is available
if ! python -c "import SEMITONES" 2>/dev/null; then
    echo "ERROR: SEMITONES not available in semitones_env"
    echo "Please install: pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip"
    exit 1
fi

echo "SEMITONES version check passed"
echo "Python: \$(which python)"
echo "=========================================="

# Set environment variables
export PROJECT_ROOT="${PROJECT_ROOT}"

# Copy data to local scratch for faster I/O
echo "Copying data to local scratch (\$TMPDIR)..."
DATA_SOURCE="${PROJECT_ROOT}/data/mcc/benchmark/30"
LOCAL_DATA="\${TMPDIR}/mcc_data"
mkdir -p \${LOCAL_DATA}

echo "  Copying adata.h5ad..."
cp \${DATA_SOURCE}/adata.h5ad \${LOCAL_DATA}/
echo "  Copying sampling indices..."
cp -r \${DATA_SOURCE}/sps \${LOCAL_DATA}/
cp -r \${DATA_SOURCE}/random \${LOCAL_DATA}/
echo "  Data copied to \${LOCAL_DATA}"
ls -lh \${LOCAL_DATA}/

echo "=========================================="

# Change to analysis directory
cd ${SCRIPT_DIR}

echo "Working directory: \$(pwd)"
echo "=========================================="

# Run the analysis with local data path
echo "Starting ${method} analysis..."
python semitones_marker_pipeline.py \\
    --method ${method} \\
    --ncpu 8 \\
    --skip-permutation \\
    --data-path \${LOCAL_DATA}

exit_code=\$?

echo "=========================================="
if [ \$exit_code -eq 0 ]; then
    echo "${method^^} analysis completed successfully"
else
    echo "${method^^} analysis failed with exit code \$exit_code"
fi
echo "=========================================="

exit \$exit_code
EOF

    echo "Job submitted for ${method}"
}

# Function to submit comparison job
submit_comparison_job() {
    echo ""
    echo "Submitting comparison job..."
    
    sbatch <<EOF
#!/bin/bash

#SBATCH -J semitones_compare
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_DIR}/compare_%j.out
#SBATCH -e ${LOG_DIR}/compare_%j.err
#SBATCH -t 01:00:00
#SBATCH --mem=32G
#SBATCH -A ohler
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# =============================================================================
# SEMITONES Marker Gene Comparison
# =============================================================================

echo "=========================================="
echo "SEMITONES Marker Gene Comparison"
echo "=========================================="
echo "Node: \$(hostname)"
echo "Date: \$(date)"
echo "Job ID: \${SLURM_JOB_ID}"
echo "=========================================="

# Load conda
source ~/.bashrc

# Activate facs_sampling environment (for plotting)
conda activate facs_sampling

# Set environment variables
export PROJECT_ROOT="${PROJECT_ROOT}"

# Change to analysis directory
cd ${SCRIPT_DIR}

echo "Working directory: \$(pwd)"
echo "=========================================="

# Check if results exist
if [ ! -f "results/sps_marker_genes.pkl" ] || [ ! -f "results/random_marker_genes.pkl" ]; then
    echo "ERROR: Results not found. Please run both SPS and Random analyses first."
    exit 1
fi

# Run the comparison
echo "Starting comparison..."
python compare_marker_genes.py

exit_code=\$?

echo "=========================================="
if [ \$exit_code -eq 0 ]; then
    echo "Comparison completed successfully"
    echo ""
    echo "Results saved to:"
    echo "  ${SCRIPT_DIR}/results/marker_gene_comparison.csv"
    echo ""
    echo "Figures saved to:"
    echo "  ${SCRIPT_DIR}/results/figures/"
else
    echo "Comparison failed with exit code \$exit_code"
fi
echo "=========================================="

exit \$exit_code
EOF

    echo "Comparison job submitted"
}

# Main logic
case "$METHOD" in
    sps)
        submit_analysis_job "sps"
        ;;
    random)
        submit_analysis_job "random"
        ;;
    both)
        submit_analysis_job "sps"
        submit_analysis_job "random"
        echo ""
        echo "NOTE: After both jobs complete, run:"
        echo "  bash run_marker_gene_analysis.sh compare"
        ;;
    compare)
        submit_comparison_job
        ;;
    *)
        echo "Unknown method: $METHOD"
        echo "Usage: bash run_marker_gene_analysis.sh [sps|random|both|compare]"
        exit 1
        ;;
esac

echo ""
echo "========================================"
echo "Jobs submitted successfully!"
echo "========================================"
echo ""
echo "Monitor jobs with:"
echo "  squeue -u \$USER"
echo ""
echo "View logs in:"
echo "  ${LOG_DIR}/"
echo ""
echo "Results will be saved to:"
echo "  ${SCRIPT_DIR}/results/"
echo "========================================"
