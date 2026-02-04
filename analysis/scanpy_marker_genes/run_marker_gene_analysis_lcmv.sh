#!/bin/bash

# ============================================================================
# Scanpy Marker Gene Analysis for LCMV Dataset - SLURM Job Submission Script
# ============================================================================
#
# This script runs the Scanpy marker gene analysis pipeline for LCMV dataset:
# - Full dataset (ground truth)
# - SPS 100k subsample
# - Random 100k subsample
#
# Using multiple DE methods: Wilcoxon, t-test, logistic regression
#
# Usage:
#   bash run_marker_gene_analysis_lcmv.sh [mode]
#
# Modes:
#   full    - Run analysis on full dataset only
#   sps     - Run analysis on SPS subsample only
#   random  - Run analysis on Random subsample only
#   all     - Run analysis on all three datasets (default)
#
# Example:
#   bash run_marker_gene_analysis_lcmv.sh all      # Run all analyses
#
# ============================================================================

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
JOBS_DIR="$PROJECT_ROOT/jobs/scanpy_marker_genes_lcmv"

# Create jobs directory if it doesn't exist
mkdir -p "$JOBS_DIR"
mkdir -p "$PROJECT_ROOT/logs"

# Parse command line argument
MODE="${1:-all}"

# SLURM configuration - LCMV is a much larger dataset
PARTITION="normal"
ACCOUNT="ohler"
TIME_FULL="24:00:00"      # Full dataset is ~34M cells, needs long time
TIME_SUBSAMPLE="4:00:00"  # Subsamples are much faster
MEM_FULL="256G"           # Full dataset needs lots of memory
MEM_SUBSAMPLE="64G"       # Subsamples need less memory
CPUS="8"
LOCAL_TMP="100G"          # Local temp storage

# DE methods to run
DE_METHODS="wilcoxon t-test logreg"

# Conda environment
CONDA_ENV="facs_sampling"

# ============================================================================
# Create SLURM job scripts
# ============================================================================

create_job_script() {
    local method=$1
    local time=$2
    local mem=$3
    local job_script="$JOBS_DIR/run_lcmv_${method}.sh"
    
    cat > "$job_script" << EOF
#!/bin/bash
#SBATCH --job-name=scanpy_marker_lcmv_${method}
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --time=${time}
#SBATCH --mem=${mem}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --gres=localtmp:${LOCAL_TMP}
#SBATCH --output=${PROJECT_ROOT}/logs/scanpy_marker_lcmv_${method}_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/scanpy_marker_lcmv_${method}_%j.err

# Use local temp directory to avoid memory issues
export TMPDIR="\${LOCAL_TMPDIR:-/tmp}"
echo "Using TMPDIR: \$TMPDIR"

# Load conda
source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Navigate to script directory
cd ${SCRIPT_DIR}

# Run the pipeline
echo "Starting Scanpy marker gene analysis (LCMV) for ${method}..."
echo "DE methods: ${DE_METHODS}"
echo "Start time: \$(date)"
echo "Node: \$(hostname)"
echo "Memory: ${mem}"

python scanpy_marker_pipeline_lcmv.py \\
    --method ${method} \\
    --de-methods ${DE_METHODS}

echo "End time: \$(date)"
echo "Done!"
EOF
    
    chmod +x "$job_script"
    echo "Created job script: $job_script"
}

# ============================================================================
# Submit jobs based on mode
# ============================================================================

submit_job() {
    local method=$1
    local time=$2
    local mem=$3
    
    create_job_script "$method" "$time" "$mem"
    
    echo "Submitting job for ${method}..."
    sbatch "$JOBS_DIR/run_lcmv_${method}.sh"
}

case "$MODE" in
    "full")
        echo "Running analysis on full LCMV dataset..."
        submit_job "full" "$TIME_FULL" "$MEM_FULL"
        ;;
    "sps")
        echo "Running analysis on SPS subsample..."
        submit_job "sps" "$TIME_SUBSAMPLE" "$MEM_SUBSAMPLE"
        ;;
    "random")
        echo "Running analysis on Random subsample..."
        submit_job "random" "$TIME_SUBSAMPLE" "$MEM_SUBSAMPLE"
        ;;
    "all")
        echo "Running analysis on all datasets..."
        # Note: Full dataset analysis takes very long and lots of memory
        # Consider running subsamples first to verify pipeline works
        submit_job "full" "$TIME_FULL" "$MEM_FULL"
        submit_job "sps" "$TIME_SUBSAMPLE" "$MEM_SUBSAMPLE"
        submit_job "random" "$TIME_SUBSAMPLE" "$MEM_SUBSAMPLE"
        echo ""
        echo "All jobs submitted!"
        echo ""
        echo "NOTE: Full dataset analysis may take 12-24 hours due to ~34M cells"
        ;;
    "subsamples")
        echo "Running analysis on subsamples only (SPS and Random)..."
        submit_job "sps" "$TIME_SUBSAMPLE" "$MEM_SUBSAMPLE"
        submit_job "random" "$TIME_SUBSAMPLE" "$MEM_SUBSAMPLE"
        echo ""
        echo "Subsample jobs submitted!"
        ;;
    "local")
        echo "Running all analyses locally (not on SLURM)..."
        echo "This may take a VERY long time!"
        echo "Consider using 'subsamples' mode first."
        
        source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh
        conda activate ${CONDA_ENV}
        
        cd "$SCRIPT_DIR"
        
        echo "Running SPS subsample..."
        python scanpy_marker_pipeline_lcmv.py --method sps --de-methods $DE_METHODS
        
        echo "Running Random subsample..."
        python scanpy_marker_pipeline_lcmv.py --method random --de-methods $DE_METHODS
        
        echo "Running full dataset (this will take a long time)..."
        python scanpy_marker_pipeline_lcmv.py --method full --de-methods $DE_METHODS
        
        echo "Done!"
        ;;
    *)
        echo "Usage: $0 [full|sps|random|all|subsamples|local]"
        echo ""
        echo "Modes:"
        echo "  full       - Run analysis on full LCMV dataset only (~34M cells, needs 256GB, ~24h)"
        echo "  sps        - Run analysis on SPS subsample only (100k cells)"
        echo "  random     - Run analysis on Random subsample only (100k cells)"
        echo "  all        - Run analysis on all three datasets (default)"
        echo "  subsamples - Run analysis on SPS and Random only (recommended first)"
        echo "  local      - Run all analyses locally (not on SLURM)"
        exit 1
        ;;
esac

echo ""
echo "Job scripts saved to: $JOBS_DIR"
echo "Logs will be saved to: $PROJECT_ROOT/logs/"
