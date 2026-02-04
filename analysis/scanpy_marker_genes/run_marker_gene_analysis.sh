#!/bin/bash

# ============================================================================
# Scanpy Marker Gene Analysis - SLURM Job Submission Script
# ============================================================================
#
# This script runs the Scanpy marker gene analysis pipeline for:
# - Full dataset (ground truth)
# - SPS 100k subsample
# - Random 100k subsample
#
# Using multiple DE methods: Wilcoxon, t-test, logistic regression
#
# Usage:
#   bash run_marker_gene_analysis.sh [mode]
#
# Modes:
#   full    - Run analysis on full dataset only
#   sps     - Run analysis on SPS subsample only
#   random  - Run analysis on Random subsample only
#   all     - Run analysis on all three datasets (default)
#   compare - Run comparison script only (after all analyses complete)
#
# Example:
#   bash run_marker_gene_analysis.sh all      # Run all analyses
#   bash run_marker_gene_analysis.sh compare  # Run comparison after analyses
#
# ============================================================================

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
JOBS_DIR="$PROJECT_ROOT/jobs/scanpy_marker_genes"

# Create jobs directory if it doesn't exist
mkdir -p "$JOBS_DIR"

# Parse command line argument
MODE="${1:-all}"

# SLURM configuration
PARTITION="normal"
ACCOUNT="ohler"
TIME_FULL="8:00:00"      # Full dataset takes longer
TIME_SUBSAMPLE="4:00:00"  # Subsamples are faster
TIME_COMPARE="1:00:00"    # Comparison is quick
MEM_FULL="128G"          # Full dataset needs more memory
MEM_SUBSAMPLE="64G"      # Subsamples need less memory
MEM_COMPARE="16G"        # Comparison needs little memory
CPUS="4"

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
    local job_script="$JOBS_DIR/run_${method}.sh"
    
    cat > "$job_script" << EOF
#!/bin/bash
#SBATCH --job-name=scanpy_marker_${method}
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --time=${time}
#SBATCH --mem=${mem}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --output=${PROJECT_ROOT}/logs/scanpy_marker_${method}_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/scanpy_marker_${method}_%j.err

# Load conda
source ~/.bashrc
conda activate ${CONDA_ENV}

# Navigate to script directory
cd ${SCRIPT_DIR}

# Run the pipeline
echo "Starting Scanpy marker gene analysis for ${method}..."
echo "DE methods: ${DE_METHODS}"
echo "Start time: \$(date)"

python scanpy_marker_pipeline.py \\
    --method ${method} \\
    --de-methods ${DE_METHODS}

echo "End time: \$(date)"
echo "Done!"
EOF
    
    chmod +x "$job_script"
    echo "Created job script: $job_script"
}

create_compare_job_script() {
    local job_script="$JOBS_DIR/run_compare.sh"
    
    cat > "$job_script" << EOF
#!/bin/bash
#SBATCH --job-name=scanpy_marker_compare
#SBATCH --partition=${PARTITION}
#SBATCH --account=${ACCOUNT}
#SBATCH --time=${TIME_COMPARE}
#SBATCH --mem=${MEM_COMPARE}
#SBATCH --cpus-per-task=2
#SBATCH --output=${PROJECT_ROOT}/logs/scanpy_marker_compare_%j.out
#SBATCH --error=${PROJECT_ROOT}/logs/scanpy_marker_compare_%j.err

# Load conda
source ~/.bashrc
conda activate ${CONDA_ENV}

# Navigate to script directory
cd ${SCRIPT_DIR}

# Run comparison
echo "Starting marker gene comparison..."
echo "Start time: \$(date)"

python compare_marker_genes.py \\
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
    sbatch "$JOBS_DIR/run_${method}.sh"
}

case "$MODE" in
    "full")
        echo "Running analysis on full dataset..."
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
        submit_job "full" "$TIME_FULL" "$MEM_FULL"
        submit_job "sps" "$TIME_SUBSAMPLE" "$MEM_SUBSAMPLE"
        submit_job "random" "$TIME_SUBSAMPLE" "$MEM_SUBSAMPLE"
        echo ""
        echo "All jobs submitted. After they complete, run:"
        echo "  bash $0 compare"
        ;;
    "compare")
        echo "Running comparison script..."
        create_compare_job_script
        sbatch "$JOBS_DIR/run_compare.sh"
        ;;
    "local")
        echo "Running all analyses locally (not on SLURM)..."
        echo "This may take a long time!"
        
        source ~/.bashrc
        conda activate ${CONDA_ENV}
        
        cd "$SCRIPT_DIR"
        
        echo "Running full dataset..."
        python scanpy_marker_pipeline.py --method full --de-methods $DE_METHODS
        
        echo "Running SPS subsample..."
        python scanpy_marker_pipeline.py --method sps --de-methods $DE_METHODS
        
        echo "Running Random subsample..."
        python scanpy_marker_pipeline.py --method random --de-methods $DE_METHODS
        
        echo "Running comparison..."
        python compare_marker_genes.py --de-methods $DE_METHODS
        
        echo "Done!"
        ;;
    *)
        echo "Usage: $0 [full|sps|random|all|compare|local]"
        echo ""
        echo "Modes:"
        echo "  full    - Run analysis on full dataset only"
        echo "  sps     - Run analysis on SPS subsample only"
        echo "  random  - Run analysis on Random subsample only"
        echo "  all     - Run analysis on all three datasets (default)"
        echo "  compare - Run comparison script only (after all analyses complete)"
        echo "  local   - Run all analyses locally (not on SLURM)"
        exit 1
        ;;
esac

echo ""
echo "Job scripts saved to: $JOBS_DIR"
echo "Logs will be saved to: $PROJECT_ROOT/logs/"
