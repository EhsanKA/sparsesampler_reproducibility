#!/bin/bash

# Feature Index Classification for MCC_01 and MCC_05
# Runs RF classification for each feature index (1-30)
#
# Usage:
#   bash run_mcc01_mcc05_slurm.sh mcc_01
#   bash run_mcc01_mcc05_slurm.sh mcc_05

DATASET=${1:-mcc_01}

# Load cluster configuration
PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
SCRIPT_DIR="${PROJECT_ROOT}/jobs/feature_index_classification"

# Create log directories
mkdir -p ${PROJECT_ROOT}/logs/feature_index_classification_${DATASET}

# Define feature indices to test
FEATURE_INDICES=($(seq 1 30))
SIZE=100000

# Calculate total combinations
total_combinations=${#FEATURE_INDICES[@]}

echo "=============================================="
echo "Feature Index Classification Analysis"
echo "Dataset: ${DATASET}"
echo "=============================================="
echo "Feature indices: 1-30"
echo "Sample size: ${SIZE}"
echo "Total jobs: ${total_combinations}"
echo "=============================================="

# Submit the job array
sbatch <<EOF
#!/bin/bash

#SBATCH -J fi_clf_${DATASET}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ekarimi@mdc-berlin.de
#SBATCH -o ${PROJECT_ROOT}/logs/feature_index_classification_${DATASET}/output_%A_%a.log
#SBATCH -e ${PROJECT_ROOT}/logs/feature_index_classification_${DATASET}/error_%A_%a.log
#SBATCH -t 06:00:00
#SBATCH --mem=256G
#SBATCH -A ohler
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=localtmp:200G
#SBATCH --array=1-${total_combinations}

# Define feature indices inside job
FEATURE_INDICES=($(seq 1 30))
SIZE=${SIZE}
DATASET=${DATASET}

# Get feature index for this task
task_id=\$SLURM_ARRAY_TASK_ID
feature_index=\${FEATURE_INDICES[\$((task_id - 1))]}

# Use local temp directory
export TMPDIR="\${LOCAL_TMPDIR:-/tmp}"
echo "Using TMPDIR: \$TMPDIR"

# Load environment
source ~/.bashrc
conda activate facs_sampling

# Set working directory
cd ${SCRIPT_DIR}

# Set output directory
if [ "\${DATASET}" == "mcc" ]; then
    OUTPUT_DIR="${SCRIPT_DIR}/results"
else
    OUTPUT_DIR="${SCRIPT_DIR}/results_\${DATASET}"
fi

# Check if results already exist
output_file="\${OUTPUT_DIR}/summary_feature_index_\${feature_index}.csv"
if [ -f "\$output_file" ]; then
    echo "Results for feature_index=\${feature_index} already exist, skipping..."
    exit 0
fi

echo "=============================================="
echo "Feature Index Classification"
echo "=============================================="
echo "Dataset: \$DATASET"
echo "Task ID: \$task_id"
echo "Feature Index: \$feature_index"
echo "Size: \$SIZE"
echo "Node: \$(hostname)"
echo "Start time: \$(date)"
echo "=============================================="

# Run analysis
python classify_by_feature_index_unified.py \
    --dataset \$DATASET \
    --feature-index \$feature_index \
    --size \$SIZE \
    --rep 0

echo "=============================================="
echo "Completed"
echo "End time: \$(date)"
echo "=============================================="
EOF

echo ""
echo "Job array submitted for ${DATASET}!"
echo "Monitor with: squeue -u \$USER"
echo "Results will be in: ${SCRIPT_DIR}/results_${DATASET}/"
