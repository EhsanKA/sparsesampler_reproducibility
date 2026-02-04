#!/bin/bash

# Feature Index Classification Analysis
# Runs RF classification for each feature index (1-30) to find the best one
#
# Usage:
#   bash run_all_feature_indices.sh

# Load cluster configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"

if [ -f "${PROJECT_ROOT}/config/cluster_config.sh" ]; then
    source ${PROJECT_ROOT}/config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found"
    exit 1
fi

# Define feature indices to test
FEATURE_INDICES=($(seq 1 30))
SIZE=100000

# Calculate total combinations
total_combinations=${#FEATURE_INDICES[@]}

echo "=============================================="
echo "Feature Index Classification Analysis"
echo "=============================================="
echo "Feature indices: 1-30"
echo "Sample size: ${SIZE}"
echo "Total jobs: ${total_combinations}"
echo "=============================================="

# Submit the job array
sbatch <<EOF
#!/bin/bash

#SBATCH -J feature_idx_clf
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${PROJECT_ROOT}/logs/feature_index_classification/output_%A_%a.log
#SBATCH -e ${PROJECT_ROOT}/logs/feature_index_classification/error_%A_%a.log
#SBATCH -t 04:00:00
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

# Check if results already exist
output_file="${PROJECT_ROOT}/jobs/feature_index_classification/results/summary_feature_index_\${feature_index}.csv"
if [ -f "\$output_file" ]; then
    echo "Results for feature_index=\${feature_index} already exist, skipping..."
    exit 0
fi

echo "=============================================="
echo "Feature Index Classification"
echo "=============================================="
echo "Task ID: \$task_id"
echo "Feature Index: \$feature_index"
echo "Size: \$SIZE"
echo "Node: \$(hostname)"
echo "Start time: \$(date)"
echo "=============================================="

# Run analysis
python classify_by_feature_index.py \
    --feature-index \$feature_index \
    --size \$SIZE \
    --rep 0

echo "=============================================="
echo "Completed"
echo "End time: \$(date)"
echo "=============================================="
EOF

echo ""
echo "Job array submitted!"
echo "Monitor with: squeue -u \$USER"
echo "Results will be in: ${PROJECT_ROOT}/jobs/feature_index_classification/results/"
