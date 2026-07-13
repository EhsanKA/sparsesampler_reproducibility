#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"

if [ -f "${PROJECT_ROOT}/config/cluster_config.sh" ]; then
    source ${PROJECT_ROOT}/config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found"
    exit 1
fi

total_combinations=30
SIZE=100000
mkdir -p ${PROJECT_ROOT}/logs/feature_index_classification

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

FEATURE_INDICES=($(seq 1 30))
feature_index=\${FEATURE_INDICES[\$((SLURM_ARRAY_TASK_ID - 1))]}
SIZE=${SIZE}

source ~/.bashrc
conda activate facs_sampling
cd ${SCRIPT_DIR}

output_file="${PROJECT_ROOT}/jobs/feature_index_classification/results/summary_feature_index_\${feature_index}.csv"
if [ -f "\$output_file" ]; then
    exit 0
fi

python classify_by_feature_index.py --feature-index \$feature_index --size \$SIZE --rep 0
EOF
