#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
SIZE=100000
total_combinations=25

if [ -f "${PROJECT_ROOT}/config/cluster_config.sh" ]; then
    source ${PROJECT_ROOT}/config/cluster_config.sh
else
    CLUSTER_EMAIL="ekarimi@mdc-berlin.de"
fi

mkdir -p ${PROJECT_ROOT}/logs/feature_index_classification_lcmv

sbatch <<EOF
#!/bin/bash
#SBATCH -J lcmv_fi_clf
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${PROJECT_ROOT}/logs/feature_index_classification_lcmv/output_%A_%a.log
#SBATCH -e ${PROJECT_ROOT}/logs/feature_index_classification_lcmv/error_%A_%a.log
#SBATCH -t 04:00:00
#SBATCH --mem=256G
#SBATCH -A ohler
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=localtmp:200G
#SBATCH --array=1-${total_combinations}

FEATURE_INDICES=($(seq 1 25))
feature_index=\${FEATURE_INDICES[\$((SLURM_ARRAY_TASK_ID - 1))]}

source ~/.bashrc
conda activate facs_sampling
cd ${SCRIPT_DIR}

output_file="${PROJECT_ROOT}/jobs/feature_index_classification_lcmv/results/summary_feature_index_\${feature_index}.csv"
if [ -f "\$output_file" ]; then
    exit 0
fi

python classify_by_feature_index_lcmv.py --feature-index \$feature_index --size ${SIZE} --rep 0
EOF
