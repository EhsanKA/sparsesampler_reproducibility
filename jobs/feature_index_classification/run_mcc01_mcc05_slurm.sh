#!/bin/bash

DATASET=${1:-mcc_01}
PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
SCRIPT_DIR="${PROJECT_ROOT}/jobs/feature_index_classification"
SIZE=100000
total_combinations=30

mkdir -p ${PROJECT_ROOT}/logs/feature_index_classification_${DATASET}

sbatch <<EOF
#!/bin/bash
#SBATCH -J fi_clf_${DATASET}
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

FEATURE_INDICES=($(seq 1 30))
feature_index=\${FEATURE_INDICES[\$((SLURM_ARRAY_TASK_ID - 1))]}

source ~/.bashrc
conda activate facs_sampling
cd ${SCRIPT_DIR}

if [ "${DATASET}" == "mcc" ]; then
    OUTPUT_DIR="${SCRIPT_DIR}/results"
else
    OUTPUT_DIR="${SCRIPT_DIR}/results_${DATASET}"
fi

if [ -f "\${OUTPUT_DIR}/summary_feature_index_\${feature_index}.csv" ]; then
    exit 0
fi

python classify_by_feature_index_unified.py \
    --dataset ${DATASET} \
    --feature-index \$feature_index \
    --size ${SIZE} \
    --rep 0
EOF
