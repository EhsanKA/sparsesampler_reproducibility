#!/bin/bash

if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found"
    exit 1
fi

total_jobs=11
max_array_index=$((total_jobs - 1))
mkdir -p ${LOG_PATH}/test_feature_index/logs

sbatch <<EOF
#!/bin/bash
#SBATCH -J create_feature_index_tables
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/test_feature_index/logs/create_table_%A_%a.stdlog
#SBATCH -e ${LOG_PATH}/test_feature_index/logs/create_table_%A_%a.stderr
#SBATCH --array=0-${max_array_index}
#SBATCH -t 01:00:00
#SBATCH --mem=200G
#SBATCH -A ohler
#SBATCH -p long
#SBATCH -C cascade-lake
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --gres=localtmp:50G
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate facs_sampling
export PROJECT_ROOT=${PROJECT_ROOT}
cd \${PROJECT_ROOT}/jobs/test_feature_index
python ../make_tables.py --type feature_index --array-index \${SLURM_ARRAY_TASK_ID}
EOF
