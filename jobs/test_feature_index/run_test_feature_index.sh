#!/bin/bash

if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found"
    exit 1
fi

mkdir -p ${LOG_PATH}/test_feature_index/logs

for dataset in lcmv mcc mcc_01 mcc_05; do
    if [ "$dataset" == "lcmv" ]; then
        total_refs=5
    else
        total_refs=5
    fi

    sbatch <<EOF
#!/bin/bash
#SBATCH -J test_feature_index_${dataset}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/test_feature_index/logs/output_${dataset}_%A_%a.stdlog
#SBATCH -e ${LOG_PATH}/test_feature_index/logs/error_${dataset}_%A_%a.stderr
#SBATCH -t 16:00:00
#SBATCH --mem=200G
#SBATCH -A ohler
#SBATCH -p long
#SBATCH -C cascade-lake
#SBATCH --nodelist=maxg11,maxg12,maxg13,maxg14,maxg15,maxg16,maxg17,maxg18,maxg21,maxg22,maxg23,maxg24,maxg25,maxg26
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=localtmp:50G
#SBATCH --array=1-${total_refs}

$(if [ "$dataset" == "lcmv" ]; then echo "REFS=(1 5 10 20 34)"; else echo "REFS=(5 10 20 25 30)"; fi)
ref=\${REFS[\$((SLURM_ARRAY_TASK_ID - 1))]}

source ~/.bashrc
conda activate facs_sampling
export PROJECT_ROOT=${PROJECT_ROOT}
cd ${PROJECT_ROOT}/jobs/test_feature_index
python test_feature_index.py --dataset ${dataset} --ref \$ref --all
EOF
    sleep 1
done
