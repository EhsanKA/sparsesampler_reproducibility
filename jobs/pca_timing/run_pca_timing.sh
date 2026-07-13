#!/bin/bash

if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found"
    exit 1
fi

mkdir -p ${LOG_PATH}/pca_timing/logs

sbatch <<EOF
#!/bin/bash
#SBATCH -J pca_timing_benchmark
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/pca_timing/logs/output_%j.stdlog
#SBATCH -e ${LOG_PATH}/pca_timing/logs/error_%j.stderr
#SBATCH -t 02:00:00
#SBATCH --mem=200G
#SBATCH -A ohler
#SBATCH -p long
#SBATCH -C cascade-lake
#SBATCH --nodelist=maxg11,maxg12,maxg13,maxg14,maxg15,maxg16,maxg17,maxg18,maxg21,maxg22,maxg23,maxg24,maxg25,maxg26
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=localtmp:50G

source ~/.bashrc
conda activate ${CONDA_ENV}
export PROJECT_ROOT=${PROJECT_ROOT}
cd ${PROJECT_ROOT}/jobs/pca_timing
python pca_timing_benchmark.py
EOF
