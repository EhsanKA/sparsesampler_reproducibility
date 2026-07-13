#!/bin/bash

if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found"
    exit 1
fi

mkdir -p ${LOG_PATH}/test_sampling_methods/logs

sbatch <<EOF
#!/bin/bash
#SBATCH -J create_sampling_methods_table
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/test_sampling_methods/logs/create_sampling_methods_table_%A.stdlog
#SBATCH -e ${LOG_PATH}/test_sampling_methods/logs/create_sampling_methods_table_%A.stderr
#SBATCH -t 02:00:00
#SBATCH --mem=100G
#SBATCH -A ohler
#SBATCH -p long
#SBATCH -C cascade-lake
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate facs_sampling
export PROJECT_ROOT=${PROJECT_ROOT}
cd \${PROJECT_ROOT}/jobs/test_sampling_methods
python ../make_tables.py --type sampling_methods
EOF
