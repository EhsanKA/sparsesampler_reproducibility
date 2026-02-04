#!/bin/bash

# Load cluster configuration
if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found. Please copy from template and configure."
    echo "cp ../../config/cluster_config.sh.template ../../config/cluster_config.sh"
    exit 1
fi

# Define datasets and their sizes
declare -A DATASET_SIZES
DATASET_SIZES[lcmv]="50000 100000 200000"
DATASET_SIZES[mcc]="50000 100000 200000 300000"
DATASET_SIZES[mcc_01]="50000 100000 200000 300000"
DATASET_SIZES[mcc_05]="50000 100000 200000 300000"

# Calculate total number of combinations
total_jobs=0
for dataset in lcmv mcc mcc_01 mcc_05; do
    sizes=${DATASET_SIZES[$dataset]}
    for size in $sizes; do
        total_jobs=$((total_jobs + 1))
    done
done

# Calculate max array index (0-indexed, so max is total_jobs - 1)
max_array_index=$((total_jobs - 1))

# Create logs directory if it doesn't exist
mkdir -p ${LOG_PATH}/test_feature_index/logs

# Submit job array
echo "Submitting job array with ${total_jobs} tasks (one task per table)..."

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

# Load environment
source ~/.bashrc
conda activate facs_sampling

# Set PROJECT_ROOT environment variable
export PROJECT_ROOT=${PROJECT_ROOT}

# Change directory to the script folder
cd \${PROJECT_ROOT}/jobs/test_feature_index

# Run the table creation script using array index
# SLURM_ARRAY_TASK_ID is 0-indexed and maps to dataset/size combinations
python create_feature_index_table.py --array-index \${SLURM_ARRAY_TASK_ID}

EOF

echo "Job array submitted successfully!"
echo "Array ID: Will be shown after submission"
echo "Each task will create one table: {dataset}_feature_index_table_size_{size}.csv"
echo "Total tasks: ${total_jobs}"
