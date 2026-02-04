#!/bin/bash

# Load cluster configuration
if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found. Please copy from template and configure."
    echo "cp ../../config/cluster_config.sh.template ../../config/cluster_config.sh"
    exit 1
fi

# Create logs directory if it doesn't exist
mkdir -p ${LOG_PATH}/test_sampling_methods/logs

# Submit a single job that processes all datasets and sizes
echo "Submitting job to create sampling methods tables for all datasets..."

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

# Load environment
source ~/.bashrc
conda activate facs_sampling

# Set PROJECT_ROOT environment variable
export PROJECT_ROOT=${PROJECT_ROOT}

# Change directory to the script folder
cd \${PROJECT_ROOT}/jobs/test_sampling_methods

# Run the table creation script (processes all datasets and sizes)
python create_sampling_methods_table.py

EOF

echo "Table creation job submitted successfully!"

