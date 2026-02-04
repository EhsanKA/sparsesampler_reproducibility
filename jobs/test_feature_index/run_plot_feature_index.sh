#!/bin/bash

# Load cluster configuration
if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found. Please copy from template and configure."
    echo "cp ../../config/cluster_config.sh.template ../../config/cluster_config.sh"
    exit 1
fi

# Submit the plotting job
sbatch <<EOF
#!/bin/bash

#SBATCH -J plot_feature_index
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/test_feature_index/logs/plot_feature_index_%A.stdlog
#SBATCH -e ${LOG_PATH}/test_feature_index/logs/plot_feature_index_%A.stderr
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
cd ${PROJECT_ROOT}/jobs/test_feature_index

# Run the plotting script
python plot_feature_index_results.py

EOF

echo "Plotting job submitted successfully!"


