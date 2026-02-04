#!/bin/bash

# Load cluster configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "${SCRIPT_DIR}/../../config/cluster_config.sh" ]; then
    source "${SCRIPT_DIR}/../../config/cluster_config.sh"
else
    echo "Error: cluster_config.sh not found"
    exit 1
fi

# Create logs directory
mkdir -p ${LOG_PATH}/test_feature_index/logs

# Submit job
sbatch <<EOF
#!/bin/bash
#SBATCH -J evr_coverage_plot
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/test_feature_index/logs/evr_coverage_plot_%j.out
#SBATCH -e ${LOG_PATH}/test_feature_index/logs/evr_coverage_plot_%j.err
#SBATCH -t 01:00:00
#SBATCH --mem=32G
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Load environment
source ~/.bashrc
conda activate ${CONDA_ENV}

# Set PROJECT_ROOT environment variable
export PROJECT_ROOT=${PROJECT_ROOT}

# Change to project directory
cd ${PROJECT_ROOT}

# Run the EVR coverage plot script
python jobs/test_feature_index/plot_evr_from_tables_with_coverage.py

echo "Job completed at \$(date)"
EOF

echo "Job submitted!"
