#!/bin/bash

# Load cluster configuration
if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found. Please copy from template and configure."
    echo "cp ../../config/cluster_config.sh.template ../../config/cluster_config.sh"
    exit 1
fi

# Create log directory
mkdir -p ${LOG_PATH}/pca_timing/logs

# Submit the job
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

# Load environment
source ~/.bashrc
conda activate ${CONDA_ENV}

# Set PROJECT_ROOT environment variable
export PROJECT_ROOT=${PROJECT_ROOT}

# Change directory to the script folder
cd ${PROJECT_ROOT}/jobs/pca_timing

# Debugging info
echo "Running PCA Timing Benchmark"
echo "Node: \$(hostname)"
echo "CPU Info: \$(lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
echo "CPU Min MHz: \$(lscpu | grep 'CPU min MHz' | cut -d: -f2 | xargs)"
echo "CPU Max MHz: \$(lscpu | grep 'CPU max MHz' | cut -d: -f2 | xargs)"
echo "CPU Current MHz: \$(lscpu | grep '^CPU MHz' | head -1 | cut -d: -f2 | xargs)"
echo ""

# Run the Python script
python pca_timing_benchmark.py

EOF
