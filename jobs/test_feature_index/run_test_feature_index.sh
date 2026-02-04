#!/bin/bash

# Load cluster configuration
if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found. Please copy from template and configure."
    echo "cp ../../config/cluster_config.sh.template ../../config/cluster_config.sh"
    exit 1
fi

# Define datasets and their references
# We'll create separate job arrays for each dataset
DATASETS=('lcmv' 'mcc' 'mcc_01' 'mcc_05')

# LCMV references
LCMV_REFS=(1 5 10 20 34)

# MCC references (same for mcc, mcc_01, mcc_05)
MCC_REFS=(5 10 20 25 30)

# Feature indices to test (1 to 30)
FEATURE_INDICES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)

# Calculate total combinations for each dataset
# For each dataset, we have: num_refs combinations
# We'll use --all flag to process all feature_index and size combinations per ref

# Function to submit jobs for a dataset
submit_dataset_jobs() {
    local dataset=$1
    local refs_array=$2
    
    # Get reference array name
    local refs_var_name="${dataset^^}_REFS"
    if [ "$dataset" == "lcmv" ]; then
        local refs=("${LCMV_REFS[@]}")
    else
        local refs=("${MCC_REFS[@]}")
    fi
    
    local total_refs=${#refs[@]}
    
    echo "Submitting jobs for dataset: $dataset with ${total_refs} references"
    
    # Submit the job array for all references
    # Each job processes one (dataset, reference) combination
    # The job loads the reference data once and processes all feature_index and size combinations
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

# Define the references array inside the job script
$(if [ "$dataset" == "lcmv" ]; then
    echo "REFS=(1 5 10 20 34)"
else
    echo "REFS=(5 10 20 25 30)"
fi)

# Compute task index based on SLURM_ARRAY_TASK_ID
task_id=\$SLURM_ARRAY_TASK_ID
ref_index=\$((task_id - 1))
ref=\${REFS[\$ref_index]}

# Load environment
source ~/.bashrc
conda activate facs_sampling

# Set PROJECT_ROOT environment variable
export PROJECT_ROOT=${PROJECT_ROOT}

# Change directory to the script folder
cd ${PROJECT_ROOT}/jobs/test_feature_index

# Debugging info
echo "Running with Dataset=${dataset}, Reference=\$ref"
echo "Node: \$(hostname)"
echo "CPU Info: \$(lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
echo "CPU Min MHz: \$(lscpu | grep 'CPU min MHz' | cut -d: -f2 | xargs)"
echo "CPU Max MHz: \$(lscpu | grep 'CPU max MHz' | cut -d: -f2 | xargs)"
echo "CPU Current MHz: \$(lscpu | grep '^CPU MHz' | head -1 | cut -d: -f2 | xargs)"

# Run the Python script with --all flag to process all feature_index and size combinations
# This loads the reference data once and processes all combinations efficiently
python test_feature_index.py --dataset ${dataset} --ref \$ref --all

EOF
}

# Submit jobs for each dataset
for dataset in "${DATASETS[@]}"; do
    submit_dataset_jobs "$dataset"
    sleep 1  # Small delay between submissions
done

echo "All jobs submitted successfully!"

