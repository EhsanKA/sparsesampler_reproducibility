#!/bin/bash

# Load cluster configuration
if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found. Please copy from template and configure."
    echo "cp ../../config/cluster_config.sh.template ../../config/cluster_config.sh"
    exit 1
fi

# Define your combinations
DATASETS=('mcc' 'lcmv')
METHODS=('sps' 'random')
SIZES=(100000 200000)
REPS=(0)

# Calculate total combinations
total_combinations=$(( ${#DATASETS[@]} * ${#METHODS[@]} * ${#SIZES[@]} * ${#REPS[@]} ))

echo "Submitting DE analysis job array with $total_combinations combinations"
echo "Datasets: ${DATASETS[@]}"
echo "Methods: ${METHODS[@]}"
echo "Sizes: ${SIZES[@]}"
echo "Reps: ${REPS[@]}"

# Submit the job array for all combinations
sbatch <<EOF
#!/bin/bash

#SBATCH -J de_analysis_array
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/de_analysis/logs/output_de_%A_%a.stdlog
#SBATCH -e ${LOG_PATH}/de_analysis/logs/output_de_%A_%a.stderr
#SBATCH -t 12:00:00
#SBATCH --mem=100G
#SBATCH -A ohler
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-${total_combinations}

# Define the arrays inside the job script
DATASETS=('mcc' 'lcmv')
METHODS=('sps' 'random')
SIZES=(100000 200000)
REPS=(0)

# Compute task index based on SLURM_ARRAY_TASK_ID
task_id=\$SLURM_ARRAY_TASK_ID

# Calculate indices
rep_index=\$(( (task_id - 1) % \${#REPS[@]} ))
count=\$(( (task_id - 1) / \${#REPS[@]} ))

size_index=\$(( (count) %  \${#SIZES[@]} ))
count=\$(( (count) / \${#SIZES[@]} ))

method_index=\$(( (count) %  \${#METHODS[@]} ))
count=\$(( (count) /  \${#METHODS[@]} ))

dataset_index=\$(( (count) %  \${#DATASETS[@]} ))
count=\$(( (count) /  \${#DATASETS[@]} ))

# Assign values based on calculated indices
dataset=\${DATASETS[\$dataset_index]}
method=\${METHODS[\$method_index]}
size=\${SIZES[\$size_index]}
rep=\${REPS[\$rep_index]}

# Load environment
source ~/.bashrc
conda activate ${CONDA_ENV}

# Set PROJECT_ROOT environment variable
export PROJECT_ROOT=${PROJECT_ROOT}

# Change directory to the analysis folder
cd ${PROJECT_ROOT}/analysis

# File path to check if job is already done
output_file="${PROJECT_ROOT}/analysis/results/de_analysis/\${dataset}_\${method}_size\${size}_rep\${rep}_de_detailed.csv"

# Check if the file exists and skip if it's already done
if [ -f "\$output_file" ]; then
    echo "File \$output_file exists, skipping task \$task_id..."
    exit 0
fi

# Debugging info
echo "=========================================="
echo "DE Analysis Task \$task_id"
echo "=========================================="
echo "Dataset: \$dataset"
echo "Method: \$method"
echo "Size: \$size"
echo "Rep: \$rep"
echo "Node: \$(hostname)"
echo "CPU Info: \$(lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
echo "=========================================="

# Run the Python script with the computed parameters
python run_de_analysis.py --dataset \$dataset --method \$method --size \$size --rep \$rep

echo "Task \$task_id completed"
EOF

echo "Job array submitted successfully!"
echo "Monitor with: squeue -u \$USER"
