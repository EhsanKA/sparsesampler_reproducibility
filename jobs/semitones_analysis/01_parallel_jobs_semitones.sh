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
# Focus on SPS samples as references (as suggested by boss)
# Using 100k samples for SEMITONES reference cells
DATASETS=('mcc_01' 'mcc_05' 'lcmv')
METHODS=('sps' 'random')  # Compare SPS vs Random as references
SIZES=(100000)  # Using 100k samples as reference cells
REPS=(0)

# Calculate total combinations
total_combinations=$(( ${#DATASETS[@]} * ${#METHODS[@]} * ${#SIZES[@]} * ${#REPS[@]} ))

echo "Submitting SEMITONES analysis job array with $total_combinations combinations"
echo "Datasets: ${DATASETS[@]}"
echo "Methods: ${METHODS[@]}"
echo "Sizes: ${SIZES[@]}"
echo "Reps: ${REPS[@]}"
echo ""
echo "Note: Using semitones_env conda environment"
echo "Note: Sampled cells will be used as REFERENCE cells in SEMITONES"

# Submit the job array for all combinations
sbatch <<EOF
#!/bin/bash

#SBATCH -J semitones_analysis
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/semitones_analysis/logs/output_semitones_%A_%a.stdlog
#SBATCH -e ${LOG_PATH}/semitones_analysis/logs/output_semitones_%A_%a.stderr
#SBATCH -t 24:00:00
#SBATCH --mem=150G
#SBATCH -A ohler
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-${total_combinations}

# Define the arrays inside the job script
DATASETS=('mcc_01' 'mcc_05' 'lcmv')
METHODS=('sps' 'random')
SIZES=(100000)  # Using 100k samples as reference cells
REPS=(0)

# Compute task index based on SLURM_ARRAY_TASK_ID
task_id=\$SLURM_ARRAY_TASK_ID

# Calculate indices
rep_index=\$(( (task_id - 1) % \${#REPS[@]} ))
count=\$(( (task_id - 1) / \${#REPS[@]} ))

size_index=\$(( (count) %  \${#SIZES[@]} ))
count=\$(( (count) /  \${#SIZES[@]} ))

method_index=\$(( (count) %  \${#METHODS[@]} ))
count=\$(( (count) /  \${#METHODS[@]} ))

dataset_index=\$(( (count) %  \${#DATASETS[@]} ))
count=\$(( (count) /  \${#DATASETS[@]} ))

# Assign values based on calculated indices
dataset=\${DATASETS[\$dataset_index]}
method=\${METHODS[\$method_index]}
size=\${SIZES[\$size_index]}
rep=\${REPS[\$rep_index]}

# Set reference size based on dataset
if [ "\$dataset" == "mcc_01" ] || [ "\$dataset" == "mcc_05" ]; then
    ref=30
elif [ "\$dataset" == "lcmv" ]; then
    ref=34
else
    echo "Unknown dataset: \$dataset"
    exit 1
fi

# Load environment
source ~/.bashrc
conda activate semitones_env

# Verify SEMITONES is available
if ! python -c "import SEMITONES" 2>/dev/null; then
    echo "ERROR: SEMITONES not available in semitones_env"
    echo "Please install: conda activate semitones_env && pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip"
    exit 1
fi

# Set PROJECT_ROOT environment variable
export PROJECT_ROOT=${PROJECT_ROOT}

# Change directory to the analysis folder
cd \${PROJECT_ROOT}/analysis

# File path to check if job is already done
output_file="\${PROJECT_ROOT}/analysis/results/semitones_downstream/\${dataset}_\${method}_size\${size}_rep\${rep}_results.pkl"

# Check if the file exists and skip if it's already done
if [ -f "\$output_file" ]; then
    echo "File \$output_file exists, skipping task \$task_id..."
    exit 0
fi

# Debugging info
echo "=========================================="
echo "SEMITONES Analysis Task \$task_id"
echo "=========================================="
echo "Dataset: \$dataset"
echo "Reference size: \$ref"
echo "Method: \$method (used as REFERENCE cells)"
echo "Sample size: \$size"
echo "Rep: \$rep"
echo "Node: \$(hostname)"
echo "CPU Info: \$(lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
echo "Conda env: semitones_env"
echo "=========================================="

# Run the Python script with the computed parameters
# This uses sampled cells as REFERENCE cells in SEMITONES
python semitones_downstream_analysis.py \\
    --dataset \$dataset \\
    --ref \$ref \\
    --method \$method \\
    --size \$size \\
    --rep \$rep

exit_code=\$?

if [ \$exit_code -eq 0 ]; then
    echo "Task \$task_id completed successfully"
else
    echo "Task \$task_id failed with exit code \$exit_code"
    exit \$exit_code
fi
EOF

echo ""
echo "Job array submitted successfully!"
echo "Monitor with: squeue -u \$USER"
echo ""
echo "Results will be saved to:"
echo "  \${PROJECT_ROOT}/analysis/results/semitones_downstream/"
echo ""
echo "Logs will be saved to:"
echo "  \${LOG_PATH}/semitones_analysis/logs/"

