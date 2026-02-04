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
REFERENCES=(5 10 20 25 30)
methods=('random')
SIZES=(50000 100000 200000 300000)
REPS=(0 1 2 3 4)

# Calculate total combinations
total_combinations=$(( ${#REFERENCES[@]} * ${#methods[@]} * ${#SIZES[@]} * ${#REPS[@]} ))

# Submit the job array for all combinations
sbatch <<EOF
#!/bin/bash

#SBATCH -J random_array_job_python
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/mcc_01_production/logs/output_random_%A_%a.stdlog
#SBATCH -t 06:00:00
#SBATCH --mem=120G
#SBATCH -A ohler
#SBATCH -p long
#SBATCH -C cascade-lake
#SBATCH --nodelist=maxg11,maxg12,maxg13,maxg14,maxg15,maxg16,maxg17,maxg18,maxg21,maxg22,maxg23,maxg24,maxg25,maxg26
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-${total_combinations}

# Define the arrays inside the job script
REFERENCES=(5 10 20 25 30)
methods=('random')
SIZES=(50000 100000 200000 300000)
REPS=(0 1 2 3 4)

# Initial seed for reproducibility
initial_seed=789

# Compute task index based on SLURM_ARRAY_TASK_ID
task_id=\$SLURM_ARRAY_TASK_ID

# Calculate indices
rep_index=\$(( (task_id - 1) % \${#REPS[@]} ))
count=\$(( (task_id - 1) / \${#REPS[@]} ))

size_index=\$(( (count) %  \${#SIZES[@]} ))
count=\$(( (count) / \${#SIZES[@]} ))

method_index=\$(( (count) %  \${#methods[@]} ))
count=\$(( (count) /  \${#methods[@]} ))

ref_index=\$(( (count) %  \${#REFERENCES[@]} ))
count=\$(( (count) /  \${#REFERENCES[@]} ))

# Assign values based on calculated indices
ref=\${REFERENCES[\$ref_index]}
method=\${methods[\$method_index]}
size=\${SIZES[\$size_index]}
rep=\${REPS[\$rep_index]}

# Compute the seed for this task
seed=\$((initial_seed + task_id))

# Load environment
source ~/.bashrc
conda activate facs_sampling

# Change directory to the working folder
cd ${NOTEBOOK_PATH}/mcc_01

# File path to check if job is already done
file_path="${DATA_PATH}/mcc_01/benchmark/\$ref/\$method/\$size/\$rep/results.pkl"

# Check if the file exists and skip if it's already done
if [ -f "\$file_path" ]; then
    echo "File \$file_path exists, skipping task \$task_id..."
    exit 0
fi

# Debugging info (optional)
echo "Running with Reference=\$ref, Method=\$method, Size=\$size, Rep=\$rep, Seed=\$seed"
echo "Node: \$(hostname)"
echo "CPU Info: \$(lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
echo "CPU Min MHz: \$(lscpu | grep 'CPU min MHz' | cut -d: -f2 | xargs)"
echo "CPU Max MHz: \$(lscpu | grep 'CPU max MHz' | cut -d: -f2 | xargs)"
echo "CPU Current MHz: \$(lscpu | grep '^CPU MHz' | head -1 | cut -d: -f2 | xargs)"

# Run the Python script with the computed parameters
python parallel.py --ref \$ref --method \$method --size \$size --rep \$rep --seed \$seed
EOF
