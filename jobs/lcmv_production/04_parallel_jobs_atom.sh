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
REFERENCES=(1 5 10 20 34)
methods=('atomic')
SIZES=(50000 100000 200000)
REPS=(0 1 2 3 4)

# Calculate total combinations
total_combinations=$(( ${#REFERENCES[@]} * ${#methods[@]} * ${#SIZES[@]} * ${#REPS[@]} ))

# Submit the job array for all combinations
qsub -t 1-$total_combinations <<EOF
#!/bin/bash

#$ -N sampling_array_job_R
#$ -m beas
#$ -M ${CLUSTER_EMAIL}
#$ -o ${LOG_PATH}/lcmv_production/logs/output_atom.stdlog
#$ -j y
#$ -l h_rt=10:00:00
#$ -l m_mem_free=100G

# Define the arrays inside the job script
REFERENCES=(1 5 10 20 34)
methods=('atomic')
SIZES=(50000 100000 200000)
REPS=(0 1 2 3 4)

# Initial seed for reproducibility
initial_seed=3956

# Compute task index based on SGE_TASK_ID
task_id=\$SGE_TASK_ID

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

# Load the Guix environment
GUIX_PROFILE=\${HOME}/.guix-profile
source \${GUIX_PROFILE}/etc/profile
source /etc/profile.d/ge_interactive-environment.sh



# Change directory to the working folder
cd ${NOTEBOOK_PATH}/lcmv

# File path to check if job is already done
file_path="${DATA_PATH}/lcmv/benchmark/\${ref}/\${method}/\${size}/\${rep}/results.csv"

# Check if the file exists and skip if it's already done
if [ -f "\$file_path" ]; then
    echo "File \$file_path exists, skipping task \$task_id..."
    exit 0
fi

# Debugging info (optional)
echo "Running with Reference=\$ref, Method=\$method, Size=\$size, Rep=\$rep, Seed=\$seed"

# Run the R script with the computed parameters
~/.guix-profile/bin/Rscript atomic_sketching_lcmv.R \$ref \$method \$size \$rep \$seed
EOF
