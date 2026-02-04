#!/bin/bash

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Load cluster configuration
CONFIG_FILE="${SCRIPT_DIR}/../../config/cluster_config.sh"
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    echo "Error: cluster_config.sh not found at: $CONFIG_FILE"
    echo "Please copy from template and configure:"
    echo "cp ${SCRIPT_DIR}/../../config/cluster_config.sh.template ${CONFIG_FILE}"
    exit 1
fi


# Define your combinations
REFERENCES=(5 10 20 25 30)
methods=('atomic')
SIZES=(50000 100000 200000 300000)
REPS=(0 1 2 3 4)

# Calculate total combinations
total_combinations=$(( ${#REFERENCES[@]} * ${#methods[@]} * ${#SIZES[@]} * ${#REPS[@]} ))

# Submit the job array for all combinations
sbatch <<EOF
#!/bin/bash

#SBATCH -J sampling_array_job_R
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/mcc_05_production/logs/output_atom_%A_%a.stdlog
#SBATCH -e ${LOG_PATH}/mcc_05_production/logs/error_atom_%A_%a.errlog
#SBATCH -t 10:00:00
#SBATCH --mem=600G
#SBATCH -A ohler
#SBATCH -p long
#SBATCH -C cascade-lake
#SBATCH --nodelist=maxg11,maxg12,maxg13,maxg14,maxg15,maxg16,maxg17,maxg18,maxg21,maxg22,maxg23,maxg24,maxg25,maxg26
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=localtmp:50G
#SBATCH --array=1-${total_combinations}

# Define the arrays inside the job script
REFERENCES=(5 10 20 25 30)
methods=('atomic')
SIZES=(50000 100000 200000 300000)
REPS=(0 1 2 3 4)

# Initial seed for reproducibility
initial_seed=3956

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

# Load the Guix environment
GUIX_PROFILE=\${HOME}/.guix-profile
source \${GUIX_PROFILE}/etc/profile
# Suppress error if file doesn't exist (non-critical)
source /etc/profile.d/ge_interactive-environment.sh 2>/dev/null || true



# Change directory to the working folder
cd ${NOTEBOOK_PATH}/mcc_05

# File path to check if job is already done
file_path="${DATA_PATH}/mcc_05/benchmark/\${ref}/\${method}/\${size}/\${rep}/results.csv"

# Error handler function (define before traps)
error_handler() {
    exit_code=\$?
    echo "=========================================" >&2
    echo "ERROR: Job failed with exit code \$exit_code" >&2
    echo "Task ID: \$task_id" >&2
    echo "Reference: \$ref, Method: \$method, Size: \$size, Rep: \$rep" >&2
    echo "Node: \$(hostname)" >&2
    echo "Time: \$(date)" >&2
    echo "=========================================" >&2
    # Also print to stdout for visibility
    echo "ERROR: Job failed with exit code \$exit_code"
    echo "Task ID: \$task_id"
    echo "Reference: \$ref, Method: \$method, Size: \$size, Rep: \$rep"
    exit \$exit_code
}

# Cleanup handler for signals
cleanup_handler() {
    exit_code=\$?
    # Capture final memory state
    final_mem=\$(free -h | grep '^Mem:')
    final_mem_used=\$(free -g | grep '^Mem:' | awk '{print \$2}')
    
    # Get actual SLURM memory limit from environment variables
    # SLURM_MEM_PER_NODE is in MB, convert to GB for display
    if [ -n "\$SLURM_MEM_PER_NODE" ]; then
        slurm_mem_mb=\$SLURM_MEM_PER_NODE
        slurm_mem_gb=\$((slurm_mem_mb / 1024))
        slurm_mem_limit="\$slurm_mem_gb G (from SLURM_MEM_PER_NODE=\$slurm_mem_mb MB)"
    elif [ -n "\$SLURM_MEM_PER_CPU" ] && [ -n "\$SLURM_CPUS_PER_TASK" ]; then
        slurm_mem_mb=\$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK))
        slurm_mem_gb=\$((slurm_mem_mb / 1024))
        slurm_mem_limit="\$slurm_mem_gb G (calculated from SLURM_MEM_PER_CPU=\$SLURM_MEM_PER_CPU MB × \$SLURM_CPUS_PER_TASK CPUs)"
    else
        slurm_mem_limit="Unknown (SLURM memory env vars not set)"
    fi
    
    if [ \$exit_code -ne 0 ] && [ \$exit_code -ne 130 ] && [ \$exit_code -ne 143 ]; then
        echo "=========================================" >&2
        echo "WARNING: Job terminated (exit code: \$exit_code)" >&2
        echo "Task ID: \$task_id" >&2
        echo "Reference: \$ref, Method: \$method, Size: \$size, Rep: \$rep" >&2
        echo "Node: \$(hostname)" >&2
        echo "Time: \$(date)" >&2
        echo "Final Memory State: \$final_mem" >&2
        echo "Memory Used (GB): \$final_mem_used" >&2
        echo "SLURM Memory Limit: \$slurm_mem_limit" >&2
        echo "" >&2
        echo "This may indicate:" >&2
        if [ -n "\$slurm_mem_gb" ]; then
            echo "  - Out of memory (OOM) - check if memory used approaches \$slurm_mem_gb G" >&2
        else
            echo "  - Out of memory (OOM) - check memory usage" >&2
        fi
        echo "  - Time limit exceeded" >&2
        echo "  - Manual cancellation" >&2
        echo "  - System signal (SIGTERM, SIGKILL)" >&2
        echo "" >&2
        echo "To check job cancellation reason, run:" >&2
        echo "  sacct -j \$SLURM_JOB_ID --format=JobID,State,ExitCode,MaxRSS,ReqMem" >&2
        echo "=========================================" >&2
        # Also print to stdout
        echo "WARNING: Job terminated (exit code: \$exit_code)"
        echo "Task ID: \$task_id"
        echo "Reference: \$ref, Method: \$method, Size: \$size, Rep: \$rep"
        echo "Final Memory: \$final_mem"
        echo "SLURM Memory Limit: \$slurm_mem_limit"
    fi
}

# Set up error handling
trap 'error_handler' ERR
trap 'cleanup_handler' EXIT TERM INT

# Check if the file exists and skip if it's already done
if [ -f "\$file_path" ]; then
    echo "File \$file_path exists, skipping task \$task_id..."
    exit 0
fi

# Enable error trapping (after file check to avoid false positives)
set -e

# Debugging info (optional)
echo "Running with Reference=\$ref, Method=\$method, Size=\$size, Rep=\$rep, Seed=\$seed"
echo "Node: \$(hostname)"
echo "CPU Info: \$(lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
echo "CPU Min MHz: \$(lscpu | grep 'CPU min MHz' | cut -d: -f2 | xargs)"
echo "CPU Max MHz: \$(lscpu | grep 'CPU max MHz' | cut -d: -f2 | xargs)"
echo "CPU Current MHz: \$(lscpu | grep '^CPU MHz' | head -1 | cut -d: -f2 | xargs)"
mem_info=\$(free -h | grep '^Mem:')
echo "Memory Info: \$mem_info"

# Display SLURM memory allocation
if [ -n "\$SLURM_MEM_PER_NODE" ]; then
    slurm_mem_mb=\$SLURM_MEM_PER_NODE
    slurm_mem_gb=\$((slurm_mem_mb / 1024))
    echo "SLURM Allocated Memory: \$slurm_mem_gb G (\$slurm_mem_mb MB from SLURM_MEM_PER_NODE)"
elif [ -n "\$SLURM_MEM_PER_CPU" ] && [ -n "\$SLURM_CPUS_PER_TASK" ]; then
    slurm_mem_mb=\$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK))
    slurm_mem_gb=\$((slurm_mem_mb / 1024))
    echo "SLURM Allocated Memory: \$slurm_mem_gb G (\$slurm_mem_mb MB = \$SLURM_MEM_PER_CPU MB × \$SLURM_CPUS_PER_TASK CPUs)"
else
    echo "SLURM Memory Environment Variables:"
    echo "  SLURM_MEM_PER_NODE: \${SLURM_MEM_PER_NODE:-not set}"
    echo "  SLURM_MEM_PER_CPU: \${SLURM_MEM_PER_CPU:-not set}"
    echo "  SLURM_CPUS_PER_TASK: \${SLURM_CPUS_PER_TASK:-not set}"
fi

# Run the R script with the computed parameters
# Errors will be captured in both stderr log and stdout log
echo "Starting R script execution at \$(date)"
echo "Command: ~/.guix-profile/bin/Rscript atomic_sketching_mcc.R \$ref \$method \$size \$rep \$seed"

# Monitor memory usage in background (log every 30 seconds)
(
    while true; do
        mem_usage=\$(free -h | grep '^Mem:' | awk '{print "Used: " \$3 "/" \$2 ", Available: " \$7}')
        echo "[Memory Monitor] \$(date): \$mem_usage" >&2
        sleep 30
    done
) &
monitor_pid=\$!

# Run R script and capture exit code
# Temporarily disable exit-on-error to capture exit code manually
set +e
~/.guix-profile/bin/Rscript atomic_sketching_mcc.R \$ref \$method \$size \$rep \$seed 2>&1
rscript_exit_code=\$?
set -e  # Re-enable exit-on-error

# Stop memory monitor
kill \$monitor_pid 2>/dev/null || true

if [ \$rscript_exit_code -ne 0 ]; then
    echo "=========================================" >&2
    echo "ERROR: R script exited with code \$rscript_exit_code" >&2
    echo "Task ID: \$task_id" >&2
    echo "Reference: \$ref, Method: \$method, Size: \$size, Rep: \$rep" >&2
    echo "Node: \$(hostname)" >&2
    echo "Time: \$(date)" >&2
    echo "=========================================" >&2
    # Also print to stdout
    echo "ERROR: R script exited with code \$rscript_exit_code"
    echo "Task ID: \$task_id"
    echo "Reference: \$ref, Method: \$method, Size: \$size, Rep: \$rep"
    exit \$rscript_exit_code
fi

echo "R script completed successfully at \$(date)"
EOF
