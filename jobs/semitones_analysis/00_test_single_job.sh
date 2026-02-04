#!/bin/bash

# Test script for running a single SEMITONES analysis job
# Useful for testing before submitting the full job array

# Load cluster configuration
if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found. Please copy from template and configure."
    echo "cp ../../config/cluster_config.sh.template ../../config/cluster_config.sh"
    exit 1
fi

# Test configuration (single job)
# Using 100k samples for SEMITONES reference cells
DATASET="mcc"
REF=30
METHOD="sps"
SIZE=100000
REP=0

echo "Submitting single SEMITONES test job"
echo "Dataset: $DATASET"
echo "Reference: $REF"
echo "Method: $METHOD"
echo "Size: $SIZE"
echo "Rep: $REP"
echo ""

# Submit the job
sbatch <<EOF
#!/bin/bash

#SBATCH -J semitones_test
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/semitones_analysis/logs/test_semitones_%j.stdlog
#SBATCH -e ${LOG_PATH}/semitones_analysis/logs/test_semitones_%j.stderr
#SBATCH -t 24:00:00
#SBATCH --mem=300G
#SBATCH -A ohler
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

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

# Debugging info
echo "=========================================="
echo "SEMITONES Test Job"
echo "=========================================="
echo "Dataset: ${DATASET}"
echo "Reference size: ${REF}"
echo "Method: ${METHOD} (used as REFERENCE cells)"
echo "Sample size: ${SIZE}"
echo "Rep: ${REP}"
echo "Node: \$(hostname)"
echo "CPU Info: \$(lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
echo "Conda env: semitones_env"
echo "=========================================="

# Run the Python script
python semitones_downstream_analysis.py \\
    --dataset ${DATASET} \\
    --ref ${REF} \\
    --method ${METHOD} \\
    --size ${SIZE} \\
    --rep ${REP}

exit_code=\$?

if [ \$exit_code -eq 0 ]; then
    echo "Test job completed successfully"
else
    echo "Test job failed with exit code \$exit_code"
    exit \$exit_code
fi
EOF

echo ""
echo "Test job submitted successfully!"
echo "Monitor with: squeue -u \$USER"
echo ""
echo "Results will be saved to:"
echo "  \${PROJECT_ROOT}/analysis/results/semitones_downstream/${DATASET}_${METHOD}_size${SIZE}_rep${REP}_*"

