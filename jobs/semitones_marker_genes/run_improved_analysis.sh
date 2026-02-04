#!/bin/bash

# Improved SEMITONES Marker Gene Analysis Job Script
# Tests different SEMITONES distance metrics to improve SPS performance

# Load cluster configuration
if [ -f "../../config/cluster_config.sh" ]; then
    source ../../config/cluster_config.sh
else
    echo "Error: cluster_config.sh not found. Please copy from template and configure."
    echo "cp ../../config/cluster_config.sh.template ../../config/cluster_config.sh"
    exit 1
fi

echo "========================================"
echo "Improved SEMITONES Marker Gene Analysis"
echo "========================================"
echo "Dataset: lcmv (ref=34)"
echo "Sample size: 100K"
echo "Methods: SPS, Random"
echo "Marker methods: logfc, semitones_cosine, semitones_euclidean, semitones_correlation"
echo "Environment: semitones_env"
echo "========================================"

# Submit the job
sbatch <<EOF
#!/bin/bash

#SBATCH -J improved_marker_analysis
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${CLUSTER_EMAIL}
#SBATCH -o ${LOG_PATH}/semitones_marker_genes/logs/output_improved_%j.stdlog
#SBATCH -e ${LOG_PATH}/semitones_marker_genes/logs/output_improved_%j.stderr
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -A ohler
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Load environment - use semitones_env which has SEMITONES installed
source ~/.bashrc
conda activate semitones_env

# Verify SEMITONES is available
if ! python -c "import SEMITONES" 2>/dev/null; then
    echo "ERROR: SEMITONES not available in semitones_env"
    echo "Please install: pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip"
    exit 1
fi

# Set PROJECT_ROOT environment variable
export PROJECT_ROOT=${PROJECT_ROOT}

# Change directory to the analysis folder
cd \${PROJECT_ROOT}/analysis/semitones_marker_genes

# Debugging info
echo "=========================================="
echo "Improved SEMITONES Marker Gene Analysis"
echo "=========================================="
echo "Node: \$(hostname)"
echo "CPU Info: \$(lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
echo "Memory: \$(free -h | grep Mem | awk '{print \$2}')"
echo "Conda env: semitones_env"
echo "Working directory: \$(pwd)"
echo "=========================================="

# Run the Python script
python improved_semitones_analysis.py

exit_code=\$?

if [ \$exit_code -eq 0 ]; then
    echo "Improved analysis completed successfully"
else
    echo "Improved analysis failed with exit code \$exit_code"
    exit \$exit_code
fi
EOF

echo ""
echo "Job submitted successfully!"
echo "Monitor with: squeue -u \$USER"
echo ""
echo "Results will be saved to:"
echo "  ${PROJECT_ROOT}/analysis/semitones_marker_genes/results/improved_semitones_analysis_results.csv"
echo ""
echo "Logs will be saved to:"
echo "  ${LOG_PATH}/semitones_marker_genes/logs/"