#!/bin/bash
# Script to run downstream analysis comparison for reviewer response

# Set project root if not already set
if [ -z "$PROJECT_ROOT" ]; then
    export PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
fi

# Change to analysis directory
cd "$PROJECT_ROOT/analysis"

# Activate conda environment if needed (adjust as necessary)
# conda activate facs_sampling

echo "Starting downstream analysis comparison..."
echo "This may take a while as it processes multiple datasets and configurations."

# Run the analysis script
python downstream_analysis_comparison.py

echo "Analysis complete. Generating figures..."

# Generate figures
python plot_downstream_analysis.py

echo "All done! Check the results in:"
echo "  - Results: $PROJECT_ROOT/analysis/results/downstream_analysis/"
echo "  - Figures: $PROJECT_ROOT/figures/revision/downstream_analysis/"







