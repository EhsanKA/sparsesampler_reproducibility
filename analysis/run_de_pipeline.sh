#!/bin/bash
# Script to run DE analysis and visualization pipeline

# Set project root if not already set
if [ -z "$PROJECT_ROOT" ]; then
    export PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
fi

# Change to analysis directory
cd "$PROJECT_ROOT/analysis"

echo "=========================================="
echo "Step 1: Running DE Analysis"
echo "=========================================="
echo "This will analyze DE genes for:"
echo "  - MCC full dataset (ref 30)"
echo "  - LCMV full dataset (ref 34)"
echo "  - Comparing with subsampled data (100K, 200K cells)"
echo "  - Methods: sps, random"
echo ""
echo "This may take a while (DE analysis is computationally intensive)..."
echo ""

python run_de_analysis.py

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Step 2: Generating Visualizations"
    echo "=========================================="
    echo ""
    
    python plot_de_analysis.py
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "=========================================="
        echo "Pipeline completed successfully!"
        echo "=========================================="
        echo ""
        echo "Results:"
        echo "  - DE analysis results: $PROJECT_ROOT/analysis/results/de_analysis/"
        echo "  - Figures: $PROJECT_ROOT/figures/revision/de_analysis/"
        echo ""
    else
        echo "Error generating visualizations"
        exit 1
    fi
else
    echo "Error running DE analysis"
    exit 1
fi






