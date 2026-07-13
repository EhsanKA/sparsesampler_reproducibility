#!/bin/bash
#SBATCH --job-name=gen_figs
#SBATCH --output=logs/gen_figs_%A_%a.out
#SBATCH --error=logs/gen_figs_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=normal
#SBATCH --array=1-5

# Job array for figure generation
# Task 1: Main coverage figures (MCC + LCMV with 1% rarity)
# Task 2: Main time performance figures (MCC + LCMV)
# Task 3: Revision coverage figures (MCC_01 with 0.1% + MCC_05 with 0.5%)
# Task 4: Revision time performance figures (MCC_01 + MCC_05)
# Task 5: Main UMAP figures (MCC + LCMV)

set -e

# Set project paths
PROJECT_ROOT="/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility"
FIGURES_DIR="$PROJECT_ROOT/figures"

echo "Job array task ID: $SLURM_ARRAY_TASK_ID"
echo "Project root: $PROJECT_ROOT"
export PROJECT_ROOT

# Create logs directory
mkdir -p "$PROJECT_ROOT/logs"

# Activate conda environment
source /fast/AG_Ohler/ekarimi/miniforge/etc/profile.d/conda.sh
conda activate sps
echo "Using Python: $(which python)"

echo "============================================"
echo "Starting task $SLURM_ARRAY_TASK_ID at $(date)"
echo "============================================"

case $SLURM_ARRAY_TASK_ID in
    1)
        echo "--- Task 1: Main coverage figures (MCC + LCMV, 1% rarity) ---"
        echo "Rare cell definition: distance > 75th percentile AND frequency < 1%"
        cd "$FIGURES_DIR/main/coverage"
        python coverage_combined.py
        echo "Generated:"
        echo "  - figures/figure1_combined_coverage.jpg"
        echo "  - figures/supp_figure3_mcc_coverage.jpg"
        echo "  - figures/supp_figure4_lcmv_coverage.jpg"
        ;;
    2)
        echo "--- Task 2: Main time performance figures (MCC + LCMV) ---"
        cd "$FIGURES_DIR/main/time_performance"
        python combined_time_plot.py
        echo "Generated:"
        echo "  - figures/figure1_combined_time.jpg"
        echo "  - figures/supp_figure5-6_mcc_time.jpg"
        echo "  - figures/supp_figure5-6_lcmv_time.jpg"
        ;;
    3)
        echo "--- Task 3: Revision coverage figures (MCC_01 + MCC_05) ---"
        echo "Rare cell definition: distance > 75th percentile AND frequency < threshold"
        echo "  - MCC_01: frequency < 0.1%"
        echo "  - MCC_05: frequency < 0.5%"
        cd "$FIGURES_DIR/revision/coverage"
        python coverage_combined.py
        echo "Generated:"
        echo "  - figures/revision/figure1_combined_coverage.jpg"
        echo "  - figures/revision/supp_figure3_mcc_01_coverage.jpg"
        echo "  - figures/revision/supp_figure3_mcc_05_coverage.jpg"
        ;;
    4)
        echo "--- Task 4: Revision time performance figures (MCC_01 + MCC_05) ---"
        cd "$FIGURES_DIR/revision/time_performance"
        python combined_time_plot.py
        echo "Generated:"
        echo "  - figures/revision/figure1_combined_time.jpg"
        echo "  - figures/revision/supp_figure5-6_mcc_01_time.jpg"
        echo "  - figures/revision/supp_figure5-6_mcc_05_time.jpg"
        ;;
    5)
        echo "--- Task 5: Main UMAP figures (MCC + LCMV) ---"
        cd "$FIGURES_DIR/main/umaps"
        python combined_umaps.py
        echo "Generated:"
        echo "  - figures/main/umaps/ UMAP figures"
        ;;
    *)
        echo "Unknown task ID: $SLURM_ARRAY_TASK_ID"
        exit 1
        ;;
esac

echo ""
echo "============================================"
echo "Task $SLURM_ARRAY_TASK_ID completed at $(date)"
echo "============================================"
