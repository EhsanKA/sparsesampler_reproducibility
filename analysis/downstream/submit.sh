#!/bin/bash
cd "$(dirname "$0")"
JOB=$1
WRAPPER="bash worker.sh $JOB"
PROJECT_ROOT="${PROJECT_ROOT:-/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility}"

case "$JOB" in
  bins|quality|proportions|adjacent)
    sbatch --job-name="ds_${JOB}" \
      --output="results/${JOB}_%j.out" --error="results/${JOB}_%j.err" \
      --time=06:00:00 --mem=200G --cpus-per-task=1 --partition=normal -A ohler \
      --wrap="$WRAPPER"
    ;;
  leiden|within_pop)
    sbatch --job-name="ds_${JOB}" \
      --output="results/${JOB}_%A_%a.out" --error="results/${JOB}_%A_%a.err" \
      --time=03:00:00 --mem=200G --cpus-per-task=8 --partition=normal -A ohler \
      --array=0-35 --wrap="$WRAPPER"
    ;;
  within_pop_precompute)
    sbatch --job-name="ds_wpop_pre" \
      --output="results/wpop_pre_%a.out" --error="results/wpop_pre_%a.err" \
      --time=06:00:00 --mem=200G --cpus-per-task=4 --partition=normal -A ohler \
      --array=0-1 --wrap="$WRAPPER"
    ;;
  *)
    echo "Usage: bash submit.sh {leiden|bins|quality|proportions|adjacent|within_pop_precompute|within_pop}"
    exit 1
    ;;
esac
