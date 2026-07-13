#!/bin/bash
set -e
cd "$(dirname "$0")"
mkdir -p results

source ~/.bashrc
conda activate facs_sampling
export PROJECT_ROOT="${PROJECT_ROOT:-/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility}"

JOB=$1
DATASETS=(lcmv mcc)
METHODS=(sps random)
SIZES=(50000 100000 200000)
REPS=(0 1 2)

decode_array() {
    local idx=$1
    local n=$(( ${#METHODS[@]} * ${#SIZES[@]} * ${#REPS[@]} ))
    local dataset_idx=$((idx / n))
    local r=$((idx % n))
    local method_idx=$((r / (${#SIZES[@]} * ${#REPS[@]})))
    r=$((r % (${#SIZES[@]} * ${#REPS[@]})))
    local size_idx=$((r / ${#REPS[@]}))
    local rep_idx=$((r % ${#REPS[@]}))
    DS="${DATASETS[$dataset_idx]}"
    METHOD="${METHODS[$method_idx]}"
    SIZE="${SIZES[$size_idx]}"
    REP="${REPS[$rep_idx]}"
}

case "$JOB" in
  leiden)
    decode_array "${SLURM_ARRAY_TASK_ID:-0}"
    python leiden_discovery.py --dataset "$DS" --method "$METHOD" --size "$SIZE" --rep "$REP"
    ;;
  within_pop)
    decode_array "${SLURM_ARRAY_TASK_ID:-0}"
    python within_pop.py --dataset "$DS" --method "$METHOD" --size "$SIZE" --rep "$REP"
    ;;
  bins)
    python bins.py --task characterize --dataset all
    ;;
  quality)
    python bins.py --task quality --dataset all
    ;;
  proportions)
    python proportions.py --dataset all
    ;;
  adjacent)
    python adjacent_bins.py --dataset all
    ;;
  within_pop_precompute)
    python within_pop.py --precompute --dataset "${DATASETS[$SLURM_ARRAY_TASK_ID]}"
    ;;
  *)
    echo "Unknown job: $JOB"
    exit 1
    ;;
esac
