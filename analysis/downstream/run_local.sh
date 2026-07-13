#!/bin/bash
# Local driver to regenerate the full Leiden downstream block on a single big node.
# Runs all modes x conditions with a bounded concurrency to stay within RAM.
cd "$(dirname "$0")"
PY=/fast/AG_Ohler/ekarimi/miniforge/envs/facs_sampling/bin/python
export PROJECT_ROOT="${PROJECT_ROOT:-/fast/AG_Ohler/ekarimi/projects/sparsesampler_reproducibility}"
mkdir -p results logs

# Cap BLAS threads per job so concurrent jobs don't oversubscribe the 112 cores.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-12}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export NUMBA_NUM_THREADS="$OMP_NUM_THREADS"

MAXJOBS="${MAXJOBS:-6}"
MODES=("${MODES:-leiden}")
DATASETS=(lcmv mcc)
METHODS=(sps random)
SIZES=(50000 100000 200000)
REPS=(0 1 2)

run_one() {
    local mode=$1 ds=$2 method=$3 size=$4 rep=$5
    local tag="${mode}_${ds}_${method}_${size}_rep${rep}"
    local log="logs/${tag}.log"
    local t0=$(date +%s)
    "$PY" leiden_discovery.py --dataset "$ds" --method "$method" --size "$size" --rep "$rep" > "$log" 2>&1
    local rc=$?
    echo "[$(date +%H:%M:%S)] DONE $tag rc=$rc $(( $(date +%s) - t0 ))s"
}
export -f run_one
export PY PROJECT_ROOT

JOBS=()
for mode in $MODES; do
    for ds in "${DATASETS[@]}"; do
        for method in "${METHODS[@]}"; do
            for size in "${SIZES[@]}"; do
                for rep in "${REPS[@]}"; do
                    JOBS+=("$mode $ds $method $size $rep")
                done
            done
        done
    done
done

echo "Total jobs: ${#JOBS[@]} | MAXJOBS=$MAXJOBS | threads/job=$OMP_NUM_THREADS"
printf "%s\n" "${JOBS[@]}" | xargs -P "$MAXJOBS" -L1 bash -c 'run_one "$@"' _
echo "ALL_DONE"
