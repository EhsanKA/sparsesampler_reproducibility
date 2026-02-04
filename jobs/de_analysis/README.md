# DE Analysis Parallel Jobs

This directory contains SLURM job scripts for running DE gene analysis in parallel.

## Files

1. **`01_parallel_jobs_de_analysis.sh`** - Main SLURM job array script
   - Submits parallel jobs for all combinations of:
     - Datasets: mcc, lcmv
     - Methods: sps, random
     - Sample sizes: 100000, 200000
     - Replicates: 0
   - Total: 2 × 2 × 2 × 1 = 8 jobs

2. **`02_combine_results.sh`** - Script to combine results after all jobs complete
   - Combines all detailed result files into a summary CSV
   - Run this after all SLURM jobs finish

## Usage

### Step 1: Submit parallel jobs
```bash
cd jobs/de_analysis
bash 01_parallel_jobs_de_analysis.sh
```

This will submit a SLURM job array with 8 jobs (one for each combination).

### Step 2: Monitor jobs
```bash
squeue -u $USER
```

### Step 3: After all jobs complete, combine results
```bash
bash 02_combine_results.sh
```

### Step 4: Generate visualizations
```bash
cd ../../analysis
python plot_de_analysis.py
```

## Output

- **Detailed results**: `analysis/results/de_analysis/{dataset}_{method}_size{size}_rep{rep}_de_detailed.csv`
- **Summary results**: `analysis/results/de_analysis/de_analysis_summary.csv` (created after combining)
- **Logs**: `jobs/de_analysis/logs/output_de_{job_id}_{task_id}.stdlog`

## Job Configuration

- **Time limit**: 12 hours
- **Memory**: 200GB
- **CPUs**: 4 per task
- **Partition**: long
- **Nodes**: cascade-lake nodes (maxg11-maxg26)

## Notes

- Each job analyzes ALL cell types (not limited to 10)
- Each job compares top 50 DE genes between full and subsampled data
- Jobs will skip if output file already exists (idempotent)






