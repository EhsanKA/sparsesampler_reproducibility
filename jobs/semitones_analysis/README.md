# SEMITONES Analysis Jobs

This directory contains SLURM job scripts for running SEMITONES analysis using sampled cells as references.

## Overview

SEMITONES (Similarity Enrichment for Multiple Interrogation Techniques in Omics via Neighborhoods of Exemplars) is used to evaluate how well sampled cells represent the full dataset. The sampled cells are used as **reference cells**, and all other cells are **query cells** that get compared to these references.

## Files

- `01_parallel_jobs_semitones.sh` - Main SLURM job array script for running SEMITONES analysis

## Usage

### Submit Jobs

```bash
cd jobs/semitones_analysis
bash 01_parallel_jobs_semitones.sh
```

### What It Does

The script runs SEMITONES analysis for:
- **Datasets**: mcc_01, mcc_05, lcmv
- **Methods**: sps, random (comparing SPS vs Random as references)
- **Sample sizes**: 100000 (100k samples used as reference cells)
- **Replications**: 0

For each combination, it:
1. Loads the full dataset
2. Loads the sampled cell indices (these become REFERENCE cells)
3. Uses SEMITONES to calculate enrichment scores for all query cells
4. Analyzes enrichment by cell type
5. Creates visualizations
6. Saves results to `analysis/results/semitones_downstream/`

### Job Configuration

- **Time limit**: 24 hours
- **Memory**: 150GB
- **CPUs**: 8
- **Conda environment**: `semitones_env`
- **Account**: ohler
- **Partition**: normal

## Output

### Results Files
Located in: `analysis/results/semitones_downstream/`

- `{dataset}_{method}_size{size}_rep{rep}_results.pkl` - Full results (pickle)
- `{dataset}_{method}_size{size}_rep{rep}_summary.csv` - Summary by cell type (CSV)
- `{dataset}_{method}_size{size}_enrichment_analysis.png` - Visualizations

### Log Files
Located in: `jobs/semitones_analysis/logs/`

- `output_semitones_{job_id}_{array_id}.stdlog` - Standard output
- `output_semitones_{job_id}_{array_id}.stderr` - Standard error

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# Check specific job
squeue -j <job_id>

# View log file
tail -f jobs/semitones_analysis/logs/output_semitones_<job_id>_<array_id>.stdlog
```

## Understanding the Results

### Enrichment Scores
- **High scores**: Query cell is well-represented by reference cells ✅
- **Low scores**: Query cell is poorly represented (might be rare/novel type) ⚠️

### Cell Type Coverage
- Shows which cell types are well/poorly represented in the reference set
- SPS samples should have better coverage than random samples

### Key Metrics
- **Mean enrichment**: Overall quality of reference set
- **Per-cell-type enrichment**: Which types are well represented
- **Poorly represented cells**: Cells with very low enrichment

## Comparison: SPS vs Random

The script runs both SPS and Random sampling as references. Expected results:
- **SPS samples** should have:
  - Higher mean enrichment scores
  - Better coverage of rare cell types
  - More consistent enrichment across cell types

- **Random samples** should have:
  - Lower mean enrichment scores
  - Poorer coverage of rare cell types
  - More variable enrichment

## Notes

- The script uses `semitones_env` conda environment (not `facs_sampling`)
- SEMITONES requires Python < 3.10
- Jobs will skip if output files already exist (idempotent)
- Large datasets may require more memory/time

## Troubleshooting

### SEMITONES Not Found
If you see "SEMITONES not available":
```bash
conda activate semitones_env
pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip
```

### Memory Issues
If jobs fail due to memory:
- Increase `--mem` in the SLURM script
- Process smaller datasets first
- Use smaller sample sizes

### Timeout Issues
If jobs timeout:
- Increase `-t` time limit in SLURM script
- Process smaller datasets first

## Related Files

- `analysis/semitones_downstream_analysis.py` - Main analysis script
- `analysis/SEMITONES_EXPLANATION.md` - Detailed explanation of SEMITONES
- `analysis/SEMITONES_SUMMARY.md` - Quick reference guide

