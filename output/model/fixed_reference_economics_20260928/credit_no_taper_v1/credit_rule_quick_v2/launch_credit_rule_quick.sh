#!/usr/bin/env bash
# Prepared only.  Do not run after the stated absolute deadline.
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:10:00
#SBATCH --array=0-1
set -euo pipefail
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
remote=/scratch/td2248/projects/fixed_reference_credit_rule_quick_20260929
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/credit_no_taper_v1/credit_rule_quick_v1
source=$remote/source_credit_rule_quick_v2
case=(ours author); selected=${case[$SLURM_ARRAY_TASK_ID]}
[[ $(date -u +%s) -lt $(date -u -d '2026-09-30T01:00:25Z' +%s) ]] || { echo 'absolute deadline expired' >&2; exit 2; }
[[ -f "$source/run_credit_rule_quick.py" && -f "$source/strict_tenure.py" && -f "$source/plan.json" ]] || { echo 'stage pinned packet first' >&2; exit 2; }
mkdir -p "$remote/results_v2" "$remote/cache"
[[ ! -e "$remote/results_v2/$selected" ]] || { echo 'refusing repeat' >&2; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 NUMBA_CACHE_DIR=/work/numba_cache MPLBACKEND=Agg
apptainer exec --bind "$reference:$original:ro,$source:/work/credit_rule_source:ro,$remote/results_v2:/work/results,$remote/cache:/work/numba_cache" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/credit_rule_source/run_credit_rule_quick.py --case "$selected" --output "/work/results/$selected"
