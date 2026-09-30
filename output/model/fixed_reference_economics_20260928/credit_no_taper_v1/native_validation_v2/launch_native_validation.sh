#!/usr/bin/env bash
#SBATCH --job-name=e5f_no_taper_native_validation
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:20:00
set -euo pipefail

# Prepared only: lead must pin plan.json fields and submit.  The mock controller
# tests run first and perform zero lifecycle solves; no retry is present.
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
remote=/scratch/td2248/projects/fixed_reference_credit_no_taper_validation_20260929
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/credit_no_taper_v1/native_validation_v2
source=$remote/source_native_validation_v2
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'submit with sbatch after lead review' >&2; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export NATIVE_VALIDATION_STARTED_EPOCH="$(date +%s)"
export NATIVE_VALIDATION_DEADLINE_EPOCH="$((NATIVE_VALIDATION_STARTED_EPOCH + 1200))"
export NUMBA_CACHE_DIR="/work/numba_cache"
mkdir -p "$remote/results" "$NUMBA_CACHE_DIR"
[[ -f "$source/plan.json" && -f "$source/run_native_validation.py" ]] || { echo 'immutable validation source was not staged' >&2; exit 2; }
[[ ! -e "$remote/results/native_validation_v2" ]] || { echo 'refusing duplicate validation output' >&2; exit 2; }
[[ ! -e "$remote/results/controller_mock_tests_v2" ]] || { echo 'refusing duplicate mock-test output' >&2; exit 2; }
cd "$source"
python - <<'PY'
import hashlib, json, pathlib
p=pathlib.Path('plan.json'); d=pathlib.Path('run_native_validation.py'); plan=json.loads(p.read_text())
if plan['driver_sha256'].startswith('TO_BE_'): raise SystemExit('lead must pin driver and file hashes before submission')
if plan['driver_sha256'] != hashlib.sha256(d.read_bytes()).hexdigest(): raise SystemExit('driver pin mismatch')
if plan['schema'] != 'block0506_renter_no_taper_native_validation_v2': raise SystemExit('wrong plan schema')
if plan['total_seconds'] != 1200 or plan['case_seconds'] != 360: raise SystemExit('wrong time budget')
PY
python=(apptainer exec --bind "$reference:$original:ro,$source:/work/no_taper_source:ro,$remote/results:/work/results,$remote/cache:/work/numba_cache" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python)
"${python[@]}" /work/no_taper_source/run_native_validation.py --plan /work/no_taper_source/plan.json --output /work/results/controller_mock_tests_v2 --self-test
"${python[@]}" /work/no_taper_source/run_native_validation.py --plan /work/no_taper_source/plan.json --output /work/results/native_validation_v2
