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
export NATIVE_VALIDATION_STARTED_EPOCH="$(date +%s)"
export NATIVE_VALIDATION_DEADLINE_EPOCH="$((NATIVE_VALIDATION_STARTED_EPOCH + 1200))"

# Prepared only: lead must pin plan.json fields and submit.  The mock controller
# tests run first and perform zero lifecycle solves; no retry is present.
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
remote=/scratch/td2248/projects/fixed_reference_credit_no_taper_validation_20260929
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/credit_no_taper_v1/native_validation_v3
source=$remote/source_native_validation_v3
mode=${1:-}
preflight_sha=${2:-}
if [[ -n "$mode" ]]; then
  if [[ "$mode" == "--self-test-only" ]]; then
    [[ -z "$preflight_sha" ]] || { echo 'self-test-only takes no receipt hash' >&2; exit 2; }
  elif [[ "$mode" == "--preflight-receipt-sha" && -n "$preflight_sha" ]]; then
    mode=""
  else
    echo 'usage: sbatch launch_native_validation.sh [--self-test-only | --preflight-receipt-sha SHA256]' >&2; exit 2
  fi
else
  echo 'production requires --preflight-receipt-sha SHA256' >&2; exit 2
fi
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'submit with sbatch after lead review' >&2; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export NUMBA_DISABLE_JIT=0
export NUMBA_CACHE_DIR="/work/numba_cache"
mkdir -p "$remote/results" "$remote/cache"
[[ -f "$source/plan.json" && -f "$source/run_native_validation.py" ]] || { echo 'immutable validation source was not staged' >&2; exit 2; }
[[ ! -e "$remote/results/native_validation_v3" ]] || { echo 'refusing duplicate validation output' >&2; exit 2; }
if [[ "$mode" == "--self-test-only" ]]; then
  [[ ! -e "$remote/results/controller_mock_tests_v3" ]] || { echo 'refusing duplicate mock-test output' >&2; exit 2; }
fi
cd "$source"
python - <<'PY'
import hashlib, json, pathlib
p=pathlib.Path('plan.json'); d=pathlib.Path('run_native_validation.py'); plan=json.loads(p.read_text())
if plan['driver_sha256'].startswith('TO_BE_'): raise SystemExit('lead must pin driver and file hashes before submission')
if plan['driver_sha256'] != hashlib.sha256(d.read_bytes()).hexdigest(): raise SystemExit('driver pin mismatch')
if plan['schema'] != 'block0506_renter_no_taper_native_validation_v3': raise SystemExit('wrong plan schema')
if plan['total_seconds'] != 1200 or plan['case_seconds'] != 360: raise SystemExit('wrong time budget')
if plan['preflight_contract'] != {'required': True, 'receipt': 'controller_mock_tests.json', 'requires_matching_plan_and_driver': True, 'requires_matching_scientific_fields': True}: raise SystemExit('wrong preflight contract')
PY
python=(apptainer exec --bind "$reference:$original:ro,$source:/work/no_taper_source:ro,$remote/results:/work/results,$remote/cache:/work/numba_cache" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python)
if [[ "$mode" == "--self-test-only" ]]; then
  "${python[@]}" /work/no_taper_source/run_native_validation.py --plan /work/no_taper_source/plan.json --output /work/results/controller_mock_tests_v3 --self-test-only
  exit 0
fi
PREFLIGHT_SHA="$preflight_sha" "${python[@]}" - <<'PY'
import hashlib, json, os, pathlib
receipt=pathlib.Path('/work/results/controller_mock_tests_v3/controller_mock_tests.json')
plan=pathlib.Path('/work/no_taper_source/plan.json')
if not receipt.is_file(): raise SystemExit('successful Torch preflight receipt is required before lifecycle cases')
actual=hashlib.sha256(receipt.read_bytes()).hexdigest()
if actual != os.environ['PREFLIGHT_SHA']: raise SystemExit('preflight receipt SHA differs from explicit production pin')
r=json.loads(receipt.read_text()); p=json.loads(plan.read_text())
fields=('schema','reference_label','reference_manifest_sha256','cases','maximum_lifecycle_solves','case_seconds','total_seconds','threads','memory_gib','grid_nodes','lambda_d','natural_credit','renormalize_births','fixed_prices_rents_psi_fiscal','compiled_mode','required_readout','overlay_parameters_path','overlay_parameters_sha256','base_driver_path','pinned_files')
scientific={k:p[k] for k in fields}
digest=hashlib.sha256(json.dumps(scientific,sort_keys=True,separators=(',',':')).encode()).hexdigest()
if r.get('status') != 'passed' or r.get('model_solves') != 0: raise SystemExit('preflight did not pass as a zero-solve test')
if r.get('plan_sha256') != hashlib.sha256(plan.read_bytes()).hexdigest() or r.get('driver_sha256') != p['driver_sha256']: raise SystemExit('preflight plan/driver differs')
if r.get('scientific_plan_fields') != scientific or r.get('scientific_plan_sha256') != digest: raise SystemExit('preflight scientific fields differ')
PY
"${python[@]}" /work/no_taper_source/run_native_validation.py --plan /work/no_taper_source/plan.json --output /work/results/native_validation_v3
