#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_fit_sensitivity_v1/numba_cache_production
export PYTHONPATH=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1/code/model/tools
python - /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_fit_sensitivity_v1 <<'CHECK'
from pathlib import Path
import hashlib,json,sys,datetime
r=Path(sys.argv[1]);m=json.loads((r/'launch_manifest.json').read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(r/'run_e5f_income_fit_sensitivity.py')==m['driver_sha256'],'new controller hash mismatch'
assert sha(m['prior_plan'])==m['prior_plan_sha256'],'parent plan hash mismatch'
assert sha(m['selected_summary'])==m['selected_summary_sha256'],'selected reference hash mismatch'
assert datetime.datetime.now(datetime.timezone.utc)<datetime.datetime.fromisoformat(m['finish_by_utc']),'absolute finish deadline already passed'
CHECK
remaining=$(python -c "import time,datetime;print(max(0,min(10200,int(datetime.datetime(2026,9,21,11,30,tzinfo=datetime.timezone.utc).timestamp()-time.time()))))")
[[ "$remaining" -gt 0 ]]
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "stage=production exit=%s\n" "$rc" > /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_fit_sensitivity_v1/production_failure.txt; fi' EXIT
timeout --signal=TERM --kill-after=20s "${remaining}s" python -B /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_fit_sensitivity_v1/run_e5f_income_fit_sensitivity.py --mode production --plan /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1/plan.json --selected-summary /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_overnight_v1/production/summary.json --output /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_fit_sensitivity_v1/production --workers 8 --case-timeout 900 --global-timeout "$remaining" --smoke-summary /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/income_fit_sensitivity_v1/smoke/smoke_summary.json
