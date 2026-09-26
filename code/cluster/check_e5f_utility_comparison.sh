#!/usr/bin/env bash
# Torch only. Bounded unit and native zero-solve checks; no calibration launch.
set -euo pipefail
bundle="${1:?path to isolated preparation bundle}"
result_name="${2:?fresh result directory name}"
reference=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$bundle/tools:$reference/source/code/model/tools:$reference/source/code/model:/scratch/td2248/commute_pdf_qa_deps"
result="$bundle/$result_name"
test ! -e "$result"
mkdir "$result"
python3 -m unittest -v test_e5f_utility_comparison_runtime test_e5f_utility_comparison_design \
  test_run_e5f_utility_comparison test_e5f_utility_comparison_search \
  test_collect_e5f_utility_comparison > "$result/tests.log" 2>&1
python3 "$bundle/tools/run_e5f_utility_comparison.py" --stage prepare \
  --reference-root "$reference" --output "$result/preparation" \
  --provenance "$bundle/inputs/proposed_common_contract.json" \
  --pension-receipt "$bundle/inputs/pension_receipt.json" \
  --pension-contract "$bundle/inputs/pension_measurement_contract.json" \
  --timing-receipt "$bundle/inputs/observed_timing.json" \
  --workers-per-arm 10 --total-seconds 28800 --repeat-reserve-seconds 4500 \
  --export-reserve-seconds 900 --overhead-seconds 120 > "$result/prepare.log"
EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256=$(sha256sum "$result/preparation/contract.json" | cut -d ' ' -f 1)
export EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256
for arm in floor_linear floor_concave shares_linear shares_concave; do
  python3 "$bundle/tools/run_e5f_utility_comparison.py" --stage preflight \
    --contract "$result/preparation/contract.json" --arm "$arm" \
    --output "$result/native_$arm" > "$result/preflight_$arm.log" 2>&1
done
python3 - "$result" <<'PY'
import json, pathlib, sys
root=pathlib.Path(sys.argv[1])
arms={p.parent.name:json.loads(p.read_text()) for p in root.glob('native_*/preflight.json')}
assert len(arms)==4 and all(v['native_solve_count']==0 for v in arms.values())
result={'status':'focused_tests_and_four_native_preflights_passed','native_solve_count':0,
        'arms':{k:{'status':v['status'],'tax':v['fiscal_rule']['payroll_tax'],
                   'bound_seeds':len(v['records'])} for k,v in arms.items()},
        'tests_log':str(root/'tests.log'),'production_launch_authorized':False}
(root/'validation_receipt.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,sort_keys=True))
PY
