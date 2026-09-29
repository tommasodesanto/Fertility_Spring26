#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=05:00:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
overlay=/scratch/td2248/projects/fixed_reference_transition_20260928
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
relative=output/model/fixed_reference_transition_20260928
packet="$relative/four_shock_v1/budget_diagnostic_v1"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR="$overlay/four_shock_v1/numba_cache"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
export E5F_DIAGNOSTIC_PACKET="$packet"
export PYTHONPATH="$original/$packet/source/code/model/tools:$original/code/model/tools"
apptainer exec --bind "$stage:$original:ro" --bind "$overlay:$original/$relative:rw" --pwd "$original" \
 /share/apps/images/ubuntu-24.04.4.sif bash -s <<'PAYLOAD'
set -euo pipefail
python=/share/apps/anaconda3/2025.06/bin/python
packet="$E5F_DIAGNOSTIC_PACKET"
"$python" -m unittest -v test_e5f_preference_estimation test_e5f_preference_estimation_batch test_e5f_preference_transition test_e5f_four_shock_acceleration test_e5f_exact_policy_cache test_e5f_preference_budget_diagnostic > "$packet/tests.log" 2>&1
"$python" - <<'PREFLIGHT'
import importlib.util,json,os,pathlib
p=pathlib.Path(os.environ['E5F_DIAGNOSTIC_PACKET'])
spec=importlib.util.spec_from_file_location('diagnostic_preflight',p/'source/code/cluster/run_e5f_preference_budget_diagnostic.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
for mode in ('smoke','full'):
 c=m.read(p/'configs'/f'{mode}.json');m.validate_config(c);m.check_sources(c)
 if mode=='full':
  m.validate_baseline(m.read(m.pinned(c['baseline_mapping'])),104,dict(q=.7898695017462086,psi=.1355551166583114,pension=.917784047463731))
m.write(p/'preflight.json',dict(status='PASS',both_actual_configs=True,full_saved_baseline_valid=True,tests_passed=True))
PREFLIGHT
for mode in smoke full; do
 config="$packet/configs/$mode.json"
 config_sha="$(cat "$packet/configs/$mode.sha256")"
 "$python" "$packet/source/code/cluster/run_e5f_preference_budget_diagnostic.py" --config "$config" --config-sha256 "$config_sha" --output "$packet/$mode"
 # Full runs only after smoke returns success; the driver also writes its result.
 "$python" -c 'import json,sys; r=json.load(open(sys.argv[1])); assert r["numerical_certified"] is True and r["fitted_shocks"] is False and r["production_ready"] is False' "$packet/$mode/complete.json"
done
PAYLOAD
