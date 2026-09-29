#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
overlay=/scratch/td2248/projects/fixed_reference_transition_20260928
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
relative=output/model/fixed_reference_transition_20260928
packet="$relative/four_shock_v1/budget_diagnostic_v2"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$original/$packet/source/code/model/tools:$original/code/model/tools"
export E5F_DIAGNOSTIC_PACKET="$packet"
apptainer exec --bind "$stage:$original:ro" --bind "$overlay:$original/$relative:rw" --pwd "$original" \
 /share/apps/images/ubuntu-24.04.4.sif bash -s <<'PAYLOAD'
set -euo pipefail
python=/share/apps/anaconda3/2025.06/bin/python
packet="$E5F_DIAGNOSTIC_PACKET"
"$python" -m unittest -v test_e5f_preference_estimation test_e5f_preference_estimation_batch test_e5f_preference_transition test_e5f_four_shock_acceleration test_e5f_exact_policy_cache test_e5f_preference_budget_diagnostic > "$packet/verification_tests.log" 2>&1
"$python" - <<'REGRESSION'
import copy,importlib.util,json,os,pathlib
import numpy as np
from run_e5f_preference_transition import plain
p=pathlib.Path(os.environ['E5F_DIAGNOSTIC_PACKET']);old=p.parent/'budget_diagnostic_v1'
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
before=load('before',old/'source/code/cluster/run_e5f_preference_budget_diagnostic.py')
after=load('after',p/'source/code/cluster/run_e5f_preference_budget_diagnostic.py')
config=after.read(p/'configs/full.json')
endpoint=after.read(after.pinned(config['endpoint_receipt']))
baseline_path=old/'smoke/baseline/mapping.json';record=after.read(baseline_path)
terminal=copy.deepcopy(endpoint['terminal'])
for key in ('terminal_birth_queue','stationary_birth_queue'):terminal[key]=np.asarray(terminal[key],dtype=float)
receipt=after.summarize(record,terminal)
try:before.write(p/'old_writer_regression.json',receipt)
except TypeError as exc:
 assert 'ndarray' in str(exc),str(exc)
else:raise AssertionError('Old writer did not reproduce the observed failure')
after.write(p/'saved_data_roundtrip.json',receipt)
assert after.read(p/'saved_data_roundtrip.json')==plain(receipt)
assert plain(terminal)==plain(copy.deepcopy(terminal))
for mode in ('smoke','full'):
 c=after.read(p/'configs'/f'{mode}.json');after.validate_config(c);after.check_sources(c)
 if mode=='full':after.validate_baseline(after.read(after.pinned(c['baseline_mapping'])),104,dict(q=.7898695017462086,psi=.1355551166583114,pension=.917784047463731))
after.write(p/'reporting_regression.json',dict(status='PASS',old_ndarray_failure_reproduced=True,new_saved_data_roundtrip_exact=True,terminal_arrays_restored_to_native_types=True,new_terminal_comparison_passes=True,both_actual_configs_pass=True,model_solves=0,input_pins=dict(smoke_baseline=dict(path=str(baseline_path),sha256=after.sha(baseline_path)),endpoint=config['endpoint_receipt']),driver_sha256=after.sha(p/'source/code/cluster/run_e5f_preference_budget_diagnostic.py')))
print('Saved-data reporting regression and both actual-config preflights PASS; zero model solves')
REGRESSION
PAYLOAD
