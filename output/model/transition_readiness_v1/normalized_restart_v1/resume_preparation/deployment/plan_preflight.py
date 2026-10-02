"""Invoke the existing pinned plan/preparation validators; zero model calls."""
import argparse,importlib.util,json,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--plan',required=True);p.add_argument('--output',required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[6];code=root/'code/model/experiments/transition_readiness';sys.path.insert(0,str(code))
import one_shock_floor as c
plan=json.loads(Path(a.plan).read_text());result=c.preflight(plan)
prepared=json.loads(c.pinned(plan['prepared_native_inputs']).read_text());measured=c.validate_preparation(plan,prepared)
import numpy as np
matrix=np.load(c.pinned(measured['matrix']));assert matrix.shape==(24,24) and np.isfinite(matrix).all()
result.update(status='exact_mounted_preparation_preflight_passed',native_calls=0,reference_checkpoint_pin_verified=True,measured_matrix_pin_verified=True,approved_controller_critical_AST_verified=True,original_preparation_sha256=plan['prepared_native_inputs']['sha256'])
if 'diagnostic_measurement_reuse' in plan:
    reuse=json.loads(c.pinned(plan['diagnostic_measurement_reuse']).read_text())
    cases=[c.reused_diagnostic_measurement(plan,case['psi']) for case in reuse['cases']]
    assert len(cases)==4 and all(case['accepted_for_relaxed_diagnostic'] for case in cases)
    result.update(authenticated_completed_measurements_verified=4,cached_native_calls=0)
c.write(Path(a.output)/'fit_preflight.json',result);print(json.dumps(result))
