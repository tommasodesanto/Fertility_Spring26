"""Zero-model-solve checks for the broad initial-vector contract."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
design=json.loads((HERE/'design.json').read_text())
plan=json.loads((PACKET/'plan.json').read_text())
assert hashlib.sha256((PACKET/'plan.json').read_bytes()).hexdigest()==design['original_plan_sha256']
assert hashlib.sha256((PACKET/'restart_controller_v2/controller.py').read_bytes()).hexdigest()==design['reviewed_censor_source_sha256']
assert len(design['starts'])==16
assert design['maximum_objective_calls_per_chain']==80
assert design['wall_seconds_per_chain']==4500 and design['final_reserve_seconds']==900
assert design['absolute_deadline_epoch']==1790932500
for arm_index,arm in enumerate(('hard','quarter')):
 rows=design['starts'][8*arm_index:8*(arm_index+1)]
 assert [r['slot'] for r in rows]==list(range(8*arm_index,8*(arm_index+1)))
 assert all(r['arm']==arm and r['weight_profile']=='base_control' for r in rows)
 for k,(lo,hi) in design['ranges'].items():
  vals=[r['parameters'][k] for r in rows]
  bins=[int((v-lo)/(hi-lo)*8) for v in vals]
  assert sorted(bins)==list(range(8)),(arm,k,bins)
  assert all(plan['bounds'][k][0]<=v<=plan['bounds'][k][1] for v in vals)
assert all(set(r['parameters'])==set(plan['free_coordinates']) for r in design['starts'])
print(json.dumps(dict(status='passed_zero_model_solves',arms=2,starts_per_arm=8,
                      dimensions=10,distinct_strata_per_dimension=8)))
