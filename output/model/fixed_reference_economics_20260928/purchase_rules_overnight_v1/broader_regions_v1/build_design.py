"""Create the 16 fixed, broad starting vectors; no model solve."""
import hashlib
import json
import random
from pathlib import Path

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
plan_bytes=(PACKET/'plan.json').read_bytes()
plan=json.loads(plan_bytes)
ranges={
 'beta_annual':(.955,.98),'chi':(.85,1.3),
 'first_birth_fixed_cost':(.2,.65),'kappa_fert':(.06,.22),
 'kappa_fert_continuation':(.18,.65),'theta0':(.065,.17),
 'child_benefit_curvature':(.05,.2),'tenure_choice_kappa':(.007,.03),
 'psi_child':(.1,.25),'h_P':(1.7,2.6),
}
assert set(ranges)==set(plan['free_coordinates'])
for k,(lo,hi) in ranges.items():
 blo,bhi=plan['bounds'][k]
 assert blo<=lo<hi<=bhi,(k,lo,hi,blo,bhi)
starts=[]
for arm_index,arm in enumerate(('hard','quarter')):
 rng=random.Random(20261002+arm_index)
 cols={}
 for k,(lo,hi) in ranges.items():
  bins=list(range(8));rng.shuffle(bins)
  cols[k]=[lo+(hi-lo)*(j+.5)/8 for j in bins]
  assert sorted(bins)==list(range(8))
 for row in range(8):
  slot=8*arm_index+row
  starts.append(dict(slot=slot,arm=arm,original_chain_index=row if arm=='hard' else 24+row,
                     parameters={k:cols[k][row] for k in ranges},weight_profile='base_control',
                     initial_price=plan['initial_price']))
design=dict(status='experimental_initialization_only',method='deterministic_stratified_midpoint_latin_hypercube',
            seed_hard=20261002,seed_quarter=20261003,original_plan_sha256=hashlib.sha256(plan_bytes).hexdigest(),
            original_target_contract_sha256=plan['target_fingerprint'],
            original_weight_fingerprint=plan['weight_fingerprint'],ranges=ranges,
            reviewed_censor_source_sha256=hashlib.sha256((PACKET/'restart_controller_v2/controller.py').read_bytes()).hexdigest(),
            full_approved_bounds=plan['bounds'],free_coordinates=plan['free_coordinates'],
            scored_targets=10,total_targets=14,starts=starts,
            economic_changes=['initial optimizer vectors only; no model, target, weight or bound change'],
            maximum_objective_calls_per_chain=80,wall_seconds_per_chain=4500,
            final_reserve_seconds=900,absolute_deadline_epoch=1790932500)
(HERE/'design.json').write_text(json.dumps(design,sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(starts=len(starts),plan_sha256=design['original_plan_sha256'],
                      first=starts[0]['parameters'],last=starts[-1]['parameters'])))
