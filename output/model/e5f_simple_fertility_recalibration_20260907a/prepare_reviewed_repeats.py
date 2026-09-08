"""Prepare the two already-budgeted exact repeats after reviewing rejected case5."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'code/model/tools'))
import run_e5f_bounded_calibration_refinement as adapter
import build_e5f_bounded_refinement_plan as planner
base=Path('output/model/e5f_simple_fertility_recalibration_20260907a').resolve()
contract=adapter.read_json(base/'contract.json')
for path,sha in contract['code_sha256'].items():adapter.verify(path,sha)
path=base/'search/joint/plan.json'
sha='565705df82c0fba636accc33fc9675e2133219d2400bcaf7784914915cb036a6'
plan=adapter.load_plan(path,sha)
status=adapter.read_json(path.parent/'report/summary.json')
assert status['plan_sha256']==sha and status['completed_cases']==11 and status['missing_cases']==[5]
failure=adapter.read_json(path.parent/'task_005/adapter_failure.json')
assert failure['error']=='Housing market did not clear: residual=3.588e-04.'
best=adapter.read_json(base/'search/best_so_far.json')['best']
assert best['id']==12 and best['loss']==26.249682727266702
rows=adapter.read_csv(path.parent/'report/all_candidates.csv')
planner.repeats(plan,status,rows,base/'search/reviewed_repeats')
