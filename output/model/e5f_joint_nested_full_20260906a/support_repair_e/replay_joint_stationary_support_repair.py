"""Replay original failed wide candidate with the reviewed support-status repair."""
from pathlib import Path
from types import SimpleNamespace
import copy,json,sys,time
root=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907e')
old=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907d')
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_joint_overnight_case as a
from intergen_eqscale_seq_optimized import solver
original=old/'output/model/joint_nested_overnight/search/initial_population_part_01'
plan=a.read_json(original/'plan.json')
a.verify(original/'plan.json','be5f74af91913c7649fb1505e56ad27ea9796d611ace0afc9c53528e29fb9462')
case=next(c for c in plan['cases'] if c['id']==26)
a.verify(original/case['center'],case['center_sha256'])
out=root/'output/model/joint_nested_overnight/support_replay';out.mkdir(exist_ok=False)
(out/case['center']).write_bytes((original/case['center']).read_bytes())
new=copy.deepcopy(plan);new.update(cases=[case],code_bundle_sha256=a.BUNDLE,adapter_sha256=a.digest(a.__file__),launch_deadline_epoch=time.time()+3600)
a.write_json(out/'plan.json',new)
source_function=solver.forward_distribution_markov_income
missing=[]
def observe(*args,**kwargs):
    result=source_function(*args,**kwargs)
    stats=result[1]
    if getattr(stats,'stationary_eventstudy_status',None)=='undefined_first_birth_support':
        missing.append(dict(status=stats.stationary_eventstudy_status,masses=stats.stationary_eventstudy_branch_masses))
        a.write_json(out/'unavailable_intermediate_support.json',dict(count=len(missing),trials=missing))
    return result
solver.forward_distribution_markov_income=observe
try:
    a.run_case(SimpleNamespace(plan=out/'plan.json',plan_sha256=a.digest(out/'plan.json'),case_id=26))
finally:
    solver.forward_distribution_markov_income=source_function
    a.write_json(out/'replay_provenance.json',dict(original_plan_sha256=a.digest(original/'plan.json'),
        original_center_sha256=case['center_sha256'],new_plan_sha256=a.digest(out/'plan.json'),
        new_bundle=a.BUNDLE,instrumentation='Records missing-support status; unchanged returned arrays and statistics',
        unavailable_trial_count=len(missing),complete=(out/'task_026/case_receipt.json').exists()))
