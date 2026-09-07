"""Instrumentation-only replay of failed wide case26; original exception retained."""
from pathlib import Path
from types import SimpleNamespace
import copy, json, sys, time
import numpy as np

root=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907d')
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_joint_overnight_case as a
from intergen_eqscale_seq_optimized import joint_nested as joint

origin=root/'output/model/joint_nested_overnight/search/initial_population_part_01'
execution=a.read_json(origin/'task_026/execution_contract.json')
plan=a.load_plan(origin/'plan.json',execution['plan_sha256'])
case=next(c for c in plan['cases'] if c['id']==26)
a.verify(origin/case['center'],case['center_sha256'])
out=root/'output/model/joint_nested_overnight/stationary_support_diagnosis'
out.mkdir(exist_ok=False)
(out/case['center']).write_bytes((origin/case['center']).read_bytes())
new=copy.deepcopy(plan);new.update(cases=[case],launch_deadline_epoch=time.time()+540)
a.write_json(out/'plan.json',new)
original=joint.stationary_first_birth_response

def observed(*args,**kwargs):
    try:return original(*args,**kwargs)
    except RuntimeError as error:
        if str(error)!='Invalid stationary matched joint branch mass':raise
        tb=error.__traceback__;values=None
        while tb:
            if tb.tb_frame.f_code.co_name=='stationary_first_birth_response':values=tb.tb_frame.f_locals
            tb=tb.tb_next
        if values is None:raise
        g,j,P=args[:3];masses=values['masses'];means=values['means']
        result=dict(status='original_failure_reproduced',instrumentation_only=True,
            original_plan_sha256=execution['plan_sha256'],original_center_sha256=case['center_sha256'],
            source_bundle=a.BUNDLE,stationary_function_source_sha256=a.digest(joint.__file__),
            masses=masses,means=means,mass_difference=abs(masses[0]-masses[1]),
            unequal_mass_condition=bool(abs(masses[0]-masses[1])>2e-10),
            low_support_condition=bool(min(masses)<=1e-14),
            pre_choice_mass=float(g.sum()),psi_child=float(P.psi_child),
            kappa=float(P.tenure_choice_kappa),nest_lambda=float(P.joint_nest_lambda),
            probability_min=float(np.min(j.probabilities)),probability_max=float(np.max(j.probabilities)),
            probability_nonfinite=int(np.count_nonzero(~np.isfinite(j.probabilities))),
            original_exception=str(error))
        a.write_json(out/'support_trace.json',result)
        print(json.dumps(result),flush=True)
        raise

joint.stationary_first_birth_response=observed
try:
    a.run_case(SimpleNamespace(plan=out/'plan.json',plan_sha256=a.digest(out/'plan.json'),case_id=26))
except RuntimeError as error:
    if str(error)!='Invalid stationary matched joint branch mass' or not (out/'support_trace.json').exists():raise
    print('DIAGNOSTIC_COMPLETE: original failure retained; no repaired model or calibrated result',flush=True)
else:
    raise RuntimeError('Expected original failure did not reproduce')
