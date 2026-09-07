"""Replay original failed case18 and measure pruning loss without changing its result."""
from pathlib import Path
from types import SimpleNamespace
import ast, copy, gzip, inspect, json, pickle, sys, textwrap, time
import numpy as np

source_root=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907g')
root=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_mass_diagnosis_20260907')
sys.path[:0]=[str(source_root/'code/model'),str(source_root/'code/model/tools')]
import run_e5f_joint_overnight_case as a
import run_e5f_transition_calibration as calibration
from intergen_eqscale_seq_optimized import solver
origin=source_root/'output/model/joint_nested_overnight/search/initial_population_part_01'
execution=a.read_json(origin/'task_018/execution_contract.json')
plan=a.load_plan(origin/'plan.json',execution['plan_sha256'])
case=next(c for c in plan['cases'] if c['id']==18)
a.verify(origin/case['center'],case['center_sha256'])
out=root/'output/model/joint_nested_overnight/branch_mass_diagnosis'

# Construct an observational alternative of the same linear transition.
# Only its three absolute-mass pruning comparisons change; never return it
# to the original case or treat the replay as a completed calibration.
tree=ast.parse(textwrap.dedent(inspect.getsource(solver.advance_cohort_one_period_markov_income)))
changed=[]
for node in ast.walk(tree):
    if isinstance(node,ast.If) and isinstance(node.test,ast.Compare):
        t=node.test
        if len(t.ops)==1 and isinstance(t.ops[0],ast.Lt) and len(t.comparators)==1 and isinstance(t.comparators[0],ast.Constant) and t.comparators[0].value==1e-15:
            changed.append(ast.unparse(t));t.ops=[ast.LtE()];t.comparators=[ast.Constant(value=0.0)]
assert len(changed)==3,changed
scope=dict(solver.__dict__)
exec(compile(ast.fix_missing_locations(tree),'<observational_zero_only_transition>','exec'),scope)
strict=scope['advance_cohort_one_period_markov_income']
if '--preflight' in sys.argv:
    print(json.dumps(dict(status='preflight_passed',original_plan_sha256=execution['plan_sha256'],original_center_sha256=case['center_sha256'],original_comparisons=changed,source_bundle=a.BUNDLE,full_history_cases=1,case_cap_seconds=3600,outputs_created=False)),flush=True)
    sys.exit(0)
out.mkdir(parents=True,exist_ok=False)
(out/case['center']).write_bytes((origin/case['center']).read_bytes())
new=copy.deepcopy(plan);new.update(cases=[case],launch_deadline_epoch=time.time()+3600)
a.write_json(out/'plan.json',new)
original=calibration.begin_dated_first_birth_housing_branch

def observed(*args,**kwargs):
    try:return original(*args,**kwargs)
    except RuntimeError as error:
        if 'treated_branch_advancement mass gate failed:' not in str(error):raise
        tb=error.__traceback__;frame=None
        while tb:
            if tb.tb_frame.f_code.co_name=='advance':frame=dict(tb.tb_frame.f_locals)
            tb=tb.tb_next
        if frame is None:raise RuntimeError('Missing original advance frame') from error
        evaluation,P,b_grid,shared=args[:4]
        branch=frame['branch'];policy=frame['selected_policy'];expected=float(frame['expected_survivor_mass'])
        assert bool(P.use_numba_scatter) and solver.NUMBA_AVAILABLE
        _,_,income=solver.income_transition_values(P)
        ust=bool(getattr(P,'use_stochastic_aging',False) and hasattr(P,'Pi_child'))
        Pia=P.Pi_child if ust else None
        unpruned=np.zeros_like(branch);rescaled=np.zeros_like(branch);rows=[]
        for j in range(int(P.J)-1):
            survival=float(P.survival_probs[j]) if bool(getattr(P,'use_age_survival',False)) else 1.0
            cohort=survival*branch[:,:,:,j,:,:,:];mass=float(cohort.sum())
            if mass==0:continue
            tail=(j,policy.loc_probs,policy.tenure_choice,policy.tenure_probs,policy.bp_pol,P,b_grid,shared,policy.maps.lmm_idx,policy.maps.lmm_wt,policy.maps.tmx_idx,policy.maps.tmx_wt,ust,Pia,income)
            original_next=solver.advance_cohort_one_period_markov_income(cohort,*tail)
            strict_next=strict(cohort,*tail)
            scaled_next=mass*solver.advance_cohort_one_period_markov_income(cohort/mass,*tail)
            unpruned[:,:,:,j+1,:,:,:]=strict_next;rescaled[:,:,:,j+1,:,:,:]=scaled_next
            rows.append(dict(age_index=j,input_mass=mass,original_mass=float(original_next.sum()),zero_only_mass=float(strict_next.sum()),rescaled_mass=float(scaled_next.sum()),original_lost_mass=float((strict_next-original_next).sum()),strict_vs_scaled_l1=float(np.abs(strict_next-scaled_next).sum())))
        strict_gap=abs(float(unpruned.sum())-expected)/expected
        scaled_gap=abs(float(rescaled.sum())-expected)/expected
        packet=dict(branch=branch,selected_policy=policy,parameters=P,b_grid=b_grid,shared=shared,expected_mass=expected,origin_period=kwargs['origin_period'],evaluation=evaluation)
        with gzip.open(out/'failed_branch_state.pkl.gz','wb') as f:pickle.dump(packet,f,protocol=pickle.HIGHEST_PROTOCOL)
        trace=dict(status='original_failure_reproduced',instrumentation_only=True,original_error=str(error),original_plan_sha256=execution['plan_sha256'],original_center_sha256=case['center_sha256'],source_bundle=a.BUNDLE,solver_sha256=a.digest(solver.__file__),calibration_sha256=a.digest(calibration.__file__),original_comparisons=changed,origin_period=kwargs['origin_period'],original_mass=float(frame['next_pre'].sum()),expected_mass=expected,zero_only_mass=float(unpruned.sum()),zero_only_relative_gap=strict_gap,rescaled_mass=float(rescaled.sum()),rescaled_relative_gap=scaled_gap,original_gate_tolerance=5e-9,zero_only_passes_unchanged_gate=bool(strict_gap<=5e-9),rescaled_passes_unchanged_gate=bool(scaled_gap<=5e-9),cohorts=rows,saved_state_sha256=a.digest(out/'failed_branch_state.pkl.gz'),production_promoted=False)
        a.write_json(out/'branch_mass_trace.json',trace)
        print(json.dumps(trace),flush=True)
        raise

calibration.begin_dated_first_birth_housing_branch=observed
try:
    a.run_case(SimpleNamespace(plan=out/'plan.json',plan_sha256=a.digest(out/'plan.json'),case_id=18))
except RuntimeError as error:
    if 'treated_branch_advancement mass gate failed:' not in str(error) or not (out/'branch_mass_trace.json').exists():raise
    print('DIAGNOSTIC_COMPLETE: original failure retained; no repaired model or calibration',flush=True)
else:raise RuntimeError('Original mass failure did not reproduce')
finally:
    calibration.begin_dated_first_birth_housing_branch=original
