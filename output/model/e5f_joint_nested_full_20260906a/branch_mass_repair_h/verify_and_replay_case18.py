"""Verify the precise repair on the original failed state, then replay its full case."""
from pathlib import Path
from types import SimpleNamespace
import copy,gc,gzip,json,pickle,sys,time
import numpy as np
root=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907h')
old=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907g')
trace_root=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_mass_diagnosis_20260907/output/model/joint_nested_overnight/branch_mass_diagnosis')
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_joint_overnight_case as a
import run_e5f_transition_calibration as cal
from intergen_eqscale_seq_optimized import solver
for rel,sha in a.read_json(root/'source_manifest.json').items():a.verify(root/rel,sha)
assert cal.code_fingerprint_contract(solver)['bundle_sha256']==a.BUNDLE
out=root/'output/model/joint_nested_overnight/mass_repair_replay';out.mkdir(parents=True,exist_ok=False)

def verify_state():
    trace=a.read_json(trace_root/'branch_mass_trace.json')
    assert trace['source_bundle']=='50b4342797eb271e71b15651c60f4e45c6c740eb6205dd609fcd32501c7428eb'
    assert trace['zero_only_passes_unchanged_gate'] and trace['origin_period']==2
    a.verify(trace_root/'failed_branch_state.pkl.gz',trace['saved_state_sha256'])
    with gzip.open(trace_root/'failed_branch_state.pkl.gz','rb') as f:packet=pickle.load(f)
    cal.transition.configure_sequential_model()
    branch=cal.begin_dated_first_birth_housing_branch(packet['evaluation'],packet['parameters'],packet['b_grid'],packet['shared'],origin_period=packet['origin_period'])
    gates=branch['branch_mass_roundoff_gates'];assert len(gates)==2
    for gate in gates:
        assert gate['relative_tolerance']==5e-9 and gate['relative_gap']<=5e-9
    treated=gates[0]
    assert treated['initial_actual_mass']==trace['original_mass']
    assert treated['actual_mass_before_normalization']==trace['zero_only_mass']
    assert treated['initial_relative_gap']>5e-9
    assert treated['transport_precision']=='positive_mass_retained_after_original_gate_failure'
    assert np.isfinite(branch['treated_next_pre']).all() and np.min(branch['treated_next_pre'])>=0
    assert np.isfinite(branch['control_next_pre']).all() and np.min(branch['control_next_pre'])>=0
    proof=dict(status='pass_saved_original_failure_state',trace_sha256=a.digest(trace_root/'branch_mass_trace.json'),saved_state_sha256=trace['saved_state_sha256'],new_bundle=a.BUNDLE,solver_sha256=a.digest(solver.__file__),calibration_sha256=a.digest(cal.__file__),origin_period=branch['origin_period'],origin_mass=branch['origin_mass'],survivor_mass=branch['survivor_mass_before_destination_gating'],gates=gates,full_replay_complete=False,production_promoted=False)
    a.write_json(out/'saved_state_repair_verification.json',proof)
    print(json.dumps(proof),flush=True)
verify_state();gc.collect()
origin=old/'output/model/joint_nested_overnight/search/initial_population_part_01'
a.verify(origin/'plan.json','3884c6202ecefdc018ae2a53d4f4f506731708524dc79132a7b3877ab6fc9a71')
plan=a.read_json(origin/'plan.json');case=next(c for c in plan['cases'] if c['id']==18)
a.verify(origin/case['center'],'40eb3f4df91b5de25925112c12a5aee94464ca2220c9e5042939a3543d75a0c4')
for name,sha in plan['helper_sha256'].items():a.verify(root/'code/model/tools'/name,sha)
(out/case['center']).write_bytes((origin/case['center']).read_bytes())
new=copy.deepcopy(plan);new.update(cases=[case],code_bundle_sha256=a.BUNDLE,adapter_sha256=a.digest(a.__file__),controller_sha256=a.digest(root/'code/model/tools/run_e5f_joint_nested_long_search.py'),planner_sha256=a.digest(root/'code/model/tools/build_e5f_bounded_refinement_plan.py'),launch_deadline_epoch=time.time()+3600)
a.write_json(out/'plan.json',new)
try:a.run_case(SimpleNamespace(plan=out/'plan.json',plan_sha256=a.digest(out/'plan.json'),case_id=18))
finally:
    a.write_json(out/'replay_provenance.json',dict(original_plan_sha256=a.digest(origin/'plan.json'),original_center_sha256=case['center_sha256'],new_plan_sha256=a.digest(out/'plan.json'),new_bundle=a.BUNDLE,complete=(out/'task_018/case_receipt.json').exists(),production_promoted=False))
