"""Reproduce the saved budget failure, verify its correction, replay original case17."""
from pathlib import Path
from types import SimpleNamespace
import copy,gc,json,sys,time
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_joint_overnight_case as a
import run_e5f_independent_numerical_audit as audit
import run_e5f_transition_calibration as cal
from intergen_eqscale_seq_optimized import solver
for rel,sha in a.read_json(root/'source_manifest.json').items():a.verify(root/rel,sha)
assert cal.code_fingerprint_contract(solver)['bundle_sha256']==a.BUNDLE
origin=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907i/output/model/joint_nested_overnight/search/initial_population')
a.verify(origin/'plan.json','d228918f2344f5e09c585baf8a428cae26facde50d2824e7fc141ef58c917075')
out=root/'output/model/joint_nested_overnight/renter_repair_replay';out.mkdir(parents=True,exist_ok=False)
def verify_state():
    checkpoint=origin/'task_017/dated_state.pkl.gz'
    packet=audit.load_checkpoint(checkpoint);P=copy.deepcopy(packet['parameters']);old=packet['evaluation'];bg=packet['b_grid']
    before=out/'before';before.mkdir();after=out/'after';after.mkdir()
    budget0=audit.budget_audit(packet,before)
    assert budget0['budget_excess_mass']==2.426372670781678e-10
    audit.transition.configure_sequential_model()
    audit.calendar.apply_fertility=audit.transition.apply_sequential_fertility
    audit.calendar.advance_calendar_distribution=audit.transition.advance_sequential_calendar_distribution
    shared=solver.precompute_shared(P,bg);start=time.time()
    ev=audit.calendar.evaluate_period(old.policy.price,old.inherited_g_pre.copy(),P,bg,shared,audit.calendar.SolveCounter(),packet['supply_rule'])
    current=dict(packet,parameters=P,shared=shared,evaluation=ev)
    names=('V','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value','fert2_probs')
    exact={name:bool(np.array_equal(getattr(old.policy,name),getattr(ev.policy,name))) for name in names}
    exact.update({name:bool(np.array_equal(getattr(old,name),getattr(ev,name))) for name in ('g_pre','g_post_fertility','g_current')})
    exact.update({name:bool(np.array_equal(getattr(old.policy.joint_choice,name),getattr(ev.policy.joint_choice,name))) for name in ('probabilities','products','wait_probabilities')})
    assert all(exact.values()),exact
    budget1=audit.budget_audit(current,after);arrays=audit.policy_array_audit(current,after)
    assert budget1['budget_excess_mass']<=2e-10 and arrays['occupied_negative_steps']==0
    # Independently reconstruct the original worst renter state's allocation.
    row=a.read_csv(origin/'task_017/worst_occupied_budget_states.csv')[0]
    j=int(round((float(row['age'])-P.age_start)/P.da));b=int(np.argmin(abs(bg-float(row['wealth']))));nn=int(row['parity']);cs=int(row['dependent_children']);flat=nn+P.n_parity*cs
    idx=(b,0,0,j,int(row['income_index']),nn,cs)
    cb=float(shared.cb_flat.reshape(-1)[flat]);hb=float(shared.hb_flat.reshape(-1)[flat]);al=float(shared.alpha_flat.reshape(-1)[flat])
    rent=float(P.user_cost_rate*old.policy.price[0]);resources=float(row['resources']);saving=float(row['saving'])
    surplus=resources-cb-rent*hb-saving;h=hb+min((1-al)*surplus/rent,P.hR_max-hb);c=resources-rent*h-saving
    assert abs(c-float(ev.policy.c_pol[idx]))<1e-14 and abs(h-float(ev.policy.hR_pol[idx]))<1e-14
    dh=ev.policy.hR_pol-old.policy.hR_pol
    proof=dict(status='pass_original_state_fixed_prices',old_checkpoint=str(checkpoint),old_checkpoint_sha256=a.digest(checkpoint),new_bundle=a.BUNDLE,elapsed_seconds=time.time()-start,budget_before=budget0,budget_after=budget1,policy_arrays=arrays,exact_unchanged_arrays=exact,
      original_worst_state=dict(row,optimizer_consumption=c,optimizer_housing=h,reported_housing_before=float(old.policy.hR_pol[idx]),reported_consumption_after=float(ev.policy.c_pol[idx])),
      housing_changed_cells=int(np.count_nonzero(dh)),housing_changed_current_mass=float(ev.g_current[dh!=0].sum()),weighted_absolute_housing_change=float(np.sum(ev.g_current*np.abs(dh))),demand_before=old.demand_by_loc.tolist(),demand_after=ev.demand_by_loc.tolist(),
      scope='Fixed-price diagnostic only. Fresh histories and recleared policies required; no production promotion.')
    a.write_json(out/'saved_state_repair_verification.json',proof);print(json.dumps(proof),flush=True)
verify_state();gc.collect()
plan=a.read_json(origin/'plan.json');case=next(c for c in plan['cases'] if c['id']==17)
a.verify(origin/case['center'],case['center_sha256'])
(out/case['center']).write_bytes((origin/case['center']).read_bytes())
new=copy.deepcopy(plan);new.update(cases=[case],code_bundle_sha256=a.BUNDLE,adapter_sha256=a.digest(a.__file__),controller_sha256=a.digest(root/'code/model/tools/run_e5f_joint_nested_long_search.py'),planner_sha256=a.digest(root/'code/model/tools/build_e5f_bounded_refinement_plan.py'),launch_deadline_epoch=min(time.time()+3600,1788788100))
a.write_json(out/'plan.json',new)
try:a.run_case(SimpleNamespace(plan=out/'plan.json',plan_sha256=a.digest(out/'plan.json'),case_id=17))
finally:a.write_json(out/'replay_provenance.json',dict(original_plan_sha256=a.digest(origin/'plan.json'),original_center_sha256=case['center_sha256'],new_plan_sha256=a.digest(out/'plan.json'),new_bundle=a.BUNDLE,complete=(out/'task_017/case_receipt.json').exists(),production_promoted=False))
