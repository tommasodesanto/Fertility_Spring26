from pathlib import Path
import copy,sys,time
import numpy as np
root=Path(__file__).resolve().parent
sys.path[:0]=[str(root/'code/model'),str(root/'code/model/tools')]
import run_e5f_independent_numerical_audit as audit
import run_e5f_joint_overnight_case as a
import run_e5f_transition_calibration as cal
from intergen_eqscale_seq_optimized import solver
for name,sha in a.read_json(root/'source_manifest.json').items():a.verify(root/name,sha)
assert cal.code_fingerprint_contract(solver)['bundle_sha256']==a.BUNDLE
source=Path('/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907j/output/model/joint_nested_overnight/policy_loop_smoke/supply-plus-20/date_2023/dated_state.pkl.gz')
a.verify(source,'ffdd4ee394db635f0278a99306177e01eee9ee6d42312670c22c7e0a4b4e6639')
packet=audit.load_checkpoint(source);P=copy.deepcopy(packet['parameters']);old=packet['evaluation'];bg=packet['b_grid']
out=root/'output/model/joint_nested_overnight/support_repair';out.mkdir(parents=True,exist_ok=False)
before=out/'before';before.mkdir();after=out/'after';after.mkdir()
budget0=audit.budget_audit(packet,before);assert budget0['budget_excess_mass']==3.727996554261621e-10
cal.transition.configure_sequential_model();audit.calendar.apply_fertility=cal.transition.apply_sequential_fertility;audit.calendar.advance_calendar_distribution=cal.transition.advance_sequential_calendar_distribution
shared=solver.precompute_shared(P,bg);start=time.time()
ev=audit.calendar.evaluate_period(old.policy.price,old.inherited_g_pre.copy(),P,bg,shared,audit.calendar.SolveCounter(),packet['supply_rule'])
current=dict(packet,parameters=P,evaluation=ev,shared=shared)
budget1=audit.budget_audit(current,after);arrays=audit.policy_array_audit(current,after)
audit.standard_diagnostics(current,after,validate_production_young=False)
assert budget1['budget_excess_mass']<=2e-10 and arrays['occupied_negative_steps']==0
assert ev.feasibility_projection_mass<=1e-6
assert abs(float(ev.g_current.sum())-float(old.g_current.sum()))<=2e-10
for age in (58,62):
 j=int(round((age-P.age_start)/P.da));b=int(np.argmin(abs(bg+.25581395348837255)))
 assert ev.g_current[b,0,0,j,0,3,3]==0
proof=dict(status='pass_original_failed_policy_state',checkpoint_sha256=a.digest(source),new_bundle=a.BUNDLE,elapsed_seconds=time.time()-start,budget_before=budget0,budget_after=budget1,policy_arrays=arrays,projection_mass=ev.feasibility_projection_mass,price=ev.policy.price.tolist(),demand_before=old.demand_by_loc.tolist(),demand_after=ev.demand_by_loc.tolist(),births_before=old.births,births_after=ev.births,scope='Original inherited population and price; requires fresh recleared histories/policy paths before calibration',production_promoted=False)
a.write_json(out/'original_state_verification.json',proof);print(proof,flush=True)
