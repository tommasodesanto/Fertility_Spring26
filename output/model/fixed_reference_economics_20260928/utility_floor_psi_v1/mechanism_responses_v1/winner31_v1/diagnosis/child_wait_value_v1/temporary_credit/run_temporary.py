#!/usr/bin/env python3
"""Two native policy evaluations: baseline replay and one-period credit surprise."""
import os,sys,copy,time,signal,json,importlib.util,hashlib
from pathlib import Path
from types import SimpleNamespace
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).resolve().parent;WIN=HERE.parents[2];sys.path.insert(0,str(WIN))
import fixed_price_responses as d
import numpy as np
BASE=WIN/'purchase_ltv_v1/local_run/retry5/results/baseline_80_80';PERM=WIN/'purchase_ltv_v1/local_run/retry9/results/both_100_100';PATCH=WIN/'purchase_ltv_v1/local_run/household_buyer_estate.py'
FIELDS=['V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value']
def load(folder):
 with np.load(folder/'solution_arrays.npz',allow_pickle=False)as a:return SimpleNamespace(**{k:a[k].copy()for k in FIELDS},fert2_probs=a['fert2_probs'].copy(),bp_pol_stay=a['bp_pol_stay'].copy(),c_pol_stay=a['c_pol_stay'].copy())
def alarm(*_):raise TimeoutError('Native policy evaluation time budget')
def mean(x,w):return float(np.sum(np.nan_to_num(x)*w)/w.sum()) if w.sum()>0 else None
def invert(sol,P,pi,w):
 p=sol.fert_probs[:,0,0,:,:,:2];F=sol.fert_value[:,0,0,:,:];mask=(w>0)&(p[...,0]>0)&(p[...,1]>0)&(pi>0)
 W=np.full_like(F,np.nan);A=W.copy();G=W.copy();C=W.copy();W[mask]=F[mask]+P.kappa_fert*np.log(p[...,0][mask]);A[mask]=F[mask]+P.kappa_fert*np.log(p[...,1][mask]);G[mask]=A[mask]-W[mask];C[mask]=W[mask]+G[mask]/np.broadcast_to(pi,W.shape)[mask]
 tp=sol.tenure_probs[:,0,0,:,:,:,:,:];return dict(p=p[...,1],W=W,A=A,G=G,C=C,mask=mask,wait_size=tp[:,:,:,0,0,1:],child_size=tp[:,:,:,1,1,1:])
def main():
 start=time.time();deadline=start+300;auth=d.authenticate_candidate(HERE/'runtime_preparation');P0=auth['P'];grid=auth['grid'];q=float(d.BINDING['candidate_price']);base=load(BASE);perm=load(PERM)
 from refactor_lab.engine import household as h
 from refactor_lab.engine.parameters import get_fecundity_by_age
 d.require(d.sha(h.__file__)==json.loads((BASE.parent/'source_receipt.json').read_text())['original_household_sha256'],'Baseline native source drift')
 d.require(d.sha(PATCH)==json.loads((PERM.parent/'source_receipt.json').read_text())['override_sha256'],'Permanent estate-patched source drift')
 spec=importlib.util.spec_from_file_location('refactor_lab.engine.temporary_credit_household',PATCH);patch=importlib.util.module_from_spec(spec);sys.modules[spec.name]=patch;spec.loader.exec_module(patch)
 pre=np.load(BASE.parent/'q0_reference_inherited_states.npz')['g_pre'];w=pre[:,0,0,:,:,0,0].copy();w[:,7:]=0;pi=get_fecundity_by_age(P0)[None,:,None];occup=pre>0
 baseline_V=base.V.copy();original_next_hash=hashlib.sha256(baseline_V.tobytes()).hexdigest();attempts=[];saved={'baseline':base,'permanent_both100':perm}
 for label,phi,fn in [('baseline_replay',.8,h.solve_bellman_full_markov_income),('temporary_both100',1.,patch.solve_bellman_full_markov_income)]:
  P=copy.deepcopy(P0);P.phi=np.full_like(np.asarray(P0.phi),phi)
  if hasattr(P,'_purchase_ltv_override'):del P._purchase_ltv_override
  d.require(not getattr(P,'native_solvency_credit',False) and P.native_due_stayer_credit,'Credit contract drift')
  sd=auth['solver'].precompute_shared(P,grid);end=min(deadline,time.time()+150);d.require(end>time.time(),'Global solve deadline')
  d.write(HERE/'latest.json',dict(status='policy_evaluation_claimed',label=label,attempts=len(attempts)+1,deadline_epoch=end,price=q,continuation_V_sha256=original_next_hash))
  before=time.monotonic();old=signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,end-time.time())
  try:objects=fn(np.asarray([P.user_cost_rate*q]),np.asarray([q]),P,grid,sd,continuation_V=baseline_V)
  finally:signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
  sol=SimpleNamespace(**dict(zip(FIELDS,objects[:9])),fert2_probs=P._fert2_probs.copy(),bp_pol_stay=P._bp_pol_stay.copy(),c_pol_stay=P._c_pol_stay.copy());elapsed=time.monotonic()-before
  d.require(hashlib.sha256(baseline_V.tobytes()).hexdigest()==original_next_hash,'Baseline continuation array modified')
  record=dict(label=label,current_financed_share=phi,seconds=elapsed,future_continuation='exact saved baseline V; native DEAD regions and transition/child-aging unchanged',native_function_source=d.sha(fn.__code__.co_filename))
  if label=='baseline_replay':
   errors={}
   for key in ('V','fert_value','fert_probs'):
    a=getattr(sol,key);b=getattr(base,key);mask=occup if key=='V' else (pre[:,:,:,:,:,:,0].sum(-1)>0) if key=='fert_value' else (pre[:,:,:,:,:,:,0]>0)
    errors[key]=float(np.max(np.abs(a-b)))
   d.require(max(errors.values())<2e-9,'Baseline continuation replay differs from certified saved policy')
   record['full_array_max_errors']=errors
  else:saved[label]=sol
  attempts.append(record);d.write(HERE/'latest.json',dict(status='policy_completed',attempts=attempts))
 # Saved-policy algebra only, no new forward lifecycle or equilibrium root.
 result={name:invert(sol,P0,pi,w)for name,sol in saved.items()};common=np.logical_and.reduce([x['mask']for x in result.values()]);ww=w*common;reference=result['baseline'];rows=[];sizes=[];ages=[]
 for name,x in result.items():
  if name=='baseline':continue
  db=pi*(x['p']-reference['p']);dg=x['G']-reference['G'];sens=pi*reference['p']*(1-reference['p'])/P0.kappa_fert
  rows.append(dict(arm=name,wait_gain=mean(x['W']-reference['W'],ww),attempt_gain=mean(x['A']-reference['A'],ww),success_net_cost_gain=mean(x['C']-reference['C'],ww),attempt_gap_change=mean(dg,ww),birth_change_mass=float(np.sum(w*db)),positive_birth_change_mass=float(np.sum(w*np.maximum(db,0))),negative_birth_change_mass=float(np.sum(w*np.minimum(db,0))),birth_change_per_young_n0_renter=float(np.sum(w*db)/w.sum()),sensitivity_weighted_gap_change=mean(dg,ww*sens),excluded_birth_change_mass=float(np.sum(w[~common]*db[~common]))))
  for branch in ('wait_size','child_size'):
   for hi,H in enumerate(P0.H_own):sizes.append(dict(arm=name,branch=branch,owned_rooms=float(H),baseline_conditional_probability=mean(reference[branch][...,hi],ww),new_conditional_probability=mean(x[branch][...,hi],ww),probability_change=mean(x[branch][...,hi]-reference[branch][...,hi],ww)))
  for j in range(7):ages.append(dict(arm=name,age_cell_start=float(P0.age_start+j*P0.da),initial_mass=float(w[:,j].sum()),birth_change_mass=float(np.sum(w[:,j]*db[:,j]))))
 d.table(HERE/'value_comparison.csv',rows);d.table(HERE/'ownership_by_size.csv',sizes);d.table(HERE/'birth_by_age.csv',ages)
 d.write(HERE/'summary.json',dict(status='two_native_policy_evaluations_completed',policy_evaluations=attempts,lifecycle_evaluations=0,equilibrium_evaluations=0,total_seconds=time.time()-start,young_n0_renter_mass=float(w.sum()),common_interior_mass=float(ww.sum()),excluded_mass=float(w[~common].sum()),continuation_V_sha256=original_next_hash,continuation_dead_node_count=int((baseline_V<=-1e9).sum()),original_household_sha256=d.sha(h.__file__),patched_household_sha256=d.sha(PATCH),baseline_solution_sha256=d.sha(BASE/'solution_arrays.npz'),permanent_solution_sha256=d.sha(PERM/'solution_arrays.npz'),driver_sha256=d.sha(__file__),price=q,model_period_years=float(P0.da),results=rows,interpretation='Four-year current-decision credit relaxation with future values anticipated under baseline financing; permanent minus temporary includes duration/expectation interactions, not additive primitive decomposition; no dead continuation value projected or relaxed',remaining_feasibility_limitation='Uses native expected continuation and interpolation unchanged; no additional strict support operator or artificial next-period bailout. Occupied-state lifecycle/estate ledger not independently regenerated for this policy-only diagnostic.'))
 print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
