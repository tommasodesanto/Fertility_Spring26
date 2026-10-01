#!/usr/bin/env python3
"""Four native backward-policy calls for permanent psi sensitivity; no KFE/GE."""
import os,sys,copy,time,json,signal,hashlib,importlib.util
from pathlib import Path
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).resolve().parent;WIN=HERE.parents[3];sys.path.insert(0,str(WIN))
import fixed_price_responses as d
import numpy as np
BASE=WIN/'purchase_ltv_v1/local_run/retry5/results/baseline_80_80';CREDIT=WIN/'purchase_ltv_v1/local_run/retry9/results/both_100_100';PATCH=WIN/'purchase_ltv_v1/local_run/household_buyer_estate.py'
def alarm(*_):raise TimeoutError('120-second global native policy-call budget')
def mean(a,w):return float(np.sum(np.nan_to_num(a)*w)/w.sum())if w.sum()>0 else None
def invert(fp,fv,P,pi,w):
 p=fp[:,0,0,:,:,:2];F=fv[:,0,0,:,:];v=(w>0)&(p[...,0]>0)&(p[...,1]>0)&(pi>0)
 W=np.full_like(F,np.nan);A=W.copy();G=W.copy();C=W.copy();W[v]=F[v]+P.kappa_fert*np.log(p[...,0][v]);A[v]=F[v]+P.kappa_fert*np.log(p[...,1][v]);G[v]=A[v]-W[v];C[v]=W[v]+G[v]/np.broadcast_to(pi,w.shape)[v]
 return dict(p=p[...,1],W=W,A=A,G=G,C=C,v=v)
def main():
 started=time.time();auth=d.authenticate_candidate(HERE/'runtime_preparation');P0=auth['P'];grid=auth['grid'];q=float(d.BINDING['candidate_price'])
 from refactor_lab.engine import household as h
 from refactor_lab.engine.parameters import get_fecundity_by_age
 d.require(d.sha(h.__file__)==json.loads((BASE.parent/'source_receipt.json').read_text())['original_household_sha256'],'Native source changed');d.require(d.sha(PATCH)==json.loads((CREDIT.parent/'source_receipt.json').read_text())['override_sha256'],'Estate-patched source changed')
 spec=importlib.util.spec_from_file_location('refactor_lab.engine.permanent_psi_credit_household',PATCH);patch=importlib.util.module_from_spec(spec);sys.modules[spec.name]=patch;spec.loader.exec_module(patch)
 pre=np.load(BASE.parent/'q0_reference_inherited_states.npz')['g_pre'];w=pre[:,0,0,:,:,0,0].copy();w[:,7:]=0;pi=get_fecundity_by_age(P0)[None,:,None]
 data={};attempts=[];deadline=time.time()+120
 for psi in (.10,.25):
  for phi,fn in [(.8,h.solve_bellman_full_markov_income),(1.,patch.solve_bellman_full_markov_income)]:
   label=f'psi{psi:.2f}_phi{phi:.1f}';P=copy.deepcopy(P0);P.psi_child=psi;P.phi=np.full_like(np.asarray(P0.phi),phi)
   if hasattr(P,'_purchase_ltv_override'):del P._purchase_ltv_override
   sd=auth['solver'].precompute_shared(P,grid);actual=auth['context']['fp'].actual_parameters(auth['context']['prepared'],P,grid)
   expected=dict(auth['actual_parameters'],psi_child=psi,child_benefit_CRRA_coefficient=(1-P.child_benefit_curvature)*psi,financed_share=phi);d.require(actual==expected,'Unapproved native scalar parameter changed')
   d.write(HERE/'latest.json',dict(status='policy_call_claimed',label=label,policy_calls=len(attempts)+1,deadline_epoch=deadline))
   old=signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.001,deadline-time.time()));t=time.monotonic()
   try:objects=fn(np.asarray([P.user_cost_rate*q]),np.asarray([q]),P,grid,sd)
   finally:signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
   d.require(time.time()<deadline,'Global120-second budget exceeded');data[psi,phi]=invert(objects[7],objects[8],P,pi,w);attempts.append(dict(label=label,seconds=time.monotonic()-t,continuation_V='omitted: fully reoptimized native lifetime recursion',actual31_parameters=actual));del objects
   d.write(HERE/'latest.json',dict(status='policy_call_completed',attempts=attempts))
 # Reference psi uses certified existing arrays, without another policy call.
 original=float(P0.psi_child)
 for phi,folder in [(.8,BASE),(1.,CREDIT)]:
  with np.load(folder/'solution_arrays.npz')as a:data[original,phi]=invert(a['fert_probs'],a['fert_value'],P0,pi,w)
 rows=[];comparisons=[]
 for psi in (.10,original,.25):
  base=data[psi,.8];credit=data[psi,1.];common=base['v']&credit['v'];ww=w*common;delta=pi*(credit['p']-base['p']);sens=pi*base['p']*(1-base['p'])/P0.kappa_fert
  comparisons.append(dict(psi_child=psi,baseline_first_birth_probability=float(np.sum(w*pi*base['p'])/w.sum()),credit_first_birth_probability=float(np.sum(w*pi*credit['p'])/w.sum()),credit_effect_pp=float(100*np.sum(w*delta)/w.sum()),positive_birth_change_mass=float(np.sum(w*np.maximum(delta,0))),negative_birth_change_mass=float(np.sum(w*np.minimum(delta,0))),birth_change_mass=float(np.sum(w*delta)),wait_value_credit_gain=mean(credit['W']-base['W'],ww),success_net_cost_credit_gain=mean(credit['C']-base['C'],ww),attempt_gap_credit_change=mean(credit['G']-base['G'],ww),sensitivity_weighted_gap_credit_change=mean(credit['G']-base['G'],ww*sens),common_interior_mass=float(ww.sum()),excluded_birth_change_mass=float(np.sum(w[~common]*delta[~common])),source='reused certified saved policies'if psi==original else 'new native fully reoptimized policies'))
  for phi,x in [(.8,base),(1.,credit)]:
   rows.append(dict(psi_child=psi,financed_share=phi,first_birth_probability=float(np.sum(w*pi*x['p'])/w.sum()),mean_wait_value=mean(x['W'],ww),mean_success_net_cost_value=mean(x['C'],ww),mean_attempt_gap=mean(x['G'],ww),credit_effect_pp=0. if phi==.8 else comparisons[-1]['credit_effect_pp'],common_interior_mass=float(ww.sum()),source=comparisons[-1]['source']))
 d.table(HERE/'four_policy_calls.csv',[x for x in rows if x['psi_child'] in (.10,.25)]);d.table(HERE/'comparison_with_saved_reference.csv',comparisons)
 d.write(HERE/'summary.json',dict(status='four_policy_calls_completed',policy_calls=attempts,policy_count=4,KFE_calls=0,GE_calls=0,total_seconds=time.time()-started,solve_budget_seconds=120,price=q,original_psi_child=original,common_PRE_sha256=hashlib.sha256(pre.tobytes()).hexdigest(),young_n0_renter_mass=float(w.sum()),original_household_sha256=d.sha(h.__file__),patched_household_sha256=d.sha(PATCH),baseline_solution_sha256=d.sha(BASE/'solution_arrays.npz'),credit_solution_sha256=d.sha(CREDIT/'solution_arrays.npz'),driver_sha256=d.sha(__file__),comparisons=comparisons,changed_scalar_fields_baseline=['psi_child','child_benefit_CRRA_coefficient'],changed_scalar_fields_credit=['psi_child','child_benefit_CRRA_coefficient','financed_share'],interpretation='Permanent preference sensitivity at prescribed q/common PRE with fully reoptimized future household policies; not recalibration, new GE, stationary population comparison or global sign frontier',limit='No new forward lifecycle/estate ledger or17 diagnostic plots regenerated; unchanged native budget/death constraints retained; source identity and31effective scalar parameters pinned.'))
 print(json.dumps(comparisons,indent=2))
if __name__=='__main__':main()
