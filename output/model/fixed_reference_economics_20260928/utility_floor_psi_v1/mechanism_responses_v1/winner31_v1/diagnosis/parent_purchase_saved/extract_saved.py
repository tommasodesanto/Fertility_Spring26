#!/usr/bin/env python3
"""Saved-policy algebra only; no lifecycle/Bellman/equilibrium solve."""
import os,sys,json,copy,hashlib
from pathlib import Path
from types import SimpleNamespace
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).resolve().parent;WINNER=HERE.parents[1]
sys.path.insert(0,str(WINNER))
import fixed_price_responses as d
import numpy as np
sys.path.insert(0,str(d.ROOT/"code/model"))
from refactor_lab.engine.parameters import children_at_home_count
BASE=WINNER/'purchase_ltv_v1/local_run/retry5/results/baseline_80_80'
BUY=WINNER/'purchase_ltv_v1/local_run/retry8/results/purchase_100_stayer_80'
PRE=BASE.parent/'q0_reference_inherited_states.npz'
def main():
 auth=d.authenticate_candidate(HERE/'runtime_preparation_final');P=auth['P'];grid=auth['grid'];sd=auth['solver'].precompute_shared(P,grid)
 cal=auth['context']['prepared'].rt['primitive'].pf.calendar
 def forbidden(*a,**k):raise RuntimeError('Numerical solve forbidden: saved-policy extraction only')
 cal.solve_policy=forbidden
 with np.load(PRE,allow_pickle=False) as a:g=a['g_pre'].copy()
 ph=hashlib.sha256(g.tobytes()).hexdigest();ages=P.age_start+np.arange(P.J)*P.da
 young=(ages>=18)&(ages<=42)
 parent=np.asarray([[children_at_home_count(n,c,P)>0 for c in range(P.n_child_states)]for n in range(P.n_parity)])
 masks={}
 for name in ['initial_renter_young_n0','initial_renter_young_current_parent','initial_renter_young_all','entire_PRE']:
  mask=np.ones(g.shape,dtype=bool) if name=='entire_PRE' else np.zeros(g.shape,dtype=bool)
  if name!='entire_PRE':
   mask[:,0,:,young,:,:,:]=True
   if name.endswith('_n0'):mask[:,:,:,:,:,1:,:]=False
   if name.endswith('_current_parent'):mask &= parent[None,None,None,None,None,:,:]
  masks[name]=mask
 rows=[];source=[]
 for label,folder in [('baseline80_80',BASE),('purchase100_stayer80',BUY)]:
  closure=json.loads((folder/'closure.json').read_text());d.require(closure['price']==float(d.BINDING['candidate_price']),'q differs');d.require(closure['baseline_state_impact']['inherited_distribution_sha256']==ph,'common PRE differs')
  source.append(dict(case=label,path=str(folder),price=closure['price'],common_PRE_sha256=ph,receipt=json.loads((folder/'receipt.json').read_text())))
  needed=['V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value','fert2_probs','bp_pol_stay','c_pol_stay']
  with np.load(folder/'solution_arrays.npz',allow_pickle=False) as a:sol=SimpleNamespace(**{k:a[k].copy() for k in needed})
  PP=copy.deepcopy(P);PP._fert2_probs=sol.fert2_probs.copy()
  policy=cal.policy_from_solution(sol,np.array([closure['price']]),PP,grid,sd)
  for group,mask in masks.items():
   inherited=g*mask;counter=cal.SolveCounter();ev=cal.evaluate_period(np.array([closure['price']]),inherited,PP,grid,sd,counter,supplied_policy=policy)
   d.require(counter.bellman==0,'Unexpected Bellman solve')
   mass=float(inherited.sum());d.require(abs(float(ev.g_pre.sum())-mass)<2e-10 and abs(float(ev.g_current.sum())-mass)<2e-10,'Group mass/feasibility changed')
   owner=float(ev.g_current[:,1:].sum());row=dict(case=label,group=group,initial_PRE_mass=mass,realized_current_owner_mass=owner,ownership_rate=owner/mass,birth_children_mass=float(ev.births),birth_children_per_initial_household=float(ev.births)/mass,projection_mass=float(ev.feasibility_projection_mass),bellman_calls=counter.bellman)
   if group=='entire_PRE':
    d.require(abs(row['ownership_rate']-closure['baseline_state_impact']['ownership_rate'])<2e-10,'Replay ownership not saved impact')
    d.require(abs(row['birth_children_mass']-closure['baseline_state_impact']['births'])<2e-10,'Replay births not saved impact')
   rows.append(row)
 d.table(HERE/'paired_group_outcomes.csv',rows)
 deltas=[]
 for group in masks:
  b=next(x for x in rows if x['case']=='baseline80_80' and x['group']==group);a=next(x for x in rows if x['case']=='purchase100_stayer80' and x['group']==group)
  deltas.append(dict(group=group,initial_PRE_mass=b['initial_PRE_mass'],baseline_ownership_rate=b['ownership_rate'],purchase_ownership_rate=a['ownership_rate'],ownership_rate_change=a['ownership_rate']-b['ownership_rate'],new_owner_mass_net=a['realized_current_owner_mass']-b['realized_current_owner_mass'],baseline_birth_children_per_household=b['birth_children_per_initial_household'],purchase_birth_children_per_household=a['birth_children_per_initial_household'],birth_children_per_household_change=a['birth_children_per_initial_household']-b['birth_children_per_initial_household']))
 d.write(HERE/'summary.json',dict(status='saved_policy_algebra_completed_no_solves',source=source,common_PRE_sha256=ph,group_outcomes=deltas,age_definition='Model cell starts18,22,26,30,34,38,42; not annual age overlap or new target clock',state_timing='Inherited PRE-fertility renter tenure; fertility then location/tenure choices realized under respective saved policy',limitations=['No pre-mask buyer-versus-renter option values or shadow prices saved; zero choice probability does not prove infeasibility','Net ownership response is not fraction of individually blocked purchasers; no joint taste-shock coupling','Birth-child flows include continuation births for current parents; n0 group is first-birth mass','Purchase-only experiment changes lifetime policies and continuation values; fixed PRE is common but policies are optimized separately'],lifecycle_calls=0,Bellman_calls=0,GE_calls=0,driver_sha256=d.sha(__file__)))
 print(json.dumps(deltas,indent=2))
if __name__=='__main__':main()
