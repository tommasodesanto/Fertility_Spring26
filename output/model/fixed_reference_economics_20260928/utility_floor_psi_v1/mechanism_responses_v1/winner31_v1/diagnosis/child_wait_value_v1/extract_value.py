#!/usr/bin/env python3
"""Exact saved first-birth logit inversion; no lifecycle/Bellman/GE calls."""
import os,sys,json,hashlib,csv
from pathlib import Path
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).resolve().parent;WIN=HERE.parents[1];sys.path.insert(0,str(WIN))
import fixed_price_responses as d
import numpy as np
ARMS={'baseline80_80':WIN/'purchase_ltv_v1/local_run/retry5/results/baseline_80_80','purchase100_stayer80':WIN/'purchase_ltv_v1/local_run/retry8/results/purchase_100_stayer_80','both100_100':WIN/'purchase_ltv_v1/local_run/retry9/results/both_100_100'}
def mean(x,w):return float(np.sum(x*w)/np.sum(w)) if np.sum(w)>0 else None
def quant(x,w,q):
 v=(w>0)&np.isfinite(x);ix=np.argsort(x[v]);xx=x[v][ix];ww=w[v][ix];return float(xx[np.searchsorted(np.cumsum(ww),q*ww.sum())])
def stats(x,w):return dict(mean=mean(x,w),p10=quant(x,w,.1),p50=quant(x,w,.5),p90=quant(x,w,.9))
def main():
 auth=d.authenticate_candidate(HERE/'runtime_preparation');P=auth['P'];grid=auth['grid']
 from refactor_lab.engine.parameters import get_fecundity_by_age,parent_age_maturation_active,independent_child_maturation_active
 d.require(P.sequential_births and not P.joint_nested_choice and P.I==1,'First-birth two-choice contract differs')
 pi=get_fecundity_by_age(P)[None,:,None];k=float(P.kappa_fert)
 pre=np.load(ARMS['baseline80_80'].parent/'q0_reference_inherited_states.npz')['g_pre'];w=pre[:,0,0,:,:,0,0].copy();ages=P.age_start+np.arange(P.J)*P.da;w[:,ages>42]=0;w[:,ages<18]=0
 fertile=((np.arange(P.J)+1>=P.A_f_start)&(np.arange(P.J)+1<=P.A_f_end))[None,:,None]
 source={};data={};pins={}
 hs=d.ROOT/'code/model/refactor_lab/engine/household.py';bs=WIN/'purchase_ltv_v1/local_run/household_buyer_estate.py'
 htext=hs.read_text();btext=bs.read_text();left='            if in_fert:\n                pi_j = float(fec[j])';right='                    # Entry (childless wait/try)'
 block=htext[htext.index(left):htext.index(right,htext.index(left))];blockb=btext[btext.index(left):btext.index(right,btext.index(left))];d.require(block==blockb,'Patched fertility block changed')
 for arm,folder in ARMS.items():
  cl=json.loads((folder/'closure.json').read_text());sr=json.loads((folder.parent/'source_receipt.json').read_text());d.require(d.sha(hs)==sr['original_household_sha256'],'Baseline household source drift')
  if arm!='baseline80_80':d.require(d.sha(bs)==sr['override_sha256'],'Buyer/death-floor override source drift')
  d.require(cl['price']==float(d.BINDING['candidate_price']),'q differs');d.require(cl['baseline_state_impact']['inherited_distribution_sha256']==hashlib.sha256(pre.tobytes()).hexdigest(),'PRE differs')
  vals={x['parameter']:float(x['estimate'])for x in csv.DictReader((folder/'parameters.csv').open())};expected=dict(auth['actual_parameters']);expected['financed_share']=1. if arm=='both100_100' else .8;d.require(vals==expected,'Other31scalar changes')
  with np.load(folder/'solution_arrays.npz',allow_pickle=False) as a:
   probs=a['fert_probs'][:,0,0,:,:,:2].copy();F=a['fert_value'][:,0,0,:,:].copy();hr=a['hR_pol'][:,0,0,:,:,0,0].copy();tp=a['tenure_probs'][:,0,0,:,:,:,:,:]
   ow=tp[:,:, :,0,0,1:].sum(-1);oc=tp[:,:, :,1,1,1:].sum(-1)
   lp=a['loc_probs'][:,0,0,0,:,:,:,:];d.require(np.max(abs(lp[:,:, :,0,0][(w>0)&fertile]-1))<1e-6,'Occupied single-location wait probability differs')
  p0=probs[...,0];p1=probs[...,1];interior=(w>0)&fertile&(p0>0)&(p1>0)&(pi>0)
  d.require(np.max(abs(lp[:,:,:,1,1][interior]-1))<1e-6,'Occupied single-location successful-child probability differs')
  W=np.full_like(F,np.nan);A=W.copy();G=W.copy();C=W.copy();W[interior]=F[interior]+k*np.log(p0[interior]);A[interior]=F[interior]+k*np.log(p1[interior]);G[interior]=A[interior]-W[interior];C[interior]=W[interior]+G[interior]/np.broadcast_to(pi,G.shape)[interior]
  d.require(np.max(abs((p0+p1)[interior]-1))<2e-12,'Interior probabilities do not sum1')
  expected_owner=(1-pi*p1)*ow+pi*p1*oc
  data[arm]=dict(p0=p0,p1=p1,F=F,W=W,A=A,G=G,C=C,interior=interior,own_wait=ow,own_child=oc,own=expected_owner,hr=hr)
  pins[arm]=dict(folder=str(folder),source_receipt=sr,solution_arrays_sha256=d.sha(folder/'solution_arrays.npz'),closure_sha256=d.sha(folder/'closure.json'),parameters_sha256=d.sha(folder/'parameters.csv'),interior_mass=float(w[interior].sum()),nonfertile_mass=float(w[(w>0)&~np.broadcast_to(fertile,w.shape)].sum()),fertile_noninterior_mass=float(w[(w>0)&fertile&~interior].sum()))
 common=np.logical_and.reduce([v['interior']for v in data.values()]);ww=np.where(common,w,0);base=data['baseline80_80'];cap=base['hr']>=float(P.hR_max)-1e-9
 summary=[];byage=[];bysplit=[];dist=[]
 for arm,x in data.items():
  for name in ('p1','W','A','G','C','own_wait','own_child'):dist.append(dict(arm=arm,object=name,**stats(np.nan_to_num(x[name]),ww)))
  if arm=='baseline80_80':continue
  db=pi*(x['p1']-base['p1']);dg=x['G']-base['G'];sens=pi*base['p1']*base['p0']/k;do=x['own']-base['own'];positive=(do>0)
  summary.append(dict(arm=arm,interior_mass=float(ww.sum()),wait_gain=mean(np.nan_to_num(x['W']-base['W']),ww),attempt_gain=mean(np.nan_to_num(x['A']-base['A']),ww),success_net_cost_gain=mean(np.nan_to_num(x['C']-base['C']),ww),attempt_gap_change=mean(np.nan_to_num(dg),ww),first_birth_change_mass=float(np.sum(w*db)),positive_birth_change_mass=float(np.sum(w*np.maximum(db,0))),negative_birth_change_mass=float(np.sum(w*np.minimum(db,0))),birth_change_per_young_renter=float(np.sum(w*db)/w.sum()),birth_linear_gap_approximation_mass=float(np.sum(ww*sens*np.nan_to_num(dg))),sensitivity_weighted_gap_change=mean(np.nan_to_num(dg),ww*sens),own_wait_change=mean(x['own_wait']-base['own_wait'],ww),own_child_change=mean(x['own_child']-base['own_child'],ww),net_expected_owner_mass_change=float(np.sum(w*do)),positive_owner_change_weighted_gap_change=mean(np.nan_to_num(dg),ww*np.maximum(do,0)),signed_owner_change_weighted_gap_change=mean(np.nan_to_num(dg),ww*do)))
  for j,age in enumerate(ages):
   if age>42:continue
   agew=np.zeros_like(w);agew[:,j]=w[:,j];iv=agew*common
   byage.append(dict(arm=arm,age_cell_start=float(age),initial_mass=float(agew.sum()),first_birth_change_mass=float(np.sum(agew*db)),positive_birth_change_mass=float(np.sum(agew*np.maximum(db,0))),negative_birth_change_mass=float(np.sum(agew*np.minimum(db,0))),wait_gain=mean(np.nan_to_num(x['W']-base['W']),iv),success_net_cost_gain=mean(np.nan_to_num(x['C']-base['C']),iv),attempt_gap_change=mean(np.nan_to_num(dg),iv),net_expected_owner_mass_change=float(np.sum(agew*do))))
  splits={'baseline_wait_renter_at_cap':cap,'baseline_wait_renter_below_cap':~cap,'baseline_try_p_lt05':base['p1']<.05,'baseline_try_p_05_25':(base['p1']>=.05)&(base['p1']<.25),'baseline_try_p_25_75':(base['p1']>=.25)&(base['p1']<.75),'baseline_try_p_ge75':base['p1']>=.75,'expected_ownership_increases':do>0,'expected_ownership_decreases':do<0}
  for name,mask in splits.items():
   sw=np.where(mask,w,0);iv=sw*common
   bysplit.append(dict(arm=arm,split=name,initial_mass=float(sw.sum()),interior_mass=float(iv.sum()),baseline_try_probability=mean(base['p1'],sw),first_birth_change_mass=float(np.sum(sw*db)),positive_birth_change_mass=float(np.sum(sw*np.maximum(db,0))),negative_birth_change_mass=float(np.sum(sw*np.minimum(db,0))),attempt_gap_change=mean(np.nan_to_num(dg),iv),wait_gain=mean(np.nan_to_num(x['W']-base['W']),iv),success_net_cost_gain=mean(np.nan_to_num(x['C']-base['C']),iv),net_expected_owner_mass_change=float(np.sum(sw*do))))
 excluded=(w>0)&~common
 checks=dict(excluded_group_mass=float(w[excluded].sum()),excluded_group_fraction=float(w[excluded].sum()/w.sum()),excluded_group_birth_change_mass={arm:float(np.sum(w[excluded]*np.broadcast_to(pi,w.shape)[excluded]*(x['p1']-base['p1'])[excluded])) for arm,x in data.items()},child_location_assertion='Each occupied interior child-success single-location probability checked within1e-6; same check as wait',no_clipping_or_population_renormalization=True)
 d.table(HERE/'weighted_distributions.csv',dist);d.table(HERE/'by_age.csv',byage);d.table(HERE/'by_baseline_pressure.csv',bysplit)
 d.write(HERE/'summary.json',dict(small_invariant_checks=checks,status='saved_value_inversion_complete_no_solves',summaries=summary,source_pins=pins,household_sha256=d.sha(hs),buyer_override_sha256=d.sha(bs),first_birth_block_sha256=hashlib.sha256(block.encode()).hexdigest(),kappa_fert=k,first_birth_fixed_cost=P.first_birth_fixed_cost,age_cells=ages.tolist(),fecundity_by_age=pi.ravel().tolist(),fertile_age_indices_1based=[P.A_f_start,P.A_f_end],young_n0_renter_PRE_mass=float(w.sum()),common_interior_mass=float(ww.sum()),newborn_exempt_success_value_active=bool(parent_age_maturation_active(P)and independent_child_maturation_active(P)),rental_cap=float(P.hR_max),conditional_ownership_scope='Realized saved post-fertility branch policies, not pre-constraint latent willingness; if newborn-exempt optimization active, success-value VI_ex policies are not separately saved',purchase_feasibility_split='Not measured: no pre-mask value/desired size object; rental-cap split is conditional WAIT renter hR policy, not realized rental status',lifecycle_calls=0,Bellman_calls=0,GE_calls=0,driver_sha256=d.sha(__file__)))
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
