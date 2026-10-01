#!/usr/bin/env python3
"""Saved unnormalized-logit removal of current two-room ownership only."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;WIN=HERE.parents[2];BASE=WIN/'purchase_ltv_v1/local_run/retry5/results/baseline_80_80';PERM=WIN/'purchase_ltv_v1/local_run/retry9/results/both_100_100'
def main():
 meta=json.loads((HERE.parent/'summary.json').read_text());params={x['parameter']:float(x['estimate'])for x in csv.DictReader((BASE/'parameters.csv').open())};k=params['kappa_fert'];kt=params['tenure_choice_kappa'];pi=np.asarray(meta['fecundity_by_age'])[None,:,None]
 w=np.load(BASE.parent/'q0_reference_inherited_states.npz')['g_pre'][:,0,0,:,:,0,0].copy();w[:,7:]=0;arm={};records=[]
 for name,folder in [('baseline',BASE),('permanent_both100',PERM)]:
  with np.load(folder/'solution_arrays.npz')as a:
   p=a['fert_probs'][:,0,0,:,:,:2];p2=a['tenure_probs'][:,0,0,:,:,0,0,1].astype(float);child2=a['tenure_probs'][:,0,0,:,:,1,1,1];F=a['fert_value'][:,0,0,:,:]
  occupied=w>0;assert np.max(p2[occupied])<1 and np.min(p2[occupied])>=0;assert np.max(abs(child2[occupied]))==0;assert np.max(abs(p.sum(-1)[occupied]-1))<2e-12
  dW=kt*np.log1p(-p2);shift=-pi*dW/k;mult=np.exp(shift);den=p[...,0]+p[...,1]*mult;assert np.all(den[occupied]>0);pnew=np.divide(p[...,1]*mult,den,out=np.zeros_like(den),where=den>0)
  original=p[...,1];iv=occupied&(p[...,0]>0)&(p[...,1]>0);records.append(dict(arm=name,p2_max=float(p2[occupied].max()),rounded_p2_one_mass=float(w[occupied&(p2>=1)].sum()),attempt_zero_mass=float(w[occupied&(p[...,1]==0)].sum()),mean_wait_value_change_interior=float(np.sum(w[iv]*dW[iv])/w[iv].sum()),mean_attempt_gap_change_from_ablation_interior=float(np.sum(w[iv]*(-np.broadcast_to(pi,w.shape)[iv]*dW[iv]))/w[iv].sum()),baseline_population_first_birth_mass=float(np.sum(w*pi*original)),ablated_first_birth_mass=float(np.sum(w*pi*pnew)),largest_odds_multiplier=float(mult[occupied].max()),zero_probability_handling='Saved zero probabilities remain zero under the finite odds multiplier; underlying underflow uncertainty is unchanged and negligible, branch levels not invented'))
  arm[name]=dict(original=original,pnew=pnew)
 raw=pi*(arm['permanent_both100']['original']-arm['baseline']['original']);new=pi*(arm['permanent_both100']['pnew']-arm['baseline']['pnew'])
 result=dict(status='zero_solve_current_option_ablation',formulas='DeltaW=kappa_tenure*log(1-p_owned2_wait); C unchanged because child-owned2 prob exactly0; Delta(A-W)=-pi*DeltaW; new_p_try=p_try*exp(-pi*DeltaW/kappa_fert)/(p_wait+p_try*exp(-pi*DeltaW/kappa_fert))',current_menu_intervention='Remove two-room ownership from CURRENT decision only; retained original future value functions and all other choices/preferences/financing',population_mass=float(w.sum()),arms=records,original_permanent_minus_baseline_birth_change_mass=float(np.sum(w*raw)),ablated_permanent_minus_baseline_birth_change_mass=float(np.sum(w*new)),original_probability_change=float(np.sum(w*raw)/w.sum()),ablated_probability_change=float(np.sum(w*new)/w.sum()),positive_ablated_change_mass=float(np.sum(w*np.maximum(new,0))),negative_ablated_change_mass=float(np.sum(w*np.minimum(new,0))),new_model_solves=0,numerical_caveat='Analytic unnormalized-logit identity evaluated from saved float32 tenure probabilities; exact theory, finite saved precision. No clipping, probability renormalization or population dropping.',interpretation='Diagnostic removal, not proposed adoption or full-menu recalibration; changing current choice availability while future menus remain original. Baseline and permanent-credit policy branch levels remain separately optimized.',driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (HERE/'two_room_ablation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2,allow_nan=False))
if __name__=='__main__':main()
