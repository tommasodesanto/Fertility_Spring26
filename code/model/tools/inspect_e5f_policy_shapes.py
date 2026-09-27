#!/usr/bin/env python3
"""Supplemental saved-policy slices; zero solves, no policy modifications."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from inspect_e5f_saved_households import load,weighted_quantile

def write_occupied_declines(rows, output):
 """Reduce saved slices; upper-node mass is descriptive, not causal incidence."""
 results=[]
 for age in (30,42):
  age_rows=[r for r in rows if int(r['age'])==age]
  total=sum(float(r['post_fertility_pre_tenure_mass']) for r in age_rows)
  for key in ('owner_probability','standard_modal_housing','expected_realized_housing','expected_realized_consumption','expected_next_assets'):
   drops=[]
   for z in sorted({int(r['income_index']) for r in age_rows}):
    line=sorted([r for r in age_rows if int(r['income_index'])==z],key=lambda r:float(r['wealth']))
    for left,right in zip(line,line[1:]):
     if not all(str(v['valid']).lower()=='true' and float(v['post_fertility_pre_tenure_mass'])>1e-10 for v in (left,right)):continue
     step=float(right[key])-float(left[key])
     if step < -1e-5:drops.append((step,z,float(left['wealth']),float(right['wealth']),float(right['post_fertility_pre_tenure_mass'])))
   worst=min(drops) if drops else (0,None,None,None,0)
   mass=sum(d[4] for d in drops)
   results.append(dict(age=age,object=key,adjacent_mass_threshold=1e-10,slice_mass=total,drop_pair_count=len(drops),upper_node_mass=mass,share_of_slice=mass/total if total else 0,largest_drop=worst[0],largest_drop_income_index=worst[1],wealth_left=worst[2],wealth_right=worst[3]))
 with (output/'occupied_declines.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
 return results

def main():
 a=argparse.ArgumentParser();a.add_argument('--contract',type=Path,required=True);a.add_argument('--case',type=Path,required=True);a.add_argument('--output',type=Path,required=True);x=a.parse_args()
 packet,rt,r=load(x.contract,x.case,x.output)
 P=packet['parameters'];e=packet['evaluation'];p=e.policy;b=packet['b_grid']; model=rt['model']; z,_,_=model.income_transition_values(P)
 assert P.I==1 and p.joint_choice is None
 assert np.isfinite(p.tenure_probs).all()
 assert p.tenure_probs.min()>=0 and p.tenure_probs.max()<=1
 assert np.max(abs(p.tenure_probs.sum(axis=-1)[e.g_post_fertility>1e-12]-1))<1e-6
 print('shapes',p.tenure_probs.shape,p.maps.tmx_idx.shape,'z',z,'childstate',p.V.shape[-1],flush=True)
 print('runtime',model.__file__,flush=True)
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 cs=0
 # Diagnostic convention is readiness_settled_state; authenticate actual helper.
 import importlib
 diag=importlib.import_module(model.__package__+'.diagnostics') if model.__package__ else None
 if diag is not None:cs=int(diag.readiness_settled_state(P))
 print('cs',cs,flush=True)
 rows=[]; findings=[]; tails=[]; allmass=e.g_pre.sum(); commonmass=e.g_post_fertility.sum(axis=(1,2,3,4,5,6)); commonlow,common99=weighted_quantile(b,commonmass,[.001,.99])
 current99=weighted_quantile(b,e.g_current.sum(axis=(1,2,3,4,5,6)),[.99])[0]
 for age in (30,42):
  j=int(round((age-P.age_start)/P.da));fig,axs=plt.subplots(2,3,figsize=(13,7))
  # Three separated actual income ranks keep the figure legible; CSV retains all15.
  selected=[0,len(z)//2,len(z)-1]
  for zz,zval in enumerate(z):
   ix=(slice(None),0,0,j,zz,0,cs); mass=e.g_post_fertility[ix];valid=p.V[ix]>-1e9
   own=p.tenure_probs[ix][:,1:].sum(axis=1);tc=p.tenure_choice[ix];hr=p.hR_pol[ix]; hh=np.where(tc>0,np.asarray(P.H_own)[np.maximum(tc-1,0)],hr)
   annual=float(model.annual_gross_income_at_state(P,0,j,zval))
   for k in [int(np.argmin(abs(b-1000))),len(b)-2,len(b)-1]:
    probs=p.tenure_probs[ix][k];tails.append(dict(age=age,income_index=zz,wealth=float(b[k]),probabilities=probs.tolist(),max_deviation_equal=float(np.max(abs(probs-1/6))),state_mass=float(mass[k])))
   expected=np.zeros(len(b));expectc=np.zeros(len(b));expectbp=np.zeros(len(b))
   for t in range(len(P.H_own)+1):
    kt=p.maps.tmx_idx[0,0,t,0,cs,:];wt=p.maps.tmx_wt[0,0,t,0,cs,:]
    dest=(slice(None),t,0,j,zz,0,cs)
    hhcond=p.hR_pol[dest] if t==0 else np.full(len(b),P.H_own[t-1])
    pr=p.tenure_probs[ix][:,t]
    expected+=pr*((1-wt)*hhcond[kt]+wt*hhcond[kt+1])
    expectc+=pr*((1-wt)*p.c_pol[dest][kt]+wt*p.c_pol[dest][kt+1])
    expectbp+=pr*((1-wt)*p.bp_pol[dest][kt]+wt*p.bp_pol[dest][kt+1])
   own_drop=np.r_[False,(np.diff(own)<-1e-5)&valid[1:]&valid[:-1]]
   h_drop=np.r_[False,(np.diff(hh)<-1e-5)&valid[1:]&valid[:-1]]
   ex_drop=np.r_[False,(np.diff(expected)<-1e-5)&valid[1:]&valid[:-1]]
   summaries={}
   for nm,flag in [('ownership_drop',own_drop),('standard_housing_drop',h_drop),('expected_housing_drop',ex_drop)]:
    summaries[nm]=dict(first_wealth=float(b[flag][0]) if flag.any() else None,minimum_adjacent_step=float(np.min(np.diff({'ownership_drop':own,'standard_housing_drop':hh,'expected_housing_drop':expected}[nm])[valid[1:]&valid[:-1]])),upper_node_mass=float(mass[flag].sum()),upper_node_population_share=float(mass[flag].sum()/allmass))
   findings.append(dict(age=age,income_index=zz,productivity=float(zval),annual_income=annual,slice_mass=float(mass.sum()),slice_wealth99=weighted_quantile(b,mass,[.99])[0] if mass.sum()>0 else None,**summaries))
   for k in range(len(b)):
    rows.append(dict(age=age,income_index=zz,productivity=float(zval),annual_income=annual,wealth=float(b[k]),valid=bool(valid[k]),post_fertility_pre_tenure_mass=float(mass[k]),conditional_rent_consumption=float(p.c_pol[ix][k]),conditional_rent_housing=float(hr[k]),modal_tenure=int(tc[k]),standard_modal_housing=float(hh[k]),owner_probability=float(own[k]),expected_realized_housing=float(expected[k]),expected_realized_consumption=float(expectc[k]),expected_next_assets=float(expectbp[k]),**{f'tenure_probability_{t}':float(p.tenure_probs[ix][k,t]) for t in range(len(P.H_own)+1)}))
   if zz in selected:
    label=f'z rank {zz+1}/{len(z)}; annual income {annual:.3f}'
    for ax,line,title in zip(axs.ravel(),[own,hh,expected,p.c_pol[ix],expectc,np.cumsum(mass)/mass.sum() if mass.sum()>0 else mass],['Owner probability, conditional family state','Standard modal-tenure housing','Expected housing, transaction maps applied','Consumption conditional on renting','Expected consumption, transaction maps applied','CDF within this state slice']):
     ax.plot(b,np.where(valid & (b>=commonlow) & (b<=common99),line,np.nan),label=label);ax.set_title(title,fontsize=10);ax.set_xlim(commonlow,common99);ax.grid(alpha=.2);ax.set_xlabel('Beginning financial wealth')
  axs[0,0].set_ylim(-.02,1.02);axs[1,2].set_ylim(-.02,1.02)
  handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
  fig.suptitle(f'Supplemental age {age}, childless incoming renter; wealth through population p99={common99:.3f}')
  fig.tight_layout(rect=[0,.06,1,.95]);fig.savefig(x.output/f'occupied_policy_age{age}.png',dpi=150);plt.close(fig)
 with (x.output/'policy_slices.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 write_occupied_declines(rows,x.output)
 (x.output/'shape_audit.json').write_text(json.dumps(dict(native_solves=0,child_state=cs,pre_tenure_wealth_p001=commonlow,common_wealth99=common99,post_tenure_wealth99=current99,findings=findings,tail_probabilities=tails,notes=['Mass is post-fertility and before tenure choice, so family state is conditioned on rather than integrating over births.','Decline masses refer to the upper grid node of a decreasing adjacent pair, not causal incidence.','Expected controls use the saved tenure lotteries and exact transaction interpolation maps; one location.','Plots select income ranks 1,8,15; CSV includes every income state.','No theorem requires housing or ownership to be increasing in wealth.']),indent=2)+'\n')
 print(json.dumps({'done':True,'rows':len(rows),'common99':common99}),flush=True)
if __name__=='__main__':main()
