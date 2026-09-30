#!/usr/bin/env python3
"""Two fresh-process fixed-price household cases; Torch only, never GE."""
from __future__ import annotations
import argparse,copy,csv,hashlib,json,os,subprocess,sys,time,traceback
from pathlib import Path
import numpy as np
LABEL='2007 stationary reference — block0506, September 28 verified export'
ROOT=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
BASE=ROOT/'output/model/fixed_reference_economics_20260928/sources/fixed_price_v1/run_fixed_price.py'
CASES=('ours','author')
BASE_SHA='96d6923a252f57bc4d8c44fd6479b13f48ba217d74edf8ef629d120428b03b44'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x): Path(p).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def require(x,m):
 if not x: raise RuntimeError(m)
def load(p,n):
 import importlib.util
 s=importlib.util.spec_from_file_location(n,p); m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m
def weights(P):
 a=np.zeros(P.J+1); surv=np.asarray(P.survival_probs)
 for j in range(P.J): a[j+1]=0. if (j==P.J-1 or (P.use_age_survival and surv[j]<1.)) else 1.
 return a
def summary(ev,P,grid,base):
 x=base.aggregates(ev,P,grid)
 # This is current wealth mass, not a realized chosen-saving statistic.
 g=ev.g_current; renter=g[:,0].sum(axis=(1,2,3,4,5)); x['negative_current_renter_wealth_mass']=float(renter[grid<0].sum())
 return x
def apply_rule(P,case):
 P.lambda_d=0.; P.debt_taper_weights=np.zeros(P.J+1) if case=='author' else weights(P); P.debt_caps=np.zeros(P.J+1)
 return {'lambda_d':0.,'debt_taper_weights':P.debt_taper_weights.tolist(),'debt_caps':P.debt_caps.tolist(),'strict_sale_gate':case=='author'}
def child(a):
 out=a.output;out.mkdir(parents=True,exist_ok=False); require(sha(BASE)==BASE_SHA,'original helper pin differs'); base=load(BASE,'_fixed_base'); manifest,contract,obj,runtime,prepared,ref=base.authenticate(out)
 model=prepared.rt['model']; cal=prepared.rt['primitive'].pf.calendar; P=copy.deepcopy(ref['parameters']); grid=np.asarray(ref['b_grid']).copy(); before=base.serialized(vars(P))
 override=apply_rule(P,a.case)
 strict=None
 if a.case=='author':
  sp=Path(__file__).with_name('strict_tenure.py'); require(sha(sp)=='c64f98bf72484facc36be5eb4d38d581c1f813b419feae9ab28eee3e96117677','strict module pin differs'); strict=load(sp,'credit_rule_quick_strict_tenure'); model.tenure_logit_kernel=strict.tenure_logit_kernel;model.tenure_choice_kernel=strict.tenure_choice_kernel
 after=base.serialized(vars(P)); changed={k:{'reference':before.get(k),'experiment':after.get(k)} for k in set(before)|set(after) if before.get(k)!=after.get(k)}
 require(set(changed)<= {'lambda_d','debt_taper_weights','debt_caps'},'override changed undeclared primitives')
 require(np.all(P.debt_taper_weights==0) if a.case=='author' else np.array_equal(P.debt_taper_weights,weights(P)),'arrays rebuilt')
 price=np.asarray(ref['solution'].p_eq); sd=model.precompute_shared(P,grid)
 t=time.monotonic(); sol=model.solve_markov_income_at_prices(price,P,grid,SD=sd,verbose=False,fast_stats=False); elapsed=time.monotonic()-t; P._fert2_probs=sol.fert2_probs.copy(); pol=cal.policy_from_solution(sol,price,P,grid,sd); pre,recon=cal.reconstruct_stationary_pre_fertility(sol,pol,P,grid,sd)
 supply=cal.HousingSupplyRule('static-elastic',float(price[0]),float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0])); ev=cal.evaluate_period(price,pre,P,grid,sd,cal.SolveCounter(),supply_rule=supply,supplied_policy=pol)
 core=summary(ev,P,grid,base); completed=float(prepared.rt['chain'].extract_moments(sol,P)['tfr']); write(out/'preliminary_numbers.json',{'case':a.case,'reference_label':LABEL,'lifecycle_solve_seconds':elapsed,**core,'completed_fertility':completed,'actual_renewal_B':float(sol.adult_entry_adjusted_birth_children),'actual_renewal_E':float(sol.entry_rate),'q0':float(price[0]),'certified':False,'note':'Core fixed-price cohort numbers written before validation; no GE or demographic steady state.'})
 impact=cal.evaluate_period(price,ref['stationary_g_pre'],P,grid,sd,cal.SolveCounter(),supply_rule=supply,supplied_policy=pol); write(out/'baseline_g_pre_impact.json',summary(impact,P,grid,base))
 packet={'parameters':P,'b_grid':grid,'shared':sd,'solution':sol,'evaluation':ev,'stationary_g_pre':pre,'supply_rule':supply}
 failure=None
 try: gates=base.gates(packet,prepared,out,stationary=True)
 except Exception as e: failure={'type':type(e).__name__,'message':str(e),'traceback':traceback.format_exc()}; write(out/'validation_failure.json',failure); gates=None
 fertility={p:prepared.rt['observe_initial_fertility'](ev,P,age_projection=p) for p in ('uniform_birth_time','constant_post_cell')}; housing=prepared.rt['observe_initial_housing_wealth'](ev,P,grid,sd,diagnostic_enabled=True,age_projection='uniform_within_age_cell',diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True); recent=prepared.rt['observe_recent_parent_flow'](ev,P,diagnostic_enabled=True,snapshot=prepared.rt['SNAPSHOT'],age_projection=prepared.rt['AGE_PROJECTION'],diagnostic_allow_residence_proxy=True,input_provenance={'case':a.case}); fits=runtime.score_targets(obj,fertility,housing,recent['model_value'],float(prepared.rt['chain'].extract_moments(sol,P)['tfr']))
 for name,rows in [('target_fit.csv',fits),('parameters.csv',manifest['full_parameter_table'])]:
  with open(out/name,'w',newline='') as f: w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 plot_failure=None
 try:
  prepared.rt['audit'].standard_diagnostics(packet,out,validate_production_young=False)
  require(len(list((out/'standard_diagnostics').glob('*.png')))==17,'standard plot count differs')
 except Exception as e:
  plot_failure={'type':type(e).__name__,'message':str(e)}; write(out/'plot_failure.json',plot_failure)
 write(out/'receipt.json',{'status':'preliminary','case':a.case,'reference_label':LABEL,'lifecycle_solves':1,'lifecycle_solve_seconds':elapsed,'array_override':override,'strict_module':str(Path(__file__).with_name('strict_tenure.py')) if strict else None,'validation_failure':failure,'plot_failure':plot_failure,'gates_written':gates is not None,'not_GE':True})
def main():
 p=argparse.ArgumentParser();p.add_argument('--case',choices=CASES);p.add_argument('--output',type=Path,required=True);a=p.parse_args();require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm required');child(a)
if __name__=='__main__': main()
