#!/usr/bin/env python3
"""Preliminary two-rule closed stationary GE root controller (Torch only)."""
from __future__ import annotations
import argparse, copy, csv, hashlib, importlib.util, json, math, os, signal, subprocess, sys, time, traceback
from pathlib import Path
import numpy as np

ROOT=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
HERE=Path(__file__).resolve().parent; SIBLING=HERE.parent/'credit_rule_quick_v1'
BASE=ROOT/'output/model/fixed_reference_economics_20260928/sources/fixed_price_v1/run_fixed_price.py'
PARENT=HERE/'run_ge_accounting_source.py'
LABEL='2007 stationary reference — block0506, September 28 verified export'; RENEWAL=1e-6; PAYGO=1e-6
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):
 p=Path(p); t=p.with_suffix(p.suffix+'.tmp'); t.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n');t.replace(p)
def read(p): return json.loads(Path(p).read_text())
def req(x,m):
 if not x: raise RuntimeError(m)
def load(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def progress(out,phase,**kw): write(Path(out)/'progress.json',dict(phase=phase,time_epoch=time.time(),**kw))
def bracket(rows):
 r=sorted(rows,key=lambda x:x['factor'])
 return next(((a,b) for a,b in zip(r,r[1:]) if a['residual']*b['residual']<=0),None)
def logsecant(a,b):
 lo,hi=math.log(a['factor']),math.log(b['factor']); f,g=a['residual'],b['residual']; req(f*g<=0 and f!=g,'invalid bracket')
 x=(lo*g-hi*f)/(g-f); return math.exp(min(max(x,lo+.1*(hi-lo)),hi-.1*(hi-lo)))
def root_eligible(row): return bool(row['gate_passed']) and abs(row['residual'])<=RENEWAL
def verify_inputs(plan):
 s=plan['source_sha256']; req(sha(HERE/'rule_interface.py')==s['rule_interface.py'],'rule bridge identity differs')
 req(sha(SIBLING/'run_credit_rule_quick.py')==s['../credit_rule_quick_v1/run_credit_rule_quick.py'],'sibling rule driver identity differs')
 req(sha(SIBLING/'strict_tenure.py')==s['../credit_rule_quick_v1/strict_tenure.py'],'strict tenure identity differs')
 req(sha(PARENT)==s['run_ge_accounting_source.py'],'reviewed GE accounting identity differs')

def child(a):
 plan=read(a.plan); verify_inputs(plan); out=a.output; out.mkdir(parents=True,exist_ok=False); base=load(BASE,'_ge_base'); parent=load(PARENT,'_ge_parent'); api=load(HERE/'rule_interface.py','_ge_rule_api')
 manifest,contract,obj,runtime,prepared,ref=base.authenticate(out); model=prepared.rt['model'];cal=prepared.rt['primitive'].pf.calendar
 P=copy.deepcopy(ref['parameters']); grid=api.prepare(ref); req(len(grid)==160,'required original 160 grid absent')
 before=base.serialized(vars(P)); override=api.apply_rule(model,P,a.case); after=base.serialized(vars(P)); changed={k:1 for k in set(before)|set(after) if before.get(k)!=after.get(k)}
 req(set(changed) <= {'lambda_d','debt_taper_weights','debt_caps'},'debt-rule bridge mutated undeclared primitives')
 q0=float(np.asarray(ref['solution'].p_eq)[0]); price=np.asarray([q0*a.factor]); sd=model.precompute_shared(P,grid)
 progress(out,'one_lifecycle_solve',factor=a.factor,deadline=a.deadline); req(time.time()<a.deadline,'deadline before solve')
 t=time.monotonic(); sol=model.solve_markov_income_at_prices(price,P,grid,SD=sd,verbose=False,fast_stats=False); elapsed=time.monotonic()-t; req(elapsed<=180 and time.time()<a.deadline,'price solve exceeded 180 seconds/deadline')
 P._fert2_probs=sol.fert2_probs.copy(); pol=cal.policy_from_solution(sol,price,P,grid,sd); pre,recon=cal.reconstruct_stationary_pre_fertility(sol,pol,P,grid,sd)
 supply=cal.HousingSupplyRule('static-elastic',float(price[0]),float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0])); ev=cal.evaluate_period(price,pre,P,grid,sd,cal.SolveCounter(),supply_rule=supply,supplied_policy=pol)
 closure=parent.closed_accounting(float(sol.entry_rate),float(sol.adult_entry_adjusted_birth_children),float(np.asarray(ev.demand_by_loc).sum()),float(np.asarray(ev.supply_by_loc).sum()),float(price[0])); closure['supply_exponent']=float(P.xi_supply[0])
 core=base.aggregates(ev,P,grid); completed=float(prepared.rt['chain'].extract_moments(sol,P)['tfr'])
 write(out/'preliminary_numbers.json',dict(case=a.case,price=float(price[0]),rent=float(P.user_cost_rate*price[0]),lifecycle_solve_seconds=elapsed,renewal_residual=closure['renewal_residual'],population_scale=closure['population_scale'],actual_births=closure['absolute_adjusted_births'],actual_entry=closure['absolute_entry'],actual_renewal_B=float(sol.adult_entry_adjusted_birth_children),actual_renewal_E=float(sol.entry_rate),q0=q0,completed_fertility=completed,**core,certified=False,note='Written immediately after solve; gates/rendering pending.'))
 fiscal_failure=None; gate_failure=None; gates=None; fiscal=None
 try:
  packet=dict(parameters=P,b_grid=grid,shared=sd,solution=sol,evaluation=ev,stationary_g_pre=pre,supply_rule=supply); gates=base.gates(packet,prepared,out,stationary=True); fiscal=float(gates['fiscal']['scaled_pension_budget_residual']); req(abs(fiscal)<=PAYGO,'actual PAYGO failure')
 except Exception as e: gate_failure={'type':type(e).__name__,'message':str(e)}; fiscal_failure=gate_failure
 fertility={x:prepared.rt['observe_initial_fertility'](ev,P,age_projection=x) for x in ('uniform_birth_time','constant_post_cell')}; housing=prepared.rt['observe_initial_housing_wealth'](ev,P,grid,sd,diagnostic_enabled=True,age_projection='uniform_within_age_cell',diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True); recent=prepared.rt['observe_recent_parent_flow'](ev,P,diagnostic_enabled=True,snapshot=prepared.rt['SNAPSHOT'],age_projection=prepared.rt['AGE_PROJECTION'],diagnostic_allow_residence_proxy=True,input_provenance={'case':a.case})
 fits=runtime.score_targets(obj,fertility,housing,recent['model_value'],completed); params=copy.deepcopy(manifest['full_parameter_table'])
 for name,rows in [('target_fit.csv',fits),('parameters.csv',params)]:
  with open(out/name,'w',newline='') as f: w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 selected=abs(closure['renewal_residual'])<=RENEWAL and gate_failure is None
 native=None
 if selected:
  try: native=parent.scaled_native_step(packet,prepared,closure)
  except Exception as e: selected=False; gate_failure={'type':type(e).__name__,'message':str(e)}; fiscal_failure=gate_failure
 plots=[]
 if selected and time.time()+10<a.deadline:
  # Parent GE's reporting-unit convention is retained: economics clears in
  # absolute units while the diagnostic packet displays rooms per household.
  report_ev=copy.copy(ev); report_ev.supply_by_loc=np.asarray(ev.supply_by_loc)/closure['population_scale']
  prepared.rt['audit'].standard_diagnostics(dict(packet,evaluation=report_ev),out,validate_production_young=False)
  plots=sorted(x.name for x in (out/'standard_diagnostics').glob('*.png'))
  req(len(plots)==17,'standard diagnostic set incomplete')
 write(out/'receipt.json',dict(status='preliminary' if not selected else 'root_gate_passed_unverified_repeat',case=a.case,lifecycle_solves=1,price=float(price[0]),price_factor=a.factor,rent=float(P.user_cost_rate*price[0]),renewal_residual=closure['renewal_residual'],paygo_residual=fiscal,closure=closure,gate_failure=gate_failure,actual_fiscal_gate=gate_failure is None,native_16_20_queue=native,override=override,grid_nodes=len(grid),fixed_psi=float(P.psi_child),completed_fertility=completed,tenure_births_consumption_rooms=core,normalization_performed=False,repeat_verified=False,standard_plots_rendered=len(plots)==17,standard_plot_count=len(plots),deadline=a.deadline))

def controller(a):
 p=read(a.plan); verify_inputs(p); out=a.output; out.mkdir(parents=True,exist_ok=False); start=time.time(); hard=min(start+p['total_seconds'],p['absolute_deadline_epoch']-60); req(start<hard,'absolute readout reserve already reached')
 write(out/'launch.json',dict(start_epoch=start,hard_stop_epoch=hard,deadline_utc=p['absolute_deadline_utc'],case=a.case,plan_sha256=sha(a.plan))); rows=[]; write(out/'latest_completed.json',{'completed':[]});write(out/'best_so_far.json',{'status':'none'})
 def run(name,factor):
  req(len(rows)<8 and time.time()+180<hard,'solve budget/deadline reserve reached')
  cd=min(hard,time.time()+180); cmd=[sys.executable,str(Path(__file__).resolve()),'--plan',str(a.plan),'--output',str(out/name),'--case',a.case,'--factor',repr(factor),'--deadline',repr(cd)]
  with open(out/(name+'.log'),'w') as log:
   z=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
   while z.poll() is None:
    progress(out,'case_running',case=name,factor=factor,completed=len(rows),deadline=cd)
    if time.time()>=cd: os.killpg(z.pid,signal.SIGTERM); raise TimeoutError(name)
    time.sleep(1)
  req(z.returncode==0,'case failed without retry: '+name); r=read(out/name/'receipt.json'); row=dict(case=name,factor=factor,residual=r['renewal_residual'],price=r['price'],rent=r['rent'],population=r['closure']['population_scale'],gate_passed=bool(r['actual_fiscal_gate'] and r['native_16_20_queue'] is not None),receipt=str(out/name/'receipt.json'));rows.append(row);write(out/'latest_completed.json',{'completed':rows,'lifecycle_solves':len(rows)}); eligible=[x for x in rows if x['gate_passed']];write(out/'best_so_far.json',min(eligible,key=lambda x:abs(x['residual'])) if eligible else {'status':'no_gate_passed_candidate'});return row
 q=run('q0',1.); selected=q if root_eligible(q) else None; b=None
 if not selected:
  direction=1.03 if q['residual']>0 else .97
  for i,f in enumerate((direction,1.08 if q['residual']>0 else .92,1.20 if q['residual']>0 else .80)):
   r=run('bracket_'+str(i+1),f); b=bracket(rows)
   if root_eligible(r): selected=r;break
   if b:break
 if not selected and b:
  while len(rows)<8:
   r=run('root_'+str(len(rows)),logsecant(*b))
   if root_eligible(r): selected=r;break
   b=bracket(rows);write(out/'bracket.json',{'left':b[0],'right':b[1]} if b else {'status':'lost'})
 write(out/'completed.json',dict(status='preliminary_root_found_unverified_repeat' if selected else 'no_cleared_GE_root',selected=selected,records=rows,repeat_required=True,reason='No exact repeat or 17-plot rendering under the urgent eight-solve/deadline contract.'))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--case',choices=('ours','author'),required=True);ap.add_argument('--factor',type=float);ap.add_argument('--deadline',type=float);a=ap.parse_args();a.output=a.output.resolve(); req(os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm required')
 try: child(a) if a.factor is not None else controller(a)
 except BaseException as e:
  if a.output.exists(): write(a.output/'failure.json',{'status':'failed','type':type(e).__name__,'message':str(e),'traceback':traceback.format_exc(),'retries':0})
  raise
if __name__=='__main__': main()
