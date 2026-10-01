#!/usr/bin/env python3
"""Six native policy calls, same price/PRE; matched preferences and broad credit."""
import os,sys,copy,time,signal,json,csv,hashlib
from pathlib import Path
from types import SimpleNamespace
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
HERE=Path(__file__).resolve().parent;WIN=HERE.parents[1];sys.path.insert(0,str(WIN))
import fixed_price_responses as d
import numpy as np
BASE=WIN/'purchase_ltv_v1/local_run/retry5/results/baseline_80_80';BROAD=WIN/'collected/responses/04_lifetime_repayment_only_p1.00';TREAT=WIN/'collected/treatment_metadata.json';HIST=d.ROOT/'output/model/fixed_reference_economics_20260928/credit_v1/solve_v1/grid_control/parameters.csv'
FIELDS=['V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value']
def alarm(*_):raise TimeoutError('180-second six-policy-call budget')
def main():
 started=time.time();auth=d.authenticate_candidate(HERE/'runtime_preparation');P0=auth['P'];grid=auth['grid'];q=float(d.BINDING['candidate_price'])
 from refactor_lab.engine import household as h
 from refactor_lab.engine.parameters import get_fecundity_by_age
 treatment=json.loads(TREAT.read_text());historical={x['parameter']:float(x['estimate'])for x in csv.DictReader(HIST.open())};loading=historical['delta_alpha_jump'];d.require(loading==.13487659325900622,'Historical share loading differs');d.require(historical['utility_reference_rent']==P0.utility_reference_rent==.11046592704873838,'Reference rent differs')
 d.require(d.sha(h.__file__)==treatment['source']['pinned_household_sha256'],'Matched native household source differs')
 pre=np.load(BASE.parent/'q0_reference_inherited_states.npz')['g_pre'];commonhash=hashlib.sha256(pre.tobytes()).hexdigest();d.require(commonhash=='459ab9229e7ce58376bda8c933b583ec059104ede31c57627ce7368c41e5a51b','Common PRE fingerprint differs');young=np.zeros_like(pre);young[:,0,0,:7,:,0,0]=pre[:,0,0,:7,:,0,0];youngmass=float(young.sum())
 cal=auth['context']['prepared'].rt['primitive'].pf.calendar;rows=[];values={};records=[];deadline=time.time()+180
 for utility in ('floor','shares_A','shares_no_A'):
  for regime in ('reference','broad'):
   label=utility+'_'+regime;out=HERE/label;out.mkdir()
   P=copy.deepcopy(auth['P'] if regime=='reference' else auth['natural']);P.native_inherited_distribution_evidence_dir=str(out/'inherited_state_diagnostics')
   if utility!='floor':P.hbar_first_child_jump=0.;P.hbar_child_rooms=0.;P.child_room_floor=False;P.delta_alpha_jump=loading;P.compensated_child_housing_shares=(utility=='shares_A')
   expected=dict(auth['actual_parameters'])
   if utility!='floor':expected.update(h_P=0.,delta_alpha_jump=loading)
   actual=auth['context']['fp'].actual_parameters(auth['context']['prepared'],P,grid);d.require(actual==expected,'Other31economic scalars changed')
   flags={key:getattr(P,key,False)for key in ('native_due_stayer_credit','native_solvency_credit','native_purchase_income','use_pti_constraint','joint_nested_choice','fertility_nest_choice','parent_dp_waiver','exhaustive_saving_control')};flags.update(unsecured_credit_limit=getattr(P,'unsecured_credit_limit',None),stored_phi=np.asarray(P.phi).tolist(),child_room_floor=P.child_room_floor,compensated_child_housing_shares=P.compensated_child_housing_shares)
   expectedcredit=treatment['reference_credit' if regime=='reference' else 'expanded_credit'];d.require(all(flags[k]==v for k,v in expectedcredit.items()),'Broad/reference financial flags differ');d.require(not flags['use_pti_constraint'] and not flags['joint_nested_choice'] and not flags['fertility_nest_choice'] and not flags['parent_dp_waiver'],'Extra choice/financial mechanism');d.require(np.all(np.asarray(P.phi)==.8),'Storedphi changed')
   if regime=='broad':d.require(auth['solver'].validate_native_solvency_mode(P),'Broad mode validation failed')
   sd=auth['solver'].precompute_shared(P,grid);d.write(HERE/'latest.json',dict(status='policy_call_claimed',label=label,policy_calls=len(records)+1,deadline_epoch=deadline))
   old=signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.001,deadline-time.time()));t=time.monotonic()
   try:objects=h.solve_bellman_full_markov_income(np.array([P.user_cost_rate*q]),np.array([q]),P,grid,sd)
   finally:signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
   d.require(time.time()<deadline,'180-second policy budget exceeded');sol=SimpleNamespace(**dict(zip(FIELDS,objects[:9])),fert2_probs=P._fert2_probs.copy(),bp_pol_stay=getattr(P,'_bp_pol_stay',None),c_pol_stay=getattr(P,'_c_pol_stay',None));elapsed=time.monotonic()-t
   policy=cal.policy_from_solution(sol,np.array([q]),P,grid,sd);counter=cal.SolveCounter();ev=cal.evaluate_period(np.array([q]),pre,P,grid,sd,counter,supplied_policy=policy)
   d.require(counter.total==0,'Native calendar launched a solve');d.require(ev.feasibility_projection_mass==0.,'Positive occupied PRE projection');d.require(abs(float(ev.g_pre.sum())-float(pre.sum()))<2e-10 and abs(float(ev.g_current.sum())-float(pre.sum()))<2e-10,'WholePRE mass changed')
   ag=auth['report_helpers'].aggregates(auth['context']['fp'],ev,P,grid,auth['context']['prepared'].rt['model'])
   ye=cal.evaluate_period(np.array([q]),young,P,grid,sd,cal.SolveCounter(),supplied_policy=policy);d.require(ye.feasibility_projection_mass==0. and abs(float(ye.g_current.sum())-youngmass)<2e-10,'YoungPRE projected/changed')
   row=dict(utility=utility,credit=regime,whole_PRE_birth_children=float(ev.births),whole_PRE_first_births=float(ag['first_births']),whole_PRE_ownership=float(ag['ownership_rate']),whole_PRE_mean_rooms=float(ag['rooms_per_household']),young_n0_renter_first_birth_probability=float(ye.births)/youngmass,whole_projection_mass=float(ev.feasibility_projection_mass),young_projection_mass=float(ye.feasibility_projection_mass),policy_seconds=elapsed)
   if utility=='floor':
    saved=json.loads(((BASE if regime=='reference' else BROAD)/'closure.json').read_text())['baseline_state_impact']
    errors={k:abs(row[rk]-float(saved[k]))for k,rk in [('births','whole_PRE_birth_children'),('first_births','whole_PRE_first_births'),('ownership_rate','whole_PRE_ownership')]};d.require(max(errors.values())<2e-9,'Floor control not reproduced');row['control_max_error']=max(errors.values())
   # Matched value gaps, with actual conception success and first-birth scale.
   w=young[:,0,0,:,:,0,0];fp=sol.fert_probs[:,0,0,:,:,:2];F=sol.fert_value[:,0,0,:,:];pi=get_fecundity_by_age(P)[None,:,None];valid=(w>0)&(fp[...,0]>0)&(fp[...,1]>0)&(pi>0);gap=np.zeros_like(F);gap[valid]=P.kappa_fert*np.log(fp[...,1][valid]/fp[...,0][valid]);values[utility,regime]=(w.copy(),valid.copy(),gap.copy(),fp[...,1].copy(),pi.copy())
   rows.append(row);records.append(dict(label=label,actual31_parameters=actual,effective_flags=flags,source_sha256=d.sha(h.__file__),continuation_V='omitted: full reoptimized lifetime recursion',utility_specification='floor unchanged'if utility=='floor'else'no physical floor; historical first-child shares; '+('compensatedA'if utility=='shares_A'else'noA'),seconds=elapsed));d.write(out/'receipt.json',dict(outcomes=row,record=records[-1]));d.write(HERE/'latest.json',dict(status='completed_case',completed=rows))
   del sol,objects,policy,ev,ye
 comparisons=[]
 for utility in ('floor','shares_A','shares_no_A'):
  b=next(x for x in rows if x['utility']==utility and x['credit']=='reference');a=next(x for x in rows if x['utility']==utility and x['credit']=='broad');w,v0,g0,p0,pi=values[utility,'reference'];_,v1,g1,p1,_=values[utility,'broad'];mask=v0&v1;ww=w*mask;db=pi*(p1-p0)
  comparisons.append(dict(utility=utility,baseline_whole_births=b['whole_PRE_birth_children'],broad_whole_births=a['whole_PRE_birth_children'],whole_birth_change=a['whole_PRE_birth_children']-b['whole_PRE_birth_children'],whole_first_birth_change=a['whole_PRE_first_births']-b['whole_PRE_first_births'],whole_ownership_change=a['whole_PRE_ownership']-b['whole_PRE_ownership'],baseline_young_firstbirth=b['young_n0_renter_first_birth_probability'],broad_young_firstbirth=a['young_n0_renter_first_birth_probability'],young_firstbirth_credit_effect_pp=100*(a['young_n0_renter_first_birth_probability']-b['young_n0_renter_first_birth_probability']),young_positive_birth_change_mass=float(np.sum(w*np.maximum(db,0))),young_negative_birth_change_mass=float(np.sum(w*np.minimum(db,0))),mean_attempt_gap_change_interior=float(np.sum(ww*(g1-g0))/ww.sum()),value_common_interior_mass=float(ww.sum())))
 d.table(HERE/'six_policy_cases.csv',rows);d.table(HERE/'matched_credit_effects.csv',comparisons);d.write(HERE/'summary.json',dict(status='six_matched_policy_calls_complete',policy_calls=6,KFE_calls=0,GE_calls=0,policy_seconds=sum(x['seconds']for x in records),total_seconds=time.time()-started,price=q,common_PRE_sha256=commonhash,whole_PRE_mass=float(pre.sum()),young_PRE_mass=youngmass,source_pins=dict(native_household=d.sha(h.__file__),winner_driver=d.sha(WIN/'fixed_price_responses.py'),treatment_metadata=d.sha(TREAT),historical_parameter_table=d.sha(HIST),driver=d.sha(__file__)),treatment_metadata=treatment,records=records,comparisons=comparisons,interpretation='Preference packages and broad all-tenure credit at prescribed price on identical inherited PRE; not floor-only attribution, historical calibration replication, candidate, stationary population or GE comparison',limit='Original native financial/feasibility operators retained; nopositivePREprojection accepted. No full forward lifecycle/death ledger/17plots rerun; broad unoccupied support/grid convergence remains uncertified.'))
 print(json.dumps(comparisons,indent=2))
if __name__=='__main__':main()
