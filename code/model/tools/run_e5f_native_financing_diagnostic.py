"""Full-lifecycle, fixed-price financing diagnostic from the frozen native packet.

The experiment uses the saved beginning-of-period mass and the native calendar
period mapper.  It is partial equilibrium: a changed stationary distribution,
market residual, or fiscal residual is never reported as a new equilibrium.
"""
from __future__ import annotations
import argparse, copy, gzip, hashlib, importlib, json, pickle, os, signal, subprocess, sys, time
from pathlib import Path
from typing import Any, Mapping
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
DEFAULT_CHECKPOINT=ROOT/'output/model/paper_baseline_sep14/replay_20260917/native_output/raw/repetition_02/initial_state.pkl.gz'
DEFAULT_REPLAY=ROOT/'output/model/paper_baseline_sep14/replay_20260917'
DEFAULT_SOURCE_ROOT=ROOT/'tmp/paper_baseline_sep14/code/model'
DEFAULT_OUTPUT=ROOT/'output/model/native_financing_diagnostic_20260919'
CHECKPOINT_SHA256='3322a61994fb3654d67f4b1d6cf2d0f7cacbb3668d06a417e192ee363c174993'
SOLVER_SHA256='2992412586b81cef3a3e58d92191bb51f54d3f9cc600d7675bbadaed7d1682da'
CASES={'baseline':(.8,0.),'mortgage_only':(1.,0.),'unsecured_only':(.8,5.),'both':(1.,5.)}
POLICY_NAMES=('V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value','fert2_probs','price')

def sha(path:Path)->str:
 h=hashlib.sha256()
 with path.open('rb') as f:
  for x in iter(lambda:f.read(1048576),b''): h.update(x)
 return h.hexdigest()
def get(x:Any,k:str,d:Any=None)->Any: return x.get(k,d) if isinstance(x,Mapping) else getattr(x,k,d)
def write(path:Path,x:Mapping[str,Any])->None:
 path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+'.tmp'); tmp.write_text(json.dumps(x,indent=2,sort_keys=True,default=str)+'\n');tmp.replace(path)
def install_paths(root:Path)->None:
 for p in (root,root/'tools'):
  if str(p.resolve()) not in sys.path: sys.path.insert(0,str(p.resolve()))
def validate_contract(checkpoint:Path,replay:Path,source:Path)->dict[str,str]:
 selected_path=replay/'native_output/selected_checkpoint.json'
 if not selected_path.exists():raise ValueError('checkpoint contract hash mismatch: selected checkpoint record missing')
 actual=sha(checkpoint); selected=json.loads(selected_path.read_text()).get('checkpoint_sha256')
 if actual!=CHECKPOINT_SHA256 or selected!=CHECKPOINT_SHA256: raise ValueError('checkpoint contract hash mismatch')
 solver=source/'intergen_eqscale_seq_optimized/solver.py'
 if sha(solver)!=SOLVER_SHA256: raise ValueError('frozen solver hash mismatch')
 return {'checkpoint_sha256':actual,'solver_sha256':sha(solver)}
def packet(path:Path)->Any:
 with gzip.open(path,'rb') as f:return pickle.load(f)
def changed(a:Any,b:Any)->list[str]:
 def eq(x,y):
  try:return np.array_equal(x,y) if isinstance(x,np.ndarray) or isinstance(y,np.ndarray) else bool(x==y)
  except Exception:return False
 return sorted(k for k in set(vars(a))|set(vars(b)) if not eq(getattr(a,k,None),getattr(b,k,None)))
def refresh_finance(P:Any)->None:
 importlib.import_module('intergen_eqscale_seq_optimized.parameters').build_debt_caps(P)
 if not(np.isfinite(P.debt_taper_weights).all() and np.isfinite(P.debt_caps).all()):raise ValueError('nonfinite rebuilt finance arrays')
def arm_parameters(base:Any,case:str)->Any:
 if case not in CASES:raise KeyError(case)
 P=copy.deepcopy(base);phi,lam=CASES[case];P.phi=np.full_like(np.asarray(base.phi,dtype=float),phi);P.lambda_d=lam
 if lam:P.debt_taper_start_age,P.debt_taper_end_age=82.,86.
 refresh_finance(P); altered=changed(base,P);allowed={'phi','lambda_d','debt_taper_start_age','debt_taper_end_age','debt_taper_weights','debt_caps'}
 if set(altered)-allowed:raise ValueError(f'unaccounted parameter changes: {set(altered)-allowed}')
 if case=='baseline' and altered:raise ValueError(f'baseline changed fields: {altered}')
 return P
def policy_arrays(policy:Any)->dict[str,np.ndarray]:
 missing=[n for n in POLICY_NAMES if get(policy,n) is None]
 if missing:raise ValueError(f'missing mandatory policy arrays: {missing}')
 return {n:np.asarray(get(policy,n)) for n in POLICY_NAMES}
def compare_exact(saved:Any,observed:Any)->None:
 a,b=policy_arrays(saved),policy_arrays(observed)
 if set(a)!=set(b):raise ValueError('mandatory policy schema mismatch')
 for n in a:np.testing.assert_allclose(b[n],a[n],atol=1e-10,rtol=0,err_msg=n)
def gates(ev:Any)->dict[str,float]:
 out={'pre_mass':float(ev.g_pre.sum()),'post_fertility_mass':float(ev.g_post_fertility.sum()),'current_mass':float(ev.g_current.sum()),'birth_mass':float(np.asarray(ev.births).sum())}
 if any(np.min(a)<-1e-12 for a in (ev.g_pre,ev.g_post_fertility,ev.g_current)):raise ValueError('negative household mass')
 if not np.isfinite(list(out.values())).all() or min(out['pre_mass'],out['post_fertility_mass'],out['current_mass']) < -1e-12 or max(abs(out['pre_mass']-out['post_fertility_mass']),abs(out['pre_mass']-out['current_mass']))>1e-10:raise ValueError(f'period mass gate failed: {out}')
 arrays=policy_arrays(ev.policy)
 for n,a in arrays.items():
  if not np.isfinite(a).all():raise ValueError(f'nonfinite policy: {n}')
 for n in ('tenure_probs','loc_probs','fert_probs','fert2_probs'):
  if arrays[n].min()<-1e-12 or arrays[n].max()>1+1e-12:raise ValueError(f'probability gate failed: {n}')
 return out
def native_arm(x:Mapping[str,Any],P:Any)->tuple[Any,Any]:
 model=importlib.import_module('intergen_eqscale_seq_optimized.solver'); primitive=importlib.import_module('run_e5f_matched_pf_smoke');price=np.asarray(x['evaluation'].policy.price);grid=x['b_grid'];chain,model=primitive.pf.transition.configure_sequential_model();primitive.pf.calendar.apply_fertility=primitive.pf.transition.apply_sequential_fertility;primitive.pf.calendar.advance_calendar_distribution=primitive.pf.transition.advance_sequential_calendar_distribution;shared=model.precompute_shared(P,grid)
 sol=model.solve_markov_income_at_prices(price,P,grid,verbose=False,fast_stats=False);P._fert2_probs=np.asarray(sol.fert2_probs).copy()
 policy=primitive.pf.calendar.policy_from_solution(sol,price,P,grid,shared)
 ev=primitive.pf.calendar.evaluate_period(price,x['stationary_g_pre'],P,grid,shared,primitive.pf.calendar.SolveCounter(),supply_rule=x['supply_rule'],supplied_policy=policy)
 return ev,primitive.dated_budget(ev,P,shared,grid,float(P.user_cost_rate*price[0]))
def inspect(args:Any)->dict[str,Any]:
 c=validate_contract(args.checkpoint,args.replay,args.source_root);install_paths(args.source_root);primitive=importlib.import_module('run_e5f_matched_pf_smoke');
 if not(hasattr(primitive,'dated_budget') and hasattr(primitive,'pf') and hasattr(primitive.pf.calendar,'evaluate_period')):raise ValueError('native primitive binding unavailable')
 x=packet(args.checkpoint);gates(x['evaluation']);P=x['parameters'];cases=[]
 for case in CASES:
  q=arm_parameters(P,case);cases.append({'case':case,'changed_fields':changed(P,q),'phi':np.asarray(q.phi).tolist(),'lambda_d':q.lambda_d,'taper':[q.debt_taper_start_age,q.debt_taper_end_age]})
 return {'status':'inspect_only','contract':c,'source_root':str(args.source_root),'packet_keys':sorted(x),'g_pre_shape':list(x['stationary_g_pre'].shape),'g_pre_mass':float(x['stationary_g_pre'].sum()),'policy_shapes':{k:list(v.shape) for k,v in policy_arrays(x['evaluation'].policy).items()},'cases':cases,'scope':'fixed-price full-lifecycle PE; saved pre-choice mass; no GE/closure claim'}
def run_case(args:Any)->dict[str,Any]:
 c=validate_contract(args.checkpoint,args.replay,args.source_root);install_paths(args.source_root);x=packet(args.checkpoint);P=arm_parameters(x['parameters'],args.case);finance_changes=changed(x['parameters'],P);started=time.monotonic();ev,budget=native_arm(x,P)
 if float(budget.get('budget_excess_mass',np.inf))>2e-10 or float(budget.get('maximum_occupied_excess',np.inf))>1e-9:raise ValueError(f'native budget gate failed: {budget}')
 if not np.array_equal(ev.g_pre,x['stationary_g_pre']):raise ValueError('period mapper changed saved beginning-of-period mass')
 if args.case=='baseline':
  compare_exact(x['evaluation'].policy,ev.policy)
  for n,a,b in (('g_current',x['evaluation'].g_current,ev.g_current),('births',x['evaluation'].births,ev.births)):
   np.testing.assert_allclose(np.asarray(b),np.asarray(a),atol=1e-10,rtol=0,err_msg=n)
 r={'status':'completed','case':args.case,'label':args.case_label or args.case,'elapsed_seconds':time.monotonic()-started,'contract':c,'parameter_changes':finance_changes,'mass':gates(ev),'relative_market_residual':float(ev.relative_market_residual),'budget':budget,'interpretation':'fixed-price PE; residuals are diagnostics, not clearing gates'};out=args.output/'cases'/(args.case_label or args.case)
 if (out/'receipt.json').exists():raise FileExistsError(f'refusing to overwrite existing receipt: {out}')
 out.mkdir(parents=True,exist_ok=True)
 audit=importlib.import_module('run_e5f_independent_numerical_audit')
 r['policy_audit']=audit.policy_array_audit({'evaluation':ev,'parameters':P,'b_grid':x['b_grid']},out)
 if r['policy_audit']['occupied_negative_steps']:raise ValueError('occupied value monotonicity gate failed')
 np.savez_compressed(out/'arrays.npz',g_pre=ev.g_pre,g_post_fertility=ev.g_post_fertility,g_current=ev.g_current,births=ev.births,**policy_arrays(ev.policy));write(out/'receipt.json',r);return r
def sequence(args:Any)->None:
 """Two exact controls, then treatments; each arm is a killable process group."""
 if args.case_budget_seconds<=0 or args.total_budget_seconds<=0:raise ValueError('time budgets must be positive')
 started=time.monotonic(); order=('baseline','baseline','mortgage_only','unsecured_only','both')
 for i,case in enumerate(order,1):
  left=args.total_budget_seconds-(time.monotonic()-started)
  if left<=0:raise TimeoutError('total diagnostic budget exhausted before launch')
  label=f'{case}_{i:02d}';cmd=[sys.executable,str(Path(__file__).resolve()),'--mode','case','--case',case,'--case-label',label,'--checkpoint',str(args.checkpoint),'--replay',str(args.replay),'--source-root',str(args.source_root),'--output',str(args.output)]
  child=subprocess.Popen(cmd,start_new_session=True);limit=min(args.case_budget_seconds,left)
  case_start=time.monotonic()
  try:
   while True:
    remaining=limit-(time.monotonic()-case_start)
    if remaining<=0:raise subprocess.TimeoutExpired(cmd,limit)
    try:child.wait(timeout=min(30.,remaining));break
    except subprocess.TimeoutExpired:
     write(args.output/'heartbeat.json',{'case':case,'elapsed_seconds':time.monotonic()-started,'case_elapsed_seconds':time.monotonic()-case_start})
  except subprocess.TimeoutExpired:
   write(args.output/'timeout.json',{'status':'timeout','case':case,'limit_seconds':limit,'elapsed_seconds':time.monotonic()-started})
   os.killpg(child.pid,signal.SIGTERM)
   try: child.wait(timeout=10)
   except subprocess.TimeoutExpired: os.killpg(child.pid,signal.SIGKILL); child.wait(timeout=10)
   raise TimeoutError(f'{case} exceeded {limit:g}s budget')
  if child.returncode:raise RuntimeError(f'{case} child exited {child.returncode}')
  receipt=json.loads((args.output/'cases'/label/'receipt.json').read_text())
  write(args.output/'latest_completed.json',{'status':'progress','completed':i,'of':len(order),'last_case':case,'elapsed_seconds':time.monotonic()-started,'receipt':receipt})
def main(argv:list[str]|None=None)->None:
 p=argparse.ArgumentParser();p.add_argument('--checkpoint',type=Path,default=DEFAULT_CHECKPOINT);p.add_argument('--replay',type=Path,default=DEFAULT_REPLAY);p.add_argument('--source-root',type=Path,default=DEFAULT_SOURCE_ROOT);p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT);p.add_argument('--mode',choices=('inspect-only','case','sequence'),default='inspect-only');p.add_argument('--case',choices=tuple(CASES));p.add_argument('--case-label');p.add_argument('--case-budget-seconds',type=float,default=300);p.add_argument('--total-budget-seconds',type=float,default=1500);a=p.parse_args(argv)
 try:
  if a.mode=='inspect-only':write(a.output/'inspect.json',inspect(a))
  elif a.mode=='sequence':sequence(a)
  elif a.case:run_case(a)
  else:p.error('--case required')
 except Exception as exc:
  write(a.output/('failure_'+(a.case_label or a.mode)+'.json'),{'status':'failed','type':type(exc).__name__,'error':str(exc),'case':a.case})
  raise
if __name__=='__main__':main()
