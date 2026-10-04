"""Bounded CES normalized-share calibration; no implicit retries or adoption."""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, math, os, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parents[4]
PLAN=ROOT/'output/model/experiments/ces_normalized_shares/overnight_v1/start_plan.json'
RESERVE=1800; MAX_CALLS=500; MAX_LIFECYCLE=32; PENALTY=1e12
RAW_STEPS={'beta_annual':.0005,'chi':.02,'first_birth_fixed_cost':.015,
 'kappa_fert':.01,'kappa_fert_continuation':.01,'delta_alpha_jump':.01,
 'delta_alpha':.005,'tenure_choice_kappa':.001,'theta0':.005,'psi_child':.004,'child_benefit_curvature':.01}
class BudgetStop(RuntimeError): pass
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(path,value):
 path=Path(path); tmp=path.with_name(path.name+'.tmp')
 tmp.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'); os.replace(tmp,path)
def norm_bounds(value): return {str(k):[float(v[0]),float(v[1])] for k,v in value.items()}
def adapter():
 p=Path(__file__).with_name('adapter.py'); spec=importlib.util.spec_from_file_location('ces_adapter',p)
 model_root=str(ROOT/'code/model')
 if model_root not in sys.path: sys.path.insert(0,model_root)
 m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
 missing=[x for x in ('contract','load_inputs','bind_parameters','install','make_evaluator') if not callable(getattr(m,x,None))]
 if missing: raise RuntimeError('CES adapter contract missing: '+', '.join(missing))
 return m
def checked_plan(path,expected):
 if sha(path)!=expected: raise RuntimeError('start-plan SHA-256 drift')
 p=json.loads(Path(path).read_text())
 if len(p.get('starts',()))!=4 or len(p.get('coordinates',()))!=11 or norm_bounds(p['bounds']).get('delta_alpha_jump')!=[0.,.25] or norm_bounds(p['bounds']).get('delta_alpha')!=[0.,.25] or 'h_P' in p['coordinates'] or sum(r['role']=='scored' for r in p['target_contract'])!=11: raise RuntimeError('CES plan contract drift')
 return p
def reflect(z):
 z=float(z)%2.; return z if z<=1. else 2.-z
def as_point(x,coords,bounds):
 q={k:float(bounds[k][0]+reflect(x[j])*(bounds[k][1]-bounds[k][0])) for j,k in enumerate(coords)}
 if not all(math.isfinite(v) and bounds[k][0]<=v<=bounds[k][1] for k,v in q.items()): raise RuntimeError('nonfinite/out-of-bound candidate')
 return q
def simplex(seed,coords,bounds):
 lo=np.array([bounds[k][0] for k in coords]); span=np.array([bounds[k][1]-bounds[k][0] for k in coords])
 x0=np.array([(seed[k]-lo[j])/span[j] for j,k in enumerate(coords)])
 raw=np.array([RAW_STEPS[k] for k in coords])/span
 ans=np.vstack([x0]+[np.array([reflect(v) for v in x0+np.eye(len(coords))[j]*raw[j]]) for j in range(len(coords))])
 if len({tuple(np.round(v,15)) for v in ans})!=len(ans): raise RuntimeError('initial simplex has duplicate vertices')
 return lo,span,ans
def hashes(report):
 d=Path(report)/'standard_diagnostics'; h={p.name:sha(p) for p in sorted(d.glob('*.png'))}
 if len(h)!=17: raise RuntimeError('expected exactly 17 standard diagnostic PNGs')
 return h
def csv_rows(report,name,expected):
 p=Path(report)/name
 if not p.is_file(): raise RuntimeError('missing native report table: '+str(p))
 with p.open(newline='') as f: rows=list(csv.DictReader(f))
 if len(rows)!=expected: raise RuntimeError(f'{name} has {len(rows)}, expected {expected} rows')
 return rows
def utility_pins(report):
 p=Path(report)/'utility_contract.json'
 if not p.is_file(): raise RuntimeError('missing utility contract: '+str(p))
 value=json.loads(p.read_text())
 required=('immutable_contract_fingerprint','effective_input_fingerprint','experiment')
 if any(k not in value for k in required): raise RuntimeError('utility contract lacks immutable pins')
 return {k:value[k] for k in required}
def validate_selected_identity(rows,chosen,coords,bounds):
 by_name={row.get('parameter'):row for row in rows}
 if len(by_name)!=31 or set(chosen)!=set(coords): raise RuntimeError('selected coordinate identity drift')
 for key in coords:
  row=by_name.get(key)
  if row is None: raise RuntimeError('chosen coordinate absent from parameter table: '+key)
  value=float(row['estimate'])
  agrees = abs(value-float(chosen[key]))<=2e-12 if key=='beta_annual' else value==float(chosen[key])
  if not agrees or not bounds[key][0]<=value<=bounds[key][1]: raise RuntimeError('chosen coordinate/table/bound mismatch: '+key)
 for key in ('delta_alpha_jump','delta_alpha'):
  if by_name.get(key,{}).get('lower') not in ('0.0','0') or by_name.get(key,{}).get('upper') not in ('0.25','.25'):
   raise RuntimeError(key+' bounds not recorded as free [0,.25]')
 if float(by_name.get('h_P',{}).get('estimate','nan'))!=0.0: raise RuntimeError('h_P must remain fixed at zero')
def exact_postcheck(evaluate,chosen,deadline,coords,bounds):
 if set(chosen)!=set(coords): raise RuntimeError('selected coordinates drift')
 for k,v in chosen.items():
  if not bounds[k][0]<=float(v)<=bounds[k][1]: raise RuntimeError('selected coordinate outside bound: '+k)
 # One native GE already performs its internal exact selected repeat, including arrays.
 root=evaluate('selected_root',chosen,deadline)
 if root.get('status')!='passed' or len(root.get('target_fit',()))!=14 or len(root.get('parameter_table',()))!=31 or len(root.get('residual',()))!=11: raise RuntimeError('experimental postcheck lacks 14-row/11-residual/31-parameter reports')
 root_report=Path(root['report']); repeat_report=root_report.parent/'selected_repeat_final'
 root_fit=csv_rows(root_report,'target_fit.csv',14); repeat_fit=csv_rows(repeat_report,'target_fit.csv',14)
 root_parameters=csv_rows(root_report,'parameters.csv',31); repeat_parameters=csv_rows(repeat_report,'parameters.csv',31)
 if root_fit!=repeat_fit: raise RuntimeError('internal native repeat 14-target table differs')
 if root_parameters!=repeat_parameters: raise RuntimeError('internal native repeat 31-parameter table differs')
 root_exp=csv_rows(root_report,'target_fit_experimental.csv',14); repeat_exp=csv_rows(repeat_report,'target_fit_experimental.csv',14)
 if root_exp!=repeat_exp or sum(r['role']=='scored' for r in root_exp)!=11: raise RuntimeError('independently rescored experimental CSV differs or lacks eleven scored rows')
 root_hashes,repeat_hashes=hashes(root_report),hashes(repeat_report)
 if root_hashes!=repeat_hashes: raise RuntimeError('internal native repeat 17 PNG hashes differ')
 if utility_pins(root_report)!=utility_pins(repeat_report): raise RuntimeError('root/repeat utility-contract pins differ')
 validate_selected_identity(root_parameters,chosen,coords,bounds)
 return dict(status='full_native_postcheck_passed',selected_root={**root, 'report':str(root_report)},
             selected_repeat_final=dict(report=str(repeat_report),utility_contract_pins=utility_pins(repeat_report)),
             target_rows=14,parameter_rows=31,residual_count=11,standard_plot_hashes=root_hashes,
             exact_tables=dict(target_fit_csv_equal=True,target_fit_experimental_csv_equal=True,parameters_csv_equal=True),
             selected_input_identity=dict(coordinates=chosen,bounds=bounds))
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--chain',type=int,required=True); ap.add_argument('--out',type=Path,required=True); ap.add_argument('--deadline-epoch',type=float,required=True); ap.add_argument('--starts-file',type=Path,default=PLAN); ap.add_argument('--starts-file-sha256'); ap.add_argument('--smoke',action='store_true'); ap.add_argument('--mock-smoke',action='store_true'); ap.add_argument('--preflight-evaluator',action='store_true'); ap.add_argument('--postcheck-only',action='store_true'); ap.add_argument('--search-receipt',type=Path); a=ap.parse_args()
 if a.out.exists(): raise RuntimeError('Refusing existing output directory')
 if a.postcheck_only != (a.search_receipt is not None): raise RuntimeError('postcheck requires exactly one search receipt')
 a.starts_file_sha256=a.starts_file_sha256 or sha(a.starts_file); plan=checked_plan(a.starts_file,a.starts_file_sha256)
 if not 0<=a.chain<4: raise RuntimeError('chain outside 0--3')
 a.out.mkdir(parents=True); out=a.out.resolve(); deadline=min(float(a.deadline_epoch),time.time()+21600); coords=plan['coordinates']; bounds=norm_bounds(plan['bounds']); seed=plan['starts'][a.chain]
 contract=dict(chain=a.chain,seed=seed,free_coordinates=coords,bounds=bounds,starts_file_sha256=a.starts_file_sha256,starts_count=4,target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'],selected_source_sha256=plan['selected_source_sha256'],source_checkpoint_sha256=plan['source_checkpoint_sha256'],objective_calls_max=MAX_CALLS,lifecycle_solves_per_GE_max=MAX_LIFECYCLE,final_native_reserve_seconds=RESERVE,no_auto_retry=True)
 atomic(out/'start_contract.json',contract)
 ad=None
 if not a.mock_smoke:
  ad=adapter(); ac=ad.contract()
  if not isinstance(ac,dict) or ac.get('target_fingerprint')!=plan['target_fingerprint'] or ac.get('weight_fingerprint')!=plan['weight_fingerprint'] or norm_bounds(ac.get('bounds',{}))!=bounds or ac.get('total_target_rows')!=14 or ac.get('scored_count')!=11: raise RuntimeError('Adapter target, bounds, or reporting contract conflict')
  P,grid=ad.load_inputs(seed)
  if a.preflight_evaluator:
   r=ad.preflight_contexts(out/'actual_reporting_context')
   if r.get('status')!='passed_zero_solves' or r.get('actual_lifecycle_solves')!=0: raise RuntimeError('actual reporting preflight failed')
   ad.make_evaluator(out,'ces_normalized_shares',P,grid,deadline,target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'])
   atomic(out/'completed.json',dict(status='evaluator_initialized_zero_solves',lifecycle_solves=0,preflight=r,contract=contract)); return
  evaluate=ad.make_evaluator(out,'ces_normalized_shares',P,grid,deadline,target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'])
  if a.postcheck_only:
   search=json.loads(a.search_receipt.read_text())
   if search.get('status')!='provisional_search_finished' or search.get('chain')!=a.chain or search.get('starts_file_sha256')!=a.starts_file_sha256 or search.get('target_fingerprint')!=plan['target_fingerprint'] or search.get('weight_fingerprint')!=plan['weight_fingerprint'] or search.get('selected_source_sha256')!=plan['selected_source_sha256']: raise RuntimeError('postcheck receipt authentication failed')
   done=exact_postcheck(evaluate,search.get('selected',{}).get('parameters',{}),deadline,coords,bounds); done.update(search_receipt_sha256=sha(a.search_receipt),starts_file_sha256=a.starts_file_sha256,target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'],selected_source_sha256=plan['selected_source_sha256'],source_checkpoint_sha256=plan['source_checkpoint_sha256']); atomic(out/'completed.json',done); return
 lo,span,initial=simplex(seed,coords,bounds); atomic(out/'search_contract.json',dict(method='bounded Nelder-Mead',initial_simplex=initial.tolist(),physical_simplex_steps={k:RAW_STEPS[k] for k in coords},max_objective_calls=2 if (a.smoke or a.mock_smoke) else MAX_CALLS,final_reserve_seconds=RESERVE))
 rows=[]; cache={}; best=None; calls=0; limit=2 if (a.smoke or a.mock_smoke) else MAX_CALLS
 def beat(status,**kw): atomic(out/'heartbeat.json',dict(epoch=time.time(),status=status,objective_calls=calls,completed_full_ge=len(rows),**kw))
 def objective(x):
  nonlocal calls,best
  if calls>=limit: raise BudgetStop('objective_call_limit')
  if not a.mock_smoke and time.time()>=deadline-RESERVE: raise BudgetStop('final_native_reserve_reached')
  if not a.mock_smoke and os.statvfs(out).f_bavail*os.statvfs(out).f_frsize<350*1024**3: raise BudgetStop('shared_free_disk_below_350GiB')
  calls+=1; q=as_point(x,coords,bounds); key=tuple(q[k].hex() for k in coords)
  if key in cache: beat('cache_hit'); return cache[key]
  label=f'{len(rows):04d}_nm'; beat('running_full_GE',label=label)
  if a.mock_smoke: result=dict(status='passed',residual=(np.asarray(x)-initial[0]-.01).tolist(),loss=float(np.sum((np.asarray(x)-initial[0]-.01)**2)),lifecycle_solves=0,mock=True)
  else: result=evaluate(label,q,deadline-RESERVE)
  row=dict(label=label,parameters=q,**result)
  if result.get('status')=='passed':
   rr=np.asarray(result.get('residual',[]),float); row['loss']=float(result.get('loss',rr@rr)); row['objective']=row['loss']; best=row if best is None or row['loss']<best['loss'] else best
  elif result.get('status')=='inadmissible_numerical': row.update(objective=PENALTY,numerical_rejection=True)
  elif result.get('status')=='budget_exhausted': atomic(out/'latest_completed.json',dict(latest=row,completed_full_ge=len(rows),objective_calls=calls)); raise BudgetStop('native_evaluation_budget_exhausted')
  else: raise RuntimeError('Unexpected evaluator status: '+str(result.get('status')))
  if not a.mock_smoke:
   size=sum(p.stat().st_size for p in (out/label).rglob('*') if p.is_file())
   if size>1024**3: raise BudgetStop('case_storage_exceeded_1GiB_planning_cap')
   row['retained_case_bytes']=size
  rows.append(row); cache[key]=row['objective']; atomic(out/'cases.json',rows); atomic(out/'latest_completed.json',dict(latest=row,completed_full_ge=len(rows),objective_calls=calls)); atomic(out/'best_so_far.json',dict(status='provisional_until_fresh_native_postcheck',best=best,completed_full_ge=len(rows))); beat('case_completed',best_loss=best['loss'] if best else None); return row['objective']
 beat('initialized')
 stop='optimizer_return'
 try: minimize(objective,initial[0],method='Nelder-Mead',bounds=[(0.,1.)]*len(coords),options=dict(initial_simplex=initial,maxfev=limit,maxiter=limit,adaptive=True))
 except BudgetStop as e: stop=str(e)
 search=dict(status='provisional_search_finished',search_stop_reason=stop,chain=a.chain,selected=best,objective_calls=calls,completed_full_ge=len(rows),lifecycle_solves=sum(r.get('lifecycle_solves',0) for r in rows),starts_file_sha256=a.starts_file_sha256,target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'],selected_source_sha256=plan['selected_source_sha256'],source_checkpoint_sha256=plan['source_checkpoint_sha256'],no_auto_retry=True)
 atomic(out/'search_completed.json',search)
 if a.mock_smoke:
  if calls!=2 or len(rows)!=2 or len(cache)!=2: raise RuntimeError('mock smoke must use exactly two distinct optimizer evaluations')
  atomic(out/'completed.json',dict(search,status='mock_loop_passed_zero_solves',contract=contract)); beat('mock_completed'); return
 if best is None: atomic(out/'completed.json',dict(search,status='no_admissible_candidate')); return
 beat('native_postcheck_running'); child=out/'native_postcheck'; cmd=[sys.executable,str(Path(__file__).resolve()),'--chain',str(a.chain),'--out',str(child),'--deadline-epoch',str(deadline),'--starts-file',str(a.starts_file),'--starts-file-sha256',a.starts_file_sha256,'--postcheck-only','--search-receipt',str(out/'search_completed.json')]
 r=subprocess.run(cmd,capture_output=True,text=True,timeout=max(1,min(RESERVE,deadline-time.time())),check=False); (out/'postcheck_child.stdout.log').write_text(r.stdout); (out/'postcheck_child.stderr.log').write_text(r.stderr)
 if r.returncode: raise RuntimeError('full native postcheck child failed: '+str(r.returncode))
 checked=json.loads((child/'completed.json').read_text())
 if checked.get('status')!='full_native_postcheck_passed': raise RuntimeError('fresh selected postcheck failed')
 atomic(out/'completed.json',dict(search,status='selected_numerically_verified',selected_postcheck=checked)); beat('completed')
if __name__=='__main__': main()
