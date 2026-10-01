"""Isolated CES scales experiment; native frozen observers and bounded NM."""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, os, signal, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
if sys.platform=='darwin':
    import importlib,pathlib
    sys.modules.setdefault('pathlib._local',pathlib)
    sys.modules.setdefault('numpy._core',np.core)
    sys.modules.setdefault('numpy._core.multiarray',importlib.import_module('numpy.core.multiarray'))
    sys.modules.setdefault('numpy._core.numeric',importlib.import_module('numpy.core.numeric'))
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
PACKET=ROOT/'output/model/fixed_reference_economics_20260928/utility_ces_scales_v1'
RESERVE=900.; PENALTY=1e12

def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False,default=str)+'\n');temp.replace(path)
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def simplex(seed,bounds,coordinates):
    lo=np.array([bounds[k][0] for k in coordinates]);span=np.array([bounds[k][1]-bounds[k][0] for k in coordinates]);c=(np.array([seed[k] for k in coordinates])-lo)/span
    a=np.tile(c,(len(c)+1,1));steps={k:min({'beta_annual':.002,'first_birth_fixed_cost':.06,'theta0':.02,'child_benefit_curvature':.01,'tenure_choice_kappa':.002,'psi_child':.015,'lambda_housing':.05}.get(k,max(.1*abs(seed[k]),.01)),.1*span[j]) for j,k in enumerate(coordinates)}
    for j,k in enumerate(coordinates):
        d=steps[k]/span[j];a[j+1,j]+=d if c[j]+d<=1 else -d
    return lo,span,a,steps

def controller(out,seed,bounds,evaluate_full,evaluate_fast,deadline,*,toy=False):
    """At most 20 GEs: two baseline gates, seventeen proposals, one full repeat."""
    out=Path(out);out.mkdir(parents=True,exist_ok=True);coords=tuple(seed);lo,span,initial,steps=simplex(seed,bounds,coords)
    cases=[];cache={};best=None;calls=0;stop='nm_call_cap'
    def run(label,point,kind,full=False):
        nonlocal best
        require(len(cases)<20,'20-GE cap reached')
        end=deadline if kind=='selected_repeat' else deadline-RESERVE
        require(time.time()<end,'Budget exhausted before evaluation')
        write(out/'latest.json',dict(status='running_full_GE',label=label,kind=kind,completed_full_ge=len(cases),parameters=point,deadline_epoch=deadline))
        r=(evaluate_full if full else evaluate_fast)(label,point,end);row=dict(label=label,kind=kind,parameters=point,**r);cases.append(row)
        if r['status']=='passed':
            rr=np.asarray(r['residual']);require(rr.shape==(10,) and np.isfinite(rr).all(),'Ten finite residuals required');row['loss']=float(rr@rr)
            if not toy:
                row['target_fit']=list(csv.DictReader(open(Path(row['report'])/'target_fit.csv')))
                row['effective_parameters']=list(csv.DictReader(open(Path(row['report'])/'parameters.csv')))
                require(len(row['target_fit'])==14 and len(row['effective_parameters'])==31,'Complete native tables required')
                row['free_parameter_table']=[dict(parameter=k,estimate=point[k],lower=bounds[k][0],upper=bounds[k][1],near_bound=min(point[k]-bounds[k][0],bounds[k][1]-point[k])<=.01*(bounds[k][1]-bounds[k][0])) for k in coords]
            if kind!='selected_repeat' and (best is None or row['loss']<best['loss']):best=row
        else:row['loss']=PENALTY
        write(out/'cases.json',cases);write(out/'latest_completed.json',row);write(out/'best_so_far.json',dict(status='provisional',best=best));return row
    baseline=run('000_baseline',dict(seed),'baseline',True)
    require(baseline['status']=='passed','Baseline failed; search blocked')
    repeat=run('001_baseline_repeat',dict(seed),'baseline_repeat',True)
    require(repeat['status']=='passed' and baseline['residual']==repeat['residual'],'Independent baseline repeat differs')
    if not toy:
        native=load(PACKET/'runtime/utility_adapter/runner.py','ces_repeat_checker');write(out/'baseline_repeat.json',native.compare_repeated(Path(baseline['report']),Path(repeat['report'])))
    cache[tuple(float(seed[k]).hex() for k in coords)]=baseline['loss']
    def objective(x):
        nonlocal calls
        calls+=1
        if time.time()>=deadline-RESERVE:raise TimeoutError('900-second selected verification reserve reached')
        point=dict(seed) if np.array_equal(x,initial[0]) else {k:float(lo[j]+span[j]*x[j]) for j,k in enumerate(coords)};key=tuple(point[k].hex() for k in coords)
        if key in cache:return cache[key]
        if len(cases)>=18:raise TimeoutError('16 new search proposals plus cached center cap reached')
        r=run(f'{len(cases):03d}_nm',point,'nm');cache[key]=r['loss'];return r['loss']
    write(out/'search_contract.json',dict(initial_simplex=initial.tolist(),steps=steps,coordinates=list(coords),max_proposals_including_center=17,max_GE=20,reserve_seconds=RESERVE,base_weights_unchanged=True))
    try:minimize(objective,initial[0],method='Nelder-Mead',bounds=[(0,1)]*len(coords),options=dict(initial_simplex=initial,maxfev=17,maxiter=17,adaptive=True))
    except TimeoutError as e:stop=str(e)
    selected=None;verification=None
    try:
        selected=run(f'{len(cases):03d}_selected_repeat',dict(best['parameters']),'selected_repeat',True)
        require(selected['status']=='passed' and selected['residual']==best['residual'],'Selected native repeat differs')
        if not toy:
            matcher=load(PACKET/'runtime/fast_objective.py','ces_selected_matcher')
            verification=matcher.compare_saved_baseline(best,selected['report'],atol=0.)
    except (TimeoutError,RuntimeError) as exc:
        stop='selected_verification_failed: '+str(exc)
        write(out/'selected_failure.json',dict(reason=str(exc),best=best,selected=selected))
    verified=selected is not None and selected['status']=='passed' and selected['residual']==best['residual'] and not stop.startswith('selected_verification_failed')
    result=dict(status='selected_numerically_verified' if verified else 'provisional_selected_unverified',selected=best,selected_repeat=selected,table_verification=verification,completed_full_ge=len(cases),objective_calls=calls,lifecycle_solves=sum(r.get('lifecycle_solves',0) for r in cases),stop=stop,optimization_convergence_certified=False,experimental_not_adopted=True)
    write(out/'completed.json',result);return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--initialize-only',action='store_true');ap.add_argument('--toy',action='store_true');ap.add_argument('--deadline-epoch',type=float);a=ap.parse_args()
    require(not a.out.exists(),'Refusing existing output');a.out.mkdir(parents=True);start=time.time();deadline=min(start+3600,a.deadline_epoch or start+3600)
    config=json.loads((PACKET/'plan.json').read_text());seed=config['seed'];bounds=config['bounds']
    if a.toy:
        def fake(label,p,end):return dict(status='passed',residual=[(p[k]-seed[k])/(bounds[k][1]-bounds[k][0])-.03 for k in seed],lifecycle_solves=0)
        controller(a.out,seed,bounds,fake,fake,deadline,toy=True);return
    for rel,digest in json.loads((PACKET/'source_pins.json').read_text()).items():require(hashlib.sha256((PACKET/rel).read_bytes()).hexdigest()==digest,'Private source drift '+rel)
    adapter=PACKET/'runtime/utility_adapter';sys.path.insert(0,str(adapter));native=load(adapter/'runner.py','ces_native');inputs=native.inputs
    require(native.PLAN['target_contract']==config['target_contract'] and inputs.canonical(native.PLAN['target_contract'])==config['target_weight_fingerprint'],'Pinned base target weights differ')
    lane='floor_s0';inputs.LANES[lane].update(seed=seed,bounds=bounds,free_coordinates=list(seed),arm='ces');inputs.ARMS['ces']=0.
    P,grid=inputs.proposal(lane);P,entry=inputs.entry(P,grid,'nonnegative_mean');Q=inputs.bind(P,seed,bounds,'ces')
    write(a.out/'input_contract.json',dict(entry=entry,parameters=seed,bounds=bounds,fixed_CES_eta=.487,fixed_alpha=.733,compensation_off=True,floors_zero=True,target_weight_fingerprint=config['target_weight_fingerprint'],experimental_not_adopted=True))
    (a.out/'full').mkdir()
    full=native.native_evaluator(a.out/'full',lane,Q,grid,deadline,config['price_start'])
    fast=load(PACKET/'runtime/fast_objective.py','ces_fast').make_evaluator(a.out/'fast',lane,Q,grid,deadline,config['price_start'],native_runner=native)
    write(a.out/'initializer.json',dict(status='exact_initializer_passed_zero_lifecycle',lifecycle_solves=0))
    if a.initialize_only:return
    def alarm(*_):raise TimeoutError('3600-second hard stop')
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.001,deadline-time.time()))
    try:controller(a.out/'results',seed,bounds,full,fast,deadline)
    except BaseException as exc:
        write(a.out/'failure.json',dict(type=type(exc).__name__,message=str(exc),no_auto_retry=True));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':main()
