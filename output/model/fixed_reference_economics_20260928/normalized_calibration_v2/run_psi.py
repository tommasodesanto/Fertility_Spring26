"""Bounded independent Nelder--Mead on the ten-coordinate free-psi floor contract."""
from __future__ import annotations
import argparse, importlib.util, json, math, os, shutil, signal, sys, time
from pathlib import Path
import numpy as np
from types import SimpleNamespace
from pso import particle_swarm
import normalized_objective
from scipy.optimize import minimize
if sys.platform=='darwin':
    import importlib,pathlib
    sys.modules.setdefault('pathlib._local',pathlib)
    sys.modules.setdefault('numpy._core',np.core)
    sys.modules.setdefault('numpy._core.multiarray',importlib.import_module('numpy.core.multiarray'))
    sys.modules.setdefault('numpy._core.numeric',importlib.import_module('numpy.core.numeric'))

HERE = Path(__file__).resolve().parent
OLD = HERE.parent/'utility_floor_round2_v1'
sys.path.insert(0, str(OLD))
import runner as native
inputs = native.inputs
write = native.write
SEEDS = (0,)*24
MAX_EVAL = 150
RESERVE = 900.
PENALTY = 1e12
CONFIG=json.loads((HERE/'plan.json').read_text())
SCORED=[r['moment'] for r in CONFIG['base_target_contract'] if r['role']=='scored']
def weight_fingerprint(multipliers):
    return inputs.canonical(dict(base_contract=CONFIG['base_target_contract'],multipliers=multipliers))
def weighted_residual(base,multipliers):
    inputs.require(len(base)==len(SCORED)==10,'Ten base scored residuals required')
    return np.asarray(base,float)*np.sqrt([multipliers.get(k,1.) for k in SCORED])
def weighted_rows(rows,multipliers):
    inputs.require(native.target_identity(rows)==CONFIG['base_target_contract'],'Base target contract drift')
    result=[]
    for r in rows:
        x=dict(r);m=float(multipliers.get(r['moment'],1.));x['base_weight']=r['weight'];x['weight_multiplier']=str(m);x['base_loss_contribution']=r['loss_contribution']
        if r['weight']!='': x['weight']=str(float(r['weight'])*m)
        if r['loss_contribution']!='': x['loss_contribution']=str(float(r['loss_contribution'])*m)
        result.append(x)
    return result


class BudgetStop(Exception):
    pass

def steps(seed, bounds, coordinates):
    explicit = {'beta_annual': .002, 'h_P': .1, 'first_birth_fixed_cost': .06,
                'theta0': .02, 'child_benefit_curvature': .01, 'tenure_choice_kappa': .002, 'psi_child': .015}
    return {k: .25*min(explicit.get(k, max(.1*abs(seed[k]), .01)), .1*(bounds[k][1]-bounds[k][0])) for k in coordinates}

def simplex(seed, bounds, coordinates):
    lower = np.array([bounds[k][0] for k in coordinates], dtype=float)
    spans = np.array([bounds[k][1]-bounds[k][0] for k in coordinates], dtype=float)
    center = (np.array([seed[k] for k in coordinates])-lower)/spans
    result = np.tile(center, (len(center)+1, 1))
    for j,k in enumerate(coordinates):
        delta = steps(seed,bounds,coordinates)[k]/spans[j]
        result[j+1,j] += delta if center[j]+delta <= 1. else -delta
    return lower, spans, result

def compact_case(directory, keep_plots):
    """Prune only newly produced case payloads after compact receipts are parsed."""
    directory = Path(directory)
    if not directory.exists(): return
    for p in list(directory.rglob('*')):
        if p.is_file() and p.suffix.lower() in ('.npz','.npy','.pkl','.pickle'):
            p.unlink()
    if not keep_plots:
        for folder in list(directory.rglob('standard_diagnostics')):
            for p in folder.glob('*.png'): p.unlink()
    # Native price caches are task-owned case directories; preserve their JSON receipts.
    for folder in list(directory.rglob('*')):
        if folder.is_dir() and folder.name in ('exact_policy_cache','policy_cache'):
            shutil.rmtree(folder)

def optimize(out, seed, bounds, coordinates, evaluate, deadline, *, toy=False, maxeval=MAX_EVAL, multipliers=None, chain=0):
    multipliers=multipliers or {}
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    lower,spans,initial=simplex(seed,bounds,coordinates)
    cache={};cases=[];best=None;started=time.time();calls=0;stop=None
    write(out/'search_contract.json', dict(method='particle_swarm' if False else 'scipy Nelder-Mead',scaled_bounds=[0.,1.],free_coordinates=list(coordinates),seed=seed,bounds=bounds,initial_simplex=initial.tolist(),physical_simplex_steps=steps(seed,bounds,coordinates),maximum_objective_calls=maxeval,final_reserve_seconds=RESERVE,invalid_numerical_root_penalty=PENALTY,no_extra_baseline_gate=True,exact_vector_cache=True,targets_scored=10,targets_total=14))
    def objective(x):
        nonlocal best,calls
        if time.time() >= deadline-RESERVE: raise BudgetStop('three_hour_actual_start_final_reserve')
        if shutil.disk_usage(out).free < 20*1024**3: raise BudgetStop('free_disk_below_20GiB')
        if calls>=maxeval: raise BudgetStop('maximum_objective_calls')
        calls+=1
        point={k:float(lower[j]+spans[j]*x[j]) for j,k in enumerate(coordinates)}
        inputs.check_point(point,bounds)
        key=tuple(point[k].hex() for k in coordinates)
        if key in cache:
            write(out/'latest.json',dict(status='exact_vector_cache_hit',objective_calls=calls,completed_full_ge=len(cases),loss=cache[key]['objective'],parameters=point))
            return cache[key]['objective']
        if len(cases)>=maxeval: raise BudgetStop('maximum_full_ge')
        label=f'{len(cases):04d}_nm'
        write(out/'latest.json',dict(status='running_full_GE',label=label,parameters=point,objective_calls=calls,completed_full_ge=len(cases),deadline_epoch=deadline,experimental_not_adopted=True))
        t=time.time();result=evaluate(label,point,deadline-RESERVE)
        row=dict(label=label,parameters=point,seconds=time.time()-t,**result)
        if result['status']=='passed':
            base_rr=np.asarray(row['residual'],float)
            rr=weighted_residual(base_rr,multipliers)
            row['base_residual']=base_rr.tolist();row['base_loss']=float(base_rr@base_rr);row['residual']=rr.tolist();row['weight_fingerprint']=weight_fingerprint(multipliers)
            inputs.require(rr.shape==(10,) and np.isfinite(rr).all(),'Invalid scored residual dimension')
            row['objective']=row['loss']=float(rr@rr)
            if not toy:
                row['target_fit']=native.readtable(Path(row['report'])/'target_fit.csv')
                row['effective_parameters']=native.readtable(Path(row['report'])/'parameters.csv')
                inputs.require(np.array_equal(native.residual(row['target_fit']),base_rr),'Checkpoint base residual differs from full table')
                row['weighted_target_fit']=weighted_rows(row['target_fit'],multipliers)
                native.table(Path(row['report'])/'target_fit_experimental.csv',row['weighted_target_fit'])
            new_best=best is None or row['loss']<best['loss']
            if new_best:
                if best is not None and not toy: compact_case(out/best['label'],False)
                best=row
        elif result['status']=='inadmissible_numerical':
            row.update(objective=PENALTY,numerical_rejection=True,computed_valid_loss=False)
            new_best=False
            write(out/label/'numerical_rejection.json',row)
        elif result['status']=='budget_exhausted':
            cases.append(row);write(out/'cases.json',cases)
            raise BudgetStop('native_evaluation_budget_exhausted')
        else:
            raise RuntimeError('Unrecognized evaluation status: '+result['status'])
        cache[key]=row;cases.append(row)
        write(out/'latest_completed.json',dict(latest=row,completed_full_ge=len(cases),objective_calls=calls,elapsed_seconds=time.time()-started))
        write(out/'best_so_far.json',dict(status='provisional_until_selected_postcheck',best=best))
        write(out/'cases.json',cases)
        if not toy: compact_case(out/label,new_best)
        return row['objective']
    result=None
    try:
        if False:
            result=particle_swarm(objective, lower, spans, bounds, coordinates, CONFIG, chain, out, maxeval)
        else:
            result=minimize(objective,initial[0],method='Nelder-Mead',bounds=[(0.,1.)]*len(coordinates),options=dict(initial_simplex=initial,maxfev=maxeval,maxiter=maxeval,xatol=1e-4,fatol=1e-4,adaptive=True))
        stop='optimizer_success' if result.success else 'optimizer_budget_or_stop'
    except BudgetStop as exc: stop=str(exc)
    receipt=dict(status='provisional_search_finished',search_stop_reason=stop,completed_full_ge=len(cases),objective_calls=calls,cache_entries=len(cache),selected=best,lifecycle_solves=sum(r.get('lifecycle_solves',0) for r in cases),optimizer_success=bool(result.success) if result is not None else False,optimization_convergence_certified=False,experimental_not_adopted=True,elapsed_seconds=time.time()-started)
    if result is not None and hasattr(result,'final_simplex') and np.isfinite(result.final_simplex[0]).all() and np.isfinite(result.final_simplex[1]).all():
        receipt['final_simplex_scaled']=np.asarray(result.final_simplex[0],float).tolist()
        receipt['final_simplex_objective']=np.asarray(result.final_simplex[1],float).tolist()
        receipt['final_simplex_physical']=[{k:float(lower[j]+spans[j]*x[j]) for j,k in enumerate(coordinates)} for x in result.final_simplex[0]]
        receipt['final_simplex_saved_only_on_normal_scipy_return']=True
    elif result is not None and hasattr(result,'final_simplex'):
        receipt['final_simplex_not_saved_reason']='incomplete nonfinite simplex objectives'
    write(out/'search_completed.json',receipt)
    return receipt


def compare_incumbent(result):
    report=Path(result['report']);checks={}
    for name,keys in [('target_fit.csv',('target','model','gap','weight','loss_contribution')),('parameters.csv',('estimate',))]:
        actual=native.readtable(report/name);reference=native.readtable(HERE/'center'/name)
        inputs.require(len(actual)==len(reference)==(14 if name.startswith('target') else 31),'Incumbent comparison row count')
        identity='moment' if name.startswith('target') else 'parameter'
        inputs.require([r[identity] for r in actual]==[r[identity] for r in reference],'Incumbent comparison identity')
        maximum=0.
        for a,b in zip(actual,reference):
            for k in keys:
                if a[k]==b[k]=='':continue
                error=abs(float(a[k])-float(b[k]));maximum=max(maximum,error)
                inputs.require(error<=1e-10,'Incumbent economic value drift '+a[identity]+':'+k)
        checks[name]=dict(rows=len(actual),maximum_absolute_error=maximum,all_parameter_estimates_required=True)
    return checks

def alarm(*_): raise TimeoutError('Three-hour actual-start hard stop')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--chain',type=int,choices=range(24),required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--deadline-epoch',type=float,required=True);ap.add_argument('--toy',action='store_true');ap.add_argument('--initialize-only',action='store_true');ap.add_argument('--fast-objective',action='store_true');ap.add_argument('--verify-only',type=Path);ap.add_argument('--smoke-only',action='store_true');args=ap.parse_args()
    inputs.require(not args.out.exists(),'Refusing existing results directory');args.out.mkdir(parents=True)
    if not args.toy:
        native.verify_sources()
        for rel,digest in json.loads((HERE/'source_pins.json').read_text()).items():
            inputs.require(inputs.sha(native.ROOT/rel)==digest,'Multistart source drift: '+rel)
    lane=f'floor_s{SEEDS[args.chain]}';seed,bounds,_=inputs.seed_and_bounds(lane);design=CONFIG['nearby_starts'][args.chain];seed=dict(design['parameters']);profile=design['profile'];multipliers=CONFIG['profiles'][profile];bounds=dict(bounds);bounds['psi_child']=tuple(CONFIG['psi_bounds']);coordinates=tuple(inputs.parameters(lane))+('psi_child',);parameter_contract=dict(free_coordinates=list(coordinates),bounds=bounds,seed=seed,psi_bound_role=CONFIG['psi_bound_role'],target_contract_sha256=inputs.canonical(CONFIG['base_target_contract']),closure='price clears birth renewal; derived H0 clears housing at actual N0=1',H0_fixed=False,normalized_population=1.,derived_H0_bounds=[.2,80.]);parameter_fingerprint=inputs.canonical(parameter_contract);inputs.LANES[lane].update(seed=seed,bounds=bounds,free_coordinates=list(coordinates));start=time.time();deadline=min(args.deadline_epoch,start+10800)
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.001,deadline-start))
    try:
        if args.toy:
            center=np.array([seed[k] for k in coordinates]);spans=np.array([bounds[k][1]-bounds[k][0] for k in coordinates])
            def evaluate(label,point,end):
                r=(np.array([point[k] for k in coordinates])-center)/spans-.03
                return dict(status='passed',residual=r.tolist(),lifecycle_solves=0,report=str(args.out/label))
            result=optimize(args.out,seed,bounds,coordinates,evaluate,deadline,toy=True,maxeval=40,multipliers=multipliers,chain=args.chain)
            inputs.require(result['completed_full_ge']>10 and result['lifecycle_solves']==0,'Toy did not exercise NM loop')
            write(args.out/'completed.json',result);return
        inputs.require(inputs.canonical(parameter_contract)==parameter_fingerprint,'Extended parameter contract drift')
        P,grid=inputs.proposal(lane);P,entry=inputs.entry(P,grid,'nonnegative_mean');Q=native.utility_checks(P,grid,lane,args.out)
        receipt={'selected_price':float(design['initial_price'])}
        write(args.out/'input_contract.json',dict(lane=lane,entry=entry,seed=seed,bounds=bounds,target_contract=native.PLAN['target_contract'],auxiliary_trial_H0=float(P.H0[0]),H0_status='derived calibrated coefficient at N0=1',normalized_population=1.,initial_psi=float(seed['psi_child']),psi_status='free in experimental search',economic_changes=[native.PLAN['economic_changes']['floor']]+['Author-requested psi_child estimated jointly with nine floor coordinates; diagnostic bounds [.01,.5]', 'Weight profile '+profile+' multipliers '+json.dumps(multipliers,sort_keys=True)],existing_verified_adapter_reused=True,free_psi=True,psi_bounds=CONFIG['psi_bounds'],psi_bound_role=CONFIG['psi_bound_role'],target_rank_unverified=True,weight_profile=profile,weight_multipliers=multipliers,weight_contract_sha256=weight_fingerprint(multipliers),parameter_contract=parameter_contract,parameter_contract_sha256=parameter_fingerprint,base_target_contract_sha256=inputs.canonical(CONFIG['base_target_contract']),no_borrowing=True,experimental_not_adopted=True))
        if args.verify_only:
            search_result=json.loads(args.verify_only.read_text());selected=search_result['selected']
            inputs.require(selected is not None,'No computed selected point for postcheck')
            evaluate=normalized_objective.make_evaluator(args.out,lane,Q,grid,deadline,receipt['selected_price'],native_runner=native)
            verification=evaluate('selected_postcheck',selected['parameters'],deadline)
            if verification['status']=='passed':
                inputs.require(verification['residual']==selected['base_residual'],'Selected fresh-process base residual differs')
                base_rows=native.readtable(Path(verification['report'])/'target_fit.csv');rows=weighted_rows(base_rows,multipliers);native.table(Path(verification['report'])/'target_fit_experimental.csv',rows)
                rr=weighted_residual(np.asarray(verification['residual']),multipliers);inputs.require(rr.tolist()==selected['residual'],'Selected weighted residual differs');verification.update(base_loss=float(np.asarray(verification['residual'])@np.asarray(verification['residual'])),weighted_loss=float(rr@rr),weighted_target_fit=rows,base_target_fit=base_rows,weight_fingerprint=weight_fingerprint(multipliers))
            write(args.out/'completed.json',dict(status='selected_numerically_verified' if verification['status']=='passed' else 'provisional_postcheck_uncomputed',selected=selected,selected_postcheck=verification,search_result=str(args.verify_only),deadline_epoch=deadline,elapsed_seconds=time.time()-start,optimization_convergence_certified=False,no_auto_retry=True))
            return
        evaluate=normalized_objective.make_evaluator(args.out,lane,Q,grid,deadline,receipt['selected_price'],native_runner=native,exploratory=args.fast_objective)
        if args.smoke_only:
            inputs.require(not args.fast_objective,'Smoke requires full native repeat')
            verification=evaluate('normalized_incumbent_smoke',seed,deadline)
            inputs.require(verification['status']=='passed','Normalized incumbent native smoke failed')
            comparisons=compare_incumbent(verification)
            write(args.out/'completed.json',dict(status='normalized_incumbent_native_smoke_passed',result=verification,comparison=comparisons,lifecycle_solves=verification['lifecycle_solves']))
            return
        if args.initialize_only:
            write(args.out/'completed.json',dict(status='exact_initializer_passed_zero_lifecycle',parameter_contract_sha256=parameter_fingerprint,free_coordinates=list(coordinates),lifecycle_solves=0));return
        result=optimize(args.out,seed,bounds,coordinates,evaluate,deadline,multipliers=multipliers,chain=args.chain)
        result.update(selected_postcheck='pending_fresh_supervised_process',elapsed_seconds=time.time()-start,deadline_epoch=deadline,no_auto_retry=True)
        write(args.out/'completed.json',result)
    except BaseException as exc:
        write(args.out/'failure.json',dict(type=type(exc).__name__,message=str(exc),elapsed_seconds=time.time()-start,no_auto_retry=True));raise
    finally: signal.setitimer(signal.ITIMER_REAL,0)

if __name__=='__main__': main()
