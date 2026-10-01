"""Bounded independent Nelder--Mead on the unchanged nine-coordinate floor contract."""
from __future__ import annotations
import argparse, importlib.util, json, math, os, shutil, signal, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
OLD = HERE.parent/'utility_floor_round2_v1'
sys.path.insert(0, str(OLD/'local_run'))
sys.path.insert(0, str(OLD))
import runner_local as native
inputs = native.inputs
write = native.write
SEEDS = (0, 1, 5, 6, 3, 7)
MAX_EVAL = 500
RESERVE = 900.
PENALTY = 1e12

class BudgetStop(Exception):
    pass

def steps(seed, bounds, coordinates):
    explicit = {'beta_annual': .002, 'h_P': .1, 'first_birth_fixed_cost': .06,
                'theta0': .02, 'child_benefit_curvature': .01, 'tenure_choice_kappa': .002}
    return {k: min(explicit.get(k, max(.1*abs(seed[k]), .01)), .1*(bounds[k][1]-bounds[k][0])) for k in coordinates}

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

def optimize(out, seed, bounds, coordinates, evaluate, deadline, *, toy=False, maxeval=MAX_EVAL):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    lower,spans,initial=simplex(seed,bounds,coordinates)
    cache={};cases=[];best=None;started=time.time();calls=0;stop=None
    write(out/'search_contract.json', dict(method='scipy Nelder-Mead',scaled_bounds=[0.,1.],free_coordinates=list(coordinates),seed=seed,bounds=bounds,initial_simplex=initial.tolist(),physical_simplex_steps=steps(seed,bounds,coordinates),maximum_objective_calls=maxeval,final_reserve_seconds=RESERVE,invalid_numerical_root_penalty=PENALTY,no_extra_baseline_gate=True,exact_vector_cache=True,targets_scored=10,targets_total=14))
    def objective(x):
        nonlocal best,calls
        if time.time() >= deadline-RESERVE: raise BudgetStop('four_hour_budget_final_reserve')
        if shutil.disk_usage(out).free < 20*1024**3: raise BudgetStop('free_disk_below_20GiB')
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
            rr=np.asarray(row['residual'],float)
            inputs.require(rr.shape==(10,) and np.isfinite(rr).all(),'Invalid scored residual dimension')
            row['objective']=row['loss']=float(rr@rr)
            if not toy:
                row['target_fit']=native.readtable(Path(row['report'])/'target_fit.csv')
                row['effective_parameters']=native.readtable(Path(row['report'])/'parameters.csv')
                inputs.require(np.array_equal(native.residual(row['target_fit']),rr),'Checkpoint residual differs from full table')
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
        result=minimize(objective,initial[0],method='Nelder-Mead',bounds=[(0.,1.)]*len(coordinates),options=dict(initial_simplex=initial,maxfev=maxeval,maxiter=maxeval,xatol=1e-4,fatol=1e-4,adaptive=True))
        stop='optimizer_success' if result.success else 'optimizer_budget_or_stop'
    except BudgetStop as exc: stop=str(exc)
    receipt=dict(status='provisional_search_finished',search_stop_reason=stop,completed_full_ge=len(cases),objective_calls=calls,cache_entries=len(cache),selected=best,lifecycle_solves=sum(r.get('lifecycle_solves',0) for r in cases),optimizer_success=bool(result.success) if result is not None else False,optimization_convergence_certified=False,experimental_not_adopted=True,elapsed_seconds=time.time()-started)
    write(out/'search_completed.json',receipt)
    return receipt

def alarm(*_): raise TimeoutError('Four-hour hard stop')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--chain',type=int,choices=range(6),required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--deadline-epoch',type=float,required=True);ap.add_argument('--toy',action='store_true');ap.add_argument('--fast-objective',action='store_true');ap.add_argument('--verify-only',type=Path);args=ap.parse_args()
    inputs.require(not args.out.exists(),'Refusing existing results directory');args.out.mkdir(parents=True)
    lane=f'floor_s{SEEDS[args.chain]}';seed,bounds,_=inputs.seed_and_bounds(lane);coordinates=inputs.parameters(lane);start=time.time();deadline=min(args.deadline_epoch,start+14400)
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,max(.001,deadline-start))
    try:
        if args.toy:
            center=np.array([seed[k] for k in coordinates]);spans=np.array([bounds[k][1]-bounds[k][0] for k in coordinates])
            def evaluate(label,point,end):
                r=np.r_[(np.array([point[k] for k in coordinates])-center)/spans-.03,.05]
                return dict(status='passed',residual=r.tolist(),lifecycle_solves=0,report=str(args.out/label))
            result=optimize(args.out,seed,bounds,coordinates,evaluate,deadline,toy=True,maxeval=40)
            inputs.require(result['completed_full_ge']>10 and result['lifecycle_solves']==0,'Toy did not exercise NM loop')
            write(args.out/'completed.json',result);return
        native.verify_sources()
        P,grid=inputs.proposal(lane);P,entry=inputs.entry(P,grid,'nonnegative_mean');Q=native.utility_checks(P,grid,lane,args.out)
        receipt=native.validate_smoke_receipt(OLD/'local_run/fresh_gate/results/smoke_receipt.json','floor')
        write(args.out/'input_contract.json',dict(lane=lane,entry=entry,seed=seed,bounds=bounds,target_contract=native.PLAN['target_contract'],fixed_H0=float(P.H0[0]),fixed_psi=float(P.psi_child),economic_changes=native.PLAN['economic_changes']['floor'],existing_verified_smoke_reused=True,smoke_receipt_sha256=inputs.sha(OLD/'local_run/fresh_gate/results/smoke_receipt.json'),no_borrowing=True,experimental_not_adopted=True))
        if args.verify_only:
            search_result=json.loads(args.verify_only.read_text());selected=search_result['selected']
            inputs.require(selected is not None,'No computed selected point for postcheck')
            evaluate=native.native_evaluator(args.out,lane,Q,grid,deadline,receipt['selected_price'])
            verification=evaluate('selected_postcheck',selected['parameters'],deadline)
            if verification['status']=='passed':
                inputs.require(verification['residual']==selected['residual'],'Selected fresh-process postcheck residual differs')
            write(args.out/'completed.json',dict(status='selected_numerically_verified' if verification['status']=='passed' else 'provisional_postcheck_uncomputed',selected=selected,selected_postcheck=verification,search_result=str(args.verify_only),deadline_epoch=deadline,elapsed_seconds=time.time()-start,optimization_convergence_certified=False,no_auto_retry=True))
            return
        factory=native.native_evaluator
        if args.fast_objective:
            import fast_objective
            factory=lambda *a,**k: fast_objective.make_evaluator(*a,**k,native_runner=native)
        evaluate=factory(args.out,lane,Q,grid,deadline,receipt['selected_price'])
        result=optimize(args.out,seed,bounds,coordinates,evaluate,deadline)
        result.update(selected_postcheck='pending_fresh_supervised_process',elapsed_seconds=time.time()-start,deadline_epoch=deadline,no_auto_retry=True)
        write(args.out/'completed.json',result)
    except BaseException as exc:
        write(args.out/'failure.json',dict(type=type(exc).__name__,message=str(exc),elapsed_seconds=time.time()-start,no_auto_retry=True));raise
    finally: signal.setitimer(signal.ITIMER_REAL,0)

if __name__=='__main__': main()
