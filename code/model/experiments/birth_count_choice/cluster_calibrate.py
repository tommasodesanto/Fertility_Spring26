"""Matched estate-A recalibration: one birth-menu cap and one matched start."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, signal, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[4]
EXPERIMENT = Path(__file__).resolve().parent
ANCHOR = ROOT / 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json'
RESERVE, MAX_CALLS, MAX_LIFECYCLE = 1800., 500, 32
PENALTY = 1e12

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def canonical(value): return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
def write(path, value):
    temporary = Path(path).with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)
def readrows(path):
    with Path(path).open(newline='') as f: return list(csv.DictReader(f))
def require(ok, message):
    if not ok: raise RuntimeError(message)
class BudgetStop(Exception): pass

def simplex(seed, bounds, coordinates):
    explicit = {'beta_annual': .002, 'h_P': .1, 'first_birth_fixed_cost': .06, 'theta0': .02,
                'child_benefit_curvature': .01, 'tenure_choice_kappa': .002, 'psi_child': .015}
    lower = np.array([bounds[k][0] for k in coordinates])
    spans = np.array([bounds[k][1]-bounds[k][0] for k in coordinates])
    center = (np.array([seed[k] for k in coordinates])-lower)/spans
    result = np.tile(center, (len(center)+1, 1))
    for j,k in enumerate(coordinates):
        step = .25*min(explicit.get(k, max(.1*abs(seed[k]), .01)), .1*spans[j])/spans[j]
        result[j+1,j] += step if center[j]+step <= 1. else -step
    return lower, spans, result

def checked_plan(path, digest):
    require(sha(path)==digest, 'Matched starts SHA-256 drift')
    plan=json.loads(Path(path).read_text()); anchor=json.loads(ANCHOR.read_text())
    require(anchor['status']=='selected_numerically_verified' and anchor['chain']==13 and anchor['arm']=='alternative', 'Chain13 identity drift')
    require(plan['source_checkpoint_sha256']==sha(ANCHOR), 'Anchor checkpoint drift')
    require(plan['starts'][0]==anchor['selected']['parameters'], 'First start differs from chain13')
    provisional=json.loads((ROOT/plan['provisional_seed_source']).read_text())
    require(sha(ROOT/plan['provisional_seed_source'])==plan['provisional_seed_source_sha256'],'Provisional seed source drift')
    best=next(c for c in provisional['chains'] if c['chain']==6)['files']['best_so_far']['best']
    require(plan['starts'][1]==best['parameters'] and best['label']=='0060_nm' and best['exploration_unverified'] is True,'Unverified provisional seed identity drift')
    require([{k:r[k] for k in ('moment','target','weight','role')} for r in best['target_fit']]==plan['target_contract'] and best['weight_fingerprint']==plan['weight_fingerprint'],'Provisional seed target identity drift')
    require(len(plan['starts'])==5 and len(plan['bounds'])==10, 'Five matched starts / ten coordinates required')
    require(plan['bounds']['beta_annual']==[.93,.99], 'Beta bounds drift')
    require(plan['arms']=={'binary':1,'count3':3}, 'Birth menu arm contract drift')
    require(len({canonical(s) for s in plan['starts']})==5,'Duplicate starts')
    for row in plan['starts']:
        require(set(row)==set(plan['bounds']), 'Start coordinate drift')
        require(all(math.isfinite(v) and lo<=v<=hi for k,v in row.items() for lo,hi in [plan['bounds'][k]]),'Start outside bounds')
    return plan,anchor

def experimental_parameter_rows(parameters):
    revised=[dict(row) for row in parameters]
    for row in revised:
        if row['parameter']=='beta_annual':
            row['lower'],row['upper']='.93','.99'
            row['near_bound']=str(min(float(row['estimate'])-.93,.99-float(row['estimate']))<=.0006)
    return revised

def check_native(result, point, contract, bounds):
    require(result['status']=='passed', 'Native selected verification failed')
    report=Path(result['report'])
    fits=readrows(report/'target_fit_new_contract.csv'); params=readrows(report/'parameters_estate_a.csv')
    require(len(fits)==14 and len(params)==31, 'Native 14/31 report drift')
    require([{k:r[k] for k in ('moment','target','weight','role')} for r in fits]==contract, 'New target contract drift')
    require(fits==result['target_fit'], 'Evaluator/report fit drift')
    for key,value in point.items():
        rows=[r for r in params if r['parameter']==key]
        require(len(rows)==1 and float(rows[0]['estimate'])==value,'Parameter report drift: '+key)
        require([float(rows[0]['lower']),float(rows[0]['upper'])]==bounds[key], 'Parameter bound report drift: '+key)
    plots={p.name:sha(p) for p in sorted((report/'standard_diagnostics').glob('*.png'))}
    require(len(plots)==17,'Standard plot count drift')
    rr=np.asarray(result['residual']); loss=float(rr@rr)
    require(rr.shape==(10,) and np.isfinite(rr).all(),'Ten residual contract drift')
    require(abs(sum(float(r['loss_contribution'] or 0) for r in fits)-loss)<1e-8,'Native loss arithmetic drift')
    return fits,params,plots,loss

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--arm',choices=('binary','count3'),required=True);p.add_argument('--chain',type=int,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--deadline-epoch',type=float,required=True)
    p.add_argument('--starts-file',type=Path,required=True);p.add_argument('--starts-file-sha256',required=True)
    p.add_argument('--mock-smoke',action='store_true');p.add_argument('--preflight-evaluator',action='store_true')
    p.add_argument('--smoke',action='store_true');p.add_argument('--postcheck-only',action='store_true');p.add_argument('--search-receipt',type=Path)
    args=p.parse_args();require(args.postcheck_only==(args.search_receipt is not None),'Postcheck requires search receipt')
    require(not args.out.exists(),'Refusing existing output directory');args.out.mkdir(parents=True);out=args.out.resolve()
    start=time.time();deadline=min(args.deadline_epoch,start+21600)
    require(deadline>start+(0 if args.mock_smoke or args.postcheck_only else RESERVE),'Insufficient final native reserve')
    plan,anchor=checked_plan(args.starts_file,args.starts_file_sha256)
    require(0<=args.chain<5,'Matched chain outside 0..4');cap=plan['arms'][args.arm];seed=plan['starts'][args.chain]
    bounds=plan['bounds'];coordinates=tuple(seed);lower,spans,initial=simplex(seed,bounds,coordinates)
    contract=dict(arm=args.arm,birth_cap=cap,chain=args.chain,seed=seed,all_starts=plan['starts'],starts_count=5,
        bounds=bounds,free_coordinates=list(coordinates),starts_file_sha256=args.starts_file_sha256,
        selected_source_sha256=sha(ANCHOR),provisional_seed_source_sha256=plan['provisional_seed_source_sha256'],provisional_seed=plan['provisional_seed'],target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'],
        experiment_flags=dict(birth_count_choice_enabled=True,birth_count_choice_cap=cap,bequest_net_of_selling_cost=True,estate_flow_net_of_selling_cost=True),
        estate='bp + (1-psi)*P*h_prime; no extra R',scf_scope='provisional recipient/scope mismatch; retained target and weight',
        no_auto_retry=True,not_adopted_calibration=True,search_evaluator='full_native',selected_evaluator='full_native')
    write(out/'start_contract.json',contract)
    limit=2 if args.smoke or args.mock_smoke else MAX_CALLS
    write(out/'search_contract.json',dict(method='bounded Nelder-Mead',initial_simplex=initial.tolist(),max_objective_calls=limit,
        max_lifecycle_per_full_GE=32,final_reserve_seconds=RESERVE,wall_seconds=21600,deterministic_seed=plan['deterministic_seed']))
    def heartbeat(status,**more):write(out/'heartbeat.json',dict(epoch=time.time(),status=status,arm=args.arm,chain=args.chain,**more))
    heartbeat('initialized',objective_calls=0,completed_full_ge=0)
    if not args.mock_smoke:
        sys.path.insert(0,str(EXPERIMENT));from model import estate_contract as estate
        from model.inputs import load_inputs, DEFAULT_PRICE
        targets,targetpin,weightpin,livebounds=estate.contract()
        require(targets==plan['target_contract'] and targetpin==plan['target_fingerprint'] and weightpin==plan['weight_fingerprint'],'Estate target fingerprint drift')
        require({k:list(v) for k,v in livebounds.items()}==bounds,'Estate bounds drift')
        P,grid=load_inputs(parameters=seed)
        evaluate=estate.make_evaluator(out,'floor_s0',P,grid,deadline,anchor['selected']['price'],birth_cap=cap,
            target_fingerprint=targetpin,weight_fingerprint=weightpin)
        if args.preflight_evaluator:
            # Initialize the actual authenticated reporting stack, including frozen overlays.
            # This is the same caller P/grid/flags used by the evaluator, with no solve.
            estate.apply_experiment_flags(P,estate.experiment_flags(cap))
            from model.reporting import build_context
            context_out=out/'reporting_context';context_out.mkdir()
            context=build_context(P,grid,context_out,price_start=anchor['selected']['price'],deadline=deadline,max_lifecycle=MAX_LIFECYCLE,closure='population_one')
            require(all(getattr(context['P'],k)==v for k,v in estate.experiment_flags(cap).items()),'Context experiment flag drift')
            require(np.array_equal(context['b_grid'],grid),'Context grid drift')
            write(out/'completed.json',dict(status='evaluator_initialized_zero_solves',lifecycle_solves=0,actual_reporting_context_built=True,context_flags=estate.experiment_flags(cap),source_contract=contract))
            heartbeat('evaluator_initialized_zero_solves',objective_calls=0,completed_full_ge=0);return
        if args.postcheck_only:
            search=json.loads(args.search_receipt.read_text())
            require(search['status']=='provisional_search_finished' and all(search[k]==contract[k] for k in ('arm','chain','starts_file_sha256','target_fingerprint','weight_fingerprint')),'Postcheck search drift')
            chosen=search['selected'];require(chosen is not None and chosen['status']=='passed','Missing passed selected candidate')
            verify=evaluate('selected_postcheck',chosen['parameters'],deadline)
            fits,params,plots,loss=check_native(verify,chosen['parameters'],targets,bounds)
            report=Path(verify['report']);repeat_report=report.parent/'selected_repeat_final'
            require(repeat_report.is_dir(),'Missing native exact repeat report')
            # Native GE performs its exact solution-array repeat and acceptance gates internally.
            oldfits=readrows(repeat_report/'target_fit.csv')
            repeatfits,repeatrr,repeatparams,_=estate.rescore_report(repeat_report)
            repeatparams=experimental_parameter_rows(repeatparams)
            require(repeatfits==fits and np.array_equal(np.asarray(repeatrr),np.asarray(verify['residual'])),'Exact repeat experimental rows drift')
            repeatplots={z.name:sha(z) for z in sorted((repeat_report/'standard_diagnostics').glob('*.png'))}
            require(repeatparams==params and repeatplots==plots,'Exact repeat parameter / plots drift')
            write(out/'completed.json',dict(status='full_native_postcheck_passed',selected_postcheck=verify,native_loss=loss,
                target_fit=fits,parameters=params,repeat=dict(status='exact_full_ge_repeat_passed',target_rows=14,parameter_rows=31,
                standard_plot_hashes=plots,experimental_target_fit_exact=True),search_receipt_sha256=sha(args.search_receipt),
                starts_file_sha256=args.starts_file_sha256,target_fingerprint=targetpin,weight_fingerprint=weightpin));return
    cases=[];cache={};best=None;calls=0
    def objective(x):
        nonlocal best,calls
        if time.time()>=deadline-RESERVE:raise BudgetStop('final_native_reserve_reached')
        if calls>=limit:raise BudgetStop('objective_call_limit')
        if not args.mock_smoke and os.statvfs(out).f_bavail*os.statvfs(out).f_frsize<350*1024**3:raise BudgetStop('shared_free_disk_below_350GiB')
        calls+=1;point={k:float(lower[j]+spans[j]*x[j]) for j,k in enumerate(coordinates)}
        key=tuple(point[k].hex() for k in coordinates)
        if key in cache:heartbeat('cache_hit',objective_calls=calls,completed_full_ge=len(cases));return cache[key]
        label=f'{len(cases):04d}_nm';heartbeat('running_full_GE',label=label,objective_calls=calls,completed_full_ge=len(cases))
        if args.mock_smoke:
            rr=np.asarray(x)-initial[0]-.01
            result=dict(status='passed',residual=rr.tolist(),lifecycle_solves=0,mock=True)
        else:result=evaluate(label,point,deadline-RESERVE)
        row=dict(label=label,parameters=point,**result)
        if result['status']=='passed':
            rr=np.asarray(result['residual']);require(rr.shape==(10,) and np.isfinite(rr).all(),'Residual contract drift')
            row['loss']=row['objective']=float(rr@rr);row['weight_fingerprint']=plan['weight_fingerprint']
            if best is None or row['loss']<best['loss']:best=row
        elif result['status']=='inadmissible_numerical':row.update(objective=PENALTY,numerical_rejection=True,computed_valid_loss=False)
        elif result['status']=='budget_exhausted':
            write(out/'latest_completed.json',dict(latest=row,completed_full_ge=len(cases),objective_calls=calls));raise BudgetStop('native_evaluation_budget_exhausted')
        else:raise RuntimeError('Unexpected evaluator status: '+str(result['status']))
        case_bytes=sum(z.stat().st_size for z in (out/label).rglob('*') if z.is_file())
        if not args.mock_smoke and case_bytes>1024**3:raise BudgetStop('case_storage_exceeded_1GiB_planning_cap')
        row['retained_case_bytes']=case_bytes
        cases.append(row);cache[key]=row['objective']
        write(out/'cases.json',cases);write(out/'latest_completed.json',dict(latest=row,completed_full_ge=len(cases),objective_calls=calls))
        write(out/'best_so_far.json',dict(status='provisional_until_fresh_native_postcheck',best=best,completed_full_ge=len(cases)))
        heartbeat('case_completed',objective_calls=calls,completed_full_ge=len(cases),best_loss=best['loss'] if best else None)
        return row['objective']
    stop='optimizer_return'
    try:
        try:minimize(objective,initial[0],method='Nelder-Mead',bounds=[(0.,1.)]*10,options=dict(initial_simplex=initial,maxfev=limit,maxiter=limit,xatol=1e-4,fatol=1e-4,adaptive=True))
        except BudgetStop as exc:stop=str(exc)
        search=dict(status='provisional_search_finished',search_stop_reason=stop,arm=args.arm,birth_cap=cap,chain=args.chain,selected=best,
            objective_calls=calls,completed_full_ge=len(cases),lifecycle_solves=sum(r.get('lifecycle_solves',0) for r in cases),
            target_fingerprint=plan['target_fingerprint'],weight_fingerprint=plan['weight_fingerprint'],starts_file_sha256=args.starts_file_sha256,
            search_evaluator='full_native',selected_evaluator='full_native',optimization_convergence_certified=False,no_auto_retry=True)
        write(out/'search_completed.json',search)
        if args.mock_smoke:
            require(calls==2 and len(cases)==2,'Mock loop did not perform exact two cases')
            write(out/'completed.json',dict(search,status='mock_loop_passed_zero_solves'));heartbeat('mock_completed',objective_calls=calls,completed_full_ge=len(cases));return
        if best is None:write(out/'completed.json',dict(search,status='no_admissible_candidate'));return
        heartbeat('native_postcheck_running',objective_calls=calls,completed_full_ge=len(cases))
        command=[sys.executable,str(Path(__file__).resolve()),'--arm',args.arm,'--chain',str(args.chain),'--out',str(out/'native_postcheck'),
            '--deadline-epoch',str(deadline),'--starts-file',str(args.starts_file.resolve()),'--starts-file-sha256',args.starts_file_sha256,
            '--postcheck-only','--search-receipt',str(out/'search_completed.json')]
        child=subprocess.run(command,capture_output=True,text=True,timeout=min(RESERVE,deadline-time.time()),check=False)
        (out/'postcheck_child.stdout.log').write_text(child.stdout);(out/'postcheck_child.stderr.log').write_text(child.stderr)
        require(child.returncode==0,'Full native postcheck child failed: '+str(child.returncode))
        verified=json.loads((out/'native_postcheck/completed.json').read_text())
        require(verified['status']=='full_native_postcheck_passed' and verified['search_receipt_sha256']==sha(out/'search_completed.json'),'Native child receipt drift')
        require(all(verified[k]==search[k] for k in ('starts_file_sha256','target_fingerprint','weight_fingerprint')),'Native child fingerprints drift')
        require(np.allclose(best['residual'],verified['selected_postcheck']['residual'],rtol=0,atol=1e-10),'Fresh native selected residual drift')
        done=dict(search,status='selected_numerically_verified',elapsed_seconds=time.time()-start)
        done.update({k:verified[k] for k in ('native_loss','selected_postcheck','target_fit','parameters','repeat')})
        done['smoke_fast_full_comparison']=dict(status='search_full_new_target_exact',target_rows=14,nonwealth_rows_unchanged=True) if args.smoke else None
        write(out/'completed.json',done);heartbeat('completed',objective_calls=calls,completed_full_ge=len(cases),native_loss=done['native_loss'])
    except BaseException as exc:
        write(out/'failure.json',dict(type=type(exc).__name__,message=str(exc),elapsed_seconds=time.time()-start,no_auto_retry=True))
        heartbeat('failed',objective_calls=calls,completed_full_ge=len(cases),error=str(exc));raise

if __name__=='__main__':main()
