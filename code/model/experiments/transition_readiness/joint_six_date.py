#!/usr/bin/env python3
"""Bounded six-date joint root using an authenticated historical measured seed.

No endpoint solves, fabricated slopes, long horizons, fit or reference changes.
The unchanged two-block root gets at most one Newton step and a qualifying
fresh replay. Root certification and terminal convergence are separate objects.
"""
import argparse
import ast
import copy
import gzip
import json
import importlib.util
import math
import os
from pathlib import Path
import pickle
import signal
import threading
import time

import legacy_changed_psi as frozen
from selected_adapter import require, sha

PSI = frozen.PSI * 1.001
CASE_SECONDS = 600
RENDER_RESERVE = 120


def pinned(item):
    require(isinstance(item, dict) and set(item) == {'path', 'sha256'}, 'Exact artifact pin required')
    path = Path(item['path'])
    if not path.is_absolute():
        path = frozen.ROOT / path
    require(path.is_file() and sha(path) == item['sha256'], 'Pinned artifact missing or changed: ' + str(path))
    return path


def preflight(config):
    frozen.preflight()
    require(config.get('schema') == 'historical_joint_six_date_v1', 'Wrong joint diagnostic schema')
    require(config.get('runner_sha256') == sha(__file__), 'Joint runner source changed')
    require(config.get('total_seconds') == 1950 and config.get('horizon') == 6 and config.get('maximum_maps') == 3,
            'Exact 1950-second, six-date, three-map budget required')
    endpoint = json.loads(pinned(config['endpoint_receipt']).read_text())
    require(endpoint.get('reference_manifest_sha256') == frozen.REFERENCE_SHA and endpoint.get('psi_child') == PSI,
            'Changed-psi endpoint belongs to another reference/experiment')
    require(endpoint.get('housing') == 'fixed_stock' and endpoint.get('repeat_verified') is True
            and endpoint.get('native_one_step_verified') is True and endpoint.get('terminal', {}).get('all_checks_pass') is True,
            'Original changed endpoint repeat/native gates required')
    require(endpoint['checkpoint']['sha256']=='3d4d08e323be40f6a6c19267aeea167e79ad1133f6e78a55ecec79cc6fc00f06',
            'Expected exact verified changed endpoint checkpoint')
    pinned(endpoint['checkpoint'])
    seed = json.loads(pinned(config['seed_receipt']).read_text())
    require(seed.get('horizon') == 12 and seed.get('mapping_count') == 5 and seed.get('psi_child') == frozen.PSI
            and seed.get('reference_manifest_sha256') == frozen.REFERENCE_SHA and seed.get('housing') == 'fixed_stock'
            and seed.get('closure') == 'fixed_tax' and seed.get('expectations') == 'current_shock_permanent_until_next_surprise',
            'A matching measured 12-date original-reference seed is required; no guessed default')
    pinned(seed['matrix'])
    compatibility=seed_compatibility(config,seed)
    return dict(status='PASS', model_calls=0, endpoint_reused=True, fresh_endpoint_calls=0,
        horizon=6, maximum_maps=3, maximum_policy_calls_conservative=36, total_seconds=1950,
        economic_change=dict(psi_child=dict(reference=frozen.PSI, experiment=PSI, classification='experimental')),
        all_other_primitives='authenticated historical block0506 unchanged', current_floor_validation=False,
        production_ready=False, native_setup_verified=False, horizon_comparison_verified=False,seed_compatibility=compatibility)


def seed_compatibility(config,seed):
    """Authenticate the reviewed three source changes without changing receipt bytes."""
    current=json.loads((frozen.HERE/'legacy_source_pins.json').read_text())
    current.pop('run_e5f_preference_budget_diagnostic.py')
    prior=seed['source_pins'];require(set(prior)==set(current),'Exact original nine-source seed inventory required')
    expected={'e5f_exact_policy_cache.py','e5f_preference_shock_fit.py','run_e5f_preference_estimation.py'}
    require({name for name in current if current[name]!=prior[name]}==expected,'Unexpected measured-seed source difference')
    sources=config['seed_source_archive'];require(set(sources)==expected,'Three exact seed-era source pins required')
    texts={}
    for name,item in sources.items():
        require(item['sha256']==prior[name],'Seed-era archived source identity differs: '+name)
        texts[name]=pinned(item).read_text()
    old=texts['e5f_exact_policy_cache.py']
    for context in ('key = exact_call_key(original, args, kwargs)', 'payload = pickle.dumps(result, protocol=5)'):
        before=context+'\n        except Exception:'
        after=context+'\n        except TimeoutError:\n            raise\n        except Exception:'
        require(old.count(before)==1,'Reviewed cache transform context changed')
        old=old.replace(before,after,1)
    require(old==(frozen.HERE/'pinned_tools/e5f_exact_policy_cache.py').read_text(),
            'Cache difference exceeds the two reviewed timeout-propagation insertions')
    def functions(text):
        tree=ast.parse(text);return {node.name:node for node in ast.walk(tree) if isinstance(node,ast.FunctionDef)}
    before=functions(texts['run_e5f_preference_estimation.py'])
    after=functions((frozen.HERE/'pinned_tools/run_e5f_preference_estimation.py').read_text())
    for name in ('draft_plan','_reconstruct_seed','initial_jacobian'):
        require(ast.dump(before[name],include_attributes=False)==ast.dump(after[name],include_attributes=False),
                'Used estimator computation differs: '+name)
    # The sole reviewed initializer addition stores the existing deadline.
    initializer=copy.deepcopy(after['__init__'])
    initializer.body=[node for node in initializer.body if not(isinstance(node,ast.Assign) and
        any(isinstance(target,ast.Attribute) and target.attr=='stage_deadline' for target in node.targets))]
    require(ast.dump(before['__init__'],include_attributes=False)==ast.dump(initializer,include_attributes=False),
            'Estimator initialization differs beyond stage deadline bookkeeping')
    return dict(status='PASS',exact_measured_mapping_modules=6,source_differences={name:dict(seed=prior[name],runtime=current[name]) for name in sorted(expected)},
        cache_change='two TimeoutError propagation insertions only; successful-call semantics unchanged',
        estimator_used_computations_unchanged=True,seed_validator='original archived nine-source validator',
        shock_fit_invoked=False,estimator_endpoint_path_fit_invoked=False,derivative_refresh=False,
        transport='lead-reviewed numerical warm start; receipt and matrix preserved byte-for-byte')


class MapBudget:
    def __init__(self, deadline):
        self.deadline = deadline
        self.started = 0

    def reserve(self, physical_pass):
        if self.started >= 3:
            raise TimeoutError('Three fresh mappings exhausted')
        if self.started == 2 and not physical_pass:
            raise TimeoutError('Fresh replay refused: no candidate clears both physical root gates')
        future_maps = 2-self.started
        remaining = self.deadline-time.monotonic()-RENDER_RESERVE-future_maps*CASE_SECONDS
        if remaining <= 0:
            raise TimeoutError('No mapping budget after reserving fresh replay and plotting')
        self.started += 1
        return min(CASE_SECONDS, remaining)


def guarded(seconds, call):
    def timeout(*_):
        raise TimeoutError('Joint diagnostic native-call deadline reached')
    previous = signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        return call()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def solve_paths(*, inner, acceleration, packet, evaluator, terminal, endpoint, plan, jacobian, output, deadline, native=True):
    """Use the original root and callbacks; all state/economics owned by caller."""
    import numpy as np
    p=plan['path']; P=packet['parameters']; q0=float(packet['solution'].p_eq[0]); qT=endpoint['price']
    latest={}; best_score=math.inf; best_physical=False; receipts=[]; budget=MapBudget(deadline)
    path=np.full(6, PSI)
    def evaluate(q,b):
        nonlocal latest, best_score, best_physical
        seconds=budget.reserve(best_physical)
        name=f'map_{budget.started:03d}'
        inner.write(output/'started_mapping.json',dict(name=name,prices=q,pensions=b,seconds=seconds,epoch=time.time()))
        if native:
            from e5f_social_security_root import CandidateDomainError
            try:evaluator.rt['primitive'].pf.rents_from_asset_prices(q,qT,P)
            except ValueError as exc:raise CandidateDomainError(str(exc)) from exc
        result,record=guarded(seconds,lambda:inner.mapping(packet,evaluator,terminal,endpoint,q,b,path,
            'fixed_stock',output/name,64*1024**3,capture=True,measure_fertility=True,start_year=2007))
        valid=all(record['gates'].values())
        score=max(max(map(abs,record['market_residual']))/2e-4,max(map(abs,record['fiscal_residual']))/1e-6)
        physical=valid and score<=1
        check=inner.terminal_checks(packet,evaluator,terminal,endpoint,result,path,p)
        receipt=dict(name=name,mapping_valid=valid,physical_root_gates_pass=physical,
            max_market_error=max(map(abs,record['market_residual'])),max_fiscal_error=max(map(abs,record['fiscal_residual'])),
            terminal=inner.plain(check),terminal_interpretation=('terminal check conditional on joint root certification' if not physical else 'candidate terminal check'),
            prices=q,pensions=b,record=record)
        receipts.append(receipt); latest=dict(result=result,record=record,terminal=check,prices=q.copy(),pensions=b.copy())
        inner.write(output/'latest_completed.json',receipt)
        inner.dump_checkpoint(output/(name+'_checkpoint.pkl.gz'),dict(terminal_state=result.terminal_state,
            prices=q,pensions=b,psi_path=path,rows=result.rows))
        if valid and score<best_score:
            best_score=score;best_physical=physical;inner.write(output/'best_so_far.json',receipt)
        return dict(mapping_valid=valid,market_residual=record['market_residual'],fiscal_residual=record['fiscal_residual'],payload=dict(evaluation=budget.started))
    def progress(row):
        inner.write(output/'root_progress.json',row)
    root=acceleration.solve_joint_with_acceleration(closure='fixed_tax',
        initial_prices=np.linspace(q0,qT,6),initial_fiscal_values=np.linspace(P.pension,terminal['parameters'].pension,6),
        evaluate=evaluate,project_prices=lambda q:np.clip(q,q0*p['price_bound_ratios'][0],q0*p['price_bound_ratios'][1]),
        fiscal_bounds=[P.pension*v for v in p['pension_bound_ratios']],market_tolerance=p['market_tolerance'],
        fiscal_tolerance=p['fiscal_tolerance'],market_slope=p['market_slope'],fiscal_slope=p['fiscal_slope'],
        max_log_step=p['max_log_step'],damping=p['damping'],max_evaluations=3,
        deadline_monotonic=deadline-RENDER_RESERVE,max_condition_number=plan['fit']['max_condition_number'],
        worsening_factor=plan['fit']['worsening_factor'],final_reproduction_tolerance=p['final_reproduction_tolerance'],
        callback=progress,initial_jacobian=jacobian)
    inner.write(output/'root.json',root)
    check=latest.get('terminal')
    summary=dict(root_certified=bool(root['converged']),terminal_pass=bool(check and check['all_checks_pass']),
        terminal_failure_interpretable=bool(root['converged']),maps_started=budget.started,maps_completed=len(receipts),
        initial_prices='original linear reference-to-endpoint guess; all six prices endogenous to root',
        no_endpoint_price_overwrite=True,production_ready=False,current_floor_validation=False,full_horizon_verified=False,
        historical_repayment_limitation=True,fitted_shocks=False)
    inner.write(output/'joint_receipt.json',summary)
    return summary,latest


def setup(config, output, modules):
    import numpy as np
    inner=modules['run_e5f_preference_transition.py'];estimator=modules['run_e5f_preference_estimation.py']
    manifest,packet,evaluator=inner.load_reference(output/'reference')
    endpoint=inner.read(pinned(config['endpoint_receipt']))
    with gzip.open(pinned(endpoint['checkpoint']),'rb') as stream:
        terminal=pickle.load(stream)
    inner.check_endpoint_primitives(packet['parameters'],terminal['parameters'])
    require(endpoint['source_manifest_sha256']==manifest['source_manifest']['sha256'], 'Endpoint engine source differs')
    require(np.array_equal(packet['b_grid'],terminal['b_grid']) and np.array_equal(terminal['stationary_g_pre'],terminal['evaluation'].g_pre),
            'Endpoint grid/native distribution differs')
    require(abs(float(np.sum(terminal['stationary_g_pre']))-1)<=1e-9, 'Endpoint normalized population differs')
    require(endpoint['psi_child']==terminal['parameters'].psi_child and endpoint['pension']==terminal['parameters'].pension
            and endpoint['price']==float(terminal['solution'].p_eq[0]), 'Endpoint checkpoint coordinates differ')
    plan=estimator.draft_plan('one_permanent');plan['budget']['total_seconds']=1950
    plan['path']['max_evaluations']=3;plan['path']['cache_max_bytes']=64*1024**3
    plan['source_pins']=json.loads((frozen.HERE/'legacy_source_pins.json').read_text())
    plan['source_pins'].pop('run_e5f_preference_budget_diagnostic.py')
    seed=inner.read(pinned(config['seed_receipt']))
    # Original validation checks actual reference psi, source, closure and proof.
    compatibility=seed_compatibility(config,seed);inner.write(output/'seed_compatibility.json',compatibility)
    validation_plan=copy.deepcopy(plan);validation_plan['source_pins']=seed['source_pins']
    native=estimator.NativeEstimator(validation_plan,output,manifest,packet,evaluator);native.seed_receipt=seed
    archive=pinned(config['seed_source_archive']['run_e5f_preference_estimation.py'])
    spec=importlib.util.spec_from_file_location('authenticated_seed_era_estimator',archive)
    original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
    original.NativeEstimator._validate_seed_receipt(native,seed)
    measured=np.load(pinned(seed['matrix']),allow_pickle=False)
    require(np.array_equal(measured,native._reconstruct_seed(seed,12)), 'Measured matrix differs from exact lag reconstruction')
    jacobian=native.initial_jacobian(6)
    row_scale=np.r_[np.ones(6),np.full(6,200.)]
    require(np.isfinite(jacobian).all() and np.linalg.cond(row_scale[:,None]*jacobian)<=plan['fit']['max_condition_number'],
            'Measured six-date seed is ill-conditioned; guessed-slope fallback forbidden')
    return inner,packet,evaluator,terminal,endpoint,plan,jacobian


def run(config,output):
    started=time.monotonic();deadline=started+1950
    proof=preflight(config)
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and os.sys.platform=='linux','Torch Slurm only')
    for key in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
        require(os.environ.get(key)=='1','One numerical thread required: '+key)
    require(not output.exists(),'Output already exists');output.mkdir(parents=True)
    stop=threading.Event();last={};inner=None
    def timeout(*_):raise TimeoutError('Absolute 1950-second joint diagnostic deadline reached')
    previous=signal.signal(signal.SIGUSR1,timeout)
    def watchdog():
        if not stop.wait(max(0,deadline-time.monotonic())):os.kill(os.getpid(),signal.SIGUSR1)
    def heartbeat():
        while not stop.wait(30):
            if inner:inner.write(output/'heartbeat.json',dict(epoch=time.time(),elapsed_seconds=time.monotonic()-started,phase='active'))
    timer=threading.Thread(target=watchdog,daemon=True);timer.start()
    worker=threading.Thread(target=heartbeat,daemon=True);worker.start()
    try:
        modules=frozen.import_pinned();inner=modules['run_e5f_preference_transition.py']
        inner.write(output/'preflight.json',proof);inner.write(output/'heartbeat.json',dict(epoch=time.time(),phase='setup'))
        inner.write(output/'latest_completed.json',dict(status='no_completed_mapping'))
        inner.write(output/'best_so_far.json',dict(status='no_completed_mapping'))
        inner,packet,evaluator,terminal,endpoint,plan,jacobian=guarded(min(300,deadline-time.monotonic()),lambda:setup(config,output,modules))
        import e5f_four_shock_acceleration as acceleration
        summary,last=solve_paths(inner=inner,acceleration=acceleration,packet=packet,evaluator=evaluator,
            terminal=terminal,endpoint=endpoint,plan=plan,jacobian=jacobian,output=output,deadline=deadline)
        visual=guarded(max(.001,deadline-time.monotonic()),lambda:inner.render_diagnostics(last['record']['diagnostic_packets'],
            output/'diagnostics',evaluator.rt['audit'],inner.read(inner.MANIFEST)['standard_diagnostic_names']))
        inner.write(output/'complete.json',dict(summary,elapsed_seconds=time.monotonic()-started,diagnostics=visual))
    except BaseException as exc:
        if inner and (output/'latest_completed.json').is_file() and time.monotonic()<deadline and not (output/'diagnostics').exists():
            try:
                saved=inner.read(output/'latest_completed.json')
                if 'record' in saved:
                    visual=guarded(max(.001,deadline-time.monotonic()),lambda:inner.render_diagnostics(saved['record']['diagnostic_packets'],
                        output/'diagnostics',evaluator.rt['audit'],inner.read(inner.MANIFEST)['standard_diagnostic_names']))
                    inner.write(output/'diagnostics_receipt.json',visual)
            except BaseException as render_exc:
                inner.write(output/'diagnostics_failure.json',dict(error=str(render_exc),captured_packets_retained=True))
        if inner:inner.write(output/'failure.json',dict(error_type=type(exc).__name__,error=str(exc),elapsed_seconds=time.monotonic()-started,
            root_certified=False,terminal_failure_interpretable=False,production_ready=False))
        raise
    finally:
        stop.set();worker.join(timeout=1);timer.join(timeout=1);signal.signal(signal.SIGUSR1,previous)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--config-sha256',required=True);parser.add_argument('--preflight',action='store_true');parser.add_argument('--output',type=Path)
    args=parser.parse_args();require(sha(args.config)==args.config_sha256,'Config SHA differs')
    config=json.loads(args.config.read_text())
    if args.preflight:print(json.dumps(preflight(config),indent=2,sort_keys=True))
    else:
        require(args.output is not None,'Output required');run(config,args.output)


if __name__=='__main__':main()
