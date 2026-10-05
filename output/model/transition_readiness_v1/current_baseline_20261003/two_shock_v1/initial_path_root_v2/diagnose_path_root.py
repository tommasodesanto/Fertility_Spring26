#!/usr/bin/env python3
"""Full original fixed-psi 24-date joint-root diagnostic using frozen v5 code.

One fresh native 12-date/five-map Jacobian and one original joint root only.
This does not fit fertility, compare a 32-date path, or certify production.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import threading
import time
import traceback
from types import SimpleNamespace

ACCEL_RELATIVE = 'code/model/experiments/transition_readiness/pinned_tools/e5f_four_shock_acceleration.py'
NATIVE_GATE_KEYS = {'mass','policy_reproduction','projection','dated_audits'}
ROOT_EVALUATIONS = 12  # Includes the original root's reserved fresh replay.
NATIVE_CALLS = 720
TOTAL_SECONDS = 7080.


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def pinned(item):
    require(isinstance(item, dict) and set(item) == {'path','sha256'}, 'Exact path/SHA pin required')
    path = Path(item['path']).resolve(strict=True)
    require(path.is_file() and sha(path) == item['sha256'], 'Pinned file differs: ' + str(path))
    return path


def write_output(path, value):
    """Atomic JSON with every numeric array element retained; fail on other types."""
    def plain(x):
        if isinstance(x,dict):
            require(all(isinstance(k,str) for k in x), 'JSON receipt keys must be strings')
            return {k:plain(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):
            return [plain(v) for v in x]
        if isinstance(x,(str,bool,int)) or x is None:
            return x
        if isinstance(x,float):
            require(math.isfinite(x), 'Nonfinite receipt scalar')
            return x
        if type(x).__module__.startswith('numpy'):
            import numpy as np
            if isinstance(x,np.ndarray):
                require(x.dtype.kind in 'biuf' and np.isfinite(x).all(),
                        'Only finite numeric arrays can enter receipts')
                return plain(x.tolist())
            if isinstance(x,np.generic):
                return plain(x.item())
        raise TypeError('Unsupported receipt type: '+type(x).__name__)
    path = Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(plain(value),sort_keys=True,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def load_contract(config):
    require(set(config) == {'v2_config','v2_helper','driver'}, 'Root diagnostic config keys differ')
    require(pinned(config['driver']) == Path(__file__).resolve(), 'Executing root driver differs from self pin')
    source_config = pinned(config['v2_config'])
    require(source_config.name == 'config.json' and source_config.parent.name == 'initial_path_diagnostic_v2',
            'Approved v2 config location differs')
    base = json.loads(source_config.read_text())
    helper_path = pinned(config['v2_helper'])
    require(base.get('driver') == config['v2_helper'] and helper_path.name == 'diagnose_initial_path.py',
            'Approved v2 helper pin differs')
    spec = importlib.util.spec_from_file_location('pinned_initial_path_v2', helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    require(Path(helper.__file__).resolve() == helper_path, 'Imported v2 helper source differs')
    root, source, paths = helper.validate_config(base)
    require(root.name == 'execution_smoke_v5' and source.is_relative_to(root), 'Frozen v5 source differs')
    return helper, base, root, source, paths


def root_controls(plan):
    controls = plan['path']
    gates = plan['gates']
    require(controls['max_evaluations'] == 12 and
            controls['price_bound_ratios'] == [.05,20.] and
            controls['pension_bound_ratios'] == [.05,20.] and
            controls['max_log_step'] == .15 and controls['damping'] == .7,
            'Original joint-root controls differ')
    require(gates['market_tolerance'] == 2e-4 and gates['fiscal_tolerance'] == 2e-5 and
            gates['final_reproduction_tolerance'] == 1e-10,
            'Original physical/replay gates differ')
    require(plan['seed'] == dict(horizon=12,perturbed_date=5,log_step=1e-5) and
            plan['budget']['seed_seconds'] == 1800 and plan['budget']['mapping_seconds'] == 1800,
            'Original seed/map controls differ')
    return dict(max_evaluations=ROOT_EVALUATIONS, market_tolerance=gates['market_tolerance'],
                fiscal_tolerance=gates['fiscal_tolerance'],
                final_reproduction_tolerance=gates['final_reproduction_tolerance'],
                price_bound_ratios=controls['price_bound_ratios'],
                pension_bound_ratios=controls['pension_bound_ratios'],
                max_log_step=controls['max_log_step'],damping=controls['damping'],
                max_condition_number=plan['fit']['max_condition_number'],
                worsening_factor=plan['fit']['worsening_factor'])


def require_root_callback_slot(completed):
    require(type(completed) is int and 0 <= completed < ROOT_EVALUATIONS,
            'Twelve-call diagnostic root callback cap reached')


def require_native_slot(current_calls, deadline):
    require(type(current_calls) is int and 0 <= current_calls < NATIVE_CALLS,
            'Independent 720-native-call cap reached')
    require(time.monotonic() < deadline, 'Independent diagnostic deadline reached')


def require_accelerator_source(module, root, plan):
    expected = root/'frozen/source'/ACCEL_RELATIVE
    pins = plan['identity']['source_pins']
    require(ACCEL_RELATIVE in pins and Path(module.__file__).resolve() == expected and
            sha(expected) == pins[ACCEL_RELATIVE],
            'Original frozen joint-root implementation differs')
    return dict(path=str(expected),sha256=pins[ACCEL_RELATIVE])


def validate_seed(seed, folder, runtime, d, initial_hash, start_calls):
    require(seed.get('accounting_valid') is True and seed.get('policy_calls') ==
            runtime.rt.total_native_calls - start_calls, 'Fresh seed call accounting differs')
    require(seed.get('mapping_count') == 5 or len(list(folder.glob('map_*/native_record.json'))) == 5,
            'Exactly five native seed maps required')
    require(seed.get('identity') == dict(runtime.rt.identity(),stage_start_year=2007,
                                          inherited_state_sha256=initial_hash),
            'Fresh seed identity differs')
    evidence = seed.get('source_evidence')
    require(isinstance(evidence, list) and len(evidence) == 5,
            'Five pinned native seed map records required')
    baseline = json.loads((folder/'baseline_checks.json').read_text())
    require(baseline.get('valid') is True and baseline.get('terminal',{}).get('all_checks_pass') is True,
            'Original baseline seed/terminal gate failed')
    for i, item in enumerate(evidence,1):
        require(d.pinned(item) == folder/f'map_{i:03d}'/'native_record.json',
                'Seed map evidence path/hash differs')
        record = json.loads(d.pinned(item).read_text())
        gates = record.get('gates')
        require(record.get('accounting_valid') is True and isinstance(gates,dict) and
                set(gates) == NATIVE_GATE_KEYS and all(value is True for value in gates.values()),
                'Fresh seed map native gates failed')
    require(seed.get('horizon') == 12 and seed.get('perturbed_date') == 5,
            'Fresh seed horizon/date differs')
    return evidence


def execute(config, output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); deadline = started + TOTAL_SECONDS
    phase = 'config'; runtime = None; completed_maps = 0; helper = None
    lock = threading.Lock(); stop = threading.Event()
    def progress(name, **fields):
        with lock:
            payload = dict(phase=name,epoch=time.time(),
                elapsed_seconds=time.monotonic()-started,
                actual_native_calls=getattr(getattr(runtime,'rt',None),'total_native_calls',0),
                completed_root_maps=completed_maps,**fields)
            write_output(output/'progress.json',payload)
    def pulse():
        while not stop.wait(180.):
            progress('native_in_progress' if phase in ('seed','root') else phase)
    thread = threading.Thread(target=pulse,name='joint-root-heartbeat',daemon=True)
    thread.start()
    try:
        helper, base, root, source, endpoint_paths = load_contract(config)
        phase = 'preflight'
        progress(phase, diagnostic_only=True, empirical_fitted=False,
                 scientific_validation=False, production_ready=False)
        sys.path.insert(0,str(source))
        import two_shock as d
        import two_shock_runtime as native_module
        require(Path(d.__file__).resolve() == source/'two_shock.py' and
                Path(native_module.__file__).resolve() == source/'two_shock_runtime.py',
                'Frozen driver/runtime import path differs')
        plan = helper.read_pin(base['fit_manifest'],root)
        controls = root_controls(plan)
        import numpy as np
        with native_module.retained.watchdog(deadline):
            preflight = d.preflight(plan)
            require(preflight['native_calls'] == 0, 'Preflight made native calls')
            write_output(output/'diagnostic_profile.json',dict(
                profile='full_original_fixed_psi24_root_diagnostic_v2',
                diagnostic_only=True,empirical_fitted=False,scientific_validation=False,
                production_ready=False,horizon=24,seed_horizon=12,
                root_max_evaluations_including_replay=ROOT_EVALUATIONS,
                maximum_actual_native_calls=NATIVE_CALLS,total_seconds=TOTAL_SECONDS,
                seed_seconds=plan['budget']['seed_seconds'],
                mapping_seconds=plan['budget']['mapping_seconds'],
                root_controls=controls,source_config=config['v2_config'],
                source_helper=config['v2_helper'],source_manifest=base['fit_manifest'],
                initial_paths='explicit v2 config paths; no prior root or Jacobian reuse'))
            phase = 'constructor'; progress(phase)
            runtime = native_module.NativeRuntime(plan,output/'runtime')
            rt = runtime.rt
            require(rt.total_native_calls == 0, 'Constructor made native calls')
            phase = 'restore_reference'; progress(phase)
            helper.bind_reference_cache_outputs(rt,base['reference_receipt'],root,output)
            helper.restore_reference_with_bridge_redirect(rt,base['reference_receipt'],root,
                                                          output/'reference_restored')
            require(rt.total_native_calls == 0, 'Authenticated reference restore made native calls')
            initial = copy.deepcopy(rt.initial_state)
            require(hashlib.sha256(np.asarray(initial.g_pre).tobytes()).hexdigest() == helper.INITIAL_SHA,
                    'Original initial population hash differs')
            initial_hash = d.state_hash(initial,rt.pf.birth_queue_values)
            require(d.state_hash(rt.initial_state,rt.pf.birth_queue_values) == initial_hash,
                    'Original queues/state differ before seed')
            phase = 'endpoint'; progress(phase, initial_state_sha256=initial_hash)
            saved, proofs = helper.validate_endpoint(base,endpoint_paths,rt.identity())
            require(native_module.retained.stationary_mapping_valid(proofs['stationary']) and
                    native_module.retained.stationary_mapping_valid(proofs['latest_completed']),
                    'Original stationary endpoint gates failed')
            require(np.array_equal(saved['b_grid'],rt.grid), 'Endpoint grid differs')
            helper.same_public_parameters(saved['parameters'],rt.packet['parameters'])
            terminal_v_sha = hashlib.sha256(np.asarray(saved['solution'].V).tobytes()).hexdigest()
            require(terminal_v_sha == helper.TERMINAL_V_SHA, 'Endpoint V bytes differ')
            terminal = dict(evaluation=SimpleNamespace(policy=SimpleNamespace(V=saved['solution'].V)))
            endpoint = dict(price=helper.QT,population_scale=proofs['latest_completed']['population_scale'])
            require(rt.total_native_calls == 0, 'Endpoint extraction made native calls')
            original_guard = rt._guard_native_call
            def counted_guard():
                require_native_slot(rt.total_native_calls,deadline)
                original_guard(); progress('native_call')
            rt._guard_native_call = counted_guard
            try:
                with rt.native_budget(deadline,NATIVE_CALLS):
                    phase = 'seed'; progress(phase)
                    seed_folder = output/'fresh_seed'
                    seed_start_calls = rt.total_native_calls
                    seed_deadline = min(deadline,time.monotonic()+plan['budget']['seed_seconds'])
                    with native_module.retained.watchdog(seed_deadline):
                        seed = runtime.measure_seed(stage=0,start_year=2007,
                            inherited_state=copy.deepcopy(initial),folder=seed_folder,
                            deadline=seed_deadline,**plan['seed'])
                    write_output(output/'fresh_seed_full_receipt.json',seed)
                    seed_evidence = validate_seed(seed,seed_folder,runtime,d,initial_hash,seed_start_calls)
                    write_output(output/'fresh_seed_receipt.json',dict(status='five_native_maps_completed',
                        seed_horizon=12,perturbed_date=5,policy_calls=seed['policy_calls'],
                        actual_native_calls=rt.total_native_calls,source_evidence=seed_evidence,
                        identity=seed['identity'],diagnostic_only=True,
                        original_measurement_receipt=d.pin(seed_folder/'measured/receipt.json'),
                        original_matrix_file=d.pin(seed_folder/'measured/jacobian.npy'),
                        full_augmented_receipt=d.pin(output/'fresh_seed_full_receipt.json')))
                    require(d.state_hash(rt.initial_state,rt.pf.birth_queue_values) == initial_hash,
                            'Original queues/state changed during seed')
                    phase = 'root'; progress(phase)
                    original, _ = native_module.retained.original_modules()
                    from e5f_four_shock_acceleration import extend_measured_jacobian, solve_joint_with_acceleration
                    accelerator = sys.modules['e5f_four_shock_acceleration']
                    accelerator_pin = require_accelerator_source(accelerator,root,plan)
                    J = extend_measured_jacobian(seed,24)
                    slopes = [float(np.median(np.abs(np.diag(J)[i*24:(i+1)*24]))) for i in range(2)]
                    require(all(math.isfinite(x) and x > 0 for x in slopes),
                            'Fresh measured own slopes invalid')
                    q0 = rt.reference_price; b0 = float(rt.P.pension)
                    q = np.asarray(base['q_path'],float); b = np.asarray(base['pension_path'],float)
                    psi = np.asarray(base['psi_path'],float)
                    root_folder = output/'joint_root'; root_folder.mkdir()
                    latest = {}; best = {}; root_start_calls = rt.total_native_calls
                    root_policy_calls = 0
                    def evaluate(q_candidate,b_candidate):
                        nonlocal completed_maps, root_policy_calls
                        require_root_callback_slot(completed_maps)
                        require(time.monotonic() < deadline, 'Shared diagnostic deadline exhausted')
                        number = completed_maps + 1
                        map_deadline = min(deadline,time.monotonic()+plan['budget']['mapping_seconds'])
                        folder = root_folder/f'map_{number:03d}'
                        progress('root_mapping', map_number=number)
                        with native_module.retained.watchdog(map_deadline),rt.native_budget(map_deadline,NATIVE_CALLS-rt.total_native_calls):
                            native, record = rt.mapping(terminal,endpoint,q_candidate,b_candidate,psi,
                                folder,initial_state=copy.deepcopy(initial),start_year=2007)
                        completed_maps = number
                        require(record['policy_calls'] == rt.total_native_calls-root_start_calls-root_policy_calls,
                                'Root callback actual native-call count differs')
                        root_policy_calls += record['policy_calls']
                        require(record['accounting_valid'] is True and all(record['gates'].values()) and
                                len(record['rows']) == 24 and len(native.dated_states) == 24,
                                'Original dated map gates or length failed')
                        require(all(math.isfinite(float(x)) for key in ('market_residual','fiscal_residual')
                                    for x in record[key]), 'Nonfinite physical residual')
                        d.write(folder/'native_record.json',record)
                        pin = d.pin(folder/'native_record.json')
                        latest.update(native=native,record=record,pin=pin)
                        score = max(max(map(abs,record['market_residual']))/controls['market_tolerance'],
                                    max(map(abs,record['fiscal_residual']))/controls['fiscal_tolerance'])
                        case = dict(map_number=number,score=score,market_maximum=max(map(abs,record['market_residual'])),
                                    fiscal_maximum=max(map(abs,record['fiscal_residual'])),
                                    actual_native_calls=rt.total_native_calls,record=pin,
                                    prices=np.asarray(q_candidate).tolist(),pensions=np.asarray(b_candidate).tolist())
                        write_output(root_folder/'latest_completed.json',case)
                        if not best or score < best['score']:
                            best.update(case)
                            write_output(root_folder/'best_so_far.json',best)
                        progress('root_map_completed', map_number=number, residual_score=score)
                        return dict(mapping_valid=True,market_residual=record['market_residual'],
                                    fiscal_residual=record['fiscal_residual'])
                    def callback(row):
                        write_output(root_folder/'root_progress.json',row)
                        progress('root_iteration')
                    root_result = solve_joint_with_acceleration(
                        closure='fixed_tax',initial_prices=q,initial_fiscal_values=b,evaluate=evaluate,
                        project_prices=lambda values:np.clip(values,q0*controls['price_bound_ratios'][0],
                                                              q0*controls['price_bound_ratios'][1]),
                        fiscal_bounds=[b0*x for x in controls['pension_bound_ratios']],
                        market_tolerance=controls['market_tolerance'],fiscal_tolerance=controls['fiscal_tolerance'],
                        market_slope=slopes[0],fiscal_slope=slopes[1],
                        max_log_step=controls['max_log_step'],damping=controls['damping'],
                        max_evaluations=controls['max_evaluations'],deadline_monotonic=deadline,
                        max_condition_number=controls['max_condition_number'],
                        worsening_factor=controls['worsening_factor'],
                        final_reproduction_tolerance=controls['final_reproduction_tolerance'],
                        initial_jacobian=J,callback=callback)
                    write_output(root_folder/'root.json',root_result)
                    require(completed_maps <= ROOT_EVALUATIONS and rt.total_native_calls <= NATIVE_CALLS,
                            'Root map/native call allowance exceeded')
                    require(d.state_hash(rt.initial_state,rt.pf.birth_queue_values) == initial_hash and
                            d.state_hash(initial,rt.pf.birth_queue_values) == initial_hash,
                            'Original state or both queues changed')
                    require(root_start_calls + root_policy_calls == rt.total_native_calls and
                            seed_start_calls + seed['policy_calls'] == root_start_calls,
                            'Seed/root actual native-call ledger inconsistent')
                    gates = root_result.get('gates',{})
                    passed = bool(root_result.get('converged') is True and gates and
                                  all(gates.values()) and root_result.get('final_reproduction_max_abs') is not None and
                                  root_result['final_reproduction_max_abs'] <= controls['final_reproduction_tolerance'] and
                                  latest['record']['accounting_valid'] is True and all(latest['record']['gates'].values()))
                    status = 'root_converged' if passed else 'budget_limited_root_diagnostic_completed'
                    result = dict(status=status,root_converged=passed,
                                  profile='full_original_fixed_psi24_root_diagnostic_v2',diagnostic_only=True,
                                  empirical_fitted=False,scientific_validation=False,production_ready=False,
                                  horizon=24,seed_map_count=5,root_callback_count=completed_maps,
                                  seed_native_calls=seed['policy_calls'],root_native_calls=rt.total_native_calls-root_start_calls,
                                  actual_native_calls=rt.total_native_calls,initial_state_sha256=initial_hash,
                                  terminal_V_sha256=terminal_v_sha,root_controls=controls,
                                  root=root_result,latest_completed=latest['pin'],best=best,
                                  full_rows=latest['record']['rows'],
                                  market_residual_trajectory=latest['record']['market_residual'],
                                  fiscal_residual_trajectory=latest['record']['fiscal_residual'],
                                  source_config=config['v2_config'],source_helper=config['v2_helper'],
                                  source_manifest=base['fit_manifest'],endpoint_pins=base['endpoint_pins'],
                                  accelerator_source=accelerator_pin,
                                  fresh_seed_receipt=d.pin(output/'fresh_seed_full_receipt.json'))
                    write_output(output/'result.json',result)
                    phase='complete'; progress(phase)
                    return result
            finally:
                rt._guard_native_call = original_guard
    except BaseException as exc:
        failure = dict(status='failed',phase=phase,error=repr(exc),
            traceback=traceback.format_exc(),actual_native_calls=getattr(getattr(runtime,'rt',None),'total_native_calls',0),
            completed_root_maps=completed_maps,scientific_validation=False,production_ready=False)
        write_output(output/'failure.json',failure)
        progress('failed')
        raise
    finally:
        stop.set(); thread.join(timeout=1.)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args(argv)
    execute(json.loads(args.config.read_text()),args.output)


if __name__=='__main__':
    main()
