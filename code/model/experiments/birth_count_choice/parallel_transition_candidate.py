#!/usr/bin/env python3
"""One fixed-psi, one-horizon, fresh-native stage-one transition candidate.

The copied StageAdapter.evaluate below differs from frozen v5 only in its cold
price initializer and lossless extra receipts. Original roots, updates, bounds,
physical gates, endpoint, seed, and native mapping remain frozen dependencies.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import hashlib
import importlib
import json
import math
from pathlib import Path
import pickle
import sys
import threading
import time
import traceback
import numpy as np

RENT = .13689881249028354
MAX_EXTERNAL_SECONDS = 10800
MAX_INTERNAL_SECONDS = 10700
MAX_NATIVE_CALLS = 20000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def pinned(item):
    require(isinstance(item, dict) and set(item) == {'path', 'sha256'}, 'Exact path/SHA pin required')
    path = Path(item['path']).resolve(strict=True)
    require(path.is_file() and sha(path) == item['sha256'], 'Pinned file differs: '+str(path))
    return path


def pin(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), sha256=sha(path))


def write_full(path, value):
    """Atomic JSON retaining complete numeric arrays; nonfinite scalars are explicit strings."""
    def plain(x):
        if isinstance(x, dict):
            require(all(isinstance(k, str) for k in x), 'Receipt keys must be strings')
            return {k: plain(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [plain(v) for v in x]
        if isinstance(x, np.ndarray):
            require(x.dtype.kind in 'biuf', 'Only numeric arrays allowed in JSON')
            return plain(x.tolist())
        if isinstance(x, np.generic):
            return plain(x.item())
        if isinstance(x, float) and not math.isfinite(x):
            return str(x)
        if isinstance(x, (str, bool, int, float)) or x is None:
            return x
        raise TypeError('Unsupported receipt type: '+type(x).__name__)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(plain(value), sort_keys=True, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def no_arbitrage_prices(parameters, endpoint_price, horizon):
    require(type(horizon) is int and horizon in (24, 32), 'Candidate horizon must be 24 or 32')
    A = float(parameters.R_gross) + float(parameters.delta) + float(parameters.tau_H)
    qT = float(endpoint_price)
    require(math.isfinite(A) and A > 1 and math.isfinite(qT) and qT > 0, 'Invalid native no-arbitrage inputs')
    steady = RENT/(A-1)
    t = np.arange(horizon, dtype=float)
    q = steady + (qT-steady)*np.power(A, -(horizon-t))
    require(np.isfinite(q).all() and (q > 0).all(), 'Invalid dated price initialization')
    return q


def load_frozen(config):
    require(set(config) == {'driver','fit_manifest','frozen_package','psi','horizon','budget','task_id'}, 'Exact candidate config required')
    require(pinned(config['driver']) == Path(__file__).resolve(), 'Executing candidate driver differs from pin')
    root = Path(config['frozen_package']).resolve(strict=True)
    source = root/'code/model/experiments/birth_count_choice'
    require(root.name == 'source' and root.parent.name == 'frozen' and
            source.is_dir(), 'Frozen v5 package root required')
    manifest = pinned(config['fit_manifest'])
    require(manifest.name == 'fit_manifest.json' and manifest.parent.name == 'inputs' and
            manifest.parent.parent.name == 'execution_smoke_v5',
            'V5 fit manifest required')
    plan = json.loads(manifest.read_text())
    require(Path(plan['source_pins']['two_shock_driver']['path']).resolve() == source/'two_shock.py' and
            Path(plan['source_pins']['two_shock_runtime']['path']).resolve() == source/'two_shock_runtime.py',
            'Fit manifest frozen source path differs')
    require(sha(source/'two_shock.py') == plan['source_pins']['two_shock_driver']['sha256'] and
            sha(source/'two_shock_runtime.py') == plan['source_pins']['two_shock_runtime']['sha256'],
            'Frozen driver or runtime changed')
    require(type(config['psi']) in (int,float) and math.isfinite(config['psi']) and
            plan['stage_starts'][0]['bounds'][0] < config['psi'] < plan['stage_starts'][0]['bounds'][1],
            'Fixed candidate psi outside original absolute bounds')
    require(type(config['horizon']) is int and config['horizon'] in (24,32), 'Horizon must be 24 or 32')
    require(type(config['task_id']) is str and config['task_id'].startswith('c') and
            config['task_id'].endswith(f"_h{config['horizon']}") and len(config['task_id']) == 7,
            'Candidate task identity/horizon differs')
    require(config['budget'] == dict(external_seconds=MAX_EXTERNAL_SECONDS,internal_seconds=MAX_INTERNAL_SECONDS,
                                     maximum_policy_calls=MAX_NATIVE_CALLS), 'Candidate budget differs')
    require(plan['mode'] == 'diagnostic' and plan['smoke'] is False and
            plan['seed'] == dict(horizon=12, perturbed_date=5, log_step=1e-5), 'Original empirical mode/seed required')
    require(plan['path']['max_evaluations'] == 12 and plan['budget']['candidate_seconds'] == 7200 and
            plan['budget']['path_seconds'] == 6000 and plan['budget']['endpoint_seconds'] == 1800 and
            plan['budget']['mapping_seconds'] == 1800 and plan['budget']['seed_seconds'] == 1800,
            'Original candidate/root/endpoint/seed/map controls differ')
    sys.path.insert(0,str(source))
    driver = importlib.import_module('two_shock')
    native = importlib.import_module('two_shock_runtime')
    require(Path(driver.__file__).resolve() == source/'two_shock.py' and
            Path(native.__file__).resolve() == source/'two_shock_runtime.py', 'Frozen import differs')
    require(driver.preflight(plan)['native_calls'] == 0, 'Frozen preflight made native calls')
    # Preflight authenticates unmodified manifest. These limits only shorten its
    # available wall time; original per-stage controls and gates remain exact.
    plan = copy.deepcopy(plan)
    plan['budget']['total_seconds'] = MAX_INTERNAL_SECONDS
    plan['budget']['maximum_policy_calls'] = MAX_NATIVE_CALLS
    return plan, driver, native, manifest, source


def candidate_adapter_class(native):
    # Dynamic subclass keeps every inherited native operation in the frozen module.
    retained, write, pin = native.retained, native.write, native.pin
    class CandidateAdapter(native.StageAdapter):
        def evaluate(self, *, psi, start_year, horizon, seed, gates, budget, endpoint_controls, path_controls, deadline, folder):
            import numpy as np
            original, _ = retained.original_modules()
            from e5f_four_shock_acceleration import extend_measured_jacobian, solve_joint_with_acceleration
            folder = Path(folder)
            start = self.calls
            terminal, endpoint = self._endpoint(psi, folder/'endpoint', min(deadline,time.monotonic()+budget['endpoint_seconds']))
            p = path_controls
            psi_path = np.full(horizon,psi)
            qT,bT = endpoint['price'],terminal['parameters'].pension
            J = extend_measured_jacobian(seed,horizon)
            # Native own-lag slopes are used for the root's emergency reset too.
            slopes = [float(np.median(np.abs(np.diag(J)[i*horizon:(i+1)*horizon]))) for i in range(2)]
            require(all(math.isfinite(s) and s > 0 for s in slopes), 'Native measured own slopes absent; no default fallback')
            warm,warm_kind = self.select_warm(horizon,psi)
            initialization_kind = 'no_arbitrage_backward_from_fresh_endpoint' if warm is None else warm_kind
            write(folder/'warm_start.json',dict(kind=initialization_kind,target_psi_hex=float(psi).hex(),horizon=horizon,
                source_psi_hex=None if warm is None else warm['psi_hex'],identity=self.identity(),
                source_folder=None if warm is None else warm['source_folder'],
                fresh_native_mapping_required=True,residuals_or_fertility_reused=False))
            q = no_arbitrage_prices(self.rt.P, qT, horizon) if warm is None else warm['prices'].copy()
            b = np.full(horizon,bT) if warm is None else warm['fiscal_values'].copy()
            latest = {}
            best = {}
            attempted_count = 0
            completed_count = 0
            path_deadline = min(deadline,time.monotonic()+budget['path_seconds'])
            def evaluate(q,b):
                nonlocal attempted_count, completed_count
                attempted_count += 1
                map_number = attempted_count
                map_folder = folder/f'map_{map_number:03d}'
                native,record = self._mapping(terminal,endpoint,q,b,psi_path,map_folder,path_deadline)
                completed_count += 1
                latest.update(native=native,record=record,map_number=map_number,
                              mapping_pin=pin(map_folder/'native_record.json'))
                write_full(folder/'latest_completed_full.json',record)
                score=max(max(map(abs,record['market_residual']))/gates['market_tolerance'],
                          max(map(abs,record['fiscal_residual']))/gates['fiscal_tolerance'])
                if not best or score<best['score']:
                    best.update(score=score,map_number=map_number,record=latest['mapping_pin'])
                    write_full(folder/'best_so_far_full.json',best)
                return dict(mapping_valid=all(record['gates'].values()),market_residual=record['market_residual'],fiscal_residual=record['fiscal_residual'])
            root = solve_joint_with_acceleration(closure='fixed_tax',initial_prices=q,initial_fiscal_values=b,evaluate=evaluate,
                project_prices=lambda values:np.clip(values,self.q*p['price_bound_ratios'][0],self.q*p['price_bound_ratios'][1]),
                fiscal_bounds=[self.pension*x for x in p['pension_bound_ratios']],market_tolerance=gates['market_tolerance'],
                fiscal_tolerance=gates['fiscal_tolerance'],market_slope=slopes[0],fiscal_slope=slopes[1],
                max_log_step=p['max_log_step'],damping=p['damping'],max_evaluations=p['max_evaluations'],
                deadline_monotonic=path_deadline,max_condition_number=self.plan['fit']['max_condition_number'],
                worsening_factor=self.plan['fit']['worsening_factor'],final_reproduction_tolerance=1e-10,
                initial_jacobian=J if warm is None else warm['final_jacobian'].copy(),callback=lambda row:write_full(folder/'root_progress_full.json',row))
            original.inner.write(folder/'root.json',root)
            write_full(folder/'root_full.json',root)
            require(bool(latest), 'Original root produced no completed native mapping')
            record = latest['record']
            terminal_check = self.rt.terminal_checks(terminal,endpoint,latest['native'],psi_path,
                tolerance=1e-3,raw_queue_tolerance=1e-3)
            valid = bool(root['converged'] and terminal_check['all_checks_pass'])
            if valid or (self.plan['mode']=='diagnostic' and root['converged']):
                self.retain_warm(horizon,psi,root,folder)
            return dict(identity=self.identity(),reference_manifest_sha256=self.rt.identity()['reference_sha256'],
                source_pins=self.rt.identity()['source_pins'],housing=self.rt.housing,
                shock_contract=dict(start_year=self.start_year,psi=psi,expectations='permanent_until_next_surprise'),
                psi=psi,horizon=horizon,accounting_valid=bool(record['accounting_valid'] and all(record['gates'].values())),policy_calls=self.calls-start,
                root_and_terminal_pass=valid,stationary_pass=endpoint['stationary_pass'],
                root_pass=bool(root['converged']),terminal_pass=bool(terminal_check['all_checks_pass']),
                replay_pass=bool(root['gates']['market_replay'] and root['gates']['fiscal_replay']) if 'gates' in root else bool(root['converged']),
                stationary_renewal_gap=endpoint['stationary_renewal_gap'],
                market_maximum_residual=max(map(abs,record['market_residual'])),
                fiscal_maximum_residual=max(map(abs,record['fiscal_residual'])),
                replay_maximum_gap=(root['final_reproduction_max_abs']
                                    if root['final_reproduction_max_abs'] is not None else float('inf')),
                terminal=terminal_check,rows=record['rows'],fertility=record['fertility'],
                final_mapping_pin=latest['mapping_pin'],path_evaluations=attempted_count,
                path_mappings_completed=completed_count,latest_completed_map_number=latest['map_number'],
                native_reply=latest['native'],terminal_packet=terminal,endpoint=endpoint,psi_path=psi_path,root=root)

    return CandidateAdapter


def export_boundary(reply, runtime, driver, folder, *, accepted):
    """Save the exact own-vintage 2015 state, both queues, V and forecasts."""
    folder = Path(folder);folder.mkdir(parents=True,exist_ok=True)
    state = reply['dated_states'][2]['state']
    q = np.asarray(reply['prices'],float)
    b = np.asarray(reply['pensions'],float)
    V = np.asarray(reply['values'][2])
    require(len(q) == len(b) == reply['horizon'] and V.size > 0, 'Full dated forecast/boundary required')
    packet = dict(year=2015,start_year=2007,local_index=2,households=copy.deepcopy(state),
                  original_initial_state=copy.deepcopy(runtime.rt.initial_state),
                  boundary_price=float(q[2]),boundary_pension=float(b[2]),boundary_value=V.copy(),
                  first_two_values=[np.asarray(reply['values'][i]).copy() for i in (0,1)],
                  terminal_value=np.asarray(reply['values'][-1]).copy(),
                  forecast_prices=q.copy(),forecast_pensions=b.copy(),
                  forecast_psi=np.asarray(reply['psi_path'])[2:].copy(),
                  full_psi_path=np.asarray(reply['psi_path']).copy(),
                  reference_identity=runtime.rt.identity(),root_gate_passed=bool(accepted),
                  selected_stage_handoff=False)
    endpoint_path=folder/'fresh_endpoint_packet.pkl.gz'
    endpoint_tmp=endpoint_path.with_suffix(endpoint_path.suffix+'.tmp')
    with gzip.open(endpoint_tmp,'wb') as stream:
        pickle.dump(dict(terminal_packet=reply['terminal_packet'],endpoint=reply['endpoint'],
                         reference_identity=runtime.rt.identity()),stream,protocol=pickle.HIGHEST_PROTOCOL)
    endpoint_tmp.replace(endpoint_path)
    path = folder/'own_vintage_2015.pkl.gz'
    temporary = path.with_suffix(path.suffix+'.tmp')
    with gzip.open(temporary,'wb') as stream:
        pickle.dump(packet,stream,protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(path)
    with gzip.open(path,'rb') as stream:
        restored = pickle.load(stream)
    require(driver.state_hash(restored['households'],runtime.queue_values) ==
            driver.state_hash(state,runtime.queue_values) and
            np.array_equal(restored['boundary_value'],V) and
            np.array_equal(restored['forecast_prices'],q) and
            np.array_equal(restored['forecast_pensions'],b) and
            all(np.array_equal(restored['first_two_values'][i],reply['values'][i]) for i in (0,1)) and
            driver.state_hash(restored['original_initial_state'],runtime.queue_values)==
                driver.state_hash(runtime.rt.initial_state,runtime.queue_values),
            'Boundary checkpoint roundtrip differs')
    receipt = dict(year=2015,start_year=2007,local_index=2,source=pin(path),
                   fresh_endpoint_packet=pin(endpoint_path),
                   state_sha256=driver.state_hash(state,runtime.queue_values),
                   boundary_q=float(q[2]),boundary_b=float(b[2]),
                   boundary_V_sha256=hashlib.sha256(V.tobytes()).hexdigest(),
                   terminal_V_sha256=hashlib.sha256(np.asarray(reply['values'][-1]).tobytes()).hexdigest(),
                   population_sha256=hashlib.sha256(np.asarray(state.g_pre).tobytes()).hexdigest(),
                   initial_state_sha256=driver.state_hash(runtime.rt.initial_state,runtime.queue_values),
                   reference_identity=runtime.rt.identity(),
                   adjusted_queue=np.asarray(runtime.queue_values(state.scheduled_entries)).tolist(),
                   raw_queue=np.asarray(runtime.queue_values(state.scheduled_raw_entries)).tolist(),
                   forecast_prices=q.tolist(),forecast_pensions=b.tolist(),
                   forecast_psi=np.asarray(reply['psi_path'])[2:].tolist(),
                   root_gate_passed=bool(accepted),selected_stage_handoff=False,
                   pending_horizon_comparison=True,no_rescaling=True,nonanticipating_vintage=True)
    write_full(folder/'boundary_receipt.json',receipt)
    return receipt


def run(config, output):
    output = Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    started = time.monotonic();deadline = started + MAX_INTERNAL_SECONDS
    phase='config';runtime=None;stop=threading.Event();lock=threading.Lock()
    def progress(name,**fields):
        with lock:
            payload=dict(phase=name,epoch=time.time(),elapsed_seconds=time.monotonic()-started,
                         actual_native_calls=getattr(getattr(runtime,'rt',None),'total_native_calls',0))
            payload.update(fields)
            write_full(output/'progress.json',payload)
    def pulse():
        while not stop.wait(150):
            progress('native_in_progress' if phase in ('reference','seed','candidate') else phase)
    thread=threading.Thread(target=pulse,name='candidate-heartbeat',daemon=True)
    thread.start()
    try:
        plan,driver,native,manifest,source=load_frozen(config)
        write_full(output/'source_identity.json',dict(driver=config['driver'],fit_manifest=config['fit_manifest'],
                    frozen_driver=pin(source/'two_shock.py'),frozen_runtime=pin(source/'two_shock_runtime.py'),
                    reference_identity=plan['identity'],task_id=config['task_id'],psi=config['psi'],horizon=config['horizon'],
                    budget=config['budget'],scientific_validation=False))
        phase='constructor';progress(phase)
        native.StageAdapter=candidate_adapter_class(native)
        runtime=native.NativeRuntime(plan,output/'runtime')
        require(runtime.rt.total_native_calls == 0, 'Constructor made native calls')
        initial_calls=runtime.rt.total_native_calls
        runtime.bind_budget(deadline,lambda name,**fields: progress(name,**fields))
        with runtime.budget_context():
            phase='reference';progress(phase)
            prepared=runtime.prepare_reference(deadline,output/'fresh_reference')
            require(bool(prepared.get('reference_checkpoint')) and bool(prepared.get('reconstruction_receipt')),
                    'Fresh native reference evidence missing')
            initial=runtime.initial_state()
            initial_hash=driver.state_hash(initial,runtime.queue_values)
            write_full(output/'reference_receipt.json',prepared)
            phase='seed';progress(phase)
            seed_deadline=min(deadline,time.monotonic()+plan['budget']['seed_seconds'])
            seed=runtime.measure_seed(stage=0,start_year=2007,inherited_state=copy.deepcopy(initial),
                          folder=output/'fresh_seed',deadline=seed_deadline)
            require(seed['mapping_count'] == 5 and seed['horizon'] == 12 and
                    seed['identity'] == dict(plan['identity'],stage_start_year=2007,inherited_state_sha256=initial_hash) and
                    len(seed['source_evidence']) == 5 and seed['accounting_valid'] is True,
                    'Fresh original seed identity/count/accounting differs')
            for item in seed['source_evidence']:
                record=json.loads(driver.pinned(item).read_text())
                require(record['accounting_valid'] is True and all(record['gates'].values()), 'Native seed map gate failed')
            write_full(output/'seed_receipt.json',seed)
            require(driver.state_hash(runtime.rt.initial_state,runtime.queue_values)==initial_hash,
                    'Reference state changed during seed')
            phase='candidate';progress(phase,psi=config['psi'],horizon=config['horizon'])
            candidate_deadline=min(deadline,time.monotonic()+plan['budget']['candidate_seconds'])
            reply=runtime.evaluate_stage(stage=0,psi=float(config['psi']),horizon=config['horizon'],
                          start_year=2007,inherited_state=initial,deadline=candidate_deadline,
                          folder=output/'candidate')
            actual_candidate_calls=(runtime.rt.total_native_calls-initial_calls-seed['policy_calls']-
                                    prepared['policy_calls'])
            completed_candidate_calls=reply['policy_calls']
            unfinished_candidate_calls=actual_candidate_calls-completed_candidate_calls
            budget_interrupted=(reply['root'].get('status')=='time_or_evaluation_budget' and
                                reply['root'].get('converged') is False)
            require(unfinished_candidate_calls >= 0 and
                    (unfinished_candidate_calls == 0 or budget_interrupted),
                    'Native candidate-call ledger differs outside an unfinished budget outcome')
            require(driver.state_hash(runtime.rt.initial_state,runtime.queue_values)==initial_hash and
                    driver.state_hash(initial,runtime.queue_values)==initial_hash,
                    'Original initial state or queues changed during candidate')
            require(reply['identity'] == dict(plan['identity'],stage_start_year=2007,inherited_state_sha256=initial_hash) and
                    reply['source_pins'] == plan['identity']['source_pins'] and
                    reply['shock_contract'] == dict(start_year=2007,psi=float(config['psi']),
                                                    expectations='permanent_until_next_surprise') and
                    len(reply['rows']) == len(reply['dated_states']) == config['horizon'] and
                    all(r['calendar_year']==2007+4*i for i,r in enumerate(reply['rows'])),
                    'Native candidate source/clock/state differs')
            for evidence in reply['source_evidence']:
                driver.pinned(evidence)
            gates=plan['gates']
            accepted=bool(all(reply.get(key) is True for key in ('root_pass','replay_pass','stationary_pass','accounting_valid')) and
                reply['market_maximum_residual']<=gates['market_tolerance'] and
                reply['fiscal_maximum_residual']<=gates['fiscal_tolerance'] and
                reply['replay_maximum_gap']<=gates['final_reproduction_tolerance'] and
                reply['stationary_renewal_gap']<=gates['stationary_renewal_tolerance'])
            boundary=export_boundary(reply,runtime,driver,output/'own_vintage_2015',accepted=accepted)
            # Index four is the actual 2023 state for a 2007-vintage stage-one path.
            export=runtime.export_state(reply,index=4,year=2023,folder=output/'state_2023_checkpoint')
            require(export['exact_native_state'] is True and export['reconstructed_or_rescaled'] is False and
                    export['calendar_year']==2023 and export['period_index']==4,
                    'Exact own-vintage 2023 state export failed')
            fertility=reply['fertility']
            require(len(fertility)>=2,'First two historical fertility rows required')
            diagnostics=None
            if accepted:
                render_deadline=min(deadline,time.monotonic()+plan['budget']['render_seconds'])
                diagnostics=runtime.render_standard(reply,output/'standard_diagnostics',render_deadline)
                names=set(plan['standard_plot_names'])
                require(isinstance(diagnostics,dict) and diagnostics.get('sampled_dates'),
                        'Accepted candidate requires original standard diagnostic packet')
                for dated in diagnostics['sampled_dates']:
                    require(set(dated['plots'])==names,'Exact 17 standard plot names required')
                    directory=output/'standard_diagnostics'/f"date_{dated['period']:03d}"/'standard_diagnostics'
                    for name,value in dated['plots'].items():
                        require((directory/name).is_file() and sha(directory/name)==value,
                                'Standard diagnostic plot hash differs')
                write_full(output/'diagnostics.json',diagnostics)
            result=dict(status='root_gate_passed_pending_horizon' if accepted else 'candidate_unaccepted',
                        accepted_candidate=False,root_gate_passed=accepted,
                        pending_horizon_comparison=True,empirical_fitted=False,scalar_fitter_called=False,
                        scientific_validation=False,production_ready=False,two_shock_result=False,
                        task_id=config['task_id'],psi=float(config['psi']),horizon=config['horizon'],
                        reference_identity=plan['identity'],fit_manifest=config['fit_manifest'],
                        native_calls=runtime.rt.total_native_calls,seed_native_calls=seed['policy_calls'],
                        candidate_native_calls=actual_candidate_calls,
                        candidate_completed_operation_calls=completed_candidate_calls,
                        candidate_unfinished_native_calls=unfinished_candidate_calls,
                        root_pass=reply['root_pass'],replay_pass=reply['replay_pass'],
                        stationary_pass=reply['stationary_pass'],accounting_valid=reply['accounting_valid'],
                        terminal_diagnostic_pass=reply['terminal_pass'],terminal_diagnostic_gating=False,
                        market_maximum_residual=reply['market_maximum_residual'],
                        fiscal_maximum_residual=reply['fiscal_maximum_residual'],
                        replay_maximum_gap=(None if reply['root'].get('final_reproduction_max_abs') is None
                                            else reply['replay_maximum_gap']),
                        stationary_renewal_gap=reply['stationary_renewal_gap'],
                        root=reply['root'],full_dated_rows=reply['rows'],
                        final_mapping_pin=reply['final_mapping_pin'],
                        path_evaluations=reply['path_evaluations'],
                        path_mappings_completed=reply['path_mappings_completed'],
                        latest_completed_map_number=reply['latest_completed_map_number'],
                        market_residual_trajectory=reply['root']['final']['market_residual'] if reply['root'].get('final') else None,
                        fiscal_residual_trajectory=reply['root']['final']['fiscal_residual'] if reply['root'].get('final') else None,
                        fertility_rows=fertility,first_two_historical_fertility_rows=fertility[:2],
                        own_vintage_2015=boundary,stage_one_2023_diagnostic_continuation=export,
                        standard_diagnostics=diagnostics,
                        source_evidence=reply['source_evidence'])
            write_full(output/'result.json',result)
            write_full(output/'latest_completed.json',result)
            if accepted:
                write_full(output/'best_so_far.json',result)
            phase='complete';progress(phase,accepted_candidate=False,root_gate_passed=accepted,
                                      pending_horizon_comparison=True)
            return result
    except BaseException as exc:
        failure=dict(status='failed',phase=phase,error=repr(exc),traceback=traceback.format_exc(),
                     actual_native_calls=getattr(getattr(runtime,'rt',None),'total_native_calls',0),
                     accepted_candidate=False,scientific_validation=False,production_ready=False)
        write_full(output/'failure.json',failure)
        progress('failed')
        raise
    finally:
        stop.set();thread.join(timeout=1)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    run(json.loads(args.config.read_text()),args.output)


if __name__=='__main__':
    main()
