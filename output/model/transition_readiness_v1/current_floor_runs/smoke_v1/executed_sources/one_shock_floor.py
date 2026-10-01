#!/usr/bin/env python3
"""One permanent 2007 preference shock; current-floor native callbacks only.

This controller does not implement or substitute household or equilibrium math.
The runtime must provide native seed/evaluation/export/render methods described
by RuntimeAdapter below. A diagnostic short-horizon fit is never production.
The original scalar fitter and original annual target builder are reused.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import math
import signal
from contextlib import contextmanager
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
PINNED = HERE / 'pinned_tools'
GATES = dict(market_tolerance=2e-4, fiscal_tolerance=2e-5,
             final_reproduction_tolerance=1e-10,
             stationary_renewal_tolerance=1e-6, terminal_tolerance=1e-3,
             raw_queue_relative_tolerance=1e-3, horizon_relative_tolerance=1e-3)
IDENTITY_KEYS = {'reference_sha256', 'engine_sha256', 'entry_sha256',
                 'grid_sha256', 'effective_parameters_sha256', 'source_pins'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    def plain(x):
        if hasattr(x, 'tolist'):
            return plain(x.tolist())
        if isinstance(x, dict):
            return {k:plain(v) for k,v in x.items()}
        if isinstance(x, (tuple,list)):
            return [plain(v) for v in x]
        if isinstance(x,float) and not math.isfinite(x):
            return None
        return x
    temporary.write_text(json.dumps(plain(value), indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def pinned(item):
    require(isinstance(item, dict) and set(item) == {'path', 'sha256'}, 'Exact file pin required')
    path = Path(item['path'])
    require(path.is_file() and sha(path) == item['sha256'], 'Missing or changed pin: ' + str(path))
    return path


def original_modules():
    # These modules only import standard Python/numpy; no runtime model loading.
    sys.path.insert(0, str(PINNED))
    return (importlib.import_module('run_e5f_preference_estimation'),
            importlib.import_module('e5f_preference_shock_fit'))


def preflight(plan):
    require(plan['schema'] == 'current_floor_one_permanent_v1', 'Wrong plan schema')
    require(plan['kind'] == 'one_permanent' and plan['start_year'] == 2007,
            'Exactly one permanent preference shock begins in 2007')
    require(plan['fiscal_relaxation_authorized'] is True and plan['gates'] == GATES,
            'Only author-authorized dated fiscal relaxation to 2e-5 is permitted')
    require(plan['mode'] in ('production', 'diagnostic'), 'Explicit scientific mode required')
    horizons = plan['horizons']
    require(horizons == [104, 128] if plan['mode'] == 'production' else
            len(horizons) == 2 and all(type(h) is int for h in horizons) and 5 <= horizons[0] < horizons[1],
            'Production requires full 104/128 gates; diagnostic needs two explicit horizons')
    require(set(plan['identity']) == IDENTITY_KEYS, 'Complete current-floor identity required')
    for key in IDENTITY_KEYS - {'source_pins'}:
        require(isinstance(plan['identity'][key], str) and len(plan['identity'][key]) == 64,
                'SHA-256 identity required: ' + key)
    require(plan['identity']['source_pins'], 'Native runtime source pins required')
    pinned(plan['handoff'])
    require(plan['identity']['reference_sha256']==plan['handoff']['sha256'], 'Native reference identity must equal selected handoff pin')
    if plan.get('prepared_native_inputs') is not None:
        preparation=json.loads(pinned(plan['prepared_native_inputs']).read_text())
        validate_preparation(plan,preparation)
    for pin in plan['source_files'].values():
        pinned(pin)
    require({'controller', 'original_estimator', 'original_fitter', 'runtime'}.issubset(plan['source_files']),
            'Controller, retained target/fitter and new runtime source pins required')
    require(plan['source_files']['controller']['sha256'] == sha(__file__), 'Controller source differs')
    require(plan['source_files']['original_estimator']['sha256'] == sha(PINNED/'run_e5f_preference_estimation.py') and
            plan['source_files']['original_fitter']['sha256'] == sha(PINNED/'e5f_preference_shock_fit.py'),
            'Original measurement and scalar fitter source changed')
    for key in ('blocks', 'annual'):
        pinned(plan['target_contract'][key])
    targets = plan['target_contract']['rows']
    require(len(targets) == 4 and [r['decision_year'] for r in targets] == [2007, 2011, 2015, 2019],
            'Complete original four-window target contract required')
    budget = plan['budget']
    for key in ('total_seconds', 'seed_seconds', 'candidate_seconds', 'endpoint_seconds',
                'mapping_seconds', 'path_seconds', 'render_seconds', 'maximum_policy_calls'):
        require(type(budget[key]) in (float, int) and math.isfinite(budget[key]) and budget[key] > 0,
                'Explicit positive finite budget required: ' + key)
    require(type(budget['maximum_policy_calls']) is int, 'Integer policy-call cap required')
    require(plan['seed'] == dict(horizon=12, perturbed_date=5, log_step=1e-5),
            'Fresh current-floor 12-date two-block measured seed required')
    require(len(plan['standard_plot_names']) == 17 and len(set(plan['standard_plot_names'])) == 17,
            'Retained exact 17 standard diagnostic names required')
    require(plan['fit']['max_evaluations'] >= 5 and plan['path']['max_evaluations'] >= 2 and
            plan['endpoint']['max_evaluations'] >= 2, 'Explicit bounded iterations/replay allowance required')
    retained_fit = dict(log_difference_step=.01,fertility_tolerance=.005,max_log_step=.15,
        damping=.7,max_condition_number=1e8,worsening_factor=1.5,reproduction_tolerance=1e-8)
    require(all(plan['fit'].get(k)==v for k,v in retained_fit.items()), 'Retained scalar fit/measurement controls required')
    for name in ('path','endpoint'):
        require(plan[name].get('price_bound_ratios') == [.05,20.] and plan[name].get('max_log_step') == .15 and
                plan[name].get('damping') == .7, 'Retained numerical price domain/step controls required: '+name)
    require(plan['endpoint'].get('slope')==1. and plan['path'].get('pension_bound_ratios')==[.05,20.],
            'Original stationary root initialization and fiscal domain required')
    psi = plan['initial_psi']
    require(math.isfinite(psi) and psi > 0 and plan['psi_bound_ratios'] == [.01, 2.],
            'Use retained positive scalar preference domain')
    conservative_calls = 2+5*2*12+plan['fit']['max_evaluations']*(plan['endpoint']['max_evaluations']+1+
        sum(2*h*plan['path']['max_evaluations'] for h in horizons))
    return dict(status='PASS', native_calls=0, mode=plan['mode'],
                scientific_validation=False, production_ready=False, horizons=horizons,
                conservative_iteration_policy_calls=conservative_calls,
                policy_call_stop_cap=budget['maximum_policy_calls'], total_seconds=budget['total_seconds'],
                fitted_moments=1,free_shock_parameters=1,
                stop_conditions=['deadline','iteration cap','actual policy-call cap','failed native accounting','missing native seed'])


def validate_preparation(plan, preparation):
    require(preparation.get('schema')=='current_floor_native_reference_seed_v1' and
            preparation.get('reference_repeat_verified') is True and preparation.get('fresh_seed_baseline_verified') is True,
            'Completed verified native reference/seed preparation required; no fallback')
    for key in ('identity','source_files','target_contract','seed','gates'):
        require(preparation[key]==plan[key], 'Prepared native inputs differ: '+key)
    for key in ('reference_checkpoint','reconstruction_receipt','measured_seed_receipt','seed_baseline_checks'):
        pinned(preparation[key])
    reconstructed=json.loads(pinned(preparation['reconstruction_receipt']).read_text())
    require(reconstructed['status']=='passed' and reconstructed['identity']==plan['identity'] and
            reconstructed['checkpoint_sha256']==preparation['reference_checkpoint']['sha256'],
            'Prepared reference lacks exact native selected-repeat checkpoint provenance')
    baseline=json.loads(pinned(preparation['seed_baseline_checks']).read_text())
    require(baseline['valid'] is True, 'Native measured seed baseline failed stationary physical gates')
    measured=json.loads(pinned(preparation['measured_seed_receipt']).read_text())
    require(measured['identity']==plan['identity'] and measured['native_measured'] is True and
            measured['horizon']==12 and measured['mapping_count']==5 and measured['perturbation_log_step']==1e-5,
            'Prepared measured seed identity/measurement changed')
    pinned(measured['matrix'])
    return measured


def smoke_readiness(reply):
    gates=dict(accounting=reply['accounting_valid'] is True,stationary=reply['stationary_pass'] is True,
        root=reply['root_pass'] is True,replay=reply['replay_pass'] is True,
        housing=reply['market_maximum_residual']<=GATES['market_tolerance'],
        fiscal=reply['fiscal_maximum_residual']<=GATES['fiscal_tolerance'],
        physical_replay=reply['replay_maximum_gap']<=GATES['final_reproduction_tolerance'])
    passed=all(gates.values())
    return dict(status='PASS' if passed else 'NONPASS',native_setup_verified=passed,native_root_gates=gates,
        diagnostic_terminal_pass=reply['terminal_pass'] is True,
        root_and_terminal_pass=reply['root_and_terminal_pass'] is True,
        scientific_validation=False,production_ready=False)


def production_flags(plan):
    numerical=plan['mode']=='production' and plan['horizons']==[104,128]
    return dict(numerical_path_certified=numerical,
        production_ready=numerical and plan.get('policy_contract_closed') is True)


class RuntimeAdapter:
    """Native floor runtime callback contract; absence fails without fallback.

    identity() returns the exact six-key plan identity. measure_seed(**kwargs)
    measures baseline and +/- log price/pension maps: five native 12-date maps.
    evaluate(**kwargs) returns the original path receipt fields plus accounting,
    stationary and replay physical-gap checks, and native data needed for export.
    export_2023(reply, folder) exports the saved index-4 state and continuation.
    render_standard(reply, folder) returns the exact retained 17 plot file paths.
    Every numerical callback enforces the supplied deadline and policy-call cap.
    """


@contextmanager
def watchdog(deadline):
    remaining = deadline-time.monotonic()
    if remaining <= 0:
        raise TimeoutError('Numerical stage deadline exhausted')
    handler = signal.getsignal(signal.SIGALRM)
    timer = signal.getitimer(signal.ITIMER_REAL)
    def timeout(*_):
        raise TimeoutError('Bounded native numerical call timed out')
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, remaining)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, *timer)
        signal.signal(signal.SIGALRM, handler)


class NativeAdapter:
    """Original numerical orchestration over the new floor engine.

    FloorRuntime must expose reference (native stationary packet), q, psi,
    pension, identity(), stationary(), mapping(), terminal_checks(), export_2023()
    and render_standard(). Native mapping records carry actual policy_calls,
    accounting_valid, gates, original physical residuals, rows and fertility.
    No legacy packet, seed or wealth reconstruction enters this adapter.
    """
    def __init__(self, runtime, plan):
        self.rt, self.plan = runtime, plan
        self.reference = runtime.packet
        self.q = runtime.reference_price
        self.psi = float(runtime.P.psi_child)
        self.pension = float(runtime.P.pension)
        self.calls = 0
        self.warm = {}
        self.endpoints = {}

    def identity(self):
        return self.rt.identity()

    def prepare_native_reference(self, deadline, folder):
        with watchdog(deadline), self.rt.native_budget(deadline,self.plan['budget']['maximum_policy_calls']-self.calls):
            self.reference=self.rt.reconstruct_reference(folder)
        record=json.loads((Path(folder)/'reference_reconstruction.json').read_text())
        require(record['status']=='passed', 'Native selected-reference repeat reconstruction failed')
        require(self.rt.reference_verified is True and self.rt.initial_state is not None,
                'Actual scaled inherited reference state must be reconstructed and verified')
        self.calls+=record['policy_calls']
        checkpoint=Path(folder)/'selected_native_packet.pkl.gz'
        require(sha(checkpoint)==record['checkpoint_sha256'], 'Actual reconstructed native checkpoint differs')
        return dict(accounting_valid=True,policy_calls=record['policy_calls'],
            reference_checkpoint=dict(path=str(checkpoint),sha256=sha(checkpoint)),
            reconstruction_receipt=dict(path=str(Path(folder)/'reference_reconstruction.json'),sha256=sha(Path(folder)/'reference_reconstruction.json')))

    def restore_native_reference(self, preparation, folder):
        checkpoint=pinned(preparation['reference_checkpoint'])
        self.reference=self.rt.restore_reference(dict(preparation['reference_checkpoint'],
            receipt=preparation['reconstruction_receipt']),folder)
        require(self.rt.reference_verified is True and self.rt.initial_state is not None,
                'Pinned native reference restoration must preserve the actual inherited state')
        return dict(accounting_valid=True,policy_calls=0)

    def _account(self, record):
        require(record.get('accounting_valid') is True, 'Stop on failed native accounting')
        calls = record.get('policy_calls')
        require(type(calls) is int and calls >= 0, 'Actual native policy-call count required')
        self.calls += calls
        require(self.calls <= self.plan['budget']['maximum_policy_calls'], 'Native policy-call cap reached')

    def _mapping(self, terminal, endpoint, q, b, psi, folder, deadline, initial_state=None):
        call_deadline=min(deadline,time.monotonic()+self.plan['budget']['mapping_seconds'])
        with watchdog(call_deadline), self.rt.native_budget(call_deadline,self.plan['budget']['maximum_policy_calls']-self.calls):
            reply, record = self.rt.mapping(terminal, endpoint, q, b, psi, folder,
                initial_state=initial_state, start_year=2007)
        self._account(record)
        require(all(math.isfinite(float(v)) for k in ('market_residual', 'fiscal_residual') for v in record[k]),
                'Finite physical residual vectors required')
        write(Path(folder)/'native_record.json', record)
        return reply, record

    def measure_seed(self, *, horizon, perturbed_date, log_step, gates, budget, deadline, folder):
        import numpy as np
        original, _ = original_modules()
        folder = Path(folder)
        start = self.calls
        count = 0
        endpoint = dict(price=self.q, population_scale=self.rt.population_scale)
        def evaluate(q, b):
            nonlocal count
            count += 1
            native, record = self._mapping(self.reference, endpoint, q, b,
                np.full(horizon, self.psi), folder/f'map_{count:03d}', deadline)
            valid = all(record['gates'].values())
            if count == 1:
                terminal = self.rt.terminal_checks(self.reference, endpoint, native,
                    np.full(horizon, self.psi), tolerance=1e-6, raw_queue_tolerance=1e-6)
                valid = valid and terminal['all_checks_pass'] and max(map(abs, record['market_residual'])) <= gates['market_tolerance'] and max(map(abs, record['fiscal_residual'])) <= 1e-6
                write(folder/'baseline_checks.json', dict(valid=valid, terminal=terminal))
            return dict(mapping_valid=valid, market_residual=record['market_residual'], fiscal_residual=record['fiscal_residual'])
        matrix = original.inner.measure_jacobian(evaluate, np.full(horizon,self.q),
            np.full(horizon,self.pension), perturbed_date, log_step, folder/'measured',
            dict(identity=self.identity(), native_measured=True))
        receipt = json.loads((folder/'measured/receipt.json').read_text())
        receipt.update(matrix=matrix, policy_calls=self.calls-start, accounting_valid=True)
        return receipt

    def _endpoint(self, psi, folder, deadline):
        import numpy as np
        _, _ = original_modules()
        from e5f_ssj_scaled_step_root import solve_price_path_scaled
        key=float(psi).hex()
        if key in self.endpoints:
            return self.endpoints[key]
        c = self.plan['endpoint']
        latest = {}
        count = 0
        def evaluate(q):
            nonlocal count
            count += 1
            call_deadline=min(deadline,time.monotonic()+self.plan['budget']['mapping_seconds'])
            with watchdog(call_deadline), self.rt.native_budget(call_deadline,self.plan['budget']['maximum_policy_calls']-self.calls):
                packet, record = self.rt.stationary(psi, float(q[0]), folder/f'point_{count:03d}')
            self._account(record)
            latest.update(packet=packet, record=record)
            write(folder/'latest_completed.json', record)
            return dict(mapping_valid=all(record['gates'].values()), residual=np.array([record['renewal_residual']]))
        root = solve_price_path_scaled(initial_prices=np.array([self.q]), evaluate=evaluate,
            project=lambda q:np.clip(q, self.q*c['price_bound_ratios'][0], self.q*c['price_bound_ratios'][1]),
            slope=c['slope'], market_tolerance=1e-6, max_log_step=c['max_log_step'], damping=c['damping'],
            max_evaluations=c['max_evaluations'], deadline_monotonic=deadline,
            max_condition_number=self.plan['fit']['max_condition_number'], worsening_factor=self.plan['fit']['worsening_factor'],
            final_reproduction_tolerance=1e-10, callback=lambda row:write(folder/'root_progress.json', row))
        write(folder/'root.json', root)
        require(root['converged'], 'Current-floor stationary endpoint failed renewal/replay gates')
        packet, record = latest['packet'], latest['record']
        endpoint = dict(price=record['price'], population_scale=record['population_scale'], psi_child=psi)
        native, check = self._mapping(packet, endpoint,
            np.array([endpoint['price']]), np.array([packet['parameters'].pension]), np.array([psi]), folder/'one_step', deadline,
            initial_state=self.rt.stationary_state(packet,population_scale=endpoint['population_scale']))
        terminal = self.rt.terminal_checks(packet, endpoint, native,
            np.array([psi]), tolerance=1e-6, raw_queue_tolerance=1e-6)
        require(all(check['gates'].values()) and terminal['all_checks_pass'] and
            max(map(abs,check['market_residual'])) <= 2e-4 and max(map(abs,check['fiscal_residual'])) <= 1e-6,
            'Current-floor endpoint native constant-path gate failed')
        endpoint.update(stationary_renewal_gap=abs(record['renewal_residual']), stationary_pass=True)
        self.endpoints[key]=(packet,endpoint)
        return packet, endpoint

    def evaluate(self, *, psi, start_year, horizon, seed, gates, budget, endpoint_controls, path_controls, deadline, folder):
        import numpy as np
        original, _ = original_modules()
        from e5f_four_shock_acceleration import extend_measured_jacobian, solve_joint_with_acceleration
        folder = Path(folder)
        start = self.calls
        terminal, endpoint = self._endpoint(psi, folder/'endpoint', min(deadline,time.monotonic()+budget['endpoint_seconds']))
        p = path_controls
        psi_path = np.full(horizon,psi)
        q0,b0 = self.q,self.pension
        qT,bT = endpoint['price'],terminal['parameters'].pension
        J = extend_measured_jacobian(seed,horizon)
        # Native own-lag slopes are used for the root's emergency reset too.
        slopes = [float(np.median(np.abs(np.diag(J)[i*horizon:(i+1)*horizon]))) for i in range(2)]
        require(all(math.isfinite(s) and s > 0 for s in slopes), 'Native measured own slopes absent; no default fallback')
        warm = self.warm.get(horizon)
        q = np.linspace(q0,qT,horizon) if warm is None else warm['prices']
        b = np.linspace(b0,bT,horizon) if warm is None else warm['fiscal_values']
        latest = {}
        count = 0
        path_deadline = min(deadline,time.monotonic()+budget['path_seconds'])
        def evaluate(q,b):
            nonlocal count
            count += 1
            native,record = self._mapping(terminal,endpoint,q,b,psi_path,folder/f'map_{count:03d}',path_deadline)
            latest.update(native=native,record=record)
            write(folder/'latest_completed.json',record)
            return dict(mapping_valid=all(record['gates'].values()),market_residual=record['market_residual'],fiscal_residual=record['fiscal_residual'])
        root = solve_joint_with_acceleration(closure='fixed_tax',initial_prices=q,initial_fiscal_values=b,evaluate=evaluate,
            project_prices=lambda values:np.clip(values,q0*p['price_bound_ratios'][0],q0*p['price_bound_ratios'][1]),
            fiscal_bounds=[b0*x for x in p['pension_bound_ratios']],market_tolerance=gates['market_tolerance'],
            fiscal_tolerance=gates['fiscal_tolerance'],market_slope=slopes[0],fiscal_slope=slopes[1],
            max_log_step=p['max_log_step'],damping=p['damping'],max_evaluations=p['max_evaluations'],
            deadline_monotonic=path_deadline,max_condition_number=self.plan['fit']['max_condition_number'],
            worsening_factor=self.plan['fit']['worsening_factor'],final_reproduction_tolerance=1e-10,
            initial_jacobian=J if warm is None else warm['final_jacobian'],callback=lambda row:write(folder/'root_progress.json',row))
        original.inner.write(folder/'root.json',root)
        record = latest['record']
        terminal_check = self.rt.terminal_checks(terminal,endpoint,latest['native'],psi_path,
            tolerance=1e-3,raw_queue_tolerance=1e-3)
        valid = bool(root['converged'] and terminal_check['all_checks_pass'])
        if valid:
            self.warm[horizon] = dict(prices=root['final']['prices'],fiscal_values=root['final']['fiscal_values'],final_jacobian=root['final_jacobian'])
        return dict(identity=self.identity(),reference_manifest_sha256=self.identity()['reference_sha256'],
            source_pins=self.identity()['source_pins'],housing=self.rt.housing,
            shock_contract=dict(kind='one_permanent',start_year=2007,psi=psi,expectations='current_shock_permanent_until_next_surprise'),
            psi=psi,horizon=horizon,accounting_valid=True,policy_calls=self.calls-start,
            root_and_terminal_pass=valid,stationary_pass=endpoint['stationary_pass'],
            root_pass=bool(root['converged']),terminal_pass=bool(terminal_check['all_checks_pass']),
            replay_pass=bool(root['gates']['market_replay'] and root['gates']['fiscal_replay']) if 'gates' in root else bool(root['converged']),
            stationary_renewal_gap=endpoint['stationary_renewal_gap'],
            market_maximum_residual=max(map(abs,record['market_residual'])),
            fiscal_maximum_residual=max(map(abs,record['fiscal_residual'])),
            replay_maximum_gap=root['final_reproduction_max_abs'] if root['final_reproduction_max_abs'] is not None else float('inf'),
            terminal=terminal_check,rows=record['rows'],fertility=record['fertility'],
            native_reply=latest['native'],terminal_packet=terminal,endpoint=endpoint,psi_path=psi_path,root=root)

    def export_2023(self, reply, folder):
        return self.rt.export_2023(reply,folder)

    def render_standard(self, reply, folder, deadline):
        with watchdog(deadline):
            return self.rt.render_standard(reply,folder)


class Controller:
    def __init__(self, plan, runtime, output):
        self.plan, self.runtime, self.out = plan, runtime, Path(output)
        self.deadline = time.monotonic() + plan['budget']['total_seconds']
        self.count = 0
        self.policy_calls = 0
        self.last = None

    def remaining(self):
        seconds = self.deadline - time.monotonic()
        if seconds <= 0:
            raise TimeoutError('Total one-shock deadline exhausted')
        return seconds

    def progress(self, phase, **fields):
        self.remaining()
        write(self.out/'heartbeat.json', dict(phase=phase, epoch=time.time(),
              candidate=self.count, policy_calls=self.policy_calls, **fields))

    def account(self, reply):
        require(reply.get('accounting_valid') is True,
                'Native accounting failed; stop immediately, never score fertility')
        calls = reply.get('policy_calls')
        require(type(calls) is int and calls >= 0, 'Runtime must report actual policy-call count')
        self.policy_calls += calls
        require(self.policy_calls <= self.plan['budget']['maximum_policy_calls'], 'Policy-call budget exceeded')

    def prepare(self):
        preflight(self.plan)
        require(self.runtime.identity() == self.plan['identity'], 'Current-floor runtime identity differs')
        original, _ = original_modules()
        contract = original.target_contract(pinned(self.plan['target_contract']['blocks']), pinned(self.plan['target_contract']['annual']))
        require(contract == self.plan['target_contract'], 'Annual-builder complete target fingerprint changed before native calls')
        if self.plan.get('prepared_native_inputs') is not None:
            preparation=json.loads(pinned(self.plan['prepared_native_inputs']).read_text())
            measured=validate_preparation(self.plan,preparation)
            require(hasattr(self.runtime,'restore_native_reference'), 'Native reference restore API required; no fallback')
            self.progress('restore_authenticated_native_inputs')
            self.account(self.runtime.restore_native_reference(preparation,self.out/'native_reference'))
            import numpy as np
            from e5f_four_shock_acceleration import extend_measured_jacobian
            matrix=np.load(pinned(measured['matrix']),allow_pickle=False)
            require(matrix.shape==(24,24) and np.isfinite(matrix).all() and
                    np.array_equal(matrix,extend_measured_jacobian(measured,12)),
                    'Pinned measured matrix does not exactly reconstruct its actual lag profiles')
            self.seed=dict(measured,matrix=matrix,policy_calls=0,accounting_valid=True)
            write(self.out/'native_input_reuse.json',dict(status='PASS',preparation=self.plan['prepared_native_inputs'],
                identity=self.plan['identity'],policy_calls=0,reference_solves=0,seed_mappings=0))
            return
        reference_inputs=None
        if hasattr(self.runtime,'prepare_native_reference'):
            self.progress('selected_native_reference_reconstruction')
            reference_inputs=self.runtime.prepare_native_reference(
                min(self.deadline,time.monotonic()+self.plan['budget']['seed_seconds']),self.out/'native_reference')
            self.account(reference_inputs)
        self.progress('fresh_native_seed')
        self.seed = self.runtime.measure_seed(**self.plan['seed'], gates=self.plan['gates'],
            budget=self.plan['budget'], deadline=min(self.deadline, time.monotonic()+self.plan['budget']['seed_seconds']),
            folder=self.out/'seed')
        self.account(self.seed)
        require(self.seed.get('identity') == self.plan['identity'] and self.seed.get('horizon') == 12 and
                self.seed.get('mapping_count') == 5 and self.seed.get('unknown_blocks') == ['log_house_price', 'log_period_pension'] and
                self.seed.get('native_measured') is True and self.seed.get('matrix') is not None,
                'Fresh measured two-block seed for this exact engine/entry/grid is required')
        write(self.out/'seed_receipt.json', {k:v for k,v in self.seed.items() if k != 'matrix'})
        if reference_inputs is not None:
            def pin(path):return dict(path=str(path),sha256=sha(path))
            write(self.out/'native_preparation.json',dict(schema='current_floor_native_reference_seed_v1',
                **{key:self.plan[key] for key in ('identity','source_files','target_contract','seed','gates')},
                reference_checkpoint=reference_inputs['reference_checkpoint'],
                reconstruction_receipt=reference_inputs['reconstruction_receipt'],
                measured_seed_receipt=pin(self.out/'seed/measured/receipt.json'),
                seed_baseline_checks=pin(self.out/'seed/baseline_checks.json'),
                reference_repeat_verified=True,fresh_seed_baseline_verified=True,
                production_ready=False,shock_fit_verified=False))

    def evaluate(self, psi):
        self.count += 1
        folder = self.out/f'candidate_{self.count:04d}'
        candidate_deadline = min(self.deadline, time.monotonic()+self.plan['budget']['candidate_seconds'])
        self.progress('candidate', psi=psi)
        previous = None
        for horizon in self.plan['horizons']:
            self.remaining()
            reply = self.runtime.evaluate(psi=float(psi), start_year=2007, horizon=horizon,
                seed=self.seed, gates=self.plan['gates'], budget=self.plan['budget'],
                endpoint_controls=self.plan['endpoint'], path_controls=self.plan['path'],
                deadline=candidate_deadline, folder=folder/f'horizon_{horizon:03d}')
            self.account(reply)
            require(reply['identity'] == self.plan['identity'] and reply['psi'] == float(psi) and
                    reply['horizon'] == horizon, 'Native candidate identity/psi/horizon changed')
            valid = (reply['root_and_terminal_pass'] is True and reply['stationary_pass'] is True and
                     reply['market_maximum_residual'] <= GATES['market_tolerance'] and
                     reply['fiscal_maximum_residual'] <= GATES['fiscal_tolerance'] and
                     reply['replay_maximum_gap'] <= GATES['final_reproduction_tolerance'] and
                     reply['stationary_renewal_gap'] <= GATES['stationary_renewal_tolerance'])
            if not valid:
                summary = dict(certified=False, psi=psi, model=None, payload=dict(candidate=self.count, error='native gates failed'))
                write(folder/'failure.json', summary)
                return summary
            require(all(row['calendar_year'] == 2007+4*i for i,row in enumerate(reply['rows'])),
                    'Original four-year calendar/timing changed')
            if previous is not None:
                original, _ = original_modules()
                comparison = original.inner.compare_horizons(previous, reply, 4, GATES['horizon_relative_tolerance'])
                comparison['fertility_absolute_gaps'] = [abs(a['period_tfr_topcode_adjusted']-b['period_tfr_topcode_adjusted'])
                    for a,b in zip(previous['fertility'][:4], reply['fertility'][:4])]
                comparison['passed'] = comparison['passed'] and max(comparison['fertility_absolute_gaps']) <= self.plan['fit']['fertility_tolerance']/5
                write(folder/'horizon_comparison.json', comparison)
                if not comparison['passed']:
                    return dict(certified=False, model=None, payload=dict(candidate=self.count, error='horizon gates failed'))
            previous = reply
        models = [r['period_tfr_topcode_adjusted'] for r in reply['fertility'][:4]]
        require(len(models) == 4 and all(math.isfinite(v) for v in models), 'Four finite original fertility measurements required')
        target = self.plan['target_contract']['rows'][3]['target']
        summary = dict(certified=True, psi=float(psi), model=models[3], gap=models[3]-target,
            loss_contribution=(models[3]-target)**2, payload=dict(candidate=self.count, models=models, path=str(folder)),
            mode=self.plan['mode'], production_ready=False)
        write(folder/'complete.json', summary)
        write(self.out/'latest_completed.json', summary)
        best_path = self.out/'best_so_far.json'
        if not best_path.exists() or summary['loss_contribution'] < json.loads(best_path.read_text())['loss_contribution']:
            write(best_path, summary)
        self.last = reply
        return summary

    def run(self):
        self.prepare()
        original, fitter = original_modules()
        contract = original.target_contract(pinned(self.plan['target_contract']['blocks']), pinned(self.plan['target_contract']['annual']))
        require(contract == self.plan['target_contract'], 'Annual-builder complete target fingerprint changed')
        controls = dict(self.plan['fit'], total_seconds=self.remaining())
        psi = self.plan['initial_psi']
        result = fitter.fit_one(evaluate=self.evaluate, target=contract['rows'][3]['target'], initial_level=psi,
            bounds=[psi*x for x in self.plan['psi_bound_ratios']], controls=controls,
            callback=lambda row:self.progress('scalar_fit', root_row=row))
        write(self.out/'fit.json', result)
        require(result['converged'] and self.last is not None, 'No matched/reproduced one-shock equilibrium; stop')
        require(self.last['psi'] == result['parameter']['estimate'], 'Export candidate differs from reproduced fitted level')
        rows = fitter.fit_rows(contract['rows'], [r['period_tfr_topcode_adjusted'] for r in self.last['fertility'][:4]], 'one_permanent')
        self.out.mkdir(parents=True, exist_ok=True)
        with (self.out/'fertility_fit.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        self.progress('export_exact_2023')
        exported = self.runtime.export_2023(self.last, self.out/'state_2023')
        require(exported.get('calendar_year') == 2023 and exported.get('period_index') == 4 and
                exported.get('exact_native_state') is True and exported.get('reconstructed_or_rescaled') is False and
                exported.get('queue_lags') == [16, 20] and exported.get('forecast_and_continuation_saved') is True,
                'Exact index4 inherited wealth/distribution, both queues, forecast and continuation required')
        self.progress('standard_diagnostics')
        plots = self.runtime.render_standard(self.last, self.out/'diagnostics',
            deadline=min(self.deadline, time.monotonic()+self.plan['budget']['render_seconds']))
        if isinstance(plots,dict):
            require(bool(plots.get('sampled_dates')), 'Native dated standard diagnostic samples required')
            hashes={}
            for dated in plots['sampled_dates']:
                require(set(dated['plots'])==set(self.plan['standard_plot_names']), 'Exact 17 names at each retained native diagnostic date')
                directory=self.out/'diagnostics'/f"date_{dated['period']:03d}"/'standard_diagnostics'
                require(all((directory/name).is_file() and sha(directory/name)==value for name,value in dated['plots'].items()),
                        'Actual retained native dated diagnostic files/hashes required')
                hashes[str(dated['period'])]=dated['plots']
        else:
            require({Path(p).name for p in plots} == set(self.plan['standard_plot_names']) and len(plots) == 17 and
                    all(Path(p).is_file() for p in plots), 'Retain exact actual 17 standard diagnostics')
            hashes={Path(p).name:sha(p) for p in plots}
        receipt = dict(status='matched', mode=self.plan['mode'], scientific_validation=True,
            **production_flags(self.plan), fitted_parameter=result['parameter'],
            actual_policy_calls=self.policy_calls, state_2023=exported, plot_sha256=hashes,
            target_contract=contract, fiscal_relaxation='author-authorized dated gate 2e-5; stationary remains 1e-6')
        write(self.out/'complete.json', receipt)
        return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--native-smoke', action='store_true')
    args = parser.parse_args()
    plan = json.loads(Path(args.plan).read_text())
    receipt = preflight(plan)
    require(not (args.execute and args.native_smoke), 'Select smoke or actual shock fit, not both')
    if args.execute or args.native_smoke:
        require(plan.get('execution_enabled') is True, 'Pinned plan must explicitly enable execution')
        module = importlib.import_module('floor_runtime')
        runtime = NativeAdapter(module.FloorRuntime.from_handoff(plan['handoff'], Path(args.output)/'runtime'),plan)
        controller = Controller(plan, runtime, args.output)
        try:
            if args.native_smoke:
                require(plan['mode']=='diagnostic' and plan['horizons'][0]==6 and plan['path']['max_evaluations']==6,
                        'Exact six-date/six-map smoke must be explicitly diagnostic')
                require(plan.get('native_smoke_psi') == plan['initial_psi']*1.001,
                        'Smoke must explicitly disclose the proposed 0.1-percent preference perturbation')
                controller.prepare()
                started=time.monotonic()
                reply=runtime.evaluate(psi=plan['native_smoke_psi'],start_year=2007,horizon=6,seed=controller.seed,
                    gates=plan['gates'],budget=plan['budget'],endpoint_controls=plan['endpoint'],path_controls=plan['path'],
                    deadline=min(controller.deadline,time.monotonic()+plan['budget']['candidate_seconds']),
                    folder=Path(args.output)/'six_date_smoke')
                controller.account(reply)
                receipt=dict(smoke_readiness(reply),
                    actual_policy_calls=controller.policy_calls,six_date_elapsed_seconds=time.monotonic()-started,
                    identity=runtime.identity(),fresh_twelve_date_seed=True,maximum_maps=6)
                write(Path(args.output)/'smoke_receipt.json',receipt)
            else:
                receipt = controller.run()
        except Exception as exc:
            write(Path(args.output)/'failure.json', dict(status='FAILED', error_type=type(exc).__name__,
                error=str(exc), scientific_validation=False, production_ready=False, actual_policy_calls=controller.policy_calls))
            raise
    write(Path(args.output)/'preflight.json' if not (args.execute or args.native_smoke) else Path(args.output)/'run_receipt.json', receipt)
    print(json.dumps(receipt if not args.execute else {k:receipt[k] for k in ('status','mode','production_ready')}))


if __name__ == '__main__':
    main()
