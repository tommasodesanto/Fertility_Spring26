"""Pinned, diagnostic terminal PAYGO root from the verified revised initial state.

Required JSON contract (no launch or economic defaults): schema; source_sha256;
initial_checkpoint, initial_summary, initial_contract (each {path, sha256});
originating_source_commit; psi_change_from_initial; preference_rule;
demographic_seed_mode; terminal_demographic_year; seconds; root_seconds;
identity_tolerance; initial_fertility_tolerance; endpoint_controls;
audit_controls; root_controls; standard_graph_count.

The reviewed adapter owns the exact bounded root loop and reserves fresh replay
inside its evaluation budget. This driver never rebuilds demographics, solves
an additional GE for diagnostics, normalizes preferences, or promotes a fit.
The seed pickle is trusted only after its fixed approved SHA256 is verified.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import fields
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import pickle
import sys
import threading
import time

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
os.environ.setdefault('MPLBACKEND', 'Agg')

ROOT = Path(__file__).resolve().parents[3]
SCHEMA = 'e5f_balanced_terminal_probe_v1'
INITIAL_SHA256 = '4afc7fc6f4db32a1bb220bb3b2b30c6e2c0b96822a152c294f6c73c228320852'
ROOT_KEYS = {'price_bounds', 'pension_bounds', 'fiscal_tolerance', 'market_slope',
    'fiscal_slope', 'max_log_step', 'damping', 'max_evaluations',
    'max_condition_number', 'worsening_factor', 'final_reproduction_tolerance',
    'initial_jacobian'}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def verify(path, expected):
    if (not isinstance(expected, str) or len(expected) != 64
            or digest(path) != expected):
        raise ValueError(f'SHA256 mismatch: {path}')


def pinned_json(spec):
    if set(spec) != {'path', 'sha256'} or not Path(spec['path']).is_absolute():
        raise ValueError('Pinned inputs require explicit absolute path and SHA256')
    verify(spec['path'], spec['sha256'])
    return json.loads(Path(spec['path']).read_text())


def validate_contract(c):
    required = {'schema', 'source_sha256', 'initial_checkpoint', 'initial_summary',
        'initial_contract', 'originating_source_commit', 'psi_change_from_initial',
        'preference_rule', 'demographic_seed_mode', 'terminal_demographic_year',
        'seconds', 'root_seconds', 'identity_tolerance', 'initial_fertility_tolerance',
        'endpoint_controls', 'audit_controls', 'root_controls', 'standard_graph_count'}
    if set(c) != required or c['schema'] != SCHEMA:
        raise ValueError('Explicit terminal-probe contract fields/schema required')
    if (c['originating_source_commit'] != 'c6dd3508'
            or c['preference_rule'] != 'diagnostic_constant_terminal_intercept_not_fitted'
            or c['demographic_seed_mode'] != 'serialized_frozen_2023_primitives_without_realignment'
            or c['standard_graph_count'] != 17
            or not math.isfinite(c['psi_change_from_initial'])):
        raise ValueError('Approved diagnostic origin, preference and frozen-demography declarations required')
    for name in ('seconds', 'root_seconds', 'terminal_demographic_year'):
        if type(c[name]) is not int:
            raise ValueError(f'{name} must be an explicit integer')
    if (not 1 <= c['seconds'] <= 1800 or not 1 <= c['root_seconds'] <= c['seconds'] - 60
            or c['terminal_demographic_year'] < 2023):
        raise ValueError('First smoke is at most 30 minutes with at least 60 seconds reserved after the root')
    for name, ceiling in (('identity_tolerance', 2e-9), ('initial_fertility_tolerance', 5e-4)):
        if not math.isfinite(c[name]) or not 0 < c[name] <= ceiling:
            raise ValueError(f'Explicit {name} must be positive and no looser than {ceiling}')
    r = c['root_controls']
    if set(r) != ROOT_KEYS or type(r['max_evaluations']) is not int or not 2 <= r['max_evaluations'] <= 8:
        raise ValueError('Explicit root controls and 2–8 mappings including replay required')
    pin = c['initial_checkpoint']
    if (set(pin) != {'path', 'sha256'} or not Path(pin['path']).is_absolute()
            or pin['sha256'] != INITIAL_SHA256
            or Path(pin['path']).parts[-3:] != ('new_balanced_smoke', 'repetition_02', 'initial_state.pkl.gz')):
        raise ValueError('Require the approved c6dd3508 repetition_02 normalized initial checkpoint')


def verify_sources(c):
    required = {str(p.relative_to(ROOT)) for p in (ROOT / 'code/model').rglob('*.py')}
    if not required or not required.issubset(c['source_sha256']):
        raise ValueError('Source manifest must pin every code/model Python source including this driver')
    for relative, pin in c['source_sha256'].items():
        path = (ROOT / relative).resolve()
        if Path(relative).is_absolute() or not path.is_relative_to(ROOT):
            raise ValueError('Source paths must be snapshot-relative and stay within it')
        verify(path, pin)
    return len(required)


def validate_demographics(d, P, c, person, np):
    """Validate the serialized class and its own 2023 alignment; never rebuild it."""
    if type(d) is not person.AnnualDemographicPrimitives:
        raise ValueError('Missing or incompatible serialized AnnualDemographicPrimitives')
    if (set(vars(d)) != {f.name for f in fields(person.AnnualDemographicPrimitives)}
            or d.start_year != 2023 or d.last_empirical_year != c['terminal_demographic_year']
            or not math.isfinite(d.scale_model_units_per_person) or d.scale_model_units_per_person <= 0
            or type(d.initial_person_state) is not person.CohortState
            or d.initial_person_state.year != 2023):
        raise ValueError('Serialized demographic schema, dates or scale differ from the explicit contract')
    state = d.initial_person_state.validated(tolerance=c['identity_tolerance'])
    people, heads = np.asarray(state.persons), np.asarray(state.heads)
    rates, raw = np.asarray(d.headship_rates), np.asarray(d.raw_acs_headship_rates)
    if (people.shape[0] != 2 or people.sum() <= 0 or heads.sum() <= 0
            or rates.shape != people.shape or raw.shape != people.shape
            or not np.isfinite(rates).all() or not np.isfinite(raw).all()
            or np.any((rates < 0) | (rates > 1)) or np.any((raw < 0) | (raw > 1))):
        raise ValueError('Invalid serialized sex/age population or headship rates')
    age_start, width = int(P.age_start), int(P.da)
    if age_start != P.age_start or width != P.da or width != 4:
        raise ValueError('Demographic support must use the exact four-year household age cells')
    stop = age_start + width * int(P.J)
    outside = rates.copy()
    outside[:, age_start:stop] = 0.
    if people.shape[1] < stop or np.max(np.abs(outside)) > 1e-14:
        raise ValueError('Serialized terminal headship extends outside household-model age support')
    gap = float(np.max(np.abs(heads - people * rates)))
    age_mass = person.aggregate_heads_to_model_age_cells(state, age_start=age_start,
        cell_width=width, number_of_cells=int(P.J))
    expected_rates = raw.copy()
    for name in ('initial_household_age_mass', 'initial_person_head_age_mass',
                 'model_age_headship_alignment_factors'):
        a = np.asarray(getattr(d, name))
        if a.shape != (int(P.J),) or not np.isfinite(a).all() or np.any(a < 0):
            raise ValueError(f'Invalid serialized {name}')
    for j, factor in enumerate(d.model_age_headship_alignment_factors):
        expected_rates[:, age_start+j*width:age_start+(j+1)*width] *= factor
    gaps = dict(person_head=gap,
        initial_household_age=float(np.max(np.abs(age_mass - d.initial_household_age_mass))),
        initial_person_head_age=float(np.max(np.abs(age_mass - d.initial_person_head_age_mass))),
        headship_alignment=float(np.max(np.abs(expected_rates - rates))))
    if any(value > c['identity_tolerance'] for value in gaps.values()):
        raise ValueError(f'Serialized demographic internal identity failed: {gaps}')
    signatures = {}
    for name in ('birth_sex_shares', 'survival', 'net_migration'):
        values = getattr(d, name)
        if not isinstance(values, dict) or d.last_empirical_year not in values:
            raise ValueError(f'Serialized terminal-year {name} is missing')
        h = hashlib.sha256()
        for year, value in sorted(values.items()):
            a = np.asarray(value, dtype=float)
            shape = (2,) if name == 'birth_sex_shares' else people.shape
            if type(year) is not int or a.shape != shape or not np.isfinite(a).all():
                raise ValueError(f'Invalid serialized {name} year/shape/values')
            if name != 'net_migration' and np.any((a < 0) | (a > 1)):
                raise ValueError(f'Serialized {name} rates must lie in [0,1]')
            if name == 'birth_sex_shares' and abs(float(a.sum()) - 1.) > c['identity_tolerance']:
                raise ValueError('Birth sex shares do not sum to one')
            h.update(str(year).encode()); h.update(np.ascontiguousarray(a).tobytes())
        signatures[name] = h.hexdigest()
    sources = {'population_mid', 'births_mid', 'survival', 'vintage_2025', 'acs_headship'}
    if set(d.source_paths) != sources:
        raise ValueError('Serialized demographic provenance fields are incomplete')
    return dict(class_name=f'{type(d).__module__}.{type(d).__name__}',
        start_year=d.start_year, terminal_year=d.last_empirical_year,
        scale_model_units_per_person=d.scale_model_units_per_person,
        internal_identity_gaps=gaps, annual_array_sha256=signatures,
        serialized_source_paths={k: str(v) for k, v in d.source_paths.items()},
        source_pinning='serialized arrays are covered by the approved initial checkpoint SHA256',
        status=c['demographic_seed_mode'], initial_stationary_age_realignment_performed=False)


def run(args):
    verify(args.contract, args.contract_sha256)
    c = json.loads(args.contract.read_text())
    validate_contract(c)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    finished, lock = threading.Event(), threading.Lock()
    state = dict(phase='preflight', completed_evaluations=0)

    def jsonable(value):
        if isinstance(value, dict): return {str(k): jsonable(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)): return [jsonable(v) for v in value]
        if hasattr(value, 'tolist'): return jsonable(value.tolist())
        if isinstance(value, Path): return str(value)
        return value

    def save(name, value):
        with lock:
            path = out / name
            temporary = path.with_suffix(path.suffix + '.tmp')
            temporary.write_text(json.dumps(jsonable(value), indent=2, sort_keys=True) + '\n')
            os.replace(temporary, path)

    def watchdog():
        while not finished.wait(min(60., max(.01, c['seconds'] - (time.monotonic() - started)))):
            elapsed = time.monotonic() - started
            save('heartbeat.json', dict(state, elapsed_seconds=elapsed))
            if elapsed >= c['seconds']:
                save('failure.json', dict(state, status='timeout', elapsed_seconds=elapsed,
                    message='Explicit total wall-time budget exhausted; no automatic retry'))
                os._exit(124)

    threading.Thread(target=watchdog, daemon=True).start()
    save('contract.json', dict(c, contract_sha256=args.contract_sha256))
    save('heartbeat.json', state)
    try:
        verified_sources = verify_sources(c)
        initial_summary = pinned_json(c['initial_summary'])
        initial_contract = pinned_json(c['initial_contract'])
        verify(c['initial_checkpoint']['path'], INITIAL_SHA256)
        sys.path[:0] = [str(ROOT / 'code/model'), str(ROOT / 'code/model/tools')]
        import numpy as np
        import e5f_balanced_terminal as adapter
        import run_e5f_matched_pf_smoke as primitive
        import run_e5f_perfect_foresight_person_demography as person
        from e5f_stationary_paygo import certify_initial_pension
        from e5f_social_security import fiscal_accounts
        from e5f_parenthood_utility import validate_parenthood_utility
        with gzip.open(c['initial_checkpoint']['path'], 'rb') as stream:
            inherited = pickle.load(stream)
        required_packet = {'parameters', 'b_grid', 'evaluation', 'shared', 'supply_rule',
            'solution', 'stationary_g_pre', 'demographic_seed', 'contract_sha256'}
        if not isinstance(inherited, dict) or set(inherited) != required_packet:
            raise ValueError('Approved initial checkpoint packet fields differ from their source schema')
        if (inherited['contract_sha256'] != c['initial_contract']['sha256']
                or initial_contract['schema'] != 'e5f_parenthood_initial_probe_v1'
                or initial_contract['normalize'] is not True
                or initial_summary['status'] != 'passed_initial_diagnostic'
                or initial_summary['calibrated_smm'] is not False
                or initial_summary['checkpoint_sha256'] != INITIAL_SHA256):
            raise ValueError('Verified initial receipt, source contract or checkpoint parent identity differs')
        P = copy.deepcopy(inherited['parameters'])
        grid = np.asarray(inherited['b_grid']).copy()
        validate_parenthood_utility(P)
        if (P.I != 1 or P.J != 17 or P.Nb != 120 or len(grid) != 120
                or not np.array_equal(grid, primitive.model.make_grid(P))):
            raise ValueError('Terminal probe requires the exact 120-node, 17-age initial grid')
        normalization = initial_summary['normalization']
        if (normalization['status'] not in ('derived_intercept', 'normalized_at_initial_guess')
                or normalization['target'] != 2.1 or normalization['psi_child'] != float(P.psi_child)
                or not math.isfinite(normalization['completed_fertility'])
                or abs(normalization['completed_fertility'] - 2.1) > c['initial_fertility_tolerance']):
            raise ValueError('Initial normalization does not certify this unchanged preference intercept')
        initial_fiscal = certify_initial_pension(inherited['evaluation'].g_current, P,
            marginal_tolerance=c['identity_tolerance'], fiscal_tolerance=c['root_controls']['fiscal_tolerance'])
        demographics = inherited['demographic_seed']
        demographic_receipt = validate_demographics(demographics, P, c, person, np)
        supply = inherited['supply_rule']
        start_price = float(inherited['evaluation'].policy.price[0])
        if (start_price != float(initial_summary['price'])
                or not math.isfinite(start_price) or start_price <= 0
                or not math.isfinite(float(P.pension)) or P.pension <= 0):
            raise ValueError('Initial price/pension starting guesses do not match the verified initial state')
        start_pension = float(P.pension)
        old_psi = float(P.psi_child)
        P.psi_child = old_psi + float(c['psi_change_from_initial'])
        controls = adapter.EndpointControls(**c['endpoint_controls'])
        audit_controls = adapter.TerminalAuditControls(**c['audit_controls'])
        if (set(c['endpoint_controls']) != {f.name for f in fields(adapter.EndpointControls)}
                or set(c['audit_controls']) != {f.name for f in fields(adapter.TerminalAuditControls)}):
            raise ValueError('Every endpoint and household audit control must be explicit')
        adapter._validate_inputs(P, grid, demographics, supply, controls, audit_controls,
                                 c['root_controls']['fiscal_tolerance'])
        save('preflight.json', dict(status='passed', source_files_verified=verified_sources,
            initial_normalization=normalization, initial_fiscal=initial_fiscal,
            demographic_seed=demographic_receipt, initial_psi=old_psi, terminal_psi=P.psi_child,
            psi_change_from_initial=c['psi_change_from_initial'], preference_rule=c['preference_rule'],
            root_starting_guesses=dict(asset_price=start_price, pension_period=start_pension,
                status='verified initial equilibrium values; not a solved terminal state'),
            supply_rule=vars(supply), calibrated_smm=False, production_eligible=False))
        save('sizing.json', dict(maximum_root_mappings=c['root_controls']['max_evaluations'],
            maximum_stationary_household_solves=c['root_controls']['max_evaluations'],
            maximum_person_fixed_point_iterations_per_mapping=controls.maximum_inner_iterations,
            root_budget_seconds=c['root_seconds'], total_budget_seconds=c['seconds'],
            reporting_reserve_seconds=c['seconds']-c['root_seconds'],
            prior_initial_single_ge_seconds_approximate=70.,
            historical_terminal_reference=dict(elapsed_seconds=456.4692804738879, evaluations=5,
                source='output/model/e5f_person_demography_terminal_root_rebated-tax1-baseline_20260826a_production/result/summary.json',
                comparability='different utility, fiscal root and source; not a revised runtime forecast'),
            revised_terminal_runtime='unknown; exact-loop smoke measures it',
            extra_reporting_bellman_or_ge_solves=0,
            reporting_work='one supplied-policy terminal distribution evaluation, checkpoint reload, 17 standard graphs'))
        del inherited
        records = []
        save('latest_completed.json', dict(status='no_completed_root_mapping'))
        save('best_so_far.json', dict(status='no_valid_root_mapping'))
        def callback(record):
            if record.get('event') == 'complete':
                save('root_completion.json', record)
                return
            if record.get('event') == 'safeguard':
                save('latest_safeguard.json', record)
                return
            records.append(record)
            state.update(phase='terminal_root', completed_evaluations=record['evaluation'])
            save('latest_completed.json', record)
            save('root_evaluations.json', records)
            if record.get('new_best'): save('best_so_far.json', record)
            save('heartbeat.json', dict(state, elapsed_seconds=time.monotonic()-started))
        root_deadline = started + c['root_seconds']
        if time.monotonic() >= root_deadline:
            raise TimeoutError('Preflight exhausted the explicit root time allocation')
        state['phase'] = 'terminal_root'
        result = adapter.solve_balanced_terminal(parameters=P, b_grid=grid,
            demographic_primitives=demographics, supply_rule=supply, controls=controls,
            audit_controls=audit_controls, start_price=start_price,
            start_pension_period=start_pension, deadline_monotonic=root_deadline,
            callback=callback, **c['root_controls'])
        save('root_receipt.json', result.root_receipt)
        selected = result.endpoint
        final_checks = dict(endpoint_returned=selected is not None,
            root_fresh_replay=result.production_eligible, checkpoint_reload=False,
            final_accounts_reproduced=False, standard_graphs=False)
        checkpoint_hash = None
        if selected is not None:
            state['phase'] = 'final_checkpoint_and_diagnostics'
            final_P = selected.parameters
            shared = primitive.model.precompute_shared(final_P, grid)
            counter = primitive.calendar.SolveCounter()
            evaluation = primitive.calendar.evaluate_period(selected.policy.price,
                selected.fixed_point.g_pre, final_P, grid, shared, counter,
                supply_rule=supply, supplied_policy=selected.policy)
            if counter.total != 0:
                raise RuntimeError('Final diagnostic reconstruction unexpectedly solved households')
            actual = fiscal_accounts(evaluation.g_current, final_P)
            if actual != selected.social_security:
                raise RuntimeError('Final supplied-policy fiscal accounts differ from the returned terminal')
            final_checks['final_accounts_reproduced'] = True
            packet = dict(schema='e5f_balanced_terminal_probe_checkpoint_v1',
                parameters=final_P, b_grid=grid, policy=selected.policy,
                fixed_point=selected.fixed_point, evaluation=evaluation, shared=shared,
                supply_rule=supply, demographic_seed=demographics,
                endpoint=selected.endpoint, social_security=actual,
                root_receipt=result.root_receipt, contract_sha256=args.contract_sha256,
                production_eligible=False, calibrated_smm=False)
            checkpoint = out / 'terminal_state.pkl.gz'
            with gzip.open(checkpoint, 'wb', compresslevel=1) as stream:
                pickle.dump(packet, stream, protocol=5)
            with gzip.open(checkpoint, 'rb') as stream:
                replay = pickle.load(stream)
            for name, value in primitive.policy_arrays(packet['policy']).items():
                np.testing.assert_array_equal(value, primitive.policy_arrays(replay['policy'])[name])
            for name in ('g_pre', 'g_post_fertility', 'g_current'):
                np.testing.assert_array_equal(getattr(evaluation, name), getattr(replay['evaluation'], name))
            np.testing.assert_array_equal(replay['fixed_point'].persons.persons, selected.fixed_point.persons.persons)
            np.testing.assert_array_equal(replay['fixed_point'].persons.heads, selected.fixed_point.persons.heads)
            if (replay['contract_sha256'] != args.contract_sha256
                    or replay['parameters'].pension != final_P.pension
                    or fiscal_accounts(replay['evaluation'].g_current, replay['parameters']) != actual):
                raise RuntimeError('Terminal checkpoint reload fiscal/contract gate failed')
            del replay
            checkpoint_hash = digest(checkpoint)
            final_checks['checkpoint_reload'] = True
            save('checkpoint_receipt.json', dict(sha256=checkpoint_hash, checks=final_checks,
                actual_social_security=actual, production_eligible=False))
            import run_e5f_independent_numerical_audit as audit
            audit.standard_diagnostics(packet, out, validate_production_young=False)
            graph_paths = sorted((out / 'standard_diagnostics').glob('*.png'))
            final_checks['standard_graphs'] = len(graph_paths) == 17
            save('diagnostics_receipt.json', dict(standard_graph_count=len(graph_paths),
                graph_sha256={str(p.relative_to(out)): digest(p) for p in graph_paths},
                label='terminal stationary diagnostic; inherited lifecycle_2023.csv filename is not a 2023 state'))
        passed = all(final_checks.values())
        summary = dict(status='passed_terminal_root_diagnostic' if passed else 'failed_terminal_root_diagnostic',
            endpoint_numerically_verified=passed, final_checks=final_checks,
            root_status=result.root_receipt['status'], root_evaluations=result.root_receipt['evaluations'],
            checkpoint_sha256=checkpoint_hash, terminal_psi=P.psi_child,
            demographic_seed_mode=c['demographic_seed_mode'],
            production_eligible=False, calibrated_smm=False, historical_path_solved=False,
            interpretation='constant-preference terminal smoke from inherited frozen demographics; no fitted preference path',
            elapsed_seconds=time.monotonic()-started)
        save('summary.json', summary)
        return 0 if passed else 2
    except Exception as exc:
        save('failure.json', dict(state, status='failed_terminal_probe',
            error_type=type(exc).__name__, error=str(exc), elapsed_seconds=time.monotonic()-started,
            production_eligible=False, calibrated_smm=False))
        raise
    finally:
        finished.set()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    return run(parser.parse_args())


if __name__ == '__main__':
    raise SystemExit(main())
