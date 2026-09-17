#!/usr/bin/env python3
"""Hash-pinned PF household primitives, not an equilibrium or calibration.

Both arms use the same checkpoint parameters and wealth grid. Stationary entry
is retained solely for a constant-policy invariant test; this does not implement
the intended historical person/head closure, terminal equilibrium, or tax rebate.
"""
from __future__ import annotations

import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code/model'), str(ROOT / 'code/model/tools')]
import numpy as np
import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_perfect_foresight_transition as pf
from intergen_eqscale_seq_optimized import solver as model

SCHEMA = 'e5f_matched_pf_primitive_smoke_v1'
FIELDS = ('V', 'c_pol', 'hR_pol', 'bp_pol', 'tenure_choice', 'tenure_probs',
          'loc_probs', 'fert_probs', 'fert_value', 'fert2_probs')
JOINT_FIELDS = ('probabilities', 'wait_probabilities', 'failure_probabilities', 'products')


def same_population_reference(policy, evaluation, parameters):
    """Condition saved native choices on the exact dated inherited population.

    Effective tenure ratios depend on the conditioning pool. Stationary KFE
    and reconstructed pre-choice pools can differ at roundoff, which ratio
    normalization magnifies in almost empty states. Native choices, realized
    masses, and births still retain their separate strict checks.
    """
    arrays = {name: value.copy() for name, value in policy_arrays(policy).items()}
    if bool(getattr(parameters, 'joint_nested_choice', False)):
        post, effective, births, _, _ = calendar.factor_joint_distribution(
            evaluation.g_pre, policy, parameters)
        if (float(np.abs(post - evaluation.g_post_fertility).sum()) > 2e-10
                or abs(float(births.sum()) - float(evaluation.births)) > 2e-10):
            raise RuntimeError('Same-population native choice mass/birth reproduction failed')
        arrays['tenure_probs'] = effective
    return arrays


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def verify(path, expected):
    if digest(path) != expected:
        raise RuntimeError(f'Hash mismatch: {path}')


def load_contract(path, expected):
    verify(path, expected)
    c = json.loads(Path(path).read_text())
    if c['schema'] != SCHEMA:
        raise ValueError('Unknown contract schema')
    for key in ('checkpoint', 'selected_summary'):
        if not Path(c[key]).is_absolute():
            raise ValueError(f'{key} must be absolute')
        verify(c[key], c[key + '_sha256'])
    required = {str(p.relative_to(ROOT)) for p in (ROOT / 'code/model/intergen_eqscale_seq_optimized').glob('*.py')}
    required.update('code/model/tools/' + name for name in (
        'run_e5f_matched_pf_smoke.py', 'run_e5f_perfect_foresight_transition.py',
        'run_dynamic_population_transition.py', 'run_e5f_open_population_transition.py'))
    if not required.issubset(c['source_sha256']):
        raise RuntimeError(f'Missing source pins: {sorted(required - set(c["source_sha256"]))}')
    for relative, expected_hash in c['source_sha256'].items():
        p = Path(relative)
        if p.is_absolute() or '..' in p.parts:
            raise ValueError('Source pins must be repository-relative')
        verify(ROOT / p, expected_hash)
    if not (8 <= int(c['wealth_points']) <= 24 or int(c['wealth_points']) == 120):
        raise ValueError('Use 8–24 diagnostic nodes or the original 120-node grid')
    if not 1 <= int(c['seconds']) <= 840:
        raise ValueError('Stage budget must be at most 840 seconds (15-minute job)')
    if not 0 < float(c['price_relative_change']) <= .02:
        raise ValueError('Changed-price diagnostic must be positive and at most 2%')
    return c


def policy_arrays(policy):
    result = {name: np.asarray(getattr(policy, name)) for name in FIELDS
              if getattr(policy, name, None) is not None}
    if policy.joint_choice is not None:
        result.update({'joint_' + name: np.asarray(getattr(policy.joint_choice, name))
                       for name in JOINT_FIELDS})
    return result


def verify_selected_parameters(P, selected, supply_rule):
    """Check the checkpoint against every saved coordinate and dated psi."""
    theta = selected['best_candidate']['theta']
    checked = {}
    for name, expected in theta.items():
        if name == 'psi_child':
            continue  # The normalized theta intercept is replaced on history.
        actual = np.asarray(getattr(P, name), dtype=float)
        if not np.allclose(actual, expected, rtol=0, atol=1e-12):
            raise RuntimeError(f'Selected/checkpoint parameter mismatch: {name}')
        checked[name] = actual
    expected_psi = float(selected['best_candidate']['new_psi_child'])
    if abs(float(P.psi_child) - expected_psi) > 1e-12:
        raise RuntimeError('Checkpoint is not at the selected end-of-history preference')
    if not np.isclose(float(supply_rule.elasticity), .63, rtol=0, atol=1e-15):
        raise RuntimeError('Expected saved housing-supply rule elasticity=.63')
    return dict(checked_coordinates=checked, dated_psi_child=expected_psi,
                normalized_theta_psi_child=theta.get('psi_child'), supply_elasticity=.63)


def compare_arrays(left, right, tolerance=2e-10):
    if set(left) != set(right):
        raise RuntimeError('Policy array field sets differ')
    result = {}
    for name, value in left.items():
        other = right[name]
        if value.shape != other.shape:
            raise RuntimeError(f'Policy reproduction failed: {name}; shapes {value.shape} vs {other.shape}')
        if not np.allclose(value, other, rtol=0, atol=tolerance, equal_nan=False):
            difference = np.abs(value - other)
            raise RuntimeError(f'Policy reproduction failed: {name}; max_abs={difference.max():.17g}; '
                               f'cells_above_tolerance={int((difference > tolerance).sum())}; tolerance={tolerance}')
        finite = np.isfinite(value) & np.isfinite(other)
        result[name] = float(np.max(np.abs(value[finite] - other[finite]))) if np.any(finite) else 0.
    return result


def save_arrays(path, **arrays):
    np.savez_compressed(path, **arrays)
    with np.load(path, allow_pickle=False) as saved:
        compare_arrays(arrays, {name: saved[name] for name in saved.files}, tolerance=0.)


def dated_budget(evaluation, P, shared, grid, rent):
    """Post-tenure household spending at the actual dated rental price."""
    p, g = evaluation.policy, evaluation.g_current
    bad_mass, largest = 0., 0.
    for age in range(P.J):
        for tenure in range(g.shape[1]):
            for zz, z in enumerate(P.z_grid):
                for parity in range(P.n_parity):
                    for child in range(P.n_child_states):
                        index = (slice(None), tenure, 0, age, zz, parity, child)
                        mass = g[index]
                        if mass.sum() <= 0:
                            continue
                        flat = parity + P.n_parity * child
                        income = model.income_at_state(P, 0, age, float(z))
                        resources = P.R_gross * grid + income
                        grant = float(shared.gb_flat.reshape(-1)[flat])
                        resources += np.clip(grant - (P.R_gross * np.maximum(grid, 0) + income), 0, grant)
                        if tenure == 0:
                            cost = rent * p.hR_pol[index]
                        else:
                            h = P.H_own[tenure - 1]
                            cost = (P.delta + P.tau_H) * p.price[0] * h
                            cost += getattr(P, 'owner_size_cost', 0) * p.price[0] * max(h - getattr(P, 'owner_size_cost_ref', 6), 0) ** getattr(P, 'owner_size_cost_power', 2)
                        gap = p.c_pol[index] + cost + p.bp_pol[index] - resources
                        bad_mass += float(mass[gap > 1e-9].sum())
                        occupied = mass > 1e-12
                        if np.any(occupied):
                            largest = max(largest, float(gap[occupied].max()))
    if bad_mass > 2e-10:
        raise RuntimeError(f'Dated budget gate failed: mass={bad_mass}, excess={largest}')
    return {'budget_excess_mass': bad_mass, 'maximum_occupied_excess': largest,
            'actual_rent': rent, 'budget_tolerance': 1e-9}


def execute(args):
    c = load_contract(args.contract, args.contract_sha256)
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError(f'Refusing nonempty output directory: {out}')
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    progress = {'phase': 'loading', 'arm': args.arm, 'completed': []}
    stop = threading.Event()
    def save(name, value):
        pf.write_json(out / name, value)
    def heartbeat():
        while not stop.wait(15):
            elapsed = time.monotonic() - started
            save('heartbeat.json', dict(progress, elapsed_seconds=elapsed))
            if elapsed > c['seconds']:
                save('failure.json', dict(progress, error='contracted time limit', elapsed_seconds=elapsed))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    def stage(name):
        progress['phase'] = name
        save('heartbeat.json', dict(progress, elapsed_seconds=time.monotonic() - started))
        print(name, flush=True)
        return time.monotonic()
    def done(name, begin, **result):
        result.update(stage=name, elapsed_seconds=time.monotonic() - begin)
        progress['completed'].append(result)
        save('latest_completed.json', result)
        save('summary.json', dict(progress, status='in_progress'))
    transition.configure_sequential_model()
    calendar.model = model
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    opener = gzip.open if Path(c['checkpoint']).suffix == '.gz' else open
    with opener(c['checkpoint'], 'rb') as stream:
        packet = pickle.load(stream)
    P = copy.deepcopy(packet['parameters'])
    selected = json.loads(Path(c['selected_summary']).read_text())
    selected_match = verify_selected_parameters(P, selected, packet['supply_rule'])
    active = args.arm == 'nested'
    P.joint_nested_choice = active
    P.fertility_nest_choice = active
    P.two_shock_choice = False
    P.exhaustive_saving_control = True
    original_transfer = float(getattr(P, 'property_tax_lump_sum_transfer', 0.))
    P.property_tax_lump_sum_transfer = 0.
    original_grid = np.asarray(packet['b_grid'], dtype=float)
    if int(c['wealth_points']) == 120 and len(original_grid) != 120:
        raise RuntimeError('Original-grid pilot requires the saved 120-node grid')
    # Reduced diagnostics take an evenly spaced subset. At 120 nodes this
    # retains the selected grid exactly, including its low-wealth resolution.
    indices = np.unique(np.linspace(0, len(original_grid) - 1, int(c['wealth_points'])).round().astype(int))
    grid = original_grid[indices].copy()
    P.Nb = len(grid)
    price = float(np.asarray(packet['evaluation'].policy.price).reshape(-1)[0])
    if P.I != 1 or P.tenure_choice_kappa != .005:
        raise RuntimeError('Expected selected one-market profile with housing kappa=.005')
    if min(P.kappa_fert, P.kappa_fert_continuation) < P.tenure_choice_kappa:
        raise RuntimeError('Common point violates fertility-nest scale restriction')
    shared = model.precompute_shared(P, grid)
    save('contract.json', dict(c, contract_sha256=args.contract_sha256, arm=args.arm,
        selected_parameter_match=selected_match,
        target_fingerprint_recorded_not_refit=selected['target_fingerprint'],
        transfer_override={'checkpoint': original_transfer, 'diagnostic': 0.},
        flags={key: getattr(P, key) for key in ('joint_nested_choice', 'fertility_nest_choice', 'two_shock_choice', 'exhaustive_saving_control')},
        wealth_grid=grid, original_wealth_points=len(original_grid),
        first_fertility_scale=P.kappa_fert, continuation_fertility_scale=P.kappa_fert_continuation,
        housing_scale=P.tenure_choice_kappa, household_state_shape=[len(grid), P.J, len(P.z_grid), P.n_parity, P.n_child_states],
        expected_bellman_solves=6,
        scope='household PF primitives only; no price iteration or calibration',
        absent_gates=['person_head_closure', 'fiscal_balance', 'terminal_equilibrium', 'historical_fit', 'production_grid_precision'],
        stationary_entry='fixed stationary entrant cohort for invariant test only',
        wall_time_estimate='unmeasured updated arms; capped pilot measures stationary and dated solves'))
    begin = stage('stationary_fixed_price')
    sol = model.solve_markov_income_at_prices(np.array([price]), P, grid, SD=shared)
    stationary = calendar.policy_from_solution(sol, np.array([price]), P, grid, shared)
    reference_arrays = {k: v.copy() for k, v in policy_arrays(stationary).items()}
    pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(sol, stationary, P, grid, shared)
    entry = np.array([float(pre[:, :, :, 0, :, :, :].sum())])
    supply, _ = calendar.normalize_date0_housing_supply(pre, stationary, P, grid, shared, 'static-elastic')
    supply = calendar.HousingSupplyRule(
        mode=supply.mode, initial_price=supply.initial_price,
        initial_stock=supply.initial_stock, elasticity=float(packet['supply_rule'].elasticity))
    save_arrays(out / 'stationary_arrays.npz', **reference_arrays, g_pre=pre, g_current=sol.g, wealth_grid=grid)
    done('stationary_fixed_price', begin, timings=sol.timings, reconstruction=reconstruction)

    def evaluate(policy, state, parameters, sd, rent):
        ev = calendar.evaluate_period(policy.price, state, parameters, grid, sd, calendar.SolveCounter(), supply_rule=supply, supplied_policy=policy)
        if float(ev.feasibility_projection_mass) > 1e-6:
            raise RuntimeError(f'Feasibility projection gate failed: {ev.feasibility_projection_mass}')
        nxt, _, deaths, residual = transition.advance_sequential_calendar_distribution(ev, entry, parameters, grid, sd)
        if abs(residual) > 2e-10 or not np.isfinite(nxt).all() or nxt.min() < -1e-13:
            raise RuntimeError('Incumbent/entry mass accounting gate failed')
        budget = dated_budget(ev, parameters, sd, grid, rent)
        arrays = policy_arrays(policy)
        for name, values in arrays.items():
            if 'prob' in name and (not np.isfinite(values).all() or values.min() < -1e-12 or values.max() > 1 + 1e-12):
                raise RuntimeError(f'Probability range gate failed: {name}')
        return ev, nxt, dict(mass_accounting_error=residual, deaths=deaths,
            projection_mass=ev.feasibility_projection_mass,
            market_residual_diagnostic_only=ev.relative_market_residual, budget=budget)

    begin = stage('constant_continuation')
    rent = float(pf.rents_from_asset_prices([price], price, P)[0])
    const = pf.solve_date_policy(price=price, rent=rent, P=P, b_grid=grid, shared=shared, continuation_V=stationary.V)
    ev, nxt, metrics = evaluate(const, pre, P, shared, rent)
    # The stationary KFE stores distribution-conditioned tenure probabilities.
    # Compare the dated policy after the same conditioning, retaining checks
    # of every native joint probability array as well as realized distributions.
    gaps = compare_arrays(same_population_reference(stationary, ev, P), policy_arrays(ev.policy))
    current_gap = float(np.abs(ev.g_current - sol.g).sum())
    invariant_gap = float(np.abs(nxt - pre).sum())
    birth_gap = abs(float(ev.births) - float(sol.total_births_kfe))
    if current_gap > 2e-8 or invariant_gap > 2e-8 or birth_gap > 2e-10:
        raise RuntimeError(f'Stationary invariant failed: {current_gap}, {invariant_gap}, {birth_gap}')
    save_arrays(out / 'constant_arrays.npz', **policy_arrays(ev.policy), g_pre=pre, g_current=ev.g_current, next_pre=nxt)
    done('constant_continuation', begin, reproduction=gaps, current_l1=current_gap, invariant_l1=invariant_gap, births_absolute=birth_gap, **metrics)

    prices = price * np.array([1 + c['price_relative_change'], 1 + c['price_relative_change'] / 2])
    rents = pf.rents_from_asset_prices(prices, price, P)
    dated, dated_params, snapshots = [None] * 2, [None] * 2, [None] * 2
    continuation = stationary.V
    for t in (1, 0):
        begin = stage(f'backward_{t}')
        dated_params[t] = copy.deepcopy(P)
        sd = model.precompute_shared(dated_params[t], grid)
        dated[t] = pf.solve_date_policy(price=prices[t], rent=rents[t], P=dated_params[t], b_grid=grid, shared=sd, continuation_V=continuation)
        snapshots[t] = {k: v.copy() for k, v in policy_arrays(dated[t]).items()}
        done(f'backward_{t}', begin, price=prices[t], rent=rents[t], next_value_sha256=hashlib.sha256(np.ascontiguousarray(continuation).tobytes()).hexdigest())
        continuation = dated[t].V
    state = pre.copy()
    for t in (0, 1):
        begin = stage(f'forward_{t}')
        params = copy.deepcopy(P)
        sd = model.precompute_shared(params, grid)
        next_value = dated[t + 1].V if t == 0 else stationary.V
        replay = pf.solve_date_policy(price=prices[t], rent=rents[t], P=params, b_grid=grid, shared=sd, continuation_V=next_value)
        reproduction = compare_arrays(snapshots[t], policy_arrays(replay))
        ev, nxt, metrics = evaluate(replay, state, params, sd, rents[t])
        save_arrays(out / f'date_{t}_arrays.npz', **policy_arrays(replay), g_pre=state, g_current=ev.g_current, next_pre=nxt)
        done(f'forward_{t}', begin, reproduction=reproduction, births=ev.births, **metrics)
        state = nxt
    # Later solves must not mutate already saved policy objects.
    compare_arrays(reference_arrays, policy_arrays(stationary), tolerance=0.)
    for t in (0, 1):
        compare_arrays(snapshots[t], policy_arrays(dated[t]), tolerance=0.)
    stop.set()
    save('summary.json', dict(progress, status='passed_primitive_smoke_only',
        elapsed_seconds=time.monotonic() - started, actual_bellman_solves=6,
        market_clearing_tested=False, fiscal_balance_tested=False, person_head_closure_tested=False,
        equilibrium_verified=False, production_ready=False))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--contract', required=True, type=Path)
    p.add_argument('--contract-sha256', required=True)
    p.add_argument('--arm', required=True, choices=('sequential', 'nested'))
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    preexisting_output = args.output.exists() and any(args.output.iterdir())
    try:
        execute(args)
    except Exception as error:
        # Do not overwrite any existing run if startup rejected its directory.
        if not preexisting_output and args.output.exists() and (args.output / 'contract.json').exists():
            pf.write_json(args.output / 'failure.json', {'error': repr(error), 'arm': args.arm})
        raise


if __name__ == '__main__':
    main()
