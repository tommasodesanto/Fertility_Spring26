"""Conditional historical/person PF composition; no equilibrium or calibration.

The caller supplies the full dated price/preference/transfer paths, a terminal
value boundary, and frozen 2023 person/head primitives. No roots or empirical
targets are inferred here. The 2023 household decision is evaluated only once.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Callable, Sequence

import numpy as np

import run_e5f_perfect_foresight_transition as pf
import run_e5f_perfect_foresight_person_demography as person_pf


@dataclass
class ConditionalHistoryEvaluation:
    history: pf.PathEvaluation
    person_tail: person_pf.PersonPathEvaluation
    rows: list[dict[str, Any]]
    values: list[np.ndarray]
    bellman_solves: int
    initial_2023_age_head_gap: float
    scope: str = "Conditional PF path with supplied terminal boundary; not a converged equilibrium"


def evaluate_history_and_person_tail(
    *, years: Sequence[int], prices: Sequence[float], psi_path: Sequence[float],
    transfer_path: Sequence[float], terminal_price: float, terminal_V: np.ndarray,
    base_parameters: Any, b_grid: np.ndarray, initial_state: pf.PFInitialState,
    historical_conditioning: pf.HistoricalConditioning,
    initial_2023_persons: person_pf.CohortState,
    demographic_primitives: person_pf.AnnualDemographicPrimitives,
    supply_rule: Any, birth_to_entry_conversion: float,
    observer: Callable | None = None,
) -> ConditionalHistoryEvaluation:
    """Evaluate 2007--2019 history followed by the person tail from 2023.

    The optional observer receives the global date index, current evaluation,
    dated parameters, wealth grid, and shared inputs before each advancement.
    """
    dates = np.asarray(years)
    p = np.asarray(prices, dtype=float)
    psi = np.asarray(psi_path, dtype=float)
    transfers = np.asarray(transfer_path, dtype=float)
    if (dates.ndim != 1 or len(dates) < 5
            or not np.array_equal(dates, 2007 + 4 * np.arange(len(dates)))
            or float(base_parameters.period_years) != 4.):
        raise ValueError("Joined PF evaluation requires complete four-year dates from 2007 through at least 2023")
    if (any(a.shape != dates.shape or not np.isfinite(a).all() for a in (p, psi, transfers))
            or np.any(p <= 0) or np.any(transfers < 0)):
        raise ValueError("Explicit price, preference and nonnegative transfer paths must match all dates")
    if historical_conditioning.start_year != 2007:
        raise ValueError("Joined history must begin in 2007")
    historical_conditioning.validate(base_parameters, 4, initial_state, birth_to_entry_conversion)
    if observer is not None and not callable(observer):
        raise ValueError("Joined observer must be callable")
    if (observer is not None and historical_conditioning.observer is not None
            and historical_conditioning.observer is not observer):
        raise ValueError("Supply one common observer, not different historical and joined observers")
    callback = observer if observer is not None else historical_conditioning.observer
    people = initial_2023_persons.validated()
    frozen_people = demographic_primitives.initial_person_state
    if (people.year != 2023 or frozen_people.year != 2023
            or not np.array_equal(people.persons, frozen_people.persons)
            or not np.array_equal(people.heads, frozen_people.heads)):
        raise ValueError("Initial 2023 persons/heads must match the supplied frozen demographic primitives")

    tail_prices, tail_psi, tail_transfers = p[4:], psi[4:], transfers[4:]
    tail_rents = pf.rents_from_asset_prices(tail_prices, terminal_price, base_parameters)
    tail_values, tail_backward_solves = pf.backward_value_path(
        prices=tail_prices, rents=tail_rents, psi_path=tail_psi,
        terminal_V=terminal_V, base_parameters=base_parameters,
        b_grid=b_grid, transfer_path=tail_transfers,
    )
    history = pf.evaluate_path_at_prices(
        prices=p[:4], psi_path=psi[:4], transfer_path=transfers[:4],
        terminal_price=float(p[4]), terminal_V=tail_values[0],
        base_parameters=base_parameters, b_grid=b_grid, initial_state=initial_state,
        supply_rule=supply_rule, birth_to_entry_conversion=birth_to_entry_conversion,
        historical_conditioning=replace(historical_conditioning, observer=callback),
    )
    g_2023 = history.terminal_state.g_pre
    heads_2023 = person_pf.aggregate_heads_to_model_age_cells(
        people, age_start=int(base_parameters.age_start),
        cell_width=int(base_parameters.da), number_of_cells=int(base_parameters.J),
    )
    age_gap = float(np.max(np.abs(g_2023.sum(axis=(0, 1, 2, 4, 5, 6)) - heads_2023)))
    if not np.isfinite(age_gap) or age_gap > 2e-9:
        raise RuntimeError(f"Historical 2023 household/person head-age identity fails: {age_gap}")

    def tail_observer(period, evaluation, parameters, grid, shared):
        callback(period + 4, evaluation, parameters, grid, shared)

    tail = person_pf.evaluate_path_at_prices_person_demography(
        prices=tail_prices, psi_path=tail_psi, transfer_path=tail_transfers,
        terminal_price=terminal_price, terminal_V=terminal_V,
        base_parameters=base_parameters, b_grid=b_grid,
        initial_state=person_pf.PersonPFState(g_pre=g_2023.copy(), persons=people),
        demographic_primitives=demographic_primitives, supply_rule=supply_rule,
        precomputed_value_path=tail_values,
        observer=tail_observer if callback is not None else None,
    )
    rows = [dict(row) for row in history.rows]
    rows.extend(dict(row, period=int(row['period']) + 4) for row in tail.rows)
    if [row['calendar_year'] for row in rows] != dates.tolist():
        raise RuntimeError("Joined path duplicated or omitted a calendar date")
    count = history.bellman_solves + tail_backward_solves + tail.bellman_solves
    if count != 2 * len(dates):
        raise RuntimeError("Joined path repeated or omitted a backward/forward household solve")
    return ConditionalHistoryEvaluation(
        history=history, person_tail=tail, rows=rows,
        values=history.values[:-1] + tail.values, bellman_solves=count,
        initial_2023_age_head_gap=age_gap,
    )


HISTORY_SMOKE_SCHEMA = 'e5f_matched_history_smoke_v1'


def load_smoke_contract(path, expected_sha256, arm, *, maximum_seconds=840):
    """Validate all current source/data pins without invoking a model solve."""
    import json
    from pathlib import Path
    import run_e5f_matched_pf_smoke as primitive
    primitive.verify(path, expected_sha256)
    c = json.loads(Path(path).read_text())
    if (type(maximum_seconds) is not int or not 1 <= maximum_seconds <= 21600
            or c['schema'] != HISTORY_SMOKE_SCHEMA
            or not 1 <= int(c['seconds']) <= maximum_seconds):
        raise ValueError(f'Expected joined smoke schema and at most {maximum_seconds} seconds')
    for name in ('checkpoint', 'selected_summary', 'stationary_arrays', 'primitive_summary', 'primitive_contract'):
        if not Path(c[name]).is_absolute():
            raise ValueError(f'{name} must be absolute')
        primitive.verify(c[name], c[name + '_sha256'])
    root = Path(__file__).resolve().parents[3]
    required = {str(p.relative_to(root)) for package in
                ('intergen_eqscale_seq_optimized', 'demographic_transition')
                for p in (root / 'code/model' / package).glob('*.py')}
    required.update('code/model/tools/' + name for name in (
        'run_e5f_matched_pf_history.py', 'run_e5f_matched_pf_smoke.py',
        'run_e5f_perfect_foresight_transition.py', 'run_e5f_perfect_foresight_person_demography.py',
        'run_dynamic_population_transition.py', 'run_e5f_open_population_transition.py',
        'build_e5f_coherent_person_cohort_path.py', 'build_e5f_persons_demographic_satellite.py'))
    if not required.issubset(c['source_sha256']):
        raise ValueError(f'Missing source pins: {sorted(required - set(c["source_sha256"]))}')
    for relative, expected in c['source_sha256'].items():
        p = Path(relative)
        if p.is_absolute() or '..' in p.parts:
            raise ValueError('Source pins must be repository-relative')
        primitive.verify(root / p, expected)
    files = dict(population_mid='population_mid.csv', births_mid='births_mid.csv',
                 survival='survival.csv', vintage_2025='vintage_2025_age_sex.csv',
                 acs_headship='acs_headship_profiles.csv')
    if set(c['demographic_sources']) != set(files):
        raise ValueError('Exactly five demographic source pins are required')
    for name, filename in files.items():
        row = c['demographic_sources'][name]
        p = Path(row['path'])
        if not p.is_absolute() or p.name != filename:
            raise ValueError(f'Unexpected demographic source: {name}')
        primitive.verify(p, row['sha256'])
    census_dirs = {str(Path(c['demographic_sources'][name]['path']).parent)
                   for name in files if name != 'acs_headship'}
    if len(census_dirs) != 1:
        raise ValueError('Four Census inputs must occupy the same source directory')
    receipt = json.loads(Path(c['primitive_summary']).read_text())
    old = json.loads(Path(c['primitive_contract']).read_text())
    if old.get('schema') != primitive.SCHEMA:
        raise RuntimeError('Unexpected originating primitive contract schema')
    if receipt.get('status') != 'passed_primitive_smoke_only' or receipt.get('arm') != arm:
        raise RuntimeError('Originating primitive smoke must have passed for this arm')
    for name in ('checkpoint_sha256', 'selected_summary_sha256'):
        if old[name] != c[name]:
            raise RuntimeError(f'Primitive/current input differs: {name}')
    if old['arm'] != arm or old['original_wealth_points'] != 120 or len(old['wealth_grid']) != 120:
        raise RuntimeError('Originating primitive must use this arm and full 120-node grid')
    parent = Path(c['primitive_summary']).parent
    if (Path(c['primitive_contract']) != parent / 'contract.json'
            or Path(c['stationary_arrays']) != parent / 'stationary_arrays.npz'):
        raise RuntimeError('Stationary seed and contract must belong to the originating primitive packet')
    expected_flags = dict(joint_nested_choice=arm == 'nested', fertility_nest_choice=arm == 'nested',
                          two_shock_choice=False, exhaustive_saving_control=True)
    if old['flags'] != expected_flags:
        raise RuntimeError('Originating primitive has different choice/saving conventions')
    for source in (root / 'code/model/intergen_eqscale_seq_optimized').glob('*.py'):
        relative = str(source.relative_to(root))
        if old['source_sha256'].get(relative) != c['source_sha256'][relative]:
            raise RuntimeError(f'Stationary seed solver source differs: {relative}')
    return c, old


def check_smoke_gates(result, expected_years=(2007, 2011, 2015, 2019, 2023, 2027)):
    """Accounting gates only: prescribed prices do not imply market clearing."""
    gates = dict(historical_mass=(result.history.maximum_mass_accounting_error, 2e-8),
        historical_reproduction=(result.history.maximum_policy_reproduction_error, 2e-10),
        historical_projection=(result.history.maximum_feasibility_projection_mass, 1e-6),
        person_reproduction=(result.person_tail.maximum_policy_reproduction_error, 2e-10),
        person_identity=(result.person_tail.maximum_person_identity_error, 2e-9),
        head_identity=(result.person_tail.maximum_head_identity_error, 2e-9),
        household_head_identity=(result.person_tail.maximum_household_person_head_gap, 2e-9),
        age_head_identity=(result.person_tail.maximum_age_head_gap, 2e-9),
        initial_2023_age_identity=(result.initial_2023_age_head_gap, 2e-9),
        person_projection=(result.person_tail.maximum_feasibility_projection_mass, 1e-6))
    for row in result.person_tail.rows:
        gates['raw_household_mass_' + str(row['calendar_year'])] = (abs(row['raw_household_mass_residual']), 2e-8)
    for name, (actual, tolerance) in gates.items():
        if not np.isfinite(actual) or abs(actual) > tolerance:
            raise RuntimeError(f'Joined conditional smoke gate failed: {name}={actual}; limit={tolerance}')
    if (result.bellman_solves != 2 * len(expected_years)
            or [r['calendar_year'] for r in result.rows] != list(expected_years)):
        raise RuntimeError('Expected all supplied dates and exactly two Bellman calls per date')
    return {name: dict(value=float(value), tolerance=tolerance, passed=True)
            for name, (value, tolerance) in gates.items()}


def run_smoke(args):
    import copy
    import gzip
    import json
    import os
    from pathlib import Path
    import pickle
    import threading
    import time
    import run_e5f_matched_pf_smoke as primitive
    c, originating = load_smoke_contract(args.contract, args.contract_sha256, args.arm)
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'Refusing nonempty output: {out}')
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    stopped = threading.Event()
    progress = dict(arm=args.arm, phase='prepare_inputs', completed_dates=0)
    def save(name, value):
        pf.write_json(out / name, value)
    def heartbeat():
        while not stopped.wait(30):
            elapsed = time.monotonic() - started
            save('heartbeat.json', dict(progress, elapsed_seconds=elapsed))
            if elapsed > int(c['seconds']):
                save('failure.json', dict(progress, error='contracted wall-time cap', elapsed_seconds=elapsed))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    try:
        pf.transition.configure_sequential_model()
        pf.calendar.apply_fertility = pf.transition.apply_sequential_fertility
        pf.calendar.advance_calendar_distribution = pf.transition.advance_sequential_calendar_distribution
        opener = gzip.open if Path(c['checkpoint']).suffix == '.gz' else open
        with opener(c['checkpoint'], 'rb') as stream:
            packet = pickle.load(stream)
        P = copy.deepcopy(packet['parameters'])
        selected = json.loads(Path(c['selected_summary']).read_text())
        selected_match = primitive.verify_selected_parameters(P, selected, packet['supply_rule'])
        for name, flag in originating['flags'].items():
            setattr(P, name, flag)
        original_transfer = float(getattr(P, 'property_tax_lump_sum_transfer', 0.))
        P.property_tax_lump_sum_transfer = 0.
        grid = np.asarray(packet['b_grid'], dtype=float).copy()
        with np.load(c['stationary_arrays'], allow_pickle=False) as archive:
            terminal_V = archive['V'].copy()
            original_pre = archive['g_pre'].copy()
            seed_grid = archive['wealth_grid'].copy()
        if (len(grid) != 120 or not np.array_equal(grid, seed_grid)
                or not np.array_equal(grid, np.asarray(originating['wealth_grid']))
                or terminal_V.shape != original_pre.shape
                or original_pre.shape != packet['evaluation'].g_pre.shape
                or not np.isfinite(terminal_V).all() or not np.isfinite(original_pre).all()
                or original_pre.min() < 0):
            raise RuntimeError('Full-grid stationary seed shape/grid/health mismatch')
        P.Nb = len(grid)
        initial_mass = float(original_pre.sum())
        ages = P.age_start + np.arange(P.J) * P.da
        initial_g, initial_audit = pf.transition.reweight_distribution_to_observed_age_path(
            original_pre, ages, year=2007, initial_mass=initial_mass)
        template_2023, template_audit = pf.transition.reweight_distribution_to_observed_age_path(
            original_pre, ages, year=2023, initial_mass=initial_mass)
        sources = c['demographic_sources']
        demographics = person_pf.build_annual_demographic_primitives(
            template_2023, P, source_dir=Path(sources['population_mid']['path']).parent,
            headship_dir=Path(sources['acs_headship']['path']).parent, start_year=2023)
        renewal = selected['renewal_accounting_old_state']
        renewal_contract = selected['renewal_accounting_contract']
        slots = int(renewal_contract['birth_vintage_queue_waiting_slots'])
        conversion = float(renewal_contract['effective_birth_to_household_conversion'])
        if slots != 4 or abs(conversion - 1 / 2.1) > 1e-15:
            raise RuntimeError('Unexpected retained historical queue specification')
        initial = pf.PFInitialState(initial_g,
            [float(renewal['old_queue_mature_flow_B'])] * slots,
            [float(renewal['old_raw_birth_queue_flow_B'])] * slots)
        conditioning = pf.HistoricalConditioning(2007, initial_mass,
            {1: 2011, 2: 2015, 3: 2019, 4: 2023},
            float(renewal['outside_flow_M']), float(renewal['retention_rho']))
        price = float(np.asarray(packet['evaluation'].policy.price).reshape(-1)[0])
        years = [2007, 2011, 2015, 2019, 2023, 2027]
        prices = np.full(6, price)
        psi = np.full(6, float(P.psi_child))
        rents = pf.rents_from_asset_prices(prices, price, P)
        save('contract.json', dict(c, contract_sha256=args.contract_sha256, arm=args.arm,
            selected_parameter_match=selected_match, flags=originating['flags'],
            original_transfer=original_transfer, supplied_transfers=[0.] * 6,
            years=years, supplied_prices=prices, supplied_psi=psi,
            initial_seed='synthetic fixed-price stationary distribution at current selected psi; not old normalized 2.1 history',
            terminal_boundary='matching fixed-price stationary value; finite-boundary plumbing only',
            historical_law='retained selected outside flow, retention and birth queues; observed age bridges',
            supply='unchanged checkpoint rule and normalization; no reanchoring',
            scope='conditional joined history/person smoke; no GE, calibration fit, policy effect or terminal convergence claim',
            expected_bellman_solves=12, seconds=c['seconds']))
        save('initialization.json', dict(initial_mass=initial_mass, initial_age_bridge=initial_audit,
            template_2023_age_bridge=template_audit, renewal=renewal,
            headship_alignment_factors=demographics.model_age_headship_alignment_factors,
            model_person_scale=demographics.scale_model_units_per_person))
        primitive.save_arrays(out / 'initial_states.npz', g_2007=initial_g,
            g_2023_template=template_2023, persons_2023=demographics.initial_person_state.persons,
            heads_2023=demographics.initial_person_state.heads, terminal_V=terminal_V, wealth_grid=grid)
        observations = []
        def observer(index, evaluation, parameters, b_grid, shared):
            if float(evaluation.feasibility_projection_mass) > 1e-6:
                raise RuntimeError(f'Feasibility projection failed at {years[index]}')
            budget = primitive.dated_budget(evaluation, parameters, shared, b_grid, float(rents[index]))
            row = dict(period=index, calendar_year=years[index], births=float(evaluation.births),
                household_mass=float(evaluation.g_current.sum()), inherited_mass=float(evaluation.g_pre.sum()),
                market_residual_diagnostic_only=float(evaluation.relative_market_residual),
                budget=budget, elapsed_seconds=time.monotonic() - started)
            observations.append(row)
            progress.update(phase='forward', completed_dates=len(observations), current_year=years[index])
            save('latest_completed.json', row)
            save('observed_dates.json', observations)
            save('heartbeat.json', dict(progress, elapsed_seconds=time.monotonic() - started))
        progress['phase'] = 'backward_and_forward'
        save('heartbeat.json', dict(progress, elapsed_seconds=time.monotonic() - started))
        result = evaluate_history_and_person_tail(years=years, prices=prices, psi_path=psi,
            transfer_path=np.zeros(6), terminal_price=price, terminal_V=terminal_V,
            base_parameters=P, b_grid=grid, initial_state=initial, historical_conditioning=conditioning,
            initial_2023_persons=demographics.initial_person_state,
            demographic_primitives=demographics, supply_rule=packet['supply_rule'],
            birth_to_entry_conversion=conversion, observer=observer)
        final = result.person_tail.terminal_state
        primitive.save_arrays(out / 'final_state.npz', g_pre=final.g_pre,
            persons=final.persons.persons, heads=final.persons.heads,
            year=np.array(final.persons.year), wealth_grid=grid)
        pf.write_csv(out / 'transition_path.csv', result.rows)
        gates = check_smoke_gates(result)
        save('summary.json', dict(status='passed_conditional_history_smoke_only', arm=args.arm,
            elapsed_seconds=time.monotonic() - started, bellman_solves=result.bellman_solves,
            contract_sha256=args.contract_sha256, gates=gates, completed_dates=len(observations),
            market_residual_reported_not_gated=max(abs(r['relative_market_residual']) for r in result.rows),
            final_state_year=final.persons.year, final_household_mass=float(final.g_pre.sum()),
            final_persons=float(final.persons.persons.sum()), equilibrium_verified=False,
            calibrated_history=False, terminal_convergence_tested=False,
            artifact_sha256={name: primitive.digest(out / name) for name in
                ('final_state.npz', 'initial_states.npz', 'transition_path.csv', 'observed_dates.json', 'initialization.json')}))
    except Exception as error:
        save('failure.json', dict(progress, error=repr(error), elapsed_seconds=time.monotonic() - started))
        raise
    finally:
        stopped.set()


def main():
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--arm', choices=('sequential', 'nested'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    run_smoke(parser.parse_args())


if __name__ == '__main__':
    main()
