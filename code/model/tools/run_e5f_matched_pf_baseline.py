"""Bounded matched-PF normalization, endpoint and historical diagnostics.

Stages use explicit inherited source, target, fiscal and demographic contracts.
A conditional path receipt does not imply market or horizon convergence.
"""
from __future__ import annotations

import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import copy
import gzip
import json
from pathlib import Path
import pickle
import threading
import time
from types import SimpleNamespace

import numpy as np
import run_e5f_matched_pf_history as joined
import run_e5f_matched_pf_smoke as primitive

pf = joined.pf
person = joined.person_pf


def load_seed(contract, originating, arm):
    """Load a reproduced stationary policy and frozen empirical demographics."""
    with gzip.open(contract['checkpoint'], 'rb') as stream:
        packet = pickle.load(stream)
    P = copy.deepcopy(packet['parameters'])
    selected = json.loads(Path(contract['selected_summary']).read_text())
    primitive.verify_selected_parameters(P, selected, packet['supply_rule'])
    for name, value in originating['flags'].items():
        setattr(P, name, value)
    if (P.joint_nested_choice != (arm == 'nested') or float(P.tau_H) != .04
            or float(P.period_years) != 4.):
        raise RuntimeError('Unexpected architecture or retained historical tax')
    P.property_tax_lump_sum_transfer = 0.
    pf.transition.configure_sequential_model()
    pf.calendar.apply_fertility = pf.transition.apply_sequential_fertility
    pf.calendar.advance_calendar_distribution = pf.transition.advance_sequential_calendar_distribution
    pf.calendar.distribution_rows = pf.transition.independent_child_distribution_rows
    grid = np.asarray(packet['b_grid']).copy()
    with np.load(contract['stationary_arrays'], allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    if len(grid) != 120 or not np.array_equal(grid, arrays['wealth_grid']):
        raise RuntimeError('Stationary policy wealth grid mismatch')
    price = float(packet['evaluation'].policy.price[0])
    shared = primitive.model.precompute_shared(P, grid)
    joint = (SimpleNamespace(**{name: arrays['joint_' + name] for name in primitive.JOINT_FIELDS})
             if P.joint_nested_choice else None)
    policy = pf.calendar.PolicyBundle(
        **{name: arrays[name] for name in primitive.FIELDS}, price=np.array([price]),
        maps=pf.calendar.build_transition_maps(np.array([price]), P, grid, shared),
        joint_choice=joint)
    ages = P.age_start + np.arange(P.J) * P.da
    template, audit = pf.transition.reweight_distribution_to_observed_age_path(
        arrays['g_pre'], ages, year=2023, initial_mass=float(arrays['g_pre'].sum()))
    sources = contract['demographic_sources']
    demographics = person.build_annual_demographic_primitives(
        template, P, source_dir=Path(sources['population_mid']['path']).parent,
        headship_dir=Path(sources['acs_headship']['path']).parent, start_year=2023)
    return SimpleNamespace(parameters=P, grid=grid, policy=policy, pre=arrays['g_pre'],
        demographics=demographics, supply=packet['supply_rule'], price=price,
        selected=selected, shared=shared, age_audit=audit)


def load_normalized(contract, arm):
    """Require a completed, pinned normalized-old initialization for this arm."""
    for key in ('normalized_checkpoint', 'normalized_summary', 'normalized_contract'):
        primitive.verify(contract[key], contract[key + '_sha256'])
    receipt = json.loads(Path(contract['normalized_summary']).read_text())
    parent = json.loads(Path(contract['normalized_contract']).read_text())
    if (receipt.get('status') != 'passed_normalized_old_initialization'
            or receipt.get('arm') != arm or parent['arm'] != arm
            or receipt['checkpoint_sha256'] != contract['normalized_checkpoint_sha256']
            or parent['checkpoint_sha256'] != contract['checkpoint_sha256']
            or parent['selected_summary_sha256'] != contract['selected_summary_sha256']):
        raise RuntimeError('Normalized-old source/input/arm receipt mismatch')
    with gzip.open(contract['normalized_checkpoint'], 'rb') as stream:
        packet = pickle.load(stream)
    if packet['contract_sha256'] != parent['contract_sha256']:
        raise RuntimeError('Normalized checkpoint parent-contract identity mismatch')
    old = packet['old']
    if (bool(old.parameters.joint_nested_choice) != (arm == 'nested')
            or not np.array_equal(old.years, [2007, 2011, 2015, 2019, 2023])
            or float(old.supply_rule.elasticity) != .63):
        raise RuntimeError('Normalized-old architecture or historical contract changed')
    return old, packet['demographics']


def run_terminal_root(seed, c, args, out, progress, save, started):
    import e5f_matched_pf_endpoint as endpoint
    import e5f_matched_pf_price_root as price_root
    if c.get('terminal_preference_rule') != 'hold_normalized_2023_intercept':
        raise ValueError('Explicit diagnostic terminal-preference rule required')
    old, demographics = load_normalized(c, args.arm)
    base = copy.deepcopy(old.parameters)
    base.psi_child = float(old.psi_path[-1])
    controls = endpoint.EndpointControls(250, .5, 1e-9, 1e-10, 1e-8, 2e-4, 2.5e-5, 2e-9)
    save('contract.json', dict(c, contract_sha256=args.contract_sha256, arm=args.arm,
        terminal_psi=float(base.psi_child), historical_psi=old.psi_path,
        fiscal_regime='fixed_transfer', transfer=0., period_tax=float(base.tau_H),
        supply_rule=vars(old.supply_rule), endpoint_controls=vars(controls),
        terminal_preference_status='diagnostic no-further-shock continuation, not an estimated future path',
        demographic_tail='existing frozen 2100 inputs',
        scope='terminal stationary price/population root; historical path still requires a separate solve'))
    progress.update(phase='terminal_root', root_evaluations=0)
    def evaluate(price):
        P = copy.deepcopy(base)
        grid = old.b_grid
        shared = primitive.model.precompute_shared(P, grid)
        progress.update(phase='terminal_stationary_policy', root_evaluations=progress['root_evaluations']+1,
                        trial_price=price)
        solution = primitive.model.solve_markov_income_at_prices(np.array([price]), P, grid, SD=shared)
        policy = pf.calendar.policy_from_solution(solution, np.array([price]), P, grid, shared)
        pre, reconstruction = pf.calendar.reconstruct_stationary_pre_fertility(solution, policy, P, grid, shared)
        progress['phase'] = 'terminal_population_mapping'
        result = endpoint.evaluate_endpoint(parameters=P, b_grid=grid, policy=policy,
            asset_price=price, transfer=0., psi_child=float(base.psi_child),
            demographic_primitives=demographics, initial_g_pre=pre,
            supply_rule=old.supply_rule, fiscal_regime='fixed_transfer', controls=controls)
        # Retain only the policy belonging to this endpoint; the root retains
        # the current best endpoint rather than a cache of all large tensors.
        result.policy = policy
        result.parameters = P
        result.reconstruction = reconstruction
        check = pf.calendar.evaluate_period(np.array([price]), result.fixed_point.g_pre,
            P, grid, shared, pf.calendar.SolveCounter(), supply_rule=old.supply_rule,
            supplied_policy=policy)
        primitive.dated_budget(check, P, shared, grid, float(P.user_cost_rate)*price)
        if float(check.feasibility_projection_mass) > 1e-6:
            raise RuntimeError('Terminal root feasibility gate failed')
        return result
    def record(latest, best):
        save('latest_completed.json', latest)
        if best is not None:
            save('best_so_far.json', best)
        progress.update(phase='terminal_root', root_evaluations=latest['evaluation'])
    root = price_root.solve_price_root(fresh_evaluate=evaluate, start_price=float(seed.price),
        bound_ratios=(.6, 1.5), maximum_evaluations=20,
        deadline_monotonic=started+c['seconds'], market_tolerance=2e-4,
        replay_tolerance=2e-10, progress_callback=record)
    pf.write_csv(out / 'root_history.csv', root.records)
    selected = root.final_endpoint if root.converged else root.best_endpoint
    scalars = dict(status=root.status, arm=args.arm, converged=root.converged,
        evaluations=root.evaluations, cache_hits=root.cache_hits, best_record=root.best_record,
        replay_residual_gap=root.replay_residual_gap, elapsed_seconds=time.monotonic()-started,
        historical_equilibrium_verified=False, calibrated_history=False)
    if selected is not None:
        scalars.update(endpoint_residuals=selected.residuals, endpoint_gates=selected.gates,
                       endpoint_contract=selected.contract)
        checkpoint = out / 'terminal.pkl.gz'
        with gzip.open(checkpoint, 'wb', compresslevel=1) as stream:
            pickle.dump(dict(parameters=selected.parameters, b_grid=old.b_grid,
                policy=selected.policy, fixed_point=selected.fixed_point,
                demographics=demographics, supply_rule=old.supply_rule,
                contract_sha256=args.contract_sha256), stream, protocol=pickle.HIGHEST_PROTOCOL)
        with gzip.open(checkpoint, 'rb') as stream:
            replay = pickle.load(stream)
        if (not np.array_equal(replay['policy'].V, selected.policy.V)
                or not np.array_equal(replay['fixed_point'].g_pre, selected.fixed_point.g_pre)):
            raise RuntimeError('Terminal checkpoint reload failed')
        scalars['checkpoint_sha256'] = primitive.digest(checkpoint)
    save('summary.json', scalars)


def load_terminal(c, arm, old):
    """Use only a reproduced endpoint rooted from this exact old normalization."""
    for key in ('terminal_checkpoint', 'terminal_summary', 'terminal_contract'):
        primitive.verify(c[key], c[key + '_sha256'])
    receipt = json.loads(Path(c['terminal_summary']).read_text())
    parent = json.loads(Path(c['terminal_contract']).read_text())
    if (receipt.get('status') != 'complete_reproduced_root' or not receipt.get('converged')
            or receipt['arm'] != arm or parent['arm'] != arm
            or not all(receipt['endpoint_gates'].values())
            or receipt['replay_residual_gap'] > 2e-10
            or receipt['checkpoint_sha256'] != c['terminal_checkpoint_sha256']
            or parent['normalized_checkpoint_sha256'] != c['normalized_checkpoint_sha256']
            or parent['terminal_preference_rule'] != 'hold_normalized_2023_intercept'):
        raise RuntimeError('Terminal receipt/arm/normalization mismatch')
    with gzip.open(c['terminal_checkpoint'], 'rb') as stream:
        packet = pickle.load(stream)
    P = packet['parameters']
    if (packet['contract_sha256'] != parent['contract_sha256']
            or not np.array_equal(packet['b_grid'], old.b_grid)
            or vars(packet['supply_rule']) != vars(old.supply_rule)
            or float(P.psi_child) != float(old.psi_path[-1])
            or bool(P.joint_nested_choice) != (arm == 'nested')
            or float(P.tau_H) != .04 or float(P.property_tax_lump_sum_transfer) != 0.):
        raise RuntimeError('Terminal checkpoint contract mismatch')
    return packet


def run_history_probe(seed, c, args, out, progress, save, started, *, prices_override=None):
    """One complete conditional path, also used for independent price probes."""
    import e5f_matched_pf_moments as moments
    import run_e5f_perfect_foresight_person_demography_policy as terminal_checks
    import run_e5f_perfect_foresight_rebated_property_tax as rent_domain
    old, demographics = load_normalized(c, args.arm)
    terminal = load_terminal(c, args.arm, old)
    count = c['path_date_count']
    coordinate = c['probe_coordinate']
    if (type(count) is not int or not 6 <= count <= 40
            or type(coordinate) is not int or not -1 <= coordinate < count
            or c['probe_log_step'] != .01
            or c['terminal_preference_rule'] != 'hold_normalized_2023_intercept'
            or c['initial_price_rule'] != 'log_old_to_selected_2023_then_terminal'):
        raise ValueError('Explicit bounded historical price-probe contract required')
    years = 2007 + 4 * np.arange(count)
    psi = np.r_[old.psi_path, np.full(count - 5, old.psi_path[-1])]
    terminal_price = float(terminal['policy'].price[0])
    reference = SimpleNamespace(parameters=terminal['parameters'], asset_price=terminal_price,
        renter_price=float(terminal['parameters'].user_cost_rate)*terminal_price,
        equal_transfer=0., psi_child=float(psi[-1]),
        state=person.PersonPFState(g_pre=terminal['fixed_point'].g_pre,
                                  persons=terminal['fixed_point'].persons))
    # Numerical guesses only: every dated price is free in the later root.
    raw_prices = np.exp(np.r_[np.linspace(np.log(old.supply_rule.initial_price), np.log(seed.price), 5),
        np.linspace(np.log(seed.price), np.log(terminal_price), count-4)[1:]])
    anchor, anchor_projection = rent_domain.project_price_path_to_positive_rents(
        raw_prices, terminal=reference, minimum_rent_share=1e-6)
    prices = anchor.copy() if prices_override is None else np.asarray(prices_override, dtype=float).copy()
    if prices.shape != anchor.shape or not np.isfinite(prices).all() or np.any(prices <= 0):
        raise ValueError('Explicit root prices must be a positive finite vector matching the complete horizon')
    if prices_override is not None and coordinate != -1:
        raise ValueError('A root path cannot also be a coordinate probe')
    if coordinate >= 0:
        prices[coordinate] *= np.exp(c['probe_log_step'])
    prices, projection = rent_domain.project_price_path_to_positive_rents(
        prices, terminal=reference, minimum_rent_share=1e-6)
    if coordinate >= 0 and projection['adjusted_period_count']:
        raise RuntimeError('Coordinate probe would move other dates through rent projection')
    rents = pf.rents_from_asset_prices(prices, terminal_price, old.parameters)
    measurement = moments.measurement
    domain = tuple((r['name'], r['lower'], r['upper'], r['transform'])
                   for r in seed.selected['panel_design']['domain'])
    measurement.TRANSITION_SEARCH_DOMAIN = domain
    targets = measurement.e5_target_system_for_profile('baseline')
    if targets.fingerprint != c['target_fingerprint']:
        raise RuntimeError('Historical price probe target fingerprint changed')
    old_evaluation = pf.calendar.evaluate_period(old.policy.price, old.stationary_g_pre,
        old.parameters, old.b_grid, old.shared, pf.calendar.SolveCounter(),
        supply_rule=old.supply_rule, supplied_policy=old.policy)
    moment_observer = moments.HistoricalMomentObserver(target_system=targets,
        expected_target_fingerprint=c['target_fingerprint'], old_parameters=old.parameters,
        old_first_birth_accounting=measurement.first_birth_accounting_by_age(old_evaluation, old.parameters),
        old_normalization=old.diagnostics['normalization'], old_normalization_tolerance=5e-4)
    save('contract.json', dict(c, contract_sha256=args.contract_sha256, arm=args.arm,
        years=years, prices=prices, anchor_prices=anchor, psi_path=psi,
        transfer_path=np.zeros(count), terminal_price=terminal_price,
        supply_rule=vars(old.supply_rule), anchor_projection=anchor_projection, projection=projection,
        expected_bellman_solves=2*count,
        scope='conditional finite-horizon PF price evaluation; not market equilibrium or recalibration',
        price_origin='numerical_anchor' if prices_override is None else 'explicit_root_trial',
        outstanding='post-2023 preference continuation diagnostic; horizon convergence and dated markets required'))
    observations = []
    def observe(index, evaluation, parameters, grid, shared):
        budget = primitive.dated_budget(evaluation, parameters, shared, grid, float(rents[index]))
        if float(evaluation.feasibility_projection_mass) > 1e-6:
            raise RuntimeError('Historical price probe feasibility gate failed')
        moment_observer(index, evaluation, parameters, grid, shared)
        row = dict(period=index, calendar_year=int(years[index]), budget=budget,
            market_residual=float((evaluation.demand_by_loc[0] - evaluation.supply_by_loc[0]) / evaluation.supply_by_loc[0]),
            elapsed_seconds=time.monotonic()-started)
        observations.append(row)
        progress.update(phase='historical_forward', completed_dates=len(observations), current_year=int(years[index]))
        save('latest_completed.json', row)
    progress['phase'] = 'historical_backward_and_forward'
    result = joined.evaluate_history_and_person_tail(years=years, prices=prices, psi_path=psi,
        transfer_path=np.zeros(count), terminal_price=terminal_price, terminal_V=terminal['policy'].V,
        base_parameters=old.parameters, b_grid=old.b_grid, initial_state=old.initial_state,
        historical_conditioning=old.historical_conditioning,
        initial_2023_persons=demographics.initial_person_state,
        demographic_primitives=demographics, supply_rule=old.supply_rule,
        birth_to_entry_conversion=1./2.1, observer=observe)
    gates = joined.check_smoke_gates(result, expected_years=years.tolist())
    report = moment_observer.report('matched_pf_' + args.arm,
        theta=seed.selected['best_candidate']['theta'], parameter_domain=domain, supply_rule=old.supply_rule)
    tail_distance = terminal_checks.terminal_convergence_diagnostics(result.person_tail,
        terminal=reference, psi_path=psi[4:])
    residual = np.array([(r['housing_demand'] - r['housing_supply']) / r['housing_supply']
                         for r in result.rows])
    if not np.isfinite(residual).all():
        raise RuntimeError('Nonfinite historical price residuals')
    pf.write_csv(out/'transition_path.csv', result.rows)
    pf.write_csv(out/'target_fit.csv', report['target_fit_rows'])
    pf.write_csv(out/'parameters.csv', report['parameter_rows'])
    save('measurement.json', report)
    save('observed_dates.json', observations)
    final = result.person_tail.terminal_state
    primitive.save_arrays(out/'final_state.npz', g_pre=final.g_pre,
        persons=final.persons.persons, heads=final.persons.heads,
        wealth_grid=old.b_grid, year=np.array(final.persons.year))
    summary = dict(status='passed_conditional_historical_price_probe', arm=args.arm,
        contract_sha256=args.contract_sha256, elapsed_seconds=time.monotonic()-started,
        bellman_solves=result.bellman_solves, years=years, prices=prices, anchor_prices=anchor,
        residual=residual, maximum_market_residual=float(np.max(np.abs(residual))),
        mapping_valid=True, gates=gates, terminal_distance=tail_distance,
        target_fingerprint=targets.fingerprint, loss=report['loss'], target_count=12,
        estimated_parameter_count=11, market_equilibrium_verified=False,
        calibrated_history=False, production_promoted=False,
        artifact_sha256={name: primitive.digest(out/name) for name in
            ('target_fit.csv','parameters.csv','transition_path.csv','measurement.json','final_state.npz')})
    if c.get('save_initial_2023_state', False):
        if c['save_initial_2023_state'] is not True:
            raise ValueError('Dated-state saving flag must be explicitly Boolean')
        P2023 = copy.deepcopy(old.parameters)
        P2023.psi_child = float(old.psi_path[-1])
        initial2023 = person.PersonPFState(g_pre=result.history.terminal_state.g_pre.copy(),
                                         persons=demographics.initial_person_state)
        checkpoint = out/'initial_2023.pkl.gz'
        with gzip.open(checkpoint, 'wb', compresslevel=1) as stream:
            pickle.dump(dict(parameters=P2023, b_grid=old.b_grid, initial_state=initial2023,
                demographic_primitives=demographics, supply_rule=old.supply_rule,
                contract_sha256=args.contract_sha256,
                state_timing='before 2023 fertility/tenure choices; after observed historical age bridge',
                equilibrium_certified=False, policy_announcement_included=False),
                stream, protocol=pickle.HIGHEST_PROTOCOL)
        with gzip.open(checkpoint, 'rb') as stream:
            loaded = pickle.load(stream)
        if (not np.array_equal(loaded['initial_state'].g_pre, initial2023.g_pre)
                or not np.array_equal(loaded['initial_state'].persons.persons, initial2023.persons.persons)
                or not np.array_equal(loaded['initial_state'].persons.heads, initial2023.persons.heads)
                or loaded['initial_state'].persons.year != 2023
                or loaded['contract_sha256'] != args.contract_sha256):
            raise RuntimeError('2023 inherited-state checkpoint reload failed')
        summary['artifact_sha256']['initial_2023.pkl.gz'] = primitive.digest(checkpoint)
    save('summary.json', summary)
    save('best_so_far.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', required=True, type=Path)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--arm', required=True, choices=('sequential', 'nested'))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    primitive.verify(args.contract, args.contract_sha256)
    requested = json.loads(args.contract.read_text())
    mode = requested.get('experiment')
    caps = {'fixed_policy_person_endpoint_current_psi_no_rebate': 840,
            'normalized_old_initialization': 1680,
            'normalized_terminal_price_root': 1680,
            'normalized_historical_path_probe': 1800}
    if mode not in caps:
        raise ValueError('Explicit supported baseline diagnostic experiment required')
    c, originating = joined.load_smoke_contract(args.contract, args.contract_sha256,
                                               args.arm, maximum_seconds=caps[mode])
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    stop = threading.Event()
    progress = dict(arm=args.arm, phase='load', population_evaluations=0)
    def save(name, value):
        pf.write_json(out / name, value)
    def heartbeat():
        while not stop.wait(15):
            elapsed = time.monotonic() - started
            save('heartbeat.json', dict(progress, elapsed_seconds=elapsed))
            if elapsed > c['seconds']:
                save('failure.json', dict(progress, error='wall-time budget exhausted'))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    original_evaluate = pf.calendar.evaluate_period
    try:
        seed = load_seed(c, originating, args.arm)
        P = seed.parameters
        if mode == 'normalized_historical_path_probe':
            run_history_probe(seed, c, args, out, progress, save, started)
            return
        if mode == 'normalized_terminal_price_root':
            run_terminal_root(seed, c, args, out, progress, save, started)
            return
        if mode == 'normalized_old_initialization':
            import e5f_matched_pf_initial_state as initialization
            if c.get('old_fertility_tolerance') != 5e-4 or c.get('maximum_stationary_solves') != 8:
                raise ValueError('Expected retained old-normalization tolerance and bounded eight-solve cap')
            save('contract.json', dict(c, contract_sha256=args.contract_sha256, arm=args.arm,
                scope='restore normalized pre-announcement old state; no anticipated historical equilibrium'))
            def normalized_progress(record):
                progress.update(phase='old_normalization', stationary_solve=record['index'],
                                stationary_status=record['status'])
                save('latest_stationary.json', dict(record, elapsed_total=time.monotonic()-started))
                if record['status'] == 'complete':
                    save('latest_completed.json', dict(record, elapsed_total=time.monotonic()-started))
            old = initialization.initialize_normalized_old_state(
                parameters=P, b_grid=seed.grid, selected_summary=seed.selected,
                selected_supply_rule=seed.supply, arm=args.arm,
                fiscal_contract=initialization.CalibrationFiscalContract(.01, .04, 0., 'retained_calibration_unrebated'),
                completed_fertility_tolerance=c['old_fertility_tolerance'],
                max_stationary_solves=c['maximum_stationary_solves'],
                deadline_monotonic=started+c['seconds'], progress=normalized_progress)
            template, _ = pf.transition.reweight_distribution_to_observed_age_path(
                old.stationary_g_pre, old.parameters.age_start + np.arange(old.parameters.J)*old.parameters.da,
                year=2023, initial_mass=float(old.initial_state.g_pre.sum()))
            sources = c['demographic_sources']
            demographics = person.build_annual_demographic_primitives(template, old.parameters,
                source_dir=Path(sources['population_mid']['path']).parent,
                headship_dir=Path(sources['acs_headship']['path']).parent, start_year=2023)
            checkpoint = out / 'normalized_old.pkl.gz'
            with gzip.open(checkpoint, 'wb', compresslevel=1) as stream:
                pickle.dump(dict(old=old, demographics=demographics,
                    contract_sha256=args.contract_sha256), stream, protocol=pickle.HIGHEST_PROTOCOL)
            with gzip.open(checkpoint, 'rb') as stream:
                replay = pickle.load(stream)
            if (not np.array_equal(replay['old'].initial_state.g_pre, old.initial_state.g_pre)
                    or not np.array_equal(replay['old'].policy.V, old.policy.V)
                    or replay['contract_sha256'] != args.contract_sha256):
                raise RuntimeError('Normalized-old checkpoint reload failed')
            save('summary.json', dict(status='passed_normalized_old_initialization', arm=args.arm,
                elapsed_seconds=time.monotonic()-started, diagnostics=old.diagnostics,
                historical_years=old.years, historical_psi=old.psi_path,
                checkpoint_sha256=primitive.digest(checkpoint), historical_equilibrium_verified=False,
                calibrated_history=False))
            return
        save('contract.json', dict(c, contract_sha256=args.contract_sha256,
            arm=args.arm, fixed_price=seed.price, fixed_psi=float(P.psi_child),
            fiscal_regime='retained 1% annual tax, zero transfer; fiscal surplus reported',
            demographic_tail='existing 2100 survival, migration and headship held fixed',
            supply_rule=vars(seed.supply), maximum_inner_iterations=250,
            distribution_tolerance=1e-9, birth_rate_tolerance=1e-10, one_step_tolerance=1e-8,
            scope='stationary population at fixed household policy; no price root or calibration'))
        def observed(*pos, **kw):
            evaluation = original_evaluate(*pos, **kw)
            progress.update(phase='population_fixed_point',
                population_evaluations=progress['population_evaluations'] + 1)
            primitive.dated_budget(evaluation, P, seed.shared, seed.grid,
                                   float(P.user_cost_rate) * seed.price)
            if float(evaluation.feasibility_projection_mass) > 1e-6:
                raise RuntimeError('Endpoint feasibility projection exceeds unchanged diagnostic gate')
            save('latest_completed.json', dict(progress, elapsed_seconds=time.monotonic()-started,
                household_mass=float(evaluation.g_current.sum()), births=float(evaluation.births),
                market_residual_diagnostic_only=float(evaluation.relative_market_residual)))
            return evaluation
        pf.calendar.evaluate_period = observed
        result = person.solve_terminal_household_person_fixed_point(
            policy=seed.policy, parameters=P, b_grid=seed.grid, initial_g_pre=seed.pre,
            demographic_primitives=seed.demographics, supply_rule=seed.supply,
            maximum_iterations=250, damping=.5, distribution_tolerance=1e-9,
            birth_rate_tolerance=1e-10, one_step_tolerance=1e-8)
        pf.write_csv(out / 'inner_history.csv', result.history)
        scalars = {key: value for key, value in vars(result).items()
                   if isinstance(value, (int, float, bool, str, np.number))}
        primitive.save_arrays(out / 'endpoint_state.npz', g_pre=result.g_pre,
            persons=result.persons.persons, heads=result.persons.heads,
            wealth_grid=seed.grid)
        save('summary.json', dict(scalars, **progress, elapsed_seconds=time.monotonic()-started,
            status='passed_fixed_policy_population_endpoint' if result.converged else 'population_endpoint_not_converged',
            market_equilibrium_verified=False, calibrated_history=False,
            endpoint_state_sha256=primitive.digest(out / 'endpoint_state.npz')))
    except Exception as error:
        save('failure.json', dict(progress, error=repr(error), elapsed_seconds=time.monotonic()-started))
        raise
    finally:
        pf.calendar.evaluate_period = original_evaluate
        stop.set()


if __name__ == '__main__':
    main()
