"""Explicit-input historical PF housing/PAYGO root; no empirical objective.

The caller owns source/receipt hashes, a process watchdog, and persistence.
The announced 2007 history and person tail share one anticipated pension path.
No demographic input, housing supply curve, target or pension guess is fitted.
"""
from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
from collections.abc import Mapping
from types import SimpleNamespace
import time

import numpy as np

from e5f_social_security import fiscal_accounts
from e5f_social_security_root import solve_social_security_path
from e5f_balanced_terminal import TerminalAuditControls, _household_checks


@dataclass
class BalancedHistoryResult:
    path: object | None
    root_receipt: dict


def _runtime():
    import run_e5f_matched_pf_history as joined
    import run_e5f_matched_pf_smoke as primitive
    import run_e5f_perfect_foresight_person_demography_policy as terminal_checks
    import run_e5f_perfect_foresight_rebated_property_tax as rent_domain
    primitive.transition.configure_sequential_model()
    primitive.calendar.model = primitive.model
    primitive.calendar.apply_fertility = primitive.transition.apply_sequential_fertility
    primitive.calendar.advance_calendar_distribution = primitive.transition.advance_sequential_calendar_distribution
    primitive.calendar.distribution_rows = primitive.transition.independent_child_distribution_rows
    return joined, primitive, terminal_checks, rent_domain


def _same(left, right):
    """Compare frozen demographic primitives without copying their arrays."""
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        return np.array_equal(left, right)
    if is_dataclass(left):
        return type(left) is type(right) and all(_same(getattr(left, f.name), getattr(right, f.name)) for f in fields(left))
    if isinstance(left, Mapping):
        return isinstance(right, Mapping) and left.keys() == right.keys() and all(_same(v, right[k]) for k, v in left.items())
    if isinstance(left, SimpleNamespace):
        return isinstance(right, SimpleNamespace) and _same(vars(left), vars(right))
    return left == right


def _validate(old, terminal, receipt, demographics, terminal_demographics, count):
    from e5f_parenthood_utility import validate_parenthood_utility, PARENTHOOD_SEARCH_NAMES
    if type(count) is not int or not 6 <= count <= 100:
        raise ValueError('A complete history requires 6 to 100 four-year dates')
    d = old.diagnostics
    norm = d['normalization']
    if (d.get('schema') != 'e5f_approved_parenthood_initial_state_v1'
            or d.get('arm') != 'sequential' or norm['target'] != 2.1
            or not np.isfinite([norm['completed_fertility'], d['verified_solution_fertility'], norm['psi_child']]).all()
            or abs(norm['completed_fertility'] - 2.1) > 5e-4
            or abs(d['verified_solution_fertility'] - 2.1) > 5e-4
            or not np.isclose(norm['completed_fertility'], d['verified_solution_fertility'], rtol=0, atol=1e-12)
            or norm['psi_child'] != float(old.parameters.psi_child)
            or not d['stationary_pension']['marginal_gate'] or not d['stationary_pension']['fiscal_gate']
            or abs(d['birth_to_entry_conversion'] - 1/2.1) > 1e-15):
        raise ValueError('Verified approved 2.1 initial normalization and balanced pension required')
    if (receipt.get('schema') != 'e5f_balanced_terminal_v1'
            or not receipt.get('converged') or not receipt.get('endpoint_production_eligible')
            or not receipt.get('fresh_endpoint_matches_final') or not receipt.get('gates')
            or not all(receipt['gates'].values()) or not terminal.mapping_valid):
        raise ValueError('A converged, independently replayed balanced terminal endpoint is required')
    P, Q = old.parameters, terminal.parameters
    for parameters in (P, Q):
        validate_parenthood_utility(parameters)
        if (parameters.I != 1 or float(parameters.period_years) != 4.
                or not parameters.scale_flows_to_period or float(parameters.tau_H) != .04
                or float(parameters.tau_pay) != .179 or float(parameters.property_tax_lump_sum_transfer) != 0.
                or not parameters.exhaustive_saving_control
                or any(bool(getattr(parameters, flag, False)) for flag in
                       ('joint_nested_choice', 'fertility_nest_choice', 'two_shock_choice'))):
            raise ValueError('Approved exhaustive sequential, fixed payroll .179 and unrebated tax .04 required')
    names = [('beta' if n == 'beta_annual' else 'hbar_first_child_jump' if n == 'h_P' else n)
             for n in PARENTHOOD_SEARCH_NAMES]
    names += ['J', 'J_R', 'age_start', 'da', 'q', 'delta', 'user_cost_rate', 'R_gross',
              'w_hat', 'income_age_profile', 'z_grid', 'retirement_income_z_scale', 'r_bar', 'xi_supply']
    if any(not _same(getattr(P, name), getattr(Q, name)) for name in names):
        raise ValueError('Historical and terminal structural parameters differ')
    for name in ('phi', 'H_own', 'Pi_z', 'Pi_child', 'survival_probs', 'use_age_survival',
                 'entry_shares', 'z_transition', 'first_birth_fixed_cost', 'kappa_fert_continuation'):
        if (hasattr(P, name) or hasattr(Q, name)) and not _same(getattr(P, name, None), getattr(Q, name, None)):
            raise ValueError(f'Historical and terminal household primitive differs: {name}')
    if not np.isfinite(Q.pension) or Q.pension <= 0:
        raise ValueError('Balanced terminal pension must be positive and finite')
    grid = np.asarray(old.b_grid)
    if (grid.ndim != 1 or grid.size < 2 or not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0)
            or not np.array_equal(grid, terminal.b_grid) or int(P.Nb) != grid.size):
        raise ValueError('Historical and terminal wealth grids must coincide')
    psi = np.asarray(old.psi_path, dtype=float)
    if (not np.array_equal(old.years, [2007, 2011, 2015, 2019, 2023]) or psi.shape != (5,)
            or not np.isfinite(psi).all() or psi[0] != float(P.psi_child)
            or float(Q.psi_child) != psi[-1]):
        raise ValueError('Terminal psi must equal the explicit 2023 historical intercept')
    if not np.allclose(psi, psi[0] + float(d['preference_change_2023'])*np.arange(5)/4,
                       rtol=0, atol=1e-14):
        raise ValueError('Preserve the approved announced linear historical preference path')
    rule = old.supply_rule
    if rule.mode != 'static-elastic' or rule.elasticity != .63:
        raise ValueError('Inherited initial supply curve must use eta .63')
    if terminal.endpoint.contract['supply'] != {n: getattr(rule, n) for n in
                                               ('mode', 'initial_price', 'initial_stock', 'elasticity')}:
        raise ValueError('Terminal supply must be the inherited initial curve')
    for price in (float(rule.initial_price), float(rule.initial_price)*1.1):
        actual = np.asarray(rule.quantity(np.array([price])), dtype=float)
        original = P.H0[0]*(P.user_cost_rate*price/P.r_bar[0])**P.xi_supply[0]
        if actual.shape != (1,) or not np.isfinite(actual).all() or not np.isclose(actual[0], original, rtol=1e-12, atol=0):
            raise ValueError('Explicit supply differs from the approved initial parameter curve')
    if not _same(demographics, terminal_demographics):
        raise ValueError('Terminal and historical tail require identical frozen demographic primitives')
    if (demographics.start_year != 2023 or demographics.initial_person_state.year != 2023
            or not np.isfinite(demographics.scale_model_units_per_person)
            or demographics.scale_model_units_per_person <= 0):
        raise ValueError('Frozen demographics require a positive common scale and 2023 anchor')
    g = np.asarray(old.initial_state.g_pre)
    if (not np.isfinite(g).all() or np.any(g < 0) or g.sum() <= 0
            or not np.isclose(g.sum(), old.historical_conditioning.initial_mass, rtol=0, atol=2e-10)):
        raise ValueError('Historical demographic bridge must preserve its verified initial 2007 scale')
    final = receipt.get('final')
    if (final is None or not np.array_equal(np.asarray(final['prices']), terminal.policy.price)
            or not np.array_equal(np.asarray(final['fiscal_values']), [float(Q.pension)])):
        raise ValueError('Supplied terminal does not match the fresh final root coordinates')
    return grid, np.r_[psi, np.full(count-5, psi[-1])]


def solve_balanced_history(*, old_state, terminal, terminal_root_receipt,
                          demographic_primitives, terminal_demographic_primitives,
                          count, initial_prices, initial_pensions, price_bounds,
                          pension_bounds, audit_controls, market_tolerance, fiscal_tolerance,
                          market_slope, fiscal_slope, max_log_step, damping,
                          max_evaluations, deadline_monotonic, max_condition_number,
                          worsening_factor, final_reproduction_tolerance,
                          callback, observer=None, initial_jacobian=None):
    """Jointly solve dated housing and pensions with fixed payroll tax .179.

    Initial pension guesses are required for every date, including 2007.
    The terminal demographic object is supplied explicitly because the large
    endpoint does not store it; callers pin both objects to their source receipt.
    A process watchdog must bound a running evaluation. No targets are measured.
    """
    grid, psi = _validate(old_state, terminal, terminal_root_receipt,
        demographic_primitives, terminal_demographic_primitives, count)
    if not isinstance(audit_controls, TerminalAuditControls):
        raise ValueError('Explicit TerminalAuditControls required for every historical date')
    ceilings = dict(reconstruction_tolerance=5e-9,
        feasibility_projection_tolerance=1e-6, probability_tolerance=1e-12,
        occupied_mass_tolerance=1e-12, value_drop_tolerance=1e-7)
    for name, ceiling in ceilings.items():
        value = float(getattr(audit_controls, name))
        if not np.isfinite(value) or not 0 <= value <= ceiling:
            raise ValueError(f'{name} must be explicit, nonnegative and no looser than {ceiling}')
    if (type(max_evaluations) is not int or not 2 <= max_evaluations <= 8
            or not np.isfinite(deadline_monotonic) or deadline_monotonic <= time.monotonic()
            or not 0 < market_tolerance <= 2e-4 or not 0 < fiscal_tolerance <= 1e-6
            or not np.isfinite(final_reproduction_tolerance) or not 0 <= final_reproduction_tolerance <= 2e-10):
        raise ValueError('Require explicit bounded evaluations/deadline and unchanged tight market/fiscal/replay gates')
    if ((callback is not None and not callable(callback)) or (observer is not None and not callable(observer))):
        raise ValueError('Callbacks must be callable or None')
    prices, pensions = np.asarray(initial_prices, dtype=float), np.asarray(initial_pensions, dtype=float)
    lo, hi = map(float, price_bounds)
    if (prices.shape != (count,) or pensions.shape != (count,) or not np.isfinite([lo, hi]).all()
            or not 0 < lo < hi or not np.isfinite(prices).all() or np.any(prices < lo) or np.any(prices > hi)):
        raise ValueError('Explicit complete initial vectors and positive finite price bounds required')
    joined, primitive, checks, rent_domain = _runtime()
    P = old_state.parameters
    reference = SimpleNamespace(parameters=terminal.parameters, asset_price=float(terminal.policy.price[0]),
        renter_price=float(terminal.parameters.user_cost_rate)*float(terminal.policy.price[0]), equal_transfer=0.,
        psi_child=float(terminal.parameters.psi_child), state=joined.person_pf.PersonPFState(
            terminal.fixed_point.g_pre, terminal.fixed_point.persons))
    years = 2007+4*np.arange(count)
    last, best, last_trial, trial_count = None, None, None, 0

    def project(p):
        projected = rent_domain.project_price_path_to_positive_rents(np.clip(p, lo, hi),
            terminal=reference, minimum_rent_share=1e-6)[0]
        if np.any(projected > hi):
            raise ValueError('Declared price bounds cannot accommodate positive PF rents')
        return projected

    def evaluate(p, benefits):
        nonlocal last, last_trial, trial_count
        if time.monotonic() >= deadline_monotonic:
            raise TimeoutError('Historical path evaluation deadline reached')
        trial_count += 1
        last = None
        rents = joined.pf.rents_from_asset_prices(p, reference.asset_price, P)
        accounts, budget_checks, household_audits = [], [], []
        def observe(index, evaluation, dated, wealth_grid, shared):
            if (index != len(accounts) or float(dated.pension) != float(benefits[index])
                    or float(dated.tau_pay) != .179 or float(dated.property_tax_lump_sum_transfer) != 0.):
                raise RuntimeError('Dated observer did not receive the proposed pension/payroll path exactly once')
            diagnostics, household_gates = _household_checks(evaluation, dated, shared,
                wealth_grid, float(rents[index]), primitive, audit_controls)
            household_audits.append(dict(calendar_year=int(years[index]),
                diagnostics=diagnostics, gates=household_gates))
            if not all(household_gates.values()):
                raise RuntimeError(f'Historical household audit failed at {int(years[index])}: {household_gates}')
            budget_checks.append(diagnostics['budget'])
            accounts.append(fiscal_accounts(evaluation.g_current, dated))
            if observer is not None:
                observer(index, evaluation, dated, wealth_grid, shared)
        last = joined.evaluate_history_and_person_tail(years=years, prices=p, psi_path=psi,
            transfer_path=np.zeros(count), terminal_price=reference.asset_price, terminal_V=terminal.policy.V,
            base_parameters=P, b_grid=grid, initial_state=old_state.initial_state,
            historical_conditioning=old_state.historical_conditioning,
            initial_2023_persons=demographic_primitives.initial_person_state,
            demographic_primitives=demographic_primitives, supply_rule=old_state.supply_rule,
            birth_to_entry_conversion=1/2.1, pension_path=benefits,
            payroll_tax_path=np.full(count, .179), observer=observe)
        gates = joined.check_smoke_gates(last, expected_years=years.tolist())
        if len(accounts) != count:
            raise RuntimeError('Every date requires an actual post-choice fiscal ledger')
        market, fiscal = [], []
        for row, ledger in zip(last.rows, accounts):
            for name in ('payroll_tax_revenue', 'pension_outlays', 'pension_period_units', 'payroll_tax_rate'):
                if not np.isclose(row[name], ledger[name], rtol=0, atol=2e-10):
                    raise RuntimeError('Saved dated fiscal row differs from actual current household ledger')
            revenue, outlays = ledger['payroll_tax_revenue'], ledger['pension_outlays']
            fiscal.append(-1. if revenue == 0. and outlays != 0.
                          else (revenue-outlays)/max(abs(revenue), abs(outlays), 1e-12))
            if not np.isfinite(row['housing_supply']) or row['housing_supply'] <= 0:
                raise RuntimeError('Positive finite housing supply required')
            market.append((row['housing_demand']-row['housing_supply'])/row['housing_supply'])
        if not np.isfinite(market+fiscal).all():
            raise RuntimeError('Nonfinite historical housing/pension residual')
        distance = checks.terminal_convergence_diagnostics(last.person_tail, terminal=reference, psi_path=psi[4:])
        distance['last_pension_relative_gap'] = abs(float(benefits[-1])-float(terminal.parameters.pension))/float(terminal.parameters.pension)
        distance['pension_tail_tolerance_status'] = 'reported; explicit horizon certification remains caller-owned'
        last_trial = dict(trial=trial_count, prices=p.copy(), pensions=benefits.copy())
        return dict(market_residual=market, fiscal_residual=fiscal, mapping_valid=True,
            payload=dict(trial=trial_count, bellman_solves=last.bellman_solves, mapping_gates=gates,
                         terminal_distance=distance, fiscal_accounts=accounts, dated_budget_checks=budget_checks,
                         dated_household_audits=household_audits))

    def progress(record):
        nonlocal best
        if record.get('new_best'):
            best = last
        if callback is not None:
            callback(record)

    receipt = solve_social_security_path(closure='fixed_tax', initial_prices=prices,
        initial_fiscal_values=pensions, evaluate=evaluate, project_prices=project,
        fiscal_bounds=pension_bounds, market_tolerance=market_tolerance, fiscal_tolerance=fiscal_tolerance,
        market_slope=market_slope, fiscal_slope=fiscal_slope, max_log_step=max_log_step,
        damping=damping, max_evaluations=max_evaluations, deadline_monotonic=deadline_monotonic,
        max_condition_number=max_condition_number, worsening_factor=worsening_factor,
        final_reproduction_tolerance=final_reproduction_tolerance, callback=progress,
        initial_jacobian=initial_jacobian)
    final = receipt['final']
    matched = bool(final is not None and last is not None and last_trial is not None
        and final['payload']['trial'] == last_trial['trial']
        and np.array_equal(final['prices'], last_trial['prices'])
        and np.array_equal(final['fiscal_values'], last_trial['pensions']))
    selected = final if matched else receipt['best']
    distance = None if selected is None else selected['payload']['terminal_distance']
    receipt.update(schema='e5f_balanced_history_v1', years=years.tolist(), psi_path=psi.tolist(),
        audit_controls=vars(audit_controls).copy(),
        payroll_tax=.179, property_tax_period=.04, property_tax_rebate=0.,
        fresh_path_matches_final=matched, returned_path='fresh_final' if matched else 'best_diagnostic',
        finite_horizon_market_fiscal_converged=bool(receipt['converged'] and matched),
        terminal_distance=distance, terminal_distance_passed=bool(distance and distance['all_checks_pass']),
        horizon_verified=False, historical_equilibrium_certified=False, calibrated_smm=False,
        fiscal_population='actual dated post-choice household heads',
        fiscal_residual_definition='(payroll revenue - pension outlays)/max(abs(revenue),abs(outlays),1e-12)',
        fiscal_zero_budget_rule='both zero: zero; zero revenue with positive outlays: -1',
        initial_state='verified stationary conditional states reweighted to observed 2007 ages',
        demographic_scale_validation='same frozen terminal primitives; exact 2023 household/head age gate',
        demographic_continuation='serialized 2023 anchor with annual projections through the last empirical year, held fixed afterward',
        preference_continuation='supplied announced historical path, held constant after 2023')
    return BalancedHistoryResult(last if matched else best, receipt)
