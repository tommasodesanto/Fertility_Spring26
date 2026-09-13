"""Isolated rebated-tax adapter for successive permanent-surprise forecasts.

The caller supplies the terminal endpoint and all numerical controls.  This
module changes no calibration, target, demographic primitive, or core model
code.  A finite dated root is useful provisional evidence even when the supplied
terminal boundary is not close enough to certify the horizon; those statuses
are therefore recorded separately and neither one implies production use.
"""
from __future__ import annotations

import copy
from dataclasses import replace
from types import SimpleNamespace
import time

import numpy as np

import e5f_successive_surprises as base


InheritedState = base.InheritedState
SurpriseResult = base.SurpriseResult
local_conditioning = base.local_conditioning

PAYROLL_TAX = 0.179
FISCAL_RESIDUAL_SCALE = 200.0


def _runtime():
    return base._runtime()


def _path_root_solver():
    # The launcher's pinned sys.path selects tmp/e5f_matched_pf.
    from e5f_matched_pf_path_root import solve_price_path
    return solve_price_path


def _terminal_parts(terminal):
    """Accept the retained terminal result or its endpoint-style namespace."""
    parameters = getattr(terminal, "parameters", getattr(terminal, "P", None))
    policy = getattr(terminal, "policy", None)
    if policy is None and hasattr(terminal, "policyV"):
        policy = SimpleNamespace(V=terminal.policyV)
    if parameters is None or policy is None or not hasattr(policy, "V"):
        raise ValueError("Terminal must supply parameters and a terminal policy/V")
    asset_price = getattr(terminal, "asset_price", None)
    if asset_price is None:
        if hasattr(terminal, "price"):
            asset_price = float(np.asarray(terminal.price, dtype=float).reshape(-1)[0])
        elif not hasattr(policy, "price"):
            raise ValueError("Terminal must supply its asset price")
        else:
            asset_price = float(np.asarray(policy.price, dtype=float).reshape(-1)[0])
    fixed = getattr(terminal, "fixed_point", None)
    state = getattr(terminal, "state", None)
    if state is None and fixed is not None:
        state = fixed
    if state is not None and (not hasattr(state, "g_pre") or not hasattr(state, "persons")):
        raise ValueError("Optional terminal state must contain households and persons")
    return parameters, policy, float(asset_price), state


def rebated_tax_accounts(*, property_tax_revenue, transfer_per_head, head_mass):
    """Return the equal per-head property-tax rebate ledger.

    ``head_mass`` is the mass of all current household heads, renters and owners.
    Property-tax revenue already includes tax paid on rental and owner housing in
    the core period evaluator; this adapter returns all of that revenue equally.
    """
    revenue = float(property_tax_revenue)
    transfer = float(transfer_per_head)
    heads = float(head_mass)
    if (not np.isfinite([revenue, transfer, heads]).all() or revenue < 0
            or transfer < 0 or heads <= 0):
        raise ValueError("Finite nonnegative revenue/transfer and positive head mass required")
    outlays = transfer * heads
    residual = revenue - outlays
    scale = max(abs(revenue), abs(outlays), 1e-12)
    return dict(property_tax_revenue=revenue,
                equal_transfer_period_units=transfer,
                household_head_mass=heads,
                equal_transfer_outlays=outlays,
                government_budget_residual=residual,
                scaled_government_budget_residual=residual / scale,
                implied_equal_transfer=revenue / heads)


def dated_residual(*, demand, supply, payroll_accounts, tax_accounts):
    """Housing, PAYGO, and equal-rebate residuals in root-coordinate order."""
    demand, supply = float(demand), float(supply)
    if not np.isfinite([demand, supply]).all() or supply <= 0:
        raise ValueError("Finite demand and positive finite housing supply required")
    payroll_revenue = float(payroll_accounts["payroll_tax_revenue"])
    pension_outlays = float(payroll_accounts["pension_outlays"])
    tax_revenue = float(tax_accounts["property_tax_revenue"])
    transfer_outlays = float(tax_accounts["equal_transfer_outlays"])
    if not np.isfinite([payroll_revenue, pension_outlays,
                        tax_revenue, transfer_outlays]).all():
        raise ValueError("Finite dated government accounts required")
    return np.array([
        (demand - supply) / supply,
        FISCAL_RESIDUAL_SCALE * (payroll_revenue - pension_outlays)
        / max(abs(payroll_revenue), abs(pension_outlays), 1e-12),
        FISCAL_RESIDUAL_SCALE * (tax_revenue - transfer_outlays)
        / max(abs(tax_revenue), abs(transfer_outlays), 1e-12),
    ], dtype=float)


def stack_dated_residuals(blocks):
    """Stack date rows into [all housing, all PAYGO, all rebates]."""
    values = np.asarray(blocks, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3 or not np.isfinite(values).all():
        raise ValueError("Dated residual blocks must be a finite N-by-3 array")
    return values.T.reshape(-1)


def evaluate_forecast(*, inherited, old_state, demographics, prices, pensions,
                      transfers, psi, terminal, observer=None,
                      demographic_evaluator=None):
    """Evaluate one constant-current-preference forecast with explicit rebates."""
    _, joined, _, _, _ = _runtime()
    pf, person = joined.pf, joined.person_pf
    P, grid = old_state.parameters, old_state.b_grid
    terminal_parameters, terminal_policy, terminal_price, _ = _terminal_parts(terminal)
    if (not hasattr(terminal_parameters, "psi_child")
            or not np.isclose(float(terminal_parameters.psi_child), float(psi),
                              rtol=0, atol=1e-14)):
        raise ValueError("Forecast terminal must use the current constant preference")
    start = inherited.year
    if start not in (2007, 2011, 2015, 2019, 2023):
        raise ValueError("No restart beyond the fixed historical observation window")
    p = np.asarray(prices, dtype=float)
    benefits = np.asarray(pensions, dtype=float)
    rebates = np.asarray(transfers, dtype=float)
    h = (2023 - start) // 4
    if (p.ndim != 1 or len(p) < h + 2 or benefits.shape != p.shape
            or rebates.shape != p.shape):
        raise ValueError("Forecast must contain equally sized dated price, pension, and rebate paths")
    if (not np.isfinite(psi) or not np.isfinite(p).all() or np.any(p <= 0)
            or not np.isfinite(benefits).all() or np.any(benefits < 0)
            or not np.isfinite(rebates).all() or np.any(rebates < 0)):
        raise ValueError("Finite preference and nonnegative fiscal paths with positive prices required")
    psi_path = np.full(len(p), float(psi))
    taxes = np.full(len(p), PAYROLL_TAX)
    years = start + 4 * np.arange(len(p))
    tail_values, backward = pf.backward_value_path(
        prices=p[h:], rents=pf.rents_from_asset_prices(p[h:], terminal_price, P),
        psi_path=psi_path[h:], terminal_V=terminal_policy.V,
        base_parameters=P, b_grid=grid, transfer_path=rebates[h:],
        pension_path=benefits[h:], payroll_tax_path=taxes[h:])
    if h:
        conditioning = replace(
            local_conditioning(old_state.historical_conditioning, start),
            observer=observer)
        history = pf.evaluate_path_at_prices(
            prices=p[:h], psi_path=psi_path[:h], transfer_path=rebates[:h],
            terminal_price=float(p[h]), terminal_V=tail_values[0],
            base_parameters=P, b_grid=grid, initial_state=inherited.households,
            supply_rule=old_state.supply_rule, birth_to_entry_conversion=1 / 2.1,
            historical_conditioning=conditioning, pension_path=benefits[:h],
            payroll_tax_path=taxes[:h])
        g = history.terminal_state.g_pre
    else:
        g = inherited.households.g_pre
        history = SimpleNamespace(rows=[], values=[], bellman_solves=0,
            maximum_mass_accounting_error=0., maximum_policy_reproduction_error=0.,
            maximum_feasibility_projection_mass=0.)
    people = demographics.initial_person_state.validated()
    if people.year != 2023:
        raise ValueError("Retained demographic anchor must remain 2023")
    if not h and hasattr(inherited.households, "persons"):
        inherited_people = inherited.households.persons
        if (inherited_people.year != 2023
                or not np.array_equal(inherited_people.persons, people.persons)
                or not np.array_equal(inherited_people.heads, people.heads)):
            raise ValueError("Inherited 2023 person anchor differs; no silent reset allowed")
        people = inherited_people
    heads = person.aggregate_heads_to_model_age_cells(
        people, age_start=int(P.age_start), cell_width=int(P.da),
        number_of_cells=int(P.J))
    gap = float(np.max(np.abs(g.sum(axis=(0, 1, 2, 4, 5, 6)) - heads)))
    if not np.isfinite(gap) or gap > 2e-9:
        raise RuntimeError("Inherited 2023 head-age bridge failed")

    def tail_observer(i, *args):
        if observer is not None:
            observer(h + i, *args)

    evaluator = demographic_evaluator or person.evaluate_path_at_prices_person_demography
    tail = evaluator(
        prices=p[h:], psi_path=psi_path[h:], transfer_path=rebates[h:],
        terminal_price=terminal_price, terminal_V=terminal_policy.V,
        base_parameters=P, b_grid=grid,
        initial_state=person.PersonPFState(g_pre=g.copy(), persons=people),
        demographic_primitives=demographics, supply_rule=old_state.supply_rule,
        precomputed_value_path=tail_values,
        observer=tail_observer if observer is not None else None,
        pension_path=benefits[h:], payroll_tax_path=taxes[h:])
    result = joined.ConditionalHistoryEvaluation(
        history=history, person_tail=tail,
        rows=[dict(r) for r in history.rows]
             + [dict(r, period=int(r["period"]) + h) for r in tail.rows],
        values=(history.values[:-1] if h else []) + tail.values,
        bellman_solves=history.bellman_solves + backward + tail.bellman_solves,
        initial_2023_age_head_gap=gap)
    joined.check_smoke_gates(result, expected_years=years.tolist())
    return result


def first_period_state(*, inherited, old_state, demographics, path, prices,
                       pensions, transfers, psi, demographic_evaluator=None):
    """Replay the realized period under its accepted expected next price and V."""
    _, joined, _, _, _ = _runtime()
    pf, person = joined.pf, joined.person_pf
    P = old_state.parameters
    common = dict(
        prices=[float(prices[0])], psi_path=[float(psi)],
        transfer_path=[float(transfers[0])], terminal_price=float(prices[1]),
        terminal_V=path.values[1], base_parameters=P, b_grid=old_state.b_grid,
        initial_state=inherited.households, supply_rule=old_state.supply_rule,
        pension_path=[float(pensions[0])], payroll_tax_path=[PAYROLL_TAX])
    if inherited.year < 2023:
        replay = pf.evaluate_path_at_prices(
            **common, birth_to_entry_conversion=1 / 2.1,
            historical_conditioning=local_conditioning(
                old_state.historical_conditioning, inherited.year))
        if (replay.maximum_mass_accounting_error > 2e-8
                or replay.maximum_policy_reproduction_error > 2e-10
                or replay.maximum_feasibility_projection_mass > 1e-6):
            raise RuntimeError("First-period historical replay gate failed")
        state = replay.terminal_state
        if inherited.year == 2019:
            state = person.PersonPFState(
                g_pre=state.g_pre,
                persons=copy.deepcopy(demographics.initial_person_state))
    else:
        if not hasattr(inherited.households, "persons"):
            raise ValueError("2023 surprise requires inherited household/person state")
        evaluator = demographic_evaluator or person.evaluate_path_at_prices_person_demography
        replay = evaluator(**common, demographic_primitives=demographics,
                           precomputed_value_path=path.values[:2])
        state = replay.terminal_state
        limits = dict(maximum_person_identity_error=2e-9,
            maximum_head_identity_error=2e-9,
            maximum_household_person_head_gap=2e-9,
            maximum_age_head_gap=2e-9,
            maximum_policy_reproduction_error=2e-10,
            maximum_feasibility_projection_mass=1e-6)
        for name, limit in limits.items():
            value = getattr(replay, name)
            if not np.isfinite(value) or value > limit:
                raise RuntimeError("First-period person replay gate failed: " + name)
    if not np.allclose(replay.values[0], path.values[0], rtol=0, atol=2e-10):
        raise RuntimeError("First-period continuation value differs from accepted forecast")
    keys = ("asset_price", "renter_price", "housing_demand", "owner_rate",
            "birth_children_topcode_adjusted", "pension_period_units",
            "payroll_tax_revenue", "pension_outlays", "property_tax_revenue",
            "equal_transfer_outlays")
    for key in keys:
        if key in replay.rows[0] or key in path.rows[0]:
            if key not in replay.rows[0] or key not in path.rows[0] or not np.isclose(
                    float(replay.rows[0][key]), float(path.rows[0][key]),
                    rtol=0, atol=2e-10):
                raise RuntimeError("First-period replay differs: " + key)
    return InheritedState(inherited.year + 4, state)


def solve_rebated_forecast(*, inherited, psi, old_state, terminal,
                           demographic_primitives, count, initial_prices,
                           initial_pensions, initial_transfers, audit_controls,
                           root_controls, deadline_monotonic, callback=None,
                           observer=None, demographic_evaluator=None):
    """Jointly root every dated price, pension, and equal rebate.

    The terminal is caller-supplied.  Terminal-distance failure is reported as a
    boundary diagnostic and does not discard a finite, exactly reproduced dated
    root or its first-period replay.
    """
    _, joined, primitive, checks, rent_domain = _runtime()
    from e5f_social_security import fiscal_accounts
    from e5f_balanced_terminal import TerminalAuditControls, _household_checks

    terminal_parameters, terminal_policy, terminal_price, terminal_state = _terminal_parts(terminal)
    if (not np.isfinite(psi)
            or not np.isclose(float(terminal_parameters.psi_child), float(psi),
                              rtol=0, atol=1e-14)):
        raise ValueError("Forecast terminal must use the current constant preference")
    if not isinstance(audit_controls, TerminalAuditControls):
        raise ValueError("Explicit retained household audit controls required")
    audit_ceilings = dict(reconstruction_tolerance=5e-9,
        feasibility_projection_tolerance=1e-6, probability_tolerance=1e-12,
        occupied_mass_tolerance=1e-12, value_drop_tolerance=1e-7)
    for name, ceiling in audit_ceilings.items():
        value = float(getattr(audit_controls, name))
        if not np.isfinite(value) or not 0 <= value <= ceiling:
            raise ValueError("Invalid retained household audit: " + name)
    if type(count) is not int or count < 2:
        raise ValueError("Forecast count must include at least two dates")
    if not np.isfinite(deadline_monotonic) or deadline_monotonic <= time.monotonic():
        raise ValueError("Expired forecast deadline")
    p0 = np.asarray(initial_prices, dtype=float)
    b0 = np.asarray(initial_pensions, dtype=float)
    t0 = np.asarray(initial_transfers, dtype=float)
    if p0.shape != (count,) or b0.shape != (count,) or t0.shape != (count,):
        raise ValueError("Complete dated price, pension, and rebate guesses required")
    controls = dict(root_controls)
    required = {"price_bounds", "pension_bounds", "transfer_bounds", "slope",
        "market_tolerance", "max_log_step", "damping", "max_evaluations",
        "max_condition_number", "worsening_factor", "final_reproduction_tolerance"}
    missing = required - controls.keys()
    if missing:
        raise ValueError("Missing explicit root controls: " + ", ".join(sorted(missing)))
    bounds = []
    for name in ("price_bounds", "pension_bounds", "transfer_bounds"):
        pair = tuple(controls.pop(name))
        if len(pair) != 2 or not np.isfinite(pair).all() or not 0 < pair[0] < pair[1]:
            raise ValueError(name + " must be explicit finite positive bounds")
        bounds.append(pair)
    if type(controls["max_evaluations"]) is not int or not 2 <= controls["max_evaluations"] <= 8:
        raise ValueError("Retain the bounded 2--8 mapping budget")
    if not 0 < float(controls["market_tolerance"]) <= 2e-4:
        raise ValueError("Joint root tolerance cannot exceed 2e-4")
    if not 0 <= float(controls["final_reproduction_tolerance"]) <= 2e-10:
        raise ValueError("Final reproduction tolerance cannot exceed 2e-10")
    optional = {"initial_jacobian", "default_jacobian"}
    unknown = controls.keys() - (required - {"price_bounds", "pension_bounds", "transfer_bounds"}) - optional
    if unknown:
        raise ValueError("Unknown root controls: " + ", ".join(sorted(unknown)))
    for values, pair, label in ((p0, bounds[0], "price"),
                                (b0, bounds[1], "pension"),
                                (t0, bounds[2], "transfer")):
        if (not np.isfinite(values).all() or np.any(values < pair[0])
                or np.any(values > pair[1])):
            raise ValueError(f"Initial {label} guesses must lie inside explicit bounds")
    endpoint_fields = dict(
        parameters=terminal_parameters, asset_price=terminal_price,
        renter_price=float(terminal_parameters.user_cost_rate) * terminal_price,
        equal_transfer=float(getattr(terminal_parameters,
                                     "property_tax_lump_sum_transfer", 0.)),
        psi_child=float(psi))
    if terminal_state is not None:
        endpoint_fields["state"] = joined.person_pf.PersonPFState(
            terminal_state.g_pre, terminal_state.persons)
    endpoint = SimpleNamespace(**endpoint_fields)
    final_path = coordinates = None
    trial = 0

    def project(x):
        x = np.asarray(x, dtype=float).copy()
        if x.shape != (3 * count,):
            raise ValueError("Joint root projection received wrong coordinate shape")
        for j, pair in enumerate(bounds):
            x[j * count:(j + 1) * count] = np.clip(
                x[j * count:(j + 1) * count], pair[0], pair[1])
        x[:count] = rent_domain.project_price_path_to_positive_rents(
            x[:count], terminal=endpoint, minimum_rent_share=1e-6)[0]
        if np.any(x[:count] > bounds[0][1]) or np.any(x[:count] < bounds[0][0]):
            raise ValueError("Positive-rent projection exceeds explicit price bounds")
        return x

    def evaluate(x):
        nonlocal final_path, coordinates, trial
        final_path = coordinates = None
        trial += 1
        prices = np.asarray(x[:count], dtype=float)
        pensions = np.asarray(x[count:2 * count], dtype=float)
        transfers = np.asarray(x[2 * count:], dtype=float)
        accounts, head_masses, tax_ledgers, audits = [], [], [], []
        rents = joined.pf.rents_from_asset_prices(
            prices, endpoint.asset_price, old_state.parameters)

        def observe(i, evaluation, parameters, grid, shared):
            if (i != len(accounts) or parameters.pension != float(pensions[i])
                    or parameters.tau_pay != PAYROLL_TAX
                    or parameters.property_tax_lump_sum_transfer != float(transfers[i])):
                raise RuntimeError("Wrong dated fiscal path")
            diagnostics, gates = _household_checks(
                evaluation, parameters, shared, grid, float(rents[i]),
                primitive, audit_controls)
            if not all(gates.values()):
                raise RuntimeError("Dated household audit failed")
            payroll = fiscal_accounts(evaluation.g_current, parameters)
            accounts.append(payroll)
            head_masses.append(float(np.sum(evaluation.g_current)))
            audits.append(dict(year=inherited.year + 4 * i,
                               diagnostics=diagnostics, gates=gates))
            if observer is not None:
                observer(i, evaluation, parameters, grid, shared)

        result = evaluate_forecast(
            inherited=inherited, old_state=old_state,
            demographics=demographic_primitives, prices=prices,
            pensions=pensions, transfers=transfers, psi=psi, terminal=terminal,
            observer=observe, demographic_evaluator=demographic_evaluator)
        if len(accounts) != count:
            raise RuntimeError("Missing dated fiscal ledger")
        blocks = []
        for i, (row, payroll) in enumerate(zip(result.rows, accounts)):
            tax = rebated_tax_accounts(
                property_tax_revenue=row["property_tax_revenue"],
                transfer_per_head=transfers[i],
                head_mass=head_masses[i])
            # Native paths report exact outlays.  Require equality if present so
            # the root cannot balance a ledger different from household budgets.
            if "equal_transfer_outlays" in row and not np.isclose(
                    float(row["equal_transfer_outlays"]), tax["equal_transfer_outlays"],
                    rtol=0, atol=2e-10):
                raise RuntimeError("Native equal-transfer ledger differs from all-head accounting")
            tax_ledgers.append(tax)
            blocks.append(dated_residual(
                demand=row["housing_demand"], supply=row["housing_supply"],
                payroll_accounts=payroll, tax_accounts=tax))
        if terminal_state is None:
            boundary = dict(status="unverified", all_checks_pass=False,
                reason="Caller supplied continuation value/price without a terminal state")
        else:
            boundary = copy.deepcopy(checks.terminal_convergence_diagnostics(
                result.person_tail, terminal=endpoint,
                psi_path=np.full(len(result.person_tail.rows), float(psi))))
            boundary_passed = bool(boundary.get("all_checks_pass", False))
            boundary["status"] = "passed" if boundary_passed else "not_converged"
        final_path = result
        coordinates = (np.asarray(x, dtype=float).copy(), trial)
        return dict(
            residual=stack_dated_residuals(blocks),
            mapping_valid=True,
            payload=dict(trial=trial, boundary_diagnostics=boundary,
                         boundary_status=boundary["status"],
                         dated_household_audits=audits,
                         payroll_accounts=accounts,
                         rebated_tax_accounts=tax_ledgers))

    initial = np.concatenate([p0, b0, t0])
    if "default_jacobian" not in controls:
        controls["default_jacobian"] = np.diag(np.concatenate([
            np.full(count, -float(controls["slope"])),
            np.full(2 * count, -FISCAL_RESIDUAL_SCALE)]))
    receipt = _path_root_solver()(
        initial_prices=initial, evaluate=evaluate, project=project,
        deadline_monotonic=deadline_monotonic, callback=callback, **controls)
    final = receipt.get("final")
    matched = bool(final is not None and final_path is not None and coordinates is not None
                   and final["payload"]["trial"] == coordinates[1]
                   and np.array_equal(final["prices"], coordinates[0]))
    finite = bool(receipt.get("converged") and matched)
    boundary_status = (final["payload"]["boundary_status"] if finite
                       else "not_evaluated_at_accepted_root")
    receipt.update(
        schema="e5f_rebated_permanent_surprise_v1",
        start_year=inherited.year, psi=float(psi),
        expected_psi_path=[float(psi)] * count,
        payroll_tax_path=[PAYROLL_TAX] * count,
        finite_horizon_market_fiscal_converged=finite,
        finite_root_status="passed" if finite else "failed",
        boundary_status=boundary_status,
        terminal_distance_passed=bool(finite and boundary_status == "passed"),
        horizon_verified=False, production_eligible=False,
        information="Current preference persists; later preference surprises are not anticipated",
        fiscal_information="Payroll tax is 0.179; all property-tax revenue is rebated equally to current household heads")
    if not finite:
        return SurpriseResult(final_path, receipt, None, None)
    accepted = final["prices"]
    prices = accepted[:count]
    pensions = accepted[count:2 * count]
    transfers = accepted[2 * count:]
    next_state = first_period_state(
        inherited=inherited, old_state=old_state,
        demographics=demographic_primitives, path=final_path, prices=prices,
        pensions=pensions, transfers=transfers, psi=psi,
        demographic_evaluator=demographic_evaluator)
    realized = dict(final_path.rows[0], forecast_vintage_year=inherited.year,
        expected_next_asset_price=float(prices[1]),
        expected_next_pension=float(pensions[1]),
        expected_next_equal_transfer=float(transfers[1]),
        expected_constant_psi=float(psi), boundary_status=boundary_status,
        production_eligible=False)
    return SurpriseResult(final_path, receipt, next_state, realized)
