"""Original-household-law terminal steady state with the native birth queue.

The three root coordinates are the asset price, four-year pension, and equal
property-tax rebate.  For every trial, the stationary household distribution
is first solved at unit mass.  Population is then fixed by the unchanged
housing supply curve, and all four waiting birth-queue vintages are scaled by
that same economically determined population level.

This adapter never uses the post-2023 person-demography endpoint.  It performs
no historical conditioning and saves a terminal checkpoint only after a fresh
root replay and an exact one-period application of the original queue operator.
"""
from __future__ import annotations

import copy
import gzip
import json
import math
import os
import pickle
from pathlib import Path
from types import SimpleNamespace
import time
from typing import Any

import numpy as np


PAYROLL_TAX = 0.179
REPLACEMENT_FERTILITY = 2.1
RESIDUAL_SCALE = 200.0
QUEUE_WAITING_SLOTS = 4
ROOT_TOLERANCE = 1e-5
UNSCALED_ROOT_TOLERANCE = ROOT_TOLERANCE / RESIDUAL_SCALE


def _runtime():
    """Load and configure only the launcher's pinned scientific runtime."""
    import e5f_rebated_surprises as rebated

    _, joined, primitive, _, _ = rebated._runtime()
    pf = joined.pf
    return rebated, joined, primitive, pf, pf.calendar


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def _save_json(path: Path, payload: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _save_pickle(path: Path, payload: Any) -> None:
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with gzip.open(temporary, "wb") as stream:
        pickle.dump(payload, stream, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temporary, path)


def _relative_gap(left: float, right: float) -> float:
    scale = max(abs(float(left)), abs(float(right)), 1e-12)
    return (float(left) - float(right)) / scale


def _validate_supply(old: Any) -> None:
    P = old.parameters
    supply = old.supply_rule
    if int(P.I) != 1:
        raise ValueError("Original-queue terminal requires exactly one housing market")
    if (float(P.period_years) != 4.0 or not bool(P.scale_flows_to_period)
            or not math.isclose(float(P.tau_H), 0.04, rel_tol=0.0, abs_tol=1e-15)
            or not math.isclose(float(P.tau_pay), PAYROLL_TAX,
                                rel_tol=0.0, abs_tol=1e-15)):
        raise ValueError(
            "Require four-year flows, 1% annual property tax, and payroll tax .179"
        )
    if any(bool(getattr(P, name, False)) for name in
           ("joint_nested_choice", "fertility_nest_choice", "two_shock_choice")):
        raise ValueError("Terminal requires the approved sequential household law")
    if not bool(getattr(P, "exhaustive_saving_control", False)):
        raise ValueError("Terminal requires exhaustive_saving_control=True")
    if (getattr(supply, "mode", None) != "static-elastic"
            or not math.isclose(float(supply.elasticity), 0.63,
                                rel_tol=0.0, abs_tol=1e-15)):
        raise ValueError("The unchanged static-elastic supply curve is required")
    H0 = np.asarray(P.H0, dtype=float)
    r_bar = np.asarray(P.r_bar, dtype=float)
    xi = np.asarray(P.xi_supply, dtype=float)
    if (H0.shape != (1,) or r_bar.shape != (1,) or xi.shape != (1,)
            or not np.isfinite(np.r_[H0, r_bar, xi]).all()
            or np.any(H0 <= 0) or np.any(r_bar <= 0)
            or not math.isclose(float(xi[0]), 0.63, rel_tol=0.0, abs_tol=1e-15)):
        raise ValueError("Parameter supply curve is invalid or has been changed")
    for q in (float(supply.initial_price), 1.1 * float(supply.initial_price)):
        explicit = float(np.asarray(supply.quantity(np.array([q]))).reshape(-1)[0])
        parameter = float(H0[0] * (float(P.user_cost_rate) * q / r_bar[0]) ** xi[0])
        if not math.isclose(explicit, parameter, rel_tol=1e-12, abs_tol=0.0):
            raise ValueError("Explicit supply rule differs from the calibrated curve")


def _validated_controls(controls: dict[str, Any], start: Any, old: Any):
    values = dict(controls)
    missing = {
        "price_bounds", "pension_bounds", "transfer_bounds", "max_log_step",
        "damping", "max_evaluations", "max_condition_number",
        "worsening_factor", "final_reproduction_tolerance",
    } - values.keys()
    if missing:
        raise ValueError("Missing terminal root controls: " + ", ".join(sorted(missing)))
    bounds = []
    for name in ("price_bounds", "pension_bounds", "transfer_bounds"):
        pair = tuple(float(item) for item in values.pop(name))
        if len(pair) != 2 or not np.isfinite(pair).all() or not 0 < pair[0] < pair[1]:
            raise ValueError(name + " must be an explicit finite positive pair")
        bounds.append(pair)
    maximum = values["max_evaluations"]
    if isinstance(maximum, bool) or maximum not in (16, 24):
        raise ValueError("Terminal root budget must be exactly 16 or 24 evaluations")
    supplied_tolerance = float(values.pop("market_tolerance", ROOT_TOLERANCE))
    if not math.isfinite(supplied_tolerance) or supplied_tolerance <= 0:
        raise ValueError("Root tolerance must be finite and positive")
    values.pop("fiscal_tolerance", None)
    values.pop("fiscal_slope", None)
    slope = float(values.pop("slope", values.pop("market_slope", 1.0)))
    values.pop("automatic_fiscal_polish", None)
    allowed = {
        "max_log_step", "damping", "max_evaluations", "max_condition_number",
        "worsening_factor", "final_reproduction_tolerance", "initial_jacobian",
    }
    unknown = values.keys() - allowed
    if unknown:
        raise ValueError("Unknown terminal root controls: " + ", ".join(sorted(unknown)))
    if not 0 <= float(values["final_reproduction_tolerance"]) <= 2e-10:
        raise ValueError("Fresh root reproduction tolerance cannot exceed 2e-10")
    if start is None:
        initial = np.array([
            float(np.asarray(old.policy.price, dtype=float).reshape(-1)[0]),
            float(old.parameters.pension),
            float(old.parameters.property_tax_lump_sum_transfer),
        ])
    else:
        initial = np.asarray(start, dtype=float)
    if initial.shape != (3,) or not np.isfinite(initial).all() or np.any(initial <= 0):
        raise ValueError("Terminal start must contain positive finite (price,pension,rebate)")
    for value, (lower, upper), label in zip(
            initial, bounds, ("price", "pension", "transfer")):
        if not lower <= value <= upper:
            raise ValueError(f"Initial terminal {label} lies outside inherited bounds")
    arguments = {key: values[key] for key in allowed - {"initial_jacobian"}}
    arguments.update(slope=slope, market_tolerance=min(supplied_tolerance, ROOT_TOLERANCE))
    return bounds, initial, arguments, values.get("initial_jacobian")


def _candidate_reference(candidate: Any) -> dict[str, Any]:
    return {
        "asset_price": candidate.asset_price,
        "renter_price": candidate.renter_price,
        "pension_period_units": candidate.pension_period,
        "equal_transfer_period_units": candidate.transfer,
        "population_households": candidate.population_scale,
        "housing_demand": candidate.housing_demand,
        "housing_supply": candidate.housing_supply,
        "owner_rate": candidate.owner_rate,
        "birth_children_topcode_adjusted": candidate.adjusted_births,
        "entry_flow": candidate.entry_flow,
        "raw_birth_queue_flow": candidate.raw_queue_flow,
        "renewal_ratio": candidate.renewal_ratio,
        "payroll_tax_rate": PAYROLL_TAX,
    }


def endpoint_reference(endpoint: Any) -> dict[str, Any]:
    """Return a tensor-free JSON series for plotting the terminal reference."""
    receipt = getattr(endpoint, "receipt", None)
    if isinstance(receipt, dict) and isinstance(receipt.get("endpoint_reference"), dict):
        return copy.deepcopy(receipt["endpoint_reference"])
    if hasattr(endpoint, "endpoint_reference"):
        return copy.deepcopy(endpoint.endpoint_reference)
    return _candidate_reference(endpoint)


def _evaluate_trial(*, old: Any, psi: float, coordinates: np.ndarray,
                    audit: Any, deadline: float, trial: int) -> Any:
    if time.monotonic() >= deadline:
        raise TimeoutError("Terminal evaluation deadline reached")
    rebated, _, primitive, pf, calendar = _runtime()
    from e5f_balanced_terminal import _household_checks
    from e5f_social_security import bind_social_security_income, fiscal_accounts

    q, pension, transfer = (float(item) for item in coordinates)
    P = copy.deepcopy(old.parameters)
    P.psi_child = float(psi)
    P.property_tax_lump_sum_transfer = transfer
    # This binding must precede both shared-input construction and Bellman solution.
    bind_social_security_income(P, pension_period=pension, payroll_tax=PAYROLL_TAX)
    grid = np.asarray(old.b_grid, dtype=float)
    shared = calendar.model.precompute_shared(P, grid)
    price = np.array([q], dtype=float)
    solution = calendar.model.solve_markov_income_at_prices(price, P, grid, SD=shared)
    solved_grid = np.asarray(solution.b_grid, dtype=float)
    if solved_grid.shape != grid.shape or not np.array_equal(solved_grid, grid):
        raise RuntimeError("Terminal household solve changed the inherited asset grid")
    if hasattr(solution, "fert2_probs"):
        P._fert2_probs = np.asarray(solution.fert2_probs, dtype=float).copy()
    policy = calendar.policy_from_solution(solution, price, P, grid, shared)
    unit_g, reconstruction = calendar.reconstruct_stationary_pre_fertility(
        solution, policy, P, grid, shared
    )
    unit_g = np.asarray(unit_g, dtype=float)
    unit_mass = float(unit_g.sum())
    if not math.isclose(unit_mass, 1.0, rel_tol=0.0, abs_tol=2e-10):
        raise RuntimeError(f"Stationary household solution has mass {unit_mass:.12g}, not one")
    unit = calendar.evaluate_period(
        price, unit_g, P, grid, shared, calendar.SolveCounter(),
        supply_rule=old.supply_rule, supplied_policy=policy,
    )
    unit_demand = float(np.asarray(unit.demand_by_loc).reshape(-1)[0])
    housing_supply = float(np.asarray(old.supply_rule.quantity(price)).reshape(-1)[0])
    if not np.isfinite([unit_demand, housing_supply]).all() or min(unit_demand, housing_supply) <= 0:
        raise RuntimeError("Terminal unit demand and inherited housing supply must be positive")
    population_scale = housing_supply / unit_demand
    scaled_g = population_scale * unit_g
    current = calendar.evaluate_period(
        price, scaled_g, P, grid, shared, calendar.SolveCounter(),
        supply_rule=old.supply_rule, supplied_policy=policy,
    )
    accounting = primitive.transition.calendar_topcode_birth_accounting(
        current.g_pre, current.g_post_fertility, float(current.births), P
    )
    adjusted_births = float(accounting["topcode_adjusted_birth_children"])
    raw_births = float(current.births)
    entry_flow = float(np.sum(scaled_g[:, :, :, 0, :, :, :]))
    if entry_flow <= 0:
        raise RuntimeError("Terminal stationary entrant flow must be positive")
    renewal_ratio = adjusted_births / (REPLACEMENT_FERTILITY * entry_flow)
    payroll = fiscal_accounts(current.g_current, P)
    head_mass = float(np.sum(current.g_current))
    tax_revenue = float(calendar.model.property_tax_revenue_from_distribution(
        current.g_current, current.policy.hR_pol, current.policy.price, P
    ))
    tax = rebated.rebated_tax_accounts(
        property_tax_revenue=tax_revenue,
        transfer_per_head=transfer,
        head_mass=head_mass,
    )
    pension_gap = _relative_gap(payroll["payroll_tax_revenue"], payroll["pension_outlays"])
    rebate_gap = _relative_gap(tax["property_tax_revenue"], tax["equal_transfer_outlays"])
    housing_demand = float(np.asarray(current.demand_by_loc).reshape(-1)[0])
    scaled_supply = float(np.asarray(current.supply_by_loc).reshape(-1)[0])
    housing_gap = _relative_gap(housing_demand, scaled_supply)
    diagnostics, household_gates = _household_checks(
        current, P, shared, grid, float(P.user_cost_rate) * q, primitive, audit
    )
    reconstruction_gates = {
        "stationary_post_fertility_nesting_l1": (
            abs(float(reconstruction["stationary_post_fertility_nesting_l1"]))
            <= float(audit.reconstruction_tolerance)
        ),
        "stationary_post_fertility_nesting_max_abs": (
            abs(float(reconstruction["stationary_post_fertility_nesting_max_abs"]))
            <= float(audit.reconstruction_tolerance)
        ),
        "stationary_feasibility_projection": (
            0.0 <= float(reconstruction["stationary_feasibility_projection_mass"])
            <= float(audit.feasibility_projection_tolerance)
        ),
    }
    residual = RESIDUAL_SCALE * np.array(
        [renewal_ratio - 1.0, pension_gap, rebate_gap], dtype=float
    )
    mapping_valid = bool(
        np.isfinite(residual).all()
        and np.isfinite([population_scale, housing_gap]).all()
        and population_scale > 0
        and abs(housing_gap) <= UNSCALED_ROOT_TOLERANCE
        and all(household_gates.values())
        and all(reconstruction_gates.values())
    )
    raw_queue_flow = raw_births / REPLACEMENT_FERTILITY
    state = pf.PFInitialState(
        g_pre=scaled_g,
        scheduled_entries=[entry_flow] * QUEUE_WAITING_SLOTS,
        scheduled_raw_entries=[raw_queue_flow] * QUEUE_WAITING_SLOTS,
    )
    owner_rate = float(np.sum(current.g_current[:, 1:, :, :, :, :, :])) / head_mass
    return SimpleNamespace(
        trial=trial, parameters=P, b_grid=grid, policy=policy,
        asset_price=q, renter_price=float(P.user_cost_rate) * q,
        pension_period=pension, transfer=transfer,
        population_scale=population_scale, unit_mass=unit_mass,
        unit_g_pre=unit_g, state=state, reconstruction=dict(reconstruction),
        adjusted_births=adjusted_births, raw_births=raw_births,
        entry_flow=entry_flow, raw_queue_flow=raw_queue_flow,
        renewal_ratio=renewal_ratio, housing_demand=housing_demand,
        housing_supply=scaled_supply, housing_relative_gap=housing_gap,
        owner_rate=owner_rate, payroll=dict(payroll), rebate=dict(tax),
        residual=residual, pension_relative_gap=pension_gap,
        rebate_relative_gap=rebate_gap, diagnostics=diagnostics,
        household_gates=dict(household_gates),
        reconstruction_gates=reconstruction_gates, mapping_valid=mapping_valid,
        feasibility_projection_tolerance=float(audit.feasibility_projection_tolerance),
    )


def _trial_payload(candidate: Any) -> dict[str, Any]:
    return {
        "trial": candidate.trial,
        "coordinates": [candidate.asset_price, candidate.pension_period, candidate.transfer],
        "residual": candidate.residual,
        "mapping_valid": candidate.mapping_valid,
        "population_scale": candidate.population_scale,
        "unit_stationary_mass": candidate.unit_mass,
        "renewal_ratio": candidate.renewal_ratio,
        "housing_relative_gap": candidate.housing_relative_gap,
        "pension_relative_gap": candidate.pension_relative_gap,
        "rebate_relative_gap": candidate.rebate_relative_gap,
        "payroll_accounts": candidate.payroll,
        "rebate_accounts": candidate.rebate,
        "household_gates": candidate.household_gates,
        "reconstruction_gates": candidate.reconstruction_gates,
    }


def _one_step_audit(candidate: Any, old: Any, psi: float) -> dict[str, Any]:
    _, _, _, pf, _ = _runtime()
    state = candidate.state
    replay = pf.evaluate_path_at_prices(
        prices=np.array([candidate.asset_price]),
        psi_path=np.array([float(psi)]),
        transfer_path=np.array([candidate.transfer]),
        terminal_price=candidate.asset_price,
        terminal_V=candidate.policy.V,
        base_parameters=candidate.parameters,
        b_grid=candidate.b_grid,
        initial_state=state,
        supply_rule=old.supply_rule,
        birth_to_entry_conversion=1.0 / REPLACEMENT_FERTILITY,
        pension_path=np.array([candidate.pension_period]),
        payroll_tax_path=np.array([PAYROLL_TAX]),
    )
    next_state = replay.terminal_state
    next_g = np.asarray(next_state.g_pre, dtype=float)
    original_g = np.asarray(state.g_pre, dtype=float)
    original_mass = float(original_g.sum())
    next_mass = float(next_g.sum())
    level_l1 = float(np.sum(np.abs(next_g - original_g)))
    normalized_l1 = float(np.sum(np.abs(
        next_g / max(next_mass, 1e-15) - original_g / max(original_mass, 1e-15)
    )))
    entry_queue_gap = float(np.max(np.abs(
        np.asarray(next_state.scheduled_entries, dtype=float)
        - np.asarray(state.scheduled_entries, dtype=float)
    )))
    raw_queue_gap = float(np.max(np.abs(
        np.asarray(next_state.scheduled_raw_entries, dtype=float)
        - np.asarray(state.scheduled_raw_entries, dtype=float)
    )))
    checks = {
        "root_renewal": abs(candidate.renewal_ratio - 1.0) <= UNSCALED_ROOT_TOLERANCE,
        "root_paygo": abs(candidate.pension_relative_gap) <= UNSCALED_ROOT_TOLERANCE,
        "root_equal_rebate": abs(candidate.rebate_relative_gap) <= UNSCALED_ROOT_TOLERANCE,
        "housing_level": abs(candidate.housing_demand - candidate.housing_supply)
            <= 2e-10 * max(1.0, abs(candidate.housing_supply)),
        "one_step_population_level": abs(next_mass - original_mass)
            <= 2e-8 * max(1.0, original_mass),
        "one_step_distribution_level_l1": level_l1 <= 2e-8 * max(1.0, original_mass),
        "one_step_distribution_normalized_l1": normalized_l1 <= 2e-8,
        "one_step_entry_queue": entry_queue_gap <= 2e-8 * max(1.0, candidate.entry_flow),
        "one_step_raw_queue": raw_queue_gap <= 2e-8 * max(1.0, candidate.raw_queue_flow),
        "one_step_policy_reproduction": replay.maximum_policy_reproduction_error <= 2e-8,
        "one_step_mass_accounting": replay.maximum_mass_accounting_error <= 2e-8,
        "one_step_feasibility_projection": replay.maximum_feasibility_projection_mass
            <= float(getattr(candidate, "feasibility_projection_tolerance", 1e-6)),
        "household_audit": all(candidate.household_gates.values()),
        "stationary_reconstruction": all(candidate.reconstruction_gates.values()),
        "four_waiting_entry_vintages": len(state.scheduled_entries) == QUEUE_WAITING_SLOTS,
        "four_waiting_raw_vintages": len(state.scheduled_raw_entries) == QUEUE_WAITING_SLOTS,
    }
    return {
        "status": "passed" if all(checks.values()) else "failed",
        "checks": checks,
        "initial_population_mass": original_mass,
        "next_population_mass": next_mass,
        "population_absolute_gap": abs(next_mass - original_mass),
        "distribution_level_l1": level_l1,
        "distribution_normalized_l1": normalized_l1,
        "entry_queue_maximum_absolute_gap": entry_queue_gap,
        "raw_queue_maximum_absolute_gap": raw_queue_gap,
        "maximum_policy_reproduction_error": replay.maximum_policy_reproduction_error,
        "maximum_mass_accounting_error": replay.maximum_mass_accounting_error,
        "maximum_feasibility_projection_mass": replay.maximum_feasibility_projection_mass,
    }


def solve_terminal(*, old: Any, psi: float, audit: Any,
                   controls: dict[str, Any], deadline: float, folder: Path,
                   start: Any = None) -> SimpleNamespace:
    """Solve and certify the original-law stationary terminal endpoint.

    Failure is returned as ``verified=False`` with a tensor-free receipt.  When
    available, the best admissible trial is returned as an explicitly
    diagnostic policy/state, but no ``terminal.pkl.gz`` is written.
    """
    from e5f_balanced_terminal import TerminalAuditControls
    from e5f_matched_pf_path_root import solve_price_path

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    if not isinstance(audit, TerminalAuditControls):
        raise ValueError("Explicit TerminalAuditControls are required")
    audit_ceilings = {
        "reconstruction_tolerance": 5e-9,
        "feasibility_projection_tolerance": 1e-6,
        "probability_tolerance": 1e-12,
        "occupied_mass_tolerance": 1e-12,
        "value_drop_tolerance": 1e-7,
    }
    for name, ceiling in audit_ceilings.items():
        value = float(getattr(audit, name))
        if not math.isfinite(value) or not 0 <= value <= ceiling:
            raise ValueError(f"Invalid retained terminal household audit: {name}")
    if not np.isfinite(psi):
        raise ValueError("Terminal preference must be finite")
    if not np.isfinite(deadline) or deadline <= time.monotonic():
        raise ValueError("Terminal deadline has expired")
    _validate_supply(old)
    bounds, initial, root_controls, initial_jacobian = _validated_controls(
        controls, start, old
    )
    latest = None
    best_candidate = None
    trial = 0
    started = time.monotonic()

    def project(raw):
        values = np.asarray(raw, dtype=float).copy()
        if values.shape != (3,):
            raise ValueError("Terminal root must retain exactly three coordinates")
        for index, (lower, upper) in enumerate(bounds):
            values[index] = np.clip(values[index], lower, upper)
        return values

    def evaluate(raw):
        nonlocal latest, trial
        if time.monotonic() >= deadline:
            raise TimeoutError("Terminal deadline reached between model evaluations")
        trial += 1
        latest = _evaluate_trial(
            old=old, psi=float(psi), coordinates=np.asarray(raw, dtype=float),
            audit=audit, deadline=deadline, trial=trial,
        )
        _save_json(folder / "latest_phase.json", {
            "status": "completed_mapping", "trial": trial,
            "elapsed_seconds": time.monotonic() - started,
            **_trial_payload(latest),
        })
        return {
            "residual": latest.residual,
            "mapping_valid": bool(latest.mapping_valid),
            "payload": _trial_payload(latest),
        }

    def progress(record):
        nonlocal best_candidate
        enriched = dict(record, coordinate_labels=["asset_price", "pension", "rebate"])
        if record.get("evaluation") is not None:
            _save_json(folder / "latest_completed.json", enriched)
        elif record.get("event") == "complete":
            _save_json(folder / "root_status.json", enriched)
        if record.get("new_best"):
            best_candidate = latest
            _save_json(folder / "best_so_far.json", enriched)

    default_jacobian = np.diag([-RESIDUAL_SCALE * float(root_controls["slope"]),
                                -RESIDUAL_SCALE, -RESIDUAL_SCALE])
    root_arguments = dict(
        initial_prices=initial, evaluate=evaluate, project=project,
        deadline_monotonic=float(deadline), callback=progress,
        default_jacobian=default_jacobian, initial_jacobian=initial_jacobian,
        **root_controls,
    )
    try:
        root = solve_price_path(**root_arguments)
    except (RuntimeError, ValueError, TimeoutError, FloatingPointError) as exc:
        receipt = {
            "schema": "e5f_original_queue_terminal_v1",
            "status": "controlled_failure",
            "verified": False,
            "error_type": type(exc).__name__, "error": str(exc),
            "trials_completed": trial,
            "elapsed_seconds": time.monotonic() - started,
            "coordinate_labels": ["asset_price", "pension", "rebate"],
            "root_residual_definition": [
                "200*(topcode_adjusted_births/(2.1*age0_entries)-1)",
                "200*PAYGO_relative_gap", "200*equal_rebate_relative_gap",
            ],
            "best_diagnostic_available": best_candidate is not None,
            "production_eligible": False,
        }
        _save_json(folder / "root_receipt.json", receipt)
        chosen = best_candidate
        return SimpleNamespace(
            parameters=None if chosen is None else chosen.parameters,
            policy=None if chosen is None else chosen.policy,
            asset_price=None if chosen is None else chosen.asset_price,
            state=None if chosen is None else chosen.state,
            coordinates=None if chosen is None else np.array([
                chosen.asset_price, chosen.pension_period, chosen.transfer]),
            receipt=receipt, verified=False,
        )

    final = root.get("final")
    matched_fresh = bool(
        root.get("converged") and final is not None and latest is not None
        and final.get("payload", {}).get("trial") == latest.trial
        and np.array_equal(np.asarray(final["prices"], dtype=float), np.array([
            latest.asset_price, latest.pension_period, latest.transfer]))
    )
    one_step = None
    one_step_error = None
    if matched_fresh and latest.mapping_valid:
        try:
            one_step = _one_step_audit(latest, old, float(psi))
        except (RuntimeError, ValueError, TimeoutError, FloatingPointError) as exc:
            one_step_error = {"error_type": type(exc).__name__, "error": str(exc)}
    verified = bool(
        matched_fresh and latest.mapping_valid and one_step is not None
        and one_step["status"] == "passed"
        and float(root.get("final_reproduction_max_abs", math.inf)) <= 2e-10
    )
    chosen = latest if matched_fresh else best_candidate
    reference = None if chosen is None else _candidate_reference(chosen)
    receipt = _jsonable(root)
    receipt.update(
        schema="e5f_original_queue_terminal_v1",
        verified=verified,
        production_eligible=False,
        stationary_endpoint_verified=verified,
        fresh_final_mapping_matches_endpoint=matched_fresh,
        coordinate_labels=["asset_price", "pension", "rebate"],
        root_residual_definition=[
            "200*(topcode_adjusted_births/(2.1*age0_entries)-1)",
            "200*PAYGO_relative_gap", "200*equal_rebate_relative_gap",
        ],
        root_tolerance=float(root_controls["market_tolerance"]),
        unscaled_root_tolerance=UNSCALED_ROOT_TOLERANCE,
        population_scale_rule="unchanged_supply_quantity(asset_price)/unit_household_demand",
        initial_distribution_rule="stationary original household law; no historical conditioning",
        queue_rule="four waiting vintages; adjusted entry flow and raw births both divided by 2.1",
        one_step_audit=one_step,
        one_step_error=one_step_error,
        endpoint_reference=reference,
        elapsed_seconds_total=time.monotonic() - started,
    )
    _save_json(folder / "root_receipt.json", receipt)
    result = SimpleNamespace(
        parameters=None if chosen is None else chosen.parameters,
        policy=None if chosen is None else chosen.policy,
        asset_price=None if chosen is None else chosen.asset_price,
        state=None if chosen is None else chosen.state,
        coordinates=None if chosen is None else np.array([
            chosen.asset_price, chosen.pension_period, chosen.transfer]),
        receipt=receipt, verified=verified,
        endpoint_reference=reference,
    )
    if verified:
        _save_pickle(folder / "terminal.pkl.gz", result)
    return result
