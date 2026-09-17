"""Bounded adapter for the original household birth-vintage queue.

This experiment keeps the retained calibration, fiscal checks, household
kernel, grid, and housing-supply rule.  It removes the historical age bridge,
the 2023 person-state reset, outside migration, and all person-demography
updates from the state law actually evaluated.  Entrants instead follow the
original four-vintage household queue for every date.

The adapter is intentionally installed only by :func:`original_queue_adapter`.
It does not change the frozen helper or model source on disk.
"""
from __future__ import annotations

import copy
from contextlib import contextmanager
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np


PAYROLL_TAX = 0.179
BIRTH_TO_ENTRY_CONVERSION = 1.0 / 2.1
QUEUE_SLOTS = 4


def annotate_original_queue_metadata(metadata=None):
    """Return metadata that names the experimental population law exactly."""
    result = dict(metadata or {})
    result.update(
        population_law="original_household_birth_vintage_queue",
        population_closure="original_queue",
        original_queue=True,
        uses_person_demography=False,
        person_demography_case=None,
        person_A0=False,
        demographic_case="original_queue",
        historical_age_reweighting_in_realized_state=False,
        outside_migration=False,
        production_eligible=False,
    )
    return result


def initialize_original(old, packet):
    """Copy an approved old state and install its stationary household state.

    The old builder is still responsible for all parameter, normalization,
    pension, fiscal, supply, grid, and kernel checks.  Only the realized initial
    distribution and its queues change here.  The adjusted queue is initialized
    at the actual age-zero household mass, while the raw queue is initialized at
    raw births divided by 2.1.
    """
    if not hasattr(old, "initial_state"):
        raise ValueError("Approved old state must contain an initial_state")
    packet_pre = np.asarray(packet["stationary_g_pre"], dtype=float)
    retained_pre = np.asarray(getattr(old, "stationary_g_pre", packet_pre), dtype=float)
    if (packet_pre.shape != retained_pre.shape
            or not np.array_equal(packet_pre, retained_pre)):
        raise ValueError("Approved and packet stationary household states differ")
    stationary = retained_pre.copy()
    if (stationary.ndim != 7 or not np.isfinite(stationary).all()
            or float(stationary.min()) < -1e-13):
        raise ValueError("Stationary pre-fertility distribution is invalid")

    entry_flow = float(np.sum(stationary[:, :, :, 0, :, :, :]))
    solution = packet.get("solution", getattr(old, "solution", None))
    raw_births = float(getattr(solution, "total_births_kfe", np.nan))
    raw_entry_flow = raw_births / 2.1
    if (not np.isfinite([entry_flow, raw_births, raw_entry_flow]).all()
            or entry_flow <= 0.0 or raw_births < 0.0):
        raise ValueError("Finite positive actual entry mass and raw births required")

    original = old.initial_state
    if (len(original.scheduled_entries) != QUEUE_SLOTS
            or len(original.scheduled_raw_entries) != QUEUE_SLOTS):
        raise ValueError("Original household law requires four waiting vintages")
    state_type = type(original)
    initial = state_type(
        g_pre=stationary.copy(),
        scheduled_entries=[entry_flow] * QUEUE_SLOTS,
        scheduled_raw_entries=[raw_entry_flow] * QUEUE_SLOTS,
    )

    result = copy.copy(old)
    result.stationary_g_pre = stationary.copy()
    result.initial_state = initial
    result.diagnostics = annotate_original_queue_metadata(
        getattr(old, "diagnostics", {}))
    result.diagnostics.update(
        original_queue_initial_age_zero_mass=entry_flow,
        original_queue_initial_raw_births=raw_births,
        original_queue_initial_raw_entry_flow=raw_entry_flow,
        original_queue_birth_to_entry_conversion=BIRTH_TO_ENTRY_CONVERSION,
        original_queue_waiting_slots=QUEUE_SLOTS,
    )
    if (result.parameters is not old.parameters
            or result.b_grid is not old.b_grid
            or result.supply_rule is not old.supply_rule):
        raise RuntimeError("Original-queue initialization changed retained model objects")
    return result


def _validated_paths(prices, pensions, transfers):
    p = np.asarray(prices, dtype=float)
    b = np.asarray(pensions, dtype=float)
    t = np.asarray(transfers, dtype=float)
    if p.ndim != 1 or len(p) < 1 or b.shape != p.shape or t.shape != p.shape:
        raise ValueError("Equally sized one-dimensional dated paths required")
    if (not np.isfinite(p).all() or np.any(p <= 0.0)
            or not np.isfinite(b).all() or np.any(b < 0.0)
            or not np.isfinite(t).all() or np.any(t < 0.0)):
        raise ValueError("Finite positive prices and nonnegative fiscal paths required")
    return p, b, t


def queue_path(*, inherited, old_state, prices, pensions, transfers, psi,
               terminal, observer=None, demographics=None,
               demographic_evaluator=None):
    """Evaluate the original household queue over the entire supplied path.

    ``demographics`` and ``demographic_evaluator`` are accepted only for frozen
    caller compatibility.  Neither enters the state transition.
    """
    del demographics, demographic_evaluator
    import e5f_rebated_surprises as rebated

    _, joined, _, _, _ = rebated._runtime()
    pf = joined.pf
    P, grid = old_state.parameters, old_state.b_grid
    terminal_parameters, terminal_policy, terminal_price, _ = (
        rebated._terminal_parts(terminal))
    if (not hasattr(terminal_parameters, "psi_child")
            or not np.isclose(float(terminal_parameters.psi_child), float(psi),
                              rtol=0.0, atol=1e-14)):
        raise ValueError("Forecast terminal must use the current preference")
    if (not np.isfinite(psi) or not isinstance(inherited.year, (int, np.integer))
            or (int(inherited.year) - 2007) % 4 != 0):
        raise ValueError("Finite preference and a four-year model date required")
    p, benefits, rebates = _validated_paths(prices, pensions, transfers)

    native_evaluate_period = pf.calendar.evaluate_period
    observed = 0

    def dated_evaluate_period(*args, **kwargs):
        nonlocal observed
        evaluation = native_evaluate_period(*args, **kwargs)
        if observer is not None:
            # Native signature: price, g_pre, P, b_grid, shared, counter, ...
            observer(observed, evaluation, args[2], args[3], args[4])
        observed += 1
        return evaluation

    with patch.object(pf.calendar, "evaluate_period", dated_evaluate_period):
        native = pf.evaluate_path_at_prices(
            prices=p,
            psi_path=np.full(len(p), float(psi)),
            transfer_path=rebates,
            terminal_price=float(terminal_price),
            terminal_V=terminal_policy.V,
            base_parameters=P,
            b_grid=grid,
            initial_state=inherited.households,
            supply_rule=old_state.supply_rule,
            birth_to_entry_conversion=BIRTH_TO_ENTRY_CONVERSION,
            historical_conditioning=None,
            pension_path=benefits,
            payroll_tax_path=np.full(len(p), PAYROLL_TAX),
        )
    if observed != len(p):
        raise RuntimeError("Original queue did not evaluate every dated period exactly once")
    if (native.maximum_mass_accounting_error > 2e-8
            or native.maximum_policy_reproduction_error > 2e-10
            or native.maximum_feasibility_projection_mass > 1e-6):
        raise RuntimeError("Full original-queue path failed mass, policy replay, or feasibility gates")

    for i, row in enumerate(native.rows):
        row["period"] = i
        row["calendar_year"] = int(inherited.year) + 4 * i
        row["annual_net_migration_over_period"] = 0.0
        row["net_migrant_heads_over_period"] = 0.0
        row["pension_period"] = float(benefits[i])
        row["pension_period_units"] = float(benefits[i])
        row.update(annotate_original_queue_metadata())

    history = SimpleNamespace(
        rows=[], values=[], bellman_solves=0,
        maximum_mass_accounting_error=0.0,
        maximum_policy_reproduction_error=0.0,
        maximum_feasibility_projection_mass=0.0,
    )
    return SimpleNamespace(
        history=history,
        person_tail=native,
        rows=native.rows,
        values=native.values,
        bellman_solves=native.bellman_solves,
        maximum_market_residual=native.maximum_market_residual,
        maximum_mass_accounting_error=native.maximum_mass_accounting_error,
        maximum_policy_reproduction_error=native.maximum_policy_reproduction_error,
        maximum_feasibility_projection_mass=native.maximum_feasibility_projection_mass,
        elapsed_seconds=native.elapsed_seconds,
        pension_period=benefits.copy(),
        initial_2023_age_head_gap=0.0,
        metadata=annotate_original_queue_metadata(),
    )


def first_period_state(*, inherited, old_state, demographics, path, prices,
                       pensions, transfers, psi, demographic_evaluator=None):
    """Replay the first accepted period under its own next price and value."""
    import e5f_rebated_surprises as rebated

    p, benefits, rebates = _validated_paths(prices, pensions, transfers)
    if len(p) < 2 or len(path.values) < 2 or not path.rows:
        raise ValueError("First-period replay needs two prices and path values")
    terminal_parameters = SimpleNamespace(psi_child=float(psi))
    replay_terminal = SimpleNamespace(
        parameters=terminal_parameters,
        policy=SimpleNamespace(V=np.asarray(path.values[1], dtype=float)),
        asset_price=float(p[1]),
    )
    replay = queue_path(
        inherited=inherited,
        old_state=old_state,
        prices=p[:1],
        pensions=benefits[:1],
        transfers=rebates[:1],
        psi=psi,
        terminal=replay_terminal,
        observer=None,
        demographics=demographics,
        demographic_evaluator=demographic_evaluator,
    )
    native = replay.person_tail
    if (native.maximum_mass_accounting_error > 2e-8
            or native.maximum_policy_reproduction_error > 2e-10
            or native.maximum_feasibility_projection_mass > 1e-6):
        raise RuntimeError("First-period original-queue replay gate failed")
    if not np.allclose(replay.values[0], path.values[0], rtol=0.0, atol=2e-10):
        raise RuntimeError("First-period continuation value differs from accepted forecast")
    keys = (
        "asset_price", "renter_price", "housing_demand", "housing_supply",
        "owner_rate", "birth_children_topcode_adjusted",
        "pension_period_units", "payroll_tax_revenue", "pension_outlays",
        "property_tax_revenue", "equal_transfer_outlays",
        "effective_mature_entrant_flow_B", "raw_state_scheduled_mature_entrant_flow_B",
        "entrant_flow_next", "mass_accounting_residual",
    )
    for key in keys:
        expected, actual = path.rows[0].get(key), replay.rows[0].get(key)
        if ((expected is None) != (actual is None)
                or (expected is not None and not np.isclose(
                    float(actual), float(expected), rtol=0.0, atol=2e-10))):
            raise RuntimeError("First-period original-queue replay differs: " + key)
    return rebated.InheritedState(
        int(inherited.year) + 4, native.terminal_state)


@contextmanager
def original_queue_adapter():
    """Temporarily install the original queue in the frozen rebated workflow."""
    import e5f_rebated_initial_bridge as bridge
    import e5f_rebated_surprises as rebated

    native_builder = bridge.build_rebated_initial_state

    def build_original_initial_state(**kwargs):
        old = native_builder(**kwargs)
        return initialize_original(old, kwargs["packet"])

    with patch.object(bridge, "build_rebated_initial_state",
                      build_original_initial_state), patch.object(
                          rebated, "evaluate_forecast", queue_path), patch.object(
                              rebated, "first_period_state", first_period_state):
        yield SimpleNamespace(
            bridge=bridge,
            rebated=rebated,
            metadata=annotate_original_queue_metadata(),
        )
