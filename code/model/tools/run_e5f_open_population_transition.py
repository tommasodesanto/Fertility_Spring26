#!/usr/bin/env python3
"""Cohort-accounting transition for the sequential child-room-floor model.

The script makes no change to the household state space or stationary solver.
It starts from a stationary sequential-model distribution, applies a permanent
change in the child-preference shifter, and carries the full age, wealth,
income, tenure, parity, and children-at-home distribution through calendar
time. At each date, households treat the current price and primitives as
permanent and the contemporaneous housing market clears against the initial
housing stock. This is therefore a transparent temporary-equilibrium
transition diagnostic, not a perfect-foresight asset-pricing exercise.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import shlex
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for path in (MODEL_ROOT, TOOLS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import audit_closed_reproductive_closure as closure_audit
import run_dynamic_population_transition as calendar


DEFAULT_SOURCE = closure_audit.E5F_FLOOR_SOURCE
DEFAULT_OUTDIR = ROOT / "output/model/e5f_floor_open_population_transition"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--periods", type=int, default=25)
    parser.add_argument(
        "--old-psi-child",
        type=float,
        default=0.10,
        help="Child-preference shifter in the initial stationary equilibrium.",
    )
    parser.add_argument(
        "--new-psi-child",
        type=float,
        default=0.0,
        help="Permanent post-shock child-preference shifter.",
    )
    parser.add_argument(
        "--outside-origin-entry-share",
        type=float,
        default=0.169,
        help="Initial steady-state share of entrant households supplied from outside.",
    )
    parser.add_argument(
        "--market-tol",
        type=float,
        default=2e-4,
        help=(
            "Relative clearing tolerance. The transition distribution is atomic "
            "on the deterministic tenure/wealth grid, so its demand schedule can "
            "jump across the exact root."
        ),
    )
    parser.add_argument("--market-max-iter", type=int, default=30)
    parser.add_argument(
        "--renewal-clock",
        choices=("birth-vintage", "household-state"),
        default="birth-vintage",
        help=(
            "Birth-vintage is the paper transition: births wait four complete "
            "model periods before supplying entrant households. Household-state "
            "reuses the model's memoryless child-exit flow as a diagnostic."
        ),
    )
    parser.add_argument(
        "--birth-to-entry-delay-periods",
        type=int,
        default=4,
        help="Completed four-year child stages between birth and adult entry.",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run at most four calendar periods after the exact baseline gates.",
    )
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def configure_sequential_model() -> Any:
    chain = closure_audit.load_chain(profile="e5f-floor")
    from intergen_eqscale_seq_optimized import solver as sequential_model

    # Reuse the already-tested calendar-time scaffolding, but route every
    # household and transition operation through the sequential model.
    calendar.model = sequential_model
    return chain, sequential_model


def apply_sequential_fertility(
    g_pre: np.ndarray,
    fert_probs: np.ndarray,
    P: SimpleNamespace,
) -> tuple[np.ndarray, float, np.ndarray]:
    """Apply at most one sequential birth per household during the period."""
    model = calendar.model
    if not bool(getattr(P, "sequential_births", False)):
        raise RuntimeError("The E5F transition requires sequential_births=True")
    if not model.independent_child_maturation_active(P):
        raise RuntimeError("The E5F transition requires independent child maturation")

    out = np.asarray(g_pre, dtype=float).copy()
    fecundity = model.get_fecundity_by_age(P)
    continuation = getattr(P, "_fert2_probs", None)
    if continuation is None:
        raise RuntimeError("Sequential continuation-birth probabilities are unavailable")
    births = 0.0
    births_by_loc = np.zeros(int(P.I))
    for j in range(int(P.J)):
        if not (int(P.A_f_start) <= j + 1 <= int(P.A_f_end)):
            continue
        pi_j = float(fecundity[j])
        for zz in range(g_pre.shape[4]):
            # Take every continuation at-risk pool before any birth lands, as
            # in the stationary KFE. This prevents two births in one period.
            at_risk = {
                (parity, child_state): g_pre[
                    :, :, :, j, zz, parity, child_state
                ].copy()
                for parity in range(1, int(P.n_parity) - 1)
                for child_state in range(parity + 1)
            }

            settled = model.readiness_settled_state(P)
            childless = g_pre[:, :, :, j, zz, 0, settled]
            try_mass = childless * fert_probs[:, :, :, j, zz, 1]
            realized = pi_j * try_mass
            out[:, :, :, j, zz, 0, settled] -= realized
            out[:, :, :, j, zz, 1, 1] += realized
            births += float(np.sum(realized))
            births_by_loc += np.sum(realized, axis=(0, 1))

            for parity in range(1, int(P.n_parity) - 1):
                for child_state in range(parity + 1):
                    pool = at_risk[(parity, child_state)]
                    attempt_probability = continuation[
                        :, :, :, j, zz, 1, parity - 1, child_state
                    ]
                    realized = pi_j * pool * attempt_probability
                    destination = model.birth_destination_child_state(P, child_state)
                    out[:, :, :, j, zz, parity, child_state] -= realized
                    out[:, :, :, j, zz, parity + 1, destination] += realized
                    births += float(np.sum(realized))
                    births_by_loc += np.sum(realized, axis=(0, 1))

    minimum = float(np.min(out))
    if minimum < -1e-13:
        raise RuntimeError(f"Sequential fertility produced negative mass: {minimum:.3e}")
    out[out < 0.0] = 0.0
    mass_error = float(np.sum(out) - np.sum(g_pre))
    if abs(mass_error) > 2e-10:
        raise RuntimeError(f"Sequential fertility failed mass conservation: {mass_error:.3e}")
    return out, births, births_by_loc


def children_at_home_units(distribution: np.ndarray, P: SimpleNamespace) -> float:
    total = 0.0
    for parity in range(1, int(P.n_parity)):
        for child_state in range(parity + 1):
            total += child_state * float(
                np.sum(distribution[:, :, :, :, parity, child_state])
            )
    return total


def advance_sequential_calendar_distribution(
    evaluation: calendar.PeriodEvaluation,
    entry_by_loc_next: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Advance incumbents one age cell and measure independently matured children."""
    model = calendar.model
    policy = evaluation.policy
    g_post = evaluation.g_post_fertility
    next_pre = np.zeros_like(g_post)
    mature_by_loc = np.zeros(int(P.I))
    deaths = 0.0
    _, _, Pi_z = model.income_transition_values(P)
    stochastic = bool(P.use_stochastic_aging and hasattr(P, "Pi_child"))
    if int(P.I) != 1:
        raise NotImplementedError("The sequential transition is currently one-market only")

    for j in range(int(P.J) - 1):
        survival = float(P.survival_probs[j]) if bool(P.use_age_survival) else 1.0
        cohort = g_post[:, :, :, j, :, :, :]
        deaths += (1.0 - survival) * float(np.sum(cohort))
        survivors = survival * cohort
        advanced = model.advance_cohort_one_period_markov_income(
            survivors,
            j,
            policy.loc_probs,
            policy.tenure_choice,
            policy.tenure_probs,
            policy.bp_pol,
            P,
            b_grid,
            shared,
            policy.maps.lmm_idx,
            policy.maps.lmm_wt,
            policy.maps.tmx_idx,
            policy.maps.tmx_wt,
            stochastic,
            P.Pi_child if stochastic else None,
            Pi_z,
        )
        next_pre[:, :, :, j + 1, :, :, :] = advanced
        matured_children = children_at_home_units(survivors, P) - children_at_home_units(
            advanced, P
        )
        if matured_children < -2e-10:
            raise RuntimeError(
                f"Children-at-home units rose during maturation at age index {j}: "
                f"{matured_children:.3e}"
            )
        mature_by_loc[0] += max(matured_children, 0.0)

    deaths += float(np.sum(g_post[:, :, :, int(P.J) - 1, :, :, :]))
    next_pre[:, :, :, 0, :, :, :] = calendar.entrant_cohort(
        np.asarray(entry_by_loc_next, dtype=float), P, b_grid
    )
    expected_mass = float(np.sum(g_post)) - deaths + float(np.sum(entry_by_loc_next))
    mass_residual = float(np.sum(next_pre)) - expected_mass
    return next_pre, mature_by_loc, deaths, mass_residual


def independent_child_distribution_rows(
    scenario: str,
    period: int,
    g_post: np.ndarray,
    P: SimpleNamespace,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    total = float(np.sum(g_post))
    age_rows: list[dict[str, Any]] = []
    for j in range(int(P.J)):
        mass = float(np.sum(g_post[:, :, :, j, :, :, :]))
        age_rows.append(
            {
                "scenario": scenario,
                "period": period,
                "years_from_start": period * float(P.period_years),
                "age": float(P.age_start + j * P.da),
                "adult_mass": mass,
                "adult_share": mass / max(total, 1e-15),
            }
        )
    child_rows: list[dict[str, Any]] = []
    for parity in range(int(P.n_parity)):
        for child_state in range(int(P.n_child_states)):
            mass = float(np.sum(g_post[:, :, :, :, :, parity, child_state]))
            valid_child_count = child_state if parity >= 1 and child_state <= parity else 0
            child_rows.append(
                {
                    "scenario": scenario,
                    "period": period,
                    "years_from_start": period * float(P.period_years),
                    "parity_state": parity,
                    "child_state": child_state,
                    "child_state_label": (
                        "no_children_at_home"
                        if valid_child_count == 0
                        else f"active_child_stage_{valid_child_count}"
                    ),
                    "household_mass": mass,
                    "household_share": mass / max(total, 1e-15),
                    "child_units": valid_child_count * mass,
                }
            )
    return age_rows, child_rows


def operator_gates(
    solution: SimpleNamespace,
    policy: calendar.PolicyBundle,
    g_pre: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> dict[str, Any]:
    P._fert2_probs = np.asarray(solution.fert2_probs, dtype=float).copy()
    evaluation = calendar.evaluate_period(
        policy.price,
        g_pre,
        P,
        b_grid,
        shared,
        calendar.SolveCounter(),
        supplied_policy=policy,
    )
    empty_next, mature_raw_by_loc, _, zero_entry_mass_residual = (
        advance_sequential_calendar_distribution(
            evaluation, np.zeros(int(P.I)), P, b_grid, shared
        )
    )
    empty_next[:, :, :, 0, :, :, :] = calendar.entrant_cohort(
        np.asarray(solution.entry_by_loc, dtype=float), P, b_grid
    )
    effective_mature = float(P.entrant_conversion_factor) * float(
        np.sum(mature_raw_by_loc)
    )
    return {
        "stationary_post_fertility_nesting_l1": float(
            np.sum(
                np.abs(
                    evaluation.g_post_fertility
                    - np.asarray(solution.g_beginning_distribution, dtype=float)
                )
            )
        ),
        "one_step_constant_path_nesting_l1": float(np.sum(np.abs(empty_next - g_pre))),
        "zero_entry_mass_accounting_residual": float(zero_entry_mass_residual),
        "entry_flow_E": float(solution.entry_rate),
        "effective_mature_flow_B": float(solution.entrants_mature_total),
        "reconstructed_effective_mature_flow_B": effective_mature,
        "mature_flow_abs_error": abs(effective_mature - float(solution.entrants_mature_total)),
        "birth_flow": float(evaluation.births),
        "stationary_birth_flow": float(solution.total_births_kfe),
        "birth_flow_abs_error": abs(float(evaluation.births) - float(solution.total_births_kfe)),
    }


def run_birth_vintage_scenario(
    *,
    label: str,
    initial_g_pre: np.ndarray,
    baseline_price: float,
    baseline_B: float,
    baseline_births: float,
    parameters: SimpleNamespace,
    b_grid: np.ndarray,
    periods: int,
    outside_flow: float,
    retention: float,
    conversion: float,
    delay_periods: int,
    supply_rule: calendar.HousingSupplyRule,
    market_tol: float,
    market_max_iter: int,
    outdir: Path,
    counter: calendar.SolveCounter,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Run the open transition with a dated birth-to-entry queue.

    The household child state continues to govern housing demand exactly as in
    the calibrated model. Population renewal instead uses a separate aggregate
    birth-vintage queue, so a newborn cannot become an entrant household in the
    same four-year period.
    """
    if delay_periods < 1:
        raise ValueError("The birth-to-entry delay must be at least one period")
    P = copy.deepcopy(parameters)
    shared = calendar.model.precompute_shared(P, b_grid)
    g_pre = np.asarray(initial_g_pre, dtype=float).copy()
    outside_shares = np.asarray(P.entry_shares, dtype=float).reshape(-1)
    outside_shares = outside_shares / np.sum(outside_shares)
    if int(P.I) != 1:
        raise NotImplementedError("The birth-vintage transition is currently one-market only")

    # At the old steady state this yield maps births exactly into the model's
    # effective mature-entrant flow. It combines child survival and the
    # children-to-households conversion, while the queue supplies the timing.
    maturation_survival_yield = baseline_B / max(conversion * baseline_births, 1e-15)
    if not 0.0 < maturation_survival_yield <= 1.0 + 1e-10:
        raise RuntimeError(
            "The stationary birth-to-entry yield is invalid: "
            f"{maturation_survival_yield:.12g}"
        )
    scheduled_entries = [float(baseline_B)] * int(delay_periods)

    path_rows: list[dict[str, Any]] = []
    age_rows: list[dict[str, Any]] = []
    child_rows: list[dict[str, Any]] = []
    price_guess = float(baseline_price)
    initial_policy = None
    started = time.perf_counter()
    start_solves = counter.total
    max_market_residual = 0.0
    max_mass_residual = 0.0
    max_projection = 0.0
    min_mass = math.inf
    max_nonfinite = 0
    grid_resolution_fallback_periods: list[int] = []

    for period in range(int(periods)):
        target_tolerance = min(float(market_tol), 5e-5)
        try:
            evaluation = calendar.clear_scalar_housing_market(
                g_pre,
                price_guess,
                P,
                b_grid,
                shared,
                counter,
                target_tolerance,
                min(int(market_max_iter), 18),
                supply_rule,
                initial_policy=initial_policy,
            )
        except RuntimeError as error:
            if (
                float(market_tol) <= target_tolerance
                or "Housing market did not clear" not in str(error)
            ):
                raise
            # The deterministic tenure/wealth grid can jump across the exact
            # market root. Retry only the failed date at the declared grid-
            # resolution acceptance tolerance; all smooth dates retain the
            # tighter 5e-5 target.
            evaluation = calendar.clear_scalar_housing_market(
                g_pre,
                price_guess,
                P,
                b_grid,
                shared,
                counter,
                float(market_tol),
                int(market_max_iter),
                supply_rule,
                initial_policy=initial_policy,
            )
            grid_resolution_fallback_periods.append(period)
        initial_policy = evaluation.policy
        price_guess = float(evaluation.policy.price[0])
        entry_flow = float(np.sum(evaluation.g_pre[:, :, :, 0, :, :, :]))

        empty_next, model_mature_raw_by_loc, deaths, _ = (
            advance_sequential_calendar_distribution(
                evaluation, np.zeros(int(P.I)), P, b_grid, shared
            )
        )
        scheduled_B = float(scheduled_entries.pop(0))
        new_potential_B = (
            conversion * maturation_survival_yield * float(evaluation.births)
        )
        scheduled_entries.append(new_potential_B)
        retained_B = retention * scheduled_B
        entrants_next_by_loc = outside_flow * outside_shares
        entrants_next_by_loc[0] += retained_B
        empty_next[:, :, :, 0, :, :, :] = calendar.entrant_cohort(
            entrants_next_by_loc, P, b_grid
        )
        next_pre = empty_next

        health = calendar.distribution_health(
            {
                "pre_fertility": evaluation.g_pre,
                "post_fertility": evaluation.g_post_fertility,
                "current_choice": evaluation.g_current,
                "next_pre_fertility": next_pre,
            }
        )
        period_min = health["min_distribution_mass"]
        nonfinite = int(health["nonfinite_distribution_count"])
        if period_min is None or nonfinite > 0 or float(period_min) < -1e-14:
            raise RuntimeError(
                f"{label} t={period}: distribution-health gate failed: {health}"
            )
        min_mass = min(min_mass, float(period_min))
        max_nonfinite = max(max_nonfinite, nonfinite)

        expected_next_mass = (
            float(np.sum(evaluation.g_post_fertility))
            - deaths
            + float(np.sum(entrants_next_by_loc))
        )
        mass_residual = float(np.sum(next_pre)) - expected_next_mass
        max_mass_residual = max(max_mass_residual, abs(mass_residual))
        max_market_residual = max(
            max_market_residual, float(evaluation.relative_market_residual)
        )
        max_projection = max(
            max_projection, float(evaluation.feasibility_projection_mass)
        )

        adult_population = float(np.sum(evaluation.g_post_fertility))
        current_mass = float(np.sum(evaluation.g_current))
        owner_rate = float(
            np.sum(evaluation.g_current[:, 1:, :, :, :, :, :])
        ) / max(current_mass, 1e-15)
        ages = P.age_start + np.arange(int(P.J)) * P.da
        age_mass = np.sum(
            evaluation.g_post_fertility, axis=(0, 1, 2, 4, 5, 6)
        )
        mean_age = float(np.sum(ages * age_mass) / max(np.sum(age_mass), 1e-15))
        demand = float(np.sum(evaluation.demand_by_loc))
        supply = float(np.sum(evaluation.supply_by_loc))
        model_mature_effective = conversion * float(np.sum(model_mature_raw_by_loc))
        row = {
            "scenario": label,
            "period": period,
            "years_from_start": period * float(P.period_years),
            "asset_price": price_guess,
            "asset_price_index": price_guess / baseline_price,
            "housing_user_cost": float(P.user_cost_rate * price_guess),
            "adult_population": adult_population,
            "population_index": adult_population / max(float(np.sum(initial_g_pre)), 1e-15),
            "mean_adult_age": mean_age,
            "entry_flow_E": entry_flow,
            "birth_children": float(evaluation.births),
            "births_over_entry": float(evaluation.births) / max(entry_flow, 1e-15),
            "mature_children_raw": scheduled_B / max(conversion, 1e-15),
            "effective_mature_entrant_flow_B": scheduled_B,
            "mature_to_current_entry_flow_ratio_diagnostic": scheduled_B / max(entry_flow, 1e-15),
            "model_state_same_period_mature_flow_B": model_mature_effective,
            "birth_queue_new_potential_flow_B": new_potential_B,
            "birth_queue_scheduled_flows": list(scheduled_entries),
            "closure": "open_birth_vintage",
            "housing_supply_mode": supply_rule.mode,
            "outside_entry_flow_M": outside_flow,
            "retained_mature_entrants": retained_B,
            "entrant_flow_next": float(np.sum(entrants_next_by_loc)),
            "adult_deaths": deaths,
            "housing_demand": demand,
            "housing_demand_per_adult": demand / max(adult_population, 1e-15),
            "housing_supply": supply,
            "relative_market_residual": float(evaluation.relative_market_residual),
            "owner_rate": owner_rate,
            "mass_accounting_residual": mass_residual,
            "feasibility_frontier_projection_mass": float(
                evaluation.feasibility_projection_mass
            ),
            "min_distribution_mass": period_min,
            "nonfinite_distribution_count": nonfinite,
        }
        path_rows.append(row)
        period_ages, period_children = independent_child_distribution_rows(
            label, period, evaluation.g_post_fertility, P
        )
        age_rows.extend(period_ages)
        child_rows.extend(period_children)
        calendar.write_json_atomic(
            outdir / "latest_completed_period.json",
            {
                "status": "running",
                "scenario": label,
                "completed_period": period,
                "periods_requested": periods,
                "model_solve_count": counter.total,
                "latest": row,
            },
        )
        print(
            f"{label} t={period:02d} p={price_guess:.6f} "
            f"pop={adult_population:.6f} births={evaluation.births:.6f} "
            f"queue_B={scheduled_B:.6f} market={evaluation.relative_market_residual:.2e}",
            flush=True,
        )
        g_pre = next_pre

    scenario_summary = {
        "scenario": label,
        "closure": "open_birth_vintage",
        "housing_supply_mode": supply_rule.mode,
        "outside_entry_flow_M": outside_flow,
        "renewal_retention": retention,
        "entrant_conversion_factor": conversion,
        "maturation_survival_yield": maturation_survival_yield,
        "birth_to_entry_delay_periods": delay_periods,
        "birth_to_entry_delay_years": delay_periods * float(P.period_years) + float(P.period_years),
        "periods": periods,
        "model_solve_count": counter.total - start_solves,
        "elapsed_seconds": time.perf_counter() - started,
        "max_market_residual": max_market_residual,
        "market_tolerance": market_tol,
        "smooth_date_target_tolerance": min(float(market_tol), 5e-5),
        "grid_resolution_fallback_periods": grid_resolution_fallback_periods,
        "market_resolution_note": (
            "Deterministic tenure choices on the finite wealth grid make the "
            "off-stationary demand schedule locally discontinuous; the reported "
            "price is the first bracketed grid equilibrium within this tolerance."
        ),
        "max_abs_mass_accounting_residual": max_mass_residual,
        "max_feasibility_frontier_projection_mass": max_projection,
        "min_distribution_mass": min_mass if math.isfinite(min_mass) else None,
        "max_nonfinite_distribution_count": max_nonfinite,
    }
    if max_market_residual > market_tol:
        raise RuntimeError(
            f"{label}: market gate failed, max residual={max_market_residual:.3e}"
        )
    if max_mass_residual > 2e-8:
        raise RuntimeError(
            f"{label}: mass gate failed, max residual={max_mass_residual:.3e}"
        )
    if max_projection > 1e-6:
        raise RuntimeError(
            f"{label}: feasibility projection gate failed, mass={max_projection:.3e}"
        )
    return path_rows, age_rows, child_rows, scenario_summary


def solve_stationary_open_endpoint(
    *,
    chain: Any,
    base_overrides: dict[str, Any],
    old_price: float,
    new_psi_child: float,
    outside_flow: float,
    retention: float,
    conversion: float,
    maturation_survival_yield: float,
    fixed_housing_stock: float,
    outdir: Path,
) -> dict[str, Any]:
    """Solve the stationary endpoint consistent with the birth-vintage closure."""
    cache: dict[float, tuple[float, dict[str, Any]]] = {}

    def evaluate(asset_price: float) -> tuple[float, dict[str, Any]]:
        key = round(float(asset_price), 12)
        if key not in cache:
            solution, parameters, price, elapsed = closure_audit.solve_pe(
                chain,
                base_overrides,
                asset_price=float(asset_price),
                psi_child=float(new_psi_child),
            )
            readout = closure_audit.readout(
                chain,
                solution,
                parameters,
                price,
                label="postshock_open_endpoint_search",
                price_ratio=float(asset_price / old_price),
                psi_child=float(new_psi_child),
                elapsed=elapsed,
            )
            queue_B = (
                conversion
                * maturation_survival_yield
                * float(solution.total_births_kfe)
            )
            denominator = float(solution.entry_rate) - retention * queue_B
            if denominator <= 0.0:
                population_scale = math.inf
                market_gap = math.inf
            else:
                population_scale = outside_flow / denominator
                market_gap = (
                    population_scale * float(readout["housing_demand_per_adult"])
                    - fixed_housing_stock
                )
            row = dict(readout)
            row.update(
                queue_mature_entrant_flow_B=queue_B,
                renewal_denominator=denominator,
                stationary_population_scale=population_scale,
                fixed_stock_market_gap=market_gap,
                fixed_stock_relative_market_gap=(
                    market_gap / fixed_housing_stock
                    if math.isfinite(market_gap)
                    else math.inf
                ),
            )
            cache[key] = float(market_gap), row
            calendar.write_csv(outdir / "stationary_endpoint_search.csv", [item[1] for item in cache.values()])
        return cache[key]

    lower = max(1e-4, 0.05 * float(old_price))
    upper = 1.25 * float(old_price)
    f_lower, _ = evaluate(lower)
    f_upper, _ = evaluate(upper)
    if not math.isfinite(f_lower) or not math.isfinite(f_upper) or f_lower * f_upper > 0.0:
        raise RuntimeError(
            "Could not bracket the postshock open stationary endpoint: "
            f"F({lower:.6g})={f_lower:.6g}, F({upper:.6g})={f_upper:.6g}"
        )
    midpoint = 0.5 * (lower + upper)
    for _ in range(28):
        midpoint = 0.5 * (lower + upper)
        f_midpoint, row = evaluate(midpoint)
        if abs(float(row["fixed_stock_relative_market_gap"])) <= 2.5e-5:
            break
        if f_lower * f_midpoint <= 0.0:
            upper = midpoint
            f_upper = f_midpoint
        else:
            lower = midpoint
            f_lower = f_midpoint
    _, root = evaluate(midpoint)
    root["status"] = "complete"
    calendar.write_csv(outdir / "stationary_open_endpoint.csv", [root])
    return root


def make_paper_transition_plot(
    paths: list[dict[str, Any]],
    *,
    old_births: float,
    old_B: float,
    stationary_endpoint: dict[str, Any] | None,
    delay_years: float,
    outdir: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    years = np.asarray([row["years_from_start"] for row in paths], dtype=float)
    population = np.asarray([row["population_index"] for row in paths], dtype=float)
    prices = np.asarray([row["asset_price_index"] for row in paths], dtype=float)
    births = np.asarray([row["birth_children"] for row in paths], dtype=float) / old_births
    mature = np.asarray(
        [row["effective_mature_entrant_flow_B"] for row in paths], dtype=float
    ) / old_B
    period_years = float(years[1] - years[0]) if len(years) > 1 else 0.0

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), constrained_layout=True)
    axes[0].plot(years, population, lw=2.2, color="#21476b", label="Adult households")
    axes[0].plot(years, prices, lw=2.2, color="#b14f32", label="House price")
    if stationary_endpoint is not None:
        axes[0].axhline(
            float(stationary_endpoint["stationary_population_scale"]),
            color="#21476b",
            lw=1.1,
            ls=":",
            label="New population steady state",
        )
        axes[0].axhline(
            float(stationary_endpoint["price_ratio"]),
            color="#b14f32",
            lw=1.1,
            ls=":",
            label="New price steady state",
        )
    axes[0].axhline(1.0, color="black", lw=0.8, alpha=0.45)
    axes[0].set(
        title="Movement toward the new steady state",
        xlabel="Years after the preference change",
        ylabel="Index (old steady state = 1)",
    )
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2)

    axes[1].plot(years, births, lw=2.2, color="#4f7f3f", label="Births")
    axes[1].plot(
        years + period_years,
        mature,
        lw=2.2,
        color="#6a4c93",
        label="Locally born entrants",
    )
    axes[1].axvline(delay_years, color="black", lw=1.0, ls="--", alpha=0.65)
    axes[1].axhline(1.0, color="black", lw=0.8, alpha=0.45)
    axes[1].set(
        title="Births reach entry with a delay",
        xlabel="Years after the preference change",
        ylabel="Flow index (old steady state = 1)",
    )
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.2)
    fig.savefig(outdir / "steady_state_transition.png", dpi=220)
    fig.savefig(outdir / "steady_state_transition.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.periods < 1:
        raise ValueError("--periods must be positive")
    if not 0.0 < args.outside_origin_entry_share < 1.0:
        raise ValueError("--outside-origin-entry-share must lie strictly between zero and one")
    if args.smoke:
        args.periods = min(int(args.periods), 4)

    started = time.perf_counter()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    source = args.source.resolve()
    chain, model = configure_sequential_model()
    winner, metadata = closure_audit.load_winner(source)
    theta = closure_audit.theta_from_winner(winner)
    base = closure_audit.make_overrides(
        chain, theta, nb=120, profile="e5f-floor"
    )
    old_overrides = dict(base)
    old_overrides["psi_child"] = float(args.old_psi_child)

    print("OLD_STEADY_STATE_START", flush=True)
    old_solution, old_parameters, old_price, old_seconds = closure_audit.solve_ge(
        chain, old_overrides
    )
    if int(old_parameters.J) != 17 or int(old_parameters.Nb) != 120:
        raise RuntimeError(
            f"Production dimension gate failed: J={old_parameters.J}, Nb={old_parameters.Nb}"
        )
    if not model.independent_child_maturation_active(old_parameters):
        raise RuntimeError("Source routing did not activate independent child maturation")
    if not bool(getattr(old_parameters, "child_room_floor", False)):
        raise RuntimeError("Source routing did not activate the child-room floor")

    b_grid = np.asarray(old_solution.b_grid, dtype=float)
    old_shared = model.precompute_shared(old_parameters, b_grid)
    old_parameters._fert2_probs = np.asarray(old_solution.fert2_probs, dtype=float).copy()
    old_policy = calendar.policy_from_solution(
        old_solution, old_price, old_parameters, b_grid, old_shared
    )

    # Install the sequential operators before calling any generic calendar helper.
    calendar.apply_fertility = apply_sequential_fertility
    calendar.advance_calendar_distribution = advance_sequential_calendar_distribution
    calendar.distribution_rows = independent_child_distribution_rows

    initial_g_pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(
        old_solution,
        old_policy,
        old_parameters,
        b_grid,
        old_shared,
    )
    gates = operator_gates(
        old_solution,
        old_policy,
        initial_g_pre,
        old_parameters,
        b_grid,
        old_shared,
    )
    gates.update(reconstruction)
    gate_tolerance = 5e-9
    gates["passed"] = bool(
        gates["stationary_post_fertility_nesting_l1"] <= gate_tolerance
        and gates["one_step_constant_path_nesting_l1"] <= gate_tolerance
        and gates["mature_flow_abs_error"] <= gate_tolerance
        and gates["birth_flow_abs_error"] <= gate_tolerance
    )
    if not gates["passed"]:
        write_json(
            outdir / "summary.json",
            {
                "status": "failed_operator_gates",
                "source": source,
                "old_psi_child": args.old_psi_child,
                "gates": gates,
            },
        )
        raise RuntimeError(f"Sequential calendar operator failed exact nesting: {gates}")
    print(
        "OPERATOR_GATES_PASSED "
        f"one_step_l1={gates['one_step_constant_path_nesting_l1']:.3e} "
        f"mature_error={gates['mature_flow_abs_error']:.3e}",
        flush=True,
    )

    supply_rule, supply_normalization = calendar.normalize_date0_housing_supply(
        initial_g_pre,
        old_policy,
        old_parameters,
        b_grid,
        old_shared,
        "fixed-stock",
    )

    E_old = float(old_solution.entry_rate)
    B_old = float(old_solution.entrants_mature_total)
    outside_share = float(args.outside_origin_entry_share)
    retention = (1.0 - outside_share) * E_old / B_old
    if not 0.0 <= retention <= 1.0:
        raise RuntimeError(
            "Initial outside-origin share implies an invalid retention rate: "
            f"rho={retention:.8g}"
        )
    conversion = float(old_parameters.entrant_conversion_factor)
    B_old_raw = B_old / conversion

    postshock_parameters = copy.deepcopy(old_parameters)
    postshock_parameters.psi_child = float(args.new_psi_child)
    postshock_parameters._fert2_probs = np.asarray(old_solution.fert2_probs, dtype=float).copy()
    counter = calendar.SolveCounter()
    if args.renewal_clock == "birth-vintage":
        scenario = "preference_decline_open_birth_vintage_fixed_stock"
        paths, ages, children, scenario_summary = run_birth_vintage_scenario(
            label=scenario,
            initial_g_pre=initial_g_pre,
            baseline_price=float(np.asarray(old_price).reshape(-1)[0]),
            baseline_B=B_old,
            baseline_births=float(old_solution.total_births_kfe),
            parameters=postshock_parameters,
            b_grid=b_grid,
            periods=int(args.periods),
            outside_flow=outside_share * E_old,
            retention=retention,
            conversion=conversion,
            delay_periods=int(args.birth_to_entry_delay_periods),
            supply_rule=supply_rule,
            market_tol=float(args.market_tol),
            market_max_iter=int(args.market_max_iter),
            outdir=outdir,
            counter=counter,
        )
    else:
        scenario = "preference_decline_open_household_state_fixed_stock"
        paths, ages, children, scenario_summary = calendar.run_scenario(
            scenario,
            "nesting",
            initial_g_pre,
            None,
            float(np.asarray(old_price).reshape(-1)[0]),
            E_old,
            B_old_raw,
            postshock_parameters,
            b_grid,
            int(args.periods),
            retention,
            conversion,
            supply_rule,
            float(args.market_tol),
            int(args.market_max_iter),
            outdir,
            counter,
        )
    calendar.write_csv(outdir / "transition_path.csv", paths)
    calendar.write_csv(outdir / "adult_age_distribution.csv", ages)
    calendar.write_csv(outdir / "child_pipeline.csv", children)
    calendar.make_plots(paths, ages, children, outdir)

    stationary_endpoint = None
    if args.renewal_clock == "birth-vintage" and not args.smoke:
        print("STATIONARY_ENDPOINT_START", flush=True)
        stationary_endpoint = solve_stationary_open_endpoint(
            chain=chain,
            base_overrides=base,
            old_price=float(np.asarray(old_price).reshape(-1)[0]),
            new_psi_child=float(args.new_psi_child),
            outside_flow=outside_share * E_old,
            retention=retention,
            conversion=conversion,
            maturation_survival_yield=float(
                scenario_summary["maturation_survival_yield"]
            ),
            fixed_housing_stock=float(supply_rule.initial_stock),
            outdir=outdir,
        )
        print(
            "STATIONARY_ENDPOINT_DONE "
            f"p={stationary_endpoint['asset_price']:.6f} "
            f"S={stationary_endpoint['stationary_population_scale']:.6f}",
            flush=True,
        )
    make_paper_transition_plot(
        paths,
        old_births=float(old_solution.total_births_kfe),
        old_B=B_old,
        stationary_endpoint=stationary_endpoint,
        delay_years=(
            float(scenario_summary.get("birth_to_entry_delay_years", 0.0))
            if args.renewal_clock == "birth-vintage"
            else 0.0
        ),
        outdir=outdir,
    )

    old_row = closure_audit.readout(
        chain,
        old_solution,
        old_parameters,
        old_price,
        label="old_stationary_equilibrium",
        price_ratio=1.0,
        psi_child=float(args.old_psi_child),
        elapsed=old_seconds,
    )
    final_row = paths[-1]
    peak_population = max(paths, key=lambda row: float(row["adult_population"]))
    peak_price = max(paths, key=lambda row: float(row["asset_price"]))
    outside_flow = outside_share * E_old
    summary = {
        "status": "complete",
        "method": "exact cohort accounting with temporary-equilibrium household policies",
        "not_perfect_foresight": True,
        "source": source,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_arm": metadata.get("arm", metadata.get("winner_metadata", {}).get("arm")),
        "theta": theta,
        "old_psi_child": float(args.old_psi_child),
        "new_psi_child": float(args.new_psi_child),
        "renewal_clock": args.renewal_clock,
        "birth_to_entry_delay_periods": int(args.birth_to_entry_delay_periods),
        "period_years": float(old_parameters.period_years),
        "periods": int(args.periods),
        "old_steady_state": old_row,
        "renewal_closure": {
            "equation": "E_(t+1) = M + rho * B_t",
            "outside_origin_entry_share_at_old_steady_state": outside_share,
            "outside_entry_flow_M": outside_flow,
            "retention_rho": retention,
            "entrant_conversion_factor": conversion,
            "old_entry_flow_E": E_old,
            "old_effective_mature_flow_B": B_old,
            "identity_residual": E_old - outside_flow - retention * B_old,
        },
        "operator_gates": gates,
        "supply_normalization": supply_normalization,
        "scenario_summary": scenario_summary,
        "stationary_open_endpoint": stationary_endpoint,
        "peak_population": peak_population,
        "peak_asset_price": peak_price,
        "last_simulated_period": final_row,
        "model_solve_count": counter.total + 1,
        "elapsed_seconds": time.perf_counter() - started,
        "interpretation": {
            "household_expectations": "current price and post-shock primitives are treated as permanent at each date",
            "housing_stock": "fixed at the old steady-state occupied stock",
            "population_normalization": "none after the old steady state",
            "purpose": "test whether inherited cohorts can generate short-run population and price momentum after a fertility-preference decline",
        },
    }
    write_json(outdir / "summary.json", summary)

    command = " ".join(shlex.quote(item) for item in [sys.executable, *sys.argv])
    readme = f"""# Sequential-model open population transition

This packet starts from the paper's sequential child-room-floor model at
`psi_child={args.old_psi_child:g}` and permanently moves to
`psi_child={args.new_psi_child:g}`. The outside entrant flow and retention rate
are fixed at the old steady state. The full household distribution is then
carried forward without population renormalization.

The renewal clock is `{args.renewal_clock}`. The paper-facing default uses a
separate birth-vintage queue: household children-at-home states continue to
govern housing choices, but locally born entrant households arrive only after
the declared child-to-entry delay.

The calculation is a temporary-equilibrium transition diagnostic: households
treat the current price as permanent at each date. It is not yet the paper's
perfect-foresight asset-price transition and should not be used for welfare.

Exact command:

```bash
{command}
```

Key acceptance checks are in `summary.json`; the exact path is in
`transition_path.csv` and the standard visual packet is regenerated by the
same command.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")
    print("TRANSITION_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
