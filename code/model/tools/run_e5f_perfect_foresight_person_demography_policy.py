#!/usr/bin/env python3
"""Solve one rebated property-tax path with coherent person demography.

The terminal price and rebate are read from a separately fingerprinted joint
terminal-root packet. The driver reconstructs that terminal household/person
fixed point under the current source bundle, then iterates the complete price
and equal-rebate paths while births feed the annual person cohort law.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

import run_dynamic_population_transition as calendar  # noqa: E402
import e5f_post2023_tenure_sensitivity as tenure_sensitivity  # noqa: E402
import run_e5f_open_population_transition as transition  # noqa: E402
import run_e5f_perfect_foresight_person_demography as person_pf  # noqa: E402
import run_e5f_perfect_foresight_rebated_property_tax as funded  # noqa: E402
import run_e5f_perfect_foresight_transition as pf  # noqa: E402
import run_e5f_post2023_coven_property_tax_smoke as impact  # noqa: E402
import run_e5f_post2023_no_policy_continuations as continuation  # noqa: E402


DEFAULT_OUTPUT = ROOT / "output/model/e5f_pf_person_demography_policy_20260826a"
MARKET_TOLERANCE = 2e-4
FISCAL_ABSOLUTE_TOLERANCE = 2.5e-5
TERMINAL_TOLERANCES = {
    "resident_persons_relative_gap": 0.01,
    "household_heads_relative_gap": 0.01,
    "normalized_household_distribution_l1": 0.02,
    "normalized_person_age_sex_l1": 0.02,
    "normalized_head_age_sex_l1": 0.02,
    "asset_price_relative_gap": 0.01,
    "renter_price_relative_gap": 0.01,
    "equal_transfer_relative_gap": 0.01,
    "psi_absolute_gap": 1e-3,
}


@dataclass
class PersonPolicyTerminalState:
    case: Any
    psi_child: float
    asset_price: float
    renter_price: float
    equal_transfer: float
    policy: calendar.PolicyBundle
    parameters: SimpleNamespace
    state: person_pf.PersonPFState
    fixed_point: person_pf.TerminalHouseholdPersonFixedPoint
    root_summary: dict[str, Any]


@dataclass
class PersonFundedPFSolution:
    case: str
    status: str
    path_converged: bool
    iterations: int
    prices: np.ndarray
    rents: np.ndarray
    transfers: np.ndarray
    psi_path: np.ndarray
    rows: list[dict[str, Any]]
    iteration_history: list[dict[str, Any]]
    maximum_market_residual: float
    maximum_fiscal_residual: float
    maximum_log_price_residual: float
    maximum_transfer_residual: float
    maximum_policy_reproduction_error: float
    maximum_person_identity_error: float
    maximum_head_identity_error: float
    maximum_household_person_head_gap: float
    maximum_age_head_gap: float
    maximum_feasibility_projection_mass: float
    terminal_convergence: dict[str, Any]
    bellman_solves: int
    elapsed_seconds: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=sorted(impact.CASES), required=True)
    parser.add_argument("--terminal-summary", type=Path, required=True)
    parser.add_argument("--horizon", type=int, required=True)
    parser.add_argument("--selected-report", type=Path, default=pf.DEFAULT_REPORT)
    parser.add_argument(
        "--selected-transition", type=Path, default=pf.DEFAULT_SELECTED_TRANSITION
    )
    parser.add_argument("--source", type=Path, default=pf.DEFAULT_SOURCE)
    parser.add_argument("--source-dir", type=Path, default=person_pf.DEFAULT_SOURCE_DIR)
    parser.add_argument(
        "--headship-dir", type=Path, default=person_pf.DEFAULT_HEADSHIP_DIR
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--psi-persistence", type=float, default=pf.DEFAULT_PSI_PERSISTENCE)
    parser.add_argument("--maximum-path-iterations", type=int, default=35)
    parser.add_argument("--price-damping", type=float, default=0.25)
    parser.add_argument("--transfer-damping", type=float, default=0.50)
    parser.add_argument("--maximum-log-price-step", type=float, default=0.12)
    parser.add_argument("--maximum-transfer-step", type=float, default=0.08)
    parser.add_argument("--initial-path", type=Path)
    parser.add_argument("--post2023-tenure-choice-kappa", type=float)
    parser.add_argument(
        "--soft-time-limit-seconds",
        type=float,
        help=(
            "Optional wall-clock budget. After the first completed iteration "
            "that reaches this budget, write the normal best-path packet and "
            "exit with a restartable soft-time status."
        ),
    )
    return parser.parse_args()


def write_json(path: Path, payload: Any) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(pf.jsonable(payload), indent=2, sort_keys=True, default=str)
        + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def normalized_l1(candidate: np.ndarray, reference: np.ndarray) -> float:
    left = np.asarray(candidate, dtype=float)
    right = np.asarray(reference, dtype=float)
    if left.shape != right.shape:
        raise ValueError(f"Shape mismatch: {left.shape} versus {right.shape}")
    left_mass = float(np.sum(left))
    right_mass = float(np.sum(right))
    if left_mass <= 0.0 or right_mass <= 0.0:
        raise ValueError("Normalized distributions must have positive mass")
    return float(np.sum(np.abs(left / left_mass - right / right_mass)))


def terminal_supply_elasticity(root_summary: dict[str, Any]) -> float:
    """Return the terminal packet's explicit, validated direct elasticity."""
    contract = root_summary.get("supply_contract")
    if not isinstance(contract, dict):
        raise RuntimeError("Terminal root lacks a supply contract")
    try:
        elasticity = float(contract["supply_elasticity"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError(
            "Terminal root lacks a finite positive supply elasticity"
        ) from exc
    if not math.isfinite(elasticity) or elasticity <= 0.0:
        raise RuntimeError(
            "Terminal root lacks a finite positive supply elasticity"
        )
    return elasticity


def terminal_convergence_diagnostics(
    evaluation: person_pf.PersonPathEvaluation,
    *,
    terminal: PersonPolicyTerminalState,
    psi_path: Sequence[float],
) -> dict[str, Any]:
    candidate = evaluation.terminal_state
    reference = terminal.state
    candidate_people = float(np.sum(candidate.persons.persons))
    reference_people = float(np.sum(reference.persons.persons))
    candidate_heads = float(np.sum(candidate.g_pre))
    reference_heads = float(np.sum(reference.g_pre))
    last_price = float(evaluation.prices[-1])
    last_rent = float(evaluation.rents[-1])
    last_transfer = float(evaluation.rows[-1]["equal_transfer_period_units"])
    metrics = {
        "resident_persons_relative_gap": abs(candidate_people - reference_people)
        / reference_people,
        "household_heads_relative_gap": abs(candidate_heads - reference_heads)
        / reference_heads,
        "normalized_household_distribution_l1": normalized_l1(
            candidate.g_pre, reference.g_pre
        ),
        "normalized_person_age_sex_l1": normalized_l1(
            candidate.persons.persons, reference.persons.persons
        ),
        "normalized_head_age_sex_l1": normalized_l1(
            candidate.persons.heads, reference.persons.heads
        ),
        "asset_price_relative_gap": abs(last_price - terminal.asset_price)
        / terminal.asset_price,
        "renter_price_relative_gap": abs(last_rent - terminal.renter_price)
        / terminal.renter_price,
        "equal_transfer_relative_gap": abs(last_transfer - terminal.equal_transfer)
        / max(abs(terminal.equal_transfer), 1e-15),
        "psi_absolute_gap": abs(float(psi_path[-1]) - terminal.psi_child),
    }
    checks = {
        name: metrics[name] <= tolerance
        for name, tolerance in TERMINAL_TOLERANCES.items()
    }
    return {
        "status": "passed" if all(checks.values()) else "not_converged",
        "all_checks_pass": all(checks.values()),
        "metrics": metrics,
        "checks": checks,
        "tolerances": TERMINAL_TOLERANCES,
        "candidate_terminal_year": int(candidate.persons.year),
        "reference_terminal_year_label": int(reference.persons.year),
        "candidate_resident_persons_model_units": candidate_people,
        "reference_resident_persons_model_units": reference_people,
        "candidate_household_heads_model_units": candidate_heads,
        "reference_household_heads_model_units": reference_heads,
    }


def solve_person_funded_path(
    *,
    initial_prices: Sequence[float],
    initial_transfers: Sequence[float],
    psi_path: Sequence[float],
    terminal: PersonPolicyTerminalState,
    b_grid: np.ndarray,
    initial_state: person_pf.PersonPFState,
    primitives: person_pf.AnnualDemographicPrimitives,
    supply_rule: calendar.HousingSupplyRule,
    maximum_iterations: int,
    price_damping: float,
    transfer_damping: float,
    maximum_log_price_step: float,
    maximum_transfer_step: float,
    progress_callback: Callable[
        [dict[str, Any], person_pf.PersonPathEvaluation], None
    ]
    | None = None,
    stop_after_iteration: Callable[[dict[str, Any]], bool] | None = None,
) -> PersonFundedPFSolution:
    started = time.perf_counter()
    prices = np.asarray(initial_prices, dtype=float).reshape(-1).copy()
    transfers = np.asarray(initial_transfers, dtype=float).reshape(-1).copy()
    psi_values = np.asarray(psi_path, dtype=float).reshape(-1)
    if not (prices.shape == transfers.shape == psi_values.shape):
        raise ValueError("Price, transfer, and preference paths must have equal length")
    current_price_damping = float(price_damping)
    current_transfer_damping = float(transfer_damping)
    previous_score = math.inf
    history: list[dict[str, Any]] = []
    best: tuple[float, float, float, np.ndarray, np.ndarray, int] | None = None
    total_bellman = 0
    path_converged = False
    soft_time_budget_reached = False

    for iteration in range(1, int(maximum_iterations) + 1):
        evaluation = person_pf.evaluate_path_at_prices_person_demography(
            prices=prices,
            psi_path=psi_values,
            terminal_price=terminal.asset_price,
            terminal_V=terminal.policy.V,
            base_parameters=terminal.parameters,
            b_grid=b_grid,
            initial_state=initial_state,
            demographic_primitives=primitives,
            supply_rule=supply_rule,
            transfer_path=transfers,
        )
        total_bellman += evaluation.bellman_solves
        demand = np.asarray(
            [float(row["housing_demand"]) for row in evaluation.rows]
        )
        target_prices = pf.inverse_supply_price(demand, supply_rule)
        target_transfers = np.asarray(
            [
                float(row["property_tax_revenue"])
                / max(float(row["household_heads"]), 1e-15)
                for row in evaluation.rows
            ]
        )
        log_residual = np.log(target_prices) - np.log(prices)
        transfer_residual = target_transfers - transfers
        market_error = float(evaluation.maximum_market_residual)
        fiscal_error = max(
            abs(float(row["government_budget_residual"]))
            for row in evaluation.rows
        )
        max_log = float(np.max(np.abs(log_residual)))
        max_transfer = float(np.max(np.abs(transfer_residual)))
        score = max(
            market_error / MARKET_TOLERANCE,
            fiscal_error / FISCAL_ABSOLUTE_TOLERANCE,
        )
        if best is None or (score, market_error, fiscal_error) < best[:3]:
            best = (
                score,
                market_error,
                fiscal_error,
                evaluation.prices.copy(),
                transfers.copy(),
                iteration,
            )
        record = {
            "iteration": iteration,
            "maximum_market_residual": market_error,
            "maximum_fiscal_residual": fiscal_error,
            "maximum_log_price_residual": max_log,
            "maximum_transfer_residual": max_transfer,
            "price_damping": current_price_damping,
            "transfer_damping": current_transfer_damping,
            "minimum_asset_price": float(np.min(prices)),
            "maximum_asset_price": float(np.max(prices)),
            "minimum_equal_transfer": float(np.min(transfers)),
            "maximum_equal_transfer": float(np.max(transfers)),
            "bellman_solves_this_iteration": evaluation.bellman_solves,
            "bellman_solves_cumulative": total_bellman,
            "elapsed_seconds": time.perf_counter() - started,
        }
        history.append(record)
        if progress_callback is not None:
            progress_callback(record, evaluation)
        if market_error <= MARKET_TOLERANCE and fiscal_error <= FISCAL_ABSOLUTE_TOLERANCE:
            path_converged = True
            best = (
                score,
                market_error,
                fiscal_error,
                evaluation.prices.copy(),
                transfers.copy(),
                iteration,
            )
            break
        if stop_after_iteration is not None and stop_after_iteration(record):
            soft_time_budget_reached = True
            break
        current_price_damping, current_transfer_damping = adapt_path_damping(
            current_price_damping=current_price_damping,
            current_transfer_damping=current_transfer_damping,
            declared_price_damping=float(price_damping),
            declared_transfer_damping=float(transfer_damping),
            score=score,
            previous_score=previous_score,
        )
        price_step = np.clip(
            log_residual,
            -float(maximum_log_price_step),
            float(maximum_log_price_step),
        )
        transfer_step = np.clip(
            transfer_residual,
            -float(maximum_transfer_step),
            float(maximum_transfer_step),
        )
        candidate_prices = np.exp(
            np.log(prices) + current_price_damping * price_step
        )
        candidate_transfers = np.maximum(
            transfers + current_transfer_damping * transfer_step, 0.0
        )
        prices, feasibility = funded.project_price_path_to_positive_rents(
            candidate_prices, terminal=terminal, minimum_rent_share=1e-6
        )
        transfers = candidate_transfers
        record.update(
            next_price_feasibility_adjusted_period_count=feasibility[
                "adjusted_period_count"
            ],
            next_price_feasibility_maximum_relative_change=feasibility[
                "maximum_price_relative_change"
            ],
        )
        previous_score = score
        if iteration < int(maximum_iterations):
            del evaluation

    if best is None:
        raise RuntimeError("Person-demography path iteration produced no evaluation")
    _, best_market, best_fiscal, best_prices, best_transfers, best_iteration = best
    if best_iteration == len(history):
        best_evaluation = evaluation
    else:
        best_evaluation = person_pf.evaluate_path_at_prices_person_demography(
            prices=best_prices,
            psi_path=psi_values,
            terminal_price=terminal.asset_price,
            terminal_V=terminal.policy.V,
            base_parameters=terminal.parameters,
            b_grid=b_grid,
            initial_state=initial_state,
            demographic_primitives=primitives,
            supply_rule=supply_rule,
            transfer_path=best_transfers,
        )
        total_bellman += best_evaluation.bellman_solves
    terminal_diagnostics = terminal_convergence_diagnostics(
        best_evaluation, terminal=terminal, psi_path=psi_values
    )
    for row in best_evaluation.rows:
        row.update(
            policy_case=terminal.case.name,
            policy_label=terminal.case.label,
            annual_property_tax_rate=float(terminal.case.annual_tax_rate),
            equal_rebate=True,
            housing_supply_elasticity=float(supply_rule.elasticity),
            resident_persons_actual_scale=float(row["resident_persons"])
            / primitives.scale_model_units_per_person,
        )
    return PersonFundedPFSolution(
        case=terminal.case.name,
        status=(
            "converged"
            if path_converged
            else (
                "soft_time_budget_reached"
                if soft_time_budget_reached
                else "maximum_iterations_reached"
            )
        ),
        path_converged=path_converged,
        iterations=len(history),
        prices=best_evaluation.prices,
        rents=best_evaluation.rents,
        transfers=best_transfers,
        psi_path=psi_values,
        rows=best_evaluation.rows,
        iteration_history=history,
        maximum_market_residual=best_market,
        maximum_fiscal_residual=best_fiscal,
        maximum_log_price_residual=max(
            abs(
                math.log(float(row["housing_demand"]) / float(row["housing_supply"]))
            )
            / float(supply_rule.elasticity)
            for row in best_evaluation.rows
        ),
        maximum_transfer_residual=max(
            abs(
                float(row["property_tax_revenue"])
                / max(float(row["household_heads"]), 1e-15)
                - float(row["equal_transfer_period_units"])
            )
            for row in best_evaluation.rows
        ),
        maximum_policy_reproduction_error=(
            best_evaluation.maximum_policy_reproduction_error
        ),
        maximum_person_identity_error=best_evaluation.maximum_person_identity_error,
        maximum_head_identity_error=best_evaluation.maximum_head_identity_error,
        maximum_household_person_head_gap=(
            best_evaluation.maximum_household_person_head_gap
        ),
        maximum_age_head_gap=best_evaluation.maximum_age_head_gap,
        maximum_feasibility_projection_mass=(
            best_evaluation.maximum_feasibility_projection_mass
        ),
        terminal_convergence=terminal_diagnostics,
        bellman_solves=total_bellman,
        elapsed_seconds=time.perf_counter() - started,
    )


def adapt_path_damping(
    *,
    current_price_damping: float,
    current_transfer_damping: float,
    declared_price_damping: float,
    declared_transfer_damping: float,
    score: float,
    previous_score: float,
) -> tuple[float, float]:
    """Adapt damping without exceeding a deliberately lower declared value."""
    price_floor = min(0.04, float(declared_price_damping))
    transfer_floor = min(0.08, float(declared_transfer_damping))
    if score > previous_score * 1.02:
        return (
            max(price_floor, 0.5 * float(current_price_damping)),
            max(transfer_floor, 0.5 * float(current_transfer_damping)),
        )
    if score < previous_score * 0.80:
        return (
            min(
                float(declared_price_damping),
                1.08 * float(current_price_damping),
            ),
            min(
                float(declared_transfer_damping),
                1.08 * float(current_transfer_damping),
            ),
        )
    return float(current_price_damping), float(current_transfer_damping)


def initial_paths(
    *,
    horizon: int,
    impact_price: float,
    impact_transfer: float,
    terminal: PersonPolicyTerminalState,
    persistence: float,
    source: Path | None,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if source is None:
        prices = pf.initial_price_path(
            impact_price, terminal.asset_price, float(persistence), int(horizon)
        )
        weights = float(persistence) ** np.arange(int(horizon), dtype=float)
        transfers = terminal.equal_transfer + weights * (
            impact_transfer - terminal.equal_transfer
        )
        contract: dict[str, Any] = {
            "mode": "analytic_impact_to_terminal_decay",
            "persistence": float(persistence),
        }
    else:
        with source.resolve().open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        required = {"period", "asset_price", "equal_transfer_period_units"}
        if not rows or not required.issubset(rows[0]) or len(rows) > int(horizon):
            raise ValueError(f"Invalid initial path: {source}")
        if [int(row["period"]) for row in rows] != list(range(len(rows))):
            raise ValueError("Initial-path periods must be contiguous from zero")
        prices = np.full(int(horizon), terminal.asset_price)
        transfers = np.full(int(horizon), terminal.equal_transfer)
        prices[: len(rows)] = [float(row["asset_price"]) for row in rows]
        transfers[: len(rows)] = [
            float(row["equal_transfer_period_units"]) for row in rows
        ]
        contract = {
            "mode": "recorded_path_then_terminal_extension",
            "source": str(source.resolve()),
            "source_sha256": pf.file_sha256(source.resolve()),
            "source_rows": len(rows),
        }
    prices, feasibility = funded.project_price_path_to_positive_rents(
        prices, terminal=terminal, minimum_rent_share=1e-6
    )
    contract.update(
        horizon=int(horizon),
        impact_price=float(impact_price),
        impact_transfer=float(impact_transfer),
        terminal_price=terminal.asset_price,
        terminal_transfer=terminal.equal_transfer,
        feasibility=feasibility,
    )
    return prices, transfers, contract


def make_figure(
    solution: PersonFundedPFSolution,
    *,
    primitives: person_pf.AnnualDemographicPrimitives,
    terminal: PersonPolicyTerminalState,
    output_dir: Path,
) -> None:
    rows = solution.rows
    years = np.asarray([int(row["calendar_year"]) for row in rows])
    persons = np.asarray([float(row["resident_persons"]) for row in rows])
    heads = np.asarray([float(row["household_heads"]) for row in rows])
    births = np.asarray(
        [float(row["birth_children_topcode_adjusted"]) for row in rows]
    )
    market = np.asarray([abs(float(row["relative_market_residual"])) for row in rows])
    fiscal = np.asarray([abs(float(row["government_budget_residual"])) for row in rows])
    figure, axes = plt.subplots(3, 2, figsize=(11, 10), constrained_layout=True)
    axes[0, 0].plot(years, persons / primitives.scale_model_units_per_person / 1e6)
    axes[0, 0].axhline(
        float(np.sum(terminal.state.persons.persons))
        / primitives.scale_model_units_per_person
        / 1e6,
        color="0.5",
        linestyle="--",
    )
    axes[0, 0].set(title="Resident persons", ylabel="Millions")
    axes[0, 1].plot(years, heads)
    axes[0, 1].axhline(float(np.sum(terminal.state.g_pre)), color="0.5", linestyle="--")
    axes[0, 1].set(title="Household heads", ylabel="Model units")
    axes[1, 0].plot(years, births / np.maximum(heads * 4.0, 1e-15))
    axes[1, 0].axhline(
        terminal.fixed_point.annual_births_per_head, color="0.5", linestyle="--"
    )
    axes[1, 0].set(title="Annual births per household head")
    axes[1, 1].plot(years, solution.prices, label="Asset price")
    axes[1, 1].plot(years, solution.rents, label="Rent")
    axes[1, 1].axhline(terminal.asset_price, color="0.5", linestyle="--")
    axes[1, 1].legend(frameon=False)
    axes[1, 1].set(title="Housing prices")
    axes[2, 0].plot(years, solution.transfers)
    axes[2, 0].axhline(terminal.equal_transfer, color="0.5", linestyle="--")
    axes[2, 0].set(title="Equal property-tax rebate")
    axes[2, 1].semilogy(years, np.maximum(market, 1e-14), label="Housing")
    axes[2, 1].semilogy(years, np.maximum(fiscal, 1e-14), label="Fiscal")
    axes[2, 1].axhline(MARKET_TOLERANCE, color="0.55", linestyle="--")
    axes[2, 1].axhline(FISCAL_ABSOLUTE_TOLERANCE, color="0.7", linestyle=":")
    axes[2, 1].legend(frameon=False)
    axes[2, 1].set(title="Equilibrium residuals")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.set_xlabel("Year")
    figure.suptitle(terminal.case.label)
    figure.savefig(output_dir / "diagnostics.png", dpi=200)
    figure.savefig(output_dir / "diagnostics.pdf")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if int(args.horizon) < 2 or int(args.maximum_path_iterations) < 1:
        raise ValueError("Horizon must be at least two and iterations positive")
    if args.soft_time_limit_seconds is not None and (
        not math.isfinite(float(args.soft_time_limit_seconds))
        or float(args.soft_time_limit_seconds) <= 0.0
    ):
        raise ValueError("The soft time limit must be finite and positive")
    output_dir = args.output_dir.resolve()
    launcher_files = {"launch_contract.json", "driver.log", "heartbeat.json"}
    unexpected = (
        sorted(path.name for path in output_dir.iterdir() if path.name not in launcher_files)
        if output_dir.exists()
        else []
    )
    if unexpected:
        raise RuntimeError(f"Refusing to overwrite non-launcher output: {unexpected}")
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    chain, model = transition.configure_sequential_model()
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    calendar.distribution_rows = transition.independent_child_distribution_rows
    contracts, provenance = pf.load_diagnostic_contracts(
        args.selected_report.resolve(),
        args.selected_transition.resolve(),
        args.source.resolve(),
    )
    prepared = continuation.prepare_model(
        contracts, args.source.resolve(), chain, model
    )
    initial_households, history, history_gate, impact_price = pf.reconstruct_2023_state(
        prepared,
        contracts,
        args.selected_transition.resolve(),
        market_tolerance=MARKET_TOLERANCE,
    )
    tenure_sensitivity_contract = (
        tenure_sensitivity.apply_post2023_tenure_choice_kappa(
            prepared, args.post2023_tenure_choice_kappa
        )
    )
    pf.write_csv(output_dir / "reconstructed_2007_2023_history.csv", history)
    primitives = person_pf.build_annual_demographic_primitives(
        initial_households.g_pre,
        prepared.parameters,
        source_dir=args.source_dir,
        headship_dir=args.headship_dir,
    )
    root_summary_path = args.terminal_summary.resolve()
    root_summary = json.loads(root_summary_path.read_text(encoding="utf-8"))
    if root_summary.get("status") != "complete_corrected_terminal_root":
        raise RuntimeError("Terminal root is not complete")
    if root_summary.get("case") != args.case or not all(
        (root_summary.get("gates") or {}).values()
    ):
        raise RuntimeError("Terminal root case or gates do not match")
    if root_summary.get("tenure_choice_sensitivity") != tenure_sensitivity_contract:
        raise RuntimeError(
            "Terminal root and transition use different tenure-choice sensitivities"
        )
    root_provenance = dict(root_summary.get("source_provenance") or {})
    for key in ("selected_report_sha256", "selected_transition_sha256"):
        if root_provenance.get(key) != provenance.get(key):
            raise RuntimeError(
                f"Terminal root and transition use different calibration {key}."
            )
    expected_hashes = {
        "tenure_sensitivity": pf.file_sha256(
            TOOLS_ROOT / "e5f_post2023_tenure_sensitivity.py"
        ),
        "person_demography": pf.file_sha256(
            TOOLS_ROOT / "run_e5f_perfect_foresight_person_demography.py"
        ),
        "funded_policy": pf.file_sha256(
            TOOLS_ROOT / "run_e5f_perfect_foresight_rebated_property_tax.py"
        ),
        "perfect_foresight": pf.file_sha256(
            TOOLS_ROOT / "run_e5f_perfect_foresight_transition.py"
        ),
        "solver": pf.file_sha256(
            MODEL_ROOT / "intergen_eqscale_seq_optimized/solver.py"
        ),
    }
    for name, digest in expected_hashes.items():
        if (root_summary.get("source_hashes") or {}).get(name) != digest:
            raise RuntimeError(f"Terminal root used a different {name} hash")

    new_psi = float(contracts["report_best"]["new_psi_child"])
    if not math.isclose(
        float(root_summary["psi_child"]), new_psi, rel_tol=0.0, abs_tol=1e-13
    ):
        raise RuntimeError("Terminal root psi_child differs from the fitted 2023 value")
    anchor_state = continuation.DynamicState(
        g_pre=initial_households.g_pre.copy(),
        scheduled_entries=list(initial_households.scheduled_entries),
        scheduled_raw_entries=list(initial_households.scheduled_raw_entries),
        price_guess=float(impact_price),
        initial_policy=None,
    )
    anchor_parameters = copy.deepcopy(prepared.parameters)
    anchor_parameters.psi_child = new_psi
    tax1_anchor = impact.solve_fixed_price_equal_rebate(
        state=anchor_state,
        P=anchor_parameters,
        b_grid=prepared.b_grid,
        price=float(impact_price),
        case=impact.CASES["rebated-tax1-baseline"],
        counter=calendar.SolveCounter(),
    )
    supply_rule, supply_contract = impact.reanchor_supply_rule(
        inherited_price=float(impact_price),
        baseline_demand=float(tax1_anchor.evaluation.demand_by_loc[0]),
        elasticity=terminal_supply_elasticity(root_summary),
    )
    if root_summary.get("supply_contract") != supply_contract:
        raise RuntimeError("Terminal root and transition use different supply contracts")
    case = impact.CASES[args.case]
    case_impact = (
        tax1_anchor
        if args.case == "rebated-tax1-baseline"
        else impact.solve_fixed_price_equal_rebate(
            state=anchor_state,
            P=anchor_parameters,
            b_grid=prepared.b_grid,
            price=float(impact_price),
            case=case,
            counter=calendar.SolveCounter(),
        )
    )
    stationary = funded._stationary_evaluation(
        prepared=prepared,
        case=case,
        asset_price=float(root_summary["asset_price"]),
        transfer=float(root_summary["equal_transfer_period_units"]),
        psi_child=new_psi,
    )
    terminal_fixed_point = person_pf.solve_terminal_household_person_fixed_point(
        policy=stationary.policy,
        parameters=stationary.parameters,
        b_grid=stationary.b_grid,
        initial_g_pre=stationary.unit_g_pre,
        demographic_primitives=primitives,
        supply_rule=supply_rule,
    )
    reproduction = {
        "resident_persons_relative_gap": abs(
            float(np.sum(terminal_fixed_point.persons.persons))
            - float(root_summary["resident_persons_model_units"])
        )
        / float(root_summary["resident_persons_model_units"]),
        "household_heads_relative_gap": abs(
            float(np.sum(terminal_fixed_point.g_pre))
            - float(root_summary["household_heads_model_units"])
        )
        / float(root_summary["household_heads_model_units"]),
        "market_residual": abs(terminal_fixed_point.housing_excess_demand_relative),
        "fiscal_residual": abs(terminal_fixed_point.government_budget_residual),
        "person_one_step": terminal_fixed_point.person_one_step_relative_l1,
        "head_one_step": terminal_fixed_point.head_one_step_relative_l1,
    }
    reproduction_gates = {
        "resident_persons": reproduction["resident_persons_relative_gap"] <= 1e-10,
        "household_heads": reproduction["household_heads_relative_gap"] <= 1e-10,
        "housing_market": reproduction["market_residual"] <= MARKET_TOLERANCE,
        "government_budget": reproduction["fiscal_residual"]
        <= FISCAL_ABSOLUTE_TOLERANCE,
        "person_one_step": reproduction["person_one_step"] <= 1e-8,
        "head_one_step": reproduction["head_one_step"] <= 1e-8,
    }
    write_json(
        output_dir / "terminal_root_reproduction.json",
        {
            "terminal_summary": str(root_summary_path),
            "terminal_summary_sha256": pf.file_sha256(root_summary_path),
            "metrics": reproduction,
            "gates": reproduction_gates,
        },
    )
    if not all(reproduction_gates.values()):
        raise RuntimeError(f"Terminal root reproduction failed: {reproduction_gates}")
    terminal = PersonPolicyTerminalState(
        case=case,
        psi_child=new_psi,
        asset_price=float(root_summary["asset_price"]),
        renter_price=float(root_summary["renter_price"]),
        equal_transfer=float(root_summary["equal_transfer_period_units"]),
        policy=stationary.policy,
        parameters=stationary.parameters,
        state=person_pf.PersonPFState(
            g_pre=terminal_fixed_point.g_pre,
            persons=terminal_fixed_point.persons,
        ),
        fixed_point=terminal_fixed_point,
        root_summary=root_summary,
    )
    psi_values = pf.preference_path(
        new_psi, terminal.psi_child, float(args.psi_persistence), int(args.horizon)
    )
    seed_prices, seed_transfers, seed_contract = initial_paths(
        horizon=int(args.horizon),
        impact_price=float(impact_price),
        impact_transfer=float(case_impact.transfer),
        terminal=terminal,
        persistence=float(args.psi_persistence),
        source=args.initial_path,
    )
    write_json(output_dir / "initial_path_seed_contract.json", seed_contract)

    def progress(
        record: dict[str, Any], evaluation: person_pf.PersonPathEvaluation
    ) -> None:
        write_json(
            output_dir / "latest_iteration.json",
            {
                "status": "running",
                **record,
                "price_path": evaluation.prices,
                "rent_path": evaluation.rents,
                "transfer_path": [
                    float(row["equal_transfer_period_units"])
                    for row in evaluation.rows
                ],
            },
        )
        pf.write_csv(output_dir / "latest_path.csv", evaluation.rows)

    solution = solve_person_funded_path(
        initial_prices=seed_prices,
        initial_transfers=seed_transfers,
        psi_path=psi_values,
        terminal=terminal,
        b_grid=prepared.b_grid,
        initial_state=person_pf.PersonPFState(
            g_pre=initial_households.g_pre,
            persons=primitives.initial_person_state,
        ),
        primitives=primitives,
        supply_rule=supply_rule,
        maximum_iterations=int(args.maximum_path_iterations),
        price_damping=float(args.price_damping),
        transfer_damping=float(args.transfer_damping),
        maximum_log_price_step=float(args.maximum_log_price_step),
        maximum_transfer_step=float(args.maximum_transfer_step),
        progress_callback=progress,
        stop_after_iteration=(
            (
                lambda _record: time.perf_counter() - started
                >= float(args.soft_time_limit_seconds)
            )
            if args.soft_time_limit_seconds is not None
            else None
        ),
    )
    pf.write_csv(output_dir / "transition_path.csv", solution.rows)
    pf.write_csv(output_dir / "joint_iteration_history.csv", solution.iteration_history)
    make_figure(solution, primitives=primitives, terminal=terminal, output_dir=output_dir)
    accounting_gates = {
        "policy_reproduction": solution.maximum_policy_reproduction_error <= 2e-8,
        "person_identity": solution.maximum_person_identity_error <= 2e-8,
        "head_identity": solution.maximum_head_identity_error <= 2e-8,
        "household_person_head_identity": (
            solution.maximum_household_person_head_gap <= 2e-9
        ),
        "age_head_identity": solution.maximum_age_head_gap <= 2e-9,
        "distribution_feasibility": (
            solution.maximum_feasibility_projection_mass <= 2e-8
        ),
    }
    payload = {
        "status": "complete_unpromoted_person_demography_policy_path",
        "promotion_status": "not_promoted",
        "case": args.case,
        "horizon": int(args.horizon),
        "calendar_end_year": 2023 + 4 * int(args.horizon),
        "path_status": solution.status,
        "path_converged": solution.path_converged,
        "terminal_convergence": solution.terminal_convergence,
        "accounting_gates": accounting_gates,
        "maximum_market_residual": solution.maximum_market_residual,
        "maximum_fiscal_residual": solution.maximum_fiscal_residual,
        "maximum_log_price_residual": solution.maximum_log_price_residual,
        "maximum_transfer_residual": solution.maximum_transfer_residual,
        "iterations": solution.iterations,
        "soft_time_limit_seconds": args.soft_time_limit_seconds,
        "bellman_solves": solution.bellman_solves,
        "terminal_root": root_summary,
        "history_reproduction_status": history_gate.get("status"),
        "tenure_choice_sensitivity": tenure_sensitivity_contract,
        "source_provenance": provenance,
        "source_hashes": {
            "driver": pf.file_sha256(Path(__file__).resolve()),
            **expected_hashes,
            **{
                f"input_{name}": pf.file_sha256(path)
                for name, path in primitives.source_paths.items()
            },
        },
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json(output_dir / "summary.json", payload)
    if not all(accounting_gates.values()):
        raise RuntimeError(f"Accounting gates failed: {accounting_gates}")


if __name__ == "__main__":
    main()
