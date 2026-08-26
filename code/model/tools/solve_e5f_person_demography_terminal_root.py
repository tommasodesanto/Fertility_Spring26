#!/usr/bin/env python3
"""Jointly solve a corrected terminal housing and rebate equilibrium.

The fertility-preference intercept is a declared input. At every candidate
asset price and equal transfer, the household Bellman problem is solved, the
coherent demographic-household stationary mapping is closed, and the signed
housing and government-budget residuals are returned to a two-dimensional
root finder. All evaluations are checkpointed.
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
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

import build_e5f_person_demography_terminal_point as point  # noqa: E402
import run_dynamic_population_transition as calendar  # noqa: E402
import run_e5f_open_population_transition as transition  # noqa: E402
import run_e5f_perfect_foresight_person_demography as person_pf  # noqa: E402
import run_e5f_perfect_foresight_rebated_property_tax as funded  # noqa: E402
import run_e5f_perfect_foresight_transition as pf  # noqa: E402
import run_e5f_post2023_coven_property_tax_smoke as impact  # noqa: E402
import run_e5f_post2023_no_policy_continuations as continuation  # noqa: E402


DEFAULT_OUTPUT = ROOT / "output/model/e5f_person_demography_terminal_root_20260826a"
MARKET_TOLERANCE = 2e-4
FISCAL_ABSOLUTE_TOLERANCE = 2.5e-5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=sorted(impact.CASES), required=True)
    parser.add_argument("--psi-child", type=float, required=True)
    parser.add_argument("--initial-asset-price", type=float, required=True)
    parser.add_argument("--initial-equal-transfer", type=float, required=True)
    parser.add_argument("--selected-report", type=Path, default=pf.DEFAULT_REPORT)
    parser.add_argument(
        "--selected-transition", type=Path, default=pf.DEFAULT_SELECTED_TRANSITION
    )
    parser.add_argument("--source", type=Path, default=pf.DEFAULT_SOURCE)
    parser.add_argument(
        "--source-dir", type=Path, default=person_pf.DEFAULT_SOURCE_DIR
    )
    parser.add_argument(
        "--headship-dir", type=Path, default=person_pf.DEFAULT_HEADSHIP_DIR
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--maximum-root-evaluations", type=int, default=24)
    parser.add_argument("--maximum-root-iterations", type=int, default=16)
    parser.add_argument("--log-price-difference-step", type=float, default=0.01)
    parser.add_argument("--transfer-difference-step", type=float, default=0.01)
    parser.add_argument("--maximum-inner-iterations", type=int, default=250)
    parser.add_argument("--inner-damping", type=float, default=0.50)
    return parser.parse_args()


def atomic_json(path: Path, payload: Any) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(pf.jsonable(payload), indent=2, sort_keys=True, default=str)
        + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def atomic_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def main() -> None:
    args = parse_args()
    if not math.isfinite(float(args.psi_child)):
        raise ValueError("psi_child must be finite")
    if (
        not math.isfinite(float(args.initial_asset_price))
        or float(args.initial_asset_price) <= 0.0
    ):
        raise ValueError("Initial asset price must be finite and positive")
    if (
        not math.isfinite(float(args.initial_equal_transfer))
        or float(args.initial_equal_transfer) <= 0.0
    ):
        raise ValueError("Initial equal transfer must be finite and positive")
    if int(args.maximum_root_evaluations) < 4:
        raise ValueError("At least four root evaluations are required")
    if int(args.maximum_root_iterations) < 1:
        raise ValueError("At least one root iteration is required")
    if (
        float(args.log_price_difference_step) <= 0.0
        or float(args.transfer_difference_step) <= 0.0
    ):
        raise ValueError("Finite-difference steps must be positive")
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output_dir}")
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
    initial_state, _, history_gate, impact_price = pf.reconstruct_2023_state(
        prepared,
        contracts,
        args.selected_transition.resolve(),
        market_tolerance=MARKET_TOLERANCE,
    )
    primitives = person_pf.build_annual_demographic_primitives(
        initial_state.g_pre,
        prepared.parameters,
        source_dir=args.source_dir,
        headship_dir=args.headship_dir,
    )
    anchor_state = continuation.DynamicState(
        g_pre=initial_state.g_pre.copy(),
        scheduled_entries=list(initial_state.scheduled_entries),
        scheduled_raw_entries=list(initial_state.scheduled_raw_entries),
        price_guess=float(impact_price),
        initial_policy=None,
    )
    anchor_parameters = copy.deepcopy(prepared.parameters)
    anchor_parameters.psi_child = float(contracts["report_best"]["new_psi_child"])
    anchor = impact.solve_fixed_price_equal_rebate(
        state=anchor_state,
        P=anchor_parameters,
        b_grid=prepared.b_grid,
        price=float(impact_price),
        case=impact.CASES["rebated-tax1-baseline"],
        counter=calendar.SolveCounter(),
    )
    supply_rule, supply_contract = impact.reanchor_supply_rule(
        inherited_price=float(impact_price),
        baseline_demand=float(anchor.evaluation.demand_by_loc[0]),
        elasticity=1.75,
    )

    case = impact.CASES[args.case]
    evaluations: list[dict[str, Any]] = []
    cache: dict[tuple[float, float], SimpleNamespace] = {}

    class EvaluationBudgetExhausted(RuntimeError):
        pass

    def evaluate(point_value: np.ndarray) -> SimpleNamespace:
        log_point = np.asarray(point_value, dtype=float)
        price = float(math.exp(float(log_point[0])))
        transfer = float(log_point[1])
        if not math.isfinite(transfer) or transfer < 0.0:
            raise ValueError(f"Transfer must be finite and nonnegative, got {transfer}")
        key = (round(price, 12), round(transfer, 12))
        if key not in cache:
            if len(cache) >= int(args.maximum_root_evaluations):
                raise EvaluationBudgetExhausted(
                    "Terminal-root model-evaluation budget exhausted: "
                    f"{len(cache)}/{args.maximum_root_evaluations}"
                )
            evaluation_started = time.perf_counter()
            stationary = funded._stationary_evaluation(
                prepared=prepared,
                case=case,
                asset_price=price,
                transfer=transfer,
                psi_child=float(args.psi_child),
            )
            inner = person_pf.solve_terminal_household_person_fixed_point(
                policy=stationary.policy,
                parameters=stationary.parameters,
                b_grid=stationary.b_grid,
                initial_g_pre=stationary.unit_g_pre,
                demographic_primitives=primitives,
                supply_rule=supply_rule,
                maximum_iterations=int(args.maximum_inner_iterations),
                damping=float(args.inner_damping),
            )
            if not inner.converged:
                raise RuntimeError(
                    "Inner terminal mapping failed at "
                    f"price={price}, transfer={transfer}: {inner}"
                )
            head_mass = float(np.sum(inner.g_pre))
            transfer_outlays = transfer * head_mass
            property_tax_revenue = (
                float(inner.government_budget_residual) + transfer_outlays
            )
            fiscal_scale = max(
                abs(property_tax_revenue), abs(transfer_outlays), 0.1
            )
            residual = np.asarray(
                [
                    inner.housing_excess_demand_relative,
                    inner.government_budget_residual / fiscal_scale,
                ],
                dtype=float,
            )
            cache[key] = SimpleNamespace(
                point=log_point.copy(),
                price=price,
                transfer=transfer,
                stationary=stationary,
                inner=inner,
                residual=residual,
            )
            row = {
                "evaluation": len(evaluations) + 1,
                "asset_price": price,
                "equal_transfer_period_units": transfer,
                "housing_excess_demand_relative": (
                    inner.housing_excess_demand_relative
                ),
                "housing_market_relative_residual": (
                    inner.housing_market_relative_residual
                ),
                "equal_transfer_gap": inner.equal_transfer_gap,
                "government_budget_residual": inner.government_budget_residual,
                "scaled_government_budget_residual": (
                    residual[1]
                ),
                "annual_births_per_head": inner.annual_births_per_head,
                "renewal_ratio": inner.renewal_ratio,
                "resident_persons_model_units": float(
                    np.sum(inner.persons.persons)
                ),
                "household_heads_model_units": float(np.sum(inner.g_pre)),
                "inner_iterations": inner.iterations,
                "inner_distribution_mapping_relative_l1": (
                    inner.distribution_mapping_relative_l1
                ),
                "inner_person_one_step_relative_l1": (
                    inner.person_one_step_relative_l1
                ),
                "evaluation_seconds": time.perf_counter() - evaluation_started,
                "elapsed_seconds": time.perf_counter() - started,
            }
            evaluations.append(row)
            atomic_csv(output_dir / "root_evaluations.csv", evaluations)
            atomic_json(
                output_dir / "latest_evaluation.json",
                {
                    "status": "running",
                    "case": case.name,
                    "psi_child": float(args.psi_child),
                    "latest": row,
                    "evaluations": len(evaluations),
                },
            )
        return cache[key]

    point_value = np.asarray(
        [math.log(float(args.initial_asset_price)), float(args.initial_equal_transfer)],
        dtype=float,
    )
    result = evaluate(point_value)
    difference_steps = np.asarray(
        [
            float(args.log_price_difference_step),
            float(args.transfer_difference_step),
        ],
        dtype=float,
    )

    def finite_difference_jacobian(
        at: np.ndarray, base_result: SimpleNamespace
    ) -> np.ndarray:
        columns = []
        for dimension in range(2):
            perturbation = np.zeros(2)
            perturbation[dimension] = difference_steps[dimension]
            trial = evaluate(at + perturbation)
            columns.append(
                (trial.residual - base_result.residual)
                / difference_steps[dimension]
            )
        return np.column_stack(columns)

    root_status = "maximum_iterations_reached"
    root_iterations = 0
    try:
        jacobian = finite_difference_jacobian(point_value, result)
        for iteration in range(int(args.maximum_root_iterations)):
            inner = result.inner
            current_gates = (
                abs(float(result.residual[0])) <= MARKET_TOLERANCE
                and abs(float(inner.government_budget_residual))
                <= FISCAL_ABSOLUTE_TOLERANCE
                and abs(float(inner.equal_transfer_gap))
                <= FISCAL_ABSOLUTE_TOLERANCE
            )
            if current_gates:
                root_status = "converged_declared_market_and_fiscal_gates"
                break
            try:
                direction = np.linalg.solve(jacobian, -result.residual)
            except np.linalg.LinAlgError:
                direction = np.linalg.lstsq(
                    jacobian, -result.residual, rcond=None
                )[0]
            direction = np.clip(direction, [-0.25, -0.15], [0.25, 0.15])
            current_norm = float(np.linalg.norm(result.residual))
            accepted: SimpleNamespace | None = None
            accepted_point: np.ndarray | None = None
            for line_search_damping in (1.0, 0.5, 0.25, 0.125):
                trial_point = point_value + line_search_damping * direction
                trial_point[0] = float(
                    np.clip(
                        trial_point[0],
                        math.log(
                            max(
                                float(getattr(prepared.parameters, "p_min", 1e-4)),
                                1e-8,
                            )
                        ),
                        math.log(float(getattr(prepared.parameters, "p_max", 100.0))),
                    )
                )
                trial_point[1] = float(np.clip(trial_point[1], 0.0, 8.0))
                trial = evaluate(trial_point)
                if float(np.linalg.norm(trial.residual)) < current_norm:
                    accepted = trial
                    accepted_point = trial_point
                    break
            if accepted is None or accepted_point is None:
                difference_steps *= 0.5
                jacobian = finite_difference_jacobian(point_value, result)
                continue
            displacement = accepted_point - point_value
            residual_change = accepted.residual - result.residual
            denominator = float(displacement @ displacement)
            if denominator > 1e-14:
                jacobian += np.outer(
                    residual_change - jacobian @ displacement,
                    displacement,
                ) / denominator
            point_value = accepted_point
            result = accepted
            root_iterations = iteration + 1
        else:
            root_status = "maximum_iterations_reached"
    except EvaluationBudgetExhausted:
        root_status = "maximum_model_evaluations_reached"

    if (
        abs(float(result.residual[0])) <= MARKET_TOLERANCE
        and abs(float(result.inner.government_budget_residual))
        <= FISCAL_ABSOLUTE_TOLERANCE
        and abs(float(result.inner.equal_transfer_gap))
        <= FISCAL_ABSOLUTE_TOLERANCE
    ):
        root_status = "converged_declared_market_and_fiscal_gates"

    solved_price = float(result.price)
    solved_transfer = float(result.transfer)
    final_residual = np.asarray(result.residual, dtype=float)
    final = result.inner
    gates = {
        "inner_terminal_mapping": bool(final.converged),
        "finite_demographic_steady_state": final.renewal_ratio < 1.0,
        "housing_market": abs(float(final_residual[0])) <= MARKET_TOLERANCE,
        "government_budget_absolute": (
            abs(float(final.government_budget_residual))
            <= FISCAL_ABSOLUTE_TOLERANCE
        ),
        "equal_transfer_absolute_gap": (
            abs(float(final.equal_transfer_gap))
            <= FISCAL_ABSOLUTE_TOLERANCE
        ),
        "person_one_step": final.person_one_step_relative_l1 <= 1e-8,
        "head_one_step": final.head_one_step_relative_l1 <= 1e-8,
    }
    point.write_terminal_figure(
        final,
        result.stationary.parameters,
        output_dir / "terminal_diagnostics.png",
    )
    summary = {
        "status": (
            "complete_corrected_terminal_root"
            if all(gates.values())
            else "failed_corrected_terminal_root_gates"
        ),
        "promotion_status": "not_promoted",
        "case": case.name,
        "annual_property_tax_rate": float(case.annual_tax_rate),
        "psi_child": float(args.psi_child),
        "asset_price": float(solved_price),
        "renter_price": (
            float(prepared.parameters.q)
            + float(prepared.parameters.delta)
            + float(case.annual_tax_rate) * float(prepared.parameters.period_years)
        )
        * float(solved_price),
        "equal_transfer_period_units": float(solved_transfer),
        "resident_persons_model_units": float(np.sum(final.persons.persons)),
        "resident_persons_actual_scale": float(np.sum(final.persons.persons))
        / primitives.scale_model_units_per_person,
        "household_heads_model_units": float(np.sum(final.g_pre)),
        "annual_births_per_head": final.annual_births_per_head,
        "renewal_ratio": final.renewal_ratio,
        "housing_excess_demand_relative": final.housing_excess_demand_relative,
        "government_budget_residual": final.government_budget_residual,
        "scaled_government_budget_residual": (
            final.scaled_government_budget_residual
        ),
        "implied_equal_transfer": final.implied_equal_transfer,
        "equal_transfer_gap": final.equal_transfer_gap,
        "root_solver": {
            "success": bool(all(gates.values())),
            "status": root_status,
            "iterations": int(root_iterations),
            "distinct_model_evaluations": len(evaluations),
            "maximum_model_evaluations": int(args.maximum_root_evaluations),
            "maximum_iterations": int(args.maximum_root_iterations),
            "final_difference_steps": difference_steps.tolist(),
        },
        "declared_tolerances": {
            "housing_market_relative": MARKET_TOLERANCE,
            "government_budget_absolute": FISCAL_ABSOLUTE_TOLERANCE,
            "equal_transfer_absolute_gap": FISCAL_ABSOLUTE_TOLERANCE,
            "person_one_step_relative_l1": 1e-8,
            "head_one_step_relative_l1": 1e-8,
        },
        "gates": gates,
        "inner_fixed_point": {
            "iterations": final.iterations,
            "distribution_mapping_relative_l1": (
                final.distribution_mapping_relative_l1
            ),
            "annual_birth_rate_relative_gap": (
                final.annual_birth_rate_relative_gap
            ),
            "person_one_step_relative_l1": final.person_one_step_relative_l1,
            "head_one_step_relative_l1": final.head_one_step_relative_l1,
            "age_head_one_step_max_abs": final.age_head_one_step_max_abs,
            "household_person_head_gap": final.household_person_head_gap,
        },
        "supply_contract": supply_contract,
        "terminal_demographic_closure": {
            "terminal_calendar_year": int(primitives.last_empirical_year),
            "survival_after_terminal_year": "held fixed at terminal-year schedule",
            "net_migration_after_terminal_year": (
                "held fixed at terminal-year inferred net-migration flow"
            ),
            "birth_sex_shares_after_terminal_year": (
                "held fixed at terminal-year shares"
            ),
            "headship_rates": (
                "fixed 2023 ACS sex-age profile with model-age-cell alignment"
            ),
            "fertility_preference": "declared psi_child input, not estimated here",
        },
        "history_reproduction_status": history_gate.get("status"),
        "headship_alignment_factors": (
            primitives.model_age_headship_alignment_factors.tolist()
        ),
        "source_provenance": provenance,
        "source_hashes": {
            "driver": pf.file_sha256(Path(__file__).resolve()),
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
            **{
                f"input_{name}": pf.file_sha256(path)
                for name, path in primitives.source_paths.items()
            },
        },
        "elapsed_seconds": time.perf_counter() - started,
    }
    atomic_json(output_dir / "summary.json", summary)
    if not all(gates.values()):
        raise RuntimeError(f"Corrected terminal root gates failed: {gates}")


if __name__ == "__main__":
    main()
