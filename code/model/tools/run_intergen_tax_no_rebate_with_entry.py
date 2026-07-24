#!/usr/bin/env python3
"""Compare a 1% and 2% property tax with no rebate and endogenous entry."""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from intergen_housing_fertility_optimized.calibration import extract_moments
from intergen_housing_fertility_optimized.solver import (
    precompute_shared,
    run_model_cp_dt,
    solve_markov_income_at_prices,
    upgrade_fast_markov_solution,
)
from run_intergen_funded_policy_with_entry import (
    ENTRY_TASTE_SCALE,
    TARGET_QBAR,
    calibrate_outside_objects,
    population_scale,
)
from run_intergen_funded_property_tax_test import (
    case_overrides,
    jsonable,
    load_profile,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = (
    ROOT
    / "output/model/intergen_new_moment_unrestricted_overnight_20260723_w2/report/results.json"
)
DEFAULT_OUTDIR = (
    ROOT
    / "output/model/intergen_current_m_figures_policy_20260724/tax_no_rebate_entry"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def tax_overrides(
    base_overrides: dict[str, Any],
    annual_tax: float,
    smoke: bool,
) -> dict[str, Any]:
    overrides = case_overrides(
        base_overrides,
        annual_tax,
        False,
        0.0,
        smoke,
    )
    overrides["property_tax_lump_sum_transfer"] = 0.0
    overrides["birth_entry_grant"] = False
    overrides["birth_entry_grant_amount"] = 0.0
    overrides["solve_mode"] = "ge"
    overrides.pop("p_fixed", None)
    return overrides


def solve_fixed_scale(
    base_overrides: dict[str, Any],
    annual_tax: float,
    smoke: bool,
) -> tuple[Any, Any]:
    solution, parameters, _ = run_model_cp_dt(
        tax_overrides(base_overrides, annual_tax, smoke),
        verbose=False,
    )
    return solution, parameters


def solve_scaled_price(
    fixed_solution: Any,
    parameters: Any,
    outside: dict[str, float],
    *,
    smoke: bool,
) -> SimpleNamespace:
    initial_price = float(np.asarray(fixed_solution.p_eq, dtype=float).reshape(-1)[0])
    b_grid = np.asarray(fixed_solution.b_grid, dtype=float)
    shared = precompute_shared(parameters, b_grid)
    cache: dict[float, SimpleNamespace] = {}
    initial_key = float(round(initial_price, 12))

    def evaluate(price: float) -> SimpleNamespace:
        key = float(round(price, 12))
        if key not in cache:
            if key == initial_key:
                solution = fixed_solution
                value_function = np.asarray(solution.V, dtype=float)
            else:
                solution = solve_markov_income_at_prices(
                    np.array([float(price)], dtype=float),
                    parameters,
                    b_grid,
                    verbose=False,
                    fast_stats=True,
                    SD=shared,
                    retain_payload=True,
                )
                value_function = np.asarray(solution._model_payload[0], dtype=float)
            scale, _ = population_scale(
                solution,
                parameters,
                outside,
                value_function=value_function,
            )
            if not scale.finite:
                residual = math.inf
            else:
                demand = (
                    float(scale.scale_factor)
                    * float(np.asarray(solution.housing_demand, dtype=float).reshape(-1)[0])
                )
                supply = float(parameters.H0[0]) * (
                    float(parameters.user_cost_rate)
                    * float(price)
                    / float(parameters.r_bar[0])
                ) ** float(parameters.xi_supply[0])
                residual = (demand - supply) / max(abs(supply), 1e-14)
            cache[key] = SimpleNamespace(
                price=float(price),
                solution=solution,
                scale=scale,
                residual=float(residual),
            )
            print(
                f"  p={price:.6f}, S={scale.scale_factor:.5f}, "
                f"market={residual:.2e}",
                flush=True,
            )
        return cache[key]

    p_min = float(getattr(parameters, "p_min", 0.01))
    p_max = float(getattr(parameters, "p_max", 30.0))
    lower = max(p_min, 0.75 * initial_price)
    upper = min(p_max, 1.25 * initial_price)
    lower_eval = evaluate(lower)
    upper_eval = evaluate(upper)
    expansions = 0
    while lower_eval.residual * upper_eval.residual > 0.0 and expansions < 8:
        expansions += 1
        lower = max(p_min, lower / 1.5)
        upper = min(p_max, upper * 1.5)
        lower_eval = evaluate(lower)
        upper_eval = evaluate(upper)
    if lower_eval.residual * upper_eval.residual > 0.0:
        raise RuntimeError(
            f"Could not bracket scaled price: F({lower:.4g})={lower_eval.residual:.4g}, "
            f"F({upper:.4g})={upper_eval.residual:.4g}."
        )
    tolerance = 1.0e-4 if smoke else 2.5e-5
    price = 0.5 * (lower + upper)
    for _ in range(30):
        price = 0.5 * (lower + upper)
        midpoint_eval = evaluate(price)
        if abs(midpoint_eval.residual) <= tolerance:
            break
        if lower_eval.residual * midpoint_eval.residual <= 0.0:
            upper = price
            upper_eval = midpoint_eval
        else:
            lower = price
            lower_eval = midpoint_eval
    result = evaluate(price)
    if hasattr(result.solution, "_model_payload"):
        result.solution = upgrade_fast_markov_solution(
            result.solution,
            parameters,
            b_grid,
            shared,
        )
        result.scale, _ = population_scale(result.solution, parameters, outside)
        demand = (
            float(result.scale.scale_factor)
            * float(
                np.asarray(result.solution.housing_demand, dtype=float).reshape(-1)[0]
            )
        )
        supply = float(parameters.H0[0]) * (
            float(parameters.user_cost_rate)
            * float(result.price)
            / float(parameters.r_bar[0])
        ) ** float(parameters.xi_supply[0])
        result.residual = (demand - supply) / max(abs(supply), 1e-14)
    if abs(float(result.residual)) > tolerance:
        raise RuntimeError(
            f"Scaled price residual {result.residual:.3e} exceeds {tolerance:.3e}."
        )
    result.model_evaluations = len(cache)
    result.bracket_expansions = expansions
    return result


def main() -> None:
    args = parse_args()
    source = args.source.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    profile, theta, base_overrides, target_system, _ = load_profile(
        "new-moment",
        source,
        args.smoke,
    )
    metadata = {
        "status": "running",
        "profile": profile,
        "source": str(source),
        "theta": theta,
        "smoke": bool(args.smoke),
        "policy": "annual property tax increases from 1% to 2%",
        "fiscal_treatment": "property-tax revenue is not rebated; no grant",
        "entry_protocol": {
            "source": "Phase 9b node-level entry/scale closure",
            "target_baseline_city_entry_probability": TARGET_QBAR,
            "entry_taste_scale": ENTRY_TASTE_SCALE,
            "outside_objects": "recovered at the 1% no-rebate baseline and fixed",
        },
    }
    (outdir / "metadata.json").write_text(
        json.dumps(jsonable(metadata), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    started = time.perf_counter()
    print("baseline: 1% property tax, no rebate", flush=True)
    baseline_solution, baseline_parameters = solve_fixed_scale(
        base_overrides,
        0.01,
        args.smoke,
    )
    outside, entry_rows = calibrate_outside_objects(
        baseline_solution,
        baseline_parameters,
        target_qbar=TARGET_QBAR,
        taste_scale=ENTRY_TASTE_SCALE,
    )
    baseline_scale, _ = population_scale(
        baseline_solution,
        baseline_parameters,
        outside,
    )
    if abs(float(baseline_scale.scale_factor) - 1.0) > 1e-9:
        raise RuntimeError(
            f"Baseline entry identity failed: S={baseline_scale.scale_factor:.12g}."
        )
    write_csv(outdir / "benchmark_entry_nodes.csv", entry_rows)
    (outdir / "benchmark_outside_objects.json").write_text(
        json.dumps(jsonable(outside), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("counterfactual: 2% property tax, no rebate", flush=True)
    fixed_tax2_solution, tax2_parameters = solve_fixed_scale(
        base_overrides,
        0.02,
        args.smoke,
    )
    tax2 = solve_scaled_price(
        fixed_tax2_solution,
        tax2_parameters,
        outside,
        smoke=args.smoke,
    )

    cases = (
        ("tax1_no_rebate_baseline", baseline_solution, baseline_parameters, baseline_scale),
        ("tax2_no_rebate", tax2.solution, tax2_parameters, tax2.scale),
    )
    rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    for label, solution, parameters, scale in cases:
        moments = extract_moments(solution, parameters)
        rows.append(
            {
                "case": label,
                "annual_property_tax_rate": (
                    float(parameters.tau_H) / float(parameters.period_years)
                ),
                "lump_sum_transfer_period_units": float(
                    parameters.property_tax_lump_sum_transfer
                ),
                "tfr": float(moments["tfr"]),
                "price": float(np.asarray(solution.p_eq).reshape(-1)[0]),
                "population_scale": float(scale.scale_factor),
                "entry_probability": float(scale.qbar),
                "normalized_total_births": float(solution.total_births_kfe),
                "scaled_total_births": (
                    float(scale.scale_factor) * float(solution.total_births_kfe)
                ),
                "property_tax_revenue_discarded": float(
                    scale.scale_factor * solution.property_tax_revenue
                ),
                "market_residual": (
                    float(solution.best_max_abs_rel_excess)
                    if label == "tax1_no_rebate_baseline"
                    else float(tax2.residual)
                ),
            }
        )
        for name, target, weight in zip(
            target_system.moment_names,
            target_system.target_values,
            target_system.weights,
        ):
            model = float(moments[name])
            gap = model - float(target)
            target_rows.append(
                {
                    "case": label,
                    "moment": name,
                    "target": float(target),
                    "model": model,
                    "gap": gap,
                    "weight": float(weight),
                    "loss_contribution": float(weight) * gap * gap,
                }
            )

    baseline = rows[0]
    for row in rows:
        row["tfr_change_level"] = row["tfr"] - baseline["tfr"]
        row["tfr_change_percent"] = 100.0 * (
            row["tfr"] / baseline["tfr"] - 1.0
        )
        row["price_change_percent"] = 100.0 * (
            row["price"] / baseline["price"] - 1.0
        )
        row["population_change_percent"] = 100.0 * (
            row["population_scale"] / baseline["population_scale"] - 1.0
        )
        row["total_births_change_percent"] = 100.0 * (
            row["scaled_total_births"] / baseline["scaled_total_births"] - 1.0
        )
    write_csv(outdir / "summary.csv", rows)
    write_csv(outdir / "target_fit_long.csv", target_rows)

    policy = rows[1]
    lines = [
        "# Property-tax increase with no rebate and endogenous entry",
        "",
        "Current-M parameters. The annual property tax rises from 1% to 2%; "
        "all property-tax revenue is discarded in both cases, no grant is paid, "
        "and entry/population scale follows the established Phase-9b closure.",
        "",
        "| Case | TFR | TFR change | Population | Total births | House price |",
        "|---|---:|---:|---:|---:|---:|",
        (
            f"| 1% tax, no rebate | {baseline['tfr']:.4f} | — | 1.0000 | — | "
            f"{baseline['price']:.4f} |"
        ),
        (
            f"| 2% tax, no rebate | {policy['tfr']:.4f} | "
            f"{policy['tfr_change_level']:+.4f} "
            f"({policy['tfr_change_percent']:+.2f}%) | "
            f"{policy['population_scale']:.4f} "
            f"({policy['population_change_percent']:+.2f}%) | "
            f"{policy['total_births_change_percent']:+.2f}% | "
            f"{policy['price']:.4f} ({policy['price_change_percent']:+.2f}%) |"
        ),
    ]
    (outdir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    metadata.update(
        status="complete",
        elapsed_seconds=time.perf_counter() - started,
        model_evaluations_scaled_tax2=int(tax2.model_evaluations),
        scaled_price_bracket_expansions=int(tax2.bracket_expansions),
        benchmark_outside_objects=outside,
    )
    (outdir / "metadata.json").write_text(
        json.dumps(jsonable(metadata), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"done: TFR {policy['tfr_change_percent']:+.3f}%, "
        f"population {policy['population_change_percent']:+.3f}%, "
        f"total births {policy['total_births_change_percent']:+.3f}%, "
        f"price {policy['price_change_percent']:+.3f}%",
        flush=True,
    )


if __name__ == "__main__":
    main()
