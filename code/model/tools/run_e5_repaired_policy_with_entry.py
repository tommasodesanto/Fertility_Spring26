#!/usr/bin/env python3
"""Funded E5 policy counterfactual with genuine entry and population closure.

The baseline is normalized to population one.  At that baseline the driver
computes the city-entry probability required to match the empirical
outside-origin entrant share, then recovers the outside value and outside
entrant flow.  It holds those outside objects fixed and jointly solves the
counterfactual house price and lump-sum transfer.  Entry determines stationary
population scale; the government budget includes property-tax revenue,
targeted purchase grants, and the universal transfer.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

MODEL_DIR = Path(__file__).resolve().parents[1]
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

import run_intergen_funded_policy_with_entry as closure
from entry_target_contract import entry_target_from_outside_origin_share
from run_intergen_funded_property_tax_test import case_overrides, jsonable, write_csv


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = ROOT / "output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"
DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_e5_maturation_repair_policy_entry"
DEFAULT_OUTSIDE_ORIGIN_ENTRANT_SHARE = 0.169


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--outside-origin-entrant-share",
        type=float,
        default=DEFAULT_OUTSIDE_ORIGIN_ENTRANT_SHARE,
        help=(
            "Empirical share of baseline entrants originating outside the city. "
            "The model entry probability is computed from the baseline flow identity."
        ),
    )
    parser.add_argument("--entry-taste-scale", type=float, default=2.0)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def load_repaired_profile(source: Path, smoke: bool) -> tuple[dict[str, float], dict[str, Any], Any, str]:
    payload = json.loads(source.read_text(encoding="utf-8"))
    winner = (payload.get("winners") or {}).get("E1")
    if not isinstance(winner, dict) or not isinstance(winner.get("theta"), dict):
        raise ValueError(f"{source} does not contain winners.E1.theta")
    repeat = payload.get("winner_repeat_check") or {}
    if not (
        repeat.get("both_strict")
        and float(repeat.get("loss_abs_difference", math.inf)) == 0.0
        and float(repeat.get("max_abs_moment_difference", math.inf)) == 0.0
    ):
        raise ValueError(f"{source} does not preserve an exact strict repeat")
    theta = {str(key): float(value) for key, value in winner["theta"].items()}
    source_arm = str(
        (payload.get("winner_metadata") or {}).get("arm", winner.get("arm", "unknown"))
    )
    if source_arm not in {"E5_MATURATION_REPAIR", "E5F"}:
        raise ValueError(
            f"{source} reports arm {source_arm!r}; expected E5_MATURATION_REPAIR or E5F"
        )

    for name in ("E5_MATURATION_REPAIR", "E5F", "E6A", "E6B", "E6C"):
        os.environ.pop(name, None)
    os.environ["E3_L4"] = "1"
    os.environ["E5"] = "1"
    os.environ["E5_MATURATION_REPAIR"] = "1"
    if source_arm == "E5F":
        os.environ["E5F"] = "1"
    from intergen_eqscale_seq_optimized.calibration import extract_moments
    from intergen_eqscale_seq_optimized.e5_profile import e5_target_system
    from intergen_eqscale_seq_optimized import run_e1_chain
    from intergen_eqscale_seq_optimized.solver import (
        entry_wealth_grid_weights,
        income_transition_values,
        run_model_cp_dt,
    )

    run_e1_chain.load_runtime()
    runtime = SimpleNamespace(
        J=8 if smoke else 17,
        Nb=40 if smoke else 120,
        max_iter_eq=10 if smoke else 40,
        tol_eq=2.0e-4 if smoke else 2.5e-5,
    )
    overrides = run_e1_chain.common_overrides(runtime)
    overrides.update(theta)
    overrides.update(
        child_state_mode="independent_count",
        property_tax_lump_sum_transfer=0.0,
    )
    floor_is_active = bool(overrides.get("child_room_floor", False))
    if source_arm == "E5F" and not floor_is_active:
        raise RuntimeError("E5F policy routing failed to activate the child-room floor")
    if source_arm == "E5_MATURATION_REPAIR" and floor_is_active:
        raise RuntimeError("repaired-E5 control unexpectedly activates the child-room floor")

    # Reuse the already-audited closure algebra while routing every solve and
    # target readout through the repaired E package.
    closure.run_model_cp_dt = run_model_cp_dt
    closure.extract_moments = extract_moments
    closure.entry_wealth_grid_weights = entry_wealth_grid_weights
    closure.income_transition_values = income_transition_values
    return theta, overrides, e5_target_system(), source_arm


def fiscal_root(function, *, smoke: bool) -> float:
    lower, upper = 0.0, 0.5
    f_lower, f_upper = float(function(lower)), float(function(upper))
    while f_upper > 0.0 and upper < 8.0:
        upper *= 2.0
        f_upper = float(function(upper))
    if f_lower < 0.0 or f_upper > 0.0:
        raise RuntimeError(
            f"could not bracket fiscal transfer: F(0)={f_lower:.6g}, F({upper:g})={f_upper:.6g}"
        )
    tolerance = 2.0e-4 if smoke else 2.5e-5
    midpoint = 0.5 * (lower + upper)
    for _ in range(40):
        midpoint = 0.5 * (lower + upper)
        f_mid = float(function(midpoint))
        if abs(f_mid) <= tolerance:
            break
        if f_mid > 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return float(midpoint)


def solve_fixed_balanced_case(
    base_overrides: dict[str, Any],
    target_system: Any,
    *,
    label: str,
    annual_tax: float,
    grant: bool,
    smoke: bool,
) -> tuple[dict[str, Any], Any, Any]:
    cache: dict[float, tuple[Any, Any, np.ndarray]] = {}

    def solve_at(transfer: float) -> tuple[Any, Any, np.ndarray]:
        key = float(round(transfer, 12))
        if key not in cache:
            overrides = case_overrides(base_overrides, annual_tax, grant, transfer, smoke)
            cache[key] = closure.run_model_cp_dt(overrides, verbose=False)
        return cache[key]

    transfer = fiscal_root(
        lambda value: float(solve_at(value)[0].property_tax_budget_residual),
        smoke=smoke,
    )
    solution, parameters, price = solve_at(transfer)
    moments = closure.extract_moments(solution, parameters)
    targets = target_system.targets_dict()
    weights = target_system.weights_dict()
    loss = sum(weights[name] * (float(moments[name]) - targets[name]) ** 2 for name in targets)
    row = {
        "case": label,
        "annual_property_tax_rate": annual_tax,
        "purchase_grant": bool(grant),
        "lump_sum_transfer_period_units": transfer,
        "property_tax_revenue": float(solution.property_tax_revenue),
        "grant_outlays": float(solution.birth_entry_grant_outlays),
        "transfer_outlays": float(solution.property_tax_transfer_outlays),
        "government_budget_residual": float(solution.property_tax_budget_residual),
        "price": float(np.asarray(price).reshape(-1)[0]),
        "tfr": float(moments["tfr"]),
        "childless_rate": float(moments["childless_rate"]),
        "fixed_theta_loss": float(loss),
        "market_residual": float(solution.best_max_abs_rel_excess),
        "strict_converged": bool(getattr(solution, "timings", {}).get("strict_converged", False)),
        "fiscal_root_evaluations": len(cache),
    }
    return row, solution, parameters


def main() -> None:
    args = parse_args()
    if not 0.0 < args.outside_origin_entrant_share < 1.0:
        raise ValueError("outside-origin entrant share must lie strictly between zero and one")
    if not math.isfinite(args.entry_taste_scale) or args.entry_taste_scale <= 0.0:
        raise ValueError("entry taste scale must be finite and positive")
    source, outdir = args.source.resolve(), args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    theta, base_overrides, target_system, source_arm = load_repaired_profile(source, args.smoke)
    cases = (
        ("rebated_tax1_baseline", 0.01, False),
        ("rebated_tax2", 0.02, False),
        ("rebated_tax2_grant0p4_Hge6", 0.02, True),
    )
    metadata = {
        "status": "running",
        "source": str(source),
        "source_arm": source_arm,
        "theta": theta,
        "child_maturation": "independent binomial thinning; expected 18 years at home",
        "outside_origin_entrant_share_target": float(args.outside_origin_entrant_share),
        "entry_target_probability": "computed from the funded baseline flow identity",
        "entry_taste_scale": float(args.entry_taste_scale),
        "entry_taste_scale_status": "external sensitivity; empirical discipline outstanding",
        "local_birth_entry_weight_status": "fixed at one; empirical retention discipline outstanding",
        "geographic_interpretation_status": "local-versus-national policy interpretation outstanding",
        "outside_objects": "recovered at funded 1 percent baseline and fixed thereafter",
        "fiscal_rule": "property-tax revenue equals universal transfers plus targeted grants",
        "fixed_population_results": "decomposition only",
        "smoke": bool(args.smoke),
    }
    (outdir / "metadata.json").write_text(json.dumps(jsonable(metadata), indent=2) + "\n")

    fixed_rows: list[dict[str, Any]] = []
    fixed_solutions: dict[str, tuple[Any, Any]] = {}
    for label, annual_tax, grant in cases:
        print(f"fixed-population closure: {label}", flush=True)
        row, solution, parameters = solve_fixed_balanced_case(
            base_overrides,
            target_system,
            label=label,
            annual_tax=annual_tax,
            grant=grant,
            smoke=args.smoke,
        )
        fixed_rows.append(row)
        fixed_solutions[label] = (solution, parameters)
        write_csv(outdir / "fixed_population_summary.csv", fixed_rows)

    baseline_solution, baseline_parameters = fixed_solutions[cases[0][0]]
    entry_target = entry_target_from_outside_origin_share(
        baseline_solution,
        baseline_parameters,
        outside_origin_share=float(args.outside_origin_entrant_share),
    )
    outside, entry_rows = closure.calibrate_outside_objects(
        baseline_solution,
        baseline_parameters,
        target_qbar=float(entry_target["target_qbar"]),
        taste_scale=float(args.entry_taste_scale),
    )
    identity, _ = closure.population_scale(baseline_solution, baseline_parameters, outside)
    outside["identity_scale_factor"] = float(identity.scale_factor)
    outside["identity_entry_residual"] = float(identity.entry_residual)
    outside.update(entry_target)
    implied_outside_share = (
        float(identity.qbar)
        * float(outside["outside_entry_flow"])
        / float(entry_target["baseline_entry_flow"])
    )
    outside["outside_origin_entrant_share_recovered"] = implied_outside_share
    if abs(implied_outside_share - float(args.outside_origin_entrant_share)) > 1.0e-10:
        raise RuntimeError(
            "baseline outside-origin entrant share does not reproduce its empirical target: "
            f"target={args.outside_origin_entrant_share:.12g}, "
            f"recovered={implied_outside_share:.12g}"
        )
    (outdir / "benchmark_outside_objects.json").write_text(
        json.dumps(jsonable(outside), indent=2) + "\n"
    )
    write_csv(outdir / "benchmark_entry_nodes.csv", entry_rows)

    results: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    for (label, annual_tax, grant), fixed in zip(cases, fixed_rows):
        print(f"entry-adjusting closure: {label}", flush=True)
        joint = closure.solve_joint_case(
            base_overrides,
            outside,
            label=label,
            annual_tax=annual_tax,
            grant=grant,
            initial_price=float(fixed["price"]),
            initial_transfer=float(fixed["lump_sum_transfer_period_units"]),
            smoke=args.smoke,
        )
        moments = closure.extract_moments(joint.solution, joint.parameters)
        row = {
            "case": label,
            "annual_property_tax_rate": annual_tax,
            "purchase_grant": bool(grant),
            "tfr": float(moments["tfr"]),
            "childless_rate": float(moments["childless_rate"]),
            "price": float(joint.price),
            "population_scale": float(joint.scale.scale_factor),
            "entry_probability": float(joint.scale.qbar),
            "normalized_births": float(joint.solution.total_births_kfe),
            "total_births": float(joint.scale.scale_factor) * float(joint.solution.total_births_kfe),
            "lump_sum_transfer_period_units": float(joint.transfer),
            "scaled_property_tax_revenue": float(joint.scale.scale_factor) * float(joint.solution.property_tax_revenue),
            "scaled_grant_outlays": float(joint.scale.scale_factor) * float(joint.solution.birth_entry_grant_outlays),
            "scaled_transfer_outlays": float(joint.scale.scale_factor) * float(joint.solution.property_tax_transfer_outlays),
            "scaled_fiscal_residual": float(joint.scaled_fiscal_residual),
            "market_residual": float(joint.relative_excess),
            "joint_model_evaluations": int(joint.joint_model_evaluations),
        }
        results.append(row)
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
        write_csv(outdir / "summary.csv", results)
        write_csv(outdir / "target_fit_long.csv", target_rows)
        (outdir / "latest_completed_case.json").write_text(json.dumps(jsonable(row), indent=2) + "\n")

    baseline = results[0]
    for row in results:
        row["tfr_change_percent"] = 100.0 * (row["tfr"] / baseline["tfr"] - 1.0)
        row["price_change_percent"] = 100.0 * (row["price"] / baseline["price"] - 1.0)
        row["population_change_percent"] = 100.0 * (row["population_scale"] - 1.0)
        row["total_births_change_percent"] = 100.0 * (row["total_births"] / baseline["total_births"] - 1.0)
    write_csv(outdir / "summary.csv", results)
    metadata.update(
        status="complete",
        elapsed_seconds=time.perf_counter() - started,
        entry_target_probability=float(entry_target["target_qbar"]),
        outside_objects=outside,
    )
    (outdir / "metadata.json").write_text(json.dumps(jsonable(metadata), indent=2) + "\n")


if __name__ == "__main__":
    main()
