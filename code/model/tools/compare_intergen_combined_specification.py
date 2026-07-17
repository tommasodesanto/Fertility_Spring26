#!/usr/bin/env python3
"""Controlled four-arm test of bequest, finance, and housing-supply revisions."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np

from intergen_housing_fertility.calibration import (
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_BOUNDS,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    REPAIRED_TIMING_SWITCHES,
    production_profile_overrides,
)
from intergen_housing_fertility.solver import run_model_cp_dt


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_THETA = ROOT / "output/model/intergen_code_repair_20260709/active_theta_available_20260709.json"
DEFAULT_TARGETS = (
    ROOT
    / "code/data/mms_center_periphery/output_intergen_one_market_targets"
    / "intergen_one_market_acs_housing_targets.csv"
)
TARGET_NAME = "aggregate_mean_occupied_rooms_18_85"
NEW_ANNUAL_INTEREST = 0.02
NEW_ANNUAL_DEPRECIATION = 0.011
NEW_SUPPLY_ELASTICITY = 1.75


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theta-json", type=Path, default=DEFAULT_THETA)
    parser.add_argument("--housing-targets-csv", type=Path, default=DEFAULT_TARGETS)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ROOT / "output/model/intergen_combined_specification_20260710",
    )
    parser.add_argument("--mode", choices=("smoke", "fixed"), default="fixed")
    parser.add_argument(
        "--arm",
        choices=("all", "current", "bequest_normalization_only", "finance_supply_revisions_only", "all_revisions"),
        default="all",
    )
    parser.add_argument("--smoke-path", type=Path, default=None)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def load_theta(path: Path) -> dict[str, float]:
    payload = json.loads(path.read_text())
    raw = payload.get("theta", payload)
    theta = {str(key): float(value) for key, value in raw.items()}
    expected = {"beta" if name == "beta_annual" else name for name, _, _ in PRODUCTION_SEARCH_BOUNDS}
    if set(theta) != expected:
        raise ValueError(f"theta keys differ from production bounds: {sorted(set(theta) ^ expected)}")
    return theta


def load_housing_target(path: Path) -> dict[str, Any]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    matches = [row for row in rows if row.get("moment") == TARGET_NAME]
    if len(matches) != 1:
        raise ValueError(f"expected one {TARGET_NAME!r} row in {path}, found {len(matches)}")
    row = matches[0]
    value_key = next((key for key in ("value", "estimate", "mean") if key in row), None)
    if value_key is None:
        raise ValueError(f"cannot locate target value column in {list(row)}")
    value = float(row[value_key])
    if not (4.0 < value < 8.0):
        raise ValueError(f"implausible aggregate room target: {value}")
    return {"name": TARGET_NAME, "value": value, "source": str(path.resolve()), "row": row}


def finance_overrides(period_years: float) -> dict[str, Any]:
    converted = {
        "q": (1.0 + NEW_ANNUAL_INTEREST) ** period_years - 1.0,
        "delta": 1.0 - (1.0 - NEW_ANNUAL_DEPRECIATION) ** period_years,
    }
    if period_years == 4.0:
        assert abs(converted["q"] - 0.08243215999999998) <= 1e-14
        assert abs(converted["delta"] - 0.04327930935900004) <= 1e-14
    return converted


def common_overrides(theta: dict[str, float], *, J: int, Nb: int) -> dict[str, Any]:
    return {
        **base_overrides(J=J, Nb=Nb, n_house=5, max_iter_eq=PRODUCTION_MAX_ITER_EQ),
        **production_profile_overrides(),
        **REPAIRED_TIMING_SWITCHES,
        **theta,
    }


def target_fit_rows(moments: dict[str, float], targets: dict[str, float], weights: dict[str, float]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, target in targets.items():
        model = float(moments.get(name, math.nan))
        gap = model - float(target)
        weight = float(weights[name])
        rows.append({
            "moment": name,
            "target": float(target),
            "model": model,
            "gap": gap,
            "weight": weight,
            "loss_contribution": weight * gap * gap,
        })
    return rows


def parameter_rows(theta: dict[str, float]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, lower, upper in PRODUCTION_SEARCH_BOUNDS:
        key = "beta" if name == "beta_annual" else name
        estimate = float(theta[key]) ** 0.25 if name == "beta_annual" else float(theta[key])
        span = float(upper - lower)
        rows.append({
            "parameter": name,
            "estimate": estimate,
            "lower": float(lower),
            "upper": float(upper),
            "near_bound": min(estimate - float(lower), float(upper) - estimate) <= 0.02 * span,
        })
    return rows


def calibrate_h0(
    overrides: dict[str, Any],
    target_rooms: float,
    quiet: bool,
    checkpoint_path: Path | None = None,
) -> tuple[Any, Any, np.ndarray, list[dict[str, float]]]:
    """Calibrate H0 directly against aggregate rooms using a few full-GE solves."""
    evaluations: list[dict[str, float]] = []
    last_price: float | None = None

    def evaluate(log_h0: float) -> tuple[float, Any, Any, np.ndarray]:
        nonlocal last_price
        h0 = float(math.exp(log_h0))
        candidate = {**overrides, "H0": np.array([h0])}
        if last_price is not None:
            candidate["p_init_override"] = np.array([last_price])
        started = time.perf_counter()
        sol, P, p_eq = run_model_cp_dt(candidate, verbose=not quiet)
        strict, residual = strict_convergence(sol, P)
        rooms = float(sol.aggregate_housing_demand)
        last_price = float(np.asarray(p_eq).reshape(-1)[0])
        record = {
            "H0": h0,
            "price": last_price,
            "aggregate_rooms": rooms,
            "relative_target_gap": rooms / target_rooms - 1.0,
            "market_residual": residual,
            "strict_converged": strict,
            "elapsed_sec": float(time.perf_counter() - started),
        }
        evaluations.append(record)
        if checkpoint_path is not None:
            checkpoint_path.write_text(json.dumps(jsonable({
                "status": "calibrating_H0_against_aggregate_rooms",
                "target_rooms": target_rooms,
                "evaluations": evaluations,
            }), indent=2))
        if not strict or not math.isfinite(rooms) or rooms <= 0.0:
            raise RuntimeError(f"H0 calibration produced an invalid equilibrium: {record}")
        return math.log(rooms / target_rooms), sol, P, p_eq

    x0 = math.log(float(np.asarray(overrides.get("H0", [4.0])).reshape(-1)[0]))
    f0, sol0, P0, p0 = evaluate(x0)
    if abs(f0) <= 5e-4:
        return sol0, P0, p0, evaluations
    x1 = float(np.clip(x0 - f0, math.log(0.25), math.log(40.0)))
    f1, sol1, P1, p1 = evaluate(x1)
    for _ in range(4):
        if abs(f1) <= 5e-4:
            return sol1, P1, p1, evaluations
        slope_den = f1 - f0
        x2 = x1 - f1 if abs(slope_den) < 1e-8 else x1 - f1 * (x1 - x0) / slope_den
        x2 = float(np.clip(x2, x1 - math.log(2.0), x1 + math.log(2.0)))
        x2 = float(np.clip(x2, math.log(0.25), math.log(40.0)))
        f2, sol2, P2, p2 = evaluate(x2)
        x0, f0 = x1, f1
        x1, f1, sol1, P1, p1 = x2, f2, sol2, P2, p2
    if abs(f1) > 5e-4:
        raise RuntimeError(f"H0 calibration missed aggregate-rooms target after {len(evaluations)} GE solves")
    return sol1, P1, p1, evaluations


def strict_convergence(sol: Any, P: Any) -> tuple[bool, float]:
    residual = float(getattr(sol, "best_max_abs_rel_excess", math.nan))
    timings = dict(getattr(sol, "timings", {}))
    strict = bool(
        timings.get("strict_converged", getattr(sol, "converged", False))
        and math.isfinite(residual)
        and residual <= float(P.tol_eq)
    )
    return strict, residual


def run_arm(
    label: str,
    theta: dict[str, float],
    normalized: bool,
    revisions: bool,
    target_rooms: float,
    targets: dict[str, float],
    weights: dict[str, float],
    quiet: bool,
    checkpoint_path: Path | None = None,
) -> dict[str, Any]:
    overrides = common_overrides(theta, J=PRODUCTION_J, Nb=PRODUCTION_SEARCH_NB)
    overrides["normalize_bequest_utility"] = bool(normalized)
    inversion_evaluations: list[dict[str, float]] = []
    inferred = False
    started = time.perf_counter()
    if revisions:
        overrides.update(finance_overrides(float(overrides["period_years"])))
        overrides["eta_supply"] = np.array([NEW_SUPPLY_ELASTICITY])
        sol, P, p_eq, inversion_evaluations = calibrate_h0(
            overrides, target_rooms, quiet, checkpoint_path
        )
        inferred = True
    else:
        sol, P, p_eq = run_model_cp_dt(overrides, verbose=not quiet)
    elapsed = time.perf_counter() - started
    moments = extract_moments(sol, P)
    strict, residual = strict_convergence(sol, P)
    demand = float(sol.aggregate_housing_demand)
    supply = float(sol.aggregate_housing_supply)
    quantity_gap = demand - target_rooms if inferred else math.nan
    if inferred and abs(quantity_gap) > 0.01:
        strict = False
    return {
        "label": label,
        "normalize_bequest_utility": bool(normalized),
        "finance_supply_revisions": bool(revisions),
        "strict_converged": strict,
        "market_residual": residual,
        "rank_loss": float(diagnostic_loss(moments, targets=targets, weights=weights)),
        "elapsed_final_ge_sec": float(elapsed),
        "price": float(np.asarray(p_eq).reshape(-1)[0]),
        "owner_user_cost": float(np.asarray(sol.owner_user_cost).reshape(-1)[0]),
        "annual_interest": float((1.0 + P.q) ** (1.0 / P.period_years) - 1.0),
        "annual_depreciation": float(1.0 - (1.0 - P.delta) ** (1.0 / P.period_years)),
        "model_period_interest": float(P.q),
        "model_period_depreciation": float(P.delta),
        "model_period_property_tax": float(P.tau_H),
        "user_cost_rate": float(P.user_cost_rate),
        "supply_elasticity": float(np.asarray(P.xi_supply).reshape(-1)[0]),
        "H0": float(np.asarray(P.H0).reshape(-1)[0]),
        "H0_inferred": inferred,
        "inversion_price": float(np.asarray(p_eq).reshape(-1)[0]) if inferred else None,
        "aggregate_rooms_demand": demand,
        "aggregate_rooms_supply": supply,
        "aggregate_rooms_target": target_rooms if inferred else None,
        "aggregate_rooms_target_gap": quantity_gap,
        "theta": dict(theta),
        "moments": jsonable(moments),
        "target_fit": target_fit_rows(moments, targets, weights),
        "inversion_evaluations": inversion_evaluations,
        "timings": jsonable(dict(getattr(sol, "timings", {}))),
    }


def run_smoke(theta: dict[str, float], quiet: bool) -> dict[str, Any]:
    overrides = common_overrides(theta, J=5, Nb=16)
    overrides.update(finance_overrides(float(overrides["period_years"])))
    overrides["eta_supply"] = np.array([NEW_SUPPLY_ELASTICITY])
    overrides["H0"] = np.array([4.0])
    baseline, P0, p0 = run_model_cp_dt(overrides, verbose=not quiet)
    baseline_strict, baseline_residual = strict_convergence(baseline, P0)
    target = float(baseline.aggregate_housing_demand)
    recovered, P1, recovered_p, evaluations = calibrate_h0(overrides, target, quiet)
    recovered_strict, recovered_residual = strict_convergence(recovered, P1)
    recovered_rooms = float(recovered.aggregate_housing_demand)
    recovered_h0 = float(np.asarray(P1.H0).reshape(-1)[0])
    supply_at_inversion = float(recovered.aggregate_housing_supply)
    passed = bool(
        baseline_strict
        and recovered_strict
        and abs(recovered_h0 - 4.0) / 4.0 <= 1e-3
        and abs(float(np.asarray(recovered_p).reshape(-1)[0]) - float(np.asarray(p0).reshape(-1)[0])) / float(np.asarray(p0).reshape(-1)[0]) <= 1e-3
        and abs(supply_at_inversion - target) / target <= 6e-4
        and abs(recovered_rooms - target) / target <= 5e-4
    )
    return {
        "passed": passed,
        "baseline_strict_converged": baseline_strict,
        "baseline_market_residual": baseline_residual,
        "baseline_price": float(np.asarray(p0).reshape(-1)[0]),
        "baseline_target_rooms": target,
        "recovered_price": float(np.asarray(recovered_p).reshape(-1)[0]),
        "recovered_H0": recovered_h0,
        "supply_identity_value": supply_at_inversion,
        "final_ge_strict_converged": recovered_strict,
        "final_ge_market_residual": recovered_residual,
        "final_ge_price": float(np.asarray(recovered_p).reshape(-1)[0]),
        "final_ge_rooms": recovered_rooms,
        "evaluations": evaluations,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    theta = load_theta(args.theta_json)
    if args.mode == "smoke":
        smoke = run_smoke(theta, args.quiet)
        (args.outdir / "inversion_smoke.json").write_text(json.dumps(jsonable(smoke), indent=2))
        print(json.dumps(jsonable(smoke), indent=2))
        if not smoke["passed"]:
            raise RuntimeError("H0 inversion smoke test failed")
        return

    housing_target = load_housing_target(args.housing_targets_csv)
    smoke_path = args.smoke_path or (args.outdir / "inversion_smoke.json")
    if not smoke_path.exists() or not json.loads(smoke_path.read_text()).get("passed", False):
        raise RuntimeError("run and pass --mode smoke in this output folder before fixed arms")
    targets, weights = get_target_set(PRODUCTION_TARGET_SET)
    designs = (
        ("current", False, False),
        ("bequest_normalization_only", True, False),
        ("finance_supply_revisions_only", False, True),
        ("all_revisions", True, True),
    )
    if args.arm != "all":
        designs = tuple(design for design in designs if design[0] == args.arm)
    arms: list[dict[str, Any]] = []
    for label, normalized, revisions in designs:
        arm = run_arm(
            label,
            theta,
            normalized,
            revisions,
            housing_target["value"],
            targets,
            weights,
            args.quiet,
            args.outdir / f"{label}_inversion_checkpoint.json" if revisions else None,
        )
        arms.append(arm)
        (args.outdir / "latest_completed_arm.json").write_text(
            json.dumps(jsonable(arm), indent=2, sort_keys=True)
        )
        (args.outdir / "four_arm_progress.json").write_text(
            json.dumps(jsonable({"completed": arms}), indent=2, sort_keys=True)
        )
    if not all(arm["strict_converged"] for arm in arms):
        raise RuntimeError("at least one fixed-theta arm failed strict convergence or the quantity target")
    if any(arm["theta"] != arms[0]["theta"] for arm in arms[1:]):
        raise RuntimeError("fixed-theta vectors differ across arms")
    payload = {
        "status": "fixed_theta_four_arm_specification_test",
        "specification_adopted": False,
        "production_profile": PRODUCTION_PROFILE_NAME,
        "target_set": PRODUCTION_TARGET_SET,
        "Nb": PRODUCTION_SEARCH_NB,
        "theta_source": str(args.theta_json.resolve()),
        "housing_target": housing_target,
        "parameters_identical_across_arms": True,
        "parameters": parameter_rows(theta),
        "arms": arms,
    }
    write_csv(args.outdir / "four_arm_target_fit.csv", [
        {"case": arm["label"], **row} for arm in arms for row in arm["target_fit"]
    ])
    write_csv(args.outdir / "fixed_theta_parameters.csv", parameter_rows(theta))
    (args.outdir / "four_arm_comparison.json").write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True)
    )
    print(json.dumps({
        "outdir": str(args.outdir),
        "arms": [{"label": a["label"], "loss": a["rank_loss"], "price": a["price"], "H0": a["H0"]} for a in arms],
    }, indent=2))


if __name__ == "__main__":
    main()
