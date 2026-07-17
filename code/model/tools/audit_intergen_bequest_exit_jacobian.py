#!/usr/bin/env python3
"""Fresh tight-evaluator Jacobian at an A3, M1, M3, M4, or M5 calibration winner.

The diagnostic perturbs every active coordinate in the same bounded transformed
coordinates used by the search. It checkpoints every solve and reports ranks
under both target scaling and the SMM weights, leave-one-moment-out ranks and
conditioning, weak singular-vector loadings, and deterministic-repeat noise.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from intergen_housing_fertility.calibration import diagnostic_loss, extract_moments
from intergen_housing_fertility.local_panel import jsonable
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    validate_production_profile,
)
from intergen_housing_fertility.solver import InfeasibleThetaError, run_model_cp_dt
from tools.run_intergen_bequest_exit_chain import (
    BASE_DOMAIN,
    M4_THETA0_DOMAIN,
    M4_THETA1_DOMAIN,
    M5_KAPPA_DOMAIN,
    THETA0_DOMAIN,
    THETA1_DOMAIN,
    THETA_N_DOMAIN,
    arm_contract,
    common_overrides,
    load_theta,
    parameter_table,
    target_fit,
    target_system,
    theta_from_unit,
    unit_from_theta,
)


ACTIVE_DOMAIN = list(BASE_DOMAIN) + [THETA0_DOMAIN, THETA1_DOMAIN, THETA_N_DOMAIN]
M5_WINNER_JSON_DEFAULT = (
    Path(__file__).resolve().parents[3]
    / "output/model/intergen_income_disciplined_recalibration_20260716/report/results.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--winner-json",
        type=Path,
        default=None,
        help="Winner record; defaults to the M5 collector report for --arm M5.",
    )
    parser.add_argument("--winner-arm", type=str, default=None)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--arm", choices=("A3", "M1", "M3", "M4", "M5"), default="A3")
    parser.add_argument("--ltv-terminal", type=float, default=0.4)
    parser.add_argument("--theta1", type=float, default=0.25)
    parser.add_argument("--unit-step", type=float, default=0.005)
    parser.add_argument(
        "--parameter",
        action="append",
        choices=[row[0] for row in ACTIVE_DOMAIN + [M5_KAPPA_DOMAIN]],
        help="Restrict the perturbation list; repeatable. Production default is every active coordinate.",
    )
    parser.add_argument("--J", type=int, default=PRODUCTION_J)
    parser.add_argument("--Nb", type=int, default=PRODUCTION_SEARCH_NB)
    parser.add_argument("--max-iter-eq", type=int, default=40)
    parser.add_argument("--tol-eq", type=float, default=2.5e-5)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.winner_json is None:
        if args.arm != "M5":
            parser.error(f"--winner-json is required for --arm {args.arm}")
        args.winner_json = M5_WINNER_JSON_DEFAULT
        if args.winner_arm is None:
            args.winner_arm = "M5"
    return args


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(jsonable(value), indent=2, sort_keys=True))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def singular_summary(matrix: np.ndarray, label: str) -> dict[str, Any]:
    singular = np.linalg.svd(matrix, compute_uv=False)
    top = float(singular[0]) if singular.size else 0.0
    relative = singular / top if top > 0.0 else np.zeros_like(singular)
    smallest = float(singular[-1]) if singular.size else 0.0
    return {
        "scaling": label,
        "singular_values": singular.tolist(),
        "relative_singular_values": relative.tolist(),
        "rank_rel_1e_2": int(np.sum(relative >= 1e-2)),
        "rank_rel_1e_3": int(np.sum(relative >= 1e-3)),
        "condition_number": float(top / smallest) if smallest > 0.0 else math.inf,
        "matrix_shape": list(matrix.shape),
    }


def weakest_direction_rows(
    matrix: np.ndarray, parameters: list[str], scaling: str
) -> list[dict[str, Any]]:
    _, singular, vh = np.linalg.svd(matrix, full_matrices=False)
    vector = vh[-1] if vh.size else np.zeros(len(parameters))
    rows = [
        {
            "scaling": scaling,
            "parameter": parameter,
            "coefficient": float(coefficient),
            "absolute_coefficient": abs(float(coefficient)),
            "smallest_singular_value": float(singular[-1]) if singular.size else 0.0,
        }
        for parameter, coefficient in zip(parameters, vector)
    ]
    return sorted(rows, key=lambda row: row["absolute_coefficient"], reverse=True)


def main() -> None:
    args = parse_args()
    if args.smoke:
        args.Nb = 60
        args.max_iter_eq = 2
        args.tol_eq = 0.25
    else:
        validate_production_profile(
            PRODUCTION_PROFILE_NAME,
            J=args.J,
            Nb=args.Nb,
            n_house=5,
            income_states=5,
            target_set=PRODUCTION_TARGET_SET,
            # Validate the structural profile, then enforce the tighter
            # evaluator used for the winner audit below.
            max_iter_eq=PRODUCTION_MAX_ITER_EQ,
            stage="search",
        )
        if int(args.max_iter_eq) != 40 or not math.isclose(float(args.tol_eq), 2.5e-5):
            raise ValueError(
                "production Jacobian requires max_iter_eq=40 and tol_eq=2.5e-5"
            )

    step = float(args.unit_step)
    if not 1e-5 <= step <= 0.05:
        raise ValueError("unit-step must lie in [1e-5, 0.05]")
    contract_args = SimpleNamespace(
        arm=args.arm,
        ltv_terminal=float(args.ltv_terminal),
        theta1=float(args.theta1),
        seed_theta0=0.30,
        seed_kappa=0.0,
        fixed_theta0=None,
        J=int(args.J),
        Nb=int(args.Nb),
        max_iter_eq=int(args.max_iter_eq),
        tol_eq=float(args.tol_eq),
    )
    active, fixed, mechanism = arm_contract(contract_args)
    if args.arm == "M3":
        expected_active = ACTIVE_DOMAIN
    elif args.arm == "M4":
        expected_active = list(BASE_DOMAIN) + [M4_THETA0_DOMAIN, M4_THETA1_DOMAIN]
    elif args.arm == "M5":
        expected_active = list(BASE_DOMAIN) + [
            M4_THETA0_DOMAIN,
            M4_THETA1_DOMAIN,
            M5_KAPPA_DOMAIN,
        ]
    elif args.arm == "A3":
        expected_active = list(BASE_DOMAIN) + [THETA0_DOMAIN]
    else:
        expected_active = list(BASE_DOMAIN)
    if active != expected_active:
        raise RuntimeError(f"{args.arm} active domain drifted from its expected coordinates")
    selected = list(args.parameter or [row[0] for row in active])
    invalid = sorted(set(selected) - {row[0] for row in active})
    if invalid:
        raise ValueError(f"parameters are not active in {args.arm}: {invalid}")
    if not args.smoke and selected != [row[0] for row in active]:
        raise ValueError("production Jacobian must perturb every active parameter")

    theta = load_theta(args.winner_json, args.winner_arm)
    theta.update(fixed)
    x0 = unit_from_theta(theta, active)
    if args.arm == "M3":
        target_set = "candidate_replacement_bequest_internal_v1"
    elif args.arm == "M4":
        target_set = "candidate_replacement_bequest_median_composition_v1"
    elif args.arm == "M5":
        target_set = "candidate_replacement_income_disciplined_v1"
    else:
        target_set = PRODUCTION_TARGET_SET
    targets, weights = target_system(target_set)
    names = list(targets)
    overrides = common_overrides(contract_args, mechanism)

    args.outdir.mkdir(parents=True, exist_ok=True)
    cases_path = args.outdir / "cases.jsonl"
    cases_path.write_text("")
    metadata = {
        "status": "smoke" if args.smoke else f"fresh_tight_{args.arm.lower()}_identification_jacobian",
        "arm": args.arm,
        "winner_json": str(args.winner_json),
        "winner_arm": args.winner_arm,
        "active_domain": [
            {"name": name, "lower": lo, "upper": hi, "transform": transform}
            for name, lo, hi, transform in active
        ],
        "parameters_perturbed": selected,
        "parameter_count": len(selected),
        "moment_count": len(names),
        "moment_order": names,
        "target_set": target_set,
        "unit_step": step,
        "fixed_parameters": fixed,
        "mechanism": mechanism,
        "J": int(args.J),
        "Nb": int(args.Nb),
        "max_iter_eq": int(args.max_iter_eq),
        "tol_eq": float(args.tol_eq),
        "target_scaling": "d moment / max(abs(target), 0.1) per unit-coordinate change",
        "smm_scaling": "sqrt(weight) * d moment per unit-coordinate change",
    }
    write_json(args.outdir / "metadata.json", metadata)

    records: list[dict[str, Any]] = []
    started = time.perf_counter()

    def solve(
        unit: np.ndarray,
        label: str,
        perturbation: dict[str, Any],
        *,
        require_strict: bool = True,
    ) -> dict[str, Any]:
        proposal = theta_from_unit(np.clip(unit, 0.0, 1.0), active, fixed)
        t0 = time.perf_counter()
        try:
            sol, P, p_eq = run_model_cp_dt({**overrides, **proposal}, verbose=False)
            moments = extract_moments(sol, P)
            residual = float(getattr(sol, "best_max_abs_rel_excess", math.inf))
            timings = dict(getattr(sol, "timings", {}))
            strict = bool(
                timings.get("strict_converged", getattr(sol, "converged", False))
                and math.isfinite(residual)
                and residual <= float(P.tol_eq)
            )
            status, error = "ok", ""
            price = float(np.asarray(p_eq).reshape(-1)[0])
        except InfeasibleThetaError as exc:
            moments, timings = {}, {}
            residual, price, strict = math.inf, math.nan, False
            status, error = "infeasible_theta", str(exc)
        except Exception as exc:  # noqa: BLE001 - persist a diagnostic failure.
            moments, timings = {}, {}
            residual, price, strict = math.inf, math.nan, False
            status, error = f"failed:{type(exc).__name__}", str(exc)
        record = {
            "case": len(records),
            "label": label,
            "status": status,
            "strict_converged": strict,
            "loss": float(diagnostic_loss(moments, targets=targets, weights=weights)) if moments else math.inf,
            "market_residual": residual,
            "price": price,
            "theta": proposal,
            "parameters": parameter_table(proposal, active, fixed),
            "moments": moments,
            "target_fit": target_fit(moments, targets, weights) if moments else [],
            "timings": timings,
            "perturbation": perturbation,
            "error": error,
            "elapsed_sec": time.perf_counter() - t0,
        }
        with cases_path.open("a") as handle:
            handle.write(json.dumps(jsonable(record), sort_keys=True) + "\n")
        records.append(record)
        write_json(args.outdir / "latest_completed_case.json", record)
        print(
            f"case {len(records)} {label}: strict={strict} loss={record['loss']:.7g} "
            f"resid={residual:.2e} elapsed={(time.perf_counter() - started) / 60:.1f}m",
            flush=True,
        )
        if require_strict and not strict:
            raise RuntimeError(f"Jacobian case {label} did not converge tightly: {status}: {error}")
        return record

    baseline = solve(x0, "baseline", {"side": "baseline"})
    repeat = solve(x0, "baseline_repeat", {"side": "determinism_repeat"})
    base_vector = np.asarray([baseline["moments"][name] for name in names], dtype=float)
    repeat_vector = np.asarray([repeat["moments"][name] for name in names], dtype=float)

    raw_columns: list[np.ndarray] = []
    long_rows: list[dict[str, Any]] = []
    effective_steps: dict[str, float] = {}
    for parameter in selected:
        idx = [row[0] for row in active].index(parameter)
        derivative = None
        method = ""
        denominator = math.nan
        for attempt in range(3):
            delta = step / (2**attempt)
            minus = plus = None
            if x0[idx] - delta >= 0.0:
                minus_unit = x0.copy()
                minus_unit[idx] -= delta
                minus = solve(
                    minus_unit,
                    f"{parameter}_minus_a{attempt}",
                    {"parameter": parameter, "side": "minus", "unit_step": -delta, "attempt": attempt},
                    require_strict=False,
                )
            if x0[idx] + delta <= 1.0:
                plus_unit = x0.copy()
                plus_unit[idx] += delta
                plus = solve(
                    plus_unit,
                    f"{parameter}_plus_a{attempt}",
                    {"parameter": parameter, "side": "plus", "unit_step": delta, "attempt": attempt},
                    require_strict=False,
                )
            minus_ok = bool(minus and minus.get("strict_converged", False))
            plus_ok = bool(plus and plus.get("strict_converged", False))
            if minus_ok and plus_ok:
                minus_vector = np.asarray([minus["moments"][name] for name in names])
                plus_vector = np.asarray([plus["moments"][name] for name in names])
                derivative = (plus_vector - minus_vector) / (2.0 * delta)
                method, denominator = "central", 2.0 * delta
                break
            if plus_ok:
                plus_vector = np.asarray([plus["moments"][name] for name in names])
                derivative = (plus_vector - base_vector) / delta
                method, denominator = "forward_strict_side", delta
                break
            if minus_ok:
                minus_vector = np.asarray([minus["moments"][name] for name in names])
                derivative = (base_vector - minus_vector) / delta
                method, denominator = "backward_strict_side", delta
                break
        if derivative is None:
            raise RuntimeError(
                f"no tightly converged finite-difference side for {parameter} "
                f"after steps {step}, {step / 2}, {step / 4}"
            )
        effective_steps[parameter] = denominator
        raw_columns.append(derivative)
        for moment, value in zip(names, derivative):
            long_rows.append(
                {
                    "moment": moment,
                    "parameter": parameter,
                    "derivative_per_unit_coordinate": float(value),
                    "method": method,
                    "unit_denominator": denominator,
                }
            )

    raw = np.column_stack(raw_columns)
    target_scale = np.asarray([max(abs(float(targets[name])), 0.1) for name in names])
    sqrt_weight = np.sqrt(np.asarray([float(weights[name]) for name in names]))
    target_scaled = raw / target_scale[:, None]
    smm_scaled = raw * sqrt_weight[:, None]
    target_svd = singular_summary(target_scaled, "target_scaled")
    smm_svd = singular_summary(smm_scaled, "smm_weighted")
    leave_one_out: list[dict[str, Any]] = []
    for idx, moment in enumerate(names):
        for label, matrix in (("target_scaled", target_scaled), ("smm_weighted", smm_scaled)):
            dropped = np.delete(matrix, idx, axis=0)
            svd = singular_summary(dropped, label)
            leave_one_out.append(
                {
                    "dropped_moment": moment,
                    "scaling": label,
                    "rank_rel_1e_2": svd["rank_rel_1e_2"],
                    "rank_rel_1e_3": svd["rank_rel_1e_3"],
                    "condition_number": svd["condition_number"],
                    "smallest_singular_value": svd["singular_values"][-1],
                    "smallest_relative_singular_value": svd["relative_singular_values"][-1],
                }
            )
    column_rows: list[dict[str, Any]] = []
    for col, parameter in enumerate(selected):
        for label, matrix in (("target_scaled", target_scaled), ("smm_weighted", smm_scaled)):
            values = matrix[:, col]
            order = np.argsort(np.abs(values))[::-1]
            column_rows.append(
                {
                    "parameter": parameter,
                    "scaling": label,
                    "column_norm": float(np.linalg.norm(values)),
                    "top_moment_1": names[int(order[0])],
                    "top_abs_loading_1": abs(float(values[int(order[0])])),
                    "top_moment_2": names[int(order[1])],
                    "top_abs_loading_2": abs(float(values[int(order[1])])),
                    "top_moment_3": names[int(order[2])],
                    "top_abs_loading_3": abs(float(values[int(order[2])])),
                }
            )
    weak_rows = weakest_direction_rows(target_scaled, selected, "target_scaled")
    weak_rows.extend(weakest_direction_rows(smm_scaled, selected, "smm_weighted"))
    repeat_delta = (repeat_vector - base_vector) / target_scale
    theta_idx = selected.index("theta0") if "theta0" in selected else None
    theta_norm = float(np.linalg.norm(target_scaled[:, theta_idx])) if theta_idx is not None else math.nan
    theta_step = effective_steps.get("theta0", step)
    noise_equivalent = float(np.linalg.norm(repeat_delta) / max(theta_step, 1e-12))
    noise_threshold = max(10.0 * noise_equivalent, 1e-8)

    write_csv(args.outdir / "jacobian_long.csv", long_rows)
    write_csv(args.outdir / "leave_one_moment_out.csv", leave_one_out)
    write_csv(args.outdir / "parameter_column_summary.csv", column_rows)
    write_csv(args.outdir / "weakest_direction.csv", weak_rows)
    write_json(
        args.outdir / "jacobian_matrices.json",
        {
            "moments": names,
            "parameters": selected,
            "raw_per_unit_coordinate": raw,
            "target_scaled": target_scaled,
            "smm_weighted": smm_scaled,
        },
    )
    summary = {
        "status": "complete",
        "metadata": metadata,
        "baseline_loss": baseline["loss"],
        "baseline_target_fit": baseline["target_fit"],
        "baseline_parameters": baseline["parameters"],
        "determinism_repeat": {
            "max_abs_moment_difference": float(np.max(np.abs(repeat_vector - base_vector))),
            "target_scaled_difference_norm": float(np.linalg.norm(repeat_delta)),
        },
        "rank": {"target_scaled": target_svd, "smm_weighted": smm_svd},
        "leave_one_moment_out": leave_one_out,
        "parameter_column_summary": column_rows,
        "weakest_direction": weak_rows,
        "theta0_column": {
            "target_scaled_norm": theta_norm,
            "repeat_noise_derivative_equivalent_norm": noise_equivalent,
            "noise_threshold": noise_threshold,
            "above_noise_floor": bool(theta_norm > noise_threshold) if theta_idx is not None else None,
        },
        "elapsed_sec": time.perf_counter() - started,
    }
    write_json(args.outdir / "summary.json", summary)
    print(json.dumps(jsonable(summary["rank"]), indent=2), flush=True)


if __name__ == "__main__":
    main()
