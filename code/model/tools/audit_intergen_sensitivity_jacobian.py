#!/usr/bin/env python3
"""Finite-difference sensitivity audit for the one-market intergen model.

This is an identification diagnostic, not a calibration routine. It takes
saved solved parameter vectors from the June frontier output, re-solves them
under a fixed diagnostic stack, perturbs the 13 internal SMM parameters, and
writes a normalized moment-by-parameter Jacobian with rank and collinearity
summaries.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_housing_fertility.calibration import (  # noqa: E402
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_housing_fertility.local_panel import (  # noqa: E402
    GLOBAL_DE_BOUNDS,
    income_process_overrides,
)
from intergen_housing_fertility.solver import run_model_cp_dt  # noqa: E402


DEFAULT_FRONTIER_ROOT = ROOT / "output/model/cluster_pulls/intergen_overnight_frontier_20260617"
DEFAULT_OUTDIR = ROOT / "output/model/intergen_sensitivity_jacobian_20260618"
DEFAULT_POINT_LABELS = ("core_feasibility_v1", "roomcost_test_v1")
TARGET_SET = "candidate_no_timing_v0"

PARAMETERS = [
    "beta",
    "alpha_cons",
    "b_entry_fixed",
    "c_bar_0",
    "c_bar_n",
    "h_bar_0",
    "h_bar_jump",
    "h_bar_n",
    "psi_child",
    "kappa_fert",
    "chi",
    "theta0",
    "theta_n",
]

EXTRA_MOMENTS = [
    "mean_age_first_birth",
    "parity_share_0",
    "parity_share_1",
    "parity_share_2plus",
    "aggregate_own_rate",
    "own_rate_2534",
    "own_rate_3544",
    "young_owner_rate",
    "old_owner_rate",
    "prime_childless_owner_minus_renter_rooms",
    "renter25_45_all_cap_share",
    "owner_neg_liquid_share_2534",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frontier-root", type=Path, default=DEFAULT_FRONTIER_ROOT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--point-label", action="append", choices=all_point_labels(), help="Frontier family to audit; repeatable")
    parser.add_argument("--target-set", default=TARGET_SET)
    parser.add_argument("--rel-step", type=float, default=0.01)
    parser.add_argument("--min-step", type=float, default=1e-4)
    parser.add_argument("--J", type=int, default=16)
    parser.add_argument("--Nb", type=int, default=60)
    parser.add_argument("--n-house", type=int, default=5)
    parser.add_argument("--max-iter-eq", type=int, default=3)
    parser.add_argument("--income-states", type=int, default=5)
    parser.add_argument("--rank-tol", type=float, default=1e-8)
    parser.add_argument("--max-points", type=int, default=0, help="Optional cap after point selection; 0 means no cap")
    args = parser.parse_args()

    os.environ.setdefault("NUMBA_NUM_THREADS", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

    point_labels = tuple(args.point_label) if args.point_label else DEFAULT_POINT_LABELS
    if int(args.max_points) > 0:
        point_labels = point_labels[: int(args.max_points)]

    targets, weights = get_target_set(args.target_set)
    target_moments = list(targets)
    moment_keys = unique_list(target_moments + EXTRA_MOMENTS)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    config = {
        "status": "finite_difference_identification_diagnostic_not_calibration",
        "frontier_root": str(args.frontier_root),
        "outdir": str(outdir),
        "point_labels": list(point_labels),
        "target_set": str(args.target_set),
        "target_moments": target_moments,
        "parameters": PARAMETERS,
        "rel_step": float(args.rel_step),
        "min_step": float(args.min_step),
        "J": int(args.J),
        "Nb": int(args.Nb),
        "n_house": int(args.n_house),
        "max_iter_eq": int(args.max_iter_eq),
        "income_states": int(args.income_states),
        "rank_tol": float(args.rank_tol),
    }
    write_json(outdir / "audit_config.json", config)

    points = [load_frontier_best(args.frontier_root, label) for label in point_labels]
    write_json(outdir / "source_points.json", {"points": points})

    all_baseline_rows: list[dict[str, Any]] = []
    all_case_rows: list[dict[str, Any]] = []
    all_jacobian_rows: list[dict[str, Any]] = []
    all_rank_rows: list[dict[str, Any]] = []
    all_collinearity_rows: list[dict[str, Any]] = []
    all_moment_top_rows: list[dict[str, Any]] = []
    all_param_top_rows: list[dict[str, Any]] = []
    summary_points: list[dict[str, Any]] = []

    for point in points:
        point_summary = run_point_audit(
            point,
            targets=targets,
            weights=weights,
            target_moments=target_moments,
            moment_keys=moment_keys,
            args=args,
            outdir=outdir,
        )
        all_baseline_rows.extend(point_summary["baseline_rows"])
        all_case_rows.extend(point_summary["case_rows"])
        all_jacobian_rows.extend(point_summary["jacobian_rows"])
        all_rank_rows.extend(point_summary["rank_rows"])
        all_collinearity_rows.extend(point_summary["collinearity_rows"])
        all_moment_top_rows.extend(point_summary["moment_top_rows"])
        all_param_top_rows.extend(point_summary["param_top_rows"])
        summary_points.append(point_summary["summary"])
        write_json(outdir / "latest_completed_point.json", point_summary["summary"])
        print(
            f"completed point {point_summary['summary']['point_label']}: "
            f"rank={point_summary['summary']['scaled_rank']}, "
            f"condition={point_summary['summary']['condition_number']:.3g}, "
            f"elapsed={point_summary['summary']['elapsed_sec']:.1f}s",
            flush=True,
        )

    write_csv(outdir / "baseline_moments.csv", all_baseline_rows)
    write_csv(outdir / "finite_difference_cases.csv", all_case_rows)
    write_csv(outdir / "jacobian_long.csv", all_jacobian_rows)
    write_csv(outdir / "scaled_singular_values.csv", all_rank_rows)
    write_csv(outdir / "scaled_column_correlations.csv", all_collinearity_rows)
    write_csv(outdir / "moment_top_sensitivities.csv", all_moment_top_rows)
    write_csv(outdir / "parameter_top_moments.csv", all_param_top_rows)

    summary = {
        "config": config,
        "points": summary_points,
        "interpretation": (
            "Rank and condition numbers are computed from the target-moment "
            "Jacobian scaled by max(abs(parameter), 1) and max(abs(moment), "
            "abs(target), 1). Discrete median-room moments can have locally "
            "zero finite differences even when they matter globally."
        ),
    }
    write_json(outdir / "audit_summary.json", summary)
    write_readme(outdir, summary, all_moment_top_rows, all_param_top_rows, all_collinearity_rows)


def run_point_audit(
    point: dict[str, Any],
    *,
    targets: dict[str, float],
    weights: dict[str, float],
    target_moments: list[str],
    moment_keys: list[str],
    args: argparse.Namespace,
    outdir: Path,
) -> dict[str, Any]:
    t_start = time.perf_counter()
    label = str(point["point_label"])
    theta0 = {k: float(v) for k, v in dict(point["theta"]).items() if k in PARAMETERS}
    point_dir = outdir / label
    point_dir.mkdir(parents=True, exist_ok=True)

    print(f"solving baseline for {label}", flush=True)
    base_record = solve_theta(
        label=label,
        case="baseline",
        theta=theta0,
        args=args,
        targets=targets,
        weights=weights,
        moment_keys=moment_keys,
    )
    write_json(point_dir / "baseline_record.json", base_record)

    baseline_rows = baseline_moment_rows(label, base_record, targets, target_moments)
    case_rows = [case_summary_row(label, base_record, parameter="", side="baseline", step=0.0)]

    plus_minus: dict[str, dict[str, dict[str, Any]]] = {}
    for parameter in PARAMETERS:
        plus_theta, minus_theta, step_info = perturbed_thetas(theta0, parameter, args.rel_step, args.min_step)
        plus_minus[parameter] = {}
        for side, theta in [("minus", minus_theta), ("plus", plus_theta)]:
            print(f"{label}: {parameter} {side}", flush=True)
            record = solve_theta(
                label=label,
                case=f"{parameter}_{side}",
                theta=theta,
                args=args,
                targets=targets,
                weights=weights,
                moment_keys=moment_keys,
            )
            record["perturbation"] = {
                "parameter": parameter,
                "side": side,
                **step_info,
                "theta_value": float(theta[parameter]),
            }
            plus_minus[parameter][side] = record
            case_rows.append(case_summary_row(label, record, parameter=parameter, side=side, step=step_info[side + "_step"]))
            write_json(point_dir / f"{parameter}_{side}.json", record)

    jacobian_rows, matrix_payload = build_jacobian_rows(
        label,
        theta0,
        base_record,
        plus_minus,
        targets,
        target_moments,
    )
    write_csv(point_dir / "jacobian_long.csv", jacobian_rows)
    write_json(point_dir / "jacobian_matrix.json", matrix_payload["json"])

    rank_rows, sv_summary = rank_summary(label, matrix_payload["scaled_matrix"], target_moments, PARAMETERS, args.rank_tol)
    collinearity_rows = column_collinearity_rows(label, matrix_payload["scaled_matrix"], PARAMETERS)
    moment_top_rows = moment_top_sensitivity_rows(label, matrix_payload["scaled_matrix"], target_moments, PARAMETERS)
    param_top_rows = parameter_top_moment_rows(label, matrix_payload["scaled_matrix"], target_moments, PARAMETERS)

    write_csv(point_dir / "scaled_singular_values.csv", rank_rows)
    write_csv(point_dir / "scaled_column_correlations.csv", collinearity_rows)
    write_csv(point_dir / "moment_top_sensitivities.csv", moment_top_rows)
    write_csv(point_dir / "parameter_top_moments.csv", param_top_rows)

    summary = {
        "point_label": label,
        "source_target_set": point.get("source_target_set"),
        "source_rank_loss": point.get("source_rank_loss"),
        "source_path": point.get("source_path"),
        "baseline_rank_loss": base_record["rank_loss"],
        "baseline_market_residual": base_record["market_residual"],
        "baseline_p_eq": base_record["p_eq"],
        "baseline_moments": base_record["moments"],
        "scaled_rank": sv_summary["rank"],
        "condition_number": sv_summary["condition_number"],
        "singular_values": sv_summary["singular_values"],
        "elapsed_sec": float(time.perf_counter() - t_start),
    }
    write_json(point_dir / "point_summary.json", summary)
    return {
        "summary": summary,
        "baseline_rows": baseline_rows,
        "case_rows": case_rows,
        "jacobian_rows": jacobian_rows,
        "rank_rows": rank_rows,
        "collinearity_rows": collinearity_rows,
        "moment_top_rows": moment_top_rows,
        "param_top_rows": param_top_rows,
    }


def solve_theta(
    *,
    label: str,
    case: str,
    theta: dict[str, float],
    args: argparse.Namespace,
    targets: dict[str, float],
    weights: dict[str, float],
    moment_keys: list[str],
) -> dict[str, Any]:
    t0 = time.perf_counter()
    overrides = {
        **base_overrides(J=int(args.J), Nb=int(args.Nb), n_house=int(args.n_house), max_iter_eq=int(args.max_iter_eq)),
        **income_process_overrides(int(args.income_states)),
        **theta,
    }
    try:
        sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
        moments_all = extract_moments(sol, P)
        moments = {k: float(moments_all.get(k, np.nan)) for k in moment_keys}
        rank_loss = diagnostic_loss(moments_all, targets=targets, weights=weights)
        status = "ok"
        market_residual = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
        p_scalar = float(np.asarray(p_eq, dtype=float).reshape(-1)[0])
        iterations = float(getattr(sol, "iterations_completed", np.nan))
    except Exception as exc:  # noqa: BLE001 - audit should preserve failed perturbations.
        moments = {k: math.nan for k in moment_keys}
        rank_loss = math.inf
        status = f"failed: {type(exc).__name__}: {exc}"
        market_residual = math.inf
        p_scalar = math.nan
        iterations = math.nan
    return {
        "point_label": label,
        "case": case,
        "status": status,
        "theta": theta,
        "moments": moments,
        "rank_loss": float(rank_loss),
        "market_residual": float(market_residual),
        "p_eq": float(p_scalar),
        "iterations_completed": float(iterations),
        "elapsed_sec": float(time.perf_counter() - t0),
    }


def perturbed_thetas(
    theta: dict[str, float],
    parameter: str,
    rel_step: float,
    min_step: float,
) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    value = float(theta[parameter])
    lo, hi = parameter_bounds()[parameter]
    scale = parameter_scale(parameter, value)
    raw_step = max(abs(float(rel_step)) * scale, float(min_step))
    plus_value = min(value + raw_step, hi)
    minus_value = max(value - raw_step, lo)
    if not plus_value > value:
        plus_value = value
    if not minus_value < value:
        minus_value = value
    plus = dict(theta)
    minus = dict(theta)
    plus[parameter] = float(plus_value)
    minus[parameter] = float(minus_value)
    return plus, minus, {
        "base_value": value,
        "scale": float(scale),
        "raw_step": float(raw_step),
        "plus_step": float(plus_value - value),
        "minus_step": float(value - minus_value),
        "lower_bound": float(lo),
        "upper_bound": float(hi),
    }


def build_jacobian_rows(
    label: str,
    theta0: dict[str, float],
    base_record: dict[str, Any],
    plus_minus: dict[str, dict[str, dict[str, Any]]],
    targets: dict[str, float],
    target_moments: list[str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    raw_matrix = np.full((len(target_moments), len(PARAMETERS)), np.nan)
    scaled_matrix = np.full_like(raw_matrix, np.nan)
    base_moments = dict(base_record["moments"])
    for j, parameter in enumerate(PARAMETERS):
        minus = plus_minus[parameter]["minus"]
        plus = plus_minus[parameter]["plus"]
        base_value = float(theta0[parameter])
        plus_value = float(plus["theta"][parameter])
        minus_value = float(minus["theta"][parameter])
        denom = plus_value - minus_value
        pscale = parameter_scale(parameter, base_value)
        for i, moment in enumerate(target_moments):
            base_moment = float(base_moments.get(moment, np.nan))
            plus_moment = float(plus["moments"].get(moment, np.nan))
            minus_moment = float(minus["moments"].get(moment, np.nan))
            if np.isfinite(denom) and abs(denom) > 0 and np.isfinite(plus_moment) and np.isfinite(minus_moment):
                deriv = (plus_moment - minus_moment) / denom
            else:
                deriv = math.nan
            mscale = moment_scale(moment, base_moment, targets.get(moment, math.nan))
            scaled = deriv * pscale / mscale if np.isfinite(deriv) else math.nan
            raw_matrix[i, j] = deriv
            scaled_matrix[i, j] = scaled
            rows.append(
                {
                    "point_label": label,
                    "moment": moment,
                    "parameter": parameter,
                    "base_parameter": base_value,
                    "minus_parameter": minus_value,
                    "plus_parameter": plus_value,
                    "base_moment": base_moment,
                    "minus_moment": minus_moment,
                    "plus_moment": plus_moment,
                    "target": float(targets.get(moment, math.nan)),
                    "raw_derivative": deriv,
                    "parameter_scale": pscale,
                    "moment_scale": mscale,
                    "scaled_derivative": scaled,
                    "minus_status": minus["status"],
                    "plus_status": plus["status"],
                }
            )
    payload = {
        "scaled_matrix": scaled_matrix,
        "raw_matrix": raw_matrix,
        "json": {
            "point_label": label,
            "target_moments": target_moments,
            "parameters": PARAMETERS,
            "raw_matrix": jsonable(raw_matrix),
            "scaled_matrix": jsonable(scaled_matrix),
        },
    }
    return rows, payload


def rank_summary(
    label: str,
    scaled_matrix: np.ndarray,
    target_moments: list[str],
    parameters: list[str],
    rank_tol: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    mat = np.asarray(scaled_matrix, dtype=float)
    finite_mask = np.isfinite(mat)
    mat_clean = np.where(finite_mask, mat, 0.0)
    singular = np.linalg.svd(mat_clean, compute_uv=False)
    max_s = float(singular[0]) if singular.size else 0.0
    tol_abs = float(rank_tol) * max(max_s, 1.0)
    rank = int(np.sum(singular > tol_abs))
    condition = float(max_s / singular[-1]) if singular.size and singular[-1] > 0 else math.inf
    rows = [
        {
            "point_label": label,
            "singular_index": int(i + 1),
            "singular_value": float(s),
            "relative_to_largest": float(s / max_s) if max_s > 0 else math.nan,
            "rank_tol_absolute": tol_abs,
            "rank": rank,
            "condition_number": condition,
            "n_target_moments": len(target_moments),
            "n_parameters": len(parameters),
            "n_nonfinite_entries": int(np.size(mat) - np.sum(finite_mask)),
        }
        for i, s in enumerate(singular)
    ]
    return rows, {
        "singular_values": [float(s) for s in singular],
        "rank": rank,
        "condition_number": condition,
        "rank_tol_absolute": tol_abs,
    }


def column_collinearity_rows(label: str, scaled_matrix: np.ndarray, parameters: list[str]) -> list[dict[str, Any]]:
    mat = np.asarray(scaled_matrix, dtype=float)
    mat = np.where(np.isfinite(mat), mat, 0.0)
    rows: list[dict[str, Any]] = []
    for a in range(len(parameters)):
        va = mat[:, a]
        na = float(np.linalg.norm(va))
        for b in range(a + 1, len(parameters)):
            vb = mat[:, b]
            nb = float(np.linalg.norm(vb))
            corr = float(np.dot(va, vb) / (na * nb)) if na > 0 and nb > 0 else math.nan
            rows.append(
                {
                    "point_label": label,
                    "parameter_a": parameters[a],
                    "parameter_b": parameters[b],
                    "cosine_similarity": corr,
                    "abs_cosine_similarity": abs(corr) if np.isfinite(corr) else math.nan,
                    "norm_a": na,
                    "norm_b": nb,
                }
            )
    rows.sort(key=lambda r: float(r["abs_cosine_similarity"]) if np.isfinite(float(r["abs_cosine_similarity"])) else -1.0, reverse=True)
    return rows


def moment_top_sensitivity_rows(
    label: str,
    scaled_matrix: np.ndarray,
    target_moments: list[str],
    parameters: list[str],
    *,
    top_n: int = 5,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    mat = np.asarray(scaled_matrix, dtype=float)
    for i, moment in enumerate(target_moments):
        vals = [(parameters[j], float(mat[i, j])) for j in range(len(parameters)) if np.isfinite(mat[i, j])]
        vals.sort(key=lambda x: abs(x[1]), reverse=True)
        for rank, (parameter, value) in enumerate(vals[:top_n], start=1):
            rows.append(
                {
                    "point_label": label,
                    "moment": moment,
                    "rank": rank,
                    "parameter": parameter,
                    "scaled_derivative": value,
                    "abs_scaled_derivative": abs(value),
                }
            )
    return rows


def parameter_top_moment_rows(
    label: str,
    scaled_matrix: np.ndarray,
    target_moments: list[str],
    parameters: list[str],
    *,
    top_n: int = 5,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    mat = np.asarray(scaled_matrix, dtype=float)
    for j, parameter in enumerate(parameters):
        vals = [(target_moments[i], float(mat[i, j])) for i in range(len(target_moments)) if np.isfinite(mat[i, j])]
        vals.sort(key=lambda x: abs(x[1]), reverse=True)
        for rank, (moment, value) in enumerate(vals[:top_n], start=1):
            rows.append(
                {
                    "point_label": label,
                    "parameter": parameter,
                    "rank": rank,
                    "moment": moment,
                    "scaled_derivative": value,
                    "abs_scaled_derivative": abs(value),
                }
            )
    return rows


def baseline_moment_rows(
    label: str,
    base_record: dict[str, Any],
    targets: dict[str, float],
    target_moments: list[str],
) -> list[dict[str, Any]]:
    rows = []
    for moment in target_moments:
        value = float(base_record["moments"].get(moment, math.nan))
        target = float(targets.get(moment, math.nan))
        rows.append(
            {
                "point_label": label,
                "moment": moment,
                "model": value,
                "target": target,
                "gap": value - target if np.isfinite(value) and np.isfinite(target) else math.nan,
                "rank_loss": base_record["rank_loss"],
                "market_residual": base_record["market_residual"],
            }
        )
    return rows


def case_summary_row(
    label: str,
    record: dict[str, Any],
    *,
    parameter: str,
    side: str,
    step: float,
) -> dict[str, Any]:
    row = {
        "point_label": label,
        "case": record["case"],
        "parameter": parameter,
        "side": side,
        "step": float(step),
        "status": record["status"],
        "rank_loss": record["rank_loss"],
        "market_residual": record["market_residual"],
        "p_eq": record["p_eq"],
        "iterations_completed": record["iterations_completed"],
        "elapsed_sec": record["elapsed_sec"],
    }
    for name in ["tfr", "childless_rate", "own_rate", "own_family_gap", "old_age_own_rate", "housing_user_cost_share"]:
        row[name] = record["moments"].get(name, math.nan)
    return row


def load_frontier_best(frontier_root: Path, point_label: str) -> dict[str, Any]:
    summaries = sorted(frontier_root.glob(f"intergen_overnight_{point_label}_w*_20260617/collected_summary.json"))
    if not summaries:
        raise FileNotFoundError(f"No collected summaries found for {point_label!r} under {frontier_root}")
    best: dict[str, Any] | None = None
    best_path: Path | None = None
    for path in summaries:
        obj = json.loads(path.read_text())
        candidate = dict(obj.get("best") or {})
        loss = float(candidate.get("rank_loss", math.inf))
        if best is None or loss < float(best.get("rank_loss", math.inf)):
            best = candidate
            best_path = path
    if best is None or "theta" not in best:
        raise ValueError(f"No best theta found for {point_label!r}")
    return {
        "point_label": point_label,
        "source_target_set": source_target_set_from_label(point_label),
        "source_rank_loss": float(best.get("rank_loss", math.nan)),
        "source_path": str(best_path),
        "source_case": best.get("case"),
        "source_label": best.get("label"),
        "theta": {k: float(v) for k, v in dict(best["theta"]).items() if k in PARAMETERS},
        "source_moments": dict(best.get("moments", {})),
    }


def source_target_set_from_label(point_label: str) -> str:
    return "candidate_no_timing_" + point_label


def all_point_labels() -> list[str]:
    return [
        "core_feasibility_v1",
        "cost_test_v1",
        "oldage_test_v1",
        "roomcost_test_v1",
    ]


def parameter_bounds() -> dict[str, tuple[float, float]]:
    bounds: dict[str, tuple[float, float]] = {}
    for name, lo, hi in GLOBAL_DE_BOUNDS:
        if name == "beta_annual":
            bounds["beta"] = (float(lo) ** 4.0, float(hi) ** 4.0)
        else:
            bounds[name] = (float(lo), float(hi))
    return bounds


def parameter_scale(parameter: str, value: float) -> float:
    if parameter == "beta":
        return max(abs(float(value)), 0.05)
    return max(abs(float(value)), 1.0)


def moment_scale(moment: str, value: float, target: float) -> float:
    if moment in {"prime_childless_renter_median_rooms", "prime_childless_owner_median_rooms"}:
        return max(abs(float(value)), abs(float(target)), 1.0)
    return max(abs(float(value)), abs(float(target)), 1.0)


def unique_list(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value not in seen:
            out.append(value)
            seen.add(value)
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(jsonable(row.get(key))) if isinstance(row.get(key), (dict, list)) else row.get(key) for key in fields})


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True))


def write_readme(
    outdir: Path,
    summary: dict[str, Any],
    moment_top_rows: list[dict[str, Any]],
    param_top_rows: list[dict[str, Any]],
    collinearity_rows: list[dict[str, Any]],
) -> None:
    lines = [
        "# Intergen Sensitivity Jacobian Audit",
        "",
        "This is a finite-difference identification diagnostic for the current one-market intergenerational model.",
        "It is not a calibration run.",
        "",
        "## Configuration",
        "",
        f"- target set: `{summary['config']['target_set']}`",
        f"- points: `{', '.join(summary['config']['point_labels'])}`",
        f"- relative perturbation: `{summary['config']['rel_step']}`",
        f"- numerical stack: `J={summary['config']['J']}`, `Nb={summary['config']['Nb']}`, `n_house={summary['config']['n_house']}`, `max_iter_eq={summary['config']['max_iter_eq']}`",
        "",
        "## Point Summaries",
        "",
    ]
    for point in summary["points"]:
        lines.extend(
            [
                f"### {point['point_label']}",
                "",
                f"- baseline rank loss: `{point['baseline_rank_loss']:.6g}`",
                f"- market residual: `{point['baseline_market_residual']:.3g}`",
                f"- scaled rank: `{point['scaled_rank']}`",
                f"- condition number: `{point['condition_number']:.6g}`",
                f"- singular values: `{', '.join(f'{s:.3g}' for s in point['singular_values'])}`",
                "",
            ]
        )

    lines.extend(
        [
            "## Files",
            "",
            "- `baseline_moments.csv`: baseline model moments and targets.",
            "- `finite_difference_cases.csv`: all perturbation solves.",
            "- `jacobian_long.csv`: raw and scaled derivatives.",
            "- `scaled_singular_values.csv`: local rank diagnostics.",
            "- `scaled_column_correlations.csv`: near-collinear parameter columns.",
            "- `moment_top_sensitivities.csv`: largest parameter effects by moment.",
            "- `parameter_top_moments.csv`: largest moment effects by parameter.",
            "",
            "## Interpretation Note",
            "",
            "The scaled Jacobian uses `raw_derivative * parameter_scale / moment_scale`, ",
            "where `parameter_scale = max(abs(theta), 1)` except for beta and ",
            "`moment_scale = max(abs(model moment), abs(target), 1)`. This is a ",
            "local conditioning diagnostic, not an economic elasticity. Discrete ",
            "median-room moments can be locally flat under small perturbations even ",
            "when they matter globally.",
            "",
            "## Top Local Collinearities",
            "",
        ]
    )
    for row in collinearity_rows[:12]:
        lines.append(
            f"- `{row['point_label']}`: `{row['parameter_a']}` vs `{row['parameter_b']}`, "
            f"cosine `{row['cosine_similarity']:.3f}`"
        )
    lines.append("")
    lines.append("## Top Moment Sensitivities")
    lines.append("")
    for row in moment_top_rows:
        if int(row["rank"]) <= 2:
            lines.append(
                f"- `{row['point_label']}` `{row['moment']}`: `{row['parameter']}` "
                f"scaled derivative `{row['scaled_derivative']:.3g}`"
            )
    lines.append("")
    lines.append("## Top Parameter Loads")
    lines.append("")
    for row in param_top_rows:
        if int(row["rank"]) == 1:
            lines.append(
                f"- `{row['point_label']}` `{row['parameter']}` mainly moves `{row['moment']}` "
                f"with scaled derivative `{row['scaled_derivative']:.3g}`"
            )
    lines.append("")
    (outdir / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
