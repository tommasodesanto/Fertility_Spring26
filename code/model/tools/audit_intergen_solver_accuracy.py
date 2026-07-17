#!/usr/bin/env python3
"""Run read-only numerical accuracy diagnostics for the intergen model.

The script intentionally does not change solver logic. It solves a saved
parameter vector under several numerical configurations, then writes tables for
grid convergence, shape/KFE checks, Euler residuals on smooth renter states,
equilibrium-quality sensitivity, and implementation-path equivalence.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import pickle
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
TOOLS_ROOT = Path(__file__).resolve().parent
for path in (MODEL_ROOT, TOOLS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from build_intergen_mechanics_packet import load_source_record  # noqa: E402
from intergen_housing_fertility.calibration import (  # noqa: E402
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_housing_fertility.local_panel import income_process_overrides  # noqa: E402
from intergen_housing_fertility.solver import (  # noqa: E402
    apply_child_aging,
    income_at_state,
    income_transition_values,
    precompute_shared,
    run_model_cp_dt,
    solve_markov_income_at_prices,
)
from intergen_housing_fertility.utils import flat_nc, make_grid  # noqa: E402


DEFAULT_SOURCE = (
    ROOT
    / "output/model/intergen_fixedstats_overnight_review_20260626/"
    / "candidates/de_w3_task10_de_g044_i022_best.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/intergen_solver_accuracy_20260626"
DEFAULT_TARGET_SET = "candidate_replacement_roomgap_14moment_v1"
NEG_INF_CUTOFF = -1.0e9


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    source = load_source_record(args.source.resolve(), row=int(args.source_row))
    targets, weights = get_target_set(args.target_set)
    nb_list = parse_int_list(args.nb_list)
    max_iter_list = parse_int_list(args.max_iter_list)
    fixed_price = resolve_fixed_price(args, source)

    metadata = {
        "created_at": dt.datetime.now().isoformat(timespec="seconds"),
        "status": "read_only_numerical_diagnostic_not_calibration",
        "source": source["source_meta"],
        "target_set": args.target_set,
        "targets": targets,
        "weights": weights,
        "nb_list": nb_list,
        "max_iter_list": max_iter_list,
        "fixed_price": fixed_price,
        "notes": [
            "fixed_price grid convergence isolates household-solver/grid error",
            "ge grid convergence adds price-feedback error",
            "Euler residuals are only computed for smooth interior renter states",
        ],
    }
    write_json(outdir / "metadata.json", metadata)

    solved_cases: list[dict[str, Any]] = []
    if not args.skip_grid_convergence:
        for mode in ("fixed_price", "ge"):
            for nb in nb_list:
                label = f"{mode}_Nb{nb}"
                extra = {}
                if mode == "fixed_price":
                    extra = {"solve_mode": "pe", "p_fixed": np.array([fixed_price])}
                case = solve_case(
                    label=label,
                    theta=source["theta"],
                    target_set=args.target_set,
                    targets=targets,
                    weights=weights,
                    J=int(args.J),
                    Nb=int(nb),
                    income_states=int(args.income_states),
                    n_house=int(args.n_house),
                    H_own=parse_float_list(args.H_own),
                    hR_max=float(args.hR_max),
                    max_iter_eq=int(args.max_iter_eq),
                    interp_method=str(args.interp_method),
                    extra_overrides=extra,
                    save_cache=bool(args.save_case_caches),
                    cache_dir=outdir / "solution_caches",
                )
                solved_cases.append(case)
                write_case_checkpoint(outdir, solved_cases, targets, weights)

    if not args.skip_equilibrium_quality:
        for refine in (True, False):
            for mi in max_iter_list:
                label = f"eq_quality_refine{int(refine)}_maxiter{mi}"
                case = solve_case(
                    label=label,
                    theta=source["theta"],
                    target_set=args.target_set,
                    targets=targets,
                    weights=weights,
                    J=int(args.J),
                    Nb=int(args.eq_nb),
                    income_states=int(args.income_states),
                    n_house=int(args.n_house),
                    H_own=parse_float_list(args.H_own),
                    hR_max=float(args.hR_max),
                    max_iter_eq=int(mi),
                    interp_method=str(args.interp_method),
                    extra_overrides={"scalar_market_refine": bool(refine)},
                    save_cache=False,
                    cache_dir=outdir / "solution_caches",
                )
                solved_cases.append(case)
                write_case_checkpoint(outdir, solved_cases, targets, weights)

    if not args.skip_implementation_equivalence:
        equivalence_rows = run_implementation_equivalence(args, source, targets, weights, fixed_price, outdir)
        write_csv(outdir / "implementation_equivalence.csv", equivalence_rows)

    if solved_cases:
        write_case_checkpoint(outdir, solved_cases, targets, weights)
        shape_rows = [shape_kfe_checks(case) for case in solved_cases]
        write_csv(outdir / "shape_kfe_checks.csv", shape_rows)
        write_json(outdir / "shape_kfe_checks.json", shape_rows)

        baseline = pick_case(solved_cases, "ge_Nb" + str(nb_list[0])) or solved_cases[0]
        euler_rows, euler_summary = euler_residuals(baseline)
        write_csv(outdir / "euler_residuals.csv", euler_rows)
        write_json(outdir / "euler_residual_summary.json", euler_summary)

    write_readme(outdir, args, source, solved_cases)
    print(outdir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--source-row", type=int, default=1)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--target-set", default=DEFAULT_TARGET_SET)
    parser.add_argument("--J", type=int, default=17)
    parser.add_argument("--income-states", type=int, default=5)
    parser.add_argument("--n-house", type=int, default=5)
    parser.add_argument("--H-own", default="2,4,6,8,10")
    parser.add_argument("--hR-max", type=float, default=6.0)
    parser.add_argument("--max-iter-eq", type=int, default=10)
    parser.add_argument("--eq-nb", type=int, default=60)
    parser.add_argument("--nb-list", default="60,120,240")
    parser.add_argument("--max-iter-list", default="3,5,10,25")
    parser.add_argument("--interp-method", default="linear", choices=["linear", "monotone_cubic"])
    parser.add_argument("--fixed-price", type=float, default=None)
    parser.add_argument("--save-case-caches", action="store_true")
    parser.add_argument("--skip-grid-convergence", action="store_true")
    parser.add_argument("--skip-equilibrium-quality", action="store_true")
    parser.add_argument("--skip-implementation-equivalence", action="store_true")
    return parser.parse_args()


def solve_case(
    *,
    label: str,
    theta: dict[str, float],
    target_set: str,
    targets: dict[str, float],
    weights: dict[str, float],
    J: int,
    Nb: int,
    income_states: int,
    n_house: int,
    H_own: list[float],
    hR_max: float,
    max_iter_eq: int,
    interp_method: str,
    extra_overrides: dict[str, Any],
    save_cache: bool,
    cache_dir: Path,
) -> dict[str, Any]:
    overrides = {
        **base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq),
        **income_process_overrides(income_states),
        **theta,
        "H_own": np.asarray(H_own, dtype=float),
        "n_house": len(H_own),
        "hR_max": float(hR_max),
        "interp_method": str(interp_method),
        **extra_overrides,
    }
    t0 = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
    elapsed = time.perf_counter() - t0
    moments = extract_moments(sol, P)
    loss = diagnostic_loss(moments, targets=targets, weights=weights)
    row = {
        "label": label,
        "target_set": target_set,
        "J": int(J),
        "Nb": int(Nb),
        "income_states": int(income_states),
        "n_house": int(len(H_own)),
        "hR_max": float(hR_max),
        "max_iter_eq": int(max_iter_eq),
        "interp_method": str(interp_method),
        "solve_mode": str(extra_overrides.get("solve_mode", "ge")),
        "scalar_market_refine": bool(extra_overrides.get("scalar_market_refine", getattr(P, "scalar_market_refine", True))),
        "p_eq": float(np.asarray(p_eq, dtype=float).reshape(-1)[0]),
        "rank_loss": float(loss),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "elapsed_sec": float(elapsed),
        "moments": moments,
        "targets": targets,
        "weights": weights,
        "timings": jsonable(getattr(sol, "timings", {})),
        "sol": sol,
        "P": P,
        "p_eq_vector": np.asarray(p_eq, dtype=float).reshape(-1),
    }
    if save_cache:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / f"{safe_label(label)}_solution_cache.pkl"
        with cache_path.open("wb") as fh:
            pickle.dump(
                {
                    "status": "diagnostic_cache_load_only_if_trusted",
                    "label": label,
                    "target_set": target_set,
                    "baseline": {
                        "sol": sol,
                        "P": P,
                        "p_eq": row["p_eq_vector"],
                        "moments": moments,
                        "rank_loss": float(loss),
                        "market_residual": row["market_residual"],
                        "elapsed_sec": float(elapsed),
                    },
                },
                fh,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        row["solution_cache"] = str(cache_path)
    print(f"{label}: loss={loss:.6g} residual={row['market_residual']:.2e} p={row['p_eq']:.6g}", flush=True)
    return row


def write_case_checkpoint(
    outdir: Path,
    cases: list[dict[str, Any]],
    targets: dict[str, float],
    weights: dict[str, float],
) -> None:
    rows = []
    target_rows = []
    for case in cases:
        moments = case["moments"]
        base = {k: v for k, v in case.items() if k not in {"sol", "P", "moments", "targets", "weights", "timings", "p_eq_vector"}}
        for name, value in moments.items():
            base[f"moment__{name}"] = value
        rows.append(base)
        for name in sorted(targets):
            model = float(moments.get(name, np.nan))
            target = float(targets[name])
            weight = float(weights.get(name, 1.0))
            target_rows.append(
                {
                    "label": case["label"],
                    "moment": name,
                    "target": target,
                    "model": model,
                    "difference": model - target if np.isfinite(model) else np.nan,
                    "weight": weight,
                    "loss_contribution": weight * (model - target) ** 2 if np.isfinite(model) else np.nan,
                }
            )
    write_csv(outdir / "case_moments_wide.csv", rows)
    write_csv(outdir / "target_fit_long.csv", target_rows)
    write_json(
        outdir / "case_summaries.json",
        [
            {k: jsonable(v) for k, v in case.items() if k not in {"sol", "P"}}
            for case in cases
        ],
    )
    write_grid_drift_tables(outdir, cases, targets)


def write_grid_drift_tables(outdir: Path, cases: list[dict[str, Any]], targets: dict[str, float]) -> None:
    rows = []
    for mode in ("fixed_price", "ge"):
        mode_cases = [case for case in cases if str(case["label"]).startswith(mode + "_Nb")]
        if not mode_cases:
            continue
        mode_cases = sorted(mode_cases, key=lambda c: int(c["Nb"]))
        ref = mode_cases[0]
        high = mode_cases[-1]
        for name in sorted(targets):
            ref_val = float(ref["moments"].get(name, np.nan))
            high_val = float(high["moments"].get(name, np.nan))
            for case in mode_cases:
                val = float(case["moments"].get(name, np.nan))
                rows.append(
                    {
                        "mode": mode,
                        "moment": name,
                        "target": float(targets[name]),
                        "label": case["label"],
                        "Nb": int(case["Nb"]),
                        "model": val,
                        "drift_from_lowest_Nb": val - ref_val if np.isfinite(val) and np.isfinite(ref_val) else np.nan,
                        "drift_from_highest_Nb": val - high_val if np.isfinite(val) and np.isfinite(high_val) else np.nan,
                    }
                )
    write_csv(outdir / "grid_convergence_target_moment_drift.csv", rows)


def shape_kfe_checks(case: dict[str, Any]) -> dict[str, Any]:
    sol = case["sol"]
    P = case["P"]
    g = np.asarray(sol.g, dtype=float)
    V = np.asarray(sol.V, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    total_mass = float(np.sum(g))
    mass_b = np.sum(g, axis=tuple(range(1, g.ndim)))
    edge_low = float(mass_b[0] / max(total_mass, 1e-14))
    edge_high = float(mass_b[-1] / max(total_mass, 1e-14))
    mass_age = np.sum(g, axis=(0, 1, 2, 4, 5, 6))
    expected_age_mass = float(np.mean(mass_age)) if mass_age.size else np.nan
    finite = np.isfinite(V) & (V > NEG_INF_CUTOFF)
    diffs = np.diff(np.where(finite, V, np.nan), axis=0)
    mono_viol = diffs < -1.0e-7
    reachable_slices = np.sum(g, axis=0) > 1e-12
    reachable_viol = 0
    if reachable_slices.shape == V.shape[1:]:
        reachable_viol = int(np.sum(mono_viol & reachable_slices[None, ...]))
    approx_contam = reachable_continuation_penalty_crossings(sol, P, b_grid)
    z_grid, z_weights, Pi_z = income_transition_values(P)
    child_row_error = np.nan
    if hasattr(P, "Pi_child"):
        child = np.asarray(P.Pi_child, dtype=float)
        child_row_error = float(np.nanmax(np.abs(np.sum(child, axis=1) - 1.0)))
    return {
        "label": case["label"],
        "Nb": int(case["Nb"]),
        "solve_mode": case["solve_mode"],
        "interp_method": case["interp_method"],
        "rank_loss": float(case["rank_loss"]),
        "market_residual": float(case["market_residual"]),
        "total_mass": total_mass,
        "edge_mass_low_share": edge_low,
        "edge_mass_high_share": edge_high,
        "max_abs_age_mass_deviation": float(np.max(np.abs(mass_age - expected_age_mass))) if mass_age.size else np.nan,
        "income_transition_max_row_error": float(np.max(np.abs(np.sum(Pi_z, axis=1) - 1.0))),
        "income_weight_sum_error": float(abs(np.sum(z_weights) - 1.0)),
        "child_transition_max_row_error": child_row_error,
        "V_monotonicity_violations_all": int(np.nansum(mono_viol)),
        "V_monotonicity_violations_reachable_slices": int(reachable_viol),
        "reachable_negative_inf_value_cells": int(np.sum((g > 1e-14) & (V <= NEG_INF_CUTOFF))),
        **approx_contam,
    }


def reachable_continuation_penalty_crossings(sol: Any, P: Any, b_grid: np.ndarray) -> dict[str, Any]:
    # Conservative diagnostic: for positive-mass cells, check whether the stored
    # same-tenure b' interpolation would touch penalty values in next-period V.
    g = np.asarray(sol.g, dtype=float)
    bp = np.asarray(sol.bp_pol, dtype=float)
    V = np.asarray(sol.V, dtype=float)
    if g.shape != bp.shape or V.shape != bp.shape:
        return {
            "reachable_same_tenure_bp_penalty_crossings": np.nan,
            "reachable_same_tenure_bp_penalty_crossing_mass": np.nan,
        }
    total = 0
    mass = 0.0
    bmin = float(b_grid[0])
    bmax = float(b_grid[-1])
    for idx in np.argwhere(g > 1e-12):
        bb, ten, i, j, zz, nn, cs = [int(x) for x in idx]
        if j >= int(P.J) - 1:
            continue
        q = float(np.clip(bp[bb, ten, i, j, zz, nn, cs], bmin, bmax))
        hi = int(np.searchsorted(b_grid, q, side="right"))
        lo = max(0, min(hi - 1, b_grid.size - 2))
        hi = lo + 1
        endpoints_bad = False
        for znext in range(V.shape[4]):
            vals = V[[lo, hi], ten, i, j + 1, znext, nn, cs]
            if np.any(vals <= NEG_INF_CUTOFF):
                endpoints_bad = True
                break
        if endpoints_bad:
            total += 1
            mass += float(g[tuple(idx)])
    return {
        "reachable_same_tenure_bp_penalty_crossings": int(total),
        "reachable_same_tenure_bp_penalty_crossing_mass": float(mass),
    }


def euler_residuals(case: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    sol = case["sol"]
    P = case["P"]
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    g = np.asarray(sol.g, dtype=float)
    V = np.asarray(sol.V, dtype=float)
    cpol = np.asarray(sol.c_pol, dtype=float)
    hpol = np.asarray(sol.hR_pol, dtype=float)
    bppol = np.asarray(sol.bp_pol, dtype=float)
    tp = getattr(sol, "tenure_probs", None)
    tp_arr = np.asarray(tp, dtype=float) if tp is not None else None
    z_grid, z_weights, Pi_z = income_transition_values(P)
    rows = []
    classified = {"smooth_interior": 0, "constraint_or_cap": 0, "logit_mixed": 0, "nonfinite": 0}
    alpha = float(P.alpha_cons)
    sigma = float(P.sigma)
    beta = float(P.beta)
    oms = 1.0 - sigma
    interior_count = 0
    for idx in np.argwhere(g[:, 0, :, :-1, :, :, :] > 1e-10):
        bb, i, j, zz, nn, cs = [int(x) for x in idx]
        if j >= int(P.J) - 1:
            continue
        if tp_arr is not None and tp_arr[bb, 0, i, j, zz, nn, cs, 0] < 0.90:
            classified["logit_mixed"] += 1
            continue
        bp = float(bppol[bb, 0, i, j, zz, nn, cs])
        c = float(cpol[bb, 0, i, j, zz, nn, cs])
        h = float(hpol[bb, 0, i, j, zz, nn, cs])
        if not all(math.isfinite(x) for x in (bp, c, h)) or c <= 0 or h <= 0:
            classified["nonfinite"] += 1
            continue
        if bp <= b_grid[0] + 1e-6 or bp >= b_grid[-1] - 1e-6 or bp <= 1e-5:
            classified["constraint_or_cap"] += 1
            continue
        if abs(h - float(P.hR_max)) < 1e-5:
            classified["constraint_or_cap"] += 1
            continue
        Vnr = np.zeros((b_grid.size, 1 + int(P.n_house), int(P.I), int(P.n_parity), int(P.n_child_states)))
        for znext in range(len(z_grid)):
            Vnr += Pi_z[zz, znext] * V[:, :, :, j + 1, znext, :, :]
        Vc = apply_child_aging(Vnr, P, b_grid.size, 1 + int(P.n_house), int(P.I), int(P.n_parity), int(P.n_child_states))
        cidx = int(nn + int(P.n_parity) * cs)
        Vbar = flat_nc(Vc[:, 0, i, :, :], b_grid.size, int(P.n_parity) * int(P.n_child_states))[:, cidx]
        dV = finite_diff_interp_derivative(b_grid, Vbar, bp)
        if not math.isfinite(dV):
            classified["nonfinite"] += 1
            continue
        muc = alpha * (c ** (alpha * oms - 1.0)) * (h ** ((1.0 - alpha) * oms))
        rhs = beta * dV
        if not math.isfinite(muc) or not math.isfinite(rhs) or muc <= 0 or rhs <= 0:
            classified["nonfinite"] += 1
            continue
        residual = math.log10(abs(rhs / muc - 1.0) + 1e-16)
        rows.append(
            {
                "label": case["label"],
                "b_index": bb,
                "b": float(b_grid[bb]),
                "age_index": j,
                "age": float(getattr(P, "age_start", 22.0) + j * getattr(P, "period_years", 4.0)),
                "income_index": zz,
                "z": float(z_grid[zz]),
                "parity": nn,
                "child_state": cs,
                "mass": float(g[bb, 0, i, j, zz, nn, cs]),
                "bp": bp,
                "c": c,
                "h": h,
                "marginal_utility_c": muc,
                "beta_expected_Vb": rhs,
                "log10_abs_relative_residual": residual,
            }
        )
        classified["smooth_interior"] += 1
        interior_count += 1
        if interior_count >= 5000:
            break
    residuals = np.asarray([r["log10_abs_relative_residual"] for r in rows], dtype=float)
    weights = np.asarray([r["mass"] for r in rows], dtype=float)
    summary = {
        "case_label": case["label"],
        "classification_counts": classified,
        "n_residual_rows": int(len(rows)),
        "weighted_median_log10_abs_relative_residual": weighted_quantile(residuals, weights, 0.5) if len(rows) else np.nan,
        "weighted_p95_log10_abs_relative_residual": weighted_quantile(residuals, weights, 0.95) if len(rows) else np.nan,
        "max_log10_abs_relative_residual": float(np.max(residuals)) if len(rows) else np.nan,
        "note": "Euler residuals only use renter, high-renter-probability, interior, non-cap states.",
    }
    return rows, summary


def finite_diff_interp_derivative(x: np.ndarray, y: np.ndarray, q: float) -> float:
    if q <= x[0] or q >= x[-1]:
        return math.nan
    hi = int(np.searchsorted(x, q, side="right"))
    lo = max(0, min(hi - 1, x.size - 2))
    y0 = float(y[lo])
    y1 = float(y[lo + 1])
    if y0 <= NEG_INF_CUTOFF or y1 <= NEG_INF_CUTOFF:
        return math.nan
    return float((y1 - y0) / (x[lo + 1] - x[lo]))


def run_implementation_equivalence(
    args: argparse.Namespace,
    source: dict[str, Any],
    targets: dict[str, float],
    weights: dict[str, float],
    fixed_price: float,
    outdir: Path,
) -> list[dict[str, Any]]:
    base_kwargs = {
        "theta": source["theta"],
        "target_set": args.target_set,
        "targets": targets,
        "weights": weights,
        "J": 8,
        "Nb": 24,
        "income_states": 3,
        "n_house": 3,
        "H_own": [2.0, 6.0, 10.0],
        "hR_max": float(args.hR_max),
        "max_iter_eq": 3,
        "interp_method": "linear",
        "save_cache": False,
        "cache_dir": outdir / "solution_caches",
    }
    cases = [
        solve_case(label="small_numba_kernel", extra_overrides={"use_full_kernel": True}, **base_kwargs),
        solve_case(label="small_python_kernel", extra_overrides={"use_full_kernel": False}, **base_kwargs),
        solve_case(
            label="small_fixed_linear",
            extra_overrides={"solve_mode": "pe", "p_fixed": np.array([fixed_price]), "use_full_kernel": True},
            **base_kwargs,
        ),
        solve_case(
            label="small_fixed_monotone_cubic",
            extra_overrides={"solve_mode": "pe", "p_fixed": np.array([fixed_price]), "use_full_kernel": False},
            **{**base_kwargs, "interp_method": "monotone_cubic"},
        ),
    ]
    rows = []
    pairs = [
        (cases[1], cases[0], "python_kernel_minus_numba_kernel"),
        (cases[3], cases[2], "monotone_cubic_minus_linear_fixed_price"),
    ]
    for case, ref, comparison in pairs:
        for moment in sorted(set(ref["moments"]) | set(case["moments"])):
            a = float(ref["moments"].get(moment, np.nan))
            b = float(case["moments"].get(moment, np.nan))
            rows.append(
                {
                    "comparison": comparison,
                    "case": case["label"],
                    "reference": ref["label"],
                    "moment": moment,
                    "case_value": b,
                    "reference_value": a,
                    "difference": b - a if np.isfinite(a) and np.isfinite(b) else np.nan,
                }
            )
        rows.extend(array_comparison_rows(case["sol"], ref["sol"], case["label"], ref["label"], comparison))
    # Fast-stats vs full at one fixed price and one small-grid parameterization.
    full_case = cases[2]
    P = full_case["P"]
    p_eq = full_case["p_eq_vector"]
    b_grid = np.asarray(full_case["sol"].b_grid, dtype=float).reshape(-1)
    SD = precompute_shared(P, b_grid)
    fast = solve_markov_income_at_prices(p_eq, P, b_grid, verbose=False, fast_stats=True, SD=SD)
    full = solve_markov_income_at_prices(p_eq, P, b_grid, verbose=False, fast_stats=False, SD=SD)
    fast_m = extract_moments(fast, P)
    full_m = extract_moments(full, P)
    for moment in sorted(set(fast_m) | set(full_m)):
        a = float(full_m.get(moment, np.nan))
        b = float(fast_m.get(moment, np.nan))
        rows.append(
            {
                "comparison": "fast_stats_minus_full_stats",
                "case": "small_fixed_linear_fast_stats",
                "reference": "small_fixed_linear_full_stats",
                "moment": moment,
                "case_value": b,
                "reference_value": a,
                "difference": b - a if np.isfinite(a) and np.isfinite(b) else np.nan,
            }
        )
    rows.extend(array_comparison_rows(fast, full, "small_fixed_linear_fast_stats", "small_fixed_linear_full_stats", "fast_stats_minus_full_stats"))
    return rows


def array_comparison_rows(case_sol: Any, ref_sol: Any, case_label: str, ref_label: str, comparison: str) -> list[dict[str, Any]]:
    out = []
    for name in ("V", "c_pol", "hR_pol", "bp_pol", "tenure_choice", "tenure_probs", "loc_probs", "fert_probs", "g"):
        if not hasattr(case_sol, name) or not hasattr(ref_sol, name):
            continue
        a = getattr(case_sol, name)
        b = getattr(ref_sol, name)
        if a is None or b is None:
            continue
        aa = np.asarray(a)
        bb = np.asarray(b)
        row = {
            "comparison": comparison,
            "case": case_label,
            "reference": ref_label,
            "metric_type": "array_max_abs_difference",
            "object": name,
            "case_shape": str(tuple(aa.shape)),
            "reference_shape": str(tuple(bb.shape)),
        }
        if aa.shape == bb.shape and aa.size:
            diff = np.asarray(aa, dtype=float) - np.asarray(bb, dtype=float)
            finite = np.isfinite(diff)
            row.update(
                {
                    "max_abs_difference": float(np.nanmax(np.abs(diff))) if np.any(finite) else np.nan,
                    "mean_abs_difference": float(np.nanmean(np.abs(diff))) if np.any(finite) else np.nan,
                    "n_different_gt_1e-8": int(np.sum(np.abs(diff[finite]) > 1e-8)) if np.any(finite) else 0,
                }
            )
        else:
            row.update({"max_abs_difference": np.nan, "mean_abs_difference": np.nan, "n_different_gt_1e-8": np.nan})
        out.append(row)
    return out


def resolve_fixed_price(args: argparse.Namespace, source: dict[str, Any]) -> float:
    if args.fixed_price is not None:
        return float(args.fixed_price)
    raw = source.get("raw", {})
    if isinstance(raw, dict) and "p_eq" in raw:
        arr = np.asarray(raw["p_eq"], dtype=float).reshape(-1)
        if arr.size and np.isfinite(arr[0]):
            return float(arr[0])
    # Fall back to the current packet if present.
    summary = (
        ROOT
        / "output/model/intergen_fixedstats_overnight_review_20260626/"
        / "best_de_g044_i022_packet/solution_summary.json"
    )
    if summary.exists():
        payload = json.loads(summary.read_text())
        arr = np.asarray(payload.get("p_eq", []), dtype=float).reshape(-1)
        if arr.size and np.isfinite(arr[0]):
            return float(arr[0])
    raise ValueError("could not infer fixed price; pass --fixed-price")


def pick_case(cases: list[dict[str, Any]], label: str) -> dict[str, Any] | None:
    for case in cases:
        if case["label"] == label:
            return case
    return None


def parse_int_list(text: str) -> list[int]:
    return [int(x.strip()) for x in str(text).split(",") if x.strip()]


def parse_float_list(text: str) -> list[float]:
    return [float(x.strip()) for x in str(text).split(",") if x.strip()]


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=float).reshape(-1)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(ok):
        return math.nan
    values = values[ok]
    weights = weights[ok]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights) / np.sum(weights)
    return float(values[min(np.searchsorted(cdf, q), values.size - 1)])


def safe_label(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in str(label))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(jsonable(row.get(key))) if isinstance(row.get(key), (dict, list)) else row.get(key) for key in keys})


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True))


def write_readme(outdir: Path, args: argparse.Namespace, source: dict[str, Any], cases: list[dict[str, Any]]) -> None:
    lines = [
        "# Intergen Solver Accuracy Diagnostics",
        "",
        "Read-only numerical diagnostics. No solver logic or economic primitives are changed.",
        "",
        f"- source: `{source['source_meta'].get('path')}`",
        f"- target set: `{args.target_set}`",
        f"- output: `{outdir}`",
        "",
        "## Files",
        "",
        "- `case_moments_wide.csv`: one row per solve, with every extracted model moment.",
        "- `target_fit_long.csv`: every calibration target next to the model moment, weight, and contribution.",
        "- `grid_convergence_target_moment_drift.csv`: target-moment drift across `Nb` for fixed-price and GE runs.",
        "- `shape_kfe_checks.csv`: monotonicity, edge-mass, transition-row, and approximate penalty-contamination checks.",
        "- `euler_residuals.csv` / `euler_residual_summary.json`: interior renter Euler residuals, with kink/constraint/logit states excluded.",
        "- `implementation_equivalence.csv`: Python/Numba, fast/full stats, and linear/PCHIP fixed-price comparisons.",
        "",
        "## Completed Cases",
        "",
    ]
    for case in cases:
        lines.append(
            f"- `{case['label']}`: loss `{float(case['rank_loss']):.6g}`, "
            f"residual `{float(case['market_residual']):.2e}`, p `{float(case['p_eq']):.6g}`"
        )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
