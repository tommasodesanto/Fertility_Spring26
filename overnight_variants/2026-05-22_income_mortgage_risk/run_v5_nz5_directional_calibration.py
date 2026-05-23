#!/usr/bin/env python3
"""Small directional calibration probes for V5 outside closure at Nz=5."""

from __future__ import annotations

import argparse
import csv
import json
import time
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from typing import Callable

import numpy as np

from dt_cp_model.direct_calibration import compute_smm_loss
from dt_cp_model.objective import extract_moments
from dt_cp_model.solver import make_grid, solve_equilibrium_hank_z

from run_income_mortgage_risk_v2_scenarios import BENCHMARK_P, MOMENT_ORDER, namespace_float_dict
from run_income_mortgage_risk_v5_hank_z_outside_closure import prepare_parameters


OUTPUT_STEM = "v5_nz5_directional_calibration"


def configure_closure(P: SimpleNamespace, args: argparse.Namespace) -> SimpleNamespace:
    P.population_closure = "outside_option_benchmark_normalized"
    P.kappa_entry = float(args.kappa_entry)
    P.local_birth_entry_weight = 1.0
    P.target_city_entry_prob = float(args.target_city_entry_prob)
    P.calibrate_outside_value_to_entry_prob = True
    P.outside_benchmark_target_total_population = float(getattr(P, "N_target", 1.0))
    P.allow_uncalibrated_outside_value = False
    P.collect_ge_trace = True
    return P


def load_baseline(root: Path) -> dict[str, object]:
    path = root / "income_mortgage_risk_v5_hank_z_outside_closure_nz5.log"
    if not path.exists():
        raise FileNotFoundError(f"missing Nz=5 baseline log: {path}")
    data = json.loads(path.read_text())
    if not bool(data.get("ok", False)):
        raise RuntimeError(f"Nz=5 baseline log is not ok: {path}")
    return data


def record_param_snapshot(P: SimpleNamespace) -> dict[str, object]:
    return {
        "beta": float(P.beta),
        "psi_child": float(P.psi_child),
        "c_bar_n": float(P.c_bar_n),
        "h_bar_jump": float(P.h_bar_jump),
        "h_bar_n": float(P.h_bar_n),
        "kappa_fert": float(P.kappa_fert),
        "phi": float(np.asarray(P.phi).reshape(-1)[0]),
        "E_C": float(P.E_loc[1]),
        "r_bar_C": float(P.r_bar[1]),
        "theta0": float(P.theta0),
        "theta_n": float(P.theta_n),
        "alpha_cons": float(P.alpha_cons),
    }


def set_phi(P: SimpleNamespace, value: float) -> None:
    P.phi = float(value) * np.ones(P.n_parity)


def case_mutators() -> dict[str, Callable[[SimpleNamespace], None]]:
    def fertility_utility_max(P: SimpleNamespace) -> None:
        P.psi_child = 0.12

    def child_cost_low(P: SimpleNamespace) -> None:
        P.c_bar_n = 0.075
        P.h_bar_jump = 0.20
        P.h_bar_n = 0.55

    def fertility_logit_high(P: SimpleNamespace) -> None:
        P.kappa_fert = 5.0
        P.eps_fert = P.kappa_fert

    def fertility_logit_mid(P: SimpleNamespace) -> None:
        P.kappa_fert = 4.0
        P.eps_fert = P.kappa_fert

    def finance_high(P: SimpleNamespace) -> None:
        set_phi(P, 0.90)

    def fertility_mid_finance_high(P: SimpleNamespace) -> None:
        P.kappa_fert = 4.0
        P.eps_fert = P.kappa_fert
        set_phi(P, 0.90)

    def center_pull(P: SimpleNamespace) -> None:
        P.E_loc = np.asarray(P.E_loc, dtype=float).copy()
        P.r_bar = np.asarray(P.r_bar, dtype=float).copy()
        P.E_loc[1] = 0.45
        P.r_bar[1] = 0.055

    def center_mild(P: SimpleNamespace) -> None:
        P.E_loc = np.asarray(P.E_loc, dtype=float).copy()
        P.r_bar = np.asarray(P.r_bar, dtype=float).copy()
        P.E_loc[1] = 0.25
        P.r_bar[1] = 0.055

    def beta_low(P: SimpleNamespace) -> None:
        P.beta = 0.965
        P.rho = 1.0 / P.beta - 1.0
        P.rho_hat = P.rho

    def beta_mid_finance_high(P: SimpleNamespace) -> None:
        P.beta = 0.972
        P.rho = 1.0 / P.beta - 1.0
        P.rho_hat = P.rho
        set_phi(P, 0.90)

    def fertility_mid_finance_beta_mid(P: SimpleNamespace) -> None:
        P.beta = 0.972
        P.rho = 1.0 / P.beta - 1.0
        P.rho_hat = P.rho
        P.kappa_fert = 4.0
        P.eps_fert = P.kappa_fert
        set_phi(P, 0.90)

    def alpha_cons_high(P: SimpleNamespace) -> None:
        P.alpha_cons = 0.80

    def alpha_high_fertility_mid_finance(P: SimpleNamespace) -> None:
        P.alpha_cons = 0.80
        P.kappa_fert = 4.0
        P.eps_fert = P.kappa_fert
        set_phi(P, 0.90)

    def bequest_child_low(P: SimpleNamespace) -> None:
        P.theta0 = 0.30
        P.theta_n = 0.30

    def combo_soft(P: SimpleNamespace) -> None:
        P.beta = 0.965
        P.rho = 1.0 / P.beta - 1.0
        P.rho_hat = P.rho
        P.psi_child = 0.12
        P.c_bar_n = 0.075
        P.h_bar_jump = 0.20
        P.h_bar_n = 0.55
        P.kappa_fert = 5.0
        P.eps_fert = P.kappa_fert
        set_phi(P, 0.90)
        P.E_loc = np.asarray(P.E_loc, dtype=float).copy()
        P.r_bar = np.asarray(P.r_bar, dtype=float).copy()
        P.E_loc[1] = 0.45
        P.r_bar[1] = 0.055
        P.theta0 = 0.30
        P.theta_n = 0.30

    return {
        "fertility_utility_max": fertility_utility_max,
        "child_cost_low": child_cost_low,
        "fertility_logit_high": fertility_logit_high,
        "fertility_logit_mid": fertility_logit_mid,
        "finance_high": finance_high,
        "fertility_mid_finance_high": fertility_mid_finance_high,
        "center_pull": center_pull,
        "center_mild": center_mild,
        "beta_low": beta_low,
        "beta_mid_finance_high": beta_mid_finance_high,
        "fertility_mid_finance_beta_mid": fertility_mid_finance_beta_mid,
        "alpha_cons_high": alpha_cons_high,
        "alpha_high_fertility_mid_finance": alpha_high_fertility_mid_finance,
        "bequest_child_low": bequest_child_low,
        "combo_soft": combo_soft,
    }


def parse_cases(raw: str) -> list[str]:
    known = case_mutators()
    cases = [token.strip() for token in raw.split(",") if token.strip()]
    bad = [case for case in cases if case not in known]
    if bad:
        raise ValueError(f"unknown cases {bad}; known cases are {sorted(known)}")
    return cases


def run_case(case: str, args: argparse.Namespace) -> dict[str, object]:
    mutators = case_mutators()
    P, targets, weights = prepare_parameters(args.nb, args.max_iter_eq, args.tol_eq)
    P = configure_closure(P, args)
    mutators[case](P)
    params = record_param_snapshot(P)
    b_grid = make_grid(P)
    t0 = time.perf_counter()
    try:
        sol, P_out, p_eq = solve_equilibrium_hank_z(
            BENCHMARK_P,
            P,
            b_grid,
            nz=args.nz,
            rho_z=args.rho_z,
            sigma_z=args.sigma_z,
            verbose=not args.quiet,
        )
        elapsed = time.perf_counter() - t0
        moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        moments_ns.inv_pop_share_C = float(sol.pop_share[1])
        moments_ns.inv_rent_ratio_C_over_P = float(
            (P_out.user_cost_rate * p_eq[1]) / max(P_out.user_cost_rate * p_eq[0], 1e-12)
        )
        moments = namespace_float_dict(moments_ns)
        loss = compute_smm_loss(moments, targets, weights)
        scale = getattr(sol, "accounting_scale", SimpleNamespace())
        return {
            "case": case,
            "ok": True,
            "elapsed_sec": float(elapsed),
            "loss": float(loss),
            "accepted": bool(sol.timings.get("accepted", False)),
            "strict_converged": bool(sol.timings.get("strict_converged", False)),
            "convergence_reason": sol.timings.get("convergence_reason"),
            "iterations_completed": int(sol.timings.get("iterations_completed", 0)),
            "best_eq_error": float(sol.timings.get("best_eq_error", np.nan)),
            "final_eq_error": float(sol.timings.get("final_eq_error", np.nan)),
            "p_eq": [float(x) for x in p_eq],
            "scale_factor": float(getattr(scale, "scale_factor", np.nan)),
            "finite_stationary_scale": bool(getattr(scale, "finite_stationary_scale", False)),
            "city_entry_prob_total": float(getattr(scale, "city_entry_prob_total", np.nan)),
            "outside_entry_prob": float(getattr(scale, "outside_entry_prob", np.nan)),
            "outside_entry_flow": float(getattr(scale, "outside_entry_flow", np.nan)),
            "outside_value": float(getattr(scale, "outside_value", np.nan)),
            "params": params,
            "moments": moments,
            "targets": targets,
            "timings": dict(getattr(sol, "timings", {})),
        }
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        return {
            "case": case,
            "ok": False,
            "elapsed_sec": float(elapsed),
            "error": repr(exc),
            "params": params,
        }


def flat_row(row: dict[str, object]) -> dict[str, object]:
    moments = row.get("moments", {}) if isinstance(row.get("moments"), dict) else {}
    params = row.get("params", {}) if isinstance(row.get("params"), dict) else {}
    out: dict[str, object] = {
        "case": row.get("case"),
        "ok": row.get("ok"),
        "elapsed_sec": row.get("elapsed_sec"),
        "loss": row.get("loss"),
        "accepted": row.get("accepted"),
        "strict_converged": row.get("strict_converged"),
        "convergence_reason": row.get("convergence_reason"),
        "iterations_completed": row.get("iterations_completed"),
        "final_eq_error": row.get("final_eq_error"),
        "p_P": row.get("p_eq", [np.nan, np.nan])[0] if isinstance(row.get("p_eq"), list) else np.nan,
        "p_C": row.get("p_eq", [np.nan, np.nan])[1] if isinstance(row.get("p_eq"), list) else np.nan,
        "scale_factor": row.get("scale_factor"),
        "finite_stationary_scale": row.get("finite_stationary_scale"),
        "city_entry_prob_total": row.get("city_entry_prob_total"),
        "outside_entry_prob": row.get("outside_entry_prob"),
        "outside_entry_flow": row.get("outside_entry_flow"),
        "outside_value": row.get("outside_value"),
    }
    for name, value in params.items():
        out[f"param_{name}"] = value
    for name in MOMENT_ORDER:
        out[name] = moments.get(name, np.nan)
    return out


def write_summary_csv(path: Path, rows: list[dict[str, object]]) -> None:
    flat = [flat_row(row) for row in rows]
    fieldnames = list(flat[0].keys()) if flat else []
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(flat)


def write_moment_csv(path: Path, baseline: dict[str, object], rows: list[dict[str, object]]) -> None:
    base_moments = deepcopy(baseline["moments"])
    rows_out: list[dict[str, object]] = []
    targets = rows[0].get("targets", {}) if rows else {}
    for row in rows:
        moments = row.get("moments", {}) if isinstance(row.get("moments"), dict) else {}
        for name in MOMENT_ORDER:
            target = float(targets.get(name, np.nan)) if isinstance(targets, dict) else np.nan
            model = float(moments.get(name, np.nan))
            base = float(base_moments.get(name, np.nan))
            rows_out.append(
                {
                    "case": row.get("case"),
                    "moment": name,
                    "target": target,
                    "baseline": base,
                    "model": model,
                    "baseline_gap": base - target,
                    "model_gap": model - target,
                    "abs_gap_improvement": abs(base - target) - abs(model - target),
                }
            )
    with path.open("w", newline="") as f:
        fieldnames = [
            "case",
            "moment",
            "target",
            "baseline",
            "model",
            "baseline_gap",
            "model_gap",
            "abs_gap_improvement",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)


def write_report(path: Path, baseline: dict[str, object], rows: list[dict[str, object]]) -> None:
    base_moments = baseline["moments"]
    base_loss = float(baseline["loss"])
    key_moments = [
        "tfr",
        "childless_rate",
        "mean_age_first_birth",
        "tfr_gradient",
        "own_rate",
        "own_gradient",
        "young_liquid_wealth_to_income",
        "center_share_nonparents",
        "old_age_own_rate",
        "old_age_parent_childless_gap",
    ]
    target_map = rows[0].get("targets", {}) if rows else {}

    summary_table = "\n".join(
        "| {case} | {ok} | {it} | {err:.3g} | {sec:.2f} | {loss:.2f} | {dloss:.2f} | {improved}/19 |".format(
            case=row.get("case"),
            ok=row.get("ok"),
            it=int(row.get("iterations_completed", 0) or 0),
            err=float(row.get("final_eq_error", np.nan)),
            sec=float(row.get("elapsed_sec", np.nan)),
            loss=float(row.get("loss", np.nan)),
            dloss=float(row.get("loss", np.nan)) - base_loss,
            improved=sum(
                1
                for name in MOMENT_ORDER
                if abs(float(row.get("moments", {}).get(name, np.nan)) - float(target_map.get(name, np.nan)))
                < abs(float(base_moments.get(name, np.nan)) - float(target_map.get(name, np.nan)))
            ),
        )
        for row in rows
    )

    moment_lines = []
    for row in rows:
        moments = row.get("moments", {})
        moment_lines.append(f"### {row.get('case')}")
        moment_lines.append("| Moment | Target | Baseline | Probe | Direction |")
        moment_lines.append("|---|---:|---:|---:|---:|")
        for name in key_moments:
            target = float(target_map.get(name, np.nan))
            base = float(base_moments.get(name, np.nan))
            probe = float(moments.get(name, np.nan))
            improve = abs(base - target) - abs(probe - target)
            moment_lines.append(
                f"| `{name}` | {target:.3f} | {base:.3f} | {probe:.3f} | {improve:+.3f} |"
            )
        moment_lines.append("")

    best = min((row for row in rows if row.get("ok")), key=lambda r: float(r.get("loss", np.inf)), default=None)
    best_text = "No successful probe."
    if best is not None:
        best_text = (
            f"Lowest-loss probe is `{best['case']}` with loss "
            f"`{float(best['loss']):.6g}` versus baseline `{base_loss:.6g}`."
        )

    text = f"""# V5 Nz=5 Directional Calibration Audit

This is a small same-grid calibration probe for the isolated V5
benchmark-normalized outside-option closure. It is not an optimizer. The goal
is to test whether the bad `Nz=5` moments are locally movable in the right
direction before running a larger search.

Baseline reference: `income_mortgage_risk_v5_hank_z_outside_closure_nz5.log`
with loss `{base_loss:.6g}`.

## Summary

| Case | OK | iters | GE err | seconds | loss | loss minus baseline | improved moments |
|---|---:|---:|---:|---:|---:|---:|---:|
{summary_table}

{best_text}

## Key Moment Directions

Positive direction entries mean the probe reduced the absolute target gap
relative to the `Nz=5` baseline.

{chr(10).join(moment_lines)}
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cases",
        default="fertility_utility_max,child_cost_low,fertility_logit_high,finance_high,center_pull,beta_low,bequest_child_low,combo_soft",
    )
    parser.add_argument("--nb", type=int, default=30)
    parser.add_argument("--nz", type=int, default=5)
    parser.add_argument("--rho-z", type=float, default=0.95)
    parser.add_argument("--sigma-z", type=float, default=0.35)
    parser.add_argument("--target-city-entry-prob", type=float, default=0.9)
    parser.add_argument("--kappa-entry", type=float, default=1_000_000.0)
    parser.add_argument("--max-iter-eq", type=int, default=35)
    parser.add_argument("--tol-eq", type=float, default=5e-4)
    parser.add_argument("--output-stem", default=OUTPUT_STEM)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    baseline = load_baseline(root)
    rows = []
    for case in parse_cases(args.cases):
        if not args.quiet:
            print(f"running {case}", flush=True)
        rows.append(run_case(case, args))

    output_stem = str(args.output_stem)
    write_summary_csv(root / f"{output_stem}.csv", rows)
    write_moment_csv(root / f"{output_stem}_moments.csv", baseline, rows)
    (root / f"{output_stem}.json").write_text(json.dumps({"baseline": baseline, "rows": rows}, indent=2, default=str))
    report_name = "REPORT_V5_NZ5_DIRECTIONAL_CALIBRATION.md"
    if output_stem != OUTPUT_STEM:
        report_name = f"{output_stem}.md"
    write_report(root / report_name, baseline, rows)


if __name__ == "__main__":
    main()
