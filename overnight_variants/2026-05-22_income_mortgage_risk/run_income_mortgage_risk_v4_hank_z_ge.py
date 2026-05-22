#!/usr/bin/env python3
"""V4 Branch 1 prototype: full equilibrium with structural HANK-z risk."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

from dt_cp_model.direct_calibration import compute_smm_loss
from dt_cp_model.objective import extract_moments
from dt_cp_model.parameters import apply_overrides, finalize_location_choice_spec, setup_parameters
from dt_cp_model.solver import make_grid, solve_equilibrium_hank_z

from run_income_mortgage_risk_v2_scenarios import (
    BENCHMARK_MODEL,
    BENCHMARK_P,
    MOMENT_ORDER,
    build_base_parameters,
    namespace_float_dict,
)
from run_income_mortgage_risk_v3_hank_z import write_diagnostics as write_hank_z_diagnostics


def write_ge_results(path: Path, targets: dict[str, float], model: dict[str, float]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["moment", "target", "model", "benchmark_model", "notes"])
        writer.writeheader()
        for name in MOMENT_ORDER:
            writer.writerow(
                {
                    "moment": name,
                    "target": targets.get(name, ""),
                    "model": model.get(name, np.nan),
                    "benchmark_model": BENCHMARK_MODEL.get(name, ""),
                    "notes": "V4 full-equilibrium structural z-state solve; mortgage account not structural",
                }
            )


def write_ge_diagnostics(path: Path, sol, P) -> None:
    write_hank_z_diagnostics(path, sol, P)


def write_ge_trace(path: Path, sol) -> None:
    if not hasattr(sol, "ge_trace"):
        return
    rows = []
    for rec in sol.ge_trace:
        rows.append(
            {
                "diagnostic": "ge_trace",
                "iter": rec.get("iter"),
                "err": rec.get("err"),
                "err_p": rec.get("err_p"),
                "err_e": rec.get("err_e"),
                "p_P": rec.get("p", [np.nan, np.nan])[0],
                "p_C": rec.get("p", [np.nan, np.nan])[1],
                "own_rate": rec.get("own_rate"),
                "tfr": rec.get("TFR"),
                "pop_P": rec.get("pop_share", [np.nan, np.nan])[0],
                "pop_C": rec.get("pop_share", [np.nan, np.nan])[1],
            }
        )
    with path.open("w", newline="") as f:
        fieldnames = [
            "diagnostic",
            "iter",
            "err",
            "err_p",
            "err_e",
            "p_P",
            "p_C",
            "own_rate",
            "tfr",
            "pop_P",
            "pop_C",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_report(path: Path, sol, P, p_eq: np.ndarray, moments: dict[str, float], targets: dict[str, float], loss: float, elapsed: float) -> None:
    hz = sol.hank_z
    diag = hz.diagnostics
    owner_vals = [float(x["own_rate"]) for x in diag.owner_by_z]
    fert_vals = [float(x["fertility_choice"]) for x in diag.fertility_by_z]
    own_spread = max(owner_vals) - min(owner_vals) if owner_vals else np.nan
    fert_spread = max(fert_vals) - min(fert_vals) if fert_vals else np.nan
    accepted = bool(sol.timings.get("accepted", False))
    verdict = "yellow" if accepted else "red"
    if accepted and moments.get("own_gradient", 0.0) < 0:
        verdict = "yellow"
    table = "\n".join(
        f"| `{name}` | {targets.get(name, np.nan):.3f} | {BENCHMARK_MODEL.get(name, np.nan):.3f} | {moments.get(name, np.nan):.3f} |"
        for name in MOMENT_ORDER
    )
    text = f"""# Income Risk V4 HANK-z Full Equilibrium Report

Verdict: **{verdict}**

## What This Is

This is the full-equilibrium correction to V3. The copied branch solves prices
and entry shares with the structural HANK earnings state \(z\) inside the
household problem and forward distribution. It is no longer a fixed-price PE
test. The household state is \((b,d,i,a,n,s,z)\), with
\(y_{{ia}}(z)=(1-\\tau_{{pay}})w_i e_a\exp(z)\) for working ages and a common
retirement mapping.

The mortgage-account state \(\mu\) is still not structural. This run answers
the narrower question: can the classic HANK \(z\)-state version clear the
copied model's equilibrium loop on a coarse grid?

## Equilibrium Status

- accepted: `{accepted}`
- strict converged: `{bool(sol.timings.get("strict_converged", False))}`
- convergence reason: `{sol.timings.get("convergence_reason")}`
- iterations completed: `{int(sol.timings.get("iterations_completed", 0))}`
- best equilibrium error: `{float(sol.timings.get("best_eq_error", np.nan)):.6g}`
- best iteration: `{int(sol.timings.get("best_eq_iter", 0))}`
- final equilibrium error: `{float(sol.timings.get("final_eq_error", np.nan)):.6g}`
- prices: `{[float(x) for x in p_eq]}`
- \(z\) states: `{len(hz.z_grid)}`
- \(b\) states: `{P.Nb}`
- elapsed seconds: `{elapsed:.2f}`
- runtime category: `{sol.timings.get("runtime_category", "expensive")}`
- SMM loss against live targets: `{loss:.6g}`

## HANK Diagnostics

| Object | Value |
|---|---:|
| ownership spread across \(z\) | {own_spread:.3f} |
| fertility-choice spread across \(z\) | {fert_spread:.3f} |
| aggregate TFR | {moments.get("tfr", np.nan):.3f} |
| aggregate ownership | {moments.get("own_rate", np.nan):.3f} |
| young liquid wealth / income | {moments.get("young_liquid_wealth_to_income", np.nan):.3f} |

## Moment Table

| Moment | Target | Benchmark | V4 GE |
|---|---:|---:|---:|
{table}

## Read

This is the equilibrium object the branch needed: \(z\) is a real Markov state
and the price loop moves with the \(z\)-state distribution. The result is still
yellow, not green, because the coarse grid is not recalibrated and the mortgage
account is absent. The next live implementation should add only this HANK-z
core first, then decide separately whether a compact account state is worth the
extra state-space cost.
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nb", type=int, default=30)
    parser.add_argument("--nz", type=int, default=3, choices=[3, 5])
    parser.add_argument("--rho-z", type=float, default=0.82)
    parser.add_argument("--sigma-z", type=float, default=0.28)
    parser.add_argument("--max-iter-eq", type=int, default=35)
    parser.add_argument("--tol-eq", type=float, default=5e-4)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    log_path = out_dir / "income_mortgage_risk_v4_hank_z_ge.log"
    t0 = time.perf_counter()
    try:
        P_override, targets, weights = build_base_parameters(args.nb, force_full=False)
        P = setup_parameters()
        P = apply_overrides(P, P_override)
        P = finalize_location_choice_spec(P)
        P.solve_mode = "ge"
        P.Nb = int(args.nb)
        P.max_iter_eq = int(args.max_iter_eq)
        P.tol_eq = float(args.tol_eq)
        P.collect_ge_trace = True
        b_grid = make_grid(P)
        sol, P_out, p_eq = solve_equilibrium_hank_z(
            BENCHMARK_P,
            P,
            b_grid,
            nz=args.nz,
            rho_z=args.rho_z,
            sigma_z=args.sigma_z,
            verbose=not args.quiet,
        )
        moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        moments_ns.inv_pop_share_C = float(sol.pop_share[1])
        moments_ns.inv_rent_ratio_C_over_P = float((P_out.user_cost_rate * p_eq[1]) / max(P_out.user_cost_rate * p_eq[0], 1e-12))
        moments = namespace_float_dict(moments_ns)
        loss = compute_smm_loss(moments, targets, weights)
        elapsed = time.perf_counter() - t0
        write_ge_results(out_dir / "results_income_mortgage_risk_v4_hank_z_ge.csv", targets, moments)
        write_ge_diagnostics(out_dir / "diagnostics_income_mortgage_risk_v4_hank_z_ge.csv", sol, P_out)
        write_ge_trace(out_dir / "diagnostics_income_mortgage_risk_v4_hank_z_ge_trace.csv", sol)
        write_report(out_dir / "REPORT_V4_HANK_Z_GE.md", sol, P_out, p_eq, moments, targets, loss, elapsed)
        log_path.write_text(
            json.dumps(
                {
                    "ok": True,
                    "elapsed_sec": elapsed,
                    "loss": loss,
                    "moments": moments,
                    "p_eq": [float(x) for x in p_eq],
                    "z_grid": [float(x) for x in sol.hank_z.z_grid],
                    "stationary_z": [float(x) for x in sol.hank_z.stationary_z],
                    "timings": sol.timings,
                    "ge_trace": getattr(sol, "ge_trace", []),
                },
                indent=2,
                default=str,
            )
        )
    except Exception as exc:
        log_path.write_text(json.dumps({"ok": False, "error": repr(exc)}, indent=2))
        raise


if __name__ == "__main__":
    main()
