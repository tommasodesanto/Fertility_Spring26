#!/usr/bin/env python3
"""V5 Branch 1 prototype: HANK-z GE with outside-option scale closure."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from dt_cp_model.direct_calibration import compute_smm_loss
from dt_cp_model.objective import extract_moments
from dt_cp_model.parameters import apply_overrides, finalize_location_choice_spec, setup_parameters
from dt_cp_model.solver import (
    make_grid,
    solve_equilibrium_hank_z,
)

from run_income_mortgage_risk_v2_scenarios import (
    BENCHMARK_MODEL,
    BENCHMARK_P,
    MOMENT_ORDER,
    build_base_parameters,
    namespace_float_dict,
)
from run_income_mortgage_risk_v3_hank_z import write_diagnostics as write_hank_z_diagnostics


def prepare_parameters(nb: int, max_iter_eq: int, tol_eq: float) -> tuple[SimpleNamespace, dict[str, float], dict[str, float]]:
    P_override, targets, weights = build_base_parameters(nb, force_full=False)
    P = setup_parameters()
    P = apply_overrides(P, P_override)
    P = finalize_location_choice_spec(P)
    P.solve_mode = "ge"
    P.Nb = int(nb)
    P.max_iter_eq = int(max_iter_eq)
    P.tol_eq = float(tol_eq)
    P.collect_ge_trace = True
    return P, targets, weights


def write_results(path: Path, targets: dict[str, float], model: dict[str, float]) -> None:
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
                    "notes": "V5 full-equilibrium HANK-z with outside-option scale closure",
                }
            )


def write_trace(path: Path, sol: SimpleNamespace) -> None:
    if not hasattr(sol, "ge_trace"):
        return
    fieldnames = [
        "iter",
        "err",
        "err_p",
        "err_e",
        "p_P",
        "p_C",
        "entry_share_P",
        "entry_share_C",
        "scale_factor",
        "implied_total_population",
        "city_entry_prob_total",
        "outside_entry_prob",
        "outside_entry_flow",
        "scale_denominator",
        "scale_residual",
        "own_rate",
        "tfr",
        "pop_P",
        "pop_C",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in sol.ge_trace:
            writer.writerow(
                {
                    "iter": rec.get("iter"),
                    "err": rec.get("err"),
                    "err_p": rec.get("err_p"),
                    "err_e": rec.get("err_e"),
                    "p_P": rec.get("p", [np.nan, np.nan])[0],
                    "p_C": rec.get("p", [np.nan, np.nan])[1],
                    "entry_share_P": rec.get("entry_shares", [np.nan, np.nan])[0],
                    "entry_share_C": rec.get("entry_shares", [np.nan, np.nan])[1],
                    "scale_factor": rec.get("scale_factor"),
                    "implied_total_population": rec.get("implied_total_population"),
                    "city_entry_prob_total": rec.get("city_entry_prob_total"),
                    "outside_entry_prob": rec.get("outside_entry_prob"),
                    "outside_entry_flow": rec.get("outside_entry_flow"),
                    "scale_denominator": rec.get("scale_denominator"),
                    "scale_residual": rec.get("scale_residual"),
                    "own_rate": rec.get("own_rate"),
                    "tfr": rec.get("TFR"),
                    "pop_P": rec.get("pop_share", [np.nan, np.nan])[0],
                    "pop_C": rec.get("pop_share", [np.nan, np.nan])[1],
                }
            )


def write_closure_diagnostics(path: Path, sol: SimpleNamespace) -> None:
    rows: list[dict[str, object]] = []
    scale = getattr(sol, "accounting_scale", None)
    if scale is not None:
        for name in [
            "benchmark_normalized",
            "finite_stationary_scale",
            "target_city_entry_prob",
            "outside_value",
            "outside_entry_flow",
            "local_birth_entry_weight",
            "kappa_entry",
            "city_entry_prob_total",
            "outside_entry_prob",
            "denominator",
            "scale_factor",
            "implied_total_population",
            "target_total_population",
            "reference_entry_total",
            "reference_mature_cityborn_flow",
            "entry_per_unit_scale",
            "mature_cityborn_per_unit_scale",
            "implied_entry_total",
            "implied_potential_entrant_mass",
            "implied_outside_entry_mass",
            "stationary_entry_residual",
            "stationary_entry_relative_residual",
        ]:
            rows.append({"stage": "benchmark_normalization", "name": name, "value": getattr(scale, name, np.nan)})
        for i, value in enumerate(scale.entry_values):
            rows.append({"stage": "benchmark_normalization", "name": f"entry_value_loc_{i}", "value": float(value)})
        for i, value in enumerate(scale.city_entry_prob):
            rows.append({"stage": "benchmark_normalization", "name": f"city_entry_prob_loc_{i}", "value": float(value)})
        for i, value in enumerate(scale.conditional_entry_shares):
            rows.append({"stage": "benchmark_normalization", "name": f"conditional_entry_share_loc_{i}", "value": float(value)})
        for i, value in enumerate(scale.implied_housing_demand):
            rows.append({"stage": "benchmark_normalization", "name": f"implied_housing_demand_loc_{i}", "value": float(value)})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["stage", "name", "value"])
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    sol: SimpleNamespace,
    P: SimpleNamespace,
    p_eq: np.ndarray,
    moments: dict[str, float],
    targets: dict[str, float],
    loss: float,
    elapsed: float,
) -> None:
    hz = sol.hank_z
    scale = getattr(sol, "accounting_scale", None)
    accepted = bool(sol.timings.get("accepted", False))
    finite_scale = bool(getattr(scale, "finite_stationary_scale", False))
    verdict = "yellow" if accepted and finite_scale else "red"
    table = "\n".join(
        f"| `{name}` | {targets.get(name, np.nan):.3f} | {BENCHMARK_MODEL.get(name, np.nan):.3f} | {moments.get(name, np.nan):.3f} |"
        for name in MOMENT_ORDER
    )
    text = f"""# Income Risk V5 HANK-z Outside-Option Closure Report

Verdict: **{verdict}**

## What Changed

This run keeps the structural HANK earnings state \(z\) and switches the copied
GE loop from the renewal-valve closure to the paper-facing outside-option scale
closure. The closure implements
\\[
S E_0(p)=q^E(p)\\left[M+S B_0(p)\\right],
\\qquad
S(p)=\\frac{{q^E(p)M}}{{E_0(p)-q^E(p)B_0(p)}}.
\\]
For the benchmark steady state, \(S=1\) is imposed inside the GE loop. At each
candidate composition, the script computes \(E_0(p)\) and \(B_0(p)\),
calibrates \(\\bar W^E\) to the target \(q^E\) at current entry values, and sets
\[
M = E_0(p)/q^E(p)-B_0(p).
\]
This replaces the earlier outer normalization pass. The reported final \(M\)
and \(\\bar W^E\) are the objects to hold fixed in counterfactuals.

## Benchmark Normalization

- target \(q^E\): `{float(getattr(scale, "target_city_entry_prob", np.nan)):.4f}`
- calibrated outside value \(\\bar W^E\): `{float(getattr(scale, "outside_value", np.nan)):.6g}`
- entry logit scale \(\kappa_E\): `{float(getattr(scale, "kappa_entry", np.nan)):.6g}`
- residual outside-born flow \(M\): `{float(getattr(scale, "outside_entry_flow", np.nan)):.6g}`
- entry flow per unit scale \(E_0\): `{float(getattr(scale, "entry_per_unit_scale", np.nan)):.6g}`
- mature city-born flow per unit scale \(B_0\): `{float(getattr(scale, "mature_cityborn_per_unit_scale", np.nan)):.6g}`
- accounting scale factor \(S\): `{float(getattr(scale, "scale_factor", np.nan)):.6g}`
- scale residual: `{float(getattr(scale, "stationary_entry_residual", np.nan)):.3e}`

## Final GE Status

- closure: `{P.population_closure}`
- accepted: `{accepted}`
- strict converged: `{bool(sol.timings.get("strict_converged", False))}`
- convergence reason: `{sol.timings.get("convergence_reason")}`
- iterations completed: `{int(sol.timings.get("iterations_completed", 0))}`
- best equilibrium error: `{float(sol.timings.get("best_eq_error", np.nan)):.6g}`
- final equilibrium error: `{float(sol.timings.get("final_eq_error", np.nan)):.6g}`
- prices: `{[float(x) for x in p_eq]}`
- \(z\) states: `{len(hz.z_grid)}`
- \(b\) states: `{P.Nb}`
- final scale factor \(S\): `{float(getattr(scale, "scale_factor", np.nan)):.6g}`
- final city-entry probability \(q^E\): `{float(getattr(scale, "city_entry_prob_total", np.nan)):.6g}`
- final outside probability: `{float(getattr(scale, "outside_entry_prob", np.nan)):.6g}`
- finite scale: `{finite_scale}`
- elapsed seconds: `{elapsed:.2f}`
- SMM loss against live targets: `{loss:.6g}`

## Normalization Passes

No outer normalization pass was used. The benchmark normalization is imposed
directly inside each GE iteration.

## Moment Table

| Moment | Target | Benchmark | V5 outside closure |
|---|---:|---:|---:|
{table}

## Read

Auditing against `latex/model_writeup.tex` and the copied solver logic, the
outer re-solve was unnecessarily complicated for the benchmark. The writeup
normalizes the baseline at \(S=1\) and sets \(M\) residually from the baseline
objects. Numerically, the cleaner benchmark implementation is to impose that
residual normalization at the current candidate equilibrium/composition and use
normalized housing demand. A fixed-\(M\) scale formula is the counterfactual
object, not the object that should move the benchmark away from its own
normalization.

The run is still yellow because the un-recalibrated economics are not
acceptable yet: fertility is too low and too late, geography and ownership
gradients flip sign, and liquid wealth is too high. The mortgage account state
is also still absent from this HANK-z-only branch.

## Closure Lessons

The entry margin needs its own scale \(\kappa_E\). Using the incumbent
within-city location scale \(\kappa_\ell\) made the outside probability jump to
zero or one because entry values are lifetime-utility objects with very large
levels. The accepted run therefore uses an explicitly supplied
\(\kappa_E={float(getattr(scale, "kappa_entry", np.nan)):.6g}\). The Rouwenhorst grid also requires the
true stationary Markov distribution; this script now reports the binomial
stationary weights for \(N_z={len(hz.z_grid)}\), not uniform smoke weights.
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nb", type=int, default=30)
    parser.add_argument("--nz", type=int, default=7)
    parser.add_argument("--rho-z", type=float, default=0.95)
    parser.add_argument("--sigma-z", type=float, default=0.35)
    parser.add_argument("--target-city-entry-prob", type=float, default=0.9)
    parser.add_argument("--kappa-entry", type=float, default=None)
    parser.add_argument("--local-birth-entry-weight", type=float, default=1.0)
    parser.add_argument("--baseline-max-iter-eq", type=int, default=35)
    parser.add_argument("--max-iter-eq", type=int, default=45)
    parser.add_argument("--tol-eq", type=float, default=5e-4)
    parser.add_argument("--normalization-passes", type=int, default=2)
    parser.add_argument("--scale-target-tol", type=float, default=1e-5)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    log_path = out_dir / "income_mortgage_risk_v5_hank_z_outside_closure.log"
    t0 = time.perf_counter()
    try:
        P, targets, weights = prepare_parameters(args.nb, args.max_iter_eq, args.tol_eq)
        P.population_closure = "outside_option_benchmark_normalized"
        P.kappa_entry = (
            float(args.kappa_entry) if args.kappa_entry is not None else float(getattr(P, "kappa_loc", 1.0))
        )
        P.local_birth_entry_weight = float(args.local_birth_entry_weight)
        P.target_city_entry_prob = float(args.target_city_entry_prob)
        P.calibrate_outside_value_to_entry_prob = True
        P.outside_benchmark_target_total_population = float(getattr(P, "N_target", 1.0))
        P.allow_uncalibrated_outside_value = False

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
        scale = getattr(sol, "accounting_scale", None)
        if scale is None or not bool(getattr(scale, "finite_stationary_scale", False)):
            raise RuntimeError("benchmark-normalized outside-option closure produced a nonfinite scale")

        moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        moments_ns.inv_pop_share_C = float(sol.pop_share[1])
        moments_ns.inv_rent_ratio_C_over_P = float((P_out.user_cost_rate * p_eq[1]) / max(P_out.user_cost_rate * p_eq[0], 1e-12))
        moments = namespace_float_dict(moments_ns)
        loss = compute_smm_loss(moments, targets, weights)
        elapsed = time.perf_counter() - t0

        write_results(out_dir / "results_income_mortgage_risk_v5_hank_z_outside_closure.csv", targets, moments)
        write_hank_z_diagnostics(out_dir / "diagnostics_income_mortgage_risk_v5_hank_z_outside_closure.csv", sol, P_out)
        write_closure_diagnostics(out_dir / "diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_closure.csv", sol)
        write_trace(out_dir / "diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_trace.csv", sol)
        write_report(
            out_dir / "REPORT_V5_HANK_Z_OUTSIDE_CLOSURE.md",
            sol,
            P_out,
            p_eq,
            moments,
            targets,
            loss,
            elapsed,
        )
        log_path.write_text(
            json.dumps(
                {
                    "ok": True,
                    "elapsed_sec": elapsed,
                    "loss": loss,
                    "normalization_implementation": "inside_ge_loop",
                    "normalization_passes_requested": int(args.normalization_passes),
                    "scale_target_tol_requested": float(args.scale_target_tol),
                    "p_eq": [float(x) for x in p_eq],
                    "z_grid": [float(x) for x in sol.hank_z.z_grid],
                    "stationary_z": [float(x) for x in sol.hank_z.stationary_z],
                    "benchmark_accounting_scale": vars(scale),
                    "moments": moments,
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
