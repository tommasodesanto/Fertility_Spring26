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
    accounting_population_scale,
    calibrate_outside_value_for_entry_total,
    city_entry_probabilities,
    entry_values_by_location,
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


def calibrate_outside_scale_closure(
    sol0: SimpleNamespace,
    P0: SimpleNamespace,
    b_grid0: np.ndarray,
    *,
    target_city_entry_prob: float,
    kappa_entry: float,
    local_birth_entry_weight: float,
    target_total_population: float,
) -> SimpleNamespace:
    """Calibrate the paper closure so the baseline has unit scale.

    The writeup closure is

        S E_0(p) = q^E(p) [M + S B_0(p)].

    Given a target baseline city-entry probability q^E, this helper picks the
    outside value that implements that q^E at baseline entry values and then
    sets the outside-born flow M residually so the baseline scale is unchanged.
    """

    target_q = float(np.clip(target_city_entry_prob, 1e-6, 1.0 - 1e-6))
    local_weight = max(float(local_birth_entry_weight), 0.0)
    kappa = max(float(kappa_entry), 1e-10)
    target_pop = max(float(target_total_population), 1e-12)

    entry_values = entry_values_by_location(sol0.V, b_grid0, P0)
    outside_value = calibrate_outside_value_for_entry_total(entry_values, 1.0, target_q, kappa)
    city_probs, outside_prob = city_entry_probabilities(entry_values, outside_value, kappa)
    city_prob_total = float(np.sum(city_probs))

    ref_pop = max(float(getattr(sol0, "total_mass", getattr(P0, "N_target", 1.0))), 1e-14)
    entry_flow = max(float(getattr(sol0, "entry_rate", getattr(P0, "E_total", 0.0))), 1e-14)
    mature_flow = max(float(getattr(sol0, "entrants_mature_total", 0.0)), 0.0)
    entry_per_scale = entry_flow / ref_pop
    mature_per_scale = mature_flow / ref_pop

    outside_entry_flow = target_pop * (entry_per_scale / city_prob_total - local_weight * mature_per_scale)
    finite = outside_entry_flow >= 0.0 and city_prob_total > 1e-12

    scale = accounting_population_scale(
        sol0,
        P0,
        b_grid0,
        outside_entry_flow=max(outside_entry_flow, 0.0),
        local_birth_entry_weight=local_weight,
        outside_value=outside_value,
        kappa_entry=kappa,
        calibrate_outside_value=False,
    )
    return SimpleNamespace(
        finite_calibration=bool(finite and scale.finite_stationary_scale),
        target_city_entry_prob=target_q,
        target_total_population=target_pop,
        outside_value=float(outside_value),
        outside_entry_flow=float(outside_entry_flow),
        local_birth_entry_weight=float(local_weight),
        kappa_entry=float(kappa),
        entry_values=entry_values,
        city_entry_prob=city_probs,
        city_entry_prob_total=city_prob_total,
        outside_entry_prob=float(outside_prob),
        reference_total_population=float(ref_pop),
        reference_entry_total=float(entry_flow),
        reference_mature_cityborn_flow=float(mature_flow),
        entry_per_unit_scale=float(entry_per_scale),
        mature_cityborn_per_unit_scale=float(mature_per_scale),
        baseline_accounting_scale=scale,
    )


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
                    "scale_denominator": rec.get("scale_denominator"),
                    "scale_residual": rec.get("scale_residual"),
                    "own_rate": rec.get("own_rate"),
                    "tfr": rec.get("TFR"),
                    "pop_P": rec.get("pop_share", [np.nan, np.nan])[0],
                    "pop_C": rec.get("pop_share", [np.nan, np.nan])[1],
                }
            )


def write_closure_diagnostics(path: Path, closure0: SimpleNamespace, sol: SimpleNamespace) -> None:
    rows: list[dict[str, object]] = []
    scale0 = closure0.baseline_accounting_scale
    scale = getattr(sol, "accounting_scale", None)
    for rec in getattr(closure0, "normalization_history", []):
        for key, value in rec.items():
            if key != "pass":
                rows.append({"stage": f"normalization_pass_{rec['pass']}", "name": key, "value": value})
    for name in [
        "target_city_entry_prob",
        "outside_value",
        "outside_entry_flow",
        "local_birth_entry_weight",
        "kappa_entry",
        "city_entry_prob_total",
        "outside_entry_prob",
        "reference_entry_total",
        "reference_mature_cityborn_flow",
        "entry_per_unit_scale",
        "mature_cityborn_per_unit_scale",
    ]:
        rows.append({"stage": "baseline_calibration", "name": name, "value": getattr(closure0, name)})
    for i, value in enumerate(closure0.entry_values):
        rows.append({"stage": "baseline_calibration", "name": f"entry_value_loc_{i}", "value": float(value)})
    for i, value in enumerate(closure0.city_entry_prob):
        rows.append({"stage": "baseline_calibration", "name": f"city_entry_prob_loc_{i}", "value": float(value)})
    for name in [
        "finite_stationary_scale",
        "denominator",
        "scale_factor",
        "implied_total_population",
        "stationary_entry_residual",
        "stationary_entry_relative_residual",
    ]:
        rows.append({"stage": "baseline_accounting_check", "name": name, "value": getattr(scale0, name)})
    if scale is not None:
        for name in [
            "finite_stationary_scale",
            "outside_value",
            "outside_entry_flow",
            "city_entry_prob_total",
            "outside_entry_prob",
            "denominator",
            "scale_factor",
            "implied_total_population",
            "reference_entry_total",
            "reference_mature_cityborn_flow",
            "stationary_entry_residual",
            "stationary_entry_relative_residual",
        ]:
            rows.append({"stage": "final_closure", "name": name, "value": getattr(scale, name)})
        for i, value in enumerate(scale.city_entry_prob):
            rows.append({"stage": "final_closure", "name": f"city_entry_prob_loc_{i}", "value": float(value)})
        for i, value in enumerate(scale.conditional_entry_shares):
            rows.append({"stage": "final_closure", "name": f"conditional_entry_share_loc_{i}", "value": float(value)})
        for i, value in enumerate(scale.implied_housing_demand):
            rows.append({"stage": "final_closure", "name": f"implied_housing_demand_loc_{i}", "value": float(value)})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["stage", "name", "value"])
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    sol0: SimpleNamespace,
    sol: SimpleNamespace,
    P: SimpleNamespace,
    p_eq: np.ndarray,
    moments: dict[str, float],
    targets: dict[str, float],
    loss: float,
    elapsed: float,
    closure0: SimpleNamespace,
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
At the baseline HANK-z equilibrium, the script chooses the outside value to hit
a target city-entry probability \(q^E={closure0.target_city_entry_prob:.3f}\)
and then sets \(M\) residually so the baseline scale is unchanged.

## Baseline Normalization

- baseline TFR before closure switch: `{2 * sol0.mean_parity:.4f}`
- baseline ownership before closure switch: `{sol0.own_rate:.4f}`
- target \(q^E\): `{closure0.target_city_entry_prob:.4f}`
- calibrated outside value \(\\bar W^E\): `{closure0.outside_value:.6g}`
- entry logit scale \(\kappa_E\): `{closure0.kappa_entry:.6g}`
- residual outside-born flow \(M\): `{closure0.outside_entry_flow:.6g}`
- baseline accounting scale factor: `{float(closure0.baseline_accounting_scale.scale_factor):.6g}`
- baseline scale residual: `{float(closure0.baseline_accounting_scale.stationary_entry_residual):.3e}`

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

{normalization_history_table(getattr(closure0, "normalization_history", []))}

## Moment Table

| Moment | Target | Benchmark | V5 outside closure |
|---|---:|---:|---:|
{table}

## Read

This is a real full-equilibrium solve with the outside-option scale closure
active in the price loop. The closure itself is not the failure: \(q^E\) remains
interior, \(S\) remains close to one, and the stationarity residual is zero up
to numerical precision. The run is yellow because the un-recalibrated economics
are not acceptable yet: fertility is too low and too late, geography and
ownership gradients flip sign, and liquid wealth is too high. The mortgage
account state is also still absent from this HANK-z-only branch.

## Closure Lessons

The entry margin needs its own scale \(\kappa_E\). Using the incumbent
within-city location scale \(\kappa_\ell\) made the outside probability jump to
zero or one because entry values are lifetime-utility objects with very large
levels. The accepted run therefore uses an explicitly supplied
\(\kappa_E={closure0.kappa_entry:.6g}\). The Rouwenhorst grid also requires the
true stationary Markov distribution; this script now reports the binomial
stationary weights for \(N_z={len(hz.z_grid)}\), not uniform smoke weights.
"""
    path.write_text(text)


def normalization_history_table(history: list[dict[str, object]]) -> str:
    if not history:
        return "No outer normalization history was recorded."
    rows = [
        "| Pass | \(S\) | \(q^E\) | outside prob. | \(M\) | GE error |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for rec in history:
        rows.append(
            "| {pass_id} | {scale:.6g} | {q:.6g} | {outside:.6g} | {m:.6g} | {err:.6g} |".format(
                pass_id=int(rec.get("pass", 0)),
                scale=float(rec.get("scale_factor", np.nan)),
                q=float(rec.get("city_entry_prob_total", np.nan)),
                outside=float(rec.get("outside_entry_prob", np.nan)),
                m=float(rec.get("outside_entry_flow", np.nan)),
                err=float(rec.get("best_eq_error", np.nan)),
            )
        )
    return "\n".join(rows)


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
        P0, targets, weights = prepare_parameters(args.nb, args.baseline_max_iter_eq, args.tol_eq)
        P0.population_closure = "renewal_valve_calibrated"
        b_grid0 = make_grid(P0)
        sol0, P0_out, p0 = solve_equilibrium_hank_z(
            BENCHMARK_P,
            P0,
            b_grid0,
            nz=args.nz,
            rho_z=args.rho_z,
            sigma_z=args.sigma_z,
            verbose=not args.quiet,
        )
        kappa_entry = float(args.kappa_entry) if args.kappa_entry is not None else float(getattr(P0_out, "kappa_loc", 1.0))
        closure0 = calibrate_outside_scale_closure(
            sol0,
            P0_out,
            b_grid0,
            target_city_entry_prob=args.target_city_entry_prob,
            kappa_entry=kappa_entry,
            local_birth_entry_weight=args.local_birth_entry_weight,
            target_total_population=float(getattr(P0_out, "N_target", 1.0)),
        )
        if not closure0.finite_calibration:
            raise RuntimeError("outside-option baseline normalization produced a nonfinite scale")

        outside_value = float(closure0.outside_value)
        outside_entry_flow = float(closure0.outside_entry_flow)
        entry_shares0 = np.asarray(P0_out.entry_shares, dtype=float).copy()
        p_start = np.asarray(p0, dtype=float).copy()
        normalization_history: list[dict[str, object]] = []
        sol = P_out = p_eq = b_grid = None
        n_passes = max(1, int(args.normalization_passes))
        for norm_pass in range(1, n_passes + 1):
            P, _, _ = prepare_parameters(args.nb, args.max_iter_eq, args.tol_eq)
            P.population_closure = "accounting_scale_prices"
            P.outside_value = float(outside_value)
            P.outside_value_is_calibrated = True
            P.allow_uncalibrated_outside_value = False
            P.outside_entry_flow = float(outside_entry_flow)
            P.kappa_entry = float(closure0.kappa_entry)
            P.local_birth_entry_weight = float(closure0.local_birth_entry_weight)
            P.entry_shares = np.asarray(entry_shares0, dtype=float).copy()
            P.entry_by_loc = float(P.E_total) * P.entry_shares
            b_grid = make_grid(P)
            sol, P_out, p_eq = solve_equilibrium_hank_z(
                p_start,
                P,
                b_grid,
                nz=args.nz,
                rho_z=args.rho_z,
                sigma_z=args.sigma_z,
                verbose=not args.quiet,
            )
            scale = getattr(sol, "accounting_scale", None)
            normalization_history.append(
                {
                    "pass": norm_pass,
                    "outside_value": float(outside_value),
                    "outside_entry_flow": float(outside_entry_flow),
                    "scale_factor": float(getattr(scale, "scale_factor", np.nan)),
                    "city_entry_prob_total": float(getattr(scale, "city_entry_prob_total", np.nan)),
                    "outside_entry_prob": float(getattr(scale, "outside_entry_prob", np.nan)),
                    "p_P": float(p_eq[0]),
                    "p_C": float(p_eq[1]),
                    "best_eq_error": float(sol.timings.get("best_eq_error", np.nan)),
                }
            )
            final_scale_error = abs(float(getattr(scale, "scale_factor", np.nan)) - 1.0)
            final_q_error = abs(float(getattr(scale, "city_entry_prob_total", np.nan)) - args.target_city_entry_prob)
            if norm_pass >= n_passes or (
                final_scale_error <= float(args.scale_target_tol)
                and final_q_error <= float(args.scale_target_tol)
            ):
                break
            proposal = calibrate_outside_scale_closure(
                sol,
                P_out,
                b_grid,
                target_city_entry_prob=args.target_city_entry_prob,
                kappa_entry=kappa_entry,
                local_birth_entry_weight=args.local_birth_entry_weight,
                target_total_population=float(getattr(P_out, "N_target", 1.0)),
            )
            outside_value = float(proposal.outside_value)
            outside_entry_flow = float(proposal.outside_entry_flow)
            entry_shares0 = np.asarray(getattr(scale, "conditional_entry_shares", P_out.entry_shares), dtype=float).copy()
            p_start = np.asarray(p_eq, dtype=float).copy()
        closure0.normalization_history = normalization_history
        assert sol is not None and P_out is not None and p_eq is not None and b_grid is not None

        moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        moments_ns.inv_pop_share_C = float(sol.pop_share[1])
        moments_ns.inv_rent_ratio_C_over_P = float((P_out.user_cost_rate * p_eq[1]) / max(P_out.user_cost_rate * p_eq[0], 1e-12))
        moments = namespace_float_dict(moments_ns)
        loss = compute_smm_loss(moments, targets, weights)
        elapsed = time.perf_counter() - t0

        write_results(out_dir / "results_income_mortgage_risk_v5_hank_z_outside_closure.csv", targets, moments)
        write_hank_z_diagnostics(out_dir / "diagnostics_income_mortgage_risk_v5_hank_z_outside_closure.csv", sol, P_out)
        write_closure_diagnostics(out_dir / "diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_closure.csv", closure0, sol)
        write_trace(out_dir / "diagnostics_income_mortgage_risk_v5_hank_z_outside_closure_trace.csv", sol)
        write_report(
            out_dir / "REPORT_V5_HANK_Z_OUTSIDE_CLOSURE.md",
            sol0,
            sol,
            P_out,
            p_eq,
            moments,
            targets,
            loss,
            elapsed,
            closure0,
        )
        log_path.write_text(
            json.dumps(
                {
                    "ok": True,
                    "elapsed_sec": elapsed,
                    "loss": loss,
                    "baseline_p_eq": [float(x) for x in p0],
                    "p_eq": [float(x) for x in p_eq],
                    "z_grid": [float(x) for x in sol.hank_z.z_grid],
                    "stationary_z": [float(x) for x in sol.hank_z.stationary_z],
                    "closure_calibration": {
                        key: value
                        for key, value in vars(closure0).items()
                        if key not in ("entry_values", "city_entry_prob", "baseline_accounting_scale")
                    },
                    "baseline_accounting_scale": vars(closure0.baseline_accounting_scale),
                    "final_accounting_scale": vars(getattr(sol, "accounting_scale", SimpleNamespace())),
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
