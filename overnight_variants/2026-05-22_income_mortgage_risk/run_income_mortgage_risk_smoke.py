#!/usr/bin/env python3
"""Smoke test for the isolated income-risk / mortgage-account branch.

This driver deliberately runs one coarse equilibrium at the current direct
benchmark theta and then attaches the finite z and mu diagnostic state grids
implemented in the copied solver. It does not run a global calibration.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from dt_cp_model.direct_calibration import (
    build_direct_calibration_setup,
    compute_smm_loss,
)
from dt_cp_model.objective import extract_moments
from dt_cp_model.parameters import asdict
from dt_cp_model.solver import (
    attach_income_mortgage_risk_diagnostics,
    make_grid,
    run_model_cp_dt,
)
from dt_cp_model.theta import apply_theta


BEST_DIRECT_THETA = np.array(
    [
        0.9805235878335584,
        0.31883369129402983,
        0.09201816432650124,
        0.24595815946252964,
        0.7067802199754748,
        0.09951054009758557,
        3.0417838844424634,
        1.0738923846894777,
        2.6815562612708366,
        8.644845291269531,
        0.41101999602420053,
        0.6219978476345638,
        3.0532061053532615,
        0.11127541027064872,
        0.0719341870623939,
    ],
    dtype=float,
)


BENCHMARK_MODEL = {
    "tfr": 1.898,
    "childless_rate": 0.145,
    "mean_age_first_birth": 33.535,
    "tfr_gradient": 0.119,
    "own_rate": 0.643,
    "own_gradient": 0.139,
    "own_family_gap": 0.114,
    "prime_childless_renter_median_rooms": 6.365,
    "prime_childless_owner_median_rooms": 6.800,
    "housing_increment_0to1": 0.441,
    "housing_increment_1to2": 0.192,
    "young_liquid_wealth_to_income": 0.527,
    "center_share_nonparents": 0.405,
    "center_share_newparents": 0.382,
    "migration_rate": 0.035,
    "old_age_own_rate": 0.947,
    "old_age_parent_childless_gap": 0.062,
    "inv_pop_share_C": 0.441,
    "inv_rent_ratio_C_over_P": 1.182,
}


def namespace_float_dict(ns: SimpleNamespace) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, value in vars(ns).items():
        arr = np.asarray(value)
        if arr.size == 1:
            out[key] = float(arr.reshape(-1)[0])
    return out


def build_smoke_parameters(max_iter_eq: int, force_full: bool) -> tuple[SimpleNamespace, dict[str, float], dict[str, float]]:
    setup = build_direct_calibration_setup(
        setup_mode="benchmark",
        geo_weight=100.0,
        population_closure="renewal_valve_calibrated",
        renewal_retention=1.0,
    )
    theta_dict = {name: float(value) for name, value in zip(setup.names, BEST_DIRECT_THETA)}
    structural_names = [name for name in setup.names if name not in ("E_C", "r_bar_C")]
    structural_theta = [theta_dict[name] for name in structural_names]

    P = SimpleNamespace(**deepcopy(asdict(setup.P_base)))
    P = apply_theta(P, structural_theta, structural_names)
    P.E_loc = np.array([float(P.E_loc[0]), theta_dict["E_C"]])
    P.r_bar = np.array([float(P.r_bar[0]), theta_dict["r_bar_C"]])
    P.population_closure = "renewal_valve_calibrated"
    P.renewal_calibrate_outside_flow = True
    P.renewal_target_total_population = 1.0
    P.renewal_retention = 1.0
    P.max_iter_eq = int(max_iter_eq)
    P.force_full_bellman = bool(force_full)
    return P, setup.targets, setup.weights


def write_results(path: Path, targets: dict[str, float], model: dict[str, float]) -> None:
    order = [
        "tfr",
        "childless_rate",
        "mean_age_first_birth",
        "tfr_gradient",
        "own_rate",
        "own_gradient",
        "own_family_gap",
        "prime_childless_renter_median_rooms",
        "prime_childless_owner_median_rooms",
        "housing_increment_0to1",
        "housing_increment_1to2",
        "young_liquid_wealth_to_income",
        "center_share_nonparents",
        "center_share_newparents",
        "migration_rate",
        "old_age_own_rate",
        "old_age_parent_childless_gap",
        "inv_pop_share_C",
        "inv_rent_ratio_C_over_P",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["moment", "target", "model", "benchmark_model", "notes"],
        )
        writer.writeheader()
        for name in order:
            writer.writerow(
                {
                    "moment": name,
                    "target": targets.get(name, ""),
                    "model": model.get(name, np.nan),
                    "benchmark_model": BENCHMARK_MODEL.get(name, ""),
                    "notes": "target from copied direct setup; benchmark model from CALIBRATION_STATUS.md",
                }
            )


def write_diagnostics(path: Path, diag: SimpleNamespace, P: SimpleNamespace) -> None:
    rows: list[dict[str, object]] = []
    for key, value in diag.state_shape.items():
        rows.append({"diagnostic": "state_size", "dimension": key, "value": value})

    age_indices = sorted(set([0, 7, 17, 27, 47, P.J - 1]))
    for j in age_indices:
        age = P.age_start + j
        for iz, z in enumerate(diag.spec.z_grid):
            rows.append(
                {
                    "diagnostic": "z_distribution_by_age",
                    "age": age,
                    "z": z,
                    "value": diag.z_by_age[j, iz],
                }
            )
            rows.append(
                {
                    "diagnostic": "mean_liquid_wealth_by_z_age",
                    "age": age,
                    "z": z,
                    "value": diag.mean_liquid_wealth_by_z_age[j, iz],
                }
            )
            rows.append(
                {
                    "diagnostic": "ownership_by_z_age",
                    "age": age,
                    "z": z,
                    "value": diag.ownership_by_z_age[j, iz],
                }
            )
        rows.append(
            {
                "diagnostic": "bad_credit_share_by_age",
                "age": age,
                "value": diag.bad_credit_share_by_age[j],
            }
        )
        for imu, label in enumerate(diag.spec.mu_labels):
            rows.append(
                {
                    "diagnostic": "mortgage_account_mass_by_age",
                    "age": age,
                    "mu": label,
                    "value": diag.mortgage_account_mass_by_age[j, imu],
                }
            )

    for item in diag.fertility_by_z_wealth_bin:
        rows.append(
            {
                "diagnostic": "fertility_by_z_liquid_wealth_bin",
                "z": item["z"],
                "wealth_bin": item["wealth_bin"],
                "value": item["expected_completed_fertility_choice"],
                "mass": item["mass"],
            }
        )

    avg_default = np.nanmean(diag.default_rate_by_age_loc_fert_z, axis=0)
    avg_cutoff = np.nanmean(diag.near_default_cutoff_mass, axis=0)
    for i in range(avg_default.shape[0]):
        for nn in range(avg_default.shape[1]):
            for iz, z in enumerate(diag.spec.z_grid):
                rows.append(
                    {
                        "diagnostic": "default_rate_by_location_fertility_z",
                        "location": i,
                        "parity": nn,
                        "z": z,
                        "value": avg_default[i, nn, iz],
                    }
                )
                rows.append(
                    {
                        "diagnostic": "mass_near_default_cutoff",
                        "location": i,
                        "parity": nn,
                        "z": z,
                        "value": avg_cutoff[i, nn, iz],
                    }
                )

    for credit, share in diag.purchase_origination_share_by_credit.items():
        rows.append(
            {
                "diagnostic": "purchase_origination_share_by_credit_status",
                "credit_status": credit,
                "value": share,
            }
        )

    fieldnames = [
        "diagnostic",
        "dimension",
        "age",
        "location",
        "parity",
        "z",
        "mu",
        "wealth_bin",
        "credit_status",
        "value",
        "mass",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    result: dict[str, float],
    loss: float,
    timings: dict[str, object],
    diag: SimpleNamespace,
    elapsed: float,
) -> None:
    accepted = bool(timings.get("accepted", False))
    default_max = float(np.nanmax(diag.default_rate_by_age_loc_fert_z))
    mu_mass = np.asarray(diag.mortgage_account_mass_by_age)
    mu_positive = int(np.sum(np.nanmean(mu_mass, axis=0) > 1e-6))
    z_nonzero = int(np.sum(np.nanmean(diag.z_by_age, axis=0) > 1e-6))
    verdict = "red"
    if accepted and z_nonzero > 1 and mu_positive > 1:
        verdict = "yellow"
    if accepted and bool(diag.structural_bellman_augmented):
        verdict = "green"

    text = f"""# Income Risk / Mortgage Account Smoke Report

Verdict: **{verdict}**

## What Ran

The copied benchmark solver ran once at the current direct benchmark theta and
then constructed finite diagnostic grids for \(z\) and \(\mu\). The grids are
non-degenerate, but the Bellman recursion and forward distribution in this
overnight prototype do **not** yet optimize over those two states. This is a
diagnostic state augmentation, not an accepted structural HANK/mortgage branch.

## Solve Status

- accepted equilibrium from copied baseline solver: `{accepted}`
- convergence reason: `{timings.get("convergence_reason", "unknown")}`
- best equilibrium error: `{timings.get("best_eq_error", np.nan)}`
- elapsed seconds: `{elapsed:.2f}`
- runtime category: `moderate` for one copied benchmark solve; full structural
  \(z,\mu\) state expansion remains `project-scale`
- SMM loss on this smoke solve: `{loss:.6g}`

## State Sizes

| Dimension | Count |
|---|---:|
| liquid wealth b | {diag.state_shape["Nb"]} |
| tenure d | {diag.state_shape["tenure_states"]} |
| locations | {diag.state_shape["locations"]} |
| ages | {diag.state_shape["ages"]} |
| completed-fertility states | {diag.state_shape["fertility_states"]} |
| child-age states | {diag.state_shape["child_age_states"]} |
| earnings z states | {diag.state_shape["z_states"]} |
| mortgage-account mu states | {diag.state_shape["mu_states"]} |

## Branch Diagnostics

- non-degenerate z states used in diagnostic distribution: `{z_nonzero}`
- non-degenerate mu states used in diagnostic distribution: `{mu_positive}`
- maximum diagnostic default hazard: `{default_max:.6g}`
- average bad-credit share at age 35: `{diag.bad_credit_share_by_age[min(17, len(diag.bad_credit_share_by_age)-1)]:.6g}`
- origination shares by credit status: `{json.dumps(diag.purchase_origination_share_by_credit)}`

## Interpretation

This fails the full Branch 1 acceptance criterion because the new states are
not yet inside the Bellman equation. It is still useful as a smoke test of the
state definitions, diagnostic accounting, and report pipeline. The next
implementation step is to expand value, policy, and distribution arrays from
`(b,d,i,a,n,s)` to `(b,d,i,a,n,s,z,mu)` and to feed \(y_a(z)\), mortgage
payments, amortization, default, and credit recovery directly into the copied
solver.
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-iter-eq", type=int, default=35)
    parser.add_argument("--force-full", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    log_path = out_dir / "income_mortgage_risk_smoke.log"
    t0 = time.perf_counter()
    try:
        P, targets, weights = build_smoke_parameters(args.max_iter_eq, args.force_full)
        sol, P_out, p_eq = run_model_cp_dt(P, verbose=not args.quiet)
        rent_ratio = float((P_out.user_cost_rate * p_eq[1]) / (P_out.user_cost_rate * p_eq[0]))
        moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        moments_ns.inv_pop_share_C = float(sol.pop_share[1])
        moments_ns.inv_rent_ratio_C_over_P = rent_ratio
        moments = namespace_float_dict(moments_ns)
        loss = compute_smm_loss(moments, targets, weights)
        b_grid = make_grid(P_out)
        diag = attach_income_mortgage_risk_diagnostics(sol, P_out, p_eq, b_grid)
        write_results(out_dir / "results_income_mortgage_risk.csv", targets, moments)
        write_diagnostics(out_dir / "diagnostics_income_mortgage_risk.csv", diag, P_out)
        write_report(out_dir / "REPORT.md", moments, loss, getattr(sol, "timings", {}), diag, time.perf_counter() - t0)
        log_path.write_text(
            json.dumps(
                {
                    "ok": True,
                    "elapsed_sec": time.perf_counter() - t0,
                    "loss": loss,
                    "timings": getattr(sol, "timings", {}),
                    "p_eq": np.asarray(p_eq).tolist(),
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
