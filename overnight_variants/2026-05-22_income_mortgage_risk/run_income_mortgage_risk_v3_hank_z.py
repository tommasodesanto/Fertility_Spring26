#!/usr/bin/env python3
"""V3 Branch 1 prototype: structural HANK-style earnings risk.

This script corrects the earlier scenario-mixture pass by putting the
idiosyncratic earnings state ``z`` directly in the copied household Bellman
arrays, policy arrays, and forward distribution. It deliberately leaves the
minimal mortgage-account state for the next layer because structural ``mu``
requires tenure-dependent purchase, amortization, sale, and default transitions.
"""

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
from dt_cp_model.solver import make_grid, solve_partial_equilibrium_hank_z

from run_income_mortgage_risk_v2_scenarios import (
    BENCHMARK_MODEL,
    BENCHMARK_P,
    MOMENT_ORDER,
    build_base_parameters,
    namespace_float_dict,
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
                    "notes": "V3 structural z-state PE solve; mortgage account not structural",
                }
            )


def write_diagnostics(path: Path, sol: SimpleNamespace, P: SimpleNamespace) -> None:
    hz = sol.hank_z
    diag = hz.diagnostics
    rows: list[dict[str, object]] = []
    state_shape = diag.state_shape
    for key, value in state_shape.items():
        rows.append({"diagnostic": "state_size", "name": key, "value": value})
    for iz, z in enumerate(hz.z_grid):
        rows.append(
            {
                "diagnostic": "z_process",
                "z": float(z),
                "stationary_mass": float(hz.stationary_z[iz]),
                "transition_row": json.dumps([float(x) for x in hz.Pi_z[iz, :]]),
            }
        )
    for age_idx in range(P.J):
        age = P.age_start + age_idx * P.da
        for iz, z in enumerate(hz.z_grid):
            rows.append(
                {
                    "diagnostic": "z_by_age",
                    "age": age,
                    "z": float(z),
                    "mass_share": float(diag.z_by_age[age_idx, iz]),
                    "mean_liquid_wealth": float(diag.mean_liquid_wealth_by_z_age[age_idx, iz]),
                    "mean_income": float(diag.mean_income_by_z_age[age_idx, iz]),
                    "own_rate": float(diag.ownership_by_z_age[age_idx, iz]),
                }
            )
    for rec in diag.owner_by_z:
        rows.append({"diagnostic": "owner_by_z", **rec})
    for rec in diag.fertility_by_z:
        rows.append({"diagnostic": "fertility_by_z", **rec})
    for rec in diag.fertility_by_z_wealth_bin:
        rows.append({"diagnostic": "fertility_by_z_wealth_bin", **rec})
    rows.append(
        {
            "diagnostic": "mortgage_account_status",
            "name": "mu_structural",
            "value": 0,
            "notes": "V3 implements z as a true HANK state; mu remains the next structural layer.",
        }
    )
    fieldnames = [
        "diagnostic",
        "name",
        "value",
        "notes",
        "age",
        "z",
        "stationary_mass",
        "transition_row",
        "mass_share",
        "mean_liquid_wealth",
        "mean_income",
        "own_rate",
        "fertility_choice",
        "wealth_bin",
        "expected_completed_fertility_choice",
        "mass",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    sol: SimpleNamespace,
    P: SimpleNamespace,
    moments: dict[str, float],
    targets: dict[str, float],
    loss: float,
    elapsed: float,
) -> None:
    hz = sol.hank_z
    diag = hz.diagnostics
    owner_vals = [float(x["own_rate"]) for x in diag.owner_by_z]
    fert_vals = [float(x["fertility_choice"]) for x in diag.fertility_by_z]
    own_spread = max(owner_vals) - min(owner_vals) if owner_vals else np.nan
    fert_spread = max(fert_vals) - min(fert_vals) if fert_vals else np.nan
    verdict = "yellow"
    if not np.isfinite(loss) or own_spread < 1e-3:
        verdict = "red"
    table = "\n".join(
        f"| `{name}` | {targets.get(name, np.nan):.3f} | {BENCHMARK_MODEL.get(name, np.nan):.3f} | {moments.get(name, np.nan):.3f} |"
        for name in MOMENT_ORDER
    )
    text = f"""# Income Risk V3 HANK-z Structural Report

Verdict: **{verdict}**

## What Changed

This pass implements the classic HANK object the earlier passes were missing:
a finite idiosyncratic earnings state \(z\) with a Markov transition matrix
\(\Pi_z\). The copied solver now carries value functions, policies, fertility
probabilities, location probabilities, tenure choices, and the forward
distribution over \((b,d,i,a,n,s,z)\). Continuation values average over
\(\Pi_z\) after child aging. Working-age income is
\(y_{{ia}}(z)=(1-\\tau_{{pay}})w_i e_a\exp(z)\); retirement uses the copied
common pension mapping.

The mortgage account \(\mu\) is not structural in V3. That is intentional for
this correction: the HANK earnings state is now real, while defaultable account
dynamics remain the next layer because \(\mu'\) must depend on purchase, sale,
amortization, and default timing.

## Status

- solve: fixed-price partial equilibrium at copied benchmark prices
- \(z\) states: `{len(hz.z_grid)}`
- \(z\) grid: `{[float(x) for x in hz.z_grid]}`
- stationary \(z\) distribution: `{[float(x) for x in hz.stationary_z]}`
- \(b\) states: `{P.Nb}`
- tenure states: `{1 + P.n_house}`
- locations: `{P.I}`
- ages: `{P.J}`
- fertility states: `{P.n_parity}`
- child-age states: `{P.n_child_states}`
- elapsed seconds: `{elapsed:.2f}`
- runtime category: `{sol.timings.get("runtime_category", "moderate")}`
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

| Moment | Target | Benchmark | V3 |
|---|---:|---:|---:|
{table}

## Read

This is the correct computational direction for Branch 1: idiosyncratic labor
income risk is now an actual discrete state, not a post-solve label or a set of
separate scenario economies. The test should still be read as yellow rather
than green because it is fixed-price PE and does not yet include the structural
mortgage-account state \(\mu\). The next Branch 1 pass should add \(\mu\) as a
second finite state with tenure-dependent account transitions, then rerun this
same target table.
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nb", type=int, default=30)
    parser.add_argument("--nz", type=int, default=3, choices=[3, 5])
    parser.add_argument("--rho-z", type=float, default=0.82)
    parser.add_argument("--sigma-z", type=float, default=0.28)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    log_path = out_dir / "income_mortgage_risk_v3_hank_z.log"
    t0 = time.perf_counter()
    try:
        P_override, targets, weights = build_base_parameters(args.nb, force_full=False)
        P = setup_parameters()
        P = apply_overrides(P, P_override)
        P = finalize_location_choice_spec(P)
        P.Nb = int(args.nb)
        b_grid = make_grid(P)
        sol, P_out, p_eq = solve_partial_equilibrium_hank_z(
            BENCHMARK_P,
            np.ones(P.I),
            np.array([0.55, 0.45]),
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
        write_results(out_dir / "results_income_mortgage_risk_v3_hank_z.csv", targets, moments)
        write_diagnostics(out_dir / "diagnostics_income_mortgage_risk_v3_hank_z.csv", sol, P_out)
        write_report(out_dir / "REPORT_V3_HANK_Z.md", sol, P_out, moments, targets, loss, elapsed)
        log_path.write_text(
            json.dumps(
                {
                    "ok": True,
                    "elapsed_sec": elapsed,
                    "loss": loss,
                    "moments": moments,
                    "z_grid": [float(x) for x in sol.hank_z.z_grid],
                    "stationary_z": [float(x) for x in sol.hank_z.stationary_z],
                    "timings": sol.timings,
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
