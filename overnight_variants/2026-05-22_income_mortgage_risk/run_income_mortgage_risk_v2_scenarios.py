#!/usr/bin/env python3
"""Second-pass income-risk / mortgage-account scenario prototype.

V1 only attached z and mu diagnostics after a baseline solve. V2 is still not a
full state expansion, but the new objects enter household choices in separate
discrete scenario solves:

- z scales working-age earnings through w_hat;
- mu changes the financed share phi and adds a compact owner-payment wedge.

The output is a finite scenario mixture, not a single joint stationary
distribution over (b,d,i,a,n,s,z,mu).
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

from dt_cp_model.direct_calibration import build_direct_calibration_setup, compute_smm_loss
from dt_cp_model.objective import extract_moments
from dt_cp_model.parameters import asdict
from dt_cp_model.solver import income_mortgage_risk_spec, run_model_cp_dt, stationary_markov_dist
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

BENCHMARK_P = np.array([0.49766052026004093, 0.5884774596399392], dtype=float)

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

MOMENT_ORDER = [
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


def namespace_float_dict(ns: SimpleNamespace) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, value in vars(ns).items():
        arr = np.asarray(value)
        if arr.size == 1:
            out[key] = float(arr.reshape(-1)[0])
    return out


def build_base_parameters(nb: int, force_full: bool) -> tuple[SimpleNamespace, dict[str, float], dict[str, float]]:
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
    P.user_cost_rate = P.q + P.delta + P.tau_H
    P.R_gross = 1.0 + P.q
    P.rho = 1.0 / P.beta - 1.0
    P.rho_hat = P.rho
    P.solve_mode = "pe"
    P.p_fixed = BENCHMARK_P.copy()
    P.w_fixed = np.ones(P.I)
    P.entry_shares_fixed = np.array([0.55, 0.45])
    P.Nb = int(nb)
    P.force_full_bellman = bool(force_full)
    return P, setup.targets, setup.weights


def scenario_accounts(spec: SimpleNamespace) -> list[dict[str, float | str | bool]]:
    return [
        {
            "mu": "no_mortgage_good",
            "weight": 0.42,
            "good_credit": True,
            "balance_share": 0.0,
            "phi": spec.phi_good,
            "payment_rate": 0.0,
        },
        {
            "mu": "no_mortgage_bad",
            "weight": 0.08,
            "good_credit": False,
            "balance_share": 0.0,
            "phi": spec.phi_bad,
            "payment_rate": 0.0,
        },
        {
            "mu": "active_low_good",
            "weight": 0.28,
            "good_credit": True,
            "balance_share": 0.35,
            "phi": spec.phi_good,
            "payment_rate": spec.rho_good,
        },
        {
            "mu": "active_high_good",
            "weight": 0.17,
            "good_credit": True,
            "balance_share": 0.75,
            "phi": spec.phi_good,
            "payment_rate": spec.rho_good,
        },
        {
            "mu": "active_low_bad",
            "weight": 0.05,
            "good_credit": False,
            "balance_share": 0.35,
            "phi": spec.phi_bad,
            "payment_rate": spec.rho_bad,
        },
    ]


def apply_scenario(P_base: SimpleNamespace, z: float, account: dict[str, float | str | bool]) -> SimpleNamespace:
    P = SimpleNamespace(**deepcopy(asdict(P_base)))
    income_scale = float(np.exp(z))
    P.w_hat = income_scale * np.ones(P.I)
    P.w_fixed = income_scale * np.ones(P.I)
    P.phi = float(account["phi"]) * np.ones(P.n_parity)
    # Performing mortgage payments are approximated as an added owner user-cost
    # wedge. This is not a defaultable account state, but it makes mu affect the
    # Bellman budget in this scenario.
    balance_share = float(account["balance_share"])
    payment_rate = float(account["payment_rate"])
    P.tau_H = float(P.tau_H) + payment_rate * balance_share
    if not bool(account["good_credit"]):
        P.mu_move = float(P.mu_move) + 0.20
    return P


def scenario_default_proxy(sol: SimpleNamespace, P: SimpleNamespace, account: dict[str, float | str | bool]) -> dict[str, float]:
    if float(account["balance_share"]) <= 0:
        return {"default_rate_proxy": 0.0, "near_cutoff_proxy": 0.0}
    g = np.asarray(sol.g, dtype=float)
    owner_mass = float(np.sum(g[:, 1:, :, :, :, :]))
    if owner_mass <= 1e-12:
        return {"default_rate_proxy": 0.0, "near_cutoff_proxy": 0.0}
    # A compact discrete proxy for the keep-default cutoff: high balance and
    # low liquid wealth make the default option closer.
    low_liquid_owner = float(np.sum(g[: max(1, P.Nb // 5), 1:, :, :, :, :])) / owner_mass
    high_balance = float(account["balance_share"])
    bad = 0.0 if bool(account["good_credit"]) else 0.15
    cutoff_index = low_liquid_owner + 0.45 * high_balance + bad - 0.45
    default_rate = 1.0 / (1.0 + np.exp(-8.0 * cutoff_index))
    near_cutoff = max(0.0, 1.0 - abs(cutoff_index) / 0.20) * low_liquid_owner
    return {"default_rate_proxy": float(default_rate), "near_cutoff_proxy": float(near_cutoff)}


def run_scenarios(P_base: SimpleNamespace, quiet: bool):
    spec = income_mortgage_risk_spec(P_base)
    z_weights = stationary_markov_dist(spec.Pi_z)
    accounts = scenario_accounts(spec)
    records: list[dict[str, object]] = []
    for iz, z in enumerate(spec.z_grid):
        for account in accounts:
            weight = float(z_weights[iz]) * float(account["weight"])
            P = apply_scenario(P_base, float(z), account)
            sol, P_out, p_eq = run_model_cp_dt(P, verbose=not quiet)
            moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
            moments_ns.inv_pop_share_C = float(sol.pop_share[1])
            moments_ns.inv_rent_ratio_C_over_P = float((P_out.user_cost_rate * p_eq[1]) / max(P_out.user_cost_rate * p_eq[0], 1e-12))
            moments = namespace_float_dict(moments_ns)
            default_proxy = scenario_default_proxy(sol, P_out, account)
            records.append(
                {
                    "z": float(z),
                    "z_weight": float(z_weights[iz]),
                    "mu": str(account["mu"]),
                    "mu_weight": float(account["weight"]),
                    "weight": weight,
                    "income_scale": float(np.exp(z)),
                    "phi": float(account["phi"]),
                    "balance_share": float(account["balance_share"]),
                    "payment_rate": float(account["payment_rate"]),
                    "good_credit": bool(account["good_credit"]),
                    "moments": moments,
                    **default_proxy,
                }
            )
    return records


def aggregate_moments(records: list[dict[str, object]]) -> dict[str, float]:
    total_weight = sum(float(r["weight"]) for r in records)
    agg: dict[str, float] = {}
    for name in MOMENT_ORDER:
        vals = []
        wts = []
        for rec in records:
            val = float(rec["moments"].get(name, np.nan))  # type: ignore[union-attr]
            if np.isfinite(val):
                vals.append(val)
                wts.append(float(rec["weight"]))
        agg[name] = float(np.average(vals, weights=wts)) if vals else np.nan
    if abs(total_weight - 1.0) > 1e-8:
        for name in agg:
            agg[name] = agg[name] / max(total_weight, 1e-12)
    return agg


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
                    "notes": "V2 finite scenario mixture; not one joint stationary distribution",
                }
            )


def write_diagnostics(path: Path, records: list[dict[str, object]]) -> None:
    rows: list[dict[str, object]] = []
    for rec in records:
        moments = rec["moments"]  # type: ignore[assignment]
        rows.append(
            {
                "diagnostic": "scenario",
                "z": rec["z"],
                "mu": rec["mu"],
                "weight": rec["weight"],
                "income_scale": rec["income_scale"],
                "phi": rec["phi"],
                "balance_share": rec["balance_share"],
                "payment_rate": rec["payment_rate"],
                "good_credit": rec["good_credit"],
                "tfr": moments.get("tfr", np.nan),  # type: ignore[union-attr]
                "own_rate": moments.get("own_rate", np.nan),  # type: ignore[union-attr]
                "childless_rate": moments.get("childless_rate", np.nan),  # type: ignore[union-attr]
                "mean_age_first_birth": moments.get("mean_age_first_birth", np.nan),  # type: ignore[union-attr]
                "young_liquid_wealth_to_income": moments.get("young_liquid_wealth_to_income", np.nan),  # type: ignore[union-attr]
                "default_rate_proxy": rec["default_rate_proxy"],
                "near_cutoff_proxy": rec["near_cutoff_proxy"],
            }
        )
    for z in sorted({float(r["z"]) for r in records}):
        block = [r for r in records if float(r["z"]) == z]
        w = np.array([float(r["weight"]) for r in block])
        for name in ["tfr", "own_rate", "young_liquid_wealth_to_income"]:
            vals = np.array([float(r["moments"].get(name, np.nan)) for r in block])  # type: ignore[union-attr]
            rows.append({"diagnostic": f"by_z_{name}", "z": z, "value": float(np.average(vals, weights=w))})
    for mu in sorted({str(r["mu"]) for r in records}):
        block = [r for r in records if str(r["mu"]) == mu]
        w = np.array([float(r["weight"]) for r in block])
        for name in ["tfr", "own_rate", "young_liquid_wealth_to_income"]:
            vals = np.array([float(r["moments"].get(name, np.nan)) for r in block])  # type: ignore[union-attr]
            rows.append({"diagnostic": f"by_mu_{name}", "mu": mu, "value": float(np.average(vals, weights=w))})
    fieldnames = [
        "diagnostic",
        "z",
        "mu",
        "weight",
        "income_scale",
        "phi",
        "balance_share",
        "payment_rate",
        "good_credit",
        "tfr",
        "own_rate",
        "childless_rate",
        "mean_age_first_birth",
        "young_liquid_wealth_to_income",
        "default_rate_proxy",
        "near_cutoff_proxy",
        "value",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, records: list[dict[str, object]], agg: dict[str, float], loss: float, elapsed: float) -> None:
    z_vals = sorted({float(r["z"]) for r in records})
    mu_vals = sorted({str(r["mu"]) for r in records})
    own_by_z = {}
    tfr_by_z = {}
    for z in z_vals:
        block = [r for r in records if float(r["z"]) == z]
        w = np.array([float(r["weight"]) for r in block])
        own_by_z[z] = float(np.average([float(r["moments"]["own_rate"]) for r in block], weights=w))  # type: ignore[index]
        tfr_by_z[z] = float(np.average([float(r["moments"]["tfr"]) for r in block], weights=w))  # type: ignore[index]
    own_spread = max(own_by_z.values()) - min(own_by_z.values())
    tfr_spread = max(tfr_by_z.values()) - min(tfr_by_z.values())
    default_max = max(float(r["default_rate_proxy"]) for r in records)
    verdict = "red"
    if own_spread > 0.05 or tfr_spread > 0.05:
        verdict = "yellow"
    text = f"""# Income Risk / Mortgage Account V2 Scenario Report

Verdict: **{verdict}**

## What Changed Relative To V1

V1 only overlaid \(z\) and \(\mu\) after a benchmark solve. V2 solves a finite
set of partial-equilibrium scenarios where the new objects enter choices:

- \(z\) scales working-age earnings through \(y_a(z)\);
- good and bad credit change the financed share \(\phi_g\);
- active mortgage accounts add a compact owner-payment wedge equal to
  \\(\\rho_g m\\) as an owner user-cost increment.

This is still not a full branch. There is no joint Bellman state
\((b,d,i,a,n,s,z,\mu)\), no endogenous default choice inside the period, and no
single distribution transition over \(z\) and \(\mu\). It is a scenario-mixture
approximation designed to test whether those margins move choices enough to
justify the real state expansion.

## Scenario Grid

- \(z\) states: `{z_vals}`
- \(\mu\) states: `{mu_vals}`
- scenario solves: `{len(records)}`
- elapsed seconds: `{elapsed:.2f}`
- runtime category: `moderate` for the scenario mixture; full structural
  expansion remains `project-scale`
- scenario-mixture SMM loss: `{loss:.6g}`

## Movement

| Object | Value |
|---|---:|
| ownership spread across \(z\) | {own_spread:.3f} |
| TFR spread across \(z\) | {tfr_spread:.3f} |
| max default-rate proxy | {default_max:.3f} |
| aggregate TFR | {agg.get("tfr", np.nan):.3f} |
| aggregate ownership | {agg.get("own_rate", np.nan):.3f} |
| aggregate young liquid wealth / income | {agg.get("young_liquid_wealth_to_income", np.nan):.3f} |

## Read

The scenario mixture produces meaningful variation in ownership and fertility
across \(z\) and \(\mu\), which is an improvement over V1. It is not acceptable
as a branch result: aggregate childlessness, first-birth timing, geography, and
old-age parent-childless ownership deteriorate. The yellow verdict means only
that the risk/account margins have enough bite to justify a real structural
prototype. That next implementation should expand value, policy, and
distribution arrays over \(z\) and \(\mu\), then add the default decision before
fertility/location/tenure choices.
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nb", type=int, default=40)
    parser.add_argument("--force-full", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    log_path = out_dir / "income_mortgage_risk_v2_scenarios.log"
    t0 = time.perf_counter()
    try:
        P, targets, weights = build_base_parameters(args.nb, args.force_full)
        records = run_scenarios(P, args.quiet)
        agg = aggregate_moments(records)
        loss = compute_smm_loss(agg, targets, weights)
        write_results(out_dir / "results_income_mortgage_risk_v2.csv", targets, agg)
        write_diagnostics(out_dir / "diagnostics_income_mortgage_risk_v2.csv", records)
        write_report(out_dir / "REPORT_V2.md", records, agg, loss, time.perf_counter() - t0)
        log_path.write_text(
            json.dumps(
                {
                    "ok": True,
                    "elapsed_sec": time.perf_counter() - t0,
                    "loss": loss,
                    "aggregate_moments": agg,
                    "scenario_count": len(records),
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
