#!/usr/bin/env python3
"""Second-pass developer supply smoke test.

V1 intentionally pinned ordinary prices and mapped nearly all owner rungs into
the middle type. That was too crude: the middle price rose until ownership
collapsed. V2 uses three unit types and calibrates type-specific supply around
the scalar-price baseline demand before moving prices.
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
from dt_cp_model.solver import (
    compute_developer_type_diagnostics,
    developer_inverse_supply_rent,
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
    P.use_full_kernel = False
    P.use_eval_kernel = False

    P.developer_supply = True
    P.developer_n_types = 3
    P.developer_type_names = np.array(["S", "M", "L"], dtype=object)
    P.developer_middle_min = 5.0
    P.developer_middle_max = 6.5
    return P, setup.targets, setup.weights


def run_with_prices(P: SimpleNamespace, p_iq: np.ndarray, quiet: bool):
    P_iter = SimpleNamespace(**deepcopy(asdict(P)))
    P_iter.developer_p_iq = np.asarray(p_iq, dtype=float)
    P_iter.developer_r_iq = P_iter.user_cost_rate * P_iter.developer_p_iq
    P_iter.p_fixed = P_iter.developer_p_iq[:, 0].copy()
    sol, P_out, p_eq = run_model_cp_dt(P_iter, verbose=not quiet)
    diag = compute_developer_type_diagnostics(sol, P_out, P_iter.developer_p_iq, P_iter.developer_r_iq)
    return sol, P_out, p_eq, diag


def initialize_supply_from_baseline(P: SimpleNamespace, baseline_diag: SimpleNamespace, middle_wedge: float, large_wedge: float):
    scalar_p = np.tile(BENCHMARK_P.reshape(P.I, 1), (1, 3))
    p_iq = scalar_p.copy()
    p_iq[:, 1] *= 1.0 + middle_wedge
    p_iq[:, 2] *= 1.0 + large_wedge
    r_iq = P.user_cost_rate * p_iq

    demand0 = np.maximum(np.asarray(baseline_diag.demand_iq, dtype=float), 1e-5)
    eta_iq = np.column_stack(
        [
            P.eta_supply,
            np.array([1.20, 0.70]),
            np.array([1.70, 0.90]),
        ]
    )
    tau_iq = np.column_stack(
        [
            np.zeros(P.I),
            0.04 * r_iq[:, 0],
            0.02 * r_iq[:, 0],
        ]
    )
    nu_iq = np.maximum((r_iq - tau_iq) / demand0 ** (1.0 / eta_iq), 1e-8)
    F_iq = np.column_stack([np.zeros(P.I), np.array([0.001, 0.002]), np.array([0.001, 0.002])])
    return p_iq, r_iq, tau_iq, nu_iq, eta_iq, F_iq


def run_type_loop(P: SimpleNamespace, iterations: int, quiet: bool, middle_wedge: float, large_wedge: float):
    baseline_p = np.tile(BENCHMARK_P.reshape(P.I, 1), (1, 3))
    sol0, P0, _, baseline_diag = run_with_prices(P, baseline_p, quiet)
    p_iq, r_iq, tau_iq, nu_iq, eta_iq, F_iq = initialize_supply_from_baseline(P, baseline_diag, middle_wedge, large_wedge)
    P.developer_tau_iq = tau_iq
    P.developer_nu_iq = nu_iq
    P.developer_eta_iq = eta_iq
    P.developer_F_iq = F_iq

    history: list[dict] = []
    sol = sol0
    P_out = P0
    p_eq = BENCHMARK_P.copy()
    diag = baseline_diag
    update_types = [1, 2]
    for it in range(1, iterations + 1):
        sol, P_out, p_eq, diag = run_with_prices(P, p_iq, quiet)
        target_r = developer_inverse_supply_rent(P_out, diag.demand_iq)
        old = r_iq.copy()
        for q in update_types:
            r_iq[:, q] = 0.75 * r_iq[:, q] + 0.25 * target_r[:, q]
        r_iq[:, 0] = P.user_cost_rate * BENCHMARK_P
        # Keep the smoke test in the local branch-neighborhood. If the loop
        # wants to leave this band, that is a failure of the partial test.
        r_iq[:, 1] = np.clip(r_iq[:, 1], 0.92 * r_iq[:, 0], 1.45 * r_iq[:, 0])
        r_iq[:, 2] = np.clip(r_iq[:, 2], 0.90 * r_iq[:, 0], 1.35 * r_iq[:, 0])
        p_iq = r_iq / P.user_cost_rate
        max_step = float(np.max(np.abs(r_iq[:, update_types] - old[:, update_types]) / np.maximum(old[:, update_types], 1e-12)))
        max_active_excess = float(
            np.max(np.abs(diag.excess_iq[:, update_types]) / np.maximum(np.abs(diag.demand_iq[:, update_types]), 1e-8))
        )
        history.append(
            {
                "iter": it,
                "p_iq": p_iq.tolist(),
                "r_iq": r_iq.tolist(),
                "demand_iq": diag.demand_iq.tolist(),
                "supply_iq": diag.supply_iq.tolist(),
                "max_rel_price_step_active": max_step,
                "max_active_excess_rel_demand": max_active_excess,
            }
        )
        if max_step < 0.015:
            break

    final_P = SimpleNamespace(**deepcopy(asdict(P_out)))
    final_P.developer_p_iq = p_iq.copy()
    final_P.developer_r_iq = r_iq.copy()
    final_P.developer_tau_iq = tau_iq
    final_P.developer_nu_iq = nu_iq
    final_P.developer_eta_iq = eta_iq
    final_P.developer_F_iq = F_iq
    final_diag = compute_developer_type_diagnostics(sol, final_P, p_iq, r_iq)
    return sol, final_P, p_eq, final_diag, history, baseline_diag


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
                    "notes": "V2 three-type partial PE loop; S pinned, M/L updated",
                }
            )


def weighted_room_bins(sol: SimpleNamespace, P: SimpleNamespace) -> list[dict[str, object]]:
    g = sol.g
    hR = sol.hR_pol
    rows: list[dict[str, object]] = []
    a25s = max(0, round(25 - P.age_start))
    a45e = min(P.J - 1, round(45 - P.age_start))
    renter_bins = [("<=4", -np.inf, 4.0), ("5", 4.0, 5.0), ("6", 5.0, 6.0), ("7-8", 6.0, 8.0)]
    owner_bins = [("5", 4.0, 5.5), ("6", 5.5, 6.5), ("7-8", 6.5, 8.5), ("9-10", 8.5, 10.5), ("11+", 10.5, np.inf)]
    renter_total = 0.0
    owner_total = 0.0
    renter_mass = {name: 0.0 for name, _, _ in renter_bins}
    owner_mass = {name: 0.0 for name, _, _ in owner_bins}
    renter_cap = 0.0
    lower_middle_owner = 0.0
    for j in range(a25s, a45e + 1):
        for i in range(P.I):
            for nn in range(P.n_parity):
                for cs in range(P.n_child_states):
                    gr = g[:, 0, i, j, nn, cs]
                    hr = hR[:, 0, i, j, nn, cs]
                    renter_total += float(np.sum(gr))
                    renter_cap += float(np.sum(gr[np.isclose(hr, P.hR_max, atol=1e-3)]))
                    for name, lo, hi in renter_bins:
                        renter_mass[name] += float(np.sum(gr[(hr > lo) & (hr <= hi)]))
                    for ten in range(1, 1 + P.n_house):
                        mass = float(np.sum(g[:, ten, i, j, nn, cs]))
                        h = float(P.H_own[ten - 1])
                        owner_total += mass
                        if h <= 6.8:
                            lower_middle_owner += mass
                        for name, lo, hi in owner_bins:
                            if h > lo and h <= hi:
                                owner_mass[name] += mass
    for name, mass in renter_mass.items():
        rows.append({"diagnostic": "prime_age_renter_room_share", "bin": name, "value": mass / max(renter_total, 1e-12)})
    for name, mass in owner_mass.items():
        rows.append({"diagnostic": "prime_age_owner_room_share", "bin": name, "value": mass / max(owner_total, 1e-12)})
    rows.append({"diagnostic": "renter_cap_mass", "value": renter_cap / max(renter_total, 1e-12)})
    rows.append({"diagnostic": "lower_middle_owner_rung_mass", "value": lower_middle_owner / max(owner_total, 1e-12)})
    return rows


def write_diagnostics(path: Path, diag: SimpleNamespace, baseline_diag: SimpleNamespace, sol: SimpleNamespace, P: SimpleNamespace, history: list[dict]) -> None:
    rows: list[dict[str, object]] = []
    for i in range(P.I):
        for q, name in enumerate(diag.type_names):
            rows.append(
                {
                    "diagnostic": "type_market_v2",
                    "location": i,
                    "type": name,
                    "p_iq": diag.p_iq[i, q],
                    "r_iq": diag.r_iq[i, q],
                    "stock_supply": diag.supply_iq[i, q],
                    "demand": diag.demand_iq[i, q],
                    "baseline_demand": baseline_diag.demand_iq[i, q],
                    "excess": diag.excess_iq[i, q],
                    "entry_threshold_r": diag.entry_threshold_r_iq[i, q],
                }
            )
    for i in range(P.I):
        rows.append({"diagnostic": "middle_price_premium", "location": i, "value": diag.p_iq[i, 1] - diag.p_iq[i, 0]})
        rows.append({"diagnostic": "large_price_premium", "location": i, "value": diag.p_iq[i, 2] - diag.p_iq[i, 0]})
        rows.append({"diagnostic": "middle_stock_share", "location": i, "value": diag.supply_iq[i, 1] / max(np.sum(diag.supply_iq[i, :]), 1e-12)})
        rows.append({"diagnostic": "middle_demand_share", "location": i, "value": diag.demand_iq[i, 1] / max(np.sum(diag.demand_iq[i, :]), 1e-12)})
    rows.extend(weighted_room_bins(sol, P))
    for item in history:
        rows.append({"diagnostic": "price_iteration_v2", "iteration": item["iter"], "value": item["max_rel_price_step_active"]})
        rows.append({"diagnostic": "active_type_excess_v2", "iteration": item["iter"], "value": item["max_active_excess_rel_demand"]})

    fieldnames = [
        "diagnostic",
        "iteration",
        "location",
        "type",
        "bin",
        "value",
        "p_iq",
        "r_iq",
        "stock_supply",
        "demand",
        "baseline_demand",
        "excess",
        "entry_threshold_r",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, loss: float, moments: dict[str, float], diag: SimpleNamespace, history: list[dict], elapsed: float) -> None:
    active_excess = float(np.max(np.abs(diag.excess_iq[:, 1:]) / np.maximum(np.abs(diag.demand_iq[:, 1:]), 1e-8)))
    last_step = float(history[-1]["max_rel_price_step_active"]) if history else np.nan
    own = moments.get("own_rate", np.nan)
    renter_rooms = moments.get("prime_childless_renter_median_rooms", np.nan)
    h01 = moments.get("housing_increment_0to1", np.nan)
    h12 = moments.get("housing_increment_1to2", np.nan)
    owner_rooms = moments.get("prime_childless_owner_median_rooms", np.nan)
    verdict = "red"
    if np.isfinite(last_step) and last_step < 0.08 and own > 0.20:
        verdict = "yellow"
    if (
        active_excess < 0.10
        and own > 0.45
        and renter_rooms < BENCHMARK_MODEL["prime_childless_renter_median_rooms"]
        and owner_rooms <= BENCHMARK_MODEL["prime_childless_owner_median_rooms"]
        and h01 >= BENCHMARK_MODEL["housing_increment_0to1"]
        and h12 >= BENCHMARK_MODEL["housing_increment_1to2"]
    ):
        verdict = "green"
    text = f"""# Developer Missing-Middle V2 Report

Verdict: **{verdict}**

## What Changed Relative To V1

V1 treated every owner rung above 5 rooms as middle housing. That made the
middle price act like a tax on almost all ownership and collapsed the owner
sector. V2 uses three unit types:

- `S`: non-middle ordinary units, below 5 rooms;
- `M`: missing-middle units in the 5--6.5 room interval;
- `L`: large units above 6.5 rooms.

The V2 loop first measures scalar-price baseline demand by type. It then
chooses \\(\\nu_{{iq}}\\) so each type's supply equals that baseline demand at the
initial type price. This avoids interpreting arbitrary initial stocks as a
developer wedge.

## Status

- active type price iterations: `{len(history)}`
- last active relative price step: `{last_step:.6g}`
- final active type excess relative to demand: `{active_excess:.6g}`
- ordinary price pinned: `True`
- \\(F_{{iq}}\\): reported in entry thresholds but not used for discrete shutdown
- elapsed seconds: `{elapsed:.2f}`
- SMM loss: `{loss:.6g}`

## Moment Movement

| Moment | Target | Benchmark | V2 |
|---|---:|---:|---:|
| childless renter median rooms | 4.000 | {BENCHMARK_MODEL["prime_childless_renter_median_rooms"]:.3f} | {renter_rooms:.3f} |
| childless owner median rooms | 6.000 | {BENCHMARK_MODEL["prime_childless_owner_median_rooms"]:.3f} | {owner_rooms:.3f} |
| H01 | 0.664 | {BENCHMARK_MODEL["housing_increment_0to1"]:.3f} | {h01:.3f} |
| H12 | 0.566 | {BENCHMARK_MODEL["housing_increment_1to2"]:.3f} | {h12:.3f} |
| ownership | 0.627 | {BENCHMARK_MODEL["own_rate"]:.3f} | {own:.3f} |

## Read

This is a better engineering test than V1 because the mapping is no longer
taxing the whole owner ladder as middle housing. It is still not green: the
type-clearing mechanics work, but owner median rooms move the wrong way and
\(H01\) falls relative to the current benchmark. Ordinary prices are also
pinned, the household state is unchanged, and the type mapping is piecewise by
room bins rather than estimated from ACS unit-type prices.
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--nb", type=int, default=50)
    parser.add_argument("--middle-wedge", type=float, default=0.05)
    parser.add_argument("--large-wedge", type=float, default=0.00)
    parser.add_argument("--force-full", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    t0 = time.perf_counter()
    log_path = out_dir / "developer_missing_middle_v2.log"
    try:
        P, targets, weights = build_base_parameters(args.nb, args.force_full)
        sol, P_out, p_eq, diag, history, baseline_diag = run_type_loop(
            P, args.iterations, args.quiet, args.middle_wedge, args.large_wedge
        )
        moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        moments_ns.inv_pop_share_C = float(sol.pop_share[1])
        moments_ns.inv_rent_ratio_C_over_P = float(diag.r_iq[1, 0] / max(diag.r_iq[0, 0], 1e-12))
        moments = namespace_float_dict(moments_ns)
        loss = compute_smm_loss(moments, targets, weights)
        write_results(out_dir / "results_developer_missing_middle_v2.csv", targets, moments)
        write_diagnostics(out_dir / "diagnostics_developer_missing_middle_v2.csv", diag, baseline_diag, sol, P_out, history)
        write_report(out_dir / "REPORT_V2.md", loss, moments, diag, history, time.perf_counter() - t0)
        log_path.write_text(
            json.dumps(
                {
                    "ok": True,
                    "elapsed_sec": time.perf_counter() - t0,
                    "loss": loss,
                    "moments": moments,
                    "history": history,
                    "p_iq": diag.p_iq.tolist(),
                    "r_iq": diag.r_iq.tolist(),
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
