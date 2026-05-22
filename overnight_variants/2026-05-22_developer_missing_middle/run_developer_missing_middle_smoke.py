#!/usr/bin/env python3
"""Smoke test for developer-based missing-middle supply.

This is a partial type-clearing prototype: ordinary prices are pinned to the
current benchmark scalar prices, while middle-unit prices are updated from the
developer inverse supply schedule. It does not run a global calibration.
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
    developer_middle_type,
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

    # Disable kernels that assume one scalar rent per location in the renter
    # block. The copied solver fallback handles type-specific rents by column.
    P.use_full_kernel = False
    P.use_eval_kernel = False

    P.developer_supply = True
    P.developer_n_types = 2
    P.developer_type_names = np.array(["S", "M"], dtype=object)
    P.developer_middle_min = 5.0
    P.developer_middle_max = 6.0
    return P, setup.targets, setup.weights


def initialize_developer_prices(P: SimpleNamespace) -> tuple[np.ndarray, np.ndarray]:
    p_iq = np.column_stack([BENCHMARK_P, BENCHMARK_P * np.array([1.12, 1.22])])
    r_iq = P.user_cost_rate * p_iq

    base_stock = P.H0 * (P.user_cost_rate * BENCHMARK_P / P.r_bar) ** P.eta_supply
    target_stock = np.column_stack([0.72 * base_stock, 0.28 * base_stock])
    eta_iq = np.column_stack([P.eta_supply, np.array([1.15, 0.55])])
    tau_iq = np.column_stack([np.zeros(P.I), 0.18 * r_iq[:, 0]])
    nu_iq = np.maximum((r_iq - tau_iq) / np.maximum(target_stock, 1e-8) ** (1.0 / eta_iq), 1e-8)
    F_iq = np.column_stack([np.zeros(P.I), np.array([0.002, 0.004])])

    P.developer_tau_iq = tau_iq
    P.developer_eta_iq = eta_iq
    P.developer_nu_iq = nu_iq
    P.developer_F_iq = F_iq
    return p_iq, r_iq


def run_partial_type_loop(P: SimpleNamespace, iterations: int, quiet: bool) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray, SimpleNamespace, list[dict]]:
    p_iq, r_iq = initialize_developer_prices(P)
    history: list[dict] = []
    sol = P_out = p_eq = diag = None
    mid = developer_middle_type(P)
    for it in range(1, iterations + 1):
        P_iter = SimpleNamespace(**deepcopy(asdict(P)))
        P_iter.developer_p_iq = p_iq.copy()
        P_iter.developer_r_iq = r_iq.copy()
        P_iter.p_fixed = p_iq[:, 0].copy()
        sol, P_out, p_eq = run_model_cp_dt(P_iter, verbose=not quiet)
        diag = compute_developer_type_diagnostics(sol, P_out, p_iq, r_iq)
        target_r = developer_inverse_supply_rent(P_out, diag.demand_iq)
        old_middle = r_iq[:, mid].copy()
        r_iq[:, mid] = 0.55 * r_iq[:, mid] + 0.45 * target_r[:, mid]
        r_iq[:, mid] = np.maximum(r_iq[:, mid], 1.01 * r_iq[:, 0])
        r_iq[:, mid] = np.minimum(r_iq[:, mid], 2.5 * r_iq[:, 0])
        p_iq[:, mid] = r_iq[:, mid] / P_out.user_cost_rate
        max_rel_step = float(np.max(np.abs(r_iq[:, mid] - old_middle) / np.maximum(old_middle, 1e-12)))
        history.append(
            {
                "iter": it,
                "p_iq": p_iq.tolist(),
                "r_iq": r_iq.tolist(),
                "demand_iq": diag.demand_iq.tolist(),
                "supply_iq": diag.supply_iq.tolist(),
                "max_rel_middle_price_step": max_rel_step,
            }
        )
        if max_rel_step < 0.02:
            break
    assert sol is not None and P_out is not None and p_eq is not None and diag is not None
    final_P = SimpleNamespace(**deepcopy(asdict(P_out)))
    final_P.developer_p_iq = p_iq.copy()
    final_P.developer_r_iq = r_iq.copy()
    final_diag = compute_developer_type_diagnostics(sol, final_P, p_iq, r_iq)
    return sol, final_P, p_eq, final_diag, history


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
        writer = csv.DictWriter(f, fieldnames=["moment", "target", "model", "benchmark_model", "notes"])
        writer.writeheader()
        for name in order:
            writer.writerow(
                {
                    "moment": name,
                    "target": targets.get(name, ""),
                    "model": model.get(name, np.nan),
                    "benchmark_model": BENCHMARK_MODEL.get(name, ""),
                    "notes": "partial type-clearing smoke; ordinary prices pinned",
                }
            )


def room_screen_rows(sol: SimpleNamespace, P: SimpleNamespace) -> list[dict[str, object]]:
    g = sol.g
    hR = sol.hR_pol
    rows: list[dict[str, object]] = []
    a25s = max(0, round(25 - P.age_start))
    a45e = min(P.J - 1, round(45 - P.age_start))
    renter_bins = [("<=4", -np.inf, 4.0), ("5", 4.0, 5.0), ("6", 5.0, 6.0), ("7-8", 6.0, 8.0)]
    owner_bins = [("5", 4.0, 5.5), ("6", 5.5, 6.5), ("7-8", 6.5, 8.5), ("9-10", 8.5, 10.5), ("11+", 10.5, np.inf)]
    renter_total = 0.0
    renter_bin_mass = {name: 0.0 for name, _, _ in renter_bins}
    cap_mass = 0.0
    owner_total = 0.0
    owner_bin_mass = {name: 0.0 for name, _, _ in owner_bins}
    lower_middle_owner = 0.0
    for j in range(a25s, a45e + 1):
        for i in range(P.I):
            for nn in range(P.n_parity):
                for cs in range(P.n_child_states):
                    gr = g[:, 0, i, j, nn, cs]
                    hr = hR[:, 0, i, j, nn, cs]
                    renter_total += float(np.sum(gr))
                    cap_mass += float(np.sum(gr[np.isclose(hr, P.hR_max, atol=1e-3)]))
                    for name, lo, hi in renter_bins:
                        mask = (hr > lo) & (hr <= hi)
                        renter_bin_mass[name] += float(np.sum(gr[mask]))
                    for ten in range(1, 1 + P.n_house):
                        mass = float(np.sum(g[:, ten, i, j, nn, cs]))
                        owner_total += mass
                        h = float(P.H_own[ten - 1])
                        if ten <= 3:
                            lower_middle_owner += mass
                        for name, lo, hi in owner_bins:
                            if h > lo and h <= hi:
                                owner_bin_mass[name] += mass
    for name, mass in renter_bin_mass.items():
        rows.append({"diagnostic": "prime_age_renter_room_share", "bin": name, "value": mass / max(renter_total, 1e-12)})
    for name, mass in owner_bin_mass.items():
        rows.append({"diagnostic": "prime_age_owner_room_share", "bin": name, "value": mass / max(owner_total, 1e-12)})
    rows.append({"diagnostic": "renter_cap_mass", "value": cap_mass / max(renter_total, 1e-12)})
    rows.append({"diagnostic": "lower_middle_owner_rung_mass", "value": lower_middle_owner / max(owner_total, 1e-12)})
    return rows


def mechanism_rows(sol: SimpleNamespace, P: SimpleNamespace) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    g = sol.g
    a22s = max(0, round(22 - P.age_start))
    a45e = min(P.J - 1, round(45 - P.age_start))
    middle_owner_mass = 0.0
    middle_rental_mass = 0.0
    total_mass = float(np.sum(g))
    h2_owner_from_renter = 0.0
    feasible_but_rent = 0.0
    feasible_mass = 0.0
    p_iq = np.asarray(P.developer_p_iq)
    b_grid = make_grid(P)
    for j in range(a22s, a45e + 1):
        for i in range(P.I):
            for nn in range(P.n_parity):
                for cs in range(P.n_child_states):
                    gr = g[:, 0, i, j, nn, cs]
                    hr = sol.hR_pol[:, 0, i, j, nn, cs]
                    middle_rental_mass += float(np.sum(gr[hr >= P.developer_middle_min]))
                    for ten in range(1, 1 + P.n_house):
                        mass = float(np.sum(g[:, ten, i, j, nn, cs]))
                        if P.H_own[ten - 1] >= P.developer_middle_min:
                            middle_owner_mass += mass
                    if nn >= 1 and cs >= 1 and cs <= P.n_child_stages:
                        tc = sol.tenure_choice[:, 0, i, j, nn, cs]
                        h2 = 2 if P.n_house >= 2 else 1
                        h2_price = p_iq[i, 1] * P.H_own[h2 - 1]
                        dp = (1.0 - P.phi[min(nn, len(P.phi) - 1)]) * h2_price
                        feas = gr > 0
                        feas &= b_grid[0 : len(gr)] >= dp
                        feasible_mass += float(np.sum(gr[feas]))
                        feasible_but_rent += float(np.sum(gr[feas & (tc == 0)]))
                        h2_owner_from_renter += float(np.sum(gr[tc == h2]))
    rows.append({"diagnostic": "mass_choosing_middle_rental", "value": middle_rental_mass / max(total_mass, 1e-12)})
    rows.append({"diagnostic": "mass_choosing_middle_owner_rung", "value": middle_owner_mass / max(total_mass, 1e-12)})
    rows.append({"diagnostic": "h2_starter_family_margin_mass", "value": h2_owner_from_renter})
    rows.append({"diagnostic": "new_parent_feasible_but_rent_share", "value": feasible_but_rent / max(feasible_mass, 1e-12)})
    return rows


def write_diagnostics(path: Path, diag: SimpleNamespace, sol: SimpleNamespace, P: SimpleNamespace, history: list[dict]) -> None:
    rows: list[dict[str, object]] = []
    for i in range(P.I):
        for q, name in enumerate(diag.type_names):
            rows.append(
                {
                    "diagnostic": "type_market",
                    "location": i,
                    "type": name,
                    "p_iq": diag.p_iq[i, q],
                    "r_iq": diag.r_iq[i, q],
                    "stock_supply": diag.supply_iq[i, q],
                    "demand": diag.demand_iq[i, q],
                    "excess": diag.excess_iq[i, q],
                    "entry_threshold_r": diag.entry_threshold_r_iq[i, q],
                }
            )
    mid = developer_middle_type(P)
    for i in range(P.I):
        rows.append({"diagnostic": "middle_rent_wedge", "location": i, "value": diag.r_iq[i, mid] - diag.r_iq[i, 0]})
        rows.append({"diagnostic": "middle_price_wedge", "location": i, "value": diag.p_iq[i, mid] - diag.p_iq[i, 0]})
        rows.append({"diagnostic": "middle_stock_share", "location": i, "value": diag.middle_stock_share[i]})
        rows.append({"diagnostic": "middle_demand_share", "location": i, "value": diag.middle_demand_share[i]})
    rows.extend(room_screen_rows(sol, P))
    rows.extend(mechanism_rows(sol, P))
    for item in history:
        rows.append({"diagnostic": "type_price_iteration", "iteration": item["iter"], "value": item["max_rel_middle_price_step"]})

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
        "excess",
        "entry_threshold_r",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    loss: float,
    moments: dict[str, float],
    diag: SimpleNamespace,
    history: list[dict],
    elapsed: float,
) -> None:
    max_excess = float(np.max(np.abs(diag.excess_iq) / np.maximum(diag.demand_iq, 1e-8)))
    last_step = float(history[-1]["max_rel_middle_price_step"]) if history else np.nan
    pinned = True
    ownership_ok = moments.get("own_rate", 0.0) > 0.05
    room_improved = (
        moments.get("prime_childless_renter_median_rooms", np.inf) < BENCHMARK_MODEL["prime_childless_renter_median_rooms"]
        or moments.get("prime_childless_owner_median_rooms", np.inf) < BENCHMARK_MODEL["prime_childless_owner_median_rooms"]
    )
    verdict = "red"
    if np.isfinite(last_step) and last_step < 0.10 and ownership_ok and max_excess < 0.50:
        verdict = "yellow"
    if (not pinned) and max_excess < 0.05 and room_improved:
        verdict = "green"
    text = f"""# Developer Missing-Middle Smoke Report

Verdict: **{verdict}**

## What Ran

The copied solver was modified so owner budgets use \(p_{{iq(k)}}\), owner
maintenance uses \(p_{{iq(k)}}\), sale/down-payment accounting uses the same
type price, and renter choices use a two-type rent approximation based on the
state's housing need. A small outer loop pinned ordinary prices to the current
benchmark and updated only middle-unit prices from the developer inverse supply
schedule.

This is therefore a **partial** type-specific market-clearing test, not a full
2-location by 2-type general-equilibrium fixed point.

## Developer Block

- active unit types: `S` ordinary and `M` family-capable middle
- middle boundary: `{getattr(diag, "middle_boundary", "5+ rooms")}`
- \(F_{{iq}}\): documented and included in entry-threshold diagnostics, but not
  used to shut types off in the first clearing loop
- operational wedges: \\(\\tau_{{iM}}\\), \\(\\nu_{{iM}}\\), and \\(\\eta_{{iM}}\\)
- ordinary prices pinned: `{pinned}`

## Clearing Status

- middle-price iterations: `{len(history)}`
- last relative middle-price step: `{last_step:.6g}`
- final max type excess relative to demand: `{max_excess:.6g}`
- elapsed seconds: `{elapsed:.2f}`
- runtime category: `moderate-to-expensive` for partial PE type loop; full GE
  type clearing remains `expensive`
- SMM loss on partial smoke moments: `{loss:.6g}`

## Key Moment Movement

| Moment | Target | Benchmark | Variant |
|---|---:|---:|---:|
| childless renter median rooms | 4.000 | {BENCHMARK_MODEL["prime_childless_renter_median_rooms"]:.3f} | {moments.get("prime_childless_renter_median_rooms", np.nan):.3f} |
| childless owner median rooms | 6.000 | {BENCHMARK_MODEL["prime_childless_owner_median_rooms"]:.3f} | {moments.get("prime_childless_owner_median_rooms", np.nan):.3f} |
| H01 | 0.664 | {BENCHMARK_MODEL["housing_increment_0to1"]:.3f} | {moments.get("housing_increment_0to1", np.nan):.3f} |
| H12 | 0.566 | {BENCHMARK_MODEL["housing_increment_1to2"]:.3f} | {moments.get("housing_increment_1to2", np.nan):.3f} |
| ownership | 0.627 | {BENCHMARK_MODEL["own_rate"]:.3f} | {moments.get("own_rate", np.nan):.3f} |

## Interpretation

The branch is red in this smoke test because the partial type loop did not
deliver a usable market allocation: ordinary prices were pinned, type excess
remained large, and ownership collapsed. The diagnostic CSV reports type
prices, stocks, demands, middle wedges, room bins, renter cap mass,
lower/middle owner mass, and whether the H2/starter-family margin is alive.
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=4)
    parser.add_argument("--nb", type=int, default=40)
    parser.add_argument("--force-full", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    log_path = out_dir / "developer_missing_middle_smoke.log"
    t0 = time.perf_counter()
    try:
        P, targets, weights = build_base_parameters(args.nb, args.force_full)
        sol, P_out, p_eq, diag, history = run_partial_type_loop(P, args.iterations, args.quiet)
        moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        moments_ns.inv_pop_share_C = float(sol.pop_share[1])
        moments_ns.inv_rent_ratio_C_over_P = float(diag.r_iq[1, 0] / max(diag.r_iq[0, 0], 1e-12))
        moments = namespace_float_dict(moments_ns)
        loss = compute_smm_loss(moments, targets, weights)
        write_results(out_dir / "results_developer_missing_middle.csv", targets, moments)
        write_diagnostics(out_dir / "diagnostics_developer_missing_middle.csv", diag, sol, P_out, history)
        write_report(out_dir / "REPORT.md", loss, moments, diag, history, time.perf_counter() - t0)
        log_path.write_text(
            json.dumps(
                {
                    "ok": True,
                    "elapsed_sec": time.perf_counter() - t0,
                    "loss": loss,
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
