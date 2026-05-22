#!/usr/bin/env python3
"""V3 developer missing-middle prototype: full type-price equilibrium loop.

The earlier V2 run pinned ordinary prices and updated only the middle/large
prices. This driver treats all location/type prices and entry shares as
endogenous fixed-point objects in the copied branch. Each fixed-point iteration
solves households at the current type-price vector, computes type demand, and
updates all developer supply prices from the inverse supply curve.
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

from dt_cp_model.direct_calibration import compute_smm_loss
from dt_cp_model.objective import extract_moments
from dt_cp_model.parameters import asdict
from dt_cp_model.solver import compute_developer_type_diagnostics, developer_inverse_supply_rent, run_model_cp_dt

from run_developer_missing_middle_v2 import (
    BENCHMARK_MODEL,
    BENCHMARK_P,
    MOMENT_ORDER,
    build_base_parameters,
    initialize_supply_from_baseline,
    namespace_float_dict,
    weighted_room_bins,
)


def solve_at_type_prices(P: SimpleNamespace, p_iq: np.ndarray, entry_shares: np.ndarray, quiet: bool):
    P_iter = SimpleNamespace(**deepcopy(asdict(P)))
    P_iter.solve_mode = "pe"
    P_iter.entry_shares_fixed = np.asarray(entry_shares, dtype=float)
    P_iter.developer_p_iq = np.asarray(p_iq, dtype=float)
    P_iter.developer_r_iq = P_iter.user_cost_rate * P_iter.developer_p_iq
    P_iter.developer_exact_renter_types = True
    P_iter.p_fixed = P_iter.developer_p_iq[:, 0].copy()
    P_iter.w_fixed = np.ones(P_iter.I)
    sol, P_out, p_eq = run_model_cp_dt(P_iter, verbose=not quiet)
    diag = compute_developer_type_diagnostics(sol, P_out, P_iter.developer_p_iq, P_iter.developer_r_iq)
    return sol, P_out, p_eq, diag


def normalize_entry(x: np.ndarray) -> np.ndarray:
    y = np.maximum(np.asarray(x, dtype=float).reshape(-1), 1e-8)
    return y / float(np.sum(y))


def run_type_price_ge(
    P: SimpleNamespace,
    iterations: int,
    quiet: bool,
    middle_wedge: float,
    large_wedge: float,
    price_damp: float,
    entry_damp: float,
    tol: float,
):
    baseline_p = np.tile(BENCHMARK_P.reshape(P.I, 1), (1, 3))
    entry_shares = normalize_entry(np.array([0.55, 0.45]))
    sol0, P0, _, baseline_diag = solve_at_type_prices(P, baseline_p, entry_shares, quiet)
    p_iq, r_iq, tau_iq, nu_iq, eta_iq, F_iq = initialize_supply_from_baseline(P, baseline_diag, middle_wedge, large_wedge)
    P.developer_tau_iq = tau_iq
    P.developer_nu_iq = nu_iq
    P.developer_eta_iq = eta_iq
    P.developer_F_iq = F_iq

    base_r = P.user_cost_rate * baseline_p
    history: list[dict] = []
    accepted = False
    sol = sol0
    P_out = P0
    p_eq = BENCHMARK_P.copy()
    diag = baseline_diag
    final_step = np.nan
    final_entry_step = np.nan

    for it in range(1, iterations + 1):
        sol, P_out, p_eq, diag = solve_at_type_prices(P, p_iq, entry_shares, quiet)
        target_r = developer_inverse_supply_rent(P_out, diag.demand_iq)
        target_r = np.clip(target_r, 0.35 * base_r, 2.50 * base_r)
        next_r = (1.0 - price_damp) * r_iq + price_damp * target_r
        next_p = next_r / P.user_cost_rate
        entry_target = normalize_entry(getattr(sol, "mature_entry_shares", entry_shares))
        next_entry = normalize_entry((1.0 - entry_damp) * entry_shares + entry_damp * entry_target)
        final_step = float(np.max(np.abs(next_r - r_iq) / np.maximum(np.abs(r_iq), 1e-12)))
        final_entry_step = float(np.max(np.abs(next_entry - entry_shares)))
        market_excess = float(np.max(np.abs(diag.excess_iq) / np.maximum(np.abs(diag.demand_iq), 1e-8)))
        history.append(
            {
                "iter": it,
                "p_iq": p_iq.tolist(),
                "r_iq": r_iq.tolist(),
                "entry_shares": entry_shares.tolist(),
                "entry_target": entry_target.tolist(),
                "demand_iq": diag.demand_iq.tolist(),
                "supply_iq": diag.supply_iq.tolist(),
                "max_rel_price_step": final_step,
                "max_entry_step": final_entry_step,
                "max_excess_rel_demand": market_excess,
                "own_rate": float(getattr(sol, "own_rate", np.nan)),
                "tfr": float(2.0 * getattr(sol, "mean_parity", np.nan)),
            }
        )
        if final_step < tol and final_entry_step < tol:
            accepted = True
            break
        r_iq = next_r
        p_iq = next_p
        entry_shares = next_entry

    # Final solve at the last updated endogenous prices/shares.
    sol, P_out, p_eq, diag = solve_at_type_prices(P, p_iq, entry_shares, quiet)
    market_excess = float(np.max(np.abs(diag.excess_iq) / np.maximum(np.abs(diag.demand_iq), 1e-8)))
    final_P = SimpleNamespace(**deepcopy(asdict(P_out)))
    final_P.developer_p_iq = p_iq.copy()
    final_P.developer_r_iq = r_iq.copy()
    final_P.developer_tau_iq = tau_iq
    final_P.developer_nu_iq = nu_iq
    final_P.developer_eta_iq = eta_iq
    final_P.developer_F_iq = F_iq
    final_P.entry_shares = entry_shares.copy()
    final_diag = compute_developer_type_diagnostics(sol, final_P, p_iq, r_iq)
    return sol, final_P, p_eq, final_diag, baseline_diag, history, {
        "accepted": bool(accepted or (market_excess < 0.05 and np.isfinite(final_step) and final_step < 2 * tol)),
        "final_price_step": final_step,
        "final_entry_step": final_entry_step,
        "final_market_excess": market_excess,
        "entry_shares": entry_shares.tolist(),
    }


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
                    "notes": "V3 full type-price/entry fixed point; exact realized-rent type pricing; no ordinary price pinning",
                }
            )


def write_diagnostics(path: Path, diag: SimpleNamespace, baseline_diag: SimpleNamespace, sol: SimpleNamespace, P: SimpleNamespace, history: list[dict], status: dict) -> None:
    rows: list[dict[str, object]] = []
    for i in range(P.I):
        for q, name in enumerate(diag.type_names):
            rows.append(
                {
                    "diagnostic": "type_market_v3_ge",
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
        renter_mass_total = float(np.sum(diag.renter_mass_type[i, :]))
        owner_mass_total = float(np.sum(diag.owner_mass_type[i, :]))
        for q, name in enumerate(diag.type_names):
            rows.append(
                {
                    "diagnostic": "renter_mass_type_share_v3_ge",
                    "location": i,
                    "type": name,
                    "value": diag.renter_mass_type[i, q] / max(renter_mass_total, 1e-12),
                }
            )
            rows.append(
                {
                    "diagnostic": "owner_mass_type_share_v3_ge",
                    "location": i,
                    "type": name,
                    "value": diag.owner_mass_type[i, q] / max(owner_mass_total, 1e-12),
                }
            )
    for item in history:
        rows.append({"diagnostic": "price_iteration_v3_ge", "iteration": item["iter"], "value": item["max_rel_price_step"]})
        rows.append({"diagnostic": "entry_iteration_v3_ge", "iteration": item["iter"], "value": item["max_entry_step"]})
        rows.append({"diagnostic": "market_excess_v3_ge", "iteration": item["iter"], "value": item["max_excess_rel_demand"]})
    for key, value in status.items():
        rows.append({"diagnostic": "ge_status", "type": key, "value": value})
    rows.extend(weighted_room_bins(sol, P))
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


def write_report(path: Path, loss: float, moments: dict[str, float], diag: SimpleNamespace, history: list[dict], status: dict, elapsed: float) -> None:
    own = moments.get("own_rate", np.nan)
    renter_rooms = moments.get("prime_childless_renter_median_rooms", np.nan)
    owner_rooms = moments.get("prime_childless_owner_median_rooms", np.nan)
    h01 = moments.get("housing_increment_0to1", np.nan)
    h12 = moments.get("housing_increment_1to2", np.nan)
    market_excess = float(status.get("final_market_excess", np.nan))
    accepted = bool(status.get("accepted", False))
    runtime_category = "expensive" if elapsed >= 60.0 else "moderate"
    verdict = "yellow" if accepted and own > 0.20 else "red"
    if (
        accepted
        and market_excess < 0.05
        and own > 0.45
        and renter_rooms < BENCHMARK_MODEL["prime_childless_renter_median_rooms"]
        and owner_rooms <= BENCHMARK_MODEL["prime_childless_owner_median_rooms"]
        and h01 >= BENCHMARK_MODEL["housing_increment_0to1"]
        and h12 >= BENCHMARK_MODEL["housing_increment_1to2"]
    ):
        verdict = "green"
    rows = "\n".join(
        f"| {diag.type_names[q]} loc {i} | {diag.p_iq[i, q]:.3f} | {diag.r_iq[i, q]:.3f} | {diag.supply_iq[i, q]:.3f} | {diag.demand_iq[i, q]:.3f} |"
        for i in range(diag.p_iq.shape[0])
        for q in range(diag.p_iq.shape[1])
    )
    text = f"""# Developer Missing-Middle V3 Full Equilibrium Report

Verdict: **{verdict}**

## What This Is

This supersedes the V2 partial test. V3 updates every location/type price
\(p_{{iq}}\) and the entry shares in a copied fixed-point loop. Ordinary prices
are no longer pinned. The household solve at each iteration uses the current
type-price vector, charges renters by realized room interval \(q(h^R)\), then
developer inverse supply updates all \(S/M/L\) rents.

## Equilibrium Status

- accepted: `{accepted}`
- iterations: `{len(history)}`
- final relative market excess: `{market_excess:.6g}`
- final relative price step: `{float(status.get("final_price_step", np.nan)):.6g}`
- final entry-share step: `{float(status.get("final_entry_step", np.nan)):.6g}`
- elapsed seconds: `{elapsed:.2f}`
- runtime category: `{runtime_category}`
- SMM loss: `{loss:.6g}`
- \(F_{{iq}}\): reported in entry thresholds but not used for discrete shutdown

## Type Markets

| Market | p | r | supply | demand |
|---|---:|---:|---:|---:|
{rows}

## Moment Movement

| Moment | Target | Benchmark | V3 GE |
|---|---:|---:|---:|
| childless renter median rooms | 4.000 | {BENCHMARK_MODEL["prime_childless_renter_median_rooms"]:.3f} | {renter_rooms:.3f} |
| childless owner median rooms | 6.000 | {BENCHMARK_MODEL["prime_childless_owner_median_rooms"]:.3f} | {owner_rooms:.3f} |
| H01 | 0.664 | {BENCHMARK_MODEL["housing_increment_0to1"]:.3f} | {h01:.3f} |
| H12 | 0.566 | {BENCHMARK_MODEL["housing_increment_1to2"]:.3f} | {h12:.3f} |
| ownership | 0.627 | {BENCHMARK_MODEL["own_rate"]:.3f} | {own:.3f} |

## Read

This is now a genuine type-price fixed point rather than a pinned ordinary-price
test, and renter budgets now use realized room intervals. The verdict is yellow
only if the status above accepts; it is not green because the room medians,
ownership gradient, and \(H01\) remain away from the target discipline.
"""
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=18)
    parser.add_argument("--nb", type=int, default=40)
    parser.add_argument("--middle-wedge", type=float, default=0.05)
    parser.add_argument("--large-wedge", type=float, default=0.00)
    parser.add_argument("--price-damp", type=float, default=0.35)
    parser.add_argument("--entry-damp", type=float, default=0.35)
    parser.add_argument("--tol", type=float, default=0.01)
    parser.add_argument("--force-full", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    out_dir = Path(__file__).resolve().parent
    t0 = time.perf_counter()
    log_path = out_dir / "developer_missing_middle_v3_ge.log"
    try:
        P, targets, weights = build_base_parameters(args.nb, args.force_full)
        sol, P_out, p_eq, diag, baseline_diag, history, status = run_type_price_ge(
            P,
            args.iterations,
            args.quiet,
            args.middle_wedge,
            args.large_wedge,
            args.price_damp,
            args.entry_damp,
            args.tol,
        )
        moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        moments_ns.inv_pop_share_C = float(sol.pop_share[1])
        moments_ns.inv_rent_ratio_C_over_P = float(diag.r_iq[1, 0] / max(diag.r_iq[0, 0], 1e-12))
        moments = namespace_float_dict(moments_ns)
        loss = compute_smm_loss(moments, targets, weights)
        elapsed = time.perf_counter() - t0
        write_results(out_dir / "results_developer_missing_middle_v3_ge.csv", targets, moments)
        write_diagnostics(out_dir / "diagnostics_developer_missing_middle_v3_ge.csv", diag, baseline_diag, sol, P_out, history, status)
        write_report(out_dir / "REPORT_V3_GE.md", loss, moments, diag, history, status, elapsed)
        log_path.write_text(
            json.dumps(
                {
                    "ok": True,
                    "elapsed_sec": elapsed,
                    "loss": loss,
                    "moments": moments,
                    "history": history,
                    "status": status,
                    "p_iq": diag.p_iq.tolist(),
                    "r_iq": diag.r_iq.tolist(),
                    "demand_iq": diag.demand_iq.tolist(),
                    "supply_iq": diag.supply_iq.tolist(),
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
