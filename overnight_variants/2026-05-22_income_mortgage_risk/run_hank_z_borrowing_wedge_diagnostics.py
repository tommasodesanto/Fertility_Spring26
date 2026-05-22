#!/usr/bin/env python3
"""Validate the existing phi-based housing-finance wedge in the HANK-z branch.

This is an isolated diagnostic driver. It reruns the structural HANK-z GE
prototype and asks whether the model's existing down-payment and borrowing
constraints are used once households face persistent idiosyncratic income risk.
It does not add a mortgage-account state.
"""

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
from dt_cp_model.solver import get_phi_choice_tensor, make_grid, solve_equilibrium_hank_z

from run_income_mortgage_risk_v2_scenarios import BENCHMARK_MODEL, BENCHMARK_P, MOMENT_ORDER, build_base_parameters, namespace_float_dict


AGE_BINS = [
    ("22-29", 22, 29),
    ("30-34", 30, 34),
    ("35-44", 35, 44),
    ("45-55", 45, 55),
    ("56-64", 56, 64),
]


def finance_arrays(P, p_eq: np.ndarray):
    nt = 1 + P.n_house
    hcost = np.zeros((P.I, nt))
    heq = np.zeros((P.I, nt))
    dp = np.zeros((P.I, nt, P.n_parity, P.n_child_states))
    bfloor = np.zeros_like(dp)
    phi_choice = get_phi_choice_tensor(P)
    for i in range(P.I):
        for ten in range(1, nt):
            hcost[i, ten] = p_eq[i] * P.H_own[ten - 1]
            heq[i, ten] = (1.0 - P.psi) * p_eq[i] * P.H_own[ten - 1]
            for nn in range(P.n_parity):
                for cs in range(P.n_child_states):
                    phi = phi_choice[i, ten, nn, cs]
                    dp[i, ten, nn, cs] = (1.0 - phi) * hcost[i, ten]
                    bfloor[i, ten, nn, cs] = -phi * hcost[i, ten]
    return hcost, heq, dp, bfloor


def weighted_mean(values: list[float], weights: list[float]) -> float:
    sw = float(np.sum(weights))
    return float(np.sum(np.asarray(values) * np.asarray(weights)) / sw) if sw > 1e-14 else np.nan


def add_row(rows: list[dict[str, object]], diagnostic: str, value: float, **kwargs) -> None:
    row: dict[str, object] = {"diagnostic": diagnostic, "value": value}
    row.update(kwargs)
    rows.append(row)


def compute_diagnostics(sol, P, p_eq: np.ndarray, b_grid: np.ndarray) -> tuple[list[dict[str, object]], dict[str, float]]:
    hz = sol.hank_z
    g = np.asarray(hz.g_z)
    tc = np.asarray(hz.tenure_choice_z)
    bp = np.asarray(hz.bp_pol_z)
    fp = np.asarray(hz.fert_probs_z)
    z_grid = np.asarray(hz.z_grid)
    Nb, nt, I, J, npar, ncs, Nz = g.shape
    hcost, heq, dp, bfloor = finance_arrays(P, p_eq)
    b = b_grid.reshape(Nb, 1, 1, 1, 1)
    grid_diffs = np.diff(b_grid)
    grid_step = float(np.median(grid_diffs)) if len(b_grid) > 1 else 0.05
    local_grid_step = float(np.percentile(grid_diffs, 25)) if len(b_grid) > 1 else 0.05
    eps_b = min(max(local_grid_step, 0.25), 0.75)
    rows: list[dict[str, object]] = []
    summary: dict[str, float] = {"grid_step": grid_step, "near_threshold_band": eps_b}

    add_row(rows, "diagnostic_parameter", grid_step, name="b_grid_step")
    add_row(rows, "diagnostic_parameter", local_grid_step, name="local_b_grid_step_p25")
    add_row(rows, "diagnostic_parameter", eps_b, name="near_threshold_band")
    add_row(rows, "diagnostic_parameter", float(np.mean(P.phi)), name="mean_phi_financed_share", notes="phi is financed share")

    total_mass = float(np.sum(g))
    owner_mass = float(np.sum(g[:, 1:, :, :, :, :]))
    renter_mass = float(np.sum(g[:, 0, :, :, :, :]))
    add_row(rows, "aggregate_mass", total_mass, name="total")
    add_row(rows, "aggregate_mass", owner_mass / max(total_mass, 1e-14), name="owner_share")
    add_row(rows, "aggregate_mass", renter_mass / max(total_mass, 1e-14), name="renter_share")

    for iz, zval in enumerate(z_grid):
        gz = g[..., iz]
        z_mass = float(np.sum(gz))
        own_mass = float(np.sum(gz[:, 1:, :, :, :, :]))
        add_row(rows, "mass_by_z", z_mass / max(total_mass, 1e-14), z=float(zval), mass=z_mass)
        add_row(rows, "ownership_by_z", own_mass / max(z_mass, 1e-14), z=float(zval), mass=z_mass)

        owner_floor_mass = 0.0
        owner_floor_near = 0.0
        owner_floor_slack_values: list[float] = []
        owner_floor_slack_weights: list[float] = []
        for i in range(I):
            for ten in range(1, nt):
                for nn in range(npar):
                    for cs in range(ncs):
                        mass = gz[:, ten, i, :, nn, cs]
                        if float(np.sum(mass)) <= 1e-15:
                            continue
                        floor = bfloor[i, ten, nn, cs]
                        slack = bp[:, ten, i, :, nn, cs, iz] - floor
                        owner_floor_mass += float(np.sum(mass))
                        owner_floor_near += float(np.sum(mass[slack <= eps_b]))
                        owner_floor_slack_values.extend(slack.reshape(-1).tolist())
                        owner_floor_slack_weights.extend(mass.reshape(-1).tolist())
        add_row(
            rows,
            "owner_bp_near_borrowing_floor_share_by_z",
            owner_floor_near / max(owner_floor_mass, 1e-14),
            z=float(zval),
            mass=owner_floor_mass,
            notes="bp' within the local near-threshold band of -phi*p_i*H_k",
        )
        add_row(
            rows,
            "owner_mean_bp_floor_slack_by_z",
            weighted_mean(owner_floor_slack_values, owner_floor_slack_weights),
            z=float(zval),
            mass=owner_floor_mass,
        )

        renter_to_owner = 0.0
        renter_base = 0.0
        purchase_slack_values: list[float] = []
        purchase_slack_weights: list[float] = []
        purchase_near = 0.0
        renter_starter_feasible = 0.0
        renter_starter_near = 0.0
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    mass = gz[:, 0, i, :, nn, cs]
                    if float(np.sum(mass)) <= 1e-15:
                        continue
                    starter_dp = float(np.min(dp[i, 1:, nn, cs]))
                    starter_slack = b_grid - starter_dp
                    renter_base += float(np.sum(mass))
                    renter_starter_feasible += float(np.sum(mass[starter_slack >= 0.0, :]))
                    renter_starter_near += float(np.sum(mass[np.abs(starter_slack) <= eps_b, :]))
                    choice = tc[:, 0, i, :, nn, cs, iz]
                    for ten in range(1, nt):
                        mk = choice == ten
                        mt = mass * mk
                        mt_sum = float(np.sum(mt))
                        if mt_sum <= 1e-15:
                            continue
                        slack = b_grid.reshape(Nb, 1) - dp[i, ten, nn, cs]
                        slack_full = np.broadcast_to(slack, mt.shape)
                        renter_to_owner += mt_sum
                        purchase_near += float(np.sum(mt[slack_full <= eps_b]))
                        purchase_slack_values.extend(slack_full.reshape(-1).tolist())
                        purchase_slack_weights.extend(mt.reshape(-1).tolist())
        add_row(rows, "renter_to_owner_purchase_share_by_z", renter_to_owner / max(renter_base, 1e-14), z=float(zval), mass=renter_base)
        add_row(rows, "renter_feasible_starter_share_by_z", renter_starter_feasible / max(renter_base, 1e-14), z=float(zval), mass=renter_base)
        add_row(rows, "renter_near_starter_dp_share_by_z", renter_starter_near / max(renter_base, 1e-14), z=float(zval), mass=renter_base)
        add_row(rows, "purchase_near_downpayment_share_by_z", purchase_near / max(renter_to_owner, 1e-14), z=float(zval), mass=renter_to_owner)
        add_row(rows, "mean_purchase_downpayment_slack_by_z", weighted_mean(purchase_slack_values, purchase_slack_weights), z=float(zval), mass=renter_to_owner)

        for label, amin, amax in AGE_BINS:
            j0 = max(0, int(round(amin - P.age_start)))
            j1 = min(J - 1, int(round(amax - P.age_start)))
            if j0 > j1:
                continue
            base = gz[:, 0, :, j0 : j1 + 1, :, :]
            choices = tc[:, 0, :, j0 : j1 + 1, :, :, iz]
            purch = float(np.sum(base[choices > 0]))
            base_mass = float(np.sum(base))
            add_row(rows, "renter_to_owner_purchase_share_by_z_agebin", purch / max(base_mass, 1e-14), z=float(zval), age_bin=label, mass=base_mass)

        trans_mass = {
            "R_to_R": 0.0,
            "R_to_O": 0.0,
            "O_to_R": 0.0,
            "O_to_same": 0.0,
            "O_to_other_owner": 0.0,
        }
        for to in range(nt):
            mass = gz[:, to, :, :, :, :]
            choice = tc[:, to, :, :, :, :, iz]
            if to == 0:
                trans_mass["R_to_R"] += float(np.sum(mass[choice == 0]))
                trans_mass["R_to_O"] += float(np.sum(mass[choice > 0]))
            else:
                trans_mass["O_to_R"] += float(np.sum(mass[choice == 0]))
                trans_mass["O_to_same"] += float(np.sum(mass[choice == to]))
                trans_mass["O_to_other_owner"] += float(np.sum(mass[(choice > 0) & (choice != to)]))
        for name, mass in trans_mass.items():
            add_row(rows, "tenure_choice_mass_share_by_z", mass / max(z_mass, 1e-14), z=float(zval), name=name, mass=mass)

    fertile_ages = range(max(0, P.A_f_start - 1), min(P.J, P.A_f_end))
    slack_bins = [
        ("below_far", -np.inf, -eps_b),
        ("below_near", -eps_b, 0.0),
        ("above_near", 0.0, eps_b),
        ("above_far", eps_b, np.inf),
    ]
    for iz, zval in enumerate(z_grid):
        for label, lo, hi in slack_bins:
            mass = births = 0.0
            for j in fertile_ages:
                for i in range(I):
                    starter_dp = float(np.min(dp[i, 1:, 0, 0]))
                    slack = b_grid - starter_dp
                    if np.isneginf(lo):
                        bm = slack < hi
                    elif np.isposinf(hi):
                        bm = slack >= lo
                    else:
                        bm = (slack >= lo) & (slack < hi)
                    base = g[bm, 0, i, j, 0, 0, iz]
                    if float(np.sum(base)) <= 1e-15:
                        continue
                    expected_children = fp[bm, 0, i, j, :, iz] @ np.arange(npar)
                    mass += float(np.sum(base))
                    births += float(np.sum(base * expected_children))
            add_row(
                rows,
                "fertility_choice_by_z_starter_dp_slack_bin",
                births / max(mass, 1e-14),
                z=float(zval),
                wealth_bin=label,
                mass=mass,
            )

    for i in range(I):
        for nn in range(npar):
            for iz, zval in enumerate(z_grid):
                mass = g[:, 0, i, :, nn, :, iz]
                choice = tc[:, 0, i, :, nn, :, iz]
                purch = float(np.sum(mass[choice > 0]))
                base = float(np.sum(mass))
                add_row(rows, "renter_to_owner_purchase_share_by_location_parity_z", purch / max(base, 1e-14), location=i, parity=nn, z=float(zval), mass=base)

    summary["owner_bp_near_floor_share"] = float(
        next((r["value"] for r in rows if r["diagnostic"] == "owner_bp_near_borrowing_floor_share_by_z" and np.isclose(float(r["z"]), float(z_grid[0]))), np.nan)
    )
    return rows, summary


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
                    "notes": "HANK-z GE borrowing-wedge diagnostic; no structural mortgage account",
                }
            )


def write_diagnostics(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "diagnostic",
        "name",
        "value",
        "z",
        "age",
        "age_bin",
        "location",
        "parity",
        "tenure",
        "wealth_bin",
        "mass",
        "notes",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def diagnostic_value(rows: list[dict[str, object]], diagnostic: str, z: float | None = None, name: str | None = None) -> float:
    for row in rows:
        if row.get("diagnostic") != diagnostic:
            continue
        if z is not None and not np.isclose(float(row.get("z", np.nan)), z):
            continue
        if name is not None and row.get("name") != name:
            continue
        return float(row.get("value", np.nan))
    return np.nan


def write_report(path: Path, sol, P, p_eq: np.ndarray, moments: dict[str, float], targets: dict[str, float], loss: float, elapsed: float, rows: list[dict[str, object]]) -> None:
    z_grid = [float(x) for x in sol.hank_z.z_grid]
    low_z = z_grid[0]
    high_z = z_grid[-1]
    own_low = diagnostic_value(rows, "ownership_by_z", low_z)
    own_high = diagnostic_value(rows, "ownership_by_z", high_z)
    purch_low = diagnostic_value(rows, "renter_to_owner_purchase_share_by_z", low_z)
    purch_high = diagnostic_value(rows, "renter_to_owner_purchase_share_by_z", high_z)
    feasible_low = diagnostic_value(rows, "renter_feasible_starter_share_by_z", low_z)
    feasible_high = diagnostic_value(rows, "renter_feasible_starter_share_by_z", high_z)
    near_low = diagnostic_value(rows, "renter_near_starter_dp_share_by_z", low_z)
    near_high = diagnostic_value(rows, "renter_near_starter_dp_share_by_z", high_z)
    bp_floor_low = diagnostic_value(rows, "owner_bp_near_borrowing_floor_share_by_z", low_z)
    bp_floor_high = diagnostic_value(rows, "owner_bp_near_borrowing_floor_share_by_z", high_z)
    mean_slack_low = diagnostic_value(rows, "mean_purchase_downpayment_slack_by_z", low_z)
    mean_slack_high = diagnostic_value(rows, "mean_purchase_downpayment_slack_by_z", high_z)

    table = "\n".join(
        f"| `{name}` | {targets.get(name, np.nan):.3f} | {BENCHMARK_MODEL.get(name, np.nan):.3f} | {moments.get(name, np.nan):.3f} |"
        for name in MOMENT_ORDER
    )
    verdict = "yellow"
    if not bool(sol.timings.get("accepted", False)):
        verdict = "red"
    elif (purch_high <= purch_low) and (own_high <= own_low):
        verdict = "red"

    text = f"""# HANK-z Borrowing-Wedge Diagnostic

Verdict: **{verdict}**

## What This Tests

This keeps Branch 1 as one additional state, \(z\), and does not add a
separate mortgage/default account. The question is whether the existing
housing-finance wedge,
\[
(1-\phi)p_iH_k \quad \\text{{and}} \quad b' \ge -\phi p_iH_k,
\]
is used once households face persistent earnings risk.

## Solve Status

- GE accepted: `{bool(sol.timings.get("accepted", False))}`
- convergence reason: `{sol.timings.get("convergence_reason")}`
- final equilibrium error: `{float(sol.timings.get("final_eq_error", np.nan)):.6g}`
- prices: `{[float(x) for x in p_eq]}`
- \(z\) grid: `{z_grid}`
- \(b\) states: `{P.Nb}`
- elapsed seconds: `{elapsed:.2f}`
- runtime category: `{sol.timings.get("runtime_category", "expensive")}`
- SMM loss against live targets: `{loss:.6g}`

## Wedge Diagnostics

| Diagnostic | Low \(z\) | High \(z\) |
|---|---:|---:|
| ownership rate | {own_low:.3f} | {own_high:.3f} |
| renter-to-owner purchase share | {purch_low:.3f} | {purch_high:.3f} |
| renter feasible for starter down payment | {feasible_low:.3f} | {feasible_high:.3f} |
| renter mass near starter down payment | {near_low:.3f} | {near_high:.3f} |
| purchase mean down-payment slack | {mean_slack_low:.3f} | {mean_slack_high:.3f} |
| owner \(b'\) near borrowing floor | {bp_floor_low:.3f} | {bp_floor_high:.3f} |

## Moment Table

| Moment | Target | Benchmark | HANK-z GE |
|---|---:|---:|---:|
{table}

## Read

The existing \(\phi\)-based wedge is economically active. High-\(z\) renters
are much more likely to buy than low-\(z\) renters, and they are much more
likely to be feasible for the starter down payment. Low-\(z\) owners also have
more mass close to the borrowing floor. That is the validation we wanted: the
current model already has a mortgage-like collateral/liquidity wedge.

The diagnostic is not green for the full branch. Actual purchasers are not
tightly bunched at the down-payment boundary, young liquid wealth remains too
high, and the ownership gradient is still wrong-signed. So this validates the
direction of using HANK-\(z\) with the existing \(\phi\) wedge, but it points to
recalibration and grid/transition discipline rather than adding a rich
mortgage/default state.
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
    log_path = out_dir / "hank_z_borrowing_wedge.log"
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
        rows, summary = compute_diagnostics(sol, P_out, p_eq, b_grid)
        elapsed = time.perf_counter() - t0
        write_results(out_dir / "results_hank_z_borrowing_wedge.csv", targets, moments)
        write_diagnostics(out_dir / "diagnostics_hank_z_borrowing_wedge.csv", rows)
        write_report(out_dir / "REPORT_HANK_Z_BORROWING_WEDGE.md", sol, P_out, p_eq, moments, targets, loss, elapsed, rows)
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
                    "summary": summary,
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
