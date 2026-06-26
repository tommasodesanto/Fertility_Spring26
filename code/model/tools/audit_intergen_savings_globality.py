#!/usr/bin/env python3
"""Dense-grid globality probe for intergen branch-level savings choices.

This is a read-only diagnostic. It loads a saved `solution_cache.pkl`, samples
reachable states, reconstructs each target tenure branch's one-dimensional
savings objective, and compares the stored policy to a dense global sweep.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
TOOLS_ROOT = Path(__file__).resolve().parent
for path in (MODEL_ROOT, TOOLS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from audit_intergen_active_tenure_values import (  # noqa: E402
    branch_info,
    interp_scalar,
    load_cache,
    owner_asset_price_vector,
    rent_vector,
    target_probabilities,
    tenure_label,
    transaction_arrays,
)
from intergen_housing_fertility.calibration import jsonable  # noqa: E402
from intergen_housing_fertility.solver import (  # noqa: E402
    apply_child_aging,
    bequest_utility_vec,
    eval_owner,
    eval_renter,
    get_completed_fertility,
    income_at_state,
    income_transition_values,
    precompute_shared,
)
from intergen_housing_fertility.utils import flat_nc  # noqa: E402


DEFAULT_CACHE = (
    ROOT
    / "output/model/intergen_fixedstats_overnight_review_20260626/"
    / "best_de_g044_i022_packet/solution_cache.pkl"
)
DEFAULT_OUTDIR = ROOT / "output/model/intergen_solver_accuracy_20260626/savings_globality"
NEG_INF = -1.0e10


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    payload = load_cache(args.cache.resolve())
    base = payload["baseline"]
    sol = base["sol"]
    P = base["P"]
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    price = owner_asset_price_vector(sol, base)
    rent = rent_vector(sol, P, price)
    shared = precompute_shared(P, b_grid)
    hcost, heq, dp_arr, bmo = transaction_arrays(P, price, b_grid, shared)

    states = select_probe_states(sol, P, b_grid, hcost, heq, dp_arr, bmo, shared, args)
    rows: list[dict[str, Any]] = []
    continuation_cache: dict[tuple[int, int, int], np.ndarray] = {}
    for sid, state in enumerate(states, start=1):
        bb, to, i, j, zz, nn, cs = state["index"]
        key = (int(i), int(j), int(zz))
        if key not in continuation_cache:
            continuation_cache[key] = continuation_values(sol, P, b_grid, int(i), int(j), int(zz), price)
        Vc = continuation_cache[key]
        probs = target_probabilities(sol, bb, to, i, j, zz, nn, cs)
        for tn in range(1 + int(P.n_house)):
            row = evaluate_branch_globality(
                sol=sol,
                P=P,
                b_grid=b_grid,
                price=price,
                rent=rent,
                shared=shared,
                hcost=hcost,
                heq=heq,
                dp_arr=dp_arr,
                bmo=bmo,
                Vc=Vc,
                state=state,
                state_id=sid,
                target_tenure=tn,
                target_probability=float(probs[tn]) if tn < probs.size else math.nan,
                n_dense=int(args.n_dense),
                tol=float(args.tol),
            )
            rows.append(row)

    summary_rows = summarize_rows(rows, tol=float(args.tol))
    failures = [
        row for row in rows
        if row.get("feasible") and float(row.get("value_gap_dense_minus_stored", math.nan)) > float(args.tol)
    ]
    relevant_failures = [
        row for row in failures
        if float(row.get("branch_weight", 0.0)) > float(args.relevance_tol)
    ]
    write_csv(outdir / "savings_globality_branch_rows.csv", rows)
    write_csv(outdir / "savings_globality_summary.csv", summary_rows)
    write_csv(outdir / "savings_globality_failures.csv", failures)
    write_csv(outdir / "savings_globality_relevant_failures.csv", relevant_failures)
    write_json(
        outdir / "savings_globality_metadata.json",
        {
            "cache": str(args.cache.resolve()),
            "outdir": str(outdir),
            "n_states": len(states),
            "n_rows": len(rows),
            "n_failures": len(failures),
            "n_relevant_failures": len(relevant_failures),
            "tol": float(args.tol),
            "relevance_tol": float(args.relevance_tol),
            "n_dense": int(args.n_dense),
            "grid": {"Nb": int(b_grid.size), "b_min": float(b_grid[0]), "b_max": float(b_grid[-1])},
            "H_own": np.asarray(P.H_own, dtype=float).tolist(),
            "hR_max": float(P.hR_max),
            "interp_method": str(getattr(P, "interp_method", "linear")),
            "tenure_choice_kappa": float(getattr(P, "tenure_choice_kappa", math.nan)),
        },
    )
    write_readme(outdir, args, rows, summary_rows, failures, relevant_failures)
    print(outdir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--top-states", type=int, default=80)
    parser.add_argument("--top-per-wealth", type=int, default=8)
    parser.add_argument("--n-dense", type=int, default=1601)
    parser.add_argument("--tol", type=float, default=5e-3)
    parser.add_argument("--relevance-tol", type=float, default=1e-10)
    parser.add_argument(
        "--wealth-targets",
        default="",
        help="Optional comma-separated wealth points; automatic mass quantiles and thresholds are always added.",
    )
    return parser.parse_args()


def select_probe_states(
    sol: Any,
    P: Any,
    b_grid: np.ndarray,
    hcost: np.ndarray,
    heq: np.ndarray,
    dp_arr: np.ndarray,
    bmo: np.ndarray,
    shared: Any,
    args: argparse.Namespace,
) -> list[dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    flat = g.reshape(-1)
    positive = np.flatnonzero(flat > 1e-12)
    order = positive[np.argsort(flat[positive])[::-1]]
    selected: dict[tuple[int, ...], set[str]] = {}

    def add(index: tuple[int, ...], reason: str) -> None:
        selected.setdefault(tuple(int(x) for x in index), set()).add(reason)

    for flat_idx in order[: int(args.top_states)]:
        add(np.unravel_index(int(flat_idx), g.shape), "top_mass")

    wealth_targets = automatic_wealth_targets(sol, P, b_grid, hcost, heq, dp_arr, bmo, args)
    for target in wealth_targets:
        bb = int(np.argmin(np.abs(b_grid - float(target))))
        for bbi in range(max(0, bb - 1), min(b_grid.size, bb + 2)):
            slice_mass = g[bbi, ...].reshape(-1)
            pos = np.flatnonzero(slice_mass > 1e-12)
            if pos.size == 0:
                continue
            local_order = pos[np.argsort(slice_mass[pos])[::-1]][: int(args.top_per_wealth)]
            for local in local_order:
                rest = np.unravel_index(int(local), g.shape[1:])
                add((bbi, *rest), f"wealth_near_{target:.4g}")

    # Explicit coverage for low/mid/high income, child states, terminal age, and
    # old/fertile ages using nearest positive-mass b cells.
    age_indices = sorted(set([0, max(0, int(P.A_f_start) - 1), max(0, int(P.A_f_end) - 1), max(0, int(P.J_R) - 1), int(P.J) - 1]))
    z_indices = sorted(set([0, max(0, len(getattr(P, "z_grid", [1])) // 2), max(0, len(getattr(P, "z_grid", [1])) - 1)]))
    for j in age_indices:
        for zz in z_indices:
            sub = g[:, :, :, j, zz, :, :]
            pos = np.argwhere(sub > 1e-12)
            if pos.size == 0:
                continue
            masses = np.array([sub[tuple(idx)] for idx in pos], dtype=float)
            for idx in pos[np.argsort(masses)[::-1]][: max(3, int(args.top_per_wealth) // 2)]:
                bb, to, i, nn, cs = [int(x) for x in idx]
                add((bb, to, i, j, zz, nn, cs), "age_z_bucket")

    states = []
    for index, reasons in selected.items():
        states.append(
            {
                "index": index,
                "selection_reason": ";".join(sorted(reasons)),
                "mass": float(g[index]),
            }
        )
    states.sort(key=lambda s: (-float(s["mass"]), s["index"]))
    return states


def automatic_wealth_targets(
    sol: Any,
    P: Any,
    b_grid: np.ndarray,
    hcost: np.ndarray,
    heq: np.ndarray,
    dp_arr: np.ndarray,
    bmo: np.ndarray,
    args: argparse.Namespace,
) -> list[float]:
    out = [0.0, float(getattr(P, "b_entry_fixed", 0.0))]
    for item in str(args.wealth_targets).split(","):
        item = item.strip()
        if item:
            out.append(float(item))
    g = np.asarray(sol.g, dtype=float)
    mass_b = np.sum(g, axis=tuple(range(1, g.ndim)))
    if np.sum(mass_b) > 0:
        cdf = np.cumsum(mass_b) / np.sum(mass_b)
        for q in (0.001, 0.01, 0.10, 0.50, 0.90, 0.99):
            out.append(float(b_grid[min(np.searchsorted(cdf, q), b_grid.size - 1)]))
    for i in range(int(P.I)):
        for to in range(1 + int(P.n_house)):
            sale = float(heq[i, to]) if to > 0 else 0.0
            if to > 0:
                out.append(-sale)
            for tn in range(1, 1 + int(P.n_house)):
                for nn in range(int(P.n_parity)):
                    for cs in range(int(P.n_child_states)):
                        out.append(float(dp_arr[i, tn, nn, cs]))
                        out.append(float(dp_arr[i, tn, nn, cs] - sale))
                        out.append(float(bmo[i, tn, nn, cs]))
    return sorted(set(float(x) for x in out if math.isfinite(float(x)) and b_grid[0] <= float(x) <= b_grid[-1]))


def continuation_values(sol: Any, P: Any, b_grid: np.ndarray, i: int, j: int, zz: int, price: np.ndarray) -> np.ndarray:
    V = np.asarray(sol.V, dtype=float)
    Nb = b_grid.size
    nt = 1 + int(P.n_house)
    npar = int(P.n_parity)
    ncs = int(P.n_child_states)
    z_grid, _, Pi_z = income_transition_values(P)
    if j == int(P.J) - 1:
        Vnr = np.zeros((Nb, nt, 1, npar, ncs))
        p = float(np.asarray(price, dtype=float).reshape(-1)[i])
        for ten in range(nt):
            hv = p * float(P.H_own[ten - 1]) if ten > 0 else 0.0
            for nn in range(npar):
                for cs in range(ncs):
                    nk = get_completed_fertility(nn, cs, P)
                    Vnr[:, ten, 0, nn, cs] = bequest_utility_vec(b_grid + hv, nk, P)
    else:
        Vnr = np.zeros((Nb, nt, 1, npar, ncs))
        for znext in range(len(z_grid)):
            Vnr[:, :, 0, :, :] += float(Pi_z[zz, znext]) * V[:, :, i, j + 1, znext, :, :]
    return apply_child_aging(Vnr, P, Nb, nt, 1, npar, ncs)[:, :, 0, :, :]


def evaluate_branch_globality(
    *,
    sol: Any,
    P: Any,
    b_grid: np.ndarray,
    price: np.ndarray,
    rent: np.ndarray,
    shared: Any,
    hcost: np.ndarray,
    heq: np.ndarray,
    dp_arr: np.ndarray,
    bmo: np.ndarray,
    Vc: np.ndarray,
    state: dict[str, Any],
    state_id: int,
    target_tenure: int,
    target_probability: float,
    n_dense: int,
    tol: float,
) -> dict[str, Any]:
    bb, to, i, j, zz, nn, cs = state["index"]
    b = float(b_grid[bb])
    info = branch_info(
        b=b,
        origin_tenure=to,
        target_tenure=target_tenure,
        market=i,
        parity=nn,
        child_state=cs,
        P=P,
        hcost=hcost,
        heq=heq,
        dp_arr=dp_arr,
        bmo=bmo,
        shared=shared,
    )
    base_row = {
        "state_id": int(state_id),
        "selection_reason": state["selection_reason"],
        "age_index": int(j),
        "age": float(getattr(P, "age_start", 22.0) + j * getattr(P, "period_years", 4.0)),
        "z_index": int(zz),
        "z": float(np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)[zz]),
        "parity": int(nn),
        "child_state": int(cs),
        "origin_tenure": int(to),
        "origin_label": tenure_label(P, to),
        "target_tenure": int(target_tenure),
        "target_label": tenure_label(P, target_tenure),
        "state_mass": float(state["mass"]),
        "target_probability": float(target_probability),
        "branch_weight": float(state["mass"]) * max(float(target_probability), 0.0),
        "current_b": b,
        "branch_b": float(info["branch_liquid_wealth"]),
        "post_tenure_total_wealth": post_tenure_total_wealth(P, price, i, target_tenure, float(info["branch_liquid_wealth"])),
        "feasible": bool(info["feasible"]),
        "feasibility_reason": str(info["feasibility_reason"]),
        "required_current_liquid": float(info.get("required_current_liquid", math.nan)),
        "downpayment": float(info.get("downpayment", math.nan)),
        "borrowing_floor": float(info.get("borrowing_floor", math.nan)),
    }
    if not bool(info["feasible"]):
        return base_row

    bp_line = np.asarray(sol.bp_pol)[:, target_tenure, i, j, zz, nn, cs]
    bp_stored = interp_scalar(b_grid, bp_line, float(info["branch_liquid_wealth"]))
    bp_lo, bp_hi, dense_points, eval_args = branch_dense_problem(
        P=P,
        b_grid=b_grid,
        price=price,
        rent=rent,
        shared=shared,
        Vc=Vc,
        i=i,
        j=j,
        zz=zz,
        nn=nn,
        cs=cs,
        branch_tenure=target_tenure,
        branch_b=float(info["branch_liquid_wealth"]),
        branch_borrowing_floor=float(info.get("borrowing_floor", math.nan)),
        bp_stored=float(bp_stored),
        n_dense=n_dense,
    )
    if dense_points.size == 0:
        return {**base_row, "bp_lo": bp_lo, "bp_hi": bp_hi, "status": "empty_feasible_interval"}
    values = evaluate_branch_values(target_tenure, dense_points, eval_args)
    stored_value = float(evaluate_branch_values(target_tenure, np.array([bp_stored], dtype=float), eval_args)[0])
    if np.all(~np.isfinite(values)):
        return {**base_row, "bp_lo": bp_lo, "bp_hi": bp_hi, "status": "all_dense_values_nonfinite"}
    best_idx = int(np.nanargmax(values))
    best_bp = float(dense_points[best_idx])
    best_value = float(values[best_idx])
    gap = best_value - stored_value
    return {
        **base_row,
        "bp_lo": float(bp_lo),
        "bp_hi": float(bp_hi),
        "bp_stored": float(bp_stored),
        "value_stored": stored_value,
        "bp_dense_best": best_bp,
        "value_dense_best": best_value,
        "value_gap_dense_minus_stored": float(gap),
        "bp_gap": best_bp - float(bp_stored),
        "n_dense": int(n_dense),
        "n_augmented_points": int(dense_points.size),
        "near_lower_bound": bool(abs(float(bp_stored) - float(bp_lo)) < 1e-5),
        "near_upper_bound": bool(abs(float(bp_stored) - float(bp_hi)) < 1e-5),
        "is_global_within_tol": bool(gap <= tol),
        "status": "ok",
    }


def branch_dense_problem(
    *,
    P: Any,
    b_grid: np.ndarray,
    price: np.ndarray,
    rent: np.ndarray,
    shared: Any,
    Vc: np.ndarray,
    i: int,
    j: int,
    zz: int,
    nn: int,
    cs: int,
    branch_tenure: int,
    branch_b: float,
    branch_borrowing_floor: float,
    bp_stored: float,
    n_dense: int,
) -> tuple[float, float, np.ndarray, dict[str, Any]]:
    z = float(np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)[zz])
    yj = income_at_state(P, i, j, z)
    Rv = float(P.R_gross) * float(branch_b) + float(yj)
    beta = float(P.beta)
    alpha = float(P.alpha_cons)
    oms = 1.0 - float(P.sigma)
    cb = np.asarray(shared.c_bar, dtype=float)
    hb = np.asarray(shared.h_bar, dtype=float)
    psi_v = np.asarray(shared.psi_v, dtype=float)
    cb_c = float(cb[nn, cs])
    hb_c = float(hb[nn, cs])
    pc = float(psi_v[nn, cs])
    cidx = int(nn + int(P.n_parity) * cs)
    Vbar = flat_nc(Vc[:, branch_tenure, :, :], b_grid.size, int(P.n_parity) * int(P.n_child_states))[:, cidx]
    if branch_tenure == 0:
        ri = float(np.asarray(rent, dtype=float).reshape(-1)[i])
        dc = cb_c + ri * hb_c
        cc = ri * (float(P.hR_max) - hb_c) / max(1.0 - alpha, 1e-12)
        ht_cap_c = max(float(P.hR_max) - hb_c, 1e-10)
        Kr = (alpha**alpha * ((1.0 - alpha) / ri) ** (1.0 - alpha)) ** oms
        lo = max(0.0, float(b_grid[0]))
        hi = max(Rv - dc - 1e-6, lo)
        kink = Rv - dc - cc
        args = {
            "Rv": Rv,
            "Vbar": Vbar,
            "dc": dc,
            "pc": pc,
            "cc": cc,
            "cb_c": cb_c,
            "ri": ri,
            "hRmax": float(P.hR_max),
            "ht_cap_c": ht_cap_c,
            "Kr": Kr,
            "alpha": alpha,
            "oms": oms,
            "beta": beta,
            "kinks": [kink],
        }
    else:
        p = float(np.asarray(price, dtype=float).reshape(-1)[i])
        H = float(np.asarray(P.H_own, dtype=float).reshape(-1)[branch_tenure - 1])
        owner_size_cost = float(getattr(P, "owner_size_cost", 0.0))
        owner_size_cost_ref = float(getattr(P, "owner_size_cost_ref", 6.0))
        owner_size_cost_power = float(getattr(P, "owner_size_cost_power", 2.0))
        extra = owner_size_cost * p * max(H - owner_size_cost_ref, 0.0) ** owner_size_cost_power
        oc = (float(P.delta) + float(P.tau_H)) * p * H + extra
        owner_h_bar_scale = float(getattr(P, "owner_h_bar_scale", 1.0))
        owner_service_premium = max(float(getattr(P, "chi", 1.0)), 1e-8)
        ht_c = owner_service_premium * max(H - owner_h_bar_scale * hb_c, 1e-10)
        Ko_c = ht_c ** ((1.0 - alpha) * oms)
        lo = max(float(b_grid[0]), float(branch_borrowing_floor))
        hi = max(Rv - oc - cb_c - 1e-6, lo)
        args = {
            "Rv": Rv,
            "Vbar": Vbar,
            "oc": oc,
            "cb_c": cb_c,
            "pc": pc,
            "Ko_c": Ko_c,
            "alpha": alpha,
            "oms": oms,
            "beta": beta,
            "kinks": [],
        }
    dense = augmented_dense_grid(b_grid, lo, hi, bp_stored, int(n_dense), args["kinks"])
    args["b_grid"] = b_grid
    return float(lo), float(hi), dense, args


def augmented_dense_grid(b_grid: np.ndarray, lo: float, hi: float, stored: float, n_dense: int, kinks: list[float]) -> np.ndarray:
    if not math.isfinite(lo) or not math.isfinite(hi) or hi < lo:
        return np.array([], dtype=float)
    base = np.linspace(lo, hi, max(2, int(n_dense)))
    knots = b_grid[(b_grid >= lo) & (b_grid <= hi)]
    extras = [lo, hi, stored, *kinks]
    vals = np.concatenate([base, knots, np.asarray([x for x in extras if math.isfinite(float(x))], dtype=float)])
    vals = vals[(vals >= lo) & (vals <= hi)]
    return np.unique(vals.astype(float))


def evaluate_branch_values(branch_tenure: int, bp: np.ndarray, args: dict[str, Any]) -> np.ndarray:
    Rv = np.full_like(bp, float(args["Rv"]), dtype=float)
    if branch_tenure == 0:
        return eval_renter(
            bp,
            Rv,
            args["Vbar"],
            args["b_grid"],
            args["dc"],
            args["pc"],
            args["cc"],
            args["cb_c"],
            args["ri"],
            args["hRmax"],
            args["ht_cap_c"],
            args["Kr"],
            args["alpha"],
            args["oms"],
            args["beta"],
        )
    return eval_owner(
        bp,
        Rv,
        args["Vbar"],
        args["b_grid"],
        args["oc"],
        args["cb_c"],
        args["pc"],
        args["Ko_c"],
        args["alpha"],
        args["oms"],
        args["beta"],
    )


def summarize_rows(rows: list[dict[str, Any]], *, tol: float) -> list[dict[str, Any]]:
    out = []
    feasible = [row for row in rows if row.get("feasible") and math.isfinite(float(row.get("value_gap_dense_minus_stored", math.nan)))]
    groups = {"all": feasible}
    for row in feasible:
        groups.setdefault(str(row["target_label"]), []).append(row)
    for group, group_rows in groups.items():
        gaps = np.asarray([float(row["value_gap_dense_minus_stored"]) for row in group_rows], dtype=float)
        weights = np.asarray([max(float(row.get("branch_weight", 0.0)), 0.0) for row in group_rows], dtype=float)
        out.append(
            {
                "group": group,
                "n_rows": int(len(group_rows)),
                "n_failures": int(np.sum(gaps > tol)) if gaps.size else 0,
                "n_relevant_failures": int(
                    np.sum((gaps > tol) & (weights > 1e-10))
                ) if gaps.size else 0,
                "max_gap": float(np.max(gaps)) if gaps.size else math.nan,
                "unweighted_p95_gap": float(np.quantile(gaps, 0.95)) if gaps.size else math.nan,
                "weighted_median_gap": weighted_quantile(gaps, weights, 0.50),
                "weighted_p95_gap": weighted_quantile(gaps, weights, 0.95),
                "total_branch_weight": float(np.sum(weights)),
            }
        )
    return out


def post_tenure_total_wealth(P: Any, price: np.ndarray, market: int, tenure: int, branch_b: float) -> float:
    if int(tenure) <= 0:
        return float(branch_b)
    p = float(np.asarray(price, dtype=float).reshape(-1)[market])
    H = float(np.asarray(P.H_own, dtype=float).reshape(-1)[int(tenure) - 1])
    return float(branch_b + (1.0 - float(P.psi)) * p * H)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(ok):
        return math.nan
    values = values[ok]
    weights = weights[ok]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights) / np.sum(weights)
    return float(values[min(np.searchsorted(cdf, q), values.size - 1)])


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True))


def write_readme(
    outdir: Path,
    args: argparse.Namespace,
    rows: list[dict[str, Any]],
    summary_rows: list[dict[str, Any]],
    failures: list[dict[str, Any]],
    relevant_failures: list[dict[str, Any]],
) -> None:
    lines = [
        "# Savings Globality Probe",
        "",
        "Read-only dense-grid check of branch-level savings choices.",
        "",
        f"- rows: `{len(rows)}`",
        f"- failures above tolerance `{float(args.tol):.3g}`: `{len(failures)}`",
        f"- KFE-relevant failures above tolerance: `{len(relevant_failures)}`",
        f"- dense points per interval: `{int(args.n_dense)}` before knot augmentation",
        "",
        "## Summary",
        "",
    ]
    for row in summary_rows:
        lines.append(
            f"- `{row['group']}`: rows `{row['n_rows']}`, failures `{row['n_failures']}`, "
            f"weighted p95 gap `{row['weighted_p95_gap']}`"
        )
    lines += [
        "",
        "## Files",
        "",
        "- `savings_globality_branch_rows.csv`: one row per sampled state x target branch.",
        "- `savings_globality_failures.csv`: subset with dense value gap above tolerance.",
        "- `savings_globality_relevant_failures.csv`: failures with positive branch relevance weight.",
        "- `savings_globality_summary.csv`: weighted and unweighted gap summaries.",
    ]
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
