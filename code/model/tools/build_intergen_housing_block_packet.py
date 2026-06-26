#!/usr/bin/env python3
"""Build a read-only housing-block diagnostic packet from a saved intergen cache.

The packet focuses on housing/wealth behavior before any further fertility
calibration. It reads an existing ``solution_cache.pkl`` and writes mass-weighted
policy objects using the same active-tenure convention as
``build_intergen_mechanics_packet.py``: policies are evaluated after the
current-state tenure transition at branch liquid wealth.
"""

from __future__ import annotations

import argparse
import csv
import math
import pickle
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_housing_fertility.solver import precompute_shared  # noqa: E402


DEFAULT_CACHE = (
    ROOT
    / "output/model/intergen_fixedstats_overnight_review_20260626/"
    / "best_de_g044_i022_packet/solution_cache.pkl"
)


def main() -> None:
    args = parse_args()
    cache_path = args.cache.resolve()
    outdir = args.outdir.resolve() if args.outdir else cache_path.parent / "housing_block_audit"
    outdir.mkdir(parents=True, exist_ok=True)

    payload = load_cache(cache_path)
    base = payload["baseline"]
    sol = base["sol"]
    P = base["P"]

    rows = active_state_rows(sol, P, state_mass_min=float(args.state_mass_min))
    liquid_rows = aggregate_rows(rows, key="current_liquid_wealth")
    total_rows = aggregate_rows(rows, key="post_tenure_total_wealth")
    current_total_rows = aggregate_rows(rows, key="current_total_wealth")
    owner_age_rows = owner_rung_shares_by_age(sol, P)
    renter_cap_rows = renter_cap_share_by_age(rows, P)
    old_transition_rows = old_owner_transition_rates(sol, P)

    write_csv(outdir / "aggregate_by_current_liquid_wealth.csv", liquid_rows)
    write_csv(outdir / "aggregate_by_post_tenure_total_wealth.csv", total_rows)
    write_csv(outdir / "aggregate_by_current_total_wealth.csv", current_total_rows)
    write_csv(outdir / "owner_rung_shares_by_age.csv", owner_age_rows)
    write_csv(outdir / "renter_cap_share_by_age.csv", renter_cap_rows)
    write_csv(outdir / "old_owner_transition_rates.csv", old_transition_rows)

    plot_consumption_line(
        liquid_rows,
        x="current_liquid_wealth",
        title="Average consumption by current liquid wealth",
        xlabel="current liquid wealth b",
        path=outdir / "average_consumption_by_current_liquid_wealth.png",
        mass_min=float(args.plot_mass_min),
    )
    plot_consumption_line(
        total_rows,
        x="post_tenure_total_wealth",
        title="Average consumption by post-tenure total wealth",
        xlabel="post-tenure total wealth",
        path=outdir / "average_consumption_by_post_tenure_total_wealth.png",
        mass_min=float(args.plot_mass_min),
    )
    plot_mass_lines(
        liquid_rows,
        total_rows,
        outdir / "stationary_mass_by_wealth.png",
        mass_min=float(args.plot_mass_min),
    )
    plot_owner_housing_lines(
        liquid_rows,
        total_rows,
        outdir / "ownership_and_housing_by_wealth.png",
        mass_min=float(args.plot_mass_min),
    )
    plot_next_liquid_wealth(
        liquid_rows,
        outdir / "next_liquid_wealth_by_current_liquid_wealth.png",
        mass_min=float(args.plot_mass_min),
    )
    plot_owner_rung_shares(owner_age_rows, outdir / "owner_rung_shares_by_age.png")
    plot_renter_cap_share(renter_cap_rows, outdir / "renter_cap_share_by_age.png")
    plot_old_transitions(old_transition_rows, outdir / "old_owner_transition_rates.png")
    write_readme(
        outdir,
        cache_path=cache_path,
        payload=payload,
        sol=sol,
        P=P,
        liquid_rows=liquid_rows,
        total_rows=total_rows,
        current_total_rows=current_total_rows,
        old_transition_rows=old_transition_rows,
        plot_mass_min=float(args.plot_mass_min),
    )
    print(outdir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE, help="Path to solution_cache.pkl.")
    parser.add_argument("--outdir", type=Path, default=None, help="Output folder.")
    parser.add_argument(
        "--state-mass-min",
        type=float,
        default=0.0,
        help="Ignore individual state cells with mass at or below this value.",
    )
    parser.add_argument(
        "--plot-mass-min",
        type=float,
        default=1e-8,
        help="Filter exact wealth points below this mass in plots only; CSVs keep all rows.",
    )
    return parser.parse_args()


def load_cache(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    if not isinstance(payload, dict) or "baseline" not in payload:
        raise ValueError(f"{path} does not look like an intergen solution cache")
    for key in ("sol", "P"):
        if key not in payload["baseline"]:
            raise ValueError(f"{path} baseline is missing {key!r}")
    return payload


def maybe_vector_value(values: np.ndarray, idx: int) -> float:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0:
        return math.nan
    return float(arr[min(max(int(idx), 0), arr.size - 1)])


def interp_policy_scalar(grid: np.ndarray, values: np.ndarray, query: float) -> float:
    x = np.asarray(grid, dtype=float).reshape(-1)
    y = np.asarray(values, dtype=float).reshape(-1)
    if x.size == 0 or y.size != x.size or not math.isfinite(query):
        return math.nan
    q = float(np.clip(query, x[0], x[-1]))
    return float(np.interp(q, x, y))


def owner_asset_price_vector(sol: Any, P: Any) -> np.ndarray:
    fallback = np.full(int(getattr(P, "I", 1)), math.nan)
    values = getattr(sol, "owner_asset_price", getattr(sol, "p_eq", fallback))
    arr = np.asarray(values, dtype=float).reshape(-1)
    return arr if arr.size else fallback


def liquidated_housing_value(P: Any, price: np.ndarray, tenure_index: int, market_index: int = 0) -> float:
    if int(tenure_index) <= 0:
        return 0.0
    h_own = np.asarray(getattr(P, "H_own", []), dtype=float).reshape(-1)
    if h_own.size == 0:
        return math.nan
    h = float(h_own[int(np.clip(int(tenure_index) - 1, 0, h_own.size - 1))])
    p = maybe_vector_value(price, market_index)
    if not math.isfinite(p):
        p = maybe_vector_value(price, 0)
    if not math.isfinite(p):
        return math.nan
    return float((1.0 - float(getattr(P, "psi", 0.0))) * p * h)


def deterministic_branch_wealth(
    *,
    b: float,
    origin_tenure: int,
    target_tenure: int,
    market: int,
    parity: int,
    child_state: int,
    P: Any,
    price: np.ndarray,
    shared: Any,
) -> float:
    to = int(origin_tenure)
    tn = int(target_tenure)
    if tn == to:
        return float(b)
    p = maybe_vector_value(price, int(market))
    if not math.isfinite(p):
        return math.nan
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    sale = 0.0 if to <= 0 else (1.0 - float(P.psi)) * p * float(h_own[to - 1])
    if tn <= 0:
        return float(max(b + sale, 0.0))
    hcost = p * float(h_own[tn - 1])
    if to <= 0:
        branch = float(b - hcost)
        if bool(shared.birth_dp[parity, child_state, to, tn]):
            phi = float(shared.phi_choice[market, tn, parity, child_state])
            branch = max(branch, -phi * hcost)
        else:
            grant = float(shared.birth_entry_grant[market, tn, parity, child_state])
            if grant > 0.0:
                branch += grant
        return branch
    return float(b + sale - hcost)


def target_probabilities(sol: Any, bb: int, to: int, i: int, j: int, zz: int, nn: int, cs: int) -> np.ndarray:
    nt = np.asarray(sol.g).shape[1]
    tp = getattr(sol, "tenure_probs", None)
    if tp is None:
        out = np.zeros(nt, dtype=float)
        out[int(sol.tenure_choice[bb, to, i, j, zz, nn, cs])] = 1.0
        return out
    return np.asarray(tp[bb, to, i, j, zz, nn, cs, :], dtype=float).reshape(-1)


def active_state_rows(sol: Any, P: Any, *, state_mass_min: float) -> list[dict[str, float]]:
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    c_pol = np.asarray(sol.c_pol, dtype=float)
    bp_pol = np.asarray(sol.bp_pol, dtype=float)
    hR_pol = np.asarray(sol.hR_pol, dtype=float)
    price = owner_asset_price_vector(sol, P)
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    shared = precompute_shared(P, b_grid)
    age_grid = np.asarray(float(P.age_start) + np.arange(int(P.J)) * float(P.da), dtype=float)

    rows: list[dict[str, float]] = []
    nz = np.argwhere(g > float(state_mass_min))
    for bb, to, i, j, zz, nn, cs in nz:
        mass = float(g[bb, to, i, j, zz, nn, cs])
        b = float(b_grid[bb])
        probs = target_probabilities(sol, int(bb), int(to), int(i), int(j), int(zz), int(nn), int(cs))
        if probs.size != g.shape[1] or float(np.sum(probs)) <= 0.0:
            continue
        current_total = b + liquidated_housing_value(P, price, int(to), int(i))
        for tn, pr_raw in enumerate(probs):
            pr = float(pr_raw)
            if pr <= 0.0:
                continue
            bw = deterministic_branch_wealth(
                b=b,
                origin_tenure=int(to),
                target_tenure=int(tn),
                market=int(i),
                parity=int(nn),
                child_state=int(cs),
                P=P,
                price=price,
                shared=shared,
            )
            c_val = interp_policy_scalar(b_grid, c_pol[:, tn, i, j, zz, nn, cs], bw)
            bp_val = interp_policy_scalar(b_grid, bp_pol[:, tn, i, j, zz, nn, cs], bw)
            if tn <= 0:
                h_val = interp_policy_scalar(b_grid, hR_pol[:, 0, i, j, zz, nn, cs], bw)
            else:
                h_val = float(h_own[tn - 1])
            post_total = bw + liquidated_housing_value(P, price, int(tn), int(i))
            rows.append(
                {
                    "mass": mass * pr,
                    "current_liquid_wealth": b,
                    "current_total_wealth": current_total,
                    "post_tenure_total_wealth": post_total,
                    "consumption": c_val,
                    "next_liquid_wealth": bp_val,
                    "housing_services": h_val,
                    "owner_probability": float(tn > 0),
                    "target_renter_probability": float(tn <= 0),
                    "target_renter_cap_probability": float(tn <= 0 and h_val >= float(P.hR_max) - 1e-8),
                    "age": float(age_grid[j]),
                }
            )
    return rows


def aggregate_rows(rows: list[dict[str, float]], *, key: str) -> list[dict[str, float]]:
    acc: dict[float, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in rows:
        x = float(row[key])
        if not math.isfinite(x):
            continue
        k = round(x, 10)
        mass = float(row["mass"])
        if mass <= 0.0:
            continue
        item = acc[k]
        item[key] += mass * x
        item["mass"] += mass
        for col in ("consumption", "next_liquid_wealth", "housing_services", "owner_probability"):
            val = float(row[col])
            if math.isfinite(val):
                item[col] += mass * val
    out: list[dict[str, float]] = []
    for k, item in acc.items():
        mass = item["mass"]
        out.append(
            {
                key: item[key] / mass,
                "mass": mass,
                "avg_consumption": item["consumption"] / mass,
                "avg_next_liquid_wealth": item["next_liquid_wealth"] / mass,
                "avg_housing_services": item["housing_services"] / mass,
                "owner_probability": item["owner_probability"] / mass,
            }
        )
    out.sort(key=lambda r: r[key])
    return out


def owner_rung_shares_by_age(sol: Any, P: Any) -> list[dict[str, float]]:
    g = np.asarray(sol.g, dtype=float)
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    age_grid = np.asarray(float(P.age_start) + np.arange(int(P.J)) * float(P.da), dtype=float)
    rows: list[dict[str, float]] = []
    for j, age in enumerate(age_grid):
        age_mass = float(np.sum(g[:, :, :, j, :, :, :]))
        owner_mass = float(np.sum(g[:, 1:, :, j, :, :, :]))
        row: dict[str, float] = {
            "age": float(age),
            "age_mass": age_mass,
            "owner_mass": owner_mass,
            "owner_rate": owner_mass / age_mass if age_mass > 0 else math.nan,
        }
        for idx, rooms in enumerate(h_own, start=1):
            m = float(np.sum(g[:, idx, :, j, :, :, :]))
            row[f"owner_H{rooms:g}_mass"] = m
            row[f"owner_H{rooms:g}_share_of_owners"] = m / owner_mass if owner_mass > 0 else math.nan
        rows.append(row)
    return rows


def renter_cap_share_by_age(rows: list[dict[str, float]], P: Any) -> list[dict[str, float]]:
    age_acc: dict[float, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for row in rows:
        renter_prob = float(row.get("target_renter_probability", 0.0))
        cap_prob = row.get("target_renter_cap_probability", math.nan)
        if not math.isfinite(cap_prob):
            continue
        age = float(row["age"])
        mass = float(row["mass"])
        age_acc[age]["target_renter_mass"] += mass * renter_prob
        age_acc[age]["target_renter_cap_mass"] += mass * float(cap_prob)
    out: list[dict[str, float]] = []
    for age in sorted(age_acc):
        d = age_acc[age]
        rm = d["target_renter_mass"]
        cm = d["target_renter_cap_mass"]
        out.append(
            {
                "age": age,
                "hR_max": float(P.hR_max),
                "target_renter_mass": rm,
                "target_renter_cap_mass": cm,
                "target_renter_cap_share": cm / rm if rm > 0 else math.nan,
            }
        )
    return out


def old_owner_transition_rates(sol: Any, P: Any) -> list[dict[str, float]]:
    g = np.asarray(sol.g, dtype=float)
    h_own = np.asarray(P.H_own, dtype=float).reshape(-1)
    age_grid = np.asarray(float(P.age_start) + np.arange(int(P.J)) * float(P.da), dtype=float)
    totals = defaultdict(float)
    for bb in range(g.shape[0]):
        for to in range(1, g.shape[1]):
            for i in range(g.shape[2]):
                for j, age in enumerate(age_grid):
                    if age < 65.0 or age > 75.0:
                        continue
                    for zz in range(g.shape[4]):
                        for nn in range(g.shape[5]):
                            for cs in range(g.shape[6]):
                                mass = float(g[bb, to, i, j, zz, nn, cs])
                                if mass <= 0.0:
                                    continue
                                probs = target_probabilities(sol, bb, to, i, j, zz, nn, cs)
                                for tn, pr_raw in enumerate(probs):
                                    pr = float(pr_raw)
                                    if pr <= 0.0:
                                        continue
                                    w = mass * pr
                                    totals["old_owner_mass"] += w
                                    if tn <= 0:
                                        totals["sell_to_rent"] += w
                                    elif tn == to:
                                        totals["same_rung"] += w
                                    elif h_own[tn - 1] < h_own[to - 1]:
                                        totals["downsize"] += w
                                    elif h_own[tn - 1] > h_own[to - 1]:
                                        totals["upsize"] += w
                                    else:
                                        totals["other_owner_change"] += w
    denom = totals["old_owner_mass"]
    rows = []
    for name in ("same_rung", "downsize", "upsize", "sell_to_rent", "other_owner_change"):
        rows.append({"transition": name, "mass": totals[name], "share": totals[name] / denom if denom > 0 else math.nan})
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def filtered_xy(rows: list[dict[str, float]], x: str, y: str, *, mass_min: float) -> tuple[np.ndarray, np.ndarray]:
    sel = [r for r in rows if float(r.get("mass", 0.0)) >= mass_min and math.isfinite(float(r.get(y, math.nan)))]
    return np.asarray([r[x] for r in sel], dtype=float), np.asarray([r[y] for r in sel], dtype=float)


def plot_consumption_line(rows: list[dict[str, float]], *, x: str, title: str, xlabel: str, path: Path, mass_min: float) -> None:
    xs, ys = filtered_xy(rows, x, "avg_consumption", mass_min=mass_min)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xs, ys, color="black", linewidth=2)
    ax.axvline(0.0, color="0.65", linewidth=1)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("average consumption")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_mass_lines(liquid_rows: list[dict[str, float]], total_rows: list[dict[str, float]], path: Path, *, mass_min: float) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, rows, x, title in [
        (axes[0], liquid_rows, "current_liquid_wealth", "Mass by current liquid wealth"),
        (axes[1], total_rows, "post_tenure_total_wealth", "Mass by post-tenure total wealth"),
    ]:
        xs, ys = filtered_xy(rows, x, "mass", mass_min=mass_min)
        ax.plot(xs, ys, color="black", linewidth=2)
        ax.axvline(0.0, color="0.65", linewidth=1)
        ax.set_title(title)
        ax.set_xlabel(x.replace("_", " "))
        ax.set_ylabel("stationary mass")
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_owner_housing_lines(liquid_rows: list[dict[str, float]], total_rows: list[dict[str, float]], path: Path, *, mass_min: float) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    specs = [
        (axes[0, 0], liquid_rows, "current_liquid_wealth", "owner_probability", "Owner probability by liquid wealth"),
        (axes[0, 1], total_rows, "post_tenure_total_wealth", "owner_probability", "Owner probability by total wealth"),
        (axes[1, 0], liquid_rows, "current_liquid_wealth", "avg_housing_services", "Housing services by liquid wealth"),
        (axes[1, 1], total_rows, "post_tenure_total_wealth", "avg_housing_services", "Housing services by total wealth"),
    ]
    for ax, rows, x, y, title in specs:
        xs, ys = filtered_xy(rows, x, y, mass_min=mass_min)
        ax.plot(xs, ys, color="black", linewidth=2)
        ax.axvline(0.0, color="0.65", linewidth=1)
        ax.set_title(title)
        ax.set_xlabel(x.replace("_", " "))
        ax.set_ylabel(y.replace("_", " "))
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_next_liquid_wealth(rows: list[dict[str, float]], path: Path, *, mass_min: float) -> None:
    xs, ys = filtered_xy(rows, "current_liquid_wealth", "avg_next_liquid_wealth", mass_min=mass_min)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xs, ys, color="black", linewidth=2, label="E[b' | b]")
    if xs.size:
        lo = min(float(np.nanmin(xs)), float(np.nanmin(ys)))
        hi = max(float(np.nanmax(xs)), float(np.nanmax(ys)))
        ax.plot([lo, hi], [lo, hi], color="0.55", linewidth=1, linestyle="--", label="45-degree")
    ax.axvline(0.0, color="0.65", linewidth=1)
    ax.set_title("Average next liquid wealth by current liquid wealth")
    ax.set_xlabel("current liquid wealth b")
    ax.set_ylabel("average next liquid wealth")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_owner_rung_shares(rows: list[dict[str, float]], path: Path) -> None:
    if not rows:
        return
    ages = np.asarray([r["age"] for r in rows], dtype=float)
    keys = [k for k in rows[0] if k.endswith("_share_of_owners")]
    fig, ax = plt.subplots(figsize=(8, 5))
    for key in keys:
        vals = np.asarray([r.get(key, math.nan) for r in rows], dtype=float)
        label = key.replace("owner_", "").replace("_share_of_owners", "")
        ax.plot(ages, vals, linewidth=2, label=label)
    ax.set_title("Owner rung shares by age")
    ax.set_xlabel("age")
    ax.set_ylabel("share of owners")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_renter_cap_share(rows: list[dict[str, float]], path: Path) -> None:
    ages = np.asarray([r["age"] for r in rows], dtype=float)
    vals = np.asarray([r["target_renter_cap_share"] for r in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ages, vals, color="black", linewidth=2)
    ax.set_title("Target-renter cap share by age")
    ax.set_xlabel("age")
    ax.set_ylabel("share at renter cap")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_old_transitions(rows: list[dict[str, float]], path: Path) -> None:
    labels = [r["transition"] for r in rows]
    vals = [r["share"] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(labels, vals, color="0.2")
    ax.set_title("Old owner transition shares, ages 65-75")
    ax.set_ylabel("share of old-owner transition mass")
    ax.set_ylim(0.0, max(1.0, max(vals) * 1.1 if vals else 1.0))
    ax.tick_params(axis="x", rotation=25)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def largest_negative_changes(rows: list[dict[str, float]], x: str, y: str, *, mass_min: float, limit: int = 8) -> list[tuple[float, float, float]]:
    sel = [r for r in rows if float(r.get("mass", 0.0)) >= mass_min]
    out: list[tuple[float, float, float]] = []
    for prev, cur in zip(sel, sel[1:]):
        dy = float(cur[y] - prev[y])
        if dy < 0:
            out.append((float(cur[x]), dy, float(cur["mass"])))
    return sorted(out, key=lambda t: t[1])[:limit]


def write_readme(
    outdir: Path,
    *,
    cache_path: Path,
    payload: dict[str, Any],
    sol: Any,
    P: Any,
    liquid_rows: list[dict[str, float]],
    total_rows: list[dict[str, float]],
    current_total_rows: list[dict[str, float]],
    old_transition_rows: list[dict[str, float]],
    plot_mass_min: float,
) -> None:
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    mass_by_b = np.sum(g, axis=tuple(range(1, g.ndim)))
    total_mass = float(np.sum(g))
    zero_idx = int(np.argmin(np.abs(b_grid)))
    neg_liq_share = float(np.sum(mass_by_b[b_grid < 0]) / max(total_mass, 1e-14))
    zero_mass = float(mass_by_b[zero_idx])
    aggregate_c = sum(r["mass"] * r["avg_consumption"] for r in liquid_rows)
    price = owner_asset_price_vector(sol, P)
    equity_shifts = [
        float((1.0 - float(P.psi)) * maybe_vector_value(price, 0) * h)
        for h in np.asarray(P.H_own, dtype=float).reshape(-1)
    ]
    liquid_drops = largest_negative_changes(
        liquid_rows,
        "current_liquid_wealth",
        "avg_consumption",
        mass_min=plot_mass_min,
    )
    total_drops = largest_negative_changes(
        total_rows,
        "post_tenure_total_wealth",
        "avg_consumption",
        mass_min=plot_mass_min,
    )
    with (outdir / "README.md").open("w", encoding="utf-8") as fh:
        fh.write("# Housing-Block Audit Packet\n\n")
        fh.write(f"- Source cache: `{cache_path}`\n")
        fh.write(f"- Target set in cache: `{payload.get('target_set', '')}`\n")
        fh.write("- Policy convention: active post-tenure policy, integrated over tenure logit probabilities when available.\n")
        fh.write("- Liquid wealth x-axis: current state liquid wealth `b`.\n")
        fh.write("- Total wealth x-axis: post-tenure total wealth `branch b + (1 - psi) pH` for owner outcomes, branch `b` for renter outcomes.\n")
        fh.write(f"- Plot point mass filter: `{plot_mass_min:g}`; CSV files keep all positive-mass exact points.\n")
        fh.write(f"- Total stationary mass: `{total_mass:.12g}`\n")
        fh.write(f"- Aggregate consumption from liquid aggregation: `{aggregate_c:.12g}`\n")
        fh.write(f"- Mass at grid node closest to zero (`b={b_grid[zero_idx]:.6g}`): `{zero_mass:.12g}`\n")
        fh.write(f"- Negative liquid-wealth mass share: `{neg_liq_share:.12g}`\n")
        fh.write(f"- Owner equity shifts `(1-psi)pH`: `{equity_shifts}`\n\n")
        fh.write("## Files\n\n")
        for name in [
            "average_consumption_by_current_liquid_wealth.png",
            "average_consumption_by_post_tenure_total_wealth.png",
            "stationary_mass_by_wealth.png",
            "ownership_and_housing_by_wealth.png",
            "next_liquid_wealth_by_current_liquid_wealth.png",
            "owner_rung_shares_by_age.png",
            "renter_cap_share_by_age.png",
            "old_owner_transition_rates.png",
        ]:
            fh.write(f"- `{name}`\n")
        fh.write("\n## Largest Negative Adjacent Changes In Average Consumption\n\n")
        fh.write("### By Current Liquid Wealth\n\n")
        fh.write("| b | delta avg c | mass at point |\n|---:|---:|---:|\n")
        for x, dy, mass in liquid_drops:
            fh.write(f"| {x:.6g} | {dy:.6g} | {mass:.6g} |\n")
        fh.write("\n### By Post-Tenure Total Wealth\n\n")
        fh.write("| W | delta avg c | mass at point |\n|---:|---:|---:|\n")
        for x, dy, mass in total_drops:
            fh.write(f"| {x:.6g} | {dy:.6g} | {mass:.6g} |\n")
        fh.write("\n## Old Owner Transition Shares\n\n")
        fh.write("| transition | share | mass |\n|---|---:|---:|\n")
        for row in old_transition_rows:
            fh.write(f"| {row['transition']} | {row['share']:.6g} | {row['mass']:.6g} |\n")
        fh.write("\n## Current Total Wealth CSV\n\n")
        fh.write(
            "`aggregate_by_current_total_wealth.csv` is included for comparison only. "
            "The main total-wealth plot uses post-tenure wealth because consumption is an active post-tenure policy object.\n"
        )


if __name__ == "__main__":
    main()
