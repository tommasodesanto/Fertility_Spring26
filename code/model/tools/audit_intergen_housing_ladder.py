#!/usr/bin/env python3
"""Audit owner housing ladder use from a saved intergen solution cache."""

from __future__ import annotations

import argparse
import csv
import math
import pickle
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RUN_DIR = ROOT / "output" / "model" / "intergen_model_run_current"
DEFAULT_CACHE = DEFAULT_RUN_DIR / "solution_cache.pkl"
DEFAULT_OUTDIR = DEFAULT_RUN_DIR / "ladder_audit"


def main() -> None:
    args = parse_args()
    cache_path = args.cache.resolve()
    run_dir = args.run_dir.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    payload = load_cache(cache_path)
    base = payload["baseline"]
    sol = base["sol"]
    P = base["P"]

    floor_rows = family_floor_rows(P)
    rung_rows = rung_mass_rows(sol, P)
    wealth_rows = downpayment_wealth_rows(sol, P)
    rep_rows = representative_choice_rows(sol, P)

    write_csv(outdir / "family_floor_by_rung.csv", floor_rows)
    write_csv(outdir / "rung_mass_by_group.csv", rung_rows)
    write_csv(outdir / "downpayment_thresholds_vs_wealth.csv", wealth_rows)
    write_csv(outdir / "representative_rung_choice_probabilities.csv", rep_rows)

    plot_family_floors(P, outdir / "family_floor_vs_rungs.png")
    plot_downpayment_thresholds(wealth_rows, outdir / "downpayment_thresholds_vs_wealth.png")
    plot_rung_usage(rung_rows, outdir / "housing_rung_usage_by_group.png")
    plot_representative_probabilities(rep_rows, outdir / "representative_rung_choice_probabilities.png")

    readme_path = outdir / "README.md"
    readme_path.write_text(build_audit_readme(payload, base, sol, P, outdir), encoding="utf-8")
    write_plot_index(run_dir, outdir, payload, base, sol, P)
    print(f"Wrote ladder audit to {outdir}")
    print(f"Updated plot index at {run_dir / 'PLOTS_INDEX.md'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE, help="Saved solution_cache.pkl path.")
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR, help="Current run plot folder.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR, help="Ladder audit output folder.")
    return parser.parse_args()


def load_cache(path: Path) -> dict[str, Any]:
    with path.open("rb") as fh:
        payload = pickle.load(fh)
    if not isinstance(payload, dict) or "baseline" not in payload:
        raise ValueError(f"{path} does not look like an intergen mechanics cache")
    base = payload["baseline"]
    for key in ("sol", "P"):
        if key not in base:
            raise ValueError(f"{path} baseline is missing {key!r}")
    return payload


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def ages(P: Any) -> np.ndarray:
    return np.asarray(float(P.age_start) + np.arange(int(P.J)) * float(P.da), dtype=float)


def child_floor(P: Any, nn: int, cs: int) -> float:
    k = int(getattr(P, "n_child_stages", 1))
    child_present = cs >= 1 and cs < k + 1
    if not child_present:
        return float(P.h_bar_0)
    if str(getattr(P, "child_housing_spec", "")).lower() == "linear_only":
        return float(P.h_bar_0) + float(P.h_bar_n) * nn
    return float(P.h_bar_0) + float(P.h_bar_jump) + float(P.h_bar_n) * nn


def completed_children(P: Any, nn: int, cs: int) -> int:
    k = int(getattr(P, "n_child_stages", 1))
    if cs == 0:
        return 0
    if cs == k + 1:
        return 1
    if cs == k + 2:
        return 2
    return int(nn)


def child_status(P: Any, nn: int, cs: int) -> str:
    k = int(getattr(P, "n_child_stages", 1))
    done = completed_children(P, nn, cs)
    if done == 0:
        return "completed_childless"
    if 1 <= cs <= k:
        return "dependent_children"
    return "empty_nest_parent"


def current_child_bin(P: Any, nn: int, cs: int) -> str:
    k = int(getattr(P, "n_child_stages", 1))
    if cs == 0 or cs > k:
        return "no_current_dependent_child"
    if nn <= 0:
        return "no_current_dependent_child"
    if nn == 1:
        return "one_current_dependent_child"
    return "two_plus_current_dependent_children"


def family_floor_rows(P: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    chi = float(getattr(P, "chi", 1.0))
    for nn in range(int(P.n_parity)):
        for cs in range(int(P.n_child_states)):
            hbar = child_floor(P, nn, cs)
            for rung, rooms in enumerate(np.asarray(P.H_own, dtype=float), start=1):
                slack = float(rooms - hbar)
                rows.append(
                    {
                        "parity_state_nn": nn,
                        "child_state_cs": cs,
                        "completed_children": completed_children(P, nn, cs),
                        "child_status": child_status(P, nn, cs),
                        "current_child_bin": current_child_bin(P, nn, cs),
                        "h_bar": hbar,
                        "rung": rung,
                        "rooms": float(rooms),
                        "rooms_minus_h_bar": slack,
                        "below_or_at_floor": bool(slack <= 1e-10),
                        "owner_service_after_floor": float(chi * max(slack, 1e-10)),
                    }
                )
    return rows


def rung_mass_rows(sol: Any, P: Any) -> list[dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    age_grid = ages(P)
    n_house = int(P.n_house)
    groups: list[tuple[str, Callable[[int, int, int, int], bool]]] = [
        ("all_owners", lambda j, z, nn, cs: True),
        ("young_25_34", lambda j, z, nn, cs: 25.0 <= age_grid[j] <= 34.0),
        ("prime_30_55", lambda j, z, nn, cs: 30.0 <= age_grid[j] <= 55.0),
        ("old_65_75", lambda j, z, nn, cs: 65.0 <= age_grid[j] <= 75.0),
        ("completed_childless", lambda j, z, nn, cs: completed_children(P, nn, cs) == 0),
        ("completed_parent", lambda j, z, nn, cs: completed_children(P, nn, cs) > 0),
        ("current_dependent_child", lambda j, z, nn, cs: current_child_bin(P, nn, cs) != "no_current_dependent_child"),
        (
            "empty_nest_parent",
            lambda j, z, nn, cs: completed_children(P, nn, cs) > 0
            and current_child_bin(P, nn, cs) == "no_current_dependent_child",
        ),
        ("low_income_state", lambda j, z, nn, cs: z == 0),
        ("high_income_state", lambda j, z, nn, cs: z == g.shape[4] - 1),
    ]
    rows: list[dict[str, Any]] = []
    pop_total = float(np.sum(g))
    for group, keep in groups:
        masses = np.zeros(n_house, dtype=float)
        for ten in range(1, n_house + 1):
            for j in range(g.shape[3]):
                for z in range(g.shape[4]):
                    for nn in range(g.shape[5]):
                        for cs in range(g.shape[6]):
                            if keep(j, z, nn, cs):
                                masses[ten - 1] += float(np.sum(g[:, ten, :, j, z, nn, cs]))
        owner_group_mass = float(np.sum(masses))
        for idx, rooms in enumerate(np.asarray(P.H_own, dtype=float)):
            mass = float(masses[idx])
            rows.append(
                {
                    "group": group,
                    "rung": idx + 1,
                    "rooms": float(rooms),
                    "mass": mass,
                    "share_of_group_owners": mass / owner_group_mass if owner_group_mass > 1e-14 else math.nan,
                    "population_share": mass / pop_total if pop_total > 1e-14 else math.nan,
                    "group_owner_mass": owner_group_mass,
                }
            )
    return rows


def phi_for_state(P: Any, rung: int, nn: int, cs: int) -> float:
    phi = np.asarray(getattr(P, "phi", [0.8]), dtype=float).reshape(-1)
    out = float(phi[min(nn, phi.size - 1)])
    if not bool(getattr(P, "parent_dp_waiver", False)):
        return out
    k = int(getattr(P, "n_child_stages", 1))
    birth_only = bool(getattr(P, "parent_dp_waiver_birth_state_only", False))
    eligible_child_state = (cs == 1) if birth_only else (1 <= cs <= k)
    rung_vals = np.asarray(getattr(P, "parent_dp_waiver_owner_rungs", []), dtype=int).reshape(-1)
    eligible_rung = rung in set(rung_vals.tolist()) if rung_vals.size else True
    if nn >= 1 and eligible_child_state and eligible_rung:
        out = max(out, float(getattr(P, "parent_dp_waiver_phi", out)))
    return out


def house_price(P: Any, sol: Any) -> float:
    p = np.asarray(getattr(sol, "p_eq", getattr(sol, "owner_asset_price", [np.nan])), dtype=float).reshape(-1)
    if p.size and np.isfinite(p[0]):
        return float(p[0])
    return float(np.asarray(getattr(sol, "owner_asset_price", [np.nan]), dtype=float).reshape(-1)[0])


def downpayment_wealth_rows(sol: Any, P: Any) -> list[dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    age_grid = ages(P)
    p = house_price(P, sol)
    group_specs: list[tuple[str, str, Callable[[int, int, int, int], bool]]] = [
        ("all_households", "all", lambda j, z, nn, cs: True),
        ("renters_all", "renter", lambda j, z, nn, cs: True),
        ("renters_young_25_34", "renter", lambda j, z, nn, cs: 25.0 <= age_grid[j] <= 34.0),
        (
            "renters_prime_childless_30_55",
            "renter",
            lambda j, z, nn, cs: 30.0 <= age_grid[j] <= 55.0 and completed_children(P, nn, cs) == 0,
        ),
        (
            "renters_prime_parents_30_55",
            "renter",
            lambda j, z, nn, cs: 30.0 <= age_grid[j] <= 55.0 and completed_children(P, nn, cs) > 0,
        ),
        ("households_old_65_75", "all", lambda j, z, nn, cs: 65.0 <= age_grid[j] <= 75.0),
        ("renters_low_income_state", "renter", lambda j, z, nn, cs: z == 0),
        ("renters_high_income_state", "renter", lambda j, z, nn, cs: z == g.shape[4] - 1),
    ]
    rows: list[dict[str, Any]] = []
    for group, tenure_scope, keep in group_specs:
        weights = wealth_weights(g, tenure_scope, keep)
        mass = float(np.sum(weights))
        q10, q25, q50, q75, q90 = weighted_quantiles(b_grid, weights, [0.10, 0.25, 0.50, 0.75, 0.90])
        support = b_grid[weights > 0]
        support_min = float(support.min()) if support.size else math.nan
        support_max = float(support.max()) if support.size else math.nan
        for rung, rooms in enumerate(np.asarray(P.H_own, dtype=float), start=1):
            # Report the baseline childless threshold. Family-state-specific phi
            # changes are visible in family_floor_by_rung.csv if policy waivers are active.
            phi = phi_for_state(P, rung, 0, 0)
            hcost = p * float(rooms)
            dp = (1.0 - phi) * hcost
            rows.append(
                {
                    "wealth_group": group,
                    "tenure_scope": tenure_scope,
                    "mass": mass,
                    "rung": rung,
                    "rooms": float(rooms),
                    "asset_price": p,
                    "financed_share_phi_childless": phi,
                    "house_cost_pH": hcost,
                    "down_payment_threshold": dp,
                    "borrowing_floor_after_purchase": -phi * hcost,
                    "p10_liquid_wealth": q10,
                    "p25_liquid_wealth": q25,
                    "p50_liquid_wealth": q50,
                    "p75_liquid_wealth": q75,
                    "p90_liquid_wealth": q90,
                    "occupied_support_min": support_min,
                    "occupied_support_max": support_max,
                    "share_with_liquid_wealth_above_dp": share_above(b_grid, weights, dp),
                }
            )
    return rows


def wealth_weights(
    g: np.ndarray,
    tenure_scope: str,
    keep: Callable[[int, int, int, int], bool],
) -> np.ndarray:
    weights = np.zeros(g.shape[0], dtype=float)
    tenures: range
    if tenure_scope == "renter":
        tenures = range(0, 1)
    elif tenure_scope == "owner":
        tenures = range(1, g.shape[1])
    else:
        tenures = range(g.shape[1])
    for ten in tenures:
        for j in range(g.shape[3]):
            for z in range(g.shape[4]):
                for nn in range(g.shape[5]):
                    for cs in range(g.shape[6]):
                        if keep(j, z, nn, cs):
                            weights += np.sum(g[:, ten, :, j, z, nn, cs], axis=1)
    return weights


def weighted_quantiles(values: np.ndarray, weights: np.ndarray, qs: list[float]) -> list[float]:
    mass = float(np.sum(weights))
    if mass <= 1e-14:
        return [math.nan for _ in qs]
    order = np.argsort(values)
    x = values[order]
    w = weights[order]
    cdf = np.cumsum(w) / mass
    return [float(x[np.searchsorted(cdf, q, side="left").clip(0, len(x) - 1)]) for q in qs]


def share_above(values: np.ndarray, weights: np.ndarray, threshold: float) -> float:
    mass = float(np.sum(weights))
    if mass <= 1e-14:
        return math.nan
    return float(np.sum(weights[values >= threshold]) / mass)


def representative_choice_rows(sol: Any, P: Any) -> list[dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    tp = np.asarray(sol.tenure_probs, dtype=float)
    choice = np.asarray(sol.tenure_choice)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    age_grid = ages(P)
    p = house_price(P, sol)
    kappa = float(getattr(P, "tenure_choice_kappa", math.nan))
    reps = [
        ("age30_childless_renter_low_z", 30.0, 0, 0, 0, 0),
        ("age30_childless_renter_high_z", 30.0, g.shape[4] - 1, 0, 0, 0),
        ("age42_one_child_renter_mid_z", 42.0, min(2, g.shape[4] - 1), 0, 1, 1),
        ("age42_two_child_renter_high_z", 42.0, g.shape[4] - 1, 0, 2, 1),
        ("age66_childless_owner_6rooms_mid_z", 66.0, min(2, g.shape[4] - 1), 3, 0, 0),
        ("age66_parent_owner_8rooms_mid_z", 66.0, min(2, g.shape[4] - 1), 4, 2, min(3, g.shape[6] - 1)),
    ]
    rows: list[dict[str, Any]] = []
    for label, target_age, z, to, nn, cs in reps:
        j = int(np.argmin(np.abs(age_grid - target_age)))
        z = int(np.clip(z, 0, g.shape[4] - 1))
        to = int(np.clip(to, 0, g.shape[1] - 1))
        nn = int(np.clip(nn, 0, g.shape[5] - 1))
        cs = int(np.clip(cs, 0, g.shape[6] - 1))
        state_weights = np.sum(g[:, to, :, j, z, nn, cs], axis=1)
        bidx = representative_wealth_index(b_grid, state_weights)
        b = float(b_grid[bidx])
        probs = tp[bidx, to, 0, j, z, nn, cs, :]
        chosen = int(choice[bidx, to, 0, j, z, nn, cs])
        max_prob = float(np.max(probs)) if probs.size else math.nan
        for tn in range(g.shape[1]):
            destination = "renter" if tn == 0 else f"owner_{float(P.H_own[tn - 1]):g}_rooms"
            prob = float(probs[tn])
            rel_value = math.nan
            if prob > 0.0 and max_prob > 0.0 and np.isfinite(kappa):
                rel_value = float(kappa * (math.log(prob) - math.log(max_prob)))
            rows.append(
                {
                    "state": label,
                    "age": float(age_grid[j]),
                    "z_index": z,
                    "current_tenure": "renter" if to == 0 else f"owner_{float(P.H_own[to - 1]):g}_rooms",
                    "parity_state_nn": nn,
                    "child_state_cs": cs,
                    "completed_children": completed_children(P, nn, cs),
                    "child_status": child_status(P, nn, cs),
                    "current_child_bin": current_child_bin(P, nn, cs),
                    "cell_mass": float(np.sum(state_weights)),
                    "liquid_wealth_cell": b,
                    "destination_tenure_index": tn,
                    "destination": destination,
                    "rooms": 0.0 if tn == 0 else float(P.H_own[tn - 1]),
                    "static_collateral_feasible": static_collateral_feasible(P, p, b, to, tn, nn, cs),
                    "choice_probability": prob,
                    "implied_relative_value_to_best": rel_value,
                    "argmax_choice": bool(tn == chosen),
                }
            )
    return rows


def representative_wealth_index(b_grid: np.ndarray, weights: np.ndarray) -> int:
    mass = float(np.sum(weights))
    if mass <= 1e-14:
        return int(np.argmin(np.abs(b_grid)))
    median = weighted_quantiles(b_grid, weights, [0.50])[0]
    return int(np.argmin(np.abs(b_grid - median)))


def static_collateral_feasible(P: Any, p: float, b: float, to: int, tn: int, nn: int, cs: int) -> bool:
    if tn == 0:
        return True
    hcost_new = p * float(P.H_own[tn - 1])
    phi = phi_for_state(P, tn, nn, cs)
    dp = (1.0 - phi) * hcost_new
    bmo = -phi * hcost_new
    if to == tn:
        return True
    sale_proceeds = 0.0 if to == 0 else (1.0 - float(P.psi)) * p * float(P.H_own[to - 1])
    if to == 0:
        post_purchase_liquid = b - hcost_new
        return bool((b >= dp) and (post_purchase_liquid >= bmo))
    required_cash_net_sale = dp - sale_proceeds
    post_purchase_liquid = b + sale_proceeds - hcost_new
    return bool((b >= required_cash_net_sale) and (post_purchase_liquid >= bmo))


def plot_family_floors(P: Any, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    k = int(getattr(P, "n_child_stages", 1))
    states = [
        ("childless", 0, 0),
        ("one dep. child", min(1, int(P.n_parity) - 1), min(1, int(P.n_child_states) - 1)),
        ("two dep. children", min(2, int(P.n_parity) - 1), min(1, int(P.n_child_states) - 1)),
        ("one child, empty nest", min(1, int(P.n_parity) - 1), min(k + 1, int(P.n_child_states) - 1)),
        ("two children, empty nest", min(2, int(P.n_parity) - 1), min(k + 2, int(P.n_child_states) - 1)),
    ]
    labels = [s[0] for s in states]
    floors = [child_floor(P, nn, cs) for _, nn, cs in states]
    x = np.arange(len(states))

    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(x, floors, color="#76B7B2", label="family-space floor")
    rung_colors = ["#4E79A7", "#F28E2B", "#E15759", "#59A14F", "#B07AA1", "#9C755F"]
    for idx, rooms in enumerate(np.asarray(P.H_own, dtype=float)):
        ax.axhline(float(rooms), color=rung_colors[idx % len(rung_colors)], lw=1.3, alpha=0.8)
        ax.text(len(states) - 0.45, float(rooms) + 0.05, f"{rooms:g}", color=rung_colors[idx % len(rung_colors)], fontsize=9)
    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylabel("rooms")
    ax.set_title("Family-space floor relative to owner rungs")
    ax.grid(axis="y", color="0.88", linewidth=0.8)
    ax.set_ylim(0.0, max(float(np.max(P.H_own)) + 1.0, max(floors) + 1.0))
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_downpayment_thresholds(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panel_groups = [
        "all_households",
        "renters_all",
        "renters_young_25_34",
        "renters_prime_childless_30_55",
        "renters_prime_parents_30_55",
        "households_old_65_75",
    ]
    title_map = {
        "all_households": "All households",
        "renters_all": "Renters",
        "renters_young_25_34": "Young renters, 25-34",
        "renters_prime_childless_30_55": "Prime childless renters",
        "renters_prime_parents_30_55": "Prime parent renters",
        "households_old_65_75": "Old households, 65-75",
    }
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey=True)
    for ax, group in zip(axes.ravel(), panel_groups):
        sub = [r for r in rows if r["wealth_group"] == group]
        sub = sorted(sub, key=lambda r: r["rung"])
        rooms = np.array([float(r["rooms"]) for r in sub])
        dp = np.array([float(r["down_payment_threshold"]) for r in sub])
        ax.plot(rooms, dp, marker="o", color="#4C78A8", label="down payment")
        if sub:
            q25 = float(sub[0]["p25_liquid_wealth"])
            q50 = float(sub[0]["p50_liquid_wealth"])
            q75 = float(sub[0]["p75_liquid_wealth"])
            q90 = float(sub[0]["p90_liquid_wealth"])
            ax.axhline(q25, color="0.65", lw=1.0, ls=":", label="p25 wealth")
            ax.axhline(q50, color="#F58518", lw=1.1, ls="--", label="p50 wealth")
            ax.axhline(q75, color="0.45", lw=1.0, ls="-.", label="p75 wealth")
            ax.axhline(q90, color="0.25", lw=1.0, ls=(0, (4, 2)), label="p90 wealth")
            mass = float(sub[0]["mass"])
        else:
            mass = math.nan
        ax.set_xticks(rooms, [f"{r:g}" for r in rooms])
        ax.set_title(f"{title_map[group]}\nmass={mass:.3f}", fontsize=10)
        ax.grid(axis="y", color="0.88", linewidth=0.8)
    for ax in axes[:, 0]:
        ax.set_ylabel("liquid wealth / cash units")
    for ax in axes[-1, :]:
        ax.set_xlabel("owner room rung")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False)
    fig.suptitle("Collateral thresholds versus liquid-wealth distribution", y=0.995)
    fig.tight_layout(rect=[0.0, 0.06, 1.0, 0.96])
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_rung_usage(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panel_groups = [
        "all_owners",
        "young_25_34",
        "prime_30_55",
        "old_65_75",
        "completed_childless",
        "completed_parent",
    ]
    title_map = {
        "all_owners": "All owners",
        "young_25_34": "Young owners, 25-34",
        "prime_30_55": "Prime owners, 30-55",
        "old_65_75": "Old owners, 65-75",
        "completed_childless": "Completed childless owners",
        "completed_parent": "Completed-parent owners",
    }
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5), sharey=True)
    for ax, group in zip(axes.ravel(), panel_groups):
        sub = [r for r in rows if r["group"] == group]
        sub = sorted(sub, key=lambda r: r["rung"])
        x = np.arange(len(sub))
        shares = [r["share_of_group_owners"] for r in sub]
        rooms = [f"{r['rooms']:g}" for r in sub]
        ax.bar(x, shares, color="#4C78A8")
        ax.set_xticks(x, rooms)
        ax.set_ylim(0.0, 1.0)
        mass = sub[0]["group_owner_mass"] if sub else math.nan
        ax.set_title(f"{title_map[group]}\nowner mass={mass:.3f}", fontsize=10)
        ax.grid(axis="y", color="0.88", linewidth=0.8)
    for ax in axes[:, 0]:
        ax.set_ylabel("share of group owners")
    for ax in axes[-1, :]:
        ax.set_xlabel("owner room rung")
    fig.suptitle("Current owner ladder usage", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_representative_probabilities(rows: list[dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    states = []
    for row in rows:
        if row["state"] not in states:
            states.append(row["state"])
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey=True)
    for ax, state in zip(axes.ravel(), states):
        sub = [r for r in rows if r["state"] == state]
        sub = sorted(sub, key=lambda r: r["destination_tenure_index"])
        labels = ["R"] + [f"{r['rooms']:g}" for r in sub if r["destination_tenure_index"] > 0]
        probs = [r["choice_probability"] for r in sub]
        feasible = [r["static_collateral_feasible"] for r in sub]
        colors = ["#59A14F" if ok else "#BAB0AC" for ok in feasible]
        x = np.arange(len(sub))
        ax.bar(x, probs, color=colors)
        ax.set_xticks(x, labels)
        ax.set_ylim(0.0, 1.0)
        first = sub[0]
        ax.set_title(
            f"{state.replace('_', ' ')}\n"
            f"b={first['liquid_wealth_cell']:.2f}, mass={first['cell_mass']:.3g}",
            fontsize=9,
        )
        ax.grid(axis="y", color="0.88", linewidth=0.8)
    for ax in axes[:, 0]:
        ax.set_ylabel("choice probability")
    for ax in axes[-1, :]:
        ax.set_xlabel("destination: renter or owner rooms")
    fig.suptitle("Representative next-tenure/rung probabilities", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def build_audit_readme(payload: dict[str, Any], base: dict[str, Any], sol: Any, P: Any, outdir: Path) -> str:
    label = base.get("label", "baseline")
    loss = base.get("rank_loss", math.nan)
    p = house_price(P, sol)
    lines = [
        "# Housing Ladder Audit",
        "",
        "Generated from the saved current-run solution cache. This does not re-solve the model.",
        "",
        f"- Generated: `{datetime.now().isoformat(timespec='seconds')}`",
        f"- Source record: `{payload.get('source_meta', {}).get('path', 'unknown')}`",
        f"- Baseline label: `{label}`",
        f"- Rank loss in cache: `{loss}`",
        f"- Asset price `p`: `{p:.6g}`",
        f"- PTI constraint active: `{bool(getattr(P, 'use_pti_constraint', False))}`",
        f"- Owner rungs: `{', '.join(f'{x:g}' for x in np.asarray(P.H_own, dtype=float))}`",
        "",
        "## Files",
        "",
        "- `housing_rung_usage_by_group.png`: six-panel mass-share plot for current owner rung occupancy.",
        "- `rung_mass_by_group.csv`: all rung shares by age, completed-fertility, current-child, and income groups.",
        "- `family_floor_vs_rungs.png`: family-space floors overlaid with the owner room rungs.",
        "- `family_floor_by_rung.csv`: room floors versus each owner rung by family state.",
        "- `downpayment_thresholds_vs_wealth.png`: down-payment thresholds versus selected liquid-wealth quantiles.",
        "- `downpayment_thresholds_vs_wealth.csv`: down-payment thresholds compared with liquid-wealth distributions.",
        "- `representative_rung_choice_probabilities.png`: representative next-rung probabilities; gray bars fail the static collateral test at that cell.",
        "- `representative_rung_choice_probabilities.csv`: source table, including implied logit value gaps relative to the best destination.",
        "",
        "Important interpretation: `family_floor_by_rung.csv` reports residual room slack, not a hard feasibility flag. "
        "The owner utility code maps nonpositive residual housing services to a near-zero service flow, so a rung below the floor can be technically feasible but economically dominated.",
        "",
    ]
    return "\n".join(lines)


def write_plot_index(run_dir: Path, outdir: Path, payload: dict[str, Any], base: dict[str, Any], sol: Any, P: Any) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    p = house_price(P, sol)
    lines = [
        "# Current Intergen Run Plots",
        "",
        "This is the stable local folder for the current one-market intergen solution plots. "
        "Use this folder for quick visual inspection; temporary experiments should write to separate dated/testing folders.",
        "",
        f"- Folder: `{run_dir}`",
        f"- Source record: `{payload.get('source_meta', {}).get('path', 'unknown')}`",
        f"- Target set: `{payload.get('target_set', 'unknown')}`",
        f"- Cache rank loss: `{base.get('rank_loss', math.nan)}`",
        f"- Market price `p`: `{p:.6g}`",
        f"- PTI constraint active: `{bool(getattr(P, 'use_pti_constraint', False))}`",
        f"- Updated: `{datetime.now().isoformat(timespec='seconds')}`",
        "",
        "## First-Look Plots",
        "",
        "- `first_look_policies_markets.png`: compact policies and prices/quantities.",
        "- `first_look_policies_markets_total_wealth.png`: same first-look view against total wealth.",
        "- `first_look_wealth_density.png`: liquid-wealth density support.",
        "- `first_look_total_wealth_density.png`: total-wealth density support.",
        "",
        "## Housing Ladder Plots",
        "",
        "- `owner_rung_shares_all_owners.png`: simple current owner rung shares.",
        "- `diagnostics/owner_rungs.png`: owner housing demand by rung in service units.",
        "- `ladder_audit/family_floor_vs_rungs.png`: family-space floors versus the owner rung grid.",
        "- `ladder_audit/downpayment_thresholds_vs_wealth.png`: collateral thresholds versus wealth quantiles.",
        "- `ladder_audit/housing_rung_usage_by_group.png`: current owner rung mass shares by lifecycle/family group.",
        "- `ladder_audit/representative_rung_choice_probabilities.png`: representative next-rung probabilities and collateral feasibility.",
        "",
        "## Regenerating Without Solving",
        "",
        "Run:",
        "",
        "```bash",
        "cd code/model",
        ".venv/bin/python tools/audit_intergen_housing_ladder.py",
        "```",
        "",
        "For first-look policy/market plot refreshes, open `code/model/run_intergen_model.py` and leave "
        "`REFRESH_PLOTS_FROM_SAVED_SOLUTION = True` and `FAST_REFRESH_FROM_SAVED_SOLUTION = True`.",
        "",
    ]
    (run_dir / "PLOTS_INDEX.md").write_text("\n".join(lines), encoding="utf-8")

    current_note = ROOT / "output" / "model" / "CURRENT_INTERGEN_RUN.md"
    current_note.write_text(
        "\n".join(
            [
                "# Current Intergen Run",
                "",
                f"Stable plot folder: `{run_dir}`",
                "",
                "Open `PLOTS_INDEX.md` inside that folder for the map of the current figures.",
                "",
                "Main inspection figures:",
                "",
                f"- `{run_dir / 'first_look_policies_markets.png'}`",
                f"- `{run_dir / 'first_look_policies_markets_total_wealth.png'}`",
                f"- `{run_dir / 'owner_rung_shares_all_owners.png'}`",
                f"- `{outdir / 'housing_rung_usage_by_group.png'}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
