#!/usr/bin/env python3
"""Collect parallel repaired-model policy tasks into tables, plots, and a readout."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CASE_ORDER = [
    "parent_ltv95_birth", "parent_ltv100_birth", "universal_ltv95",
    "birth_grant_A0p1_Hge6", "birth_grant_A0p2_Hge6", "birth_grant_A0p4_Hge6",
    "birth_grant_A0p58_Hge6", "birth_grant_A1p0_Hge6", "birth_grant_A0p4_all_rungs",
    "property_tax_2pct", "tax2_grant_A0p4_Hge6",
    "rental_hR6p5", "rental_hR7", "rental_hR7p5", "rental_hR8",
    "supply_H0_plus10", "supply_H0_plus20", "grant_A0p4_Hge6_supply_plus20",
    "debt_line_lambda0p2", "debt_line_lambda0p4", "debt_line0p4_supply_plus20",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--job-id", default="")
    parser.add_argument("--allow-incomplete", action="store_true", help="Permit a partial battery for smoke testing.")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    baseline, policies, targets, weights = load_tasks(args.results_root, allow_incomplete=args.allow_incomplete)
    rows = [summary_row(case, baseline) for case in policies]
    rows.sort(key=lambda row: CASE_ORDER.index(row["case"]) if row["case"] in CASE_ORDER else 999)
    write_csv(args.outdir / "policy_summary.csv", rows)
    write_target_fit(args.outdir / "policy_target_fit.csv", baseline, policies, targets, weights)
    plot_effects(args.outdir / "policy_effects.png", rows)
    plot_grant_dose(args.outdir / "grant_dose_response.png", rows)
    write_readme(args.outdir / "README.md", args.results_root, args.job_id, baseline, rows)


def load_tasks(root: Path, *, allow_incomplete: bool = False) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, float], dict[str, float]]:
    baselines: list[dict[str, Any]] = []
    policies: list[dict[str, Any]] = []
    targets: dict[str, float] | None = None
    weights: dict[str, float] | None = None
    for task in sorted(root.glob("task_*")):
        records_path = task / "policy_cases" / "policy_cases.json"
        metadata_path = task / "metadata.json"
        if not records_path.exists() or not metadata_path.exists():
            continue
        records = json.loads(records_path.read_text())
        metadata = json.loads(metadata_path.read_text())
        if targets is None:
            targets = {str(k): float(v) for k, v in metadata["targets"].items()}
            weights = {str(k): float(v) for k, v in metadata["weights"].items()}
        baseline_rows = [row for row in records if row.get("case") == "baseline"]
        policy_rows = [row for row in records if row.get("case") != "baseline"]
        if len(baseline_rows) != 1 or len(policy_rows) != 1:
            raise ValueError(f"unexpected policy task contents: {records_path}")
        baselines.append(baseline_rows[0])
        policies.append(policy_rows[0])
    if not baselines or targets is None or weights is None:
        raise FileNotFoundError(f"no completed policy tasks found under {root}")
    canonical = baselines[0]
    found_cases = [str(row.get("case")) for row in policies]
    if len(found_cases) != len(set(found_cases)):
        raise ValueError("duplicate policy cases found across task outputs")
    missing = sorted(set(CASE_ORDER) - set(found_cases))
    extra = sorted(set(found_cases) - set(CASE_ORDER))
    if not allow_incomplete and (missing or extra):
        raise ValueError(f"incomplete policy battery: missing={missing}, extra={extra}")
    for row in [canonical, *policies]:
        if not math.isfinite(float(row["rank_loss"])):
            raise ValueError(f"non-finite policy loss for {row.get('case')}")
        residual = float(row["market_residual"])
        if not math.isfinite(residual) or residual > 1e-4:
            raise ValueError(f"policy equilibrium failed strict residual gate for {row.get('case')}: {residual}")
    for other in baselines[1:]:
        for key in ("rank_loss", "market_residual"):
            if not math.isclose(float(canonical[key]), float(other[key]), rel_tol=0.0, abs_tol=1e-12):
                raise ValueError(f"baseline mismatch in {key}: {canonical[key]} versus {other[key]}")
        if not np.allclose(canonical["p_eq"], other["p_eq"], rtol=0.0, atol=1e-12):
            raise ValueError("baseline price mismatch across policy tasks")
    return canonical, policies, targets, weights


def summary_row(case: dict[str, Any], baseline: dict[str, Any]) -> dict[str, Any]:
    b = baseline["moments"]
    m = case["moments"]
    p0 = float(baseline["p_eq"][0])
    p1 = float(case["p_eq"][0])
    row: dict[str, Any] = {
        "case": case["case"],
        "label": case["label"],
        "note": case["note"],
        "rank_loss": float(case["rank_loss"]),
        "market_residual": float(case["market_residual"]),
        "price": p1,
        "price_change_pct": 100.0 * (p1 / p0 - 1.0),
    }
    debt = case.get("debt_diagnostics", {})
    row["post_taper_renter_negative_b_mass"] = maybe_float(
        debt.get("post_taper_renter_negative_b_mass")
    )
    row["post_taper_owner_negative_unsecured_mass"] = maybe_float(
        debt.get("post_taper_owner_negative_unsecured_mass")
    )
    keys = [
        "tfr", "childless_rate", "own_rate_2534", "own_rate", "old_age_own_rate",
        "housing_increment_0to1", "own_family_gap", "aggregate_mean_occupied_rooms_18_85",
        "prime30_55_childless_owner_minus_renter_mean_rooms",
        "prime30_55_parent_3plus_minus_1to2_mean_rooms",
        "young_childless_renter_liquid_wealth_to_annual_gross_income_2535",
    ]
    for key in keys:
        row[key] = maybe_float(m.get(key))
        row[f"delta_{key}"] = maybe_float(m.get(key)) - maybe_float(b.get(key))
    return row


def write_target_fit(path: Path, baseline: dict[str, Any], policies: list[dict[str, Any]], targets: dict[str, float], weights: dict[str, float]) -> None:
    rows = []
    for case in [baseline, *policies]:
        for key in sorted(targets):
            model = maybe_float(case["moments"].get(key))
            gap = model - targets[key]
            rows.append({
                "case": case["case"], "moment": key, "target": targets[key], "model": model,
                "gap": gap, "weight": weights[key], "loss_contribution": weights[key] * gap * gap,
            })
    write_csv(path, rows)


def plot_effects(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [short_label(row["case"]) for row in rows]
    y = np.arange(len(rows))
    series = [
        ("delta_tfr", r"$\Delta$ completed fertility", 1.0),
        ("delta_own_rate_2534", r"$\Delta$ ownership, 25--34 (pp)", 100.0),
        ("delta_housing_increment_0to1", r"$\Delta$ first-birth rooms", 1.0),
        ("price_change_pct", "House-price change (%)", 1.0),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(16, max(7.5, 0.38 * len(rows))), sharey=True)
    for ax, (key, title, scale) in zip(axes, series):
        vals = scale * np.asarray([float(row[key]) for row in rows])
        colors = np.where(vals >= 0.0, "#35618f", "#c46b3c")
        ax.barh(y, vals, color=colors)
        ax.axvline(0.0, color="0.25", lw=0.8)
        ax.set_title(title)
        ax.grid(axis="x", alpha=0.2)
        if ax is axes[0]:
            ax.set_yticks(y, labels)
        else:
            ax.tick_params(axis="y", labelleft=False)
    axes[0].invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_grant_dose(path: Path, rows: list[dict[str, Any]]) -> None:
    doses = {"birth_grant_A0p1_Hge6": 0.1, "birth_grant_A0p2_Hge6": 0.2,
             "birth_grant_A0p4_Hge6": 0.4, "birth_grant_A0p58_Hge6": 0.58,
             "birth_grant_A1p0_Hge6": 1.0}
    panel = sorted([(doses[row["case"]], row) for row in rows if row["case"] in doses])
    if not panel:
        return
    x = np.asarray([dose for dose, _ in panel])
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8))
    specs = [
        ("delta_tfr", r"$\Delta$ completed fertility"),
        ("delta_own_rate_2534", r"$\Delta$ ownership, 25--34 (pp)"),
        ("price_change_pct", "House-price change (%)"),
    ]
    for ax, (key, title) in zip(axes, specs):
        scale = 100.0 if key == "delta_own_rate_2534" else 1.0
        ax.plot(x, [scale * row[key] for _, row in panel], marker="o", color="#35618f", lw=2)
        ax.axhline(0.0, color="0.3", lw=0.8)
        ax.set_xlabel("Birth grant (model-period goods units)")
        ax.set_title(title)
        ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_readme(path: Path, results_root: Path, job_id: str, baseline: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    debt_rows = [row for row in rows if row["case"].startswith("debt_")]
    debt_taper_ok = all(
        abs(row["post_taper_renter_negative_b_mass"]) <= 1e-6
        and abs(row["post_taper_owner_negative_unsecured_mass"]) <= 1e-6
        for row in debt_rows
    )
    best_tfr = max(rows, key=lambda row: row["delta_tfr"])
    best_young_own = max(rows, key=lambda row: row["delta_own_rate_2534"])
    best_h01 = max(rows, key=lambda row: row["delta_housing_increment_0to1"])
    lowest_price = min(rows, key=lambda row: row["price_change_pct"])
    lines = [
        "# Full repaired-model policy battery", "",
        f"Torch job: `{job_id or 'not recorded'}`. Source task root: `{results_root}`.", "",
        "All cases hold the calibrated parameter vector fixed, use `J=17`, `Nb=120`, the repaired feasibility rules, and re-clear the one-market housing equilibrium.", "",
        f"Baseline: loss `{float(baseline['rank_loss']):.6f}`, residual `{float(baseline['market_residual']):.2e}`, price `{float(baseline['p_eq'][0]):.6f}`.", "",
        f"Debt-taper acceptance (renter debt and owner unsecured-debt mass from age 62 each below `1e-6`): `{'PASS' if debt_taper_ok else 'FAIL'}`.", "",
        "## Mechanical readout", "",
        f"- Largest fertility response: **{best_tfr['label']}**, delta TFR `{best_tfr['delta_tfr']:+.4f}`.",
        f"- Largest young-ownership response: **{best_young_own['label']}**, `{100*best_young_own['delta_own_rate_2534']:+.2f}` percentage points.",
        f"- Largest first-birth housing response: **{best_h01['label']}**, `{best_h01['delta_housing_increment_0to1']:+.3f}` rooms.",
        f"- Largest house-price reduction: **{lowest_price['label']}**, `{lowest_price['price_change_pct']:+.2f}%`.", "",
        "| Case | Delta TFR | Delta childless (pp) | Delta own 25--34 (pp) | Delta H01 | Price change | Residual |", "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['label']} | {row['delta_tfr']:+.4f} | {100*row['delta_childless_rate']:+.2f} | "
            f"{100*row['delta_own_rate_2534']:+.2f} | {row['delta_housing_increment_0to1']:+.3f} | "
            f"{row['price_change_pct']:+.2f}% | {row['market_residual']:.2e} |"
        )
    lines.extend([
        "", "## Interpretation boundaries", "",
        "- Grant and LTV cases have no government-budget, guarantee-cost, default, or lender-loss accounting.",
        "- The tax-plus-grant case is a package, not a verified balanced-budget policy.",
        "- `H0` expansions and renter-cap relief have no construction or resource cost in this diagnostic.",
        "- The unsecured line cannot fund the cash down payment and has no spread/default model; it is a smoothing sensitivity.",
        "- These are mechanism counterfactuals, not welfare estimates. `policy_target_fit.csv` reports all fifteen moments for every case.",
    ])
    path.write_text("\n".join(lines) + "\n")


def short_label(case: str) -> str:
    labels = {
        "parent_ltv95_birth": "Parent LTV 95", "parent_ltv100_birth": "Parent LTV 100",
        "universal_ltv95": "Universal LTV 95", "property_tax_2pct": "Tax 2%",
        "tax2_grant_A0p4_Hge6": "Tax + grant .4", "birth_grant_A0p4_all_rungs": "Grant .4, all homes",
        "supply_H0_plus10": "Supply +10%", "supply_H0_plus20": "Supply +20%",
        "grant_A0p4_Hge6_supply_plus20": "Grant .4 + supply",
        "debt_line_lambda0p2": "Debt line .2", "debt_line_lambda0p4": "Debt line .4",
        "debt_line0p4_supply_plus20": "Debt .4 + supply",
    }
    if case in labels:
        return labels[case]
    if case.startswith("birth_grant_A"):
        return "Grant " + case.removeprefix("birth_grant_A").removesuffix("_Hge6").replace("p", ".")
    if case.startswith("rental_hR"):
        return "Rent cap " + case.removeprefix("rental_hR").replace("p", ".")
    return case.replace("_", " ")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def maybe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


if __name__ == "__main__":
    main()
