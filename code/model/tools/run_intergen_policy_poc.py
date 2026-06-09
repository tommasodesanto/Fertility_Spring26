#!/usr/bin/env python3
"""Run proof-of-concept one-market policy comparisons from the final toy best.

This is not a production counterfactual engine. It keeps the saved final
global-DE parameter vector fixed, solves a small set of policy variants, and
writes a standard comparison packet: prices/user costs, headline moments,
lifecycle profiles, distributions, owner-entry policies, and owner-rung use.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_housing_fertility.calibration import (  # noqa: E402
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_housing_fertility.diagnostics import write_diagnostics  # noqa: E402
from intergen_housing_fertility.local_panel import income_process_overrides  # noqa: E402
from intergen_housing_fertility.solver import run_model_cp_dt  # noqa: E402


DEFAULT_SUMMARY = (
    ROOT
    / "output/model/intergen_shutdown_snapshots/final_manual_summary/"
    / "shutdown_snapshots/snapshot_ny0642_manual_final_20260609/combined_summary.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/intergen_policy_poc_20260609"
TARGET_SET = "candidate_no_timing_v0"

MOMENT_PANEL_KEYS = [
    "tfr",
    "childless_rate",
    "own_rate",
    "own_family_gap",
    "old_age_own_rate",
    "old_age_parent_childless_gap",
    "housing_increment_0to1",
    "housing_increment_1to2",
    "young_liquid_wealth_to_income",
    "liquid_wealth_to_income",
    "housing_user_cost_share",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--Nb", type=int, default=60)
    parser.add_argument("--max-iter-eq", type=int, default=25)
    parser.add_argument("--room-ladder", action="store_true", help="Use [2,4,6,8,9.5,11] instead of the saved cluster ladder")
    parser.add_argument("--skip-standard-diagnostics", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    source = json.loads(args.summary.read_text())["best"]
    theta = dict(source["theta"])
    targets, weights = get_target_set(TARGET_SET)

    cases = default_cases()
    records: list[dict[str, Any]] = []
    arrays: dict[str, dict[str, Any]] = {}
    for case in cases:
        record, arr = solve_case(
            case,
            theta=theta,
            targets=targets,
            weights=weights,
            outdir=args.outdir,
            Nb=int(args.Nb),
            max_iter_eq=int(args.max_iter_eq),
            room_ladder=bool(args.room_ladder),
            write_standard_diagnostics=not bool(args.skip_standard_diagnostics),
        )
        records.append(record)
        arrays[record["case"]] = arr
        write_json(args.outdir / "latest_completed_case.json", record)
        write_outputs_checkpoint(args.outdir, source, records)
        print(
            f"{record['case']}: loss={record['rank_loss']:.4f}, "
            f"resid={record['market_residual']:.2e}, elapsed={record['elapsed_sec']:.1f}s",
            flush=True,
        )

    write_json(args.outdir / "policy_summary.json", {"source_record": source, "cases": records})
    write_csv(args.outdir / "policy_summary.csv", records)
    write_comparison_figures(args.outdir / "figures", records, arrays, targets)
    write_readme(args.outdir, source, records)


def default_cases() -> list[dict[str, Any]]:
    return [
        {
            "case": "baseline",
            "label": "Baseline toy calibration",
            "overrides": {},
            "note": "Final global-DE toy calibration. Baseline estate_tax_rate is zero.",
        },
        {
            "case": "new_parent_ltv95",
            "label": "New-parent credit relief",
            "overrides": {
                "parent_dp_waiver": True,
                "parent_dp_waiver_phi": 0.95,
                "parent_dp_waiver_birth_state_only": True,
            },
            "note": "Raises financed share to 95 percent only for owner choices by new parents in the birth child-state.",
        },
        {
            "case": "property_tax_up_1pp",
            "label": "Property tax +1pp annual",
            "overrides": {"tau_H": 0.02 * 4.0},
            "note": "Raises the annual property tax from 1 percent to 2 percent. Prices re-clear under user-cost capitalization.",
        },
        {
            "case": "estate_tax_30pct",
            "label": "Estate tax 30 percent",
            "overrides": {"estate_tax_rate": 0.30},
            "note": "Proof-of-concept terminal bequest tax wedge. No estate-tax rebate or inheritance-transfer kernel is included.",
        },
    ]


def solve_case(
    case: dict[str, Any],
    *,
    theta: dict[str, float],
    targets: dict[str, float],
    weights: dict[str, float],
    outdir: Path,
    Nb: int,
    max_iter_eq: int,
    room_ladder: bool,
    write_standard_diagnostics: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    t0 = time.perf_counter()
    overrides = {
        **base_overrides(J=16, Nb=Nb, n_house=6, max_iter_eq=max_iter_eq),
        **income_process_overrides(5),
        **theta,
        **dict(case["overrides"]),
    }
    if room_ladder:
        overrides["H_own"] = np.asarray([2.0, 4.0, 6.0, 8.0, 9.5, 11.0])
        overrides["n_house"] = 6

    sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
    elapsed = time.perf_counter() - t0
    moments = extract_moments(sol, P)
    rank_loss = diagnostic_loss(moments, targets=targets, weights=weights)
    case_dir = outdir / str(case["case"])
    case_dir.mkdir(parents=True, exist_ok=True)
    if write_standard_diagnostics:
        write_diagnostics(sol, P, case_dir / "diagnostics")

    arr = extract_arrays(sol, P)
    metrics = extract_metrics(sol, P, moments, arr)
    record = {
        "case": str(case["case"]),
        "label": str(case["label"]),
        "note": str(case["note"]),
        "overrides": jsonable(case["overrides"]),
        "elapsed_sec": float(elapsed),
        "rank_loss": float(rank_loss),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "p_eq": float(np.asarray(p_eq).reshape(-1)[0]),
        "owner_asset_price": float(np.asarray(sol.owner_asset_price).reshape(-1)[0]),
        "owner_user_cost": float(np.asarray(sol.owner_user_cost).reshape(-1)[0]),
        "H_own": jsonable(np.asarray(P.H_own, dtype=float)),
        "hR_max": float(P.hR_max),
        "estate_tax_rate": float(getattr(P, "estate_tax_rate", 0.0)),
        "tau_H_period": float(P.tau_H),
        "moments": {k: float(moments.get(k, np.nan)) for k in sorted(set(MOMENT_PANEL_KEYS + list(targets)))},
        "diagnostic_moments": {
            "mean_age_first_birth": float(moments.get("mean_age_first_birth", np.nan)),
            "prime_childless_renter_median_rooms": float(moments.get("prime_childless_renter_median_rooms", np.nan)),
            "prime_childless_owner_median_rooms": float(moments.get("prime_childless_owner_median_rooms", np.nan)),
            "aggregate_own_rate": float(moments.get("aggregate_own_rate", np.nan)),
            "own_rate_2534": float(moments.get("own_rate_2534", np.nan)),
            "own_rate_3544": float(moments.get("own_rate_3544", np.nan)),
            "owner_neg_liquid_share_2534": float(moments.get("owner_neg_liquid_share_2534", np.nan)),
            "renter25_45_all_cap_share": float(moments.get("renter25_45_all_cap_share", np.nan)),
        },
        "metrics": metrics,
        "case_dir": str(case_dir),
    }
    write_json(case_dir / "record.json", record)
    write_csv(case_dir / "asset_distribution_by_age.csv", arr["asset_rows"])
    return record, arr


def extract_arrays(sol: Any, P: Any) -> dict[str, Any]:
    ages = np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float)
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    wealth_mean = []
    max_bin_share = []
    effective_bins = []
    asset_rows = []
    for j, age in enumerate(ages):
        gj = g[:, :, :, j, :, :, :]
        mass = float(np.sum(gj))
        mass_b = np.sum(gj, axis=tuple(range(1, gj.ndim)))
        share = mass_b / max(mass, 1e-14)
        wealth_mean.append(float(np.sum(share * b_grid)))
        max_bin_share.append(float(np.max(share)))
        effective_bins.append(float(1.0 / max(np.sum(share**2), 1e-14)))
        asset_rows.append(
            {
                "age": float(age),
                "mass": mass,
                "mean_liquid_wealth": wealth_mean[-1],
                "max_asset_bin_share": max_bin_share[-1],
                "effective_asset_bins": effective_bins[-1],
                "near_zero_mass_abs_le_0p5": float(np.sum(share[np.abs(b_grid) <= 0.5])),
                "negative_mass": float(np.sum(share[b_grid < 0.0])),
            }
        )
    return {
        "ages": ages,
        "b_grid": b_grid,
        "own_by_age": np.asarray(sol.own_by_age, dtype=float),
        "fert_by_age": np.asarray(getattr(sol, "fert_by_age", np.zeros(P.J)), dtype=float),
        "mean_liquid_wealth_by_age": np.asarray(wealth_mean, dtype=float),
        "max_asset_bin_share_by_age": np.asarray(max_bin_share, dtype=float),
        "effective_asset_bins_by_age": np.asarray(effective_bins, dtype=float),
        "owner_rung_services": np.asarray(getattr(sol, "owner_demand_by_size", np.zeros(P.n_house)), dtype=float),
        "H_own": np.asarray(P.H_own, dtype=float),
        "tenure_services": np.asarray([float(sol.aggregate_rental_demand), float(sol.aggregate_owner_demand)]),
        "owner_entry_age30_z1": owner_entry_line(sol, P, age=30.0, z_target=1.0),
        "owner_entry_age42_z1": owner_entry_line(sol, P, age=42.0, z_target=1.0),
        "asset_rows": asset_rows,
    }


def owner_entry_line(sol: Any, P: Any, *, age: float, z_target: float) -> np.ndarray:
    b_grid = np.asarray(sol.b_grid, dtype=float)
    j = int(np.clip(round((age - float(P.age_start)) / float(P.da)), 0, P.J - 1))
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", [1.0])), dtype=float).reshape(-1)
    zz = int(np.argmin(np.abs(z_grid - float(z_target))))
    tp = getattr(sol, "tenure_probs", None)
    if tp is None:
        line = (np.asarray(sol.tenure_choice)[:, 0, 0, j, zz, 0, 0] > 0).astype(float)
    else:
        line = np.sum(np.asarray(tp)[:, 0, 0, j, zz, 0, 0, 1:], axis=1)
    return np.column_stack([b_grid, line])


def extract_metrics(sol: Any, P: Any, moments: dict[str, float], arr: dict[str, Any]) -> dict[str, float]:
    owner_services = np.asarray(arr["owner_rung_services"], dtype=float)
    total_owner_services = float(np.sum(owner_services))
    if total_owner_services > 1e-14:
        owner_share = owner_services / total_owner_services
        top_idx = int(np.argmax(owner_share))
        owner_rung_top_share = float(owner_share[top_idx])
        owner_rung_top_h = float(np.asarray(P.H_own, dtype=float)[top_idx])
    else:
        owner_rung_top_share = np.nan
        owner_rung_top_h = np.nan
    return {
        "old_minus_prime_ownership": float(moments.get("old_age_own_rate", np.nan) - moments.get("own_rate", np.nan)),
        "owner_rung_top_share": owner_rung_top_share,
        "owner_rung_top_h": owner_rung_top_h,
        "max_asset_bin_share_age30": float(value_at_age(arr, "max_asset_bin_share_by_age", 30.0)),
        "effective_asset_bins_age42": float(value_at_age(arr, "effective_asset_bins_by_age", 42.0)),
        "owner_entry_age30_z1_max_drop": max_downward_drop(arr["owner_entry_age30_z1"][:, 1]),
        "owner_entry_age42_z1_max_drop": max_downward_drop(arr["owner_entry_age42_z1"][:, 1]),
    }


def value_at_age(arr: dict[str, Any], key: str, age: float) -> float:
    ages = np.asarray(arr["ages"], dtype=float)
    idx = int(np.argmin(np.abs(ages - float(age))))
    return float(np.asarray(arr[key], dtype=float)[idx])


def max_downward_drop(line: np.ndarray) -> float:
    diff = np.diff(np.asarray(line, dtype=float))
    return float(np.min(diff)) if diff.size else np.nan


def write_comparison_figures(outdir: Path, records: list[dict[str, Any]], arrays: dict[str, dict[str, Any]], targets: dict[str, float]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    color_map = {r["case"]: colors[i % len(colors)] for i, r in enumerate(records)}

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    labels = [r["case"] for r in records]
    x = np.arange(len(records))
    axes[0].bar(x, [r["owner_asset_price"] for r in records], color=[color_map[r["case"]] for r in records])
    axes[0].set_title("Asset price")
    axes[0].set_xticks(x, labels, rotation=25, ha="right")
    axes[1].bar(x, [r["owner_user_cost"] for r in records], color=[color_map[r["case"]] for r in records])
    axes[1].set_title("Flow user cost / rent")
    axes[1].set_xticks(x, labels, rotation=25, ha="right")
    for ax in axes:
        ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "01_prices_user_costs.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(4, 3, figsize=(12.0, 10.0))
    axes = axes.ravel()
    for ax, key in zip(axes, MOMENT_PANEL_KEYS):
        vals = [r["moments"].get(key, np.nan) for r in records]
        ax.bar(x, vals, color=[color_map[r["case"]] for r in records])
        if key in targets:
            ax.axhline(targets[key], color="0.25", linestyle="--", lw=1.0)
        ax.set_title(key)
        ax.set_xticks(x, labels, rotation=35, ha="right", fontsize=8)
        ax.grid(axis="y", alpha=0.2)
    for ax in axes[len(MOMENT_PANEL_KEYS) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(outdir / "02_headline_moments.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for r in records:
        arr = arrays[r["case"]]
        axes[0].plot(arr["ages"], arr["own_by_age"], lw=1.8, label=r["case"], color=color_map[r["case"]])
        axes[1].plot(arr["ages"], arr["fert_by_age"], lw=1.8, label=r["case"], color=color_map[r["case"]])
    axes[0].set_title("Ownership by age")
    axes[0].set_ylabel("ownership rate")
    axes[0].set_ylim(0.0, 1.05)
    axes[1].set_title("Fertility by age")
    axes[1].set_ylabel("expected births")
    for ax in axes:
        ax.set_xlabel("age")
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "03_lifecycle_profiles.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 8.0), sharex=True)
    for r in records:
        arr = arrays[r["case"]]
        axes[0].plot(arr["ages"], arr["mean_liquid_wealth_by_age"], lw=1.8, label=r["case"], color=color_map[r["case"]])
        axes[1].plot(arr["ages"], arr["max_asset_bin_share_by_age"], lw=1.8, label=r["case"], color=color_map[r["case"]])
        axes[2].plot(arr["ages"], arr["effective_asset_bins_by_age"], lw=1.8, label=r["case"], color=color_map[r["case"]])
    axes[0].set_ylabel("mean liquid wealth")
    axes[1].set_ylabel("max asset-bin share")
    axes[2].set_ylabel("effective asset bins")
    axes[2].set_xlabel("age")
    for ax in axes:
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "04_distributions.png", dpi=180)
    plt.close(fig)

    for age_key, filename, title in [
        ("owner_entry_age30_z1", "05_owner_entry_policy_age30_z1.png", "Owner-entry policy, age 30, z=1"),
        ("owner_entry_age42_z1", "06_owner_entry_policy_age42_z1.png", "Owner-entry policy, age 42, z=1"),
    ]:
        fig, ax = plt.subplots(figsize=(8.0, 4.6))
        for r in records:
            line = arrays[r["case"]][age_key]
            ax.plot(line[:, 0], line[:, 1], lw=1.8, label=r["case"], color=color_map[r["case"]])
        ax.set_title(title)
        ax.set_xlabel("liquid wealth")
        ax.set_ylabel("owner-entry probability")
        ax.set_ylim(-0.05, 1.05)
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir / filename, dpi=180)
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    width = 0.8 / max(len(records), 1)
    for i, r in enumerate(records):
        arr = arrays[r["case"]]
        owner_services = np.asarray(arr["owner_rung_services"], dtype=float)
        share = owner_services / max(float(np.sum(owner_services)), 1e-14)
        offset = (i - 0.5 * (len(records) - 1)) * width
        axes[0].bar(np.arange(len(share)) + offset, share, width=width, label=r["case"], color=color_map[r["case"]])
        axes[1].bar(i - 0.2, arr["tenure_services"][0], width=0.38, color=color_map[r["case"]], alpha=0.7)
        axes[1].bar(i + 0.2, arr["tenure_services"][1], width=0.38, color=color_map[r["case"]])
    axes[0].set_title("Owner rung service shares")
    h_grid = np.asarray(arrays[records[0]["case"]]["H_own"], dtype=float)
    axes[0].set_xticks(np.arange(len(h_grid)), [f"{h:g}" for h in h_grid])
    axes[0].set_xlabel("owner rung")
    axes[0].set_ylabel("share of owner services")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set_title("Housing services by tenure")
    axes[1].set_xticks(x, labels, rotation=25, ha="right")
    axes[1].set_ylabel("service units per adult")
    for ax in axes:
        ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "07_tenure_owner_rungs.png", dpi=180)
    plt.close(fig)


def write_readme(outdir: Path, source: dict[str, Any], records: list[dict[str, Any]]) -> None:
    lines = [
        "# One-Market Policy Proof Of Concept",
        "",
        "Fixed-theta policy comparison from the final global-DE toy calibration. This is a proof of concept, not a production counterfactual.",
        "",
        f"Source: run `{source.get('run_dir')}`, task `{source.get('task_id')}`, case `{source.get('case')}`, label `{source.get('label')}`.",
        "",
        "## Cases",
        "",
        "| Case | Loss | Price | User cost | TFR | Childless | Prime own | Old own | Housing share |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in records:
        m = r["moments"]
        lines.append(
            f"| `{r['case']}` | {r['rank_loss']:.3f} | {r['owner_asset_price']:.3f} | {r['owner_user_cost']:.3f} | "
            f"{m['tfr']:.3f} | {m['childless_rate']:.3f} | {m['own_rate']:.3f} | "
            f"{m['old_age_own_rate']:.3f} | {m['housing_user_cost_share']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Notes",
            "",
            "- `new_parent_ltv95` is the implemented parent-targeted credit-relief case: financed share rises to 95 percent only in the birth child-state.",
            "- `property_tax_up_1pp` raises the annual property tax rate from 1 percent to 2 percent and lets prices re-clear.",
            "- `estate_tax_30pct` is a proof-of-concept bequest-tax wedge inside terminal bequest utility. The model still has no estate-tax rebate or inheritance-transfer kernel. Since baseline has zero estate tax, tax removal is read as `baseline` relative to `estate_tax_30pct`.",
            "- The final toy calibration has known lifecycle and policy-function pathologies; use these comparisons to inspect mechanisms, not as quantitative policy estimates.",
            "",
            "Standard comparison panels are in `figures/`. Each case also has a full standard diagnostic packet under `<case>/diagnostics/`.",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


def write_outputs_checkpoint(outdir: Path, source: dict[str, Any], records: list[dict[str, Any]]) -> None:
    write_json(outdir / "policy_summary_partial.json", {"source_record": source, "cases": records})
    write_csv(outdir / "policy_summary_partial.csv", records)


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(obj), indent=2, sort_keys=True))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    flat = [flatten_row(row) for row in rows]
    fields: list[str] = []
    for row in flat:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat)


def flatten_row(row: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in row.items():
        name = f"{prefix}{key}"
        if isinstance(value, dict):
            out.update(flatten_row(value, f"{name}."))
        elif isinstance(value, (list, tuple, np.ndarray)):
            out[name] = json.dumps(jsonable(value))
        else:
            out[name] = value
    return out


if __name__ == "__main__":
    main()
