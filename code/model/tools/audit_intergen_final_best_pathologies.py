#!/usr/bin/env python3
"""Audit policy/distribution pathologies for the one-market final best.

This is a fixed-theta diagnostic. It does not recalibrate. It re-solves the
saved global-DE best under numerical/menu variants, writes the standard
diagnostic packet for each variant, and records targeted metrics for owner-entry
policy monotonicity, asset-distribution concentration, and owner-rung mass.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
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
DEFAULT_OUTDIR = ROOT / "output/model/intergen_final_best_pathology_audit_20260609"

TARGET_SET = "candidate_no_timing_v0"
TARGET_KEYS = [
    "tfr",
    "childless_rate",
    "own_rate",
    "own_family_gap",
    "housing_increment_0to1",
    "housing_increment_1to2",
    "young_liquid_wealth_to_income",
    "old_age_own_rate",
    "old_age_parent_childless_gap",
    "liquid_wealth_to_income",
    "housing_user_cost_share",
    "prime_childless_renter_median_rooms",
    "prime_childless_owner_median_rooms",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--quick", action="store_true", help="Run only baseline and Nb=120 same-ladder tests")
    parser.add_argument("--skip-diagnostics", action="store_true", help="Skip standard plot packet per variant")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(args.summary.read_text())
    source_record = summary["best"]
    theta = dict(source_record["theta"])
    targets, weights = get_target_set(TARGET_SET)

    variants = [
        {
            "name": "baseline_nb60_cluster_ladder",
            "Nb": 60,
            "H_own": None,
            "note": "Exact cluster numerical/menu settings: H_own = linspace(2, 10, 6).",
        },
        {
            "name": "fine_nb120_cluster_ladder",
            "Nb": 120,
            "H_own": None,
            "note": "Same owner ladder, finer liquid-asset grid.",
        },
    ]
    if not args.quick:
        variants.extend(
            [
                {
                    "name": "baseline_nb60_room_ladder",
                    "Nb": 60,
                    "H_own": [2.0, 4.0, 6.0, 8.0, 9.5, 11.0],
                    "note": "Same asset grid, deliberate owner ladder containing 6 rooms.",
                },
                {
                    "name": "fine_nb120_room_ladder",
                    "Nb": 120,
                    "H_own": [2.0, 4.0, 6.0, 8.0, 9.5, 11.0],
                    "note": "Finer asset grid plus deliberate owner ladder containing 6 rooms.",
                },
            ]
        )

    records: list[dict[str, Any]] = []
    arrays: dict[str, dict[str, Any]] = {}
    for spec in variants:
        record, array_payload = solve_variant(spec, theta, targets, weights, args.outdir, not args.skip_diagnostics)
        records.append(record)
        arrays[spec["name"]] = array_payload
        write_json(args.outdir / "latest_completed_variant.json", record)
        write_variant_checkpoint(args.outdir, records)
        print(
            f"{spec['name']}: loss={record['rank_loss']:.4f}, "
            f"resid={record['market_residual']:.2e}, elapsed={record['elapsed_sec']:.1f}s",
            flush=True,
        )

    write_csv(args.outdir / "variant_summary.csv", records)
    write_json(args.outdir / "variant_summary.json", {"source_record": source_record, "variants": records})
    write_comparison_plots(args.outdir, records, arrays)
    write_note(args.outdir, source_record, records)


def solve_variant(
    spec: dict[str, Any],
    theta: dict[str, float],
    targets: dict[str, float],
    weights: dict[str, float],
    outdir: Path,
    write_standard_diagnostics: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    t0 = time.perf_counter()
    overrides = {
        **base_overrides(J=16, Nb=int(spec["Nb"]), n_house=6, max_iter_eq=25),
        **income_process_overrides(5),
        **theta,
    }
    if spec.get("H_own") is not None:
        H_own = np.asarray(spec["H_own"], dtype=float)
        overrides["H_own"] = H_own
        overrides["n_house"] = int(len(H_own))

    sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
    elapsed = time.perf_counter() - t0
    moments = extract_moments(sol, P)
    rank_loss = diagnostic_loss(moments, targets=targets, weights=weights)
    variant_dir = outdir / str(spec["name"])
    variant_dir.mkdir(parents=True, exist_ok=True)
    if write_standard_diagnostics:
        write_diagnostics(sol, P, variant_dir / "diagnostics")

    policy_metrics, owner_policy_payload = owner_entry_policy_metrics(sol, P)
    distribution_metrics, age_distribution_payload = asset_distribution_metrics(sol, P)
    rung_metrics = owner_rung_metrics(sol, P)
    lifecycle_metrics = lifecycle_metrics_from_solution(sol, P, moments)
    metrics = {
        **policy_metrics,
        **distribution_metrics,
        **rung_metrics,
        **lifecycle_metrics,
        "renter25_45_all_cap_share": float(moments.get("renter25_45_all_cap_share", np.nan)),
        "owner_neg_liquid_share_2534": float(moments.get("owner_neg_liquid_share_2534", np.nan)),
    }
    write_json(variant_dir / "metrics.json", metrics)
    write_csv(variant_dir / "asset_concentration_by_age.csv", age_distribution_payload["rows"])
    plot_asset_concentration(variant_dir / "asset_concentration_by_age.png", age_distribution_payload)
    plot_owner_entry_audit(variant_dir / "owner_entry_probability_audit.png", owner_policy_payload)

    record = {
        "variant": str(spec["name"]),
        "note": str(spec["note"]),
        "elapsed_sec": float(elapsed),
        "rank_loss": float(rank_loss),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "p_eq": float(np.asarray(p_eq).reshape(-1)[0]),
        "Nb": int(P.Nb),
        "H_own": jsonable(np.asarray(P.H_own, dtype=float)),
        "hR_max": float(P.hR_max),
        "moments": {k: float(moments.get(k, np.nan)) for k in TARGET_KEYS},
        "extra_moments": {
            "mean_age_first_birth": float(moments.get("mean_age_first_birth", np.nan)),
            "aggregate_own_rate": float(moments.get("aggregate_own_rate", np.nan)),
            "own_rate_2534": float(moments.get("own_rate_2534", np.nan)),
            "own_rate_3544": float(moments.get("own_rate_3544", np.nan)),
            "young_owner_rate": float(moments.get("young_owner_rate", np.nan)),
            "old_owner_rate": float(moments.get("old_owner_rate", np.nan)),
        },
        "metrics": metrics,
        "variant_dir": str(variant_dir),
    }
    write_json(variant_dir / "record.json", record)
    arrays = {
        "ages": np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float),
        "own_by_age": np.asarray(sol.own_by_age, dtype=float),
        "fert_by_age": np.asarray(getattr(sol, "fert_by_age", np.zeros(P.J)), dtype=float),
        "asset_rows": age_distribution_payload["rows"],
    }
    return record, arrays


def owner_entry_policy_metrics(sol: Any, P: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    b_grid = np.asarray(sol.b_grid, dtype=float)
    z_grid = np.asarray(getattr(sol, "type_values", getattr(P, "z_grid", [1.0])), dtype=float).reshape(-1)
    tp = getattr(sol, "tenure_probs", None)
    tc = np.asarray(sol.tenure_choice)
    V = np.asarray(getattr(sol, "V", np.empty_like(tc, dtype=float)))
    rows = []
    worst = None
    total_pairs = 0
    any_negative_pairs = 0
    material_pairs = 0
    for j in range(P.J):
        age = float(P.age_start + j * P.da)
        for zz, z_value in enumerate(z_grid):
            if tp is None:
                line = (tc[:, 0, 0, j, zz, 0, 0] > 0).astype(float)
            else:
                line = np.sum(np.asarray(tp)[:, 0, 0, j, zz, 0, 0, 1:], axis=1)
            valid = np.isfinite(line)
            if V.ndim == 7:
                valid &= V[:, 0, 0, j, zz, 0, 0] > -1e9
            line = np.asarray(line[valid], dtype=float)
            b_valid = b_grid[valid]
            if line.size < 3:
                continue
            diff = np.diff(line)
            neg = diff[diff < -1e-8]
            material = diff[diff < -0.025]
            total_pairs += 1
            if neg.size:
                any_negative_pairs += 1
            if material.size:
                material_pairs += 1
            max_drop = float(np.min(diff)) if diff.size else 0.0
            total_downward = float(-np.sum(neg)) if neg.size else 0.0
            total_variation = float(np.sum(np.abs(diff)))
            net_change = float(abs(line[-1] - line[0]))
            tv_ratio = total_variation / max(net_change, 1e-8)
            row = {
                "age": age,
                "z": float(z_value),
                "max_drop": max_drop,
                "total_downward_variation": total_downward,
                "n_negative_diffs": int(neg.size),
                "n_material_drops": int(material.size),
                "total_variation": total_variation,
                "net_change": net_change,
                "tv_ratio": tv_ratio,
                "min_prob": float(np.min(line)),
                "max_prob": float(np.max(line)),
                "prob_at_lowest_wealth": float(line[0]),
                "prob_at_highest_wealth": float(line[-1]),
            }
            rows.append(row)
            if worst is None or row["max_drop"] < worst["max_drop"]:
                worst = row

    metrics = {
        "owner_entry_pairs_checked": int(total_pairs),
        "owner_entry_pairs_with_any_negative_diff": int(any_negative_pairs),
        "owner_entry_pairs_with_material_drop_gt_0p025": int(material_pairs),
        "owner_entry_share_pairs_with_material_drop": float(material_pairs / max(total_pairs, 1)),
        "owner_entry_worst_max_drop": float(worst["max_drop"]) if worst else np.nan,
        "owner_entry_worst_age": float(worst["age"]) if worst else np.nan,
        "owner_entry_worst_z": float(worst["z"]) if worst else np.nan,
        "owner_entry_max_total_downward_variation": float(max((r["total_downward_variation"] for r in rows), default=np.nan)),
        "owner_entry_max_tv_ratio": float(max((r["tv_ratio"] for r in rows), default=np.nan)),
    }
    payload = {"b_grid": b_grid, "z_grid": z_grid, "rows": rows, "sol": sol, "P": P}
    return metrics, payload


def asset_distribution_metrics(sol: Any, P: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    g = np.asarray(sol.g, dtype=float)
    b_grid = np.asarray(sol.b_grid, dtype=float).reshape(-1)
    ages = np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float)
    mass_total = float(np.sum(g))
    mass_b = np.sum(g, axis=tuple(range(1, g.ndim)))
    overall = concentration_summary(mass_b, b_grid, mass_total)
    rows = []
    childless_rows = []
    for j, age in enumerate(ages):
        gj = g[:, :, :, j, :, :, :]
        mass_age = float(np.sum(gj))
        mb = np.sum(gj, axis=tuple(range(1, gj.ndim)))
        row = {"age": float(age), **concentration_summary(mb, b_grid, mass_age)}
        rows.append(row)
        cr = g[:, 0, 0, j, :, 0, 0]
        cmass = float(np.sum(cr))
        childless_rows.append({"age": float(age), **concentration_summary(np.sum(cr, axis=1), b_grid, cmass)})
    max_age_row = max(rows, key=lambda r: finite_or_neg_inf(r["max_asset_bin_share"]))
    min_eff_row = min(rows, key=lambda r: finite_or_pos_inf(r["effective_asset_bins"]))
    max_childless_row = max(childless_rows, key=lambda r: finite_or_neg_inf(r["max_asset_bin_share"]))
    metrics = {
        "asset_overall_first_grid_mass": overall["first_grid_mass"],
        "asset_overall_near_zero_mass_abs_le_0p5": overall["near_zero_mass_abs_le_0p5"],
        "asset_overall_negative_mass": overall["negative_mass"],
        "asset_overall_max_bin_share": overall["max_asset_bin_share"],
        "asset_overall_effective_bins": overall["effective_asset_bins"],
        "asset_age_max_bin_share": max_age_row["max_asset_bin_share"],
        "asset_age_of_max_bin_share": max_age_row["age"],
        "asset_age_min_effective_bins": min_eff_row["effective_asset_bins"],
        "asset_age_of_min_effective_bins": min_eff_row["age"],
        "childless_renter_age_max_bin_share": max_childless_row["max_asset_bin_share"],
        "childless_renter_age_of_max_bin_share": max_childless_row["age"],
    }
    payload = {"b_grid": b_grid, "ages": ages, "rows": rows, "childless_rows": childless_rows}
    return metrics, payload


def concentration_summary(mass_b: np.ndarray, b_grid: np.ndarray, total_mass: float) -> dict[str, float]:
    if total_mass <= 1e-14 or not np.isfinite(total_mass):
        return {
            "mass": float(total_mass),
            "first_grid_mass": np.nan,
            "near_zero_mass_abs_le_0p5": np.nan,
            "negative_mass": np.nan,
            "max_asset_bin_share": np.nan,
            "effective_asset_bins": np.nan,
        }
    share = np.asarray(mass_b, dtype=float) / total_mass
    return {
        "mass": float(total_mass),
        "first_grid_mass": float(share[0]),
        "near_zero_mass_abs_le_0p5": float(np.sum(share[np.abs(b_grid) <= 0.5])),
        "negative_mass": float(np.sum(share[b_grid < 0.0])),
        "max_asset_bin_share": float(np.max(share)),
        "effective_asset_bins": float(1.0 / max(float(np.sum(share**2)), 1e-14)),
    }


def owner_rung_metrics(sol: Any, P: Any) -> dict[str, float]:
    demand = np.asarray(getattr(sol, "owner_demand_by_size", np.zeros(P.n_house)), dtype=float).reshape(-1)
    total = float(np.sum(demand))
    if total <= 1e-14:
        return {
            "owner_rung_max_service_share": np.nan,
            "owner_rung_hhi": np.nan,
            "owner_rung_top_h": np.nan,
        }
    share = demand / total
    idx = int(np.argmax(share))
    return {
        "owner_rung_max_service_share": float(share[idx]),
        "owner_rung_hhi": float(np.sum(share**2)),
        "owner_rung_top_h": float(np.asarray(P.H_own, dtype=float)[idx]),
    }


def lifecycle_metrics_from_solution(sol: Any, P: Any, moments: dict[str, float]) -> dict[str, float]:
    own = np.asarray(getattr(sol, "own_by_age", np.zeros(P.J)), dtype=float)
    fert = np.asarray(getattr(sol, "fert_by_age", np.zeros(P.J)), dtype=float)
    ages = np.asarray(P.age_start + np.arange(P.J) * P.da, dtype=float)
    own_diff = np.diff(own)
    fert_mass = float(np.sum(fert))
    fertility_mean_age_from_profile = float(np.sum(ages * fert) / max(fert_mass, 1e-14))
    return {
        "own_rate_2534": float(moments.get("own_rate_2534", np.nan)),
        "own_rate_3544": float(moments.get("own_rate_3544", np.nan)),
        "old_minus_prime_ownership": float(moments.get("old_age_own_rate", np.nan) - moments.get("own_rate", np.nan)),
        "max_age_to_age_ownership_jump": float(np.max(own_diff)) if own_diff.size else np.nan,
        "age_of_max_ownership_jump_start": float(ages[int(np.argmax(own_diff))]) if own_diff.size else np.nan,
        "fertility_profile_mean_age": fertility_mean_age_from_profile,
        "fertility_profile_peak_age": float(ages[int(np.argmax(fert))]) if fert.size else np.nan,
        "fertility_profile_peak_births": float(np.max(fert)) if fert.size else np.nan,
    }


def plot_asset_concentration(path: Path, payload: dict[str, Any]) -> None:
    rows = payload["rows"]
    ages = [r["age"] for r in rows]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True)
    axes[0].plot(ages, [r["max_asset_bin_share"] for r in rows], marker="o", lw=1.6)
    axes[0].set_ylabel("max asset-bin share")
    axes[0].grid(alpha=0.2)
    axes[1].plot(ages, [r["effective_asset_bins"] for r in rows], marker="o", lw=1.6)
    axes[1].set_ylabel("effective asset bins")
    axes[1].set_xlabel("age")
    axes[1].grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_owner_entry_audit(path: Path, payload: dict[str, Any]) -> None:
    sol = payload["sol"]
    P = payload["P"]
    b_grid = payload["b_grid"]
    z_grid = payload["z_grid"]
    tp = getattr(sol, "tenure_probs", None)
    tc = np.asarray(sol.tenure_choice)
    ages_to_plot = [30.0, 42.0]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharey=True)
    for ax, age in zip(axes, ages_to_plot):
        j = int(np.clip(round((age - float(P.age_start)) / float(P.da)), 0, P.J - 1))
        for zz, z_value in enumerate(z_grid):
            if tp is None:
                line = (tc[:, 0, 0, j, zz, 0, 0] > 0).astype(float)
            else:
                line = np.sum(np.asarray(tp)[:, 0, 0, j, zz, 0, 0, 1:], axis=1)
            ax.plot(b_grid, line, lw=1.7, label=f"z={z_value:g}")
        ax.set_title(f"age {P.age_start + j * P.da:g}")
        ax.set_xlabel("liquid wealth")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("owner-entry probability")
    axes[0].legend(frameon=False, ncols=1)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_comparison_plots(outdir: Path, records: list[dict[str, Any]], arrays: dict[str, dict[str, Any]]) -> None:
    figdir = outdir / "comparison_figures"
    figdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for rec in records:
        arr = arrays[rec["variant"]]
        ax.plot(arr["ages"], arr["own_by_age"], lw=1.8, label=rec["variant"])
    ax.set_xlabel("age")
    ax.set_ylabel("ownership rate")
    ax.set_ylim(0.0, 1.05)
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "ownership_by_age_compare.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for rec in records:
        arr = arrays[rec["variant"]]
        ax.plot(arr["ages"], arr["fert_by_age"], lw=1.8, label=rec["variant"])
    ax.set_xlabel("age")
    ax.set_ylabel("expected births")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "fertility_by_age_compare.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 6.0), sharex=True)
    for rec in records:
        rows = arrays[rec["variant"]]["asset_rows"]
        ages = [r["age"] for r in rows]
        axes[0].plot(ages, [r["max_asset_bin_share"] for r in rows], lw=1.7, label=rec["variant"])
        axes[1].plot(ages, [r["effective_asset_bins"] for r in rows], lw=1.7, label=rec["variant"])
    axes[0].set_ylabel("max asset-bin share")
    axes[1].set_ylabel("effective asset bins")
    axes[1].set_xlabel("age")
    for ax in axes:
        ax.grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "asset_concentration_compare.png", dpi=180)
    plt.close(fig)

    labels = [r["variant"].replace("_", "\n") for r in records]
    metrics = [
        ("rank_loss", [r["rank_loss"] for r in records]),
        ("own_rate", [r["moments"]["own_rate"] for r in records]),
        ("owner_entry material-drop share", [r["metrics"]["owner_entry_share_pairs_with_material_drop"] for r in records]),
        ("max asset-bin share", [r["metrics"]["asset_age_max_bin_share"] for r in records]),
        ("top owner-rung service share", [r["metrics"]["owner_rung_max_service_share"] for r in records]),
    ]
    fig, axes = plt.subplots(len(metrics), 1, figsize=(8.0, 11.0))
    for ax, (title, values) in zip(axes, metrics):
        ax.bar(np.arange(len(records)), values)
        ax.set_title(title)
        ax.set_xticks(np.arange(len(records)), labels, rotation=0, fontsize=7)
        ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(figdir / "pathology_metric_bars.png", dpi=180)
    plt.close(fig)


def write_note(outdir: Path, source_record: dict[str, Any], records: list[dict[str, Any]]) -> None:
    lines = [
        "# Final-Best Pathology Audit",
        "",
        "This fixed-theta audit re-solves the saved global-DE best under numerical and owner-menu variants. It is not a recalibration.",
        "",
        f"Source best: run `{source_record.get('run_dir')}`, task `{source_record.get('task_id')}`, case `{source_record.get('case')}`, label `{source_record.get('label')}`.",
        f"Source rank loss: `{float(source_record.get('rank_loss')):.6g}`.",
        "",
        "## Variants",
        "",
        "| Variant | Loss | Prime own | Old own | Housing share | Max asset-bin share | Effective bins min age | Owner-entry material-drop share | Top owner-rung share |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for rec in records:
        m = rec["moments"]
        x = rec["metrics"]
        lines.append(
            "| "
            + " | ".join(
                [
                    rec["variant"],
                    f"{rec['rank_loss']:.3f}",
                    f"{m['own_rate']:.3f}",
                    f"{m['old_age_own_rate']:.3f}",
                    f"{m['housing_user_cost_share']:.3f}",
                    f"{x['asset_age_max_bin_share']:.3f}",
                    f"{x['asset_age_min_effective_bins']:.2f}",
                    f"{x['owner_entry_share_pairs_with_material_drop']:.3f}",
                    f"{x['owner_rung_max_service_share']:.3f}",
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Reading Guide",
            "",
            "- If `Nb=120` materially smooths policies or distributions while moments remain close, the pathology is partly numerical/grid-driven.",
            "- If the deliberate owner ladder changes owner-size concentration or median owner rooms, the temporary ladder is contaminating fit and diagnostics.",
            "- If ownership timing and low-wealth spikes survive both changes, the failure is primarily economic/objective-based rather than only numerical.",
            "",
            "Generated comparison plots are under `comparison_figures/`; per-variant standard equilibrium plots are under each variant's `diagnostics/` directory.",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


def write_variant_checkpoint(outdir: Path, records: list[dict[str, Any]]) -> None:
    write_json(outdir / "variant_summary_partial.json", {"variants": records})
    write_csv(outdir / "variant_summary_partial.csv", records)


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(obj), indent=2, sort_keys=True))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    flat_rows = [flatten_row(r) for r in rows]
    fields: list[str] = []
    for row in flat_rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat_rows)


def flatten_row(row: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in row.items():
        name = f"{prefix}{key}"
        if isinstance(value, dict):
            out.update(flatten_row(value, prefix=f"{name}."))
        elif isinstance(value, (list, tuple, np.ndarray)):
            out[name] = json.dumps(jsonable(value))
        else:
            out[name] = value
    return out


def finite_or_neg_inf(value: Any) -> float:
    try:
        x = float(value)
    except Exception:
        return -math.inf
    return x if math.isfinite(x) else -math.inf


def finite_or_pos_inf(value: Any) -> float:
    try:
        x = float(value)
    except Exception:
        return math.inf
    return x if math.isfinite(x) else math.inf


if __name__ == "__main__":
    main()
