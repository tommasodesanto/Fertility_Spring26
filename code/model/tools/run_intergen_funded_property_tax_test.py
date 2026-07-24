#!/usr/bin/env python3
"""Fixed-parameter testing battery with a balanced property-tax budget."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from intergen_housing_fertility_optimized.calibration import diagnostic_loss, extract_moments
from intergen_housing_fertility_optimized.m5_profile import M5_THETA, m5_overrides, m5_target_system
from intergen_housing_fertility_optimized.new_moment_profile import (
    NEW_MOMENT_PROFILE_NAME,
    new_moment_overrides,
    new_moment_target_system,
)
from intergen_housing_fertility_optimized.solver import run_model_cp_dt


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = ROOT / "output/model/intergen_funded_property_tax_test_20260721"
DEFAULT_NEW_MOMENT_SOURCE = (
    ROOT
    / "output/model/intergen_new_moment_unrestricted_overnight_20260723_w2/report/results.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--smoke", action="store_true", help="Use Nb=40 and loose equilibrium tolerances.")
    parser.add_argument(
        "--profile",
        choices=("m5", "new-moment"),
        default="m5",
        help="Fixed parameter vector and target contract to use.",
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=DEFAULT_NEW_MOMENT_SOURCE,
        help="Results JSON containing selected.theta for --profile new-moment.",
    )
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def bracketed_root(function, lower: float, upper: float, *, tolerance: float, max_iterations: int) -> float:
    f_lower = float(function(lower))
    f_upper = float(function(upper))
    if f_lower == 0.0:
        return float(lower)
    if f_upper == 0.0:
        return float(upper)
    if f_lower * f_upper > 0.0:
        raise ValueError("root solver requires a sign-changing bracket")
    midpoint = 0.5 * (lower + upper)
    for _ in range(max_iterations):
        midpoint = upper - f_upper * (upper - lower) / (f_upper - f_lower)
        if not lower < midpoint < upper or min(midpoint - lower, upper - midpoint) < 1e-8 * (upper - lower):
            midpoint = 0.5 * (lower + upper)
        f_midpoint = float(function(midpoint))
        if abs(f_midpoint) <= tolerance or (upper - lower) <= tolerance:
            return float(midpoint)
        if f_lower * f_midpoint <= 0.0:
            upper = midpoint
            f_upper = f_midpoint
        else:
            lower = midpoint
            f_lower = f_midpoint
    return float(midpoint)


def load_profile(
    profile: str,
    source: Path,
    smoke: bool,
) -> tuple[str, dict[str, float], dict[str, Any], Any, str]:
    if profile == "m5":
        return (
            "M5",
            dict(M5_THETA),
            m5_overrides(tight=not smoke, optimized=True),
            m5_target_system(),
            "fixed-M5 funded property-tax testing battery",
        )
    payload = json.loads(source.read_text(encoding="utf-8"))
    if payload.get("profile") != NEW_MOMENT_PROFILE_NAME:
        raise ValueError(
            f"Expected profile {NEW_MOMENT_PROFILE_NAME!r}, found {payload.get('profile')!r}."
        )
    selected = payload.get("selected")
    if not isinstance(selected, dict) or selected.get("status") != "ok":
        raise ValueError(f"Source does not contain a successful selected result: {source}")
    raw_theta = selected.get("theta")
    if not isinstance(raw_theta, dict) or not raw_theta:
        raise ValueError(f"Selected result does not contain theta: {source}")
    theta = {str(key): float(value) for key, value in raw_theta.items()}
    overrides = new_moment_overrides(tight=not smoke, optimized=True)
    overrides.update(theta)
    return (
        "current M",
        theta,
        overrides,
        new_moment_target_system(),
        "current-M funded property-tax testing battery",
    )


def case_overrides(
    base_overrides: dict[str, Any],
    annual_tax: float,
    grant: bool,
    transfer: float,
    smoke: bool,
) -> dict[str, Any]:
    overrides = dict(base_overrides)
    overrides["tau_H"] = float(annual_tax) * float(overrides["period_years"])
    overrides["property_tax_lump_sum_transfer"] = float(transfer)
    if smoke:
        overrides.update(Nb=40, max_iter_eq=10, tol_eq=1e-4)
    if grant:
        overrides.update(
            birth_entry_grant=True,
            birth_entry_grant_amount=0.4,
            birth_entry_grant_locations=np.array([], dtype=int),
            birth_entry_grant_owner_rungs=np.array([3, 4, 5], dtype=int),
        )
    else:
        overrides["birth_entry_grant"] = False
        overrides["birth_entry_grant_amount"] = 0.0
    return overrides


def solve_balanced_case(
    label: str,
    annual_tax: float,
    grant: bool,
    smoke: bool,
    base_overrides: dict[str, Any],
    target_system: Any,
) -> tuple[dict[str, Any], Any, Any]:
    cache: dict[float, tuple[Any, Any, np.ndarray]] = {}

    def solve_at(transfer: float) -> tuple[Any, Any, np.ndarray]:
        key = float(round(transfer, 12))
        if key not in cache:
            overrides = case_overrides(base_overrides, annual_tax, grant, transfer, smoke)
            cache[key] = run_model_cp_dt(overrides, verbose=False)
        return cache[key]

    def residual(transfer: float) -> float:
        solution, _, _ = solve_at(transfer)
        return float(solution.property_tax_budget_residual)

    lower = 0.0
    upper = 0.5
    f_lower = residual(lower)
    f_upper = residual(upper)
    while f_upper > 0.0 and upper < 8.0:
        upper *= 2.0
        f_upper = residual(upper)
    if f_lower < 0.0 or f_upper > 0.0:
        raise RuntimeError(
            f"could not bracket fiscal transfer for {label}: "
            f"F(0)={f_lower:.6g}, F({upper:g})={f_upper:.6g}"
        )
    transfer = bracketed_root(residual, lower, upper, tolerance=2.5e-5, max_iterations=30)
    solution, parameters, price = solve_at(transfer)
    moments = extract_moments(solution, parameters)
    loss = diagnostic_loss(
        moments,
        targets=target_system.targets_dict(),
        weights=target_system.weights_dict(),
    )
    result = {
        "case": label,
        "annual_property_tax_rate": annual_tax,
        "purchase_grant": grant,
        "lump_sum_transfer_period_units": transfer,
        "lump_sum_transfer_annual_income_units": transfer / float(parameters.period_years),
        "property_tax_revenue": float(solution.property_tax_revenue),
        "grant_recipient_mass": float(solution.birth_entry_grant_recipient_mass),
        "grant_outlays": float(solution.birth_entry_grant_outlays),
        "transfer_outlays": float(solution.property_tax_transfer_outlays),
        "government_budget_residual": float(solution.property_tax_budget_residual),
        "price": float(np.asarray(price).reshape(-1)[0]),
        "tfr": float(moments["tfr"]),
        "childless_rate": float(moments["childless_rate"]),
        "own_rate": float(moments["own_rate"]),
        "old_age_own_rate": float(moments["old_age_own_rate"]),
        "housing_increment_0to1": float(moments["housing_increment_0to1"]),
        "fixed_theta_loss": float(loss),
        "market_residual": float(solution.best_max_abs_rel_excess),
        "strict_converged": bool(getattr(solution, "timings", {}).get("strict_converged", False)),
        "fiscal_root_evaluations": len(cache),
    }
    return result, solution, parameters


def make_plot(outdir: Path, rows: list[dict[str, Any]]) -> None:
    baseline = rows[0]
    policies = rows[1:]
    labels = ["Rebated 2% tax", "Rebated 2% + grant"]
    y = np.arange(len(policies))[::-1]
    price_changes = np.asarray(
        [100.0 * (row["price"] / baseline["price"] - 1.0) for row in policies]
    )
    tfr_changes = np.asarray([row["tfr"] - baseline["tfr"] for row in policies])
    rebate_shares = np.asarray(
        [100.0 * row["transfer_outlays"] / row["property_tax_revenue"] for row in policies]
    )
    grant_shares = np.asarray(
        [100.0 * row["grant_outlays"] / row["property_tax_revenue"] for row in policies]
    )

    plt.rcdefaults()
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "legend.fontsize": 10,
        }
    )
    navy = "#17365D"
    orange = "#C45A3B"
    green = "#2B6F68"
    gray = "#A8ADB3"
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(13.0, 5.2),
        gridspec_kw={"width_ratios": [1.18, 1.02, 1.02]},
    )

    axes[0].axvline(0.0, color="#4B4B4B", linewidth=0.7)
    axes[0].hlines(y, 0.0, tfr_changes, color=gray, linewidth=1.3)
    axes[0].scatter(tfr_changes, y, color=navy, s=45, zorder=3)
    axes[0].set_yticks(y, labels)
    axes[0].set_xlim(-0.001, max(0.014, 1.18 * float(np.max(tfr_changes))))
    axes[0].set_xlabel("Change in children per woman")
    axes[0].set_title("A. Completed fertility", loc="left", fontweight="bold")
    for yy, value in zip(y, tfr_changes):
        axes[0].annotate(
            f"{value:+.3f}",
            (value, yy),
            xytext=(6, 0),
            textcoords="offset points",
            va="center",
        )

    axes[1].axvline(0.0, color="#4B4B4B", linewidth=0.7)
    axes[1].hlines(y, 0.0, price_changes, color=gray, linewidth=1.3)
    axes[1].scatter(price_changes, y, color=orange, marker="s", s=45, zorder=3)
    axes[1].set_yticks(y, ["", ""])
    axes[1].set_xlim(min(-20.0, 1.15 * float(np.min(price_changes))), 1.5)
    axes[1].set_xlabel("Percent change")
    axes[1].set_title("B. Equilibrium house price", loc="left", fontweight="bold")
    for yy, value in zip(y, price_changes):
        axes[1].annotate(
            f"{value:+.1f}%",
            (value, yy),
            xytext=(-6, 0),
            textcoords="offset points",
            ha="right",
            va="center",
        )

    axes[2].barh(y, rebate_shares, color=green, height=0.42, label="Universal rebate")
    axes[2].barh(
        y,
        grant_shares,
        left=rebate_shares,
        color=orange,
        height=0.42,
        label="Targeted grant",
    )
    axes[2].set_yticks(y, ["", ""])
    axes[2].set_xlim(0.0, 100.0)
    axes[2].set_xlabel("Share of property-tax revenue (%)")
    axes[2].set_title("C. Fiscal allocation", loc="left", fontweight="bold")
    handles, legend_labels = axes[2].get_legend_handles_labels()

    for axis in axes:
        axis.grid(axis="x", color="#D9D9D9", linewidth=0.5)
        axis.set_axisbelow(True)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0)
    note = (
        "Notes: Fixed current-M parameters and normalized stationary population. "
        "All changes are relative to the fiscally closed, rebated 1% property-tax baseline. "
        "The 0.4 grant applies to realized renter-to-owner purchases of homes with 6+ rooms; "
        "residual revenue is rebated equally."
    )
    fig.subplots_adjust(left=0.14, right=0.985, top=0.82, bottom=0.23, wspace=0.25)
    fig.legend(
        handles,
        legend_labels,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.98),
        ncol=2,
    )
    fig.text(0.14, 0.07, note, ha="left", va="bottom", fontsize=8.3, color="#333333")
    for filename in (
        "policy_general_equilibrium_feedback_nb120.png",
        "funded_policy_comparison.png",
    ):
        fig.savefig(outdir / filename, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    profile_label, theta, base_overrides, target_system, purpose = load_profile(
        args.profile,
        args.source.resolve(),
        args.smoke,
    )
    cases = (
        ("rebated_tax1_baseline", 0.01, False),
        ("rebated_tax2", 0.02, False),
        ("rebated_tax2_grant0p4_Hge6", 0.02, True),
    )
    metadata = {
        "status": "running",
        "purpose": purpose,
        "profile": profile_label,
        "source": str(args.source.resolve()) if args.profile == "new-moment" else None,
        "fiscal_rule": "property-tax revenue = equal lump-sum transfers + targeted purchase-grant outlays",
        "tax_base": "all occupied housing, consistent with the implemented owner and rental user cost",
        "theta": theta,
        "smoke": bool(args.smoke),
    }
    (outdir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    summaries: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    for index, (label, annual_tax, grant) in enumerate(cases, start=1):
        case_start = time.perf_counter()
        print(f"case {index}/{len(cases)}: {label}", flush=True)
        summary, solution, parameters = solve_balanced_case(
            label,
            annual_tax,
            grant,
            args.smoke,
            base_overrides,
            target_system,
        )
        summary["elapsed_seconds"] = time.perf_counter() - case_start
        summaries.append(summary)
        moments = extract_moments(solution, parameters)
        for name, target, weight in zip(
            target_system.moment_names,
            target_system.target_values,
            target_system.weights,
        ):
            model = float(moments[name])
            gap = model - float(target)
            target_rows.append(
                {
                    "case": label,
                    "moment": name,
                    "target": float(target),
                    "model": model,
                    "gap": gap,
                    "weight": float(weight),
                    "loss_contribution": float(weight) * gap * gap,
                }
            )
        write_csv(outdir / "summary.csv", summaries)
        write_csv(outdir / "target_fit_long.csv", target_rows)
        (outdir / "latest_completed_case.json").write_text(
            json.dumps(jsonable(summary), indent=2, sort_keys=True) + "\n"
        )
        print(
            f"completed {label}: T={summary['lump_sum_transfer_period_units']:.5f}, "
            f"p={summary['price']:.6f}, TFR={summary['tfr']:.5f}, "
            f"budget={summary['government_budget_residual']:.2e}, "
            f"elapsed={summary['elapsed_seconds']:.1f}s",
            flush=True,
        )
    baseline = summaries[0]
    for row in summaries:
        row["price_change_percent_from_funded_baseline"] = 100.0 * (row["price"] / baseline["price"] - 1.0)
        row["tfr_change_percent_from_funded_baseline"] = 100.0 * (row["tfr"] / baseline["tfr"] - 1.0)
        row["tfr_level_change_from_funded_baseline"] = row["tfr"] - baseline["tfr"]
    write_csv(outdir / "summary.csv", summaries)
    make_plot(outdir, summaries)
    metadata.update(status="complete", elapsed_seconds=time.perf_counter() - started)
    (outdir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    lines = [
        "# Funded property-tax testing battery",
        "",
        f"Fixed {profile_label} parameters; the fiscal transfer is solved jointly with each stationary equilibrium.",
        "This is a testing-mode, fixed-population comparison, not a recalibration or the paper's population-adjusted estimand.",
        "",
        "| Case | TFR | TFR change | Price change | Period transfer | Tax revenue | Grant outlays | Budget residual |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summaries:
        lines.append(
            f"| {row['case']} | {row['tfr']:.4f} | {row['tfr_change_percent_from_funded_baseline']:+.2f}% "
            f"| {row['price_change_percent_from_funded_baseline']:+.2f}% | {row['lump_sum_transfer_period_units']:.4f} "
            f"| {row['property_tax_revenue']:.4f} | {row['grant_outlays']:.4f} "
            f"| {row['government_budget_residual']:.2e} |"
        )
    lines.extend(
        [
            "",
            "The property-tax base includes rental and owner housing because the implemented tax enters both user costs.",
            "The targeted grant is paid only on realized renter-to-owner purchases of six or more rooms in the birth state.",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
