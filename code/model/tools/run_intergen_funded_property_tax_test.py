#!/usr/bin/env python3
"""Fixed-M5 testing battery with a balanced property-tax budget."""

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
from intergen_housing_fertility_optimized.solver import run_model_cp_dt


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = ROOT / "output/model/intergen_funded_property_tax_test_20260721"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--smoke", action="store_true", help="Use Nb=40 and loose equilibrium tolerances.")
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


def case_overrides(annual_tax: float, grant: bool, transfer: float, smoke: bool) -> dict[str, Any]:
    overrides = m5_overrides(tight=not smoke, optimized=True)
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
    return overrides


def solve_balanced_case(label: str, annual_tax: float, grant: bool, smoke: bool) -> tuple[dict[str, Any], Any, Any]:
    cache: dict[float, tuple[Any, Any, np.ndarray]] = {}

    def solve_at(transfer: float) -> tuple[Any, Any, np.ndarray]:
        key = float(round(transfer, 12))
        if key not in cache:
            overrides = case_overrides(annual_tax, grant, transfer, smoke)
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
    target_system = m5_target_system()
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
    labels = ["Rebated\n1% baseline", "Rebated\n2% tax", "Rebated 2%\n+ grant"]
    baseline = rows[0]
    price_changes = [100.0 * (row["price"] / baseline["price"] - 1.0) for row in rows]
    tfr_changes = [100.0 * (row["tfr"] / baseline["tfr"] - 1.0) for row in rows]
    transfers = [row["lump_sum_transfer_annual_income_units"] for row in rows]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    for axis, values, title, ylabel in (
        (axes[0], price_changes, "House-price response", "Percent from funded baseline"),
        (axes[1], tfr_changes, "Completed-fertility response", "Percent from funded baseline"),
        (axes[2], transfers, "Universal residual rebate", "Annual-income units"),
    ):
        axis.bar(labels, values, color=["#5B6B7A", "#B45F4A", "#3C7D72"])
        axis.axhline(0.0, color="black", linewidth=0.7)
        axis.set_title(title)
        axis.set_ylabel(ylabel)
        axis.tick_params(axis="x", labelsize=8)
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / "funded_policy_comparison.png", dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    cases = (
        ("rebated_tax1_baseline", 0.01, False),
        ("rebated_tax2", 0.02, False),
        ("rebated_tax2_grant0p4_Hge6", 0.02, True),
    )
    metadata = {
        "status": "running",
        "purpose": "fixed-M5 funded property-tax testing battery",
        "fiscal_rule": "property-tax revenue = equal lump-sum transfers + targeted purchase-grant outlays",
        "tax_base": "all occupied housing, consistent with the implemented owner and rental user cost",
        "theta": M5_THETA,
        "smoke": bool(args.smoke),
    }
    (outdir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    summaries: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    target_system = m5_target_system()
    started = time.perf_counter()
    for index, (label, annual_tax, grant) in enumerate(cases, start=1):
        case_start = time.perf_counter()
        print(f"case {index}/{len(cases)}: {label}", flush=True)
        summary, solution, parameters = solve_balanced_case(label, annual_tax, grant, args.smoke)
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
        "Fixed M5 parameters; the fiscal transfer is solved jointly with each stationary equilibrium.",
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
