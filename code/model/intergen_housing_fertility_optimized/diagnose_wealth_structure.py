#!/usr/bin/env python3
"""Decompose the wealth target at a certified corrected one-shot winner."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .calibration import extract_moments
from .calibration_search import load_start_theta, target_fit_rows
from .new_moment_profile import new_moment_overrides, new_moment_target_system
from .solver import (
    add_aggregate_wealth_bequest_flow_moments,
    add_annual_gross_estate_wealth_moments,
    run_model_cp_dt,
)


def balance_sheet(
    distribution: np.ndarray,
    asset_grid: np.ndarray,
    prices: np.ndarray,
    owner_sizes: np.ndarray,
) -> dict[str, float]:
    """Aggregate a balance sheet using one internally consistent state timing."""

    positive_assets = 0.0
    debt = 0.0
    gross_housing = 0.0
    for location in range(distribution.shape[2]):
        for tenure in range(distribution.shape[1]):
            mass_by_asset = np.sum(
                distribution[:, tenure, location, :, :, :, :],
                axis=(1, 2, 3, 4),
            )
            positive_assets += float(
                np.sum(mass_by_asset * np.maximum(asset_grid, 0.0))
            )
            debt += float(np.sum(mass_by_asset * np.maximum(-asset_grid, 0.0)))
            if tenure > 0:
                gross_housing += (
                    float(np.sum(mass_by_asset))
                    * float(prices[location])
                    * float(owner_sizes[tenure - 1])
                )
    net_wealth = positive_assets - debt + gross_housing
    return {
        "positive_liquid_assets": positive_assets,
        "liquid_debt": debt,
        "gross_housing_wealth": gross_housing,
        "net_wealth": net_wealth,
    }


def net_wealth_per_household_by_age(
    distribution: np.ndarray,
    asset_grid: np.ndarray,
    prices: np.ndarray,
    owner_sizes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mass = np.sum(distribution, axis=(0, 1, 2, 4, 5, 6))
    wealth = np.zeros(distribution.shape[3], dtype=float)
    for age_index in range(distribution.shape[3]):
        for location in range(distribution.shape[2]):
            for tenure in range(distribution.shape[1]):
                mass_by_asset = np.sum(
                    distribution[:, tenure, location, age_index, :, :, :],
                    axis=(1, 2, 3),
                )
                housing_value = (
                    float(prices[location]) * float(owner_sizes[tenure - 1])
                    if tenure > 0
                    else 0.0
                )
                wealth[age_index] += float(
                    np.sum(mass_by_asset * (asset_grid + housing_value))
                )
    return wealth / np.maximum(mass, 1e-12), mass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    system = new_moment_target_system()
    theta = load_start_theta(
        args.source,
        expected_target_fingerprint=system.fingerprint,
    )
    solution, parameters, price = run_model_cp_dt(
        {**new_moment_overrides(tight=True, optimized=True), **theta},
        verbose=False,
    )
    moments = extract_moments(solution, parameters)
    fit = target_fit_rows(moments, system)
    wealth = float(solution.aggregate_wealth)
    liquid = float(solution.aggregate_liquid_net_worth)
    positive_liquid = float(solution.aggregate_positive_liquid_assets)
    liquid_debt = float(solution.aggregate_liquid_debt)
    housing = float(solution.aggregate_gross_housing_wealth)
    earnings = float(solution.aggregate_annual_after_tax_earnings)
    bequests = float(solution.annual_bequest_flow)
    target_ratio = float(
        system.targets_dict()["aggregate_wealth_to_annual_after_tax_earnings"]
    )
    required_wealth = target_ratio * earnings
    asset_grid = np.asarray(solution.b_grid, dtype=float)
    prices = np.asarray(price, dtype=float).reshape(-1)
    owner_sizes = np.asarray(parameters.H_own, dtype=float)
    balance_sheets = {
        "beginning_of_period_consistent": balance_sheet(
            np.asarray(solution.g_beginning_distribution, dtype=float),
            asset_grid,
            prices,
            owner_sizes,
        ),
        "hybrid_current_target": balance_sheet(
            np.asarray(solution.g_beginning_assets_by_current_choice, dtype=float),
            asset_grid,
            prices,
            owner_sizes,
        ),
        "posttransaction_consistent": balance_sheet(
            np.asarray(solution.g, dtype=float),
            asset_grid,
            prices,
            owner_sizes,
        ),
    }
    for values in balance_sheets.values():
        values["wealth_to_annual_after_tax_earnings"] = values["net_wealth"] / earnings
    timing_distributions = {
        "beginning_of_period_consistent": np.asarray(
            solution.g_beginning_distribution,
            dtype=float,
        ),
        "hybrid_current_target": np.asarray(
            solution.g_beginning_assets_by_current_choice,
            dtype=float,
        ),
        "posttransaction_consistent": np.asarray(solution.g, dtype=float),
    }
    timing_target_moments = {}
    for timing, timing_distribution in timing_distributions.items():
        timing_stats = SimpleNamespace()
        add_aggregate_wealth_bequest_flow_moments(
            timing_stats,
            timing_distribution,
            parameters,
            asset_grid,
            prices,
        )
        add_annual_gross_estate_wealth_moments(
            timing_stats,
            timing_distribution,
            parameters,
            asset_grid,
            prices,
        )
        timing_target_moments[timing] = {
            "aggregate_wealth_to_annual_after_tax_earnings": float(
                timing_stats.aggregate_wealth_to_annual_after_tax_earnings
            ),
            "annual_bequest_flow_to_aggregate_wealth": float(
                timing_stats.annual_bequest_flow_to_aggregate_wealth
            ),
            "old_total_estate_wealth_to_annual_income_p90_p50_7684": float(
                timing_stats.old_total_estate_wealth_to_annual_income_p90_p50_7684
            ),
            "old_total_estate_wealth_to_annual_income_median_7684": float(
                timing_stats.old_total_estate_wealth_to_annual_income_median_7684
            ),
        }
    final_age_owner_shares = {
        timing: float(np.sum(distribution[:, 1:, :, -1, :, :, :]))
        / max(float(np.sum(distribution[:, :, :, -1, :, :, :])), 1e-12)
        for timing, distribution in timing_distributions.items()
    }
    summary = {
        "status": "wealth_structure_diagnostic_not_a_calibration",
        "source": str(args.source.resolve()),
        "target_fingerprint": system.fingerprint,
        "theta": theta,
        "strict_loss": float(system.loss(moments)),
        "market_residual": float(solution.best_max_abs_rel_excess),
        "price": float(np.asarray(price).reshape(-1)[0]),
        "wealth_to_annual_after_tax_earnings": wealth / earnings,
        "target_wealth_to_annual_after_tax_earnings": target_ratio,
        "aggregate_wealth": wealth,
        "aggregate_liquid_net_worth": liquid,
        "aggregate_positive_liquid_assets": positive_liquid,
        "aggregate_liquid_debt": liquid_debt,
        "aggregate_gross_housing_wealth": housing,
        "liquid_share_of_wealth": liquid / wealth,
        "gross_housing_share_of_wealth": housing / wealth,
        "positive_liquid_assets_share_of_wealth": positive_liquid / wealth,
        "liquid_debt_share_of_wealth": liquid_debt / wealth,
        "aggregate_annual_after_tax_earnings": earnings,
        "annual_bequest_flow": bequests,
        "target_implied_aggregate_wealth_at_model_earnings": required_wealth,
        "aggregate_wealth_shortfall": required_wealth - wealth,
        "wealth_increase_required_percent": 100.0 * (required_wealth / wealth - 1.0),
        "accounting_check_liquid_plus_housing_minus_total": liquid + housing - wealth,
        "wealth_timing_comparison": balance_sheets,
        "target_moment_timing_comparison": timing_target_moments,
        "final_age_owner_share_by_timing": final_age_owner_shares,
    }
    (args.outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    with (args.outdir / "target_fit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fit[0]))
        writer.writeheader()
        writer.writerows(fit)

    ages = float(parameters.age_start) + float(parameters.da) * np.arange(int(parameters.J))
    mass = np.asarray(solution.population_mass_by_age, dtype=float)
    liquid_age = np.asarray(solution.aggregate_liquid_net_worth_by_age, dtype=float)
    positive_liquid_age = np.asarray(
        solution.aggregate_positive_liquid_assets_by_age,
        dtype=float,
    )
    liquid_debt_age = np.asarray(solution.aggregate_liquid_debt_by_age, dtype=float)
    housing_age = np.asarray(solution.aggregate_gross_housing_wealth_by_age, dtype=float)
    wealth_age = np.asarray(solution.aggregate_wealth_by_age, dtype=float)
    earnings_age = np.asarray(
        solution.aggregate_annual_after_tax_earnings_by_age,
        dtype=float,
    )
    age_rows = [
        {
            "age": float(age),
            "population_mass": float(mass[index]),
            "aggregate_liquid_net_worth": float(liquid_age[index]),
            "aggregate_positive_liquid_assets": float(positive_liquid_age[index]),
            "aggregate_liquid_debt": float(liquid_debt_age[index]),
            "aggregate_gross_housing_wealth": float(housing_age[index]),
            "aggregate_total_wealth": float(wealth_age[index]),
            "aggregate_annual_after_tax_earnings": float(earnings_age[index]),
            "liquid_net_worth_per_household": float(
                liquid_age[index] / max(mass[index], 1e-12)
            ),
            "positive_liquid_assets_per_household": float(
                positive_liquid_age[index] / max(mass[index], 1e-12)
            ),
            "liquid_debt_per_household": float(
                liquid_debt_age[index] / max(mass[index], 1e-12)
            ),
            "gross_housing_wealth_per_household": float(
                housing_age[index] / max(mass[index], 1e-12)
            ),
            "total_wealth_per_household": float(
                wealth_age[index] / max(mass[index], 1e-12)
            ),
            "annual_after_tax_earnings_per_household": float(
                earnings_age[index] / max(mass[index], 1e-12)
            ),
        }
        for index, age in enumerate(ages)
    ]
    with (args.outdir / "wealth_by_age.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(age_rows[0]))
        writer.writeheader()
        writer.writerows(age_rows)

    distribution = np.asarray(
        solution.g_beginning_assets_by_current_choice,
        dtype=float,
    )
    tenure_rows = []
    for j, age in enumerate(ages):
        for tenure in range(distribution.shape[1]):
            mass_by_asset = np.sum(
                distribution[:, tenure, :, j, :, :, :],
                axis=(1, 2, 3, 4),
            )
            tenure_mass = float(np.sum(mass_by_asset))
            housing_value = (
                float(np.asarray(price).reshape(-1)[0]) * float(parameters.H_own[tenure - 1])
                if tenure > 0
                else 0.0
            )
            positive_assets = float(np.sum(mass_by_asset * np.maximum(asset_grid, 0.0)))
            debt = float(np.sum(mass_by_asset * np.maximum(-asset_grid, 0.0)))
            gross_housing = tenure_mass * housing_value
            tenure_rows.append(
                {
                    "age": float(age),
                    "tenure_index": int(tenure),
                    "housing_services": float(parameters.H_own[tenure - 1]) if tenure > 0 else 0.0,
                    "population_mass": tenure_mass,
                    "population_share_within_age": tenure_mass / max(mass[j], 1e-12),
                    "positive_liquid_assets_per_household": positive_assets
                    / max(tenure_mass, 1e-12),
                    "liquid_debt_per_household": debt / max(tenure_mass, 1e-12),
                    "gross_housing_wealth_per_household": gross_housing
                    / max(tenure_mass, 1e-12),
                    "net_wealth_per_household": (
                        positive_assets - debt + gross_housing
                    )
                    / max(tenure_mass, 1e-12),
                }
            )
    with (args.outdir / "wealth_by_age_tenure.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(tenure_rows[0]))
        writer.writeheader()
        writer.writerows(tenure_rows)
    timing_rows = [
        {"timing": timing, **values}
        for timing, values in balance_sheets.items()
    ]
    with (args.outdir / "wealth_timing_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timing_rows[0]))
        writer.writeheader()
        writer.writerows(timing_rows)
    timing_by_age = {}
    timing_mass = None
    for timing, timing_distribution in timing_distributions.items():
        timing_by_age[timing], timing_mass = net_wealth_per_household_by_age(
            timing_distribution,
            asset_grid,
            prices,
            owner_sizes,
        )
    timing_age_rows = [
        {
            "age": float(age),
            "population_mass": float(timing_mass[index]),
            **{
                timing: float(values[index])
                for timing, values in timing_by_age.items()
            },
        }
        for index, age in enumerate(ages)
    ]
    with (args.outdir / "wealth_timing_by_age.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timing_age_rows[0]))
        writer.writeheader()
        writer.writerows(timing_age_rows)
    target_timing_rows = [
        {"timing": timing, **values}
        for timing, values in timing_target_moments.items()
    ]
    with (args.outdir / "target_moment_timing_comparison.csv").open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(target_timing_rows[0]))
        writer.writeheader()
        writer.writerows(target_timing_rows)

    liquid_pc = liquid_age / np.maximum(mass, 1e-12)
    positive_liquid_pc = positive_liquid_age / np.maximum(mass, 1e-12)
    liquid_debt_pc = liquid_debt_age / np.maximum(mass, 1e-12)
    housing_pc = housing_age / np.maximum(mass, 1e-12)
    earnings_pc = earnings_age / np.maximum(mass, 1e-12)
    fig, axis = plt.subplots(figsize=(8.0, 4.8))
    axis.plot(ages, liquid_pc, marker="o", label="Liquid net worth")
    axis.plot(ages, housing_pc, marker="o", label="Gross housing wealth")
    axis.plot(ages, liquid_pc + housing_pc, color="black", linewidth=2.0, label="Total wealth")
    axis.plot(ages, earnings_pc, linestyle="--", label="Annual after-tax earnings")
    axis.axhline(0.0, color="0.75", linewidth=0.8)
    axis.set_xlabel("Age")
    axis.set_ylabel("Model units per household")
    axis.set_title("Wealth structure at the certified one-shot winner")
    axis.legend(frameon=False, ncol=2)
    axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(args.outdir / "wealth_by_age.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8.0, 4.8))
    axis.plot(ages, positive_liquid_pc, marker="o", label="Positive liquid assets")
    axis.plot(ages, housing_pc, marker="o", label="Gross housing wealth")
    axis.plot(ages, -liquid_debt_pc, marker="o", label="Liquid debt (liability)")
    axis.plot(
        ages,
        positive_liquid_pc + housing_pc - liquid_debt_pc,
        color="black",
        linewidth=2.0,
        label="Net wealth",
    )
    axis.axhline(0.0, color="0.75", linewidth=0.8)
    axis.set_xlabel("Age")
    axis.set_ylabel("Model units per household")
    axis.set_title("Household balance sheet at the certified one-shot winner")
    axis.legend(frameon=False, ncol=2)
    axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(args.outdir / "balance_sheet_by_age.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8.0, 4.8))
    axis.plot(
        ages,
        timing_by_age["beginning_of_period_consistent"],
        marker="o",
        label="Beginning of period (consistent)",
    )
    axis.plot(
        ages,
        timing_by_age["hybrid_current_target"],
        marker="o",
        label="Current target (hybrid)",
    )
    axis.plot(
        ages,
        timing_by_age["posttransaction_consistent"],
        marker="o",
        label="Post transaction (consistent)",
    )
    axis.axhline(0.0, color="0.75", linewidth=0.8)
    axis.set_xlabel("Age")
    axis.set_ylabel("Net wealth per household")
    axis.set_title("The active wealth target mixes balance-sheet timings")
    axis.legend(frameon=False)
    axis.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(args.outdir / "wealth_timing_by_age.png", dpi=180)
    plt.close(fig)
    beginning = balance_sheets["beginning_of_period_consistent"]
    posttransaction = balance_sheets["posttransaction_consistent"]
    hybrid = balance_sheets["hybrid_current_target"]
    report = f"""# Wealth-structure and timing diagnostic

This is a diagnostic at the certified corrected one-shot winner. It does not
change the calibration, targets, parameters, or model solution.

## Main finding

The active wealth moment combines inherited beginning-of-period liquid wealth
with the household's newly chosen tenure. That is not an internally consistent
balance sheet: it can double-count housing for buyers and omit sale proceeds
for sellers.

| Measurement | Wealth / annual after-tax earnings | Net wealth |
|---|---:|---:|
| Beginning of period, inherited tenure and assets | {beginning["wealth_to_annual_after_tax_earnings"]:.6f} | {beginning["net_wealth"]:.6f} |
| Active target: inherited assets, newly chosen tenure | {hybrid["wealth_to_annual_after_tax_earnings"]:.6f} | {hybrid["net_wealth"]:.6f} |
| Post transaction, newly chosen tenure and assets | {posttransaction["wealth_to_annual_after_tax_earnings"]:.6f} | {posttransaction["net_wealth"]:.6f} |
| Empirical target | {target_ratio:.6f} | {required_wealth:.6f} at model earnings |

Thus the reported ratio of {hybrid["wealth_to_annual_after_tax_earnings"]:.3f}
overstates both consistent alternatives. The low-wealth problem is not resolved
by correcting the timing; it becomes larger.

The same hybrid distribution also feeds the bequest-flow/wealth and old-estate
p90/p50 target rows. Their timing variants are reported in
`target_moment_timing_comparison.csv`. The correct economic timing for the
bequest flow requires a separate decision—death should be aligned with the
model's within-period sequence—so this diagnostic does not silently replace
that row.

The hybrid asset assignment predates the July 22 target revision. Therefore
older M4/M5 total-estate moments that use the same distribution also require
re-audit; this is not confined to the new fourteen-moment exercise.

At age {float(ages[-1]):.0f}, the inherited beginning-of-period owner share is
{final_age_owner_shares["beginning_of_period_consistent"]:.3f}, while the current
chosen owner share is
{final_age_owner_shares["posttransaction_consistent"]:.3f}. The hybrid attaches
the former balance-sheet positions to the latter tenure choices, creating the
extreme terminal distortion in the age plot.

## Active-target balance-sheet decomposition

- Positive liquid assets: {positive_liquid:.6f}
- Liquid debt: {liquid_debt:.6f}
- Gross housing wealth: {housing:.6f}
- Net wealth: {wealth:.6f}
- Annual after-tax earnings: {earnings:.6f}

See `wealth_timing_by_age.png` for the age pattern and
`wealth_timing_comparison.csv` for the exact aggregate accounting.
"""
    (args.outdir / "REPORT.md").write_text(report)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
