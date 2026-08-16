#!/usr/bin/env python3
"""Build the paper-facing report for the E5F dated-transition calibration pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from intergen_eqscale_seq_optimized.e5f_floor_profile import E5F_DOMAIN


TRANSITION_DOMAIN = tuple(row for row in E5F_DOMAIN if row[0] != "psi_child") + (
    ("psi_child_2023", -1.25, 0.20, "asinh"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--round-dir", type=Path, action="append", required=True)
    parser.add_argument("--repeat-dir", type=Path, action="append", default=[])
    parser.add_argument("--policy-baseline-dir", type=Path, default=None)
    parser.add_argument("--policy-treatment-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def candidate_records(round_dirs: list[Path]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for round_index, directory in enumerate(round_dirs, start=1):
        for summary_path in sorted(directory.glob("task_*/summary.json")):
            summary = read_json(summary_path)
            candidate = dict(summary["best_candidate"])
            panel = dict(summary["panel_design"])
            candidate.update(
                {
                    "round": round_index,
                    "round_name": directory.name,
                    "task_dir": str(summary_path.parent),
                    "source_sha256": summary["source_sha256"],
                    "target_set": summary["target_set"],
                    "target_fingerprint": summary["target_fingerprint"],
                    "old_psi_child": summary["old_psi_child"],
                    "old_completed_fertility": summary["old_model_completed_fertility"],
                    "theta_json": json.dumps(candidate["theta"], sort_keys=True),
                    "task_id": int(panel["task_id"]),
                    "design": str(panel["design"]),
                }
            )
            records.append(candidate)
    return records


def fit_lookup(task_dir: Path) -> dict[str, dict[str, float]]:
    rows = read_csv(task_dir / "target_fit_long.csv")
    return {
        row["moment"]: {
            "target": float(row["target"]),
            "model": float(row["model"]),
            "gap": float(row["gap"]),
            "weight": float(row["weight"]),
            "loss_contribution": float(row["loss_contribution"]),
            "standardized_gap": float(row["standardized_gap"]),
        }
        for row in rows
    }


def parameter_table(best: dict[str, Any]) -> list[dict[str, Any]]:
    theta = dict(best["theta"])
    bounds = {name: (lower, upper) for name, lower, upper, _ in TRANSITION_DOMAIN}
    rows: list[dict[str, Any]] = []
    for name, lower, upper, _ in TRANSITION_DOMAIN:
        if name == "psi_child_2023":
            value = float(best["new_psi_child"])
        elif name == "beta_annual":
            value = float(theta["beta"]) ** 0.25
        else:
            value = float(theta[name])
        span = upper - lower
        rows.append(
            {
                "parameter": name,
                "estimate": value,
                "lower_bound": lower,
                "upper_bound": upper,
                "near_bound": min(value - lower, upper - value) <= 0.02 * span,
                "status": "free_in_bounded_transition_search",
            }
        )
    rows.extend(
        [
            {
                "parameter": "psi_child_2007",
                "estimate": float(best["old_psi_child"]),
                "lower_bound": -3.0,
                "upper_bound": 3.0,
                "near_bound": False,
                "status": "derived_each_candidate_to_match_old_completed_fertility_2.12",
            },
            {
                "parameter": "outside_origin_entry_share",
                "estimate": 0.169,
                "lower_bound": math.nan,
                "upper_bound": math.nan,
                "near_bound": False,
                "status": "externally_fixed_provisional_anchor",
            },
        ]
    )
    return rows


def turning_point_rows(path_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    series = {
        "household_index": "population_index",
        "housing_cost_index": "asset_price_index",
        "births_per_adult_household": "topcode_adjusted_births_per_adult",
    }
    rows: list[dict[str, Any]] = []
    for label, field in series.items():
        values = np.asarray([float(row[field]) for row in path_rows])
        for index in range(1, len(path_rows) - 1):
            turning_type = None
            if values[index] > values[index - 1] and values[index] > values[index + 1]:
                turning_type = "local_peak"
            elif values[index] < values[index - 1] and values[index] < values[index + 1]:
                turning_type = "local_trough"
            if turning_type is not None:
                rows.append(
                    {
                        "series": label,
                        "turning_point": turning_type,
                        "period": int(path_rows[index]["period"]),
                        "year": 2007
                        + int(float(path_rows[index]["years_from_start"])),
                        "value": float(values[index]),
                    }
                )
    return rows


def make_figures(
    output_dir: Path,
    comparison: list[dict[str, Any]],
    transition_rows: list[dict[str, str]],
    round_rows: list[dict[str, Any]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    short_labels = {
        "tfr": "Completed fertility",
        "childless_rate": "Childlessness",
        "mean_age_first_birth": "Mean age at first birth",
        "share_first_births_age30plus": "First births at 30+",
        "housing_increment_0to1": "Rooms after first birth",
        "prime30_55_parent_3plus_minus_1to2_mean_rooms": "Room gap: 3+ vs. 1–2 children",
        "own_family_gap": "Family ownership gap",
        "own_rate": "Ownership rate",
        "aggregate_mean_occupied_rooms_18_85": "Mean occupied rooms",
        "aggregate_wealth_to_annual_gross_labor_earnings": "Wealth / earnings",
        "annual_bequest_flow_to_aggregate_wealth": "Bequests / wealth",
        "old_total_wealth_to_annual_income_p90_p50_7684": "Old wealth p90/p50",
    }
    moments = [row["moment"] for row in comparison]
    positions = np.arange(len(moments))
    fig, ax = plt.subplots(figsize=(8.6, 6.2), constrained_layout=True)
    ax.scatter(
        [row["anchor_standardized_gap"] for row in comparison],
        positions - 0.13,
        label="Transition anchor",
    )
    ax.scatter(
        [row["best_standardized_gap"] for row in comparison],
        positions + 0.13,
        label="Best bounded pilot",
    )
    ax.axvline(0.0, color="black", lw=0.8)
    ax.set_yticks(positions)
    ax.set_yticklabels([short_labels.get(moment, moment) for moment in moments])
    ax.invert_yaxis()
    ax.set_xlabel("Signed square root of loss contribution")
    ax.set_title("Dated 2023 target fit")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(output_dir / "target_fit_comparison.png", dpi=220)
    fig.savefig(output_dir / "target_fit_comparison.pdf")
    plt.close(fig)

    years = np.asarray([float(row["years_from_start"]) for row in transition_rows])
    population = np.asarray([float(row["population_index"]) for row in transition_rows])
    housing = np.asarray([float(row["asset_price_index"]) for row in transition_rows])
    psi = np.asarray([float(row["psi_child"]) for row in transition_rows])
    births = np.asarray(
        [float(row["birth_children_topcode_adjusted"]) for row in transition_rows]
    )
    births /= births[0]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    axes[0].plot(years, population, marker="o", label="Households")
    axes[0].plot(years, housing, marker="o", label="Housing cost")
    axes[0].axhline(1.0, color="black", lw=0.7)
    axes[0].set(xlabel="Years from 2007", ylabel="Index, 2007 = 1")
    axes[0].set_title("Households and housing costs both rise")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)
    birth_line = axes[1].plot(years, births, marker="o", label="Birth flow index")
    preference_axis = axes[1].twinx()
    preference_line = preference_axis.plot(
        years, psi, marker="o", color="tab:orange", label=r"Preference shifter $\psi_t$"
    )
    axes[1].set(xlabel="Years from 2007", ylabel="Birth flow, 2007 = 1")
    preference_axis.set_ylabel(r"Preference shifter $\psi_t$")
    axes[1].set_title("Preference decline lowers births")
    lines = birth_line + preference_line
    axes[1].legend(lines, [line.get_label() for line in lines], frameon=False)
    axes[1].grid(alpha=0.25)
    fig.savefig(output_dir / "best_transition_path.png", dpi=220)
    fig.savefig(output_dir / "best_transition_path.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    ax.plot(
        [row["round"] for row in round_rows],
        [row["best_loss"] for row in round_rows],
        marker="o",
    )
    ax.set(xlabel="Search round", ylabel="Best transition loss", xticks=[row["round"] for row in round_rows])
    ax.grid(alpha=0.25)
    fig.savefig(output_dir / "search_progress.png", dpi=220)
    fig.savefig(output_dir / "search_progress.pdf")
    plt.close(fig)


def make_policy_figure(output_dir: Path, rows: list[dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    years = np.asarray([float(row["year"]) for row in rows])
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0), constrained_layout=True)
    axes[0].plot(
        years,
        [row["population_percent_change"] for row in rows],
        marker="o",
        label="Households",
    )
    axes[0].plot(
        years,
        [row["housing_cost_percent_change"] for row in rows],
        marker="o",
        label="Housing cost",
    )
    axes[0].axhline(0.0, color="black", lw=0.7)
    axes[0].set(xlabel="Calendar year", ylabel="Percent difference from baseline")
    axes[0].set_title("Aggregate effects")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].plot(
        years,
        [row["owner_rate_pp_change"] for row in rows],
        marker="o",
        label="All households",
    )
    axes[1].plot(
        years,
        [row["dependent_child_owner_rate_pp_change"] for row in rows],
        marker="o",
        label="Dependent-child households",
    )
    axes[1].axhline(0.0, color="black", lw=0.7)
    axes[1].set(xlabel="Calendar year", ylabel="Ownership change (percentage points)")
    axes[1].set_title("Tenure effects")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)
    axes[2].plot(
        years,
        [row["birth_rate_percent_change"] for row in rows],
        marker="o",
        color="tab:green",
    )
    axes[2].axhline(0.0, color="black", lw=0.7)
    axes[2].set(xlabel="Calendar year", ylabel="Birth-rate difference (percent)")
    axes[2].set_title("Fertility effect")
    axes[2].grid(alpha=0.25)
    fig.savefig(output_dir / "post2023_ltv95_policy.png", dpi=220)
    fig.savefig(output_dir / "post2023_ltv95_policy.pdf")
    plt.close(fig)


def make_long_run_figure(
    output_dir: Path, baseline_rows: list[dict[str, str]]
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    years = np.asarray(
        [2007 + float(row["years_from_start"]) for row in baseline_rows]
    )
    population = np.asarray(
        [float(row["population_index"]) for row in baseline_rows]
    )
    housing = np.asarray(
        [float(row["asset_price_index"]) for row in baseline_rows]
    )
    birth_rate = np.asarray(
        [float(row["topcode_adjusted_births_per_adult"]) for row in baseline_rows]
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2), constrained_layout=True)
    axes[0].plot(years, population, label="Households")
    axes[0].plot(years, housing, label="Housing cost")
    axes[0].axvline(2023, color="black", lw=0.8, ls="--")
    axes[0].axhline(1.0, color="black", lw=0.6)
    axes[0].set(xlabel="Calendar year", ylabel="Index, 2007 = 1")
    axes[0].set_title("Demographic momentum and delayed adjustment")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].plot(years, birth_rate, color="tab:green")
    axes[1].axvline(2023, color="black", lw=0.8, ls="--")
    axes[1].set(xlabel="Calendar year", ylabel="Births per adult household")
    axes[1].set_title("Birth rate recovers as housing costs fall")
    axes[1].grid(alpha=0.25)
    fig.savefig(output_dir / "long_run_no_policy_transition.png", dpi=220)
    fig.savefig(output_dir / "long_run_no_policy_transition.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    round_dirs = [path.resolve() for path in args.round_dir]
    repeat_dirs = [path.resolve() for path in args.repeat_dir]
    if (args.policy_baseline_dir is None) != (args.policy_treatment_dir is None):
        raise ValueError("Both policy path directories must be supplied together")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    records = candidate_records(round_dirs)
    if not records:
        raise RuntimeError("No completed transition candidates were found")
    contracts = {
        (row["source_sha256"], row["target_set"], row["target_fingerprint"])
        for row in records
    }
    if len(contracts) != 1:
        raise RuntimeError(f"Mixed source or target contracts: {contracts}")
    best = min(records, key=lambda row: float(row["transition_loss"]))
    anchor_candidates = [
        row for row in records if row["round"] == 1 and int(row["task_id"]) == 1
    ]
    if len(anchor_candidates) != 1:
        raise RuntimeError("The round-one anchor is missing or duplicated")
    anchor = anchor_candidates[0]
    anchor_fit = fit_lookup(Path(anchor["task_dir"]))
    best_fit = fit_lookup(Path(best["task_dir"]))
    if set(anchor_fit) != set(best_fit):
        raise RuntimeError("Anchor and best target rows differ")
    comparison = []
    for moment in anchor_fit:
        a, b = anchor_fit[moment], best_fit[moment]
        comparison.append(
            {
                "moment": moment,
                "target": b["target"],
                "anchor_model": a["model"],
                "anchor_gap": a["gap"],
                "anchor_loss_contribution": a["loss_contribution"],
                "anchor_standardized_gap": a["standardized_gap"],
                "best_model": b["model"],
                "best_gap": b["gap"],
                "weight": b["weight"],
                "best_loss_contribution": b["loss_contribution"],
                "best_standardized_gap": b["standardized_gap"],
            }
        )
    transition_path = Path(best["task_dir"]) / "cases" / best["candidate"] / "transition_path.csv"
    transition_rows = read_csv(transition_path)
    round_rows = []
    for index, directory in enumerate(round_dirs, start=1):
        subset = [row for row in records if row["round"] == index]
        round_rows.append(
            {
                "round": index,
                "round_name": directory.name,
                "completed_candidates": len(subset),
                "best_loss": min(float(row["transition_loss"]) for row in subset),
            }
        )

    repeat_rows = []
    for directory in repeat_dirs:
        path = directory / "task_001" / "summary.json"
        summary = read_json(path)
        repeat = dict(summary["best_candidate"])
        repeat_fit = fit_lookup(path.parent)
        max_moment_gap = max(
            abs(repeat_fit[name]["model"] - best_fit[name]["model"])
            for name in best_fit
        )
        repeat_rows.append(
            {
                "repeat": directory.name,
                "loss": float(repeat["transition_loss"]),
                "loss_gap_from_selected": float(repeat["transition_loss"])
                - float(best["transition_loss"]),
                "max_market_residual": float(repeat["max_market_residual"]),
                "max_mass_residual": float(repeat["max_mass_residual"]),
                "max_abs_target_moment_gap_from_selected": max_moment_gap,
                "theta_matches_selected": repeat["theta"] == best["theta"],
            }
        )
    repeat_max_loss_gap = (
        max(abs(row["loss_gap_from_selected"]) for row in repeat_rows)
        if repeat_rows
        else math.nan
    )
    repeat_max_moment_gap = (
        max(abs(row["max_abs_target_moment_gap_from_selected"]) for row in repeat_rows)
        if repeat_rows
        else math.nan
    )

    policy_result = None
    if args.policy_baseline_dir is not None:
        baseline_dir = args.policy_baseline_dir.resolve() / "task_001"
        treatment_dir = args.policy_treatment_dir.resolve() / "task_001"
        baseline_summary = read_json(baseline_dir / "summary.json")
        treatment_summary = read_json(treatment_dir / "summary.json")
        if baseline_summary["policy_case"] != "none":
            raise RuntimeError("The policy baseline is not the no-policy path")
        if treatment_summary["policy_case"] != "dependent-child-ltv95":
            raise RuntimeError("The policy treatment is not dependent-child LTV95")
        for name in ("source_sha256", "target_set", "target_fingerprint"):
            if baseline_summary[name] != treatment_summary[name]:
                raise RuntimeError(f"Policy paths differ in {name}")
        policy_contract = (
            baseline_summary["source_sha256"],
            baseline_summary["target_set"],
            baseline_summary["target_fingerprint"],
        )
        if policy_contract not in contracts:
            raise RuntimeError("Policy paths do not match the calibration contract")
        if baseline_summary["post_2023_periods"] != treatment_summary[
            "post_2023_periods"
        ]:
            raise RuntimeError("Policy paths use different post-2023 horizons")
        baseline_candidate = baseline_summary["best_candidate"]
        treatment_candidate = treatment_summary["best_candidate"]
        if baseline_candidate["theta"] != treatment_candidate["theta"]:
            raise RuntimeError("Policy paths use different structural parameters")
        if not math.isclose(
            float(baseline_candidate["new_psi_child"]),
            float(treatment_candidate["new_psi_child"]),
            rel_tol=0.0,
            abs_tol=1e-14,
        ):
            raise RuntimeError("Policy paths use different fertility preferences")
        baseline_path = read_csv(
            baseline_dir / "cases" / "task_001" / "transition_path.csv"
        )
        treatment_path = read_csv(
            treatment_dir / "cases" / "task_001" / "transition_path.csv"
        )
        if len(baseline_path) != len(treatment_path):
            raise RuntimeError("Policy paths have different horizons")
        write_csv(output_dir / "post2023_no_policy_path.csv", baseline_path)
        write_csv(output_dir / "post2023_ltv95_path.csv", treatment_path)
        make_long_run_figure(output_dir, baseline_path)
        long_run_turning_points = turning_point_rows(baseline_path)
        if long_run_turning_points:
            write_csv(
                output_dir / "long_run_turning_points.csv",
                long_run_turning_points,
            )
        policy_rows: list[dict[str, Any]] = []
        for baseline_row, treatment_row in zip(
            baseline_path, treatment_path, strict=True
        ):
            period = int(baseline_row["period"])
            if period != int(treatment_row["period"]):
                raise RuntimeError("Policy paths have misaligned dates")
            if period < 4:
                continue
            population_base = float(baseline_row["adult_population"])
            housing_base = float(baseline_row["asset_price"])
            births_base = float(baseline_row["birth_children_topcode_adjusted"])
            birth_rate_base = float(
                baseline_row["topcode_adjusted_births_per_adult"]
            )
            policy_rows.append(
                {
                    "period": period,
                    "year": 2007 + int(float(baseline_row["years_from_start"])),
                    "population_percent_change": 100.0
                    * (float(treatment_row["adult_population"]) / population_base - 1.0),
                    "housing_cost_percent_change": 100.0
                    * (float(treatment_row["asset_price"]) / housing_base - 1.0),
                    "birth_flow_percent_change": 100.0
                    * (
                        float(treatment_row["birth_children_topcode_adjusted"])
                        / births_base
                        - 1.0
                    ),
                    "birth_rate_percent_change": 100.0
                    * (
                        float(treatment_row["topcode_adjusted_births_per_adult"])
                        / birth_rate_base
                        - 1.0
                    ),
                    "owner_rate_pp_change": 100.0
                    * (float(treatment_row["owner_rate"]) - float(baseline_row["owner_rate"])),
                    "dependent_child_owner_rate_pp_change": 100.0
                    * (
                        float(treatment_row["dependent_child_owner_rate"])
                        - float(baseline_row["dependent_child_owner_rate"])
                    ),
                }
            )
        write_csv(output_dir / "post2023_ltv95_policy_path.csv", policy_rows)
        make_policy_figure(output_dir, policy_rows)
        baseline_post2023 = [row for row in baseline_path if int(row["period"]) >= 4]
        population_peak = max(
            baseline_post2023, key=lambda row: float(row["population_index"])
        )
        housing_peak = max(
            baseline_post2023, key=lambda row: float(row["asset_price_index"])
        )
        birth_rate_trough = min(
            baseline_post2023,
            key=lambda row: float(row["topcode_adjusted_births_per_adult"]),
        )
        policy_result = {
            "status": "conditional_policy_diagnostic_not_welfare_counterfactual",
            "policy": "dependent-child-households LTV95",
            "policy_start_year": 2023,
            "numerical_gates": {
                "baseline_max_market_residual": float(
                    baseline_candidate["max_market_residual"]
                ),
                "treatment_max_market_residual": float(
                    treatment_candidate["max_market_residual"]
                ),
                "baseline_max_mass_residual": float(
                    baseline_candidate["max_mass_residual"]
                ),
                "treatment_max_mass_residual": float(
                    treatment_candidate["max_mass_residual"]
                ),
            },
            "no_policy_continuation": {
                "population_peak_year": 2007
                + int(float(population_peak["years_from_start"])),
                "population_peak_index": float(population_peak["population_index"]),
                "housing_cost_peak_year": 2007
                + int(float(housing_peak["years_from_start"])),
                "housing_cost_peak_index": float(housing_peak["asset_price_index"]),
                "birth_rate_trough_year": 2007
                + int(float(birth_rate_trough["years_from_start"])),
                "birth_rate_trough": float(
                    birth_rate_trough["topcode_adjusted_births_per_adult"]
                ),
                "birth_rate_2023": float(
                    baseline_post2023[0]["topcode_adjusted_births_per_adult"]
                ),
                "birth_rate_endpoint": float(
                    baseline_post2023[-1]["topcode_adjusted_births_per_adult"]
                ),
                "endpoint_year": 2007
                + int(float(baseline_post2023[-1]["years_from_start"])),
                "endpoint_population_index": float(
                    baseline_post2023[-1]["population_index"]
                ),
                "endpoint_housing_cost_index": float(
                    baseline_post2023[-1]["asset_price_index"]
                ),
                "endpoint_population_change_percent_last_period": 100.0
                * (
                    float(baseline_post2023[-1]["population_index"])
                    / float(baseline_post2023[-2]["population_index"])
                    - 1.0
                ),
                "endpoint_housing_cost_change_percent_last_period": 100.0
                * (
                    float(baseline_post2023[-1]["asset_price_index"])
                    / float(baseline_post2023[-2]["asset_price_index"])
                    - 1.0
                ),
                "endpoint_birth_rate_change_percent_last_period": 100.0
                * (
                    float(
                        baseline_post2023[-1][
                            "topcode_adjusted_births_per_adult"
                        ]
                    )
                    / float(
                        baseline_post2023[-2][
                            "topcode_adjusted_births_per_adult"
                        ]
                    )
                    - 1.0
                ),
                "turning_points": long_run_turning_points,
            },
            "endpoint": policy_rows[-1],
            "prepolicy_path_identity_max_abs_gap": max(
                abs(
                    float(baseline_path[period][field])
                    - float(treatment_path[period][field])
                )
                for period in range(4)
                for field in (
                    "asset_price",
                    "adult_population",
                    "birth_children_topcode_adjusted",
                    "owner_rate",
                )
            ),
        }

    candidate_table = []
    for row in sorted(records, key=lambda item: float(item["transition_loss"])):
        candidate_table.append(
            {
                "round": row["round"],
                "round_name": row["round_name"],
                "task_id": row["task_id"],
                "design": row["design"],
                "transition_loss": row["transition_loss"],
                "new_psi_child": row["new_psi_child"],
                "old_psi_child": row["old_psi_child"],
                "terminal_tfr": row["terminal_tfr"],
                "terminal_childless_rate": row["terminal_childless_rate"],
                "terminal_mean_age_first_birth": row["terminal_mean_age_first_birth"],
                "terminal_share_first_births_age30plus": row[
                    "terminal_share_first_births_age30plus"
                ],
                "terminal_housing_cost_index": row["terminal_housing_cost_index"],
                "theta_json": row["theta_json"],
            }
        )
    write_csv(output_dir / "all_completed_candidates.csv", candidate_table)
    write_csv(output_dir / "search_rounds.csv", round_rows)
    write_csv(output_dir / "target_fit_anchor_vs_best.csv", comparison)
    write_csv(output_dir / "best_parameter_table.csv", parameter_table(best))
    write_csv(output_dir / "best_transition_path.csv", transition_rows)
    if repeat_rows:
        write_csv(output_dir / "repeat_checks.csv", repeat_rows)
    make_figures(output_dir, comparison, transition_rows, round_rows)
    best_loss = float(best["transition_loss"])
    anchor_loss = float(anchor["transition_loss"])
    result = {
        "status": "complete_bounded_transition_calibration_pilot",
        "claim_status": "diagnostic_not_promoted_estimate",
        "contract": list(next(iter(contracts))),
        "rounds": round_rows,
        "completed_candidate_count": len(records),
        "anchor_loss": anchor_loss,
        "best_loss": best_loss,
        "loss_improvement": anchor_loss - best_loss,
        "loss_improvement_percent": 100.0 * (anchor_loss - best_loss) / anchor_loss,
        "best_candidate": best,
        "repeat_checks": repeat_rows,
        "repeat_max_abs_loss_gap": repeat_max_loss_gap,
        "repeat_max_abs_target_moment_gap": repeat_max_moment_gap,
        "conditional_policy": policy_result,
        "population_bridge": (
            "Census HH-3 household totals and national ACS household-head age shares; "
            "externally normalized, not generated by post-2007 fertility"
        ),
        "largest_best_loss_contributions": sorted(
            [
                {"moment": row["moment"], "loss_contribution": row["best_loss_contribution"]}
                for row in comparison
            ],
            key=lambda row: float(row["loss_contribution"]),
            reverse=True,
        ),
    }
    write_json(output_dir / "summary.json", result)
    policy_readme_lines: list[str] = []
    if policy_result is not None:
        no_policy = policy_result["no_policy_continuation"]
        endpoint = policy_result["endpoint"]
        policy_readme_lines = [
            "",
            "## Post-2023 continuation and tenure policy",
            "",
            "Holding the 2023 fertility preference fixed, the no-policy path peaks "
            f"at a household index of {no_policy['population_peak_index']:.3f} in "
            f"{no_policy['population_peak_year']} and a housing-cost index of "
            f"{no_policy['housing_cost_peak_index']:.3f} in "
            f"{no_policy['housing_cost_peak_year']}, before inherited demographic "
            "momentum dissipates.",
            "",
            f"The birth rate reaches its trough in {no_policy['birth_rate_trough_year']} "
            f"and is {no_policy['birth_rate_endpoint']:.4f} at the endpoint, compared "
            f"with {no_policy['birth_rate_2023']:.4f} in 2023. Over the last four-year "
            "interval, households and housing costs change by "
            f"{no_policy['endpoint_population_change_percent_last_period']:.3f} and "
            f"{no_policy['endpoint_housing_cost_change_percent_last_period']:.3f} "
            "percent, respectively.",
            "",
            "The dependent-child LTV95 comparison changes the endpoint household "
            f"index by {endpoint['population_percent_change']:.2f} percent, housing "
            f"costs by {endpoint['housing_cost_percent_change']:.2f} percent, the birth "
            f"rate by {endpoint['birth_rate_percent_change']:.2f} percent, and "
            "dependent-child ownership by "
            f"{endpoint['dependent_child_owner_rate_pp_change']:.2f} percentage points.",
            "This is a conditional transition diagnostic using the existing tenure",
            "mechanism, not a funded welfare counterfactual.",
        ]
    (output_dir / "README.md").write_text(
        "\n".join(
            [
                "# Dated transition-calibration pilot",
                "",
                "This packet asks whether the current sequential household model can be",
                "calibrated from an old steady-state benchmark followed by a fertility-",
                "preference trend to 2023. Each candidate derives the old preference",
                "intercept so completed fertility is 2.12, then reweights age masses to",
                "the observed 2007 householder profile. It matches the subsequent",
                "household-count and age paths through a declared Census/ACS bridge and",
                "measures all twelve calibration targets on the simulated 2023 distribution.",
                "The 2007 trend origin is a reduced-form timing normalization, not an",
                "estimate of the historical date of the U.S. preference change.",
                "",
                f"The bounded search improves loss from {anchor_loss:.3f} to {best_loss:.3f} "
                f"({100.0 * (anchor_loss - best_loss) / anchor_loss:.1f} percent). This is",
                "a feasibility and local-search result, not a promoted estimate or a",
                "claim of global optimality. The complete target and parameter tables",
                "show where the model still fails, especially childlessness and old-age",
                "wealth dispersion.",
                "",
                "Primary files: `summary.json`, `target_fit_anchor_vs_best.csv`,",
                "`best_parameter_table.csv`, `best_transition_path.csv`, and the",
                "diagnostic figures. `all_completed_candidates.csv` preserves the search",
                "trace and `repeat_checks.csv` records exact repeat verification.",
                "When supplied, `post2023_ltv95_policy_path.csv` and its figure compare",
                "the existing dependent-child LTV95 experiment with the no-policy path",
                "from the same matched 2023 state. `post2023_no_policy_path.csv` and",
                "`long_run_no_policy_transition.png` show the absolute continuation;",
                "`long_run_turning_points.csv` records the cohort-cycle extrema;",
                "this is not a welfare calculation.",
                *policy_readme_lines,
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(
        "TRANSITION_CALIBRATION_REPORT_COMPLETE "
        f"candidates={len(records)} anchor={anchor_loss:.6f} best={best_loss:.6f} "
        f"repeat_gap={repeat_max_loss_gap:.3e}",
        flush=True,
    )


if __name__ == "__main__":
    main()
