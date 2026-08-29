#!/usr/bin/env python3
"""Compare exact pre/post deterministic-tenure state diagnostics."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterator

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PRE = ROOT / "output/model/e5f_pf_person_policy_tenure_kink_20260827a_h128_l0250"
DEFAULT_POST = ROOT / "output/model/e5f_pf_person_policy_tenure_kink_20260827a_h128_l0275"
DEFAULT_OUTPUT = ROOT / "output/model/e5f_pf_person_policy_tenure_kink_state_map_20260827a"
EXPECTED_PRE_PATH_SHA256 = "9abc5ca2f1e4db721bbdf64afc5d3a12c2e835bf4c9df0a9288596b8481c1e95"
EXPECTED_POST_PATH_SHA256 = "f779a593e64a08dc0d31ccaea17bc7bfa098f29beb0c9739230e4a8be1992a5c"
EXPECTED_YEARS = (2351, 2363, 2367)
STATE_KEY_FIELDS = (
    "calendar_year",
    "wealth_index",
    "origin_tenure",
    "market",
    "age_index",
    "income_state",
    "parity",
    "child_state",
)


def parse_expected_years(value: str) -> tuple[int, ...]:
    try:
        years = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "expected years must be comma-separated integers"
        ) from exc
    if not years:
        raise argparse.ArgumentTypeError("at least one expected year is required")
    if any(year <= 0 for year in years):
        raise argparse.ArgumentTypeError("expected years must be positive")
    if tuple(sorted(set(years))) != years:
        raise argparse.ArgumentTypeError(
            "expected years must be unique and strictly increasing"
        )
    return years


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pre-dir", type=Path, default=DEFAULT_PRE)
    parser.add_argument("--post-dir", type=Path, default=DEFAULT_POST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--expected-pre-path-sha256", default=EXPECTED_PRE_PATH_SHA256)
    parser.add_argument("--expected-post-path-sha256", default=EXPECTED_POST_PATH_SHA256)
    parser.add_argument(
        "--expected-years",
        type=parse_expected_years,
        default=EXPECTED_YEARS,
        help="comma-separated diagnostic years (default: 2351,2363,2367)",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected an object in {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write an empty comparison: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def state_key(row: dict[str, str]) -> tuple[int, ...]:
    return tuple(int(row[name]) for name in STATE_KEY_FIELDS)


def iter_state_rows(path: Path) -> Iterator[tuple[tuple[int, ...], dict[str, str]]]:
    previous: tuple[int, ...] | None = None
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not set(STATE_KEY_FIELDS).issubset(reader.fieldnames):
            raise ValueError(f"State table lacks its exact key: {path}")
        for row in reader:
            key = state_key(row)
            if previous is not None and key <= previous:
                raise ValueError(f"State table is not strictly sorted: {path}")
            previous = key
            yield key, row


def merge_state_rows(
    pre_path: Path, post_path: Path
) -> Iterator[tuple[tuple[int, ...], dict[str, str] | None, dict[str, str] | None]]:
    left = iter(iter_state_rows(pre_path))
    right = iter(iter_state_rows(post_path))
    pre_item = next(left, None)
    post_item = next(right, None)
    while pre_item is not None or post_item is not None:
        if post_item is None or (pre_item is not None and pre_item[0] < post_item[0]):
            yield pre_item[0], pre_item[1], None
            pre_item = next(left, None)
        elif pre_item is None or post_item[0] < pre_item[0]:
            yield post_item[0], None, post_item[1]
            post_item = next(right, None)
        else:
            yield pre_item[0], pre_item[1], post_item[1]
            pre_item = next(left, None)
            post_item = next(right, None)


def tenure_count(row: dict[str, str]) -> int:
    indices = {
        int(match.group(1))
        for name in row
        if (match := re.fullmatch(r"tenure_(\d+)_value", name)) is not None
    }
    if not indices or indices != set(range(max(indices) + 1)):
        raise ValueError("Tenure option columns are not contiguous from zero")
    return len(indices)


def branch_value(row: dict[str, str], tenure: int, suffix: str) -> float:
    return float(row[f"tenure_{int(tenure)}_{suffix}"])


def load_path_rows(directory: Path) -> dict[int, dict[str, str]]:
    with (directory / "transition_path.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if [int(row["period"]) for row in rows] != list(range(len(rows))):
        raise ValueError(f"Noncontiguous transition path: {directory}")
    return {int(row["calendar_year"]): row for row in rows}


def maximum_path_feasibility_projection(
    path_rows: dict[int, dict[str, str]],
) -> float:
    values: list[float] = []
    for row in path_rows.values():
        if "feasibility_frontier_projection_mass" not in row:
            raise ValueError("Transition path lacks feasibility projection mass")
        value = float(row["feasibility_frontier_projection_mass"])
        if not math.isfinite(value):
            raise ValueError("Transition path has nonfinite feasibility projection mass")
        values.append(abs(value))
    if not values:
        raise ValueError("Transition path is empty")
    return max(values)


def audit_packet(
    directory: Path,
    *,
    expected_initial_path_sha256: str,
    expected_years: tuple[int, ...] = EXPECTED_YEARS,
) -> tuple[dict[str, Any], dict[str, Any], dict[int, Path], dict[int, dict[str, str]]]:
    summary_path = directory / "summary.json"
    diagnostic_path = directory / "tenure_diagnostic_summary.json"
    summary = load_json(summary_path)
    diagnostic = load_json(diagnostic_path)
    path_rows = load_path_rows(directory)
    if summary.get("status") != "complete_unpromoted_person_demography_policy_path":
        raise RuntimeError(f"Incomplete production summary: {directory}")
    if diagnostic.get("status") != "complete_unpromoted_tenure_value_diagnostic":
        raise RuntimeError(f"Incomplete tenure diagnostic: {directory}")
    if diagnostic.get("initial_path_sha256") != expected_initial_path_sha256:
        raise RuntimeError(f"Initial path hash differs: {directory}")
    if sha256(summary_path) != diagnostic.get("production_summary_sha256"):
        raise RuntimeError(f"Production summary hash differs: {directory}")
    if summary.get("history_reproduction_status") != "passed" or not all(
        (summary.get("accounting_gates") or {}).values()
    ):
        raise RuntimeError(f"Accounting or history gate failed: {directory}")
    if maximum_path_feasibility_projection(path_rows) != 0.0:
        raise RuntimeError(f"Feasibility projection was nonzero: {directory}")
    if tuple(diagnostic.get("diagnostic_years") or ()) != expected_years:
        raise RuntimeError(f"Diagnostic years differ: {directory}")
    tables: dict[int, Path] = {}
    for validation in diagnostic.get("validations") or []:
        year = int(validation["calendar_year"])
        if int(validation["choice_mismatch_count"]) != 0:
            raise RuntimeError(f"Choice reconstruction failed in {year}: {directory}")
        if float(validation["maximum_compiled_value_error"]) > 2.0e-11:
            raise RuntimeError(f"Live value reconstruction failed in {year}: {directory}")
        if float(validation["maximum_one_market_location_probability_gap"]) > 2.0e-12:
            raise RuntimeError(f"One-market location gate failed in {year}: {directory}")
        table = directory / Path(validation["state_csv"]).name
        if sha256(table) != validation["state_csv_sha256"]:
            raise RuntimeError(f"State-table hash differs in {year}: {directory}")
        tables[year] = table
    if tuple(sorted(tables)) != expected_years:
        raise RuntimeError(f"Missing state tables: {directory}")
    return summary, diagnostic, tables, path_rows


def comparison_row(
    key: tuple[int, ...], pre: dict[str, str], post: dict[str, str]
) -> dict[str, Any]:
    pre_choice = int(pre["stored_tenure"])
    post_choice = int(post["stored_tenure"])
    fields = dict(zip(STATE_KEY_FIELDS, key))
    row: dict[str, Any] = {
        **fields,
        "liquid_wealth_pre": float(pre["liquid_wealth"]),
        "liquid_wealth_post": float(post["liquid_wealth"]),
        "age_years": float(pre["age_years"]),
        "income_multiplier": float(pre["income_multiplier"]),
        "pre_mass": float(pre["decision_mass"]),
        "post_mass": float(post["decision_mass"]),
        "pre_mass_share": float(pre["decision_mass_share"]),
        "post_mass_share": float(post["decision_mass_share"]),
        "pre_choice": pre_choice,
        "post_choice": post_choice,
        "pre_best_minus_second": float(pre["best_minus_second"]),
        "post_best_minus_second": float(post["best_minus_second"]),
        "pre_value_post_minus_pre_choice": branch_value(pre, post_choice, "value")
        - branch_value(pre, pre_choice, "value"),
        "post_value_post_minus_pre_choice": branch_value(post, post_choice, "value")
        - branch_value(post, pre_choice, "value"),
        "pre_choice_housing_pre_environment": branch_value(
            pre, pre_choice, "housing_services"
        ),
        "post_choice_housing_pre_environment": branch_value(
            pre, post_choice, "housing_services"
        ),
        "pre_choice_housing_post_environment": branch_value(
            post, pre_choice, "housing_services"
        ),
        "post_choice_housing_post_environment": branch_value(
            post, post_choice, "housing_services"
        ),
        "pre_choice_tax_pre_environment": branch_value(
            pre, pre_choice, "property_tax_per_head"
        ),
        "post_choice_tax_pre_environment": branch_value(
            pre, post_choice, "property_tax_per_head"
        ),
        "pre_choice_tax_post_environment": branch_value(
            post, pre_choice, "property_tax_per_head"
        ),
        "post_choice_tax_post_environment": branch_value(
            post, post_choice, "property_tax_per_head"
        ),
        "pre_branch_liquid_wealth": branch_value(
            pre, pre_choice, "branch_liquid_wealth"
        ),
        "counterfactual_post_branch_liquid_wealth_pre_environment": branch_value(
            pre, post_choice, "branch_liquid_wealth"
        ),
        "counterfactual_pre_branch_liquid_wealth_post_environment": branch_value(
            post, pre_choice, "branch_liquid_wealth"
        ),
        "post_branch_liquid_wealth": branch_value(
            post, post_choice, "branch_liquid_wealth"
        ),
        "pre_next_liquid_wealth": branch_value(pre, pre_choice, "next_liquid_wealth"),
        "post_next_liquid_wealth": branch_value(post, post_choice, "next_liquid_wealth"),
    }
    row["pre_weighted_direct_owner_effect_pp"] = (
        row["pre_mass_share"] * (float(post_choice > 0) - float(pre_choice > 0)) * 100.0
    )
    row["post_weighted_direct_owner_effect_pp"] = (
        row["post_mass_share"] * (float(post_choice > 0) - float(pre_choice > 0)) * 100.0
    )
    row["pre_weighted_direct_housing_effect_per_head"] = row["pre_mass_share"] * (
        row["post_choice_housing_pre_environment"]
        - row["pre_choice_housing_pre_environment"]
    )
    row["post_weighted_direct_housing_effect_per_head"] = row["post_mass_share"] * (
        row["post_choice_housing_post_environment"]
        - row["pre_choice_housing_post_environment"]
    )
    row["pre_weighted_direct_tax_effect_per_head"] = row["pre_mass_share"] * (
        row["post_choice_tax_pre_environment"] - row["pre_choice_tax_pre_environment"]
    )
    row["post_weighted_direct_tax_effect_per_head"] = row["post_mass_share"] * (
        row["post_choice_tax_post_environment"] - row["pre_choice_tax_post_environment"]
    )
    return row


def compare_year(
    year: int,
    pre_table: Path,
    post_table: Path,
    pre_path_row: dict[str, str],
    post_path_row: dict[str, str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    flips: list[dict[str, Any]] = []
    total_pre = total_post = 0.0
    owner_pre = owner_post = 0.0
    missing_pre = missing_post = 0.0
    common_pre_mass = common_post_mass = 0.0
    common_rule_effect_pre = common_rule_effect_post = 0.0
    row_count_pre = row_count_post = 0
    observed_tenure_count: int | None = None
    for key, pre, post in merge_state_rows(pre_table, post_table):
        if pre is not None:
            row_count_pre += 1
            mass = float(pre["decision_mass"])
            total_pre += mass
            owner_pre += mass * float(int(pre["stored_tenure"]) > 0)
            count = tenure_count(pre)
            observed_tenure_count = count if observed_tenure_count is None else observed_tenure_count
            if observed_tenure_count != count:
                raise ValueError("Tenure count changes within a state table")
        if post is not None:
            row_count_post += 1
            mass = float(post["decision_mass"])
            total_post += mass
            owner_post += mass * float(int(post["stored_tenure"]) > 0)
            count = tenure_count(post)
            observed_tenure_count = count if observed_tenure_count is None else observed_tenure_count
            if observed_tenure_count != count:
                raise ValueError("Tenure count changes within a state table")
        if pre is None:
            missing_pre += float(post["decision_mass"])
            continue
        if post is None:
            missing_post += float(pre["decision_mass"])
            continue
        pre_mass = float(pre["decision_mass"])
        post_mass = float(post["decision_mass"])
        common_pre_mass += pre_mass
        common_post_mass += post_mass
        pre_choice = int(pre["stored_tenure"])
        post_choice = int(post["stored_tenure"])
        common_rule_effect_pre += pre_mass * (
            float(post_choice > 0) - float(pre_choice > 0)
        )
        common_rule_effect_post += post_mass * (
            float(post_choice > 0) - float(pre_choice > 0)
        )
        if pre_choice != post_choice:
            flips.append(comparison_row(key, pre, post))
    decision_owner_rate_pre = owner_pre / total_pre
    decision_owner_rate_post = owner_post / total_post
    path_owner_rate_pre = float(pre_path_row["owner_rate"])
    path_owner_rate_post = float(post_path_row["owner_rate"])
    if (
        abs(decision_owner_rate_pre - path_owner_rate_pre) > 2.0e-12
        or abs(decision_owner_rate_post - path_owner_rate_post) > 2.0e-12
    ):
        raise RuntimeError(f"Decision-state ownership does not reproduce the path in {year}")
    summary = {
        "calendar_year": int(year),
        "pre_state_rows": row_count_pre,
        "post_state_rows": row_count_post,
        "tenure_option_count": int(observed_tenure_count or 0),
        "pre_decision_mass": total_pre,
        "post_decision_mass": total_post,
        "mass_present_only_post": missing_pre,
        "mass_present_only_pre": missing_post,
        "common_pre_mass": common_pre_mass,
        "common_post_mass": common_post_mass,
        "choice_flip_state_count": len(flips),
        "choice_flip_pre_mass_share": sum(float(row["pre_mass_share"]) for row in flips),
        "choice_flip_post_mass_share": sum(float(row["post_mass_share"]) for row in flips),
        "decision_owner_rate_pre": decision_owner_rate_pre,
        "decision_owner_rate_post": decision_owner_rate_post,
        "owner_rate_effect_percentage_points": 100.0
        * (decision_owner_rate_post - decision_owner_rate_pre),
        "common_support_rule_effect_pre_weight_percentage_points": 100.0
        * common_rule_effect_pre
        / total_pre,
        "common_support_rule_effect_post_weight_percentage_points": 100.0
        * common_rule_effect_post
        / total_post,
        "flip_direct_owner_effect_pre_weight_percentage_points": sum(
            float(row["pre_weighted_direct_owner_effect_pp"]) for row in flips
        ),
        "flip_direct_owner_effect_post_weight_percentage_points": sum(
            float(row["post_weighted_direct_owner_effect_pp"]) for row in flips
        ),
        "flip_direct_housing_effect_pre_weight_per_head": sum(
            float(row["pre_weighted_direct_housing_effect_per_head"]) for row in flips
        ),
        "flip_direct_housing_effect_post_weight_per_head": sum(
            float(row["post_weighted_direct_housing_effect_per_head"]) for row in flips
        ),
        "flip_direct_tax_effect_pre_weight_per_head": sum(
            float(row["pre_weighted_direct_tax_effect_per_head"]) for row in flips
        ),
        "flip_direct_tax_effect_post_weight_per_head": sum(
            float(row["post_weighted_direct_tax_effect_per_head"]) for row in flips
        ),
        "path_market_residual_pre": float(pre_path_row["relative_market_residual"]),
        "path_market_residual_post": float(post_path_row["relative_market_residual"]),
        "path_fiscal_residual_pre": float(pre_path_row["government_budget_residual"]),
        "path_fiscal_residual_post": float(post_path_row["government_budget_residual"]),
    }
    return flips, summary


def make_plot(flips: list[dict[str, Any]], summaries: list[dict[str, Any]], path: Path) -> None:
    years = [int(row["calendar_year"]) for row in summaries]
    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axes[0, 0].bar(
        [str(year) for year in years],
        [100.0 * float(row["choice_flip_pre_mass_share"]) for row in summaries],
        label="Pre-switch mass",
    )
    axes[0, 0].bar(
        [str(year) for year in years],
        [100.0 * float(row["choice_flip_post_mass_share"]) for row in summaries],
        alpha=0.55,
        label="Post-switch mass",
    )
    axes[0, 0].set(title="Mass in states whose tenure choice flips", ylabel="Percent of heads")
    axes[0, 0].legend(frameon=False)

    for year in years:
        rows = [row for row in flips if int(row["calendar_year"]) == year]
        if rows:
            axes[0, 1].scatter(
                [float(row["pre_value_post_minus_pre_choice"]) for row in rows],
                [float(row["post_value_post_minus_pre_choice"]) for row in rows],
                s=np.maximum(
                    8.0,
                    1400.0 * np.asarray([float(row["pre_mass_share"]) for row in rows]),
                ),
                alpha=0.6,
                label=str(year),
            )
    axes[0, 1].axhline(0.0, color="0.6", linewidth=0.8)
    axes[0, 1].axvline(0.0, color="0.6", linewidth=0.8)
    axes[0, 1].set(
        title="Value ranking crosses zero at flipping states",
        xlabel="Post-choice minus pre-choice value, pre path",
        ylabel="Post-choice minus pre-choice value, post path",
    )
    axes[0, 1].legend(frameon=False)

    axes[1, 0].bar(
        [str(year) for year in years],
        [float(row["owner_rate_effect_percentage_points"]) for row in summaries],
        label="Actual ownership effect",
    )
    axes[1, 0].plot(
        [str(year) for year in years],
        [
            float(row["flip_direct_owner_effect_pre_weight_percentage_points"])
            for row in summaries
        ],
        marker="o",
        color="tab:red",
        label="Direct flip effect (pre weights)",
    )
    axes[1, 0].set(title="Ownership effect and direct tenure-rule component", ylabel="Percentage points")
    axes[1, 0].legend(frameon=False)

    by_age: dict[float, float] = defaultdict(float)
    for row in flips:
        if int(row["calendar_year"]) == years[0]:
            by_age[float(row["age_years"])] += float(row["pre_mass_share"])
    axes[1, 1].bar([str(int(age)) for age in sorted(by_age)], [100.0 * by_age[age] for age in sorted(by_age)])
    axes[1, 1].set(title=f"Flipping-state mass by age in {years[0]}", xlabel="Age", ylabel="Percent of heads")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
    figure.savefig(path, dpi=200)
    figure.savefig(path.with_suffix(".pdf"))
    plt.close(figure)


def main() -> None:
    args = parse_args()
    pre_dir = args.pre_dir.resolve()
    post_dir = args.post_dir.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    pre_summary, pre_diag, pre_tables, pre_path = audit_packet(
        pre_dir,
        expected_initial_path_sha256=args.expected_pre_path_sha256,
        expected_years=args.expected_years,
    )
    post_summary, post_diag, post_tables, post_path = audit_packet(
        post_dir,
        expected_initial_path_sha256=args.expected_post_path_sha256,
        expected_years=args.expected_years,
    )
    if pre_diag.get("source_hashes") != post_diag.get("source_hashes"):
        raise RuntimeError("Pre/post diagnostics used different source hashes")
    all_flips: list[dict[str, Any]] = []
    year_summaries: list[dict[str, Any]] = []
    for year in args.expected_years:
        flips, summary = compare_year(
            year,
            pre_tables[year],
            post_tables[year],
            pre_path[year],
            post_path[year],
        )
        all_flips.extend(flips)
        year_summaries.append(summary)
    if not all_flips:
        raise RuntimeError("No positive-mass tenure-choice flips were found")
    write_csv(output_dir / "flipping_states.csv", all_flips)
    write_csv(output_dir / "year_summary.csv", year_summaries)
    make_plot(all_flips, year_summaries, output_dir / "tenure_kink_state_map.png")
    payload = {
        "status": "complete_unaccepted_tenure_kink_state_map",
        "interpretation": (
            "Exact state-level comparison of the two audited H128 path evaluations "
            "bracketing the deterministic tenure branch switch. This diagnoses the "
            "local fixed-point failure and is not a solved 2% transition."
        ),
        "pre_directory": str(pre_dir),
        "post_directory": str(post_dir),
        "pre_initial_path_sha256": args.expected_pre_path_sha256,
        "post_initial_path_sha256": args.expected_post_path_sha256,
        "source_hashes": pre_diag["source_hashes"],
        "pre_production_summary_sha256": sha256(pre_dir / "summary.json"),
        "post_production_summary_sha256": sha256(post_dir / "summary.json"),
        "pre_diagnostic_summary_sha256": sha256(pre_dir / "tenure_diagnostic_summary.json"),
        "post_diagnostic_summary_sha256": sha256(post_dir / "tenure_diagnostic_summary.json"),
        "year_summaries": year_summaries,
        "flipping_state_count": len(all_flips),
        "outputs": {
            "flipping_states_csv": sha256(output_dir / "flipping_states.csv"),
            "year_summary_csv": sha256(output_dir / "year_summary.csv"),
            "plot_png": sha256(output_dir / "tenure_kink_state_map.png"),
            "plot_pdf": sha256(output_dir / "tenure_kink_state_map.pdf"),
        },
        "pre_path_status": pre_summary.get("path_status"),
        "post_path_status": post_summary.get("path_status"),
    }
    write_json(output_dir / "summary.json", payload)


if __name__ == "__main__":
    main()
