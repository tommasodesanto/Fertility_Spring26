#!/usr/bin/env python3
"""Collect paired H128 paths and exact state-level fertility decompositions."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
from collections import defaultdict
from itertools import zip_longest
from pathlib import Path
from typing import Any, Iterator

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MARKET_GATE = 2.0e-4
FISCAL_GATE = 2.5e-5
PERIOD_YEARS = 4.0
STATE_KEY = (
    "calendar_year",
    "wealth_index",
    "origin_tenure",
    "market",
    "age_index",
    "income_state",
    "parity",
    "child_state",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--reform-dir", type=Path, required=True)
    parser.add_argument("--baseline-diagnostic-dir", type=Path, required=True)
    parser.add_argument("--reform-diagnostic-dir", type=Path, required=True)
    parser.add_argument("--calibration-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected a JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write an empty table: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def preserve_calibration_tables(directory: Path, output: Path) -> dict[str, Any]:
    """Copy the complete active target and parameter tables with basic audits."""

    summary = read_json(directory / "summary.json")
    target_rows = read_csv(directory / "target_fit_long.csv")
    parameter_rows = read_csv(directory / "parameter_table.csv")
    if int(summary["target_count"]) != len(target_rows):
        raise RuntimeError("Calibration target count does not match target_fit_long.csv")
    free_rows = [row for row in parameter_rows if row["is_free_parameter"] == "True"]
    if int(summary["transition_free_parameter_count"]) != len(free_rows):
        raise RuntimeError("Calibration free-parameter count does not match parameter_table.csv")
    loss = sum(float(row["loss_contribution"]) for row in target_rows)
    reported = float(summary["best_candidate"]["transition_loss"])
    if not math.isclose(loss, reported, rel_tol=0.0, abs_tol=2.0e-10):
        raise RuntimeError(
            "Target loss contributions do not add to the reported loss: "
            f"{loss} versus {reported}"
        )
    shutil.copy2(
        directory / "target_fit_long.csv", output / "active_calibration_target_fit.csv"
    )
    shutil.copy2(
        directory / "parameter_table.csv", output / "active_calibration_parameters.csv"
    )
    return {
        "target_set": summary["target_set"],
        "target_fingerprint": summary["target_fingerprint"],
        "identification_status": summary["identification_status"],
        "transition_loss": reported,
        "target_count": len(target_rows),
        "free_parameter_count": len(free_rows),
        "target_fit_sha256": sha256(directory / "target_fit_long.csv"),
        "parameter_table_sha256": sha256(directory / "parameter_table.csv"),
    }


def validate_path(directory: Path, expected_case: str) -> tuple[dict[str, Any], list[dict[str, str]]]:
    summary_path = directory / "summary.json"
    path = directory / "transition_path.csv"
    summary = read_json(summary_path)
    rows = read_csv(path)
    if summary.get("status") != "complete_unpromoted_person_demography_policy_path":
        raise RuntimeError(f"Incomplete path packet: {directory}")
    if summary.get("case") != expected_case or not bool(summary.get("path_converged")):
        raise RuntimeError(f"Wrong or unconverged path packet: {directory}")
    if float(summary["maximum_market_residual"]) > MARKET_GATE:
        raise RuntimeError(f"Market gate failed: {directory}")
    if float(summary["maximum_fiscal_residual"]) > FISCAL_GATE:
        raise RuntimeError(f"Fiscal gate failed: {directory}")
    if not all((summary.get("accounting_gates") or {}).values()):
        raise RuntimeError(f"Accounting gate failed: {directory}")
    if not bool((summary.get("accounting_gates") or {}).get("distribution_feasibility")):
        raise RuntimeError(f"Distribution feasibility gate failed: {directory}")
    if summary.get("history_reproduction_status") != "passed":
        raise RuntimeError(f"History gate failed: {directory}")
    if not bool((summary.get("terminal_convergence") or {}).get("all_checks_pass")):
        raise RuntimeError(f"H128 terminal gate failed: {directory}")
    if len(rows) != int(summary["horizon"]):
        raise RuntimeError(f"Path length differs from declared horizon: {directory}")
    if [int(row["period"]) for row in rows] != list(range(len(rows))):
        raise RuntimeError(f"Noncontiguous path: {directory}")
    return summary, rows


def diagnostic_tables(directory: Path, path_sha256: str) -> tuple[dict[str, Any], dict[int, Path], dict[int, Path]]:
    summary_path = directory / "fertility_mechanism_summary.json"
    diagnostic = read_json(summary_path)
    production = read_json(directory / "summary.json")
    if diagnostic.get("status") != "complete_unpromoted_fertility_mechanism_diagnostic":
        raise RuntimeError(f"Incomplete state diagnostic: {directory}")
    if diagnostic.get("initial_path_sha256") != path_sha256:
        raise RuntimeError(f"Diagnostic path hash mismatch: {directory}")
    if sha256(directory / "summary.json") != diagnostic.get("production_summary_sha256"):
        raise RuntimeError(f"Diagnostic production hash mismatch: {directory}")
    if not bool(production.get("path_converged")):
        raise RuntimeError(f"Diagnostic did not reproduce convergence: {directory}")
    fertility: dict[int, Path] = {}
    cross: dict[int, Path] = {}
    for validation in diagnostic.get("validations") or []:
        year = int(validation["calendar_year"])
        fpath = directory / Path(validation["fertility_state_csv"]).name
        cpath = directory / Path(validation["cross_section_csv"]).name
        if sha256(fpath) != validation["fertility_state_csv_sha256"]:
            raise RuntimeError(f"Fertility state hash mismatch in {year}: {directory}")
        if sha256(cpath) != validation["cross_section_csv_sha256"]:
            raise RuntimeError(f"Cross-section hash mismatch in {year}: {directory}")
        if abs(float(validation["birth_reconstruction_gap"])) > 2.0e-10:
            raise RuntimeError(f"Birth reconstruction failed in {year}: {directory}")
        if abs(float(validation["housing_reconstruction_gap"])) > 2.0e-9:
            raise RuntimeError(f"Housing reconstruction failed in {year}: {directory}")
        fertility[year] = fpath
        cross[year] = cpath
    expected = set(int(year) for year in diagnostic.get("diagnostic_years") or [])
    if set(fertility) != expected or set(cross) != expected:
        raise RuntimeError(f"Missing diagnostic tables: {directory}")
    return diagnostic, fertility, cross


def paired_path_rows(base: list[dict[str, str]], reform: list[dict[str, str]]) -> list[dict[str, Any]]:
    if len(base) != len(reform):
        raise ValueError("Paired paths have different lengths")
    rows: list[dict[str, Any]] = []
    for left, right in zip(base, reform):
        year = int(left["calendar_year"])
        if year != int(right["calendar_year"]):
            raise ValueError("Paired paths have different calendars")
        base_birth = float(left["births_per_household_head"]) / PERIOD_YEARS
        reform_birth = float(right["births_per_household_head"]) / PERIOD_YEARS
        base_housing = float(left["housing_demand"]) / float(left["household_heads"])
        reform_housing = float(right["housing_demand"]) / float(right["household_heads"])
        rows.append(
            {
                "period": int(left["period"]),
                "calendar_year": year,
                "baseline_asset_price": float(left["asset_price"]),
                "reform_asset_price": float(right["asset_price"]),
                "asset_price_percent_change": 100.0 * (float(right["asset_price"]) / float(left["asset_price"]) - 1.0),
                "baseline_renter_price": float(left["renter_price"]),
                "reform_renter_price": float(right["renter_price"]),
                "renter_price_percent_change": 100.0 * (float(right["renter_price"]) / float(left["renter_price"]) - 1.0),
                "baseline_equal_rebate": float(left["equal_transfer_period_units"]),
                "reform_equal_rebate": float(right["equal_transfer_period_units"]),
                "equal_rebate_percent_change": 100.0 * (float(right["equal_transfer_period_units"]) / float(left["equal_transfer_period_units"]) - 1.0),
                "baseline_owner_rate": float(left["owner_rate"]),
                "reform_owner_rate": float(right["owner_rate"]),
                "owner_rate_pp_change": 100.0 * (float(right["owner_rate"]) - float(left["owner_rate"])),
                "baseline_housing_per_head": base_housing,
                "reform_housing_per_head": reform_housing,
                "housing_per_head_percent_change": 100.0 * (reform_housing / base_housing - 1.0),
                "baseline_annual_births_per_head": base_birth,
                "reform_annual_births_per_head": reform_birth,
                "annual_births_per_head_percent_change": 100.0 * (reform_birth / base_birth - 1.0),
                "annual_births_per_head_level_change": reform_birth - base_birth,
                "baseline_resident_persons": float(left["resident_persons"]),
                "reform_resident_persons": float(right["resident_persons"]),
                "resident_persons_percent_change": 100.0 * (float(right["resident_persons"]) / float(left["resident_persons"]) - 1.0),
                "baseline_market_residual": abs(float(left["relative_market_residual"])),
                "reform_market_residual": abs(float(right["relative_market_residual"])),
                "baseline_fiscal_residual": abs(float(left["government_budget_residual"])),
                "reform_fiscal_residual": abs(float(right["government_budget_residual"])),
            }
        )
    return rows


def iter_fertility(path: Path) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def key(row: dict[str, str]) -> tuple[int, ...]:
    return tuple(int(row[field]) for field in STATE_KEY)


def wealth_bin(index: int, maximum: int) -> str:
    fraction = index / max(maximum, 1)
    if fraction < 0.25:
        return "q1_grid"
    if fraction < 0.50:
        return "q2_grid"
    if fraction < 0.75:
        return "q3_grid"
    return "q4_grid"


def decompose_fertility_year(
    baseline_path: Path,
    reform_path: Path,
    *,
    baseline_heads: float,
    reform_heads: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    base_rows = iter_fertility(baseline_path)
    reform_rows = iter_fertility(reform_path)
    grouped: dict[tuple[str, str], dict[str, float]] = defaultdict(lambda: defaultdict(float))
    totals = defaultdict(float)
    maximum_wealth_index = max(
        int(row["wealth_index"]) for row in iter_fertility(baseline_path)
    )
    year: int | None = None
    sentinel = object()
    for left, right in zip_longest(base_rows, reform_rows, fillvalue=sentinel):
        if left is sentinel or right is sentinel:
            raise RuntimeError("Paired fertility state tables have different row counts")
        if key(left) != key(right):
            raise RuntimeError(f"State grids differ: {key(left)} versus {key(right)}")
        year = int(left["calendar_year"])
        mass0 = float(left["pre_fertility_mass"])
        mass1 = float(right["pre_fertility_mass"])
        if mass0 < -1.0e-15 or mass1 < -1.0e-15:
            raise RuntimeError(f"Negative pre-fertility mass in state {key(left)}")
        for label, row in (("baseline", left), ("reform", right)):
            for field in (
                "attempt_probability",
                "fecundity_probability",
                "realized_birth_probability",
            ):
                probability = float(row[field])
                if not math.isfinite(probability) or probability < -1.0e-14 or probability > 1.0 + 1.0e-14:
                    raise RuntimeError(
                        f"Invalid {label} {field} in state {key(row)}: {probability}"
                    )
            expected = float(row["expected_birth_child_units"])
            reconstructed = (
                float(row["realized_birth_probability"])
                * float(row["birth_child_units_if_realized"])
            )
            if not math.isclose(expected, reconstructed, rel_tol=0.0, abs_tol=2.0e-14):
                raise RuntimeError(
                    f"Birth-probability accounting failed for {label} state {key(row)}"
                )
        p0 = float(left["expected_birth_child_units"])
        p1 = float(right["expected_birth_child_units"])
        w0 = mass0 / baseline_heads
        w1 = mass1 / reform_heads
        policy = 0.5 * (w0 + w1) * (p1 - p0) / PERIOD_YEARS
        composition = 0.5 * (w1 - w0) * (p0 + p1) / PERIOD_YEARS
        totals["baseline_rate"] += w0 * p0 / PERIOD_YEARS
        totals["reform_rate"] += w1 * p1 / PERIOD_YEARS
        totals["policy"] += policy
        totals["composition"] += composition
        totals["at_risk_share_baseline"] += w0
        totals["at_risk_share_reform"] += w1
        gap0 = float(left["birth_value_gap_attempt_minus_wait"])
        gap1 = float(right["birth_value_gap_attempt_minus_wait"])
        if math.isfinite(gap0) and abs(gap0) <= 0.05:
            totals["near_margin_share_baseline"] += w0
        if math.isfinite(gap1) and abs(gap1) <= 0.05:
            totals["near_margin_share_reform"] += w1
        if math.isfinite(gap0) and abs(gap0) <= 0.25:
            totals["near_margin_0p25_share_baseline"] += w0
        if math.isfinite(gap1) and abs(gap1) <= 0.25:
            totals["near_margin_0p25_share_reform"] += w1
        if math.isfinite(gap0) and abs(gap0) <= 0.50:
            totals["near_margin_0p50_share_baseline"] += w0
        if math.isfinite(gap1) and abs(gap1) <= 0.50:
            totals["near_margin_0p50_share_reform"] += w1
        dimensions = {
            "age_group": left["age_group"],
            "origin_tenure": "owner" if int(left["origin_tenure"]) > 0 else "renter",
            "parity": f"parity_{int(left['parity'])}",
            "wealth_grid_quartile": wealth_bin(int(left["wealth_index"]), maximum_wealth_index),
        }
        for dimension, label in dimensions.items():
            item = grouped[(dimension, label)]
            item["policy"] += policy
            item["composition"] += composition
            item["total"] += policy + composition
    exact = totals["reform_rate"] - totals["baseline_rate"]
    if not math.isclose(totals["policy"] + totals["composition"], exact, rel_tol=0.0, abs_tol=5.0e-13):
        raise RuntimeError(
            f"Fertility decomposition does not add up: {totals['policy'] + totals['composition']} versus {exact}"
        )
    if year is None:
        raise RuntimeError(f"Empty fertility state table: {baseline_path}")
    group_rows = [
        {
            "calendar_year": year,
            "dimension": dimension,
            "group": label,
            "policy_effect_annual_births_per_head": values["policy"],
            "composition_effect_annual_births_per_head": values["composition"],
            "total_effect_annual_births_per_head": values["total"],
            "policy_effect_percent_of_baseline_rate": 100.0 * values["policy"] / totals["baseline_rate"],
            "composition_effect_percent_of_baseline_rate": 100.0 * values["composition"] / totals["baseline_rate"],
            "total_effect_percent_of_baseline_rate": 100.0 * values["total"] / totals["baseline_rate"],
        }
        for (dimension, label), values in sorted(grouped.items())
    ]
    summary = {
        "calendar_year": year,
        "baseline_annual_births_per_head": totals["baseline_rate"],
        "reform_annual_births_per_head": totals["reform_rate"],
        "total_effect_annual_births_per_head": exact,
        "total_effect_percent_of_baseline": 100.0 * exact / totals["baseline_rate"],
        "within_state_policy_effect": totals["policy"],
        "within_state_policy_effect_percent_of_baseline": 100.0 * totals["policy"] / totals["baseline_rate"],
        "distribution_composition_effect": totals["composition"],
        "distribution_composition_effect_percent_of_baseline": 100.0 * totals["composition"] / totals["baseline_rate"],
        "at_risk_mass_share_baseline": totals["at_risk_share_baseline"],
        "at_risk_mass_share_reform": totals["at_risk_share_reform"],
        "near_margin_mass_share_baseline_gap_le_0p05": totals["near_margin_share_baseline"],
        "near_margin_mass_share_reform_gap_le_0p05": totals["near_margin_share_reform"],
        "near_margin_mass_share_baseline_gap_le_0p25": totals["near_margin_0p25_share_baseline"],
        "near_margin_mass_share_reform_gap_le_0p25": totals["near_margin_0p25_share_reform"],
        "near_margin_mass_share_baseline_gap_le_0p50": totals["near_margin_0p50_share_baseline"],
        "near_margin_mass_share_reform_gap_le_0p50": totals["near_margin_0p50_share_reform"],
    }
    return summary, group_rows


def compare_cross_sections(base_path: Path, reform_path: Path) -> list[dict[str, Any]]:
    base = read_csv(base_path)
    reform = read_csv(reform_path)
    base_map = {(row["age_group"], row["parent_status"], row["tenure"]): row for row in base}
    reform_map = {(row["age_group"], row["parent_status"], row["tenure"]): row for row in reform}
    keys = sorted(set(base_map) | set(reform_map))
    year = int((base or reform)[0]["calendar_year"])
    rows: list[dict[str, Any]] = []
    for group in keys:
        left = base_map.get(group)
        right = reform_map.get(group)
        mass0 = float(left["mass"]) if left else 0.0
        mass1 = float(right["mass"]) if right else 0.0
        if mass0 < -1.0e-15 or mass1 < -1.0e-15:
            raise RuntimeError(f"Negative cross-section mass in {year}, {group}")
        h0 = float(left["mean_housing_services"]) if left else math.nan
        h1 = float(right["mean_housing_services"]) if right else math.nan
        rows.append(
            {
                "calendar_year": year,
                "age_group": group[0],
                "parent_status": group[1],
                "tenure": group[2],
                "baseline_mass": mass0,
                "reform_mass": mass1,
                "mass_change": mass1 - mass0,
                "baseline_mass_share": float(left["mass_share"]) if left else 0.0,
                "reform_mass_share": float(right["mass_share"]) if right else 0.0,
                "mass_share_pp_change": 100.0 * ((float(right["mass_share"]) if right else 0.0) - (float(left["mass_share"]) if left else 0.0)),
                "baseline_mean_housing_services": h0,
                "reform_mean_housing_services": h1,
                "housing_services_percent_change": 100.0 * (h1 / h0 - 1.0) if left and right and h0 else math.nan,
                "baseline_child_floor_binding_share": float(left["child_floor_binding_share"]) if left else math.nan,
                "reform_child_floor_binding_share": float(right["child_floor_binding_share"]) if right else math.nan,
                "child_floor_binding_pp_change": 100.0 * (float(right["child_floor_binding_share"]) - float(left["child_floor_binding_share"])) if left and right else math.nan,
                "baseline_mean_liquid_wealth": float(left["mean_liquid_wealth"]) if left else math.nan,
                "reform_mean_liquid_wealth": float(right["mean_liquid_wealth"]) if right else math.nan,
            }
        )
    return rows


def make_transition_plots(rows: list[dict[str, Any]], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    years = np.asarray([row["calendar_year"] for row in rows])
    fig, axes = plt.subplots(3, 2, figsize=(10.5, 10.0), constrained_layout=True)
    panels = (
        ("asset_price", "Asset price"),
        ("renter_price", "Rent-equivalent price"),
        ("owner_rate", "Ownership rate"),
        ("housing_per_head", "Housing services per head"),
        ("annual_births_per_head", "Annual births per head"),
        ("resident_persons", "Resident persons (model units)"),
    )
    for axis, (field, title) in zip(axes.flat, panels):
        axis.plot(years, [row[f"baseline_{field}"] for row in rows], label="Rebated 1%", lw=2.0)
        axis.plot(years, [row[f"reform_{field}"] for row in rows], label="Rebated 2%", lw=2.0)
        axis.set_title(title)
        axis.set_xlabel("Calendar year")
        axis.grid(alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    axes.flat[0].legend(frameon=False)
    fig.savefig(output / "paired_transition_levels.png", dpi=220, bbox_inches="tight")
    fig.savefig(output / "paired_transition_levels.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 6.8), constrained_layout=True)
    panels = (
        (("asset_price_percent_change", "renter_price_percent_change"), ("Asset price", "Rent-equivalent price"), "Percent"),
        (("owner_rate_pp_change",), ("Ownership",), "Percentage points"),
        (("housing_per_head_percent_change", "annual_births_per_head_percent_change"), ("Housing services/head", "Annual births/head"), "Percent"),
        (("resident_persons_percent_change",), ("Resident persons",), "Percent"),
    )
    for axis, (fields, labels, ylabel) in zip(axes.flat, panels):
        for field, label in zip(fields, labels):
            axis.plot(years, [row[field] for row in rows], label=label, lw=2.0)
        axis.axhline(0.0, color="0.45", lw=0.8)
        axis.set_ylabel(ylabel)
        axis.set_xlabel("Calendar year")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)
        axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(output / "paired_transition_effects.png", dpi=220, bbox_inches="tight")
    fig.savefig(output / "paired_transition_effects.pdf", bbox_inches="tight")
    plt.close(fig)


def make_decomposition_plot(rows: list[dict[str, Any]], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    age_rows = [row for row in rows if row["dimension"] == "age_group"]
    years = sorted({int(row["calendar_year"]) for row in age_rows})
    groups = ["under_25", "25_34", "35_44", "45_64", "65_plus"]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)
    selected = [year for year in (2023, 2035, 2051, 2079) if year in years]
    for axis, year in zip(axes.flat, selected):
        mapping = {row["group"]: row for row in age_rows if int(row["calendar_year"]) == year}
        policy = [mapping.get(group, {}).get("policy_effect_percent_of_baseline_rate", 0.0) for group in groups]
        composition = [mapping.get(group, {}).get("composition_effect_percent_of_baseline_rate", 0.0) for group in groups]
        x = np.arange(len(groups))
        axis.bar(x - 0.18, policy, width=0.36, label="Within-state fertility policy")
        axis.bar(x + 0.18, composition, width=0.36, label="Distribution/composition")
        axis.axhline(0.0, color="0.45", lw=0.8)
        axis.set_xticks(x, [group.replace("_", "–") for group in groups], rotation=25, ha="right")
        axis.set_title(str(year))
        axis.set_ylabel("Contribution (% of baseline annual births/head)")
        axis.grid(axis="y", alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.savefig(output / "fertility_age_decomposition.png", dpi=220, bbox_inches="tight")
    fig.savefig(output / "fertility_age_decomposition.pdf", bbox_inches="tight")
    plt.close(fig)


def make_incidence_plot(rows: list[dict[str, Any]], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected = [row for row in rows if int(row["calendar_year"]) in (2023, 2035, 2051, 2079)]
    labels = sorted({(row["age_group"], row["parent_status"]) for row in selected})
    compact: list[dict[str, Any]] = []
    for year in (2023, 2035, 2051, 2079):
        year_rows = [row for row in selected if int(row["calendar_year"]) == year]
        for age, parent in labels:
            cell = [row for row in year_rows if row["age_group"] == age and row["parent_status"] == parent]
            if not cell:
                continue
            base_total = sum(float(row["baseline_mass"]) for row in cell)
            reform_total = sum(float(row["reform_mass"]) for row in cell)
            base_owner = sum(float(row["baseline_mass"]) for row in cell if row["tenure"] == "owner")
            reform_owner = sum(float(row["reform_mass"]) for row in cell if row["tenure"] == "owner")
            base_h = sum(float(row["baseline_mass"]) * float(row["baseline_mean_housing_services"]) for row in cell) / max(base_total, 1e-15)
            reform_h = sum(float(row["reform_mass"]) * float(row["reform_mean_housing_services"]) for row in cell) / max(reform_total, 1e-15)
            compact.append({"calendar_year": year, "group": f"{age}:{parent}", "owner_pp": 100.0 * (reform_owner / reform_total - base_owner / base_total), "housing_pct": 100.0 * (reform_h / base_h - 1.0)})
    focus = [row for row in compact if row["group"].startswith(("25_34", "35_44"))]
    groups = sorted({row["group"] for row in focus})
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)
    for axis, field, title, ylabel in (
        (axes[0], "owner_pp", "Ownership incidence", "2% minus 1% (percentage points)"),
        (axes[1], "housing_pct", "Housing-service incidence", "2% minus 1% (percent)"),
    ):
        for group in groups:
            data = [row for row in focus if row["group"] == group]
            axis.plot([row["calendar_year"] for row in data], [row[field] for row in data], marker="o", label=group.replace(":", ", ").replace("_", "–"))
        axis.axhline(0.0, color="0.45", lw=0.8)
        axis.set_title(title)
        axis.set_ylabel(ylabel)
        axis.set_xlabel("Calendar year")
        axis.grid(alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    axes[1].legend(frameon=False, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.savefig(output / "young_household_incidence.png", dpi=220, bbox_inches="tight")
    fig.savefig(output / "young_household_incidence.pdf", bbox_inches="tight")
    plt.close(fig)


def make_child_floor_plot(rows: list[dict[str, Any]], output: Path) -> None:
    """Show whether the calibrated child-space floor is economically active."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    years = sorted({int(row["calendar_year"]) for row in rows})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)
    for axis, age in zip(axes, ("25_34", "35_44")):
        for tenure, linestyle in (("renter", "-"), ("owner", "--")):
            baseline: list[float] = []
            reform: list[float] = []
            for year in years:
                cell = [
                    row
                    for row in rows
                    if int(row["calendar_year"]) == year
                    and row["age_group"] == age
                    and row["parent_status"] == "dependent_children"
                    and row["tenure"] == tenure
                ]
                if not cell:
                    baseline.append(math.nan)
                    reform.append(math.nan)
                    continue
                row = cell[0]
                baseline.append(100.0 * float(row["baseline_child_floor_binding_share"]))
                reform.append(100.0 * float(row["reform_child_floor_binding_share"]))
            axis.plot(
                years,
                baseline,
                color="C0",
                linestyle=linestyle,
                marker="o",
                label=f"1% {tenure}",
            )
            axis.plot(
                years,
                reform,
                color="C1",
                linestyle=linestyle,
                marker="o",
                label=f"2% {tenure}",
            )
        axis.set_title(f"Parents age {age.replace('_', '–')}")
        axis.set_xlabel("Calendar year")
        axis.set_ylabel("Child-space floor binding (%)")
        axis.grid(alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    axes[1].legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.savefig(output / "child_space_floor_binding.png", dpi=220, bbox_inches="tight")
    fig.savefig(output / "child_space_floor_binding.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    base_dir = args.baseline_dir.resolve()
    reform_dir = args.reform_dir.resolve()
    calibration = preserve_calibration_tables(args.calibration_dir.resolve(), output)
    base_summary, base_path = validate_path(base_dir, "rebated-tax1-baseline")
    reform_summary, reform_path = validate_path(reform_dir, "rebated-tax2-reform")
    if base_summary.get("source_hashes") != reform_summary.get("source_hashes"):
        raise RuntimeError("Paired paths use different source hashes")
    if base_summary.get("tenure_choice_sensitivity") != reform_summary.get("tenure_choice_sensitivity"):
        raise RuntimeError("Paired paths use different tenure-choice sensitivities")
    paired = paired_path_rows(base_path, reform_path)
    base_diag, base_fertility, base_cross = diagnostic_tables(
        args.baseline_diagnostic_dir.resolve(), sha256(base_dir / "transition_path.csv")
    )
    reform_diag, reform_fertility, reform_cross = diagnostic_tables(
        args.reform_diagnostic_dir.resolve(), sha256(reform_dir / "transition_path.csv")
    )
    if set(base_fertility) != set(reform_fertility):
        raise RuntimeError("Diagnostic calendars differ")
    if base_diag["post2023_tenure_choice_kappa"] != reform_diag["post2023_tenure_choice_kappa"]:
        raise RuntimeError("Tenure sensitivities differ")
    if base_diag.get("source_hashes") != reform_diag.get("source_hashes"):
        raise RuntimeError("Paired state diagnostics use different source hashes")
    path_by_year = {int(row["calendar_year"]): row for row in paired}
    fertility_summaries: list[dict[str, Any]] = []
    decomposition_rows: list[dict[str, Any]] = []
    cross_rows: list[dict[str, Any]] = []
    for year in sorted(base_fertility):
        path_row = path_by_year[year]
        summary, groups = decompose_fertility_year(
            base_fertility[year],
            reform_fertility[year],
            baseline_heads=float(base_path[int(path_row["period"])]["household_heads"]),
            reform_heads=float(reform_path[int(path_row["period"])]["household_heads"]),
        )
        if not math.isclose(
            float(summary["baseline_annual_births_per_head"]),
            float(path_row["baseline_annual_births_per_head"]),
            rel_tol=0.0,
            abs_tol=5.0e-11,
        ) or not math.isclose(
            float(summary["reform_annual_births_per_head"]),
            float(path_row["reform_annual_births_per_head"]),
            rel_tol=0.0,
            abs_tol=5.0e-11,
        ):
            raise RuntimeError(f"State and path fertility rates differ in {year}")
        fertility_summaries.append(summary)
        decomposition_rows.extend(groups)
        cross_rows.extend(compare_cross_sections(base_cross[year], reform_cross[year]))
    write_csv(output / "paired_transition_decomposition.csv", paired)
    write_csv(output / "fertility_decomposition_summary.csv", fertility_summaries)
    write_csv(output / "fertility_decomposition_by_group.csv", decomposition_rows)
    write_csv(output / "cross_section_incidence.csv", cross_rows)
    make_transition_plots(paired, output)
    make_decomposition_plot(decomposition_rows, output)
    make_incidence_plot(cross_rows, output)
    make_child_floor_plot(cross_rows, output)
    peak_birth = max(paired, key=lambda row: abs(float(row["annual_births_per_head_percent_change"])))
    peak_owner = max(paired, key=lambda row: abs(float(row["owner_rate_pp_change"])))
    summary = {
        "status": "complete_audited_h128_property_tax_mechanism_decomposition",
        "scope": (
            "Paired accepted H128 paths plus an exact two-factor Shapley decomposition "
            "of annual births per household head into within-state policy and normalized "
            "household-distribution effects at selected dates."
        ),
        "tenure_choice_kappa": base_diag["post2023_tenure_choice_kappa"],
        "active_calibration": calibration,
        "baseline_path": str(base_dir),
        "reform_path": str(reform_dir),
        "input_hashes": {
            "baseline_summary": sha256(base_dir / "summary.json"),
            "baseline_path": sha256(base_dir / "transition_path.csv"),
            "reform_summary": sha256(reform_dir / "summary.json"),
            "reform_path": sha256(reform_dir / "transition_path.csv"),
            "baseline_diagnostic_summary": sha256(args.baseline_diagnostic_dir.resolve() / "fertility_mechanism_summary.json"),
            "reform_diagnostic_summary": sha256(args.reform_diagnostic_dir.resolve() / "fertility_mechanism_summary.json"),
        },
        "path_gates": {
            "market_gate": MARKET_GATE,
            "fiscal_gate": FISCAL_GATE,
            "baseline_maximum_market_residual": base_summary["maximum_market_residual"],
            "baseline_maximum_fiscal_residual": base_summary["maximum_fiscal_residual"],
            "reform_maximum_market_residual": reform_summary["maximum_market_residual"],
            "reform_maximum_fiscal_residual": reform_summary["maximum_fiscal_residual"],
        },
        "peak_absolute_birth_effect": {
            "calendar_year": peak_birth["calendar_year"],
            "percent_change": peak_birth["annual_births_per_head_percent_change"],
        },
        "peak_absolute_ownership_effect": {
            "calendar_year": peak_owner["calendar_year"],
            "percentage_point_change": peak_owner["owner_rate_pp_change"],
        },
        "fertility_decomposition": fertility_summaries,
    }
    write_json(output / "summary.json", summary)
    artifact_hashes = {
        path.name: sha256(path)
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "artifact_hashes.json"
    }
    write_json(output / "artifact_hashes.json", artifact_hashes)


if __name__ == "__main__":
    main()
