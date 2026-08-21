#!/usr/bin/env python3
"""Build a compact diagnosis of calibration misses and housing-policy responses."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-fit", type=Path, required=True)
    parser.add_argument("--identification-summary", type=Path, required=True)
    parser.add_argument("--mechanism-effects", type=Path, required=True)
    parser.add_argument("--unrebated-tax-effects", type=Path, required=True)
    parser.add_argument("--rebated-tax-effects", type=Path, required=True)
    parser.add_argument("--child-inclusive-path", type=Path)
    parser.add_argument("--rebated-tax-shapley", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def row_at(
    rows: list[dict[str, str]], case: str, year: int
) -> dict[str, str]:
    found = [
        row
        for row in rows
        if row["policy_case"] == case and int(row["calendar_year"]) == year
    ]
    if len(found) != 1:
        raise RuntimeError(f"Expected one {case} row in {year}; found {len(found)}")
    return found[0]


def effect_row(
    *,
    case: str,
    label: str,
    year: int,
    row: dict[str, str],
    horizon: str,
) -> dict[str, Any]:
    return {
        "policy_case": case,
        "policy_label": label,
        "calendar_year": year,
        "horizon": horizon,
        "house_price_percent": float(row["percent_difference_asset_price"]),
        "housing_user_cost_percent": (
            float(row["percent_difference_housing_user_cost"])
            if row.get("percent_difference_housing_user_cost") not in (None, "")
            else math.nan
        ),
        "births_per_adult_percent": float(
            row["percent_difference_topcode_adjusted_births_per_adult"]
        ),
        "owner_rate_pp": 100.0 * float(row["difference_owner_rate"]),
        "dependent_child_owner_rate_pp": 100.0
        * float(row["difference_dependent_child_owner_rate"]),
        "adult_household_population_percent": float(
            row["percent_difference_population_index_2023"]
        ),
    }


def finite_float(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise RuntimeError(f"Nonfinite {field}: {row[field]!r}")
    return value


def child_population_summary(path: Path) -> dict[str, Any]:
    rows = read_csv(path)
    years = [int(row["calendar_year"]) for row in rows]
    expected_years = list(range(2023, 2064, 4))
    if years != expected_years:
        raise RuntimeError(
            f"Child-inclusive path must contain {expected_years}; found {years}"
        )
    for row in rows:
        adult = finite_float(row, "adult_household_units")
        children = finite_float(row, "dependent_child_units")
        total = finite_float(row, "model_population_units")
        if min(adult, children, total) < 0.0:
            raise RuntimeError("Child-inclusive population units must be nonnegative")
        if not math.isclose(total, adult + children, rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError("Model population does not equal adults plus children")
    base = rows[0]

    def indexed(year: int) -> dict[str, float]:
        row = next(item for item in rows if int(item["calendar_year"]) == year)
        return {
            "adult_households": finite_float(row, "adult_household_units")
            / finite_float(base, "adult_household_units"),
            "dependent_children": finite_float(row, "dependent_child_units")
            / finite_float(base, "dependent_child_units"),
            "model_population": finite_float(row, "model_population_units")
            / finite_float(base, "model_population_units"),
            "house_price": finite_float(row, "asset_price")
            / finite_float(base, "asset_price"),
            "births_per_adult": finite_float(
                row, "topcode_adjusted_births_per_adult"
            )
            / finite_float(base, "topcode_adjusted_births_per_adult"),
            "dependent_child_share": finite_float(row, "dependent_child_units")
            / finite_float(row, "model_population_units"),
        }

    return {
        "accounting_unit": (
            "one adult household unit plus dependent children represented in the "
            "household child-count state"
        ),
        "years": years,
        "2023_dependent_child_share": finite_float(base, "dependent_child_units")
        / finite_float(base, "model_population_units"),
        "2035_index": indexed(2035),
        "2063_index": indexed(2063),
    }


def shapley_summary(
    path: Path, rebated_impact: dict[str, Any]
) -> dict[str, dict[str, float]]:
    rows = read_csv(path)
    factors = {"tax_rate", "asset_price", "equal_rebate"}
    metrics = {"births_per_adult", "dependent_child_owner_rate"}
    selected = [row for row in rows if row["metric"] in metrics]
    if len(selected) != 6:
        raise RuntimeError("Expected six tax Shapley rows for births and family ownership")
    result: dict[str, dict[str, float]] = {}
    for metric in metrics:
        metric_rows = [row for row in selected if row["metric"] == metric]
        if {row["factor"] for row in metric_rows} != factors:
            raise RuntimeError(f"Incomplete Shapley factor set for {metric}")
        contributions = {
            row["factor"]: finite_float(row, "reported_contribution")
            for row in metric_rows
        }
        result[metric] = contributions
    if not math.isclose(
        sum(result["births_per_adult"].values()),
        rebated_impact["births_per_adult_percent"],
        rel_tol=0.0,
        abs_tol=1e-8,
    ):
        raise RuntimeError("Birth Shapley contributions do not add to the reform effect")
    if not math.isclose(
        sum(result["dependent_child_owner_rate"].values()),
        rebated_impact["dependent_child_owner_rate_pp"],
        rel_tol=0.0,
        abs_tol=1e-8,
    ):
        raise RuntimeError(
            "Ownership Shapley contributions do not add to the reform effect"
        )
    return result


def main() -> None:
    args = parse_args()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    target_rows = read_csv(args.target_fit)
    if len(target_rows) != 12 or len({row["moment"] for row in target_rows}) != 12:
        raise RuntimeError("The calibration diagnosis requires the complete 12-row target fit")
    misses = sorted(
        (
            {
                "moment": row["moment"],
                "target": float(row["target"]),
                "model": float(row["model"]),
                "standardized_gap": float(row["standardized_gap"]),
                "loss_contribution": float(row["loss_contribution"]),
                "share_of_total_loss": 0.0,
            }
            for row in target_rows
        ),
        key=lambda row: row["loss_contribution"],
        reverse=True,
    )
    total_loss = sum(row["loss_contribution"] for row in misses)
    for row in misses:
        row["share_of_total_loss"] = row["loss_contribution"] / total_loss
    write_csv(outdir / "calibration_miss_decomposition.csv", misses)

    identification = json.loads(args.identification_summary.read_text(encoding="utf-8"))
    if not math.isclose(
        float(identification["anchor_loss"]), total_loss, rel_tol=0.0, abs_tol=1e-8
    ):
        raise RuntimeError("Identification anchor loss does not match the target ledger")

    mechanism = read_csv(args.mechanism_effects)
    unrebated = read_csv(args.unrebated_tax_effects)
    rebated = read_csv(args.rebated_tax_effects)
    impact = [
        effect_row(
            case="supply-plus-20",
            label="Housing supply +20%",
            year=2023,
            row=row_at(mechanism, "supply-plus-20", 2023),
            horizon="impact",
        ),
        effect_row(
            case="dependent-child-ltv95",
            label="Family LTV 95%",
            year=2023,
            row=row_at(mechanism, "dependent-child-ltv95", 2023),
            horizon="impact",
        ),
        effect_row(
            case="combined",
            label="Supply + family credit",
            year=2023,
            row=row_at(mechanism, "combined", 2023),
            horizon="impact",
        ),
        effect_row(
            case="property-tax-2pct-no-rebate",
            label="Property tax 1% to 2%, unrebated",
            year=2023,
            row=row_at(unrebated, "property-tax-2pct-no-rebate", 2023),
            horizon="impact",
        ),
    ]
    rebated_matches = [
        row for row in rebated if int(row["calendar_year"]) == 2023
    ]
    if len(rebated_matches) != 1:
        raise RuntimeError(
            f"Expected one rebated-tax effect row in 2023; found {len(rebated_matches)}"
        )
    rebated_2023 = rebated_matches[0]
    impact.append(
        {
            "policy_case": "rebated-tax2-reform",
            "policy_label": "Property tax 1% to 2%, equal rebate",
            "calendar_year": 2023,
            "horizon": "impact",
            "house_price_percent": float(rebated_2023["asset_price_percent_change"]),
            "housing_user_cost_percent": float(
                rebated_2023["housing_user_cost_percent_change"]
            ),
            "births_per_adult_percent": float(
                rebated_2023["births_per_adult_percent_change"]
            ),
            "owner_rate_pp": float(rebated_2023["owner_rate_pp_change"]),
            "dependent_child_owner_rate_pp": float(
                rebated_2023["dependent_child_owner_rate_pp_change"]
            ),
            "adult_household_population_percent": float(
                rebated_2023["adult_population_percent_change"]
            ),
        }
    )
    terminal = [
        effect_row(
            case=case,
            label=label,
            year=2063,
            row=row_at(mechanism, case, 2063),
            horizon="2063",
        )
        for case, label in (
            ("supply-plus-20", "Housing supply +20%"),
            ("dependent-child-ltv95", "Family LTV 95%"),
            ("combined", "Supply + family credit"),
        )
    ]
    write_csv(outdir / "policy_effect_summary.csv", impact + terminal)

    child_summary = (
        child_population_summary(args.child_inclusive_path)
        if args.child_inclusive_path is not None
        else None
    )
    rebated_impact = next(
        row for row in impact if row["policy_case"] == "rebated-tax2-reform"
    )
    shapley = (
        shapley_summary(args.rebated_tax_shapley, rebated_impact)
        if args.rebated_tax_shapley is not None
        else None
    )

    supply = next(row for row in impact if row["policy_case"] == "supply-plus-20")
    tax = next(
        row for row in impact if row["policy_case"] == "property-tax-2pct-no-rebate"
    )
    supply_elasticity = supply["births_per_adult_percent"] / supply["house_price_percent"]
    tax_elasticity = tax["births_per_adult_percent"] / tax["housing_user_cost_percent"]
    top_two_share = misses[0]["share_of_total_loss"] + misses[1]["share_of_total_loss"]
    best_coordinate_loss = float(identification["best_coordinate_loss"])
    if best_coordinate_loss > total_loss + 1e-10:
        raise RuntimeError("Best coordinate loss exceeds the anchor loss")

    labels = [row["policy_label"] for row in impact]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 3.8), constrained_layout=True)
    colors = ["#2c7fb8", "#d95f0e", "#31a354", "#756bb1", "#b0493d"]
    axes[0].barh(
        labels,
        [row["births_per_adult_percent"] for row in impact],
        color=colors,
    )
    axes[0].axvline(0.0, color="#777777", linewidth=0.8)
    axes[0].set_title("Births per adult household")
    axes[0].set_xlabel("Impact effect (%)")
    axes[1].barh(
        labels,
        [row["dependent_child_owner_rate_pp"] for row in impact],
        color=colors,
    )
    axes[1].axvline(0.0, color="#777777", linewidth=0.8)
    axes[1].set_title("Ownership among dependent-child households")
    axes[1].set_xlabel("Impact effect (percentage points)")
    for axis in axes:
        axis.grid(axis="x", alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(outdir / "policy_mechanism_diagnosis.pdf", bbox_inches="tight")
    fig.savefig(outdir / "policy_mechanism_diagnosis.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    diagnosis = {
        "status": "complete",
        "strict_loss": total_loss,
        "top_two_miss_share": top_two_share,
        "top_two_misses": [row["moment"] for row in misses[:2]],
        "identification_gate_passed": bool(identification["identification_gate_passed"]),
        "full_jacobian_condition_number": float(
            identification["full_rank"]["condition_number"]
        ),
        "active_jacobian_condition_number": float(
            identification["active_after_side_gate"]["condition_number"]
        ),
        "best_one_coordinate_loss": best_coordinate_loss,
        "best_one_coordinate_loss_improvement": total_loss - best_coordinate_loss,
        "impact_fertility_elasticity_to_supply_induced_price_change": supply_elasticity,
        "impact_fertility_elasticity_to_unrebated_tax_user_cost_change": tax_elasticity,
        "child_inclusive_forward_path": child_summary,
        "rebated_property_tax_shapley": shapley,
        "interpretation": (
            "The two distinct cost experiments imply a short-run fertility elasticity "
            "near -0.06. The credit intervention primarily changes tenure, not fertility. "
            "The objective is dominated by the first-birth rooms response and mean rooms; "
            "no active target directly identifies a causal housing-cost elasticity of fertility."
        ),
        "source_hashes": {
            "target_fit": sha256(args.target_fit),
            "identification_summary": sha256(args.identification_summary),
            "mechanism_effects": sha256(args.mechanism_effects),
            "unrebated_tax_effects": sha256(args.unrebated_tax_effects),
            "rebated_tax_effects": sha256(args.rebated_tax_effects),
            **(
                {"child_inclusive_path": sha256(args.child_inclusive_path)}
                if args.child_inclusive_path is not None
                else {}
            ),
            **(
                {"rebated_tax_shapley": sha256(args.rebated_tax_shapley)}
                if args.rebated_tax_shapley is not None
                else {}
            ),
        },
    }
    (outdir / "diagnosis.json").write_text(
        json.dumps(diagnosis, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    readme = f"""# Calibration and policy diagnosis

The active twelve-moment objective has loss `{total_loss:.6f}`. The first-birth
rooms response and mean occupied rooms account for `{100.0 * top_two_share:.1f}%`
of that loss: the model produces `{misses[0]['model']:.3f}` against a target of
`{misses[0]['target']:.3f}` for the first-birth response and
`{misses[1]['model']:.3f}` against `{misses[1]['target']:.3f}` for mean rooms.
The best one-coordinate perturbation lowers the loss by only
`{total_loss - best_coordinate_loss:.3f}`. The local eleven-parameter Jacobian
fails the declared rank gate;
after freezing the inconsistent column its condition number is
`{diagnosis['active_jacobian_condition_number']:.1f}`.

Two independent cost experiments give almost the same impact elasticity of
births per adult household: `{supply_elasticity:.3f}` for the supply-induced
house-price decline and `{tax_elasticity:.3f}` for the unrebated property-tax
increase in housing user cost. The weak fertility response is therefore a
feature of the fitted household block, not an artifact of one fiscal closure.

The family-credit intervention is substantively larger for tenure: it raises
ownership among households with dependent children by
`{next(row for row in impact if row['policy_case'] == 'dependent-child-ltv95')['dependent_child_owner_rate_pp']:.2f}`
percentage points on impact while barely changing births. The combined supply
and credit intervention raises the same ownership rate by
`{next(row for row in impact if row['policy_case'] == 'combined')['dependent_child_owner_rate_pp']:.2f}`
points on impact and
`{next(row for row in terminal if row['policy_case'] == 'combined')['dependent_child_owner_rate_pp']:.2f}`
points in 2063. These are tenure-allocation results, not large demographic effects.
"""
    if child_summary is not None:
        readme += f"""

## Forward population accounting

Population is measured as one adult household unit plus the dependent children
represented in that household. Relative to 2023, this measure is
`{child_summary['2035_index']['model_population']:.3f}` in 2035 and
`{child_summary['2063_index']['model_population']:.3f}` in 2063. The corresponding
adult-household indices are `{child_summary['2035_index']['adult_households']:.3f}`
and `{child_summary['2063_index']['adult_households']:.3f}`; the child indices are
`{child_summary['2035_index']['dependent_children']:.3f}` and
`{child_summary['2063_index']['dependent_children']:.3f}`. The decline in the child
stock therefore begins before the decline in adult household mass.
"""
    if shapley is not None:
        births_shapley = shapley["births_per_adult"]
        ownership_shapley = shapley["dependent_child_owner_rate"]
        readme += f"""

## Rebated property tax

The impact increase in births under the rebated tax is a transfer result, not a
pure housing-cost effect. The exact Shapley decomposition of the
`{rebated_impact['births_per_adult_percent']:.3f}%` net change assigns
`{births_shapley['tax_rate']:.3f}%` of baseline births to the higher tax rate,
`{births_shapley['asset_price']:.3f}%` to house-price capitalization, and
`{births_shapley['equal_rebate']:.3f}%` to the equal rebate. For ownership among
dependent-child households, the corresponding contributions are
`{ownership_shapley['tax_rate']:.3f}`, `{ownership_shapley['asset_price']:.3f}`,
and `{ownership_shapley['equal_rebate']:.3f}` percentage points. The direct tax
effect is negative; redistribution more than offsets it.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")
    print(f"POLICY_CALIBRATION_DIAGNOSIS_COMPLETE output={outdir}")


if __name__ == "__main__":
    main()
