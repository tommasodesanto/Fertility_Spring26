#!/usr/bin/env python3
"""Build a Census-anchored demographic forecast with the E5F policy wedge."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from build_e5f_persons_demographic_satellite import read_csv, sha256, write_csv


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SATELLITE = ROOT / "output/model/e5f_persons_demographic_satellite_20260826b"
DEFAULT_HEADSHIP = ROOT / "output/model/e5f_headship_demographic_bridge_20260826a"
SCENARIOS = ("zero", "low", "mid", "hi")
SUMMARY_YEARS = (2025, 2050, 2080, 2100)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--headship-dir", type=Path, default=DEFAULT_HEADSHIP)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def make_figure(rows: list[dict[str, Any]], path: Path) -> None:
    by_scenario = {
        scenario: [row for row in rows if row["immigration_scenario"] == scenario]
        for scenario in SCENARIOS
    }
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.5))
    colors = {"zero": "#aaaaaa", "low": "#7aa6c2", "mid": "#26547c", "hi": "#17324d"}
    labels = {"zero": "Zero", "low": "Low", "mid": "Main", "hi": "High"}
    for scenario in SCENARIOS:
        selected = by_scenario[scenario]
        years = [int(row["calendar_year"]) for row in selected]
        axes[0].plot(
            years,
            [float(row["baseline_resident_persons_anchored"]) / 1e6 for row in selected],
            color=colors[scenario],
            linewidth=2.1 if scenario == "mid" else 1.25,
            label=labels[scenario],
        )
        axes[1].plot(
            years,
            [float(row["baseline_household_heads_anchored"]) / 1e6 for row in selected],
            color=colors[scenario],
            linewidth=2.1 if scenario == "mid" else 1.25,
            label=labels[scenario],
        )
    main = by_scenario["mid"]
    main_years = [int(row["calendar_year"]) for row in main]
    axes[0].plot(
        main_years,
        [float(row["reform_resident_persons_anchored"]) / 1e6 for row in main],
        color="#b22222",
        linewidth=1.5,
        linestyle="--",
        label="Main + policy wedge",
    )
    axes[1].plot(
        main_years,
        [float(row["reform_household_heads_anchored"]) / 1e6 for row in main],
        color="#b22222",
        linewidth=1.5,
        linestyle="--",
        label="Main + policy wedge",
    )
    axes[0].set(
        title="Resident-person forecast",
        xlabel="Calendar year",
        ylabel="Millions",
    )
    axes[1].set(
        title="Household-head forecast",
        xlabel="Calendar year",
        ylabel="Millions, ages 18--85",
    )
    axes[0].legend(frameon=False, fontsize=7.5)
    axes[1].legend(frameon=False, fontsize=7.5)

    axes[2].plot(
        main_years,
        [float(row["policy_extra_resident_persons"]) / 1e6 for row in main],
        color="#7aa6c2",
        linewidth=2.1,
        label="Resident persons",
    )
    axes[2].plot(
        main_years,
        [float(row["policy_extra_adult_persons_age_18_85"]) / 1e6 for row in main],
        color="#26547c",
        linewidth=1.8,
        label="Adult persons",
    )
    axes[2].plot(
        main_years,
        [float(row["policy_extra_household_heads"]) / 1e6 for row in main],
        color="#b22222",
        linewidth=2.0,
        label="Household heads",
    )
    axes[2].set(
        title="Rebated 2% minus rebated 1%",
        xlabel="Calendar year",
        ylabel="Millions",
    )
    axes[2].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    satellite = args.satellite_dir.resolve()
    headship = args.headship_dir.resolve()
    outdir = args.output_dir.resolve()
    for path in (satellite, headship):
        if not path.exists():
            raise FileNotFoundError(path)
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    satellite_certification = json.loads(
        (satellite / "certification.json").read_text(encoding="utf-8")
    )
    headship_certification = json.loads(
        (headship / "certification.json").read_text(encoding="utf-8")
    )
    population_anchor = float(
        satellite_certification["vintage_2025_resident_population"]
    )
    head_anchor = float(
        headship_certification[
            "vintage_2025_implied_heads_age_18_85_at_2023_rates"
        ]
    )
    person_rows = {
        (row["scenario"], int(row["calendar_year"])): row
        for row in read_csv(satellite / "annual_persons_policy_wedge.csv")
    }
    policy_head_rows = {
        (row["immigration_scenario"], int(row["calendar_year"])): row
        for row in read_csv(headship / "policy_person_to_head_wedge.csv")
    }
    baseline_head_rows = {
        (row["immigration_scenario"], int(row["calendar_year"])): row
        for row in read_csv(headship / "census_implied_head_forecasts.csv")
        if int(row["headship_vintage"]) == 2023
    }

    population_scales = {
        scenario: population_anchor
        / float(person_rows[(scenario, 2025)]["official_baseline_resident_population"])
        for scenario in SCENARIOS
    }
    head_scales = {
        scenario: head_anchor
        / float(
            baseline_head_rows[(scenario, 2025)][
                "implied_household_heads_age_18_85"
            ]
        )
        for scenario in SCENARIOS
    }
    years = sorted(
        year for scenario, year in person_rows if scenario == "mid"
    )
    rows: list[dict[str, Any]] = []
    for scenario in SCENARIOS:
        for year in years:
            persons = person_rows[(scenario, year)]
            policy_heads = policy_head_rows[(scenario, year)]
            baseline_heads = baseline_head_rows[(scenario, year)]
            baseline_person_level = (
                float(persons["official_baseline_resident_population"])
                * population_scales[scenario]
            )
            extra_persons = float(persons["extra_resident_persons"])
            baseline_head_level = (
                float(baseline_heads["implied_household_heads_age_18_85"])
                * head_scales[scenario]
            )
            extra_heads = float(policy_heads["extra_heads_2023_acs_headship"])
            rows.append(
                {
                    "immigration_scenario": scenario,
                    "calendar_year": year,
                    "baseline_resident_persons_anchored": baseline_person_level,
                    "reform_resident_persons_anchored": baseline_person_level
                    + extra_persons,
                    "policy_extra_resident_persons": extra_persons,
                    "policy_extra_children_age_0_17": float(
                        policy_heads["extra_children_age_0_17"]
                    ),
                    "policy_extra_adult_persons_age_18_85": float(
                        policy_heads["extra_adult_persons_age_18_85"]
                    ),
                    "baseline_household_heads_anchored": baseline_head_level,
                    "reform_household_heads_anchored": baseline_head_level
                    + extra_heads,
                    "policy_extra_household_heads": extra_heads,
                    "policy_extra_persons_pct_baseline": 100.0
                    * extra_persons
                    / baseline_person_level,
                    "policy_extra_heads_pct_baseline": 100.0
                    * extra_heads
                    / baseline_head_level,
                    "population_level_anchor_factor": population_scales[scenario],
                    "head_level_anchor_factor": head_scales[scenario],
                }
            )

    write_csv(outdir / "census_anchored_policy_forecast.csv", rows)
    make_figure(rows, outdir / "census_anchored_policy_forecast.png")
    start_rows = [row for row in rows if int(row["calendar_year"]) == 2025]
    gates = {
        "all_population_paths_share_2025_anchor": all(
            abs(float(row["baseline_resident_persons_anchored"]) - population_anchor)
            <= 1e-6
            for row in start_rows
        ),
        "all_head_paths_share_2025_anchor": all(
            abs(float(row["baseline_household_heads_anchored"]) - head_anchor)
            <= 1e-6
            for row in start_rows
        ),
        "reform_equals_baseline_plus_policy_wedge": all(
            abs(
                float(row["reform_resident_persons_anchored"])
                - float(row["baseline_resident_persons_anchored"])
                - float(row["policy_extra_resident_persons"])
            )
            <= 1e-6
            and abs(
                float(row["reform_household_heads_anchored"])
                - float(row["baseline_household_heads_anchored"])
                - float(row["policy_extra_household_heads"])
            )
            <= 1e-6
            for row in rows
        ),
    }
    if not all(gates.values()):
        raise RuntimeError(f"Census-anchored forecast gates failed: {gates}")

    certification = {
        "status": "complete_exploratory_census_anchored_policy_forecast_not_promoted",
        "population_anchor_2025": population_anchor,
        "household_head_anchor_2025_age_18_85": head_anchor,
        "population_anchor_factors": population_scales,
        "head_anchor_factors": head_scales,
        "accounting_gates": gates,
        "source_hashes": {
            "satellite_certification.json": sha256(satellite / "certification.json"),
            "annual_persons_policy_wedge.csv": sha256(
                satellite / "annual_persons_policy_wedge.csv"
            ),
            "headship_certification.json": sha256(headship / "certification.json"),
            "policy_person_to_head_wedge.csv": sha256(
                headship / "policy_person_to_head_wedge.csv"
            ),
            "census_implied_head_forecasts.csv": sha256(
                headship / "census_implied_head_forecasts.csv"
            ),
        },
        "limitations": [
            "Each Census 2023 projection scenario is proportionally level-anchored to Census Vintage 2025; this is not an official revised projection.",
            "The policy wedge changes births only; policy-induced migration and headship responses are zero.",
            "Frozen 2023 ACS headship rates are an accounting mapping, not a structural forecast.",
            "The demographic satellite does not feed back into housing equilibrium or household choices.",
            "The underlying perfect-foresight policy paths remain preliminary and unpromoted.",
        ],
    }
    (outdir / "certification.json").write_text(
        json.dumps(certification, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    selected = {
        (row["immigration_scenario"], int(row["calendar_year"])): row
        for row in rows
        if int(row["calendar_year"]) in SUMMARY_YEARS
    }
    readme = f"""# Exploratory Census-anchored E5F policy forecast

Status: `complete_exploratory_census_anchored_policy_forecast_not_promoted`.

All Census 2023 projection scenarios are proportionally level-anchored to the
latest Census Vintage 2025 resident-population estimate of
{population_anchor:,.0f}. Household-head paths are separately anchored to
{head_anchor / 1e6:.3f} million heads ages 18--85, obtained from the Vintage
2025 age-sex distribution and 2023 ACS headship rates. The rebated-2% policy
path adds only the model-implied fertility wedge; migration and headship are
held common across policy cases.

## Main immigration scenario

| Year | Baseline resident persons | Reform resident persons | Extra persons | Baseline heads | Reform heads | Extra heads |
|---|---:|---:|---:|---:|---:|---:|
"""
    for year in SUMMARY_YEARS:
        row = selected[("mid", year)]
        readme += (
            f"| {year} | {float(row['baseline_resident_persons_anchored']) / 1e6:.3f}m | "
            f"{float(row['reform_resident_persons_anchored']) / 1e6:.3f}m | "
            f"{float(row['policy_extra_resident_persons']):,.0f} | "
            f"{float(row['baseline_household_heads_anchored']) / 1e6:.3f}m | "
            f"{float(row['reform_household_heads_anchored']) / 1e6:.3f}m | "
            f"{float(row['policy_extra_household_heads']):,.0f} |\n"
        )
    readme += """

## Interpretation

This is the cleanest current demographic readout for the policy: use Census
for the national person and age structure, use ACS to map adults into household
heads, and use the economic model only for the proportional policy-induced
birth response. It avoids interpreting the model's adult-head mass as total
population. It remains a satellite rather than a demographic equilibrium and
does not update the policy's housing-price or rent response.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
