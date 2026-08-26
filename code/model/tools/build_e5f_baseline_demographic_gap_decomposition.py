#!/usr/bin/env python3
"""Decompose the E5F baseline household-head forecast gap.

This is an accounting diagnostic.  It anchors the active model's adult mass
to the latest ACS-implied number of household heads, compares that path with
Census population projections, and separates three sequential wedges:

1. replacing the model population law with the Census zero-immigration path
   while retaining the queue's coarse headship schedule;
2. replacing that schedule with 2023 ACS age- and sex-specific headship; and
3. moving from zero immigration to the Census main immigration scenario.

The decomposition does not change the equilibrium solver or promote the
underlying preliminary perfect-foresight paths.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from build_e5f_persons_demographic_satellite import read_csv, sha256, write_csv


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SATELLITE = ROOT / "output/model/e5f_persons_demographic_satellite_20260826b"
DEFAULT_HEADSHIP = ROOT / "output/model/e5f_headship_demographic_bridge_20260826a"
SUMMARY_YEARS = (2025, 2050, 2080, 2100)
SCENARIOS = ("zero", "low", "mid", "hi")
SEXES = (1, 2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--headship-dir", type=Path, default=DEFAULT_HEADSHIP)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=2025)
    parser.add_argument("--maximum-year", type=int, default=2100)
    return parser.parse_args()


def interpolate(rows: list[dict[str, str]], column: str, years: np.ndarray) -> np.ndarray:
    source_years = np.asarray([float(row["calendar_year"]) for row in rows])
    values = np.asarray([float(row[column]) for row in rows])
    if np.any(np.diff(source_years) <= 0.0):
        raise ValueError("Transition calendar years must be strictly increasing")
    if years[0] < source_years[0] or years[-1] > source_years[-1]:
        raise ValueError("Requested years exceed the transition-path support")
    return np.interp(years, source_years, values)


def census_total_births(path: Path) -> dict[int, float]:
    result = {
        int(row["YEAR"]): float(row["BIRTHS"])
        for row in read_csv(path)
        if row["RACE_HISP"] == "0" and row["SEX"] == "0"
    }
    if not result:
        raise RuntimeError(f"No national total-birth rows in {path}")
    return result


def census_queue_heads(path: Path) -> dict[int, float]:
    """Apply the active queue's implicit age-headship schedule to persons."""
    result: dict[int, float] = {}
    for row in read_csv(path):
        if row["ORIGIN"] != "0" or row["RACE"] != "0" or int(row["SEX"]) not in SEXES:
            continue
        year = int(row["YEAR"])
        result[year] = result.get(year, 0.0) + sum(
            float(row[f"POP_{age}"]) / 2.1 for age in range(20, 86)
        )
    if not result:
        raise RuntimeError(f"No sex-specific national population rows in {path}")
    return result


def census_acs_heads(path: Path) -> dict[str, dict[int, float]]:
    result: dict[str, dict[int, float]] = {scenario: {} for scenario in SCENARIOS}
    for row in read_csv(path):
        scenario = row["immigration_scenario"]
        if scenario not in result or int(row["headship_vintage"]) != 2023:
            continue
        result[scenario][int(row["calendar_year"])] = float(
            row["implied_household_heads_age_18_85"]
        )
    if any(not result[scenario] for scenario in SCENARIOS):
        raise RuntimeError("Census/ACS head forecast is incomplete")
    return result


def normalized_path(raw: dict[int, float], anchor: float, years: np.ndarray) -> np.ndarray:
    start = int(years[0])
    if start not in raw:
        raise RuntimeError(f"Missing {start} normalization year")
    return np.asarray([raw[int(year)] for year in years]) * anchor / raw[start]


def make_figure(rows: list[dict[str, Any]], path: Path) -> None:
    years = np.asarray([int(row["calendar_year"]) for row in rows])
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.5))

    axes[0].plot(
        years,
        [float(row["model_baseline_heads_anchored"]) / 1e6 for row in rows],
        color="#b22222",
        linewidth=2.1,
        label="Model baseline",
    )
    axes[0].plot(
        years,
        [float(row["census_zero_queue_heads_anchored"]) / 1e6 for row in rows],
        color="#777777",
        linewidth=1.5,
        linestyle="--",
        label="Census zero migration, queue mapping",
    )
    axes[0].plot(
        years,
        [float(row["census_zero_acs_heads_anchored"]) / 1e6 for row in rows],
        color="#7aa6c2",
        linewidth=1.7,
        label="Census zero migration, ACS mapping",
    )
    axes[0].plot(
        years,
        [float(row["census_mid_acs_heads_anchored"]) / 1e6 for row in rows],
        color="#26547c",
        linewidth=2.1,
        label="Census main migration, ACS mapping",
    )
    axes[0].set(
        title="Household-head paths anchored in 2025",
        xlabel="Calendar year",
        ylabel="Millions of heads ages 18--85",
    )
    axes[0].legend(frameon=False, fontsize=7.5)

    selected = [row for row in rows if int(row["calendar_year"]) in (2050, 2080, 2100)]
    x = np.arange(len(selected), dtype=float)
    width = 0.23
    for offset, column, label, color in (
        (-width, "population_law_gap", "Population law", "#999999"),
        (0.0, "headship_mapping_gap", "Headship mapping", "#7aa6c2"),
        (width, "immigration_gap", "Immigration", "#26547c"),
    ):
        axes[1].bar(
            x + offset,
            [float(row[column]) / 1e6 for row in selected],
            width=width,
            color=color,
            label=label,
        )
    axes[1].axhline(0.0, color="black", linewidth=0.7)
    axes[1].set_xticks(x, [str(int(row["calendar_year"])) for row in selected])
    axes[1].set(
        title="Sequential gap decomposition",
        xlabel="Calendar year",
        ylabel="Millions of heads",
    )
    axes[1].legend(frameon=False, fontsize=8)

    axes[2].plot(
        years,
        [
            float(row["model_topcode_adjusted_annual_births_anchored"]) / 1e6
            for row in rows
        ],
        color="#b22222",
        linewidth=2.1,
        label="Model baseline",
    )
    axes[2].plot(
        years,
        [float(row["census_zero_annual_births"]) / 1e6 for row in rows],
        color="#7aa6c2",
        linewidth=1.6,
        label="Census zero migration",
    )
    axes[2].plot(
        years,
        [float(row["census_mid_annual_births"]) / 1e6 for row in rows],
        color="#26547c",
        linewidth=2.0,
        label="Census main migration",
    )
    axes[2].set(
        title="Annual births",
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
    source_dir = satellite / "source_data"
    outdir = args.output_dir.resolve()
    for path in (source_dir, headship):
        if not path.exists():
            raise FileNotFoundError(path)
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    headship_certification = json.loads(
        (headship / "certification.json").read_text(encoding="utf-8")
    )
    expected_hashes = headship_certification["satellite_source_hashes"]
    for name, expected in expected_hashes.items():
        observed = sha256(source_dir / name)
        if observed != expected:
            raise RuntimeError(f"Source fingerprint differs for {name}")

    baseline_rows = read_csv(source_dir / "rebated_tax1_transition_path.csv")
    reform_rows = read_csv(source_dir / "rebated_tax2_transition_path.csv")
    years = np.arange(int(args.start_year), int(args.maximum_year) + 1)
    anchor_heads = float(
        headship_certification[
            "vintage_2025_implied_heads_age_18_85_at_2023_rates"
        ]
    )
    model_baseline_mass = interpolate(baseline_rows, "adult_population", years)
    model_reform_mass = interpolate(reform_rows, "adult_population", years)
    model_births = interpolate(
        baseline_rows, "birth_children_topcode_adjusted", years
    )
    scale = anchor_heads / model_baseline_mass[0]
    model_baseline_heads = model_baseline_mass * scale
    model_reform_heads = model_reform_mass * scale
    model_annual_births = model_births * scale / 4.0

    acs_raw = census_acs_heads(headship / "census_implied_head_forecasts.csv")
    acs_paths = {
        scenario: normalized_path(acs_raw[scenario], anchor_heads, years)
        for scenario in SCENARIOS
    }
    queue_zero_raw = census_queue_heads(source_dir / "population_zero.csv")
    queue_zero = normalized_path(queue_zero_raw, anchor_heads, years)
    census_births = {
        scenario: census_total_births(source_dir / f"births_{scenario}.csv")
        for scenario in ("zero", "mid")
    }

    output: list[dict[str, Any]] = []
    for index, year in enumerate(years):
        population_law_gap = queue_zero[index] - model_baseline_heads[index]
        headship_mapping_gap = acs_paths["zero"][index] - queue_zero[index]
        immigration_gap = acs_paths["mid"][index] - acs_paths["zero"][index]
        total_gap = acs_paths["mid"][index] - model_baseline_heads[index]
        output.append(
            {
                "calendar_year": int(year),
                "model_baseline_heads_anchored": model_baseline_heads[index],
                "model_reform_heads_anchored": model_reform_heads[index],
                "model_policy_wedge_heads": model_reform_heads[index]
                - model_baseline_heads[index],
                "model_topcode_adjusted_annual_births_anchored": model_annual_births[
                    index
                ],
                "census_zero_annual_births": census_births["zero"][int(year)],
                "census_mid_annual_births": census_births["mid"][int(year)],
                "census_zero_queue_heads_anchored": queue_zero[index],
                "census_zero_acs_heads_anchored": acs_paths["zero"][index],
                "census_low_acs_heads_anchored": acs_paths["low"][index],
                "census_mid_acs_heads_anchored": acs_paths["mid"][index],
                "census_high_acs_heads_anchored": acs_paths["hi"][index],
                "population_law_gap": population_law_gap,
                "headship_mapping_gap": headship_mapping_gap,
                "immigration_gap": immigration_gap,
                "total_model_to_census_main_gap": total_gap,
                "decomposition_residual": total_gap
                - population_law_gap
                - headship_mapping_gap
                - immigration_gap,
            }
        )

    write_csv(outdir / "baseline_demographic_gap_decomposition.csv", output)
    make_figure(output, outdir / "baseline_demographic_gap_decomposition.png")

    gates = {
        "model_2025_anchor_exact": bool(
            abs(model_baseline_heads[0] - anchor_heads) <= 1e-6
        ),
        "decomposition_identity": bool(max(
            abs(float(row["decomposition_residual"])) for row in output
        )
        <= 1e-6),
        "all_head_paths_positive": bool(min(
            min(float(row[column]) for row in output)
            for column in (
                "model_baseline_heads_anchored",
                "census_zero_queue_heads_anchored",
                "census_zero_acs_heads_anchored",
                "census_mid_acs_heads_anchored",
            )
        )
        > 0.0),
    }
    if not all(gates.values()):
        raise RuntimeError(f"Baseline demographic accounting gates failed: {gates}")

    selected = {int(row["calendar_year"]): row for row in output if int(row["calendar_year"]) in SUMMARY_YEARS}
    certification = {
        "status": "complete_exploratory_baseline_gap_decomposition_not_promoted",
        "anchor": {
            "year": int(args.start_year),
            "heads_age_18_85": anchor_heads,
            "definition": "Census Vintage 2025 resident persons times 2023 ACS age-sex headship rates",
        },
        "model_scale_per_adult_mass_unit": scale,
        "source_hashes": {
            "headship_certification.json": sha256(headship / "certification.json"),
            "census_implied_head_forecasts.csv": sha256(
                headship / "census_implied_head_forecasts.csv"
            ),
            **{
                name: sha256(source_dir / name)
                for name in (
                    "rebated_tax1_transition_path.csv",
                    "rebated_tax2_transition_path.csv",
                    "population_zero.csv",
                    "births_zero.csv",
                    "births_mid.csv",
                )
            },
        },
        "accounting_gates": gates,
        "limitations": [
            "The wedges are a sequential accounting decomposition and depend on the displayed ordering.",
            "The population-law wedge bundles the model's initial age state, birth flow, four-year timing, deterministic terminal exit, and zero-migration closure; it is not a mortality estimate.",
            "Frozen 2023 ACS headship rates are not a structural forecast of household formation.",
            "The Census scenarios are projection scenarios, not confidence intervals or current point forecasts.",
            "The underlying perfect-foresight paths remain preliminary and unpromoted.",
        ],
    }
    (outdir / "certification.json").write_text(
        json.dumps(certification, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    readme = f"""# Exploratory E5F baseline demographic-gap decomposition

Status: `complete_exploratory_baseline_gap_decomposition_not_promoted`.

The model's 2025 adult mass is anchored to {anchor_heads / 1e6:.3f} million
household heads ages 18--85, implied by the latest Census Vintage 2025
resident-person age distribution and 2023 ACS age- and sex-specific headship
rates. All comparison paths share that 2025 level.

## Sequential accounting decomposition

| Year | Model heads | Census zero, queue mapping | Census zero, ACS mapping | Census main, ACS mapping | Total gap |
|---|---:|---:|---:|---:|---:|
"""
    for year in SUMMARY_YEARS:
        row = selected[year]
        readme += (
            f"| {year} | {float(row['model_baseline_heads_anchored']) / 1e6:.3f}m | "
            f"{float(row['census_zero_queue_heads_anchored']) / 1e6:.3f}m | "
            f"{float(row['census_zero_acs_heads_anchored']) / 1e6:.3f}m | "
            f"{float(row['census_mid_acs_heads_anchored']) / 1e6:.3f}m | "
            f"{float(row['total_model_to_census_main_gap']) / 1e6:.3f}m |\n"
        )
    readme += """

| Year | Population-law wedge | Headship-mapping wedge | Immigration wedge |
|---|---:|---:|---:|
"""
    for year in (2050, 2080, 2100):
        row = selected[year]
        readme += (
            f"| {year} | {float(row['population_law_gap']) / 1e6:+.3f}m | "
            f"{float(row['headship_mapping_gap']) / 1e6:+.3f}m | "
            f"{float(row['immigration_gap']) / 1e6:+.3f}m |\n"
        )
    first = selected[int(args.start_year)]
    readme += f"""

## Diagnosis

The baseline miss is already visible in the flow target. Annualizing and
scaling the model's top-code-adjusted 2025 birth flow gives
{float(first['model_topcode_adjusted_annual_births_anchored']) / 1e6:.3f} million births, versus
{float(first['census_mid_annual_births']) / 1e6:.3f} million in the Census main
projection. The calibration targets completed fertility, but it does not
directly target this period birth flow.

The queue mapping is also not a headship law. It assigns zero headship through
age 19 and `1 / 2.1` at every age from 20 through 85. The ACS profile starts
much lower for young adults and rises above that value later. More importantly,
the household-only state omits adults who are nonheads in 2025 and may form
households later. The model therefore cannot generate their subsequent entry
into headship.

The population-law wedge deliberately remains bundled: it includes the initial
model age state, the untargeted period birth flow, four-year cohort timing,
deterministic exit at the terminal age, and zero migration. The next structural
step is an explicit person cohort law with survival, migration, and an
age-specific household-formation transition. This diagnostic does not alter
the active solver or the certified baseline.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
