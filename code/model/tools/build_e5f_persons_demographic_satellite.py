#!/usr/bin/env python3
"""Build a persons-based demographic satellite for an E5F policy comparison.

The active E5F transition state is a mass of adult household heads.  This
diagnostic does not change that equilibrium.  It maps the model's proportional
change in aggregate births into a separate resident-person accounting layer:

1. apply the model birth wedge to Census projected births;
2. count the additional newborns immediately; and
3. age each additional birth cohort with Census native-born survival ratios.

Net international migration follows the same Census scenario in the policy and
baseline cases, so the policy-induced migration response is fixed to zero.  The
result is a transparent policy-population wedge, not a structural U.S. baseline
forecast and not a promoted model result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import urllib.request
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCENARIOS = ("zero", "low", "mid", "hi")
CENSUS_ROOT = (
    "https://www2.census.gov/programs-surveys/popproj/datasets/2023/"
    "2023-popproj"
)
CENSUS_SOURCES = {
    "survival": f"{CENSUS_ROOT}/np2023_a4.csv",
    **{
        f"population_{scenario}": f"{CENSUS_ROOT}/np2023_d1_{scenario}.csv"
        for scenario in SCENARIOS
    },
    **{
        f"births_{scenario}": f"{CENSUS_ROOT}/np2023_d2_{scenario}.csv"
        for scenario in SCENARIOS
    },
    "vintage_2025_age_sex": (
        "https://www2.census.gov/programs-surveys/popest/datasets/2020-2025/"
        "national/asrh/nc-est2025-agesex-res.csv"
    ),
}
SUMMARY_YEARS = (2050, 2080, 2100)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-path", type=Path, required=True)
    parser.add_argument("--reform-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--policy-start-year", type=int, default=2025)
    parser.add_argument("--maximum-year", type=int, default=2100)
    parser.add_argument("--baseline-sha256")
    parser.add_argument("--reform-sha256")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def copy_input(path: Path, destination: Path) -> None:
    with path.open("rb") as source, destination.open("xb") as target:
        shutil.copyfileobj(source, target)


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "research-replication"})
    with urllib.request.urlopen(request, timeout=60) as response:
        if int(response.status) != 200:
            raise RuntimeError(f"Download failed with HTTP {response.status}: {url}")
        with destination.open("xb") as handle:
            shutil.copyfileobj(response, handle)


def validate_hash(path: Path, expected: str | None, label: str) -> str:
    actual = sha256(path)
    if expected is not None and actual != expected:
        raise RuntimeError(
            f"{label} hash mismatch: expected={expected}, actual={actual}, path={path}"
        )
    return actual


def policy_birth_knots(
    baseline_rows: list[dict[str, str]],
    reform_rows: list[dict[str, str]],
) -> list[dict[str, float | int]]:
    if len(baseline_rows) != len(reform_rows):
        raise RuntimeError("Policy paths have different lengths")
    output: list[dict[str, float | int]] = []
    for baseline, reform in zip(baseline_rows, reform_rows):
        period = int(baseline["period"])
        if period != int(reform["period"]):
            raise RuntimeError("Policy paths have different period indexes")
        baseline_births = float(baseline["birth_children_topcode_adjusted"])
        reform_births = float(reform["birth_children_topcode_adjusted"])
        if baseline_births <= 0.0 or reform_births <= 0.0:
            raise RuntimeError("Policy path contains nonpositive adjusted births")
        output.append(
            {
                "period": period,
                "elapsed_years": 4 * period,
                "source_calendar_year": int(baseline["calendar_year"]),
                "baseline_adjusted_births_model_units": baseline_births,
                "reform_adjusted_births_model_units": reform_births,
                "aggregate_birth_effect": reform_births / baseline_births - 1.0,
                "aggregate_birth_effect_pct": 100.0
                * (reform_births / baseline_births - 1.0),
            }
        )
    return output


def interpolate_effect(knots: list[dict[str, float | int]], elapsed_years: int) -> float:
    elapsed = float(elapsed_years)
    if elapsed <= float(knots[0]["elapsed_years"]):
        return float(knots[0]["aggregate_birth_effect"])
    for lower, upper in zip(knots, knots[1:]):
        x0 = float(lower["elapsed_years"])
        x1 = float(upper["elapsed_years"])
        if x0 <= elapsed <= x1:
            weight = (elapsed - x0) / (x1 - x0)
            return (1.0 - weight) * float(lower["aggregate_birth_effect"]) + weight * float(
                upper["aggregate_birth_effect"]
            )
    raise RuntimeError(
        f"Policy effect would require extrapolation beyond elapsed year {elapsed_years}"
    )


def census_population(rows: list[dict[str, str]]) -> dict[int, float]:
    selected = [
        row
        for row in rows
        if row["SEX"] == "0" and row["ORIGIN"] == "0" and row["RACE"] == "0"
    ]
    result = {int(row["YEAR"]): float(row["TOTAL_POP"]) for row in selected}
    if not result:
        raise RuntimeError("Census population source has no national-total rows")
    return result


def census_births_by_sex(rows: list[dict[str, str]]) -> dict[tuple[int, int], float]:
    result = {
        (int(row["SEX"]), int(row["YEAR"])): float(row["BIRTHS"])
        for row in rows
        if row["RACE_HISP"] == "0" and row["SEX"] in ("1", "2")
    }
    if not result:
        raise RuntimeError("Census births source has no total-group sex-specific rows")
    return result


def native_survival_by_sex(
    rows: list[dict[str, str]],
) -> dict[tuple[int, int], list[float]]:
    result = {
        (int(row["SEX"]), int(row["YEAR"])): [
            float(row[f"SRAT_{age}"]) for age in range(101)
        ]
        for row in rows
        if row["NATIVITY"] == "1"
        and row["GROUP"] == "0"
        and row["SEX"] in ("1", "2")
    }
    if not result:
        raise RuntimeError("Census survival source has no native-born sex-specific rows")
    for key, values in result.items():
        if any(not 0.0 <= value <= 1.0 for value in values):
            raise RuntimeError(f"Invalid survival ratios at {key}")
    return result


def vintage_2025_total(rows: list[dict[str, str]]) -> float:
    matches = [row for row in rows if row["SEX"] == "0" and row["AGE"] == "999"]
    if len(matches) != 1:
        raise RuntimeError("Vintage 2025 source does not have one national total row")
    return float(matches[0]["POPESTIMATE2025"])


def build_scenario(
    scenario: str,
    policy_start_year: int,
    maximum_year: int,
    knots: list[dict[str, float | int]],
    projected_population: dict[int, float],
    projected_births: dict[tuple[int, int], float],
    survival: dict[tuple[int, int], list[float]],
) -> list[dict[str, Any]]:
    cohorts: dict[int, dict[int, float]] = {1: {}, 2: {}}
    cumulative_extra_births = 0.0
    rows: list[dict[str, Any]] = []
    for year in range(policy_start_year, maximum_year + 1):
        effect = interpolate_effect(knots, year - policy_start_year)
        annual_extra_births = 0.0
        for sex in (1, 2):
            schedule = survival[(sex, year)]
            cohorts[sex] = {
                age + 1: mass * schedule[age + 1]
                for age, mass in cohorts[sex].items()
                if age < 100
            }
            extra_births = projected_births[(sex, year)] * effect
            cohorts[sex][0] = extra_births * schedule[0]
            annual_extra_births += extra_births
        cumulative_extra_births += annual_extra_births
        extra_children = sum(
            mass
            for sex in (1, 2)
            for age, mass in cohorts[sex].items()
            if age <= 17
        )
        extra_adults = sum(
            mass
            for sex in (1, 2)
            for age, mass in cohorts[sex].items()
            if age >= 18
        )
        extra_total = extra_children + extra_adults
        baseline_total = projected_population[year]
        rows.append(
            {
                "scenario": scenario,
                "calendar_year": year,
                "model_aggregate_birth_effect_pct": 100.0 * effect,
                "annual_extra_births": annual_extra_births,
                "cumulative_extra_births": cumulative_extra_births,
                "extra_children_age_0_17": extra_children,
                "extra_adults_age_18_plus": extra_adults,
                "extra_resident_persons": extra_total,
                "official_baseline_resident_population": baseline_total,
                "official_baseline_plus_policy_wedge": baseline_total + extra_total,
                "extra_resident_persons_pct_baseline": 100.0
                * extra_total
                / baseline_total,
            }
        )
    return rows


def make_figure(rows: list[dict[str, Any]], path: Path) -> None:
    by_scenario = {
        scenario: [row for row in rows if row["scenario"] == scenario]
        for scenario in SCENARIOS
    }
    main = by_scenario["mid"]
    years = [int(row["calendar_year"]) for row in main]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))
    axes[0].plot(
        years,
        [float(row["model_aggregate_birth_effect_pct"]) for row in main],
        color="#b22222",
        linewidth=2.1,
    )
    axes[0].axhline(0.0, color="#333333", linewidth=0.8)
    axes[0].set(
        title="Model-implied aggregate birth wedge",
        xlabel="Calendar year",
        ylabel="Rebated 2% vs. rebated 1% (percent)",
    )
    colors = {"zero": "#999999", "low": "#7aa6c2", "mid": "#26547c", "hi": "#d95f02"}
    labels = {"zero": "Zero migration", "low": "Low", "mid": "Main", "hi": "High"}
    for scenario in SCENARIOS:
        scenario_rows = by_scenario[scenario]
        axes[1].plot(
            years,
            [float(row["extra_resident_persons"]) / 1e6 for row in scenario_rows],
            color=colors[scenario],
            linewidth=2.1 if scenario == "mid" else 1.35,
            label=labels[scenario],
        )
    axes[1].axhline(0.0, color="#333333", linewidth=0.8)
    axes[1].set(
        title="Additional resident persons",
        xlabel="Calendar year",
        ylabel="Millions",
    )
    axes[1].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    baseline_path = args.baseline_path.resolve()
    reform_path = args.reform_path.resolve()
    outdir = args.output_dir.resolve()
    for path in (baseline_path, reform_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    sources = outdir / "source_data"
    sources.mkdir()

    baseline_hash = validate_hash(baseline_path, args.baseline_sha256, "baseline")
    reform_hash = validate_hash(reform_path, args.reform_sha256, "reform")
    copied_baseline = sources / "rebated_tax1_transition_path.csv"
    copied_reform = sources / "rebated_tax2_transition_path.csv"
    copy_input(baseline_path, copied_baseline)
    copy_input(reform_path, copied_reform)

    downloaded: dict[str, Path] = {}
    source_metadata: dict[str, Any] = {
        "policy_paths": {
            "baseline": {"source": str(baseline_path), "sha256": baseline_hash},
            "reform": {"source": str(reform_path), "sha256": reform_hash},
        },
        "census_sources": {},
    }
    for label, url in CENSUS_SOURCES.items():
        destination = sources / f"{label}.csv"
        download(url, destination)
        downloaded[label] = destination
        source_metadata["census_sources"][label] = {
            "url": url,
            "sha256": sha256(destination),
        }

    baseline_rows = read_csv(copied_baseline)
    reform_rows = read_csv(copied_reform)
    knots = policy_birth_knots(baseline_rows, reform_rows)
    required_elapsed = int(args.maximum_year) - int(args.policy_start_year)
    if required_elapsed > int(knots[-1]["elapsed_years"]):
        raise RuntimeError(
            f"Maximum year requires {required_elapsed} elapsed years, but the policy "
            f"path supplies only {knots[-1]['elapsed_years']}"
        )
    survival = native_survival_by_sex(read_csv(downloaded["survival"]))
    all_rows: list[dict[str, Any]] = []
    for scenario in SCENARIOS:
        population = census_population(read_csv(downloaded[f"population_{scenario}"]))
        births = census_births_by_sex(read_csv(downloaded[f"births_{scenario}"]))
        all_rows.extend(
            build_scenario(
                scenario,
                int(args.policy_start_year),
                int(args.maximum_year),
                knots,
                population,
                births,
                survival,
            )
        )
    if any(float(row["extra_resident_persons"]) < 0.0 for row in all_rows):
        raise RuntimeError("The policy satellite produced a negative population wedge")
    if any(
        float(row["extra_resident_persons"])
        > float(row["cumulative_extra_births"]) + 1e-8
        for row in all_rows
    ):
        raise RuntimeError("Surviving policy cohorts exceed cumulative extra births")

    summary_rows = [
        row
        for row in all_rows
        if int(row["calendar_year"]) in SUMMARY_YEARS
    ]
    write_csv(outdir / "annual_persons_policy_wedge.csv", all_rows)
    write_csv(outdir / "summary.csv", summary_rows)
    write_csv(outdir / "model_birth_effect_knots.csv", knots)
    make_figure(all_rows, outdir / "persons_policy_wedge.png")

    actual_2025 = vintage_2025_total(read_csv(downloaded["vintage_2025_age_sex"]))
    mid_2025 = next(
        float(row["official_baseline_resident_population"])
        for row in all_rows
        if row["scenario"] == "mid" and int(row["calendar_year"]) == 2025
    )
    certification = {
        "status": "complete_exploratory_persons_satellite_not_promoted",
        "policy_start_year": int(args.policy_start_year),
        "maximum_year": int(args.maximum_year),
        "policy_path_status": "unpromoted preliminary perfect-foresight comparison",
        "policy_migration_response": "fixed to zero",
        "demographic_baselines": "Census 2023 main and alternative immigration scenarios",
        "survival": "Census native-born sex-specific projected survival ratios",
        "vintage_2025_resident_population": actual_2025,
        "census_2023_projection_for_2025": mid_2025,
        "vintage_minus_projection_percent": 100.0 * (actual_2025 / mid_2025 - 1.0),
        "accounting_gates": {
            "nonnegative_policy_population_wedge": True,
            "survivors_do_not_exceed_cumulative_births": True,
            "no_policy_effect_extrapolation": True,
        },
        "limitations": [
            "This is a satellite accounting exercise and does not feed people back into household decisions or housing markets.",
            "The Census 2023 projection is older than the Vintage 2025 population estimate.",
            "The model birth wedge is linearly interpolated between four-year equilibrium dates.",
            "The policy-induced migration response is fixed to zero.",
            "Adult persons are not converted into household heads with an age-specific headship hazard.",
            "The underlying perfect-foresight policy paths remain unpromoted preliminary diagnostics.",
        ],
        "sources": source_metadata,
    }
    (outdir / "certification.json").write_text(
        json.dumps(certification, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    main_summary = {
        int(row["calendar_year"]): row
        for row in summary_rows
        if row["scenario"] == "mid"
    }
    readme = f"""# Exploratory persons-based demographic satellite

Status: `complete_exploratory_persons_satellite_not_promoted`.

This diagnostic keeps the household equilibrium unchanged. It applies the
model-implied proportional change in aggregate births from the rebated 1% to
the rebated 2% property-tax path to the Census projected birth series, counts
the additional children immediately, and ages those cohorts with Census
native-born sex-specific survival ratios.

The policy starts in {int(args.policy_start_year)}. Net international migration
is common to the two policy cases, so the policy-induced migration response is
zero. The four Census immigration scenarios change baseline births and hence
the scale of the same proportional policy wedge.

## Census main-series readout

| Year | Extra children 0--17 | Extra adults 18+ | Extra resident persons | Percent of baseline |
|---|---:|---:|---:|---:|
"""
    for year in SUMMARY_YEARS:
        row = main_summary[year]
        readme += (
            f"| {year} | {float(row['extra_children_age_0_17']):,.0f} | "
            f"{float(row['extra_adults_age_18_plus']):,.0f} | "
            f"{float(row['extra_resident_persons']):,.0f} | "
            f"{float(row['extra_resident_persons_pct_baseline']):.3f}% |\n"
        )
    readme += f"""

The latest Census Vintage 2025 resident-population estimate is
{actual_2025:,.0f}, which is {100.0 * (actual_2025 / mid_2025 - 1.0):.2f}% above
the 2023 projection's value for 2025. The person wedge is therefore more
credible than treating the Census 2023 main series as a current point forecast.

## Interpretation and limits

This fixes the narrow accounting error that policy-born children were absent
from reported population until they entered as modeled household heads. It
does **not** yet create a demographic equilibrium: it does not map persons into
household heads, add spouses and non-head adults to the model state, or allow
policy-induced migration. The active four-year birth effect is linearly
interpolated to annual Census births. The underlying perfect-foresight policy
paths are preliminary and unpromoted.

All downloaded Census files, hashes, and copied policy paths are stored under
`source_data/`; `certification.json` records the exact provenance and gates.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
