#!/usr/bin/env python3
"""Diagnose the persons-to-household-head bridge for the E5F transition.

This script combines three objects without changing the model equilibrium:

* national ACS headship rates by sex and four-year age cell;
* Census resident-population projections by age and sex; and
* the policy-born person cohorts from the E5F demographic satellite.

It measures how many adult household heads a resident-person path implies and
compares that mapping with the active queue convention, which converts births
to household entrants at exactly age 20 using ``births / 2.1``.
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
import pandas as pd

from build_e5f_persons_demographic_satellite import (
    SCENARIOS,
    census_births_by_sex,
    interpolate_effect,
    native_survival_by_sex,
    policy_birth_knots,
    read_csv,
    sha256,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_EXTRACT = (
    ROOT / "code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
)
DEFAULT_SATELLITE = (
    ROOT / "output/model/e5f_persons_demographic_satellite_20260826b"
)
EXPECTED_EXTRACT_SHA256 = (
    "edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e"
)
ACS_YEARS = (2007, 2011, 2015, 2019, 2023)
AGE_LOWERS = tuple(range(18, 83, 4))
SUMMARY_YEARS = (2025, 2050, 2080, 2100)
SEXES = (1, 2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", type=Path, default=DEFAULT_EXTRACT)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--policy-start-year", type=int, default=2025)
    parser.add_argument("--maximum-year", type=int, default=2100)
    return parser.parse_args()


def load_acs_headship(path: Path) -> dict[int, dict[str, np.ndarray]]:
    """Load national ACS resident persons, household persons, and heads."""
    reader = pd.read_stata(
        path,
        iterator=True,
        convert_categoricals=False,
        convert_dates=False,
    )
    try:
        dtype = getattr(reader, "_dtype", None)
        data_location = getattr(reader, "data_location", None)
        raw_names = getattr(dtype, "names", None)
        if dtype is None or data_location is None or raw_names is None:
            raise RuntimeError("Unsupported pandas StataReader metadata")
        field = dict(zip(reader.varlist, raw_names))
        required = {
            "year",
            "gq",
            "pernum",
            "relate",
            "hhwt",
            "perwt",
            "sex",
            "age",
        }
        missing = sorted(required - set(field))
        if missing:
            raise ValueError(f"ACS extract is missing fields: {missing}")
        record_size = int(dtype.itemsize)
        stream = reader.path_or_buf

        def read_records(start: int, count: int) -> np.ndarray:
            stream.seek(int(data_location) + int(start) * record_size)
            raw = np.frombuffer(
                stream.read(int(count) * record_size), dtype=dtype, count=int(count)
            )
            if reader.byteorder != reader._native_byteorder:
                raw = raw.byteswap().newbyteorder()
            return raw

        def lower_year_bound(target_year: int) -> int:
            lower, upper = 0, int(reader.nobs)
            while lower < upper:
                midpoint = (lower + upper) // 2
                observed = int(read_records(midpoint, 1)[field["year"]][0])
                if observed < target_year:
                    lower = midpoint + 1
                else:
                    upper = midpoint
            return lower

        result: dict[int, dict[str, np.ndarray]] = {}
        for year in ACS_YEARS:
            first, last = lower_year_bound(year), lower_year_bound(year + 1)
            if first >= last:
                raise RuntimeError(f"No ACS records for {year}")
            arrays = {
                key: np.zeros((3, len(AGE_LOWERS)), dtype=float)
                for key in ("resident_people", "household_people", "heads")
            }
            for start in range(first, last, 250_000):
                raw = read_records(start, min(250_000, last - start))
                ages = raw[field["age"]]
                in_age = (ages >= 18) & (ages <= 85)
                bins = ((ages - 18) // 4).astype(int)
                valid_person = (
                    in_age
                    & np.isfinite(raw[field["perwt"]])
                    & (raw[field["perwt"]] > 0.0)
                )
                household_person = valid_person & np.isin(raw[field["gq"]], [1, 2])
                head = (
                    in_age
                    & np.isin(raw[field["gq"]], [1, 2])
                    & (raw[field["pernum"]] == 1)
                    & (raw[field["relate"]] == 1)
                    & np.isfinite(raw[field["hhwt"]])
                    & (raw[field["hhwt"]] > 0.0)
                )
                for sex in SEXES:
                    is_sex = raw[field["sex"]] == sex
                    for key, mask, weights in (
                        ("resident_people", valid_person & is_sex, raw[field["perwt"]]),
                        (
                            "household_people",
                            household_person & is_sex,
                            raw[field["perwt"]],
                        ),
                        ("heads", head & is_sex, raw[field["hhwt"]]),
                    ):
                        arrays[key][sex] += np.bincount(
                            bins[mask],
                            weights=weights[mask].astype(float),
                            minlength=len(AGE_LOWERS),
                        )
            for key in arrays:
                arrays[key][0] = arrays[key][1] + arrays[key][2]
            if np.any(arrays["resident_people"] <= 0.0) or np.any(
                arrays["heads"] <= 0.0
            ):
                raise RuntimeError(f"Nonpositive ACS headship input in {year}")
            arrays["headship_per_resident"] = (
                arrays["heads"] / arrays["resident_people"]
            )
            arrays["headship_per_household_person"] = (
                arrays["heads"] / arrays["household_people"]
            )
            arrays["group_quarters_share"] = 1.0 - (
                arrays["household_people"] / arrays["resident_people"]
            )
            result[year] = arrays
        return result
    finally:
        reader.close()


def acs_rows(acs: dict[int, dict[str, np.ndarray]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for year in ACS_YEARS:
        for sex in (0, 1, 2):
            for index, age in enumerate(AGE_LOWERS):
                rows.append(
                    {
                        "year": year,
                        "sex": sex,
                        "age_lower": age,
                        "age_upper": age + 3,
                        "resident_people": float(acs[year]["resident_people"][sex, index]),
                        "household_people": float(
                            acs[year]["household_people"][sex, index]
                        ),
                        "household_heads": float(acs[year]["heads"][sex, index]),
                        "headship_per_resident_person": float(
                            acs[year]["headship_per_resident"][sex, index]
                        ),
                        "headship_per_household_person": float(
                            acs[year]["headship_per_household_person"][sex, index]
                        ),
                        "group_quarters_share": float(
                            acs[year]["group_quarters_share"][sex, index]
                        ),
                    }
                )
    return rows


def census_population_by_sex_age(
    rows: list[dict[str, str]],
) -> dict[tuple[int, int], np.ndarray]:
    result: dict[tuple[int, int], np.ndarray] = {}
    for row in rows:
        if row["SEX"] not in ("1", "2") or row["ORIGIN"] != "0" or row["RACE"] != "0":
            continue
        result[(int(row["SEX"]), int(row["YEAR"]))] = np.asarray(
            [float(row[f"POP_{age}"]) for age in range(101)], dtype=float
        )
    if not result:
        raise RuntimeError("Census population source has no sex-specific national rows")
    return result


def vintage_population_by_sex_age(
    rows: list[dict[str, str]],
) -> dict[int, np.ndarray]:
    result: dict[int, np.ndarray] = {}
    for sex in SEXES:
        selected = {
            int(row["AGE"]): float(row["POPESTIMATE2025"])
            for row in rows
            if row["SEX"] == str(sex) and row["AGE"] != "999"
        }
        if set(selected) != set(range(101)):
            raise RuntimeError(f"Vintage 2025 age coverage is incomplete for sex {sex}")
        result[sex] = np.asarray([selected[age] for age in range(101)], dtype=float)
    return result


def implied_heads(
    population: dict[int, np.ndarray],
    headship: np.ndarray,
) -> tuple[float, float]:
    people = 0.0
    heads = 0.0
    for sex in SEXES:
        for age in range(18, 86):
            index = (age - 18) // 4
            mass = float(population[sex][age])
            people += mass
            heads += mass * float(headship[sex, index])
    return people, heads


def baseline_head_forecasts(
    source_dir: Path,
    acs: dict[int, dict[str, np.ndarray]],
    policy_start_year: int,
    maximum_year: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for scenario in SCENARIOS:
        population = census_population_by_sex_age(
            read_csv(source_dir / f"population_{scenario}.csv")
        )
        for headship_year in ACS_YEARS:
            headship = acs[headship_year]["headship_per_resident"]
            for year in range(policy_start_year, maximum_year + 1):
                people, heads = implied_heads(
                    {sex: population[(sex, year)] for sex in SEXES}, headship
                )
                rows.append(
                    {
                        "immigration_scenario": scenario,
                        "headship_vintage": headship_year,
                        "calendar_year": year,
                        "resident_persons_age_18_85": people,
                        "implied_household_heads_age_18_85": heads,
                        "aggregate_headship_rate": heads / people,
                    }
                )
    return rows


def policy_headship_forecasts(
    source_dir: Path,
    acs: dict[int, dict[str, np.ndarray]],
    policy_start_year: int,
    maximum_year: int,
) -> list[dict[str, Any]]:
    baseline_rows = read_csv(source_dir / "rebated_tax1_transition_path.csv")
    reform_rows = read_csv(source_dir / "rebated_tax2_transition_path.csv")
    knots = policy_birth_knots(baseline_rows, reform_rows)
    survival = native_survival_by_sex(read_csv(source_dir / "survival.csv"))
    headship = acs[2023]["headship_per_resident"]
    output: list[dict[str, Any]] = []
    for scenario in SCENARIOS:
        births = census_births_by_sex(read_csv(source_dir / f"births_{scenario}.csv"))
        cohorts: dict[int, dict[int, float]] = {1: {}, 2: {}}
        birth_history: list[tuple[int, float]] = []
        for year in range(policy_start_year, maximum_year + 1):
            effect = interpolate_effect(knots, year - policy_start_year)
            annual_extra_births = 0.0
            for sex in SEXES:
                schedule = survival[(sex, year)]
                cohorts[sex] = {
                    age + 1: mass * schedule[age + 1]
                    for age, mass in cohorts[sex].items()
                    if age < 100
                }
                extra_births = float(births[(sex, year)]) * effect
                cohorts[sex][0] = extra_births * schedule[0]
                annual_extra_births += extra_births
            birth_history.append((year, annual_extra_births))
            extra_children = sum(
                mass
                for sex in SEXES
                for age, mass in cohorts[sex].items()
                if age <= 17
            )
            extra_adult_persons = sum(
                mass
                for sex in SEXES
                for age, mass in cohorts[sex].items()
                if 18 <= age <= 85
            )
            extra_heads = sum(
                mass * float(headship[sex, (age - 18) // 4])
                for sex in SEXES
                for age, mass in cohorts[sex].items()
                if 18 <= age <= 85
            )
            queue_proxy_heads = sum(
                births_mass / 2.1
                for birth_year, births_mass in birth_history
                if birth_year + 20 <= year <= birth_year + 84
            )
            output.append(
                {
                    "immigration_scenario": scenario,
                    "calendar_year": year,
                    "model_aggregate_birth_effect_pct": 100.0 * effect,
                    "annual_extra_births": annual_extra_births,
                    "extra_children_age_0_17": extra_children,
                    "extra_adult_persons_age_18_85": extra_adult_persons,
                    "extra_heads_2023_acs_headship": extra_heads,
                    "queue_proxy_extra_heads_births_divided_by_2_1_at_age_20": (
                        queue_proxy_heads
                    ),
                    "queue_proxy_minus_headship_mapped_heads": queue_proxy_heads
                    - extra_heads,
                }
            )
    return output


def make_figure(
    acs: dict[int, dict[str, np.ndarray]],
    baseline_rows: list[dict[str, Any]],
    policy_rows: list[dict[str, Any]],
    path: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.5))
    ages = np.asarray(AGE_LOWERS, dtype=float) + 1.5
    for year, color in ((2007, "#999999"), (2015, "#7aa6c2"), (2023, "#26547c")):
        axes[0].plot(
            ages,
            acs[year]["headship_per_resident"][0],
            color=color,
            linewidth=2.0 if year == 2023 else 1.35,
            label=str(year),
        )
    axes[0].axhline(1.0 / 2.1, color="#b22222", linestyle="--", linewidth=1.2, label="1 / 2.1")
    axes[0].set(
        title="Headship by age",
        xlabel="Age-cell midpoint",
        ylabel="Household heads / resident persons",
    )
    axes[0].legend(frameon=False, fontsize=8)

    main = [
        row
        for row in baseline_rows
        if row["immigration_scenario"] == "mid"
    ]
    for year, color in ((2007, "#999999"), (2015, "#7aa6c2"), (2023, "#26547c")):
        selected = [row for row in main if int(row["headship_vintage"]) == year]
        axes[1].plot(
            [int(row["calendar_year"]) for row in selected],
            [float(row["implied_household_heads_age_18_85"]) / 1e6 for row in selected],
            color=color,
            linewidth=2.0 if year == 2023 else 1.35,
            label=f"{year} rates",
        )
    axes[1].set(
        title="Heads under Census main population",
        xlabel="Calendar year",
        ylabel="Millions, ages 18--85",
    )
    axes[1].legend(frameon=False, fontsize=8)

    policy = [row for row in policy_rows if row["immigration_scenario"] == "mid"]
    years = [int(row["calendar_year"]) for row in policy]
    axes[2].plot(
        years,
        [float(row["extra_adult_persons_age_18_85"]) / 1e6 for row in policy],
        color="#7aa6c2",
        linewidth=1.8,
        label="Extra adult persons",
    )
    axes[2].plot(
        years,
        [float(row["extra_heads_2023_acs_headship"]) / 1e6 for row in policy],
        color="#26547c",
        linewidth=2.1,
        label="Heads via ACS headship",
    )
    axes[2].plot(
        years,
        [
            float(row["queue_proxy_extra_heads_births_divided_by_2_1_at_age_20"])
            / 1e6
            for row in policy
        ],
        color="#b22222",
        linewidth=1.4,
        linestyle="--",
        label="Current queue proxy",
    )
    axes[2].set(
        title="Policy-born adults and implied heads",
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
    extract = args.extract.resolve()
    satellite = args.satellite_dir.resolve()
    source_dir = satellite / "source_data"
    outdir = args.output_dir.resolve()
    for path in (extract, source_dir):
        if not path.exists():
            raise FileNotFoundError(path)
    if sha256(extract) != EXPECTED_EXTRACT_SHA256:
        raise RuntimeError("The canonical ACS extract fingerprint differs")
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    acs = load_acs_headship(extract)
    headship_rows = acs_rows(acs)
    baseline_rows = baseline_head_forecasts(
        source_dir,
        acs,
        int(args.policy_start_year),
        int(args.maximum_year),
    )
    policy_rows = policy_headship_forecasts(
        source_dir,
        acs,
        int(args.policy_start_year),
        int(args.maximum_year),
    )
    write_csv(outdir / "acs_headship_profiles.csv", headship_rows)
    write_csv(outdir / "census_implied_head_forecasts.csv", baseline_rows)
    write_csv(outdir / "policy_person_to_head_wedge.csv", policy_rows)
    make_figure(acs, baseline_rows, policy_rows, outdir / "headship_bridge.png")

    vintage = vintage_population_by_sex_age(
        read_csv(source_dir / "vintage_2025_age_sex.csv")
    )
    vintage_people, vintage_heads = implied_heads(
        vintage, acs[2023]["headship_per_resident"]
    )
    main_summary = {
        int(row["calendar_year"]): row
        for row in policy_rows
        if row["immigration_scenario"] == "mid"
        and int(row["calendar_year"]) in SUMMARY_YEARS
    }
    baseline_summary = {
        int(row["calendar_year"]): row
        for row in baseline_rows
        if row["immigration_scenario"] == "mid"
        and int(row["headship_vintage"]) == 2023
        and int(row["calendar_year"]) in SUMMARY_YEARS
    }

    rate_18_21 = float(acs[2023]["headship_per_resident"][0, 0])
    rate_22_25 = float(acs[2023]["headship_per_resident"][0, 1])
    queue_rate = 1.0 / 2.1
    gates = {
        "headship_rates_between_zero_and_one": all(
            0.0 <= float(row["headship_per_resident_person"]) <= 1.0
            for row in headship_rows
        ),
        "group_quarters_shares_between_zero_and_one": all(
            0.0 <= float(row["group_quarters_share"]) <= 1.0
            for row in headship_rows
        ),
        "policy_head_counts_do_not_exceed_adult_persons": all(
            float(row["extra_heads_2023_acs_headship"])
            <= float(row["extra_adult_persons_age_18_85"]) + 1e-8
            for row in policy_rows
        ),
    }
    if not all(gates.values()):
        raise RuntimeError(f"Headship accounting gates failed: {gates}")

    source_hashes = {
        path.name: sha256(path)
        for path in sorted(source_dir.iterdir())
        if path.is_file()
    }
    certification = {
        "status": "complete_exploratory_headship_bridge_not_promoted",
        "acs_extract": str(extract),
        "acs_extract_sha256": EXPECTED_EXTRACT_SHA256,
        "satellite_source_hashes": source_hashes,
        "policy_start_year": int(args.policy_start_year),
        "maximum_year": int(args.maximum_year),
        "vintage_2025_resident_persons_age_18_85": vintage_people,
        "vintage_2025_implied_heads_age_18_85_at_2023_rates": vintage_heads,
        "queue_birth_to_head_rate": queue_rate,
        "acs_2023_headship_rate_18_21": rate_18_21,
        "acs_2023_headship_rate_22_25": rate_22_25,
        "queue_over_acs_18_21_headship": queue_rate / rate_18_21,
        "accounting_gates": gates,
        "limitations": [
            "Headship rates are accounting mappings, not estimated structural household-formation hazards.",
            "Future headship is shown under frozen historical ACS vintages rather than forecast behavior.",
            "The policy-induced migration and headship responses are fixed to zero.",
            "The underlying perfect-foresight policy paths remain unpromoted preliminary diagnostics.",
        ],
    }
    (outdir / "certification.json").write_text(
        json.dumps(certification, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    readme = f"""# Exploratory persons-to-household-head bridge

Status: `complete_exploratory_headship_bridge_not_promoted`.

## Main diagnosis

The active renewal queue converts each birth cohort into household heads at
age 20 at rate `1 / 2.1 = {queue_rate:.3f}`. In the canonical 2023 ACS extract,
only {rate_18_21:.3f} of resident persons ages 18--21 and {rate_22_25:.3f} of
persons ages 22--25 are household heads. The queue conversion is therefore
{queue_rate / rate_18_21:.1f} times the observed 18--21 headship rate. It is a
replacement normalization, not a defensible household-formation schedule.

Applying frozen 2023 age- and sex-specific headship rates to the latest Census
Vintage 2025 age distribution implies {vintage_heads / 1e6:.2f} million
household heads ages 18--85 out of {vintage_people / 1e6:.2f} million resident
persons in that age range.

## Census-main demographic and policy readout

| Year | Baseline heads 18--85 | Extra adult persons | Extra heads via ACS rates | Current queue proxy heads |
|---|---:|---:|---:|---:|
"""
    for year in SUMMARY_YEARS:
        baseline = baseline_summary[year]
        policy = main_summary[year]
        readme += (
            f"| {year} | {float(baseline['implied_household_heads_age_18_85']) / 1e6:.3f}m | "
            f"{float(policy['extra_adult_persons_age_18_85']):,.0f} | "
            f"{float(policy['extra_heads_2023_acs_headship']):,.0f} | "
            f"{float(policy['queue_proxy_extra_heads_births_divided_by_2_1_at_age_20']):,.0f} |\n"
        )
    readme += """

## Interpretation

The persons satellite fixes population reporting immediately. This bridge
shows how policy-born adult persons could subsequently enter the household
model. It should replace the single age-20 conversion with an age-specific
headship or household-formation schedule before demographic feedback is fed
back into housing equilibrium. Immigration must remain a separate age-specific
person flow; the historical head residual cannot identify it.

The output is an accounting diagnostic. It does not alter the active solver or
promote the preliminary perfect-foresight policy comparison.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
