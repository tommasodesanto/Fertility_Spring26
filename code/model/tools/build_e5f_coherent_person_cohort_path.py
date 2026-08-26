#!/usr/bin/env python3
"""Build an isolated coherent person-to-household-head demographic path.

The driver does not call the household solver.  It replaces the diagnostic
``births / 2.1`` head-renewal queue with an annual cohort-component layer and
writes the exact four-year household-head inputs that a later equilibrium
integration can consume.

The baseline is the Census 2023 main projection, cohort-anchored to the latest
Vintage 2025 age-sex population.  Births, survival, and the inferred residual
net-migration path are kept as separate ledgers.  Fixed 2023 ACS age-sex
headship rates then generate the minimum gross nonhead-to-head formation or
head-dissolution flow required in each cell.  The policy wedge changes births
only, applies native-born survival, and holds migration fixed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demographic_transition.person_cohort_law import (
    CohortState,
    infer_net_migration_path,
    simulate_path,
)
from build_e5f_persons_demographic_satellite import (
    interpolate_effect,
    policy_birth_knots,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE_DIR = (
    ROOT
    / "output/model/e5f_persons_demographic_satellite_20260826b/source_data"
)
DEFAULT_HEADSHIP_DIR = (
    ROOT / "output/model/e5f_headship_demographic_bridge_20260826a"
)
SEXES = (1, 2)
MAX_AGE = 100
MODEL_AGE_LOWERS = tuple(range(18, 83, 4))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--headship-dir", type=Path, default=DEFAULT_HEADSHIP_DIR)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=2025)
    parser.add_argument("--maximum-year", type=int, default=2100)
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
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def population_path(rows: list[dict[str, str]]) -> dict[int, np.ndarray]:
    selected: dict[int, np.ndarray] = {}
    for row in rows:
        if row["SEX"] not in ("1", "2") or row["ORIGIN"] != "0" or row["RACE"] != "0":
            continue
        year = int(row["YEAR"])
        selected.setdefault(year, np.zeros((2, MAX_AGE + 1), dtype=float))
        selected[year][int(row["SEX"]) - 1] = np.asarray(
            [float(row[f"POP_{age}"]) for age in range(MAX_AGE + 1)],
            dtype=float,
        )
    if not selected:
        raise RuntimeError("Population file has no national sex-specific rows")
    return selected


def vintage_2025(rows: list[dict[str, str]]) -> np.ndarray:
    result = np.zeros((2, MAX_AGE + 1), dtype=float)
    for sex in SEXES:
        by_age = {
            int(row["AGE"]): float(row["POPESTIMATE2025"])
            for row in rows
            if row["SEX"] == str(sex) and row["AGE"] != "999"
        }
        if set(by_age) != set(range(MAX_AGE + 1)):
            raise RuntimeError(f"Vintage 2025 age coverage is incomplete for sex {sex}")
        result[sex - 1] = [by_age[age] for age in range(MAX_AGE + 1)]
    return result


def births_path(rows: list[dict[str, str]]) -> dict[int, np.ndarray]:
    result: dict[int, np.ndarray] = {}
    for row in rows:
        if row["RACE_HISP"] != "0" or row["SEX"] not in ("1", "2"):
            continue
        year = int(row["YEAR"])
        result.setdefault(year, np.zeros(2, dtype=float))
        result[year][int(row["SEX"]) - 1] = float(row["BIRTHS"])
    if not result:
        raise RuntimeError("Birth file has no national sex-specific rows")
    return result


def survival_path(
    rows: list[dict[str, str]], *, nativity: int
) -> dict[int, np.ndarray]:
    result: dict[int, np.ndarray] = {}
    for row in rows:
        if (
            row["NATIVITY"] != str(nativity)
            or row["GROUP"] != "0"
            or row["SEX"] not in ("1", "2")
        ):
            continue
        year = int(row["YEAR"])
        result.setdefault(year, np.zeros((2, MAX_AGE + 1), dtype=float))
        result[year][int(row["SEX"]) - 1] = np.asarray(
            [float(row[f"SRAT_{age}"]) for age in range(MAX_AGE + 1)],
            dtype=float,
        )
    if not result:
        raise RuntimeError(f"Survival file has no nativity={nativity} rows")
    return result


def fixed_2023_headship(rows: list[dict[str, str]]) -> np.ndarray:
    result = np.zeros((2, MAX_AGE + 1), dtype=float)
    selected = [
        row
        for row in rows
        if int(row["year"]) == 2023 and int(row["sex"]) in SEXES
    ]
    if len(selected) != 2 * len(MODEL_AGE_LOWERS):
        raise RuntimeError("Headship file does not contain the complete 2023 profile")
    for row in selected:
        sex = int(row["sex"]) - 1
        lower = int(row["age_lower"])
        upper = int(row["age_upper"])
        result[sex, lower : upper + 1] = float(
            row["headship_per_resident_person"]
        )
    return result


def cohort_anchor_projection(
    projection: dict[int, np.ndarray],
    vintage: np.ndarray,
    *,
    start_year: int,
    maximum_year: int,
) -> dict[int, np.ndarray]:
    """Carry the Vintage-2025 cohort correction through the projection."""

    if 2025 not in projection:
        raise RuntimeError("Projection does not contain 2025")
    projected_2025 = projection[2025]
    if np.any(projected_2025 <= 0.0):
        raise RuntimeError("Projection has nonpositive 2025 age-sex cells")
    ratio = vintage / projected_2025
    result: dict[int, np.ndarray] = {}
    for year in range(start_year, maximum_year + 1):
        if year not in projection:
            raise RuntimeError(f"Projection does not contain {year}")
        anchored = projection[year].copy()
        elapsed = year - 2025
        for age in range(MAX_AGE):
            age_in_2025 = age - elapsed
            if 0 <= age_in_2025 < MAX_AGE:
                anchored[:, age] *= ratio[:, age_in_2025]
        anchored[:, MAX_AGE] *= ratio[:, MAX_AGE]
        result[year] = anchored
    if not np.allclose(result[2025], vintage, atol=1e-8, rtol=0.0):
        raise RuntimeError("Cohort anchor does not reproduce Vintage 2025")
    return result


def aggregate_model_heads(state: CohortState) -> float:
    return float(np.sum(state.heads[:, 18:86]))


def annual_rows(
    baseline_states: dict[int, CohortState],
    baseline_ledgers: dict[int, Any],
    policy_states: dict[int, CohortState],
    policy_ledgers: dict[int, Any],
    baseline_births: dict[int, np.ndarray],
    policy_births: dict[int, np.ndarray],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for year in sorted(baseline_states):
        baseline = baseline_states[year]
        wedge = policy_states[year]
        ledger = baseline_ledgers.get(year)
        wedge_ledger = policy_ledgers.get(year)
        rows.append(
            {
                "calendar_year": year,
                "baseline_resident_persons": float(np.sum(baseline.persons)),
                "baseline_heads_age_18_85": aggregate_model_heads(baseline),
                "policy_extra_resident_persons": float(np.sum(wedge.persons)),
                "policy_extra_heads_age_18_85": aggregate_model_heads(wedge),
                "reform_resident_persons": float(
                    np.sum(baseline.persons) + np.sum(wedge.persons)
                ),
                "reform_heads_age_18_85": aggregate_model_heads(baseline)
                + aggregate_model_heads(wedge),
                "baseline_births": float(np.sum(baseline_births.get(year, np.zeros(2)))),
                "policy_extra_births": float(np.sum(policy_births.get(year, np.zeros(2)))),
                "net_migration": 0.0 if ledger is None else float(np.sum(ledger.net_migration)),
                "new_heads_from_nonheads": 0.0
                if ledger is None
                else float(np.sum(ledger.new_heads_from_nonheads[:, 18:86])),
                "head_dissolutions": 0.0
                if ledger is None
                else float(np.sum(ledger.head_dissolutions[:, 18:86])),
                "net_migrant_heads": 0.0
                if ledger is None
                else float(np.sum(ledger.net_migrant_heads[:, 18:86])),
                "existing_head_deaths": 0.0
                if ledger is None
                else float(np.sum(ledger.existing_head_deaths)),
                "policy_new_heads_from_nonheads": 0.0
                if wedge_ledger is None
                else float(np.sum(wedge_ledger.new_heads_from_nonheads[:, 18:86])),
                "baseline_person_identity_residual": 0.0
                if ledger is None
                else float(ledger.person_identity_residual),
                "baseline_head_identity_residual": 0.0
                if ledger is None
                else float(ledger.head_identity_residual),
                "policy_person_identity_residual": 0.0
                if wedge_ledger is None
                else float(wedge_ledger.person_identity_residual),
                "policy_head_identity_residual": 0.0
                if wedge_ledger is None
                else float(wedge_ledger.head_identity_residual),
            }
        )
    return rows


def age_cell_rows(
    baseline_states: dict[int, CohortState],
    policy_states: dict[int, CohortState],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    years = [
        year
        for year in sorted(baseline_states)
        if (year - min(baseline_states)) % 4 == 0
    ]
    for year in years:
        for sex in (0, 1, 2):
            for lower in MODEL_AGE_LOWERS:
                upper = lower + 3
                if sex == 0:
                    baseline_persons = float(
                        np.sum(baseline_states[year].persons[:, lower : upper + 1])
                    )
                    baseline_heads = float(
                        np.sum(baseline_states[year].heads[:, lower : upper + 1])
                    )
                    policy_persons = float(
                        np.sum(policy_states[year].persons[:, lower : upper + 1])
                    )
                    policy_heads = float(
                        np.sum(policy_states[year].heads[:, lower : upper + 1])
                    )
                else:
                    baseline_persons = float(
                        np.sum(baseline_states[year].persons[sex - 1, lower : upper + 1])
                    )
                    baseline_heads = float(
                        np.sum(baseline_states[year].heads[sex - 1, lower : upper + 1])
                    )
                    policy_persons = float(
                        np.sum(policy_states[year].persons[sex - 1, lower : upper + 1])
                    )
                    policy_heads = float(
                        np.sum(policy_states[year].heads[sex - 1, lower : upper + 1])
                    )
                rows.append(
                    {
                        "calendar_year": year,
                        "sex": sex,
                        "age_lower": lower,
                        "age_upper": upper,
                        "baseline_resident_persons": baseline_persons,
                        "baseline_household_heads": baseline_heads,
                        "policy_extra_resident_persons": policy_persons,
                        "policy_extra_household_heads": policy_heads,
                        "reform_resident_persons": baseline_persons + policy_persons,
                        "reform_household_heads": baseline_heads + policy_heads,
                    }
                )
    return rows


def four_year_flow_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_year = {int(row["calendar_year"]): row for row in rows}
    first = min(by_year)
    output: list[dict[str, Any]] = []
    for year in range(first, max(by_year) + 1, 4):
        window = [by_year[value] for value in range(max(first + 1, year - 3), year + 1)]
        row = by_year[year]
        output.append(
            {
                "calendar_year": year,
                "baseline_resident_persons": row["baseline_resident_persons"],
                "baseline_heads_age_18_85": row["baseline_heads_age_18_85"],
                "policy_extra_resident_persons": row["policy_extra_resident_persons"],
                "policy_extra_heads_age_18_85": row["policy_extra_heads_age_18_85"],
                "reform_resident_persons": row["reform_resident_persons"],
                "reform_heads_age_18_85": row["reform_heads_age_18_85"],
                "births_over_previous_four_years": sum(
                    float(value["baseline_births"]) for value in window
                ),
                "policy_extra_births_over_previous_four_years": sum(
                    float(value["policy_extra_births"]) for value in window
                ),
                "new_heads_from_nonheads_over_previous_four_years": sum(
                    float(value["new_heads_from_nonheads"]) for value in window
                ),
                "head_dissolutions_over_previous_four_years": sum(
                    float(value["head_dissolutions"]) for value in window
                ),
                "net_migration_over_previous_four_years": sum(
                    float(value["net_migration"]) for value in window
                ),
                "net_migrant_heads_over_previous_four_years": sum(
                    float(value["net_migrant_heads"]) for value in window
                ),
            }
        )
    return output


def make_figure(
    rows: list[dict[str, Any]],
    age_rows: list[dict[str, Any]],
    path: Path,
) -> None:
    years = np.asarray([int(row["calendar_year"]) for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.0))
    axes[0, 0].plot(
        years,
        [float(row["baseline_resident_persons"]) / 1e6 for row in rows],
        color="#26547c",
        linewidth=2.0,
        label="Baseline persons",
    )
    axes[0, 0].plot(
        years,
        [float(row["baseline_heads_age_18_85"]) / 1e6 for row in rows],
        color="#7aa6c2",
        linewidth=2.0,
        label="Baseline heads 18–85",
    )
    axes[0, 0].set(title="Cohort-anchored baseline", ylabel="Millions")
    axes[0, 0].legend(frameon=False, fontsize=8)

    axes[0, 1].plot(
        years,
        [float(row["policy_extra_resident_persons"]) / 1e6 for row in rows],
        color="#b22222",
        linewidth=2.0,
        label="Extra persons",
    )
    axes[0, 1].plot(
        years,
        [float(row["policy_extra_heads_age_18_85"]) / 1e6 for row in rows],
        color="#d17c4b",
        linewidth=2.0,
        label="Extra heads 18–85",
    )
    axes[0, 1].set(title="Birth-only policy wedge", ylabel="Millions")
    axes[0, 1].legend(frameon=False, fontsize=8)

    axes[1, 0].plot(
        years,
        [float(row["new_heads_from_nonheads"]) / 1e6 for row in rows],
        color="#2a9d8f",
        linewidth=1.8,
        label="Nonhead → head",
    )
    axes[1, 0].plot(
        years,
        [float(row["head_dissolutions"]) / 1e6 for row in rows],
        color="#e76f51",
        linewidth=1.8,
        label="Head dissolution",
    )
    axes[1, 0].plot(
        years,
        [float(row["net_migrant_heads"]) / 1e6 for row in rows],
        color="#777777",
        linewidth=1.3,
        linestyle="--",
        label="Net migrant heads",
    )
    axes[1, 0].set(title="Household-head flow ledger", ylabel="Millions/year")
    axes[1, 0].legend(frameon=False, fontsize=8)

    for year, color in ((2025, "#999999"), (2049, "#7aa6c2"), (2081, "#26547c")):
        selected = [
            row
            for row in age_rows
            if int(row["calendar_year"]) == year and int(row["sex"]) == 0
        ]
        if not selected:
            continue
        axes[1, 1].plot(
            [0.5 * (int(row["age_lower"]) + int(row["age_upper"])) for row in selected],
            [float(row["baseline_household_heads"]) / 1e6 for row in selected],
            color=color,
            linewidth=2.0 if year == 2080 else 1.5,
            label=str(year),
        )
    axes[1, 1].set(
        title="Baseline household heads by age cell",
        xlabel="Age-cell midpoint",
        ylabel="Millions",
    )
    axes[1, 1].legend(frameon=False, fontsize=8)
    for axis in axes.flat:
        axis.grid(alpha=0.18)
        axis.set_xlabel(axis.get_xlabel() or "Calendar year")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    source_dir = args.source_dir.resolve()
    headship_dir = args.headship_dir.resolve()
    outdir = args.output_dir.resolve()
    if not source_dir.is_dir() or not headship_dir.is_dir():
        raise FileNotFoundError("Source or headship directory is missing")
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    source_paths = {
        "population_mid": source_dir / "population_mid.csv",
        "births_mid": source_dir / "births_mid.csv",
        "survival": source_dir / "survival.csv",
        "vintage_2025": source_dir / "vintage_2025_age_sex.csv",
        "rebated_tax1": source_dir / "rebated_tax1_transition_path.csv",
        "rebated_tax2": source_dir / "rebated_tax2_transition_path.csv",
        "acs_headship": headship_dir / "acs_headship_profiles.csv",
    }
    for path in source_paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    start_year = int(args.start_year)
    maximum_year = int(args.maximum_year)
    raw_population = population_path(read_csv(source_paths["population_mid"]))
    vintage = vintage_2025(read_csv(source_paths["vintage_2025"]))
    target = cohort_anchor_projection(
        raw_population,
        vintage,
        start_year=start_year,
        maximum_year=maximum_year,
    )
    births = births_path(read_csv(source_paths["births_mid"]))
    survival_total = survival_path(read_csv(source_paths["survival"]), nativity=0)
    survival_native = survival_path(read_csv(source_paths["survival"]), nativity=1)
    headship = fixed_2023_headship(read_csv(source_paths["acs_headship"]))

    migration = infer_net_migration_path(
        target,
        births,
        survival_total,
        topcoded_last_age=True,
    )
    baseline_initial = CohortState(
        year=start_year,
        persons=target[start_year],
        heads=target[start_year] * headship,
    )
    baseline_states, baseline_ledgers = simulate_path(
        baseline_initial,
        births=births,
        survival=survival_total,
        net_migration=migration,
        headship_rates=headship,
        final_year=maximum_year,
    )

    knots = policy_birth_knots(
        read_csv(source_paths["rebated_tax1"]),
        read_csv(source_paths["rebated_tax2"]),
    )
    birth_wedge = {
        year: births[year] * interpolate_effect(knots, year - start_year)
        for year in range(start_year + 1, maximum_year + 1)
    }
    zero_migration = {
        year: np.zeros_like(target[start_year])
        for year in range(start_year + 1, maximum_year + 1)
    }
    zero = np.zeros_like(target[start_year])
    policy_initial = CohortState(start_year, zero, zero)
    policy_states, policy_ledgers = simulate_path(
        policy_initial,
        births=birth_wedge,
        survival=survival_native,
        net_migration=zero_migration,
        headship_rates=headship,
        final_year=maximum_year,
    )

    rows = annual_rows(
        baseline_states,
        baseline_ledgers,
        policy_states,
        policy_ledgers,
        births,
        birth_wedge,
    )
    cells = age_cell_rows(baseline_states, policy_states)
    bridge = four_year_flow_rows(rows)
    migration_rows = [
        {
            "calendar_year": year,
            "sex": sex + 1,
            "age": age,
            "net_migration_residual": float(migration[year][sex, age]),
        }
        for year in sorted(migration)
        for sex in range(2)
        for age in range(MAX_AGE + 1)
    ]
    write_csv(outdir / "annual_person_head_path.csv", rows)
    write_csv(outdir / "household_head_age_cells_4year.csv", cells)
    write_csv(outdir / "model_head_bridge_4year.csv", bridge)
    write_csv(outdir / "net_migration_residual_by_age_sex.csv", migration_rows)
    make_figure(rows, cells, outdir / "coherent_person_cohort_path.png")

    maximum_target_gap = max(
        float(np.max(np.abs(baseline_states[year].persons - target[year])))
        for year in baseline_states
    )
    all_ledgers = list(baseline_ledgers.values()) + list(policy_ledgers.values())
    max_person_residual = max(abs(row.person_identity_residual) for row in all_ledgers)
    max_head_residual = max(abs(row.head_identity_residual) for row in all_ledgers)
    max_mapping_residual = max(row.headship_mapping_residual for row in all_ledgers)
    vintage_total = float(np.sum(vintage))
    initial_total = float(np.sum(baseline_states[start_year].persons))
    gates = {
        "vintage_2025_age_sex_reproduced": bool(
            np.allclose(baseline_states[2025].persons, vintage, atol=1e-8, rtol=0.0)
        ),
        "baseline_target_path_reproduced": maximum_target_gap <= 1e-7,
        "person_flow_identity": max_person_residual <= 1e-6,
        "head_flow_identity": max_head_residual <= 1e-6,
        "headship_mapping_identity": max_mapping_residual <= 1e-7,
        "policy_migration_response_zero": all(
            np.count_nonzero(row.net_migration) == 0 for row in policy_ledgers.values()
        ),
        "nonnegative_baseline_and_policy_mass": all(
            np.min(state.persons) >= 0.0 and np.min(state.heads) >= 0.0
            for state in list(baseline_states.values()) + list(policy_states.values())
        ),
    }
    if not all(gates.values()):
        raise RuntimeError(f"Coherent demographic-law gates failed: {gates}")

    selected = {int(row["calendar_year"]): row for row in rows}
    certification = {
        "status": "complete_isolated_coherent_demographic_law_not_integrated",
        "source_hashes": {name: sha256(path) for name, path in source_paths.items()},
        "start_year": start_year,
        "maximum_year": maximum_year,
        "baseline": "Census 2023 main projection cohort-anchored to Vintage 2025 age-sex population",
        "headship": "fixed 2023 ACS sex-by-four-year-age-cell rates; minimum one-way formation/dissolution accounting flow",
        "policy": "rebated-2%-minus-rebated-1% model birth wedge with native-born survival and zero migration response",
        "vintage_2025_total": vintage_total,
        "initial_total": initial_total,
        "maximum_baseline_target_gap": maximum_target_gap,
        "maximum_person_identity_residual": max_person_residual,
        "maximum_head_identity_residual": max_head_residual,
        "maximum_headship_mapping_residual": max_mapping_residual,
        "accounting_gates": gates,
        "selected_policy_wedges": {
            str(year): {
                "extra_resident_persons": selected[year]["policy_extra_resident_persons"],
                "extra_heads_age_18_85": selected[year]["policy_extra_heads_age_18_85"],
            }
            for year in (2050, 2080, 2100)
        },
        "limitations": [
            "This path is isolated from the household Bellman problem and housing-market solver.",
            "Headship rates are fixed at 2023 ACS values and the minimum gross flow is an accounting closure, not an estimated formation hazard.",
            "The baseline migration residual includes every difference between the anchored Census target and the maintained total-population survival/birth law.",
            "The policy starts from the common 2025 state, so the first changed birth flow is 2026.",
            "The Census projection ends in 2100; terminal demographic primitives are not yet extrapolated or solved as a stationary fixed point.",
        ],
    }
    (outdir / "certification.json").write_text(
        json.dumps(certification, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    readme = f"""# Coherent person-to-household-head path

Status: `complete_isolated_coherent_demographic_law_not_integrated`.

This packet removes the diagnostic `births / 2.1 at age 20` renewal rule. It
starts from the exact Vintage 2025 age-sex population ({vintage_total:,.0f}
resident persons), advances births and survival by single year of age, records
the age-sex net-migration residual separately, and converts persons into
household heads through explicit nonhead-to-head formation and head-dissolution
flows.

The accounting closes exactly: the largest person-flow residual is
`{max_person_residual:.3e}`, the largest head-flow residual is
`{max_head_residual:.3e}`, and the baseline reproduces every anchored Census
age-sex cell with maximum absolute error `{maximum_target_gap:.3e}` persons.

| Year | Extra resident persons | Extra household heads, ages 18–85 |
|---|---:|---:|
| 2050 | {float(selected[2050]['policy_extra_resident_persons']):,.0f} | {float(selected[2050]['policy_extra_heads_age_18_85']):,.0f} |
| 2080 | {float(selected[2080]['policy_extra_resident_persons']):,.0f} | {float(selected[2080]['policy_extra_heads_age_18_85']):,.0f} |
| 2100 | {float(selected[2100]['policy_extra_resident_persons']):,.0f} | {float(selected[2100]['policy_extra_heads_age_18_85']):,.0f} |

`model_head_bridge_4year.csv` and `household_head_age_cells_4year.csv` are the
prospective interfaces to the four-year household model. They are not yet fed
back into household choices or housing-market clearing. No active solver,
certified baseline, or slide file is changed.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
