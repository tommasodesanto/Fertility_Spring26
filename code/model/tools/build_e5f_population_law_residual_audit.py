#!/usr/bin/env python3
"""Audit the population-law component of the E5F head forecast gap.

The diagnostic replaces Census zero-immigration births with the model's
top-code-adjusted birth flow, ages the signed cohort difference with Census
native-born survival, and applies the active queue's implicit headship mapping.
It thereby separates the population-law wedge into a birth-flow component and
a remaining initial-state/aging/exit component.  The latter is deliberately a
residual, not a mortality estimate.
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

from build_e5f_persons_demographic_satellite import (
    census_births_by_sex,
    native_survival_by_sex,
    read_csv,
    sha256,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SATELLITE = ROOT / "output/model/e5f_persons_demographic_satellite_20260826b"
DEFAULT_HEADSHIP = ROOT / "output/model/e5f_headship_demographic_bridge_20260826a"
SUMMARY_YEARS = (2025, 2050, 2080, 2100)
SEXES = (1, 2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--headship-dir", type=Path, default=DEFAULT_HEADSHIP)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=2025)
    parser.add_argument("--maximum-year", type=int, default=2100)
    return parser.parse_args()


def interpolate(
    rows: list[dict[str, str]], column: str, years: np.ndarray
) -> np.ndarray:
    source_years = np.asarray([float(row["calendar_year"]) for row in rows])
    values = np.asarray([float(row[column]) for row in rows])
    if np.any(np.diff(source_years) <= 0.0):
        raise ValueError("Transition calendar years must be strictly increasing")
    return np.interp(years, source_years, values)


def national_total_births(rows: list[dict[str, str]]) -> dict[int, float]:
    result = {
        int(row["YEAR"]): float(row["BIRTHS"])
        for row in rows
        if row["RACE_HISP"] == "0" and row["SEX"] == "0"
    }
    if not result:
        raise RuntimeError("Missing national total-birth rows")
    return result


def queue_mapped_heads(rows: list[dict[str, str]]) -> dict[int, float]:
    result: dict[int, float] = {}
    for row in rows:
        if row["ORIGIN"] != "0" or row["RACE"] != "0" or row["SEX"] not in ("1", "2"):
            continue
        year = int(row["YEAR"])
        result[year] = result.get(year, 0.0) + sum(
            float(row[f"POP_{age}"]) / 2.1 for age in range(20, 86)
        )
    if not result:
        raise RuntimeError("Missing sex-specific national population rows")
    return result


def make_figure(rows: list[dict[str, Any]], path: Path) -> None:
    years = np.asarray([int(row["calendar_year"]) for row in rows])
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.5))

    axes[0].plot(
        years,
        [float(row["model_topcode_adjusted_annual_births"]) / 1e6 for row in rows],
        color="#b22222",
        linewidth=2.1,
        label="Model baseline",
    )
    axes[0].plot(
        years,
        [float(row["census_zero_annual_births"]) / 1e6 for row in rows],
        color="#7aa6c2",
        linewidth=1.7,
        label="Census zero migration",
    )
    axes[0].plot(
        years,
        [float(row["census_mid_annual_births"]) / 1e6 for row in rows],
        color="#26547c",
        linewidth=2.0,
        label="Census main migration",
    )
    axes[0].set(title="Annual births", xlabel="Calendar year", ylabel="Millions")
    axes[0].legend(frameon=False, fontsize=8)

    selected = [row for row in rows if int(row["calendar_year"]) in (2050, 2080, 2100)]
    x = np.arange(len(selected), dtype=float)
    width = 0.28
    axes[1].bar(
        x - width / 2,
        [float(row["birth_flow_wedge"]) / 1e6 for row in selected],
        width=width,
        color="#7aa6c2",
        label="Birth flow",
    )
    axes[1].bar(
        x + width / 2,
        [float(row["initial_state_aging_exit_residual"]) / 1e6 for row in selected],
        width=width,
        color="#777777",
        label="State / aging / exit residual",
    )
    axes[1].axhline(0.0, color="black", linewidth=0.7)
    axes[1].set_xticks(x, [str(int(row["calendar_year"])) for row in selected])
    axes[1].set(
        title="Population-law wedge",
        xlabel="Calendar year",
        ylabel="Millions of heads",
    )
    axes[1].legend(frameon=False, fontsize=8)

    axes[2].plot(
        years,
        [float(row["model_annualized_mature_entrants"]) / 1e6 for row in rows],
        color="#26547c",
        linewidth=2.0,
        label="Mature entrants",
    )
    axes[2].plot(
        years,
        [float(row["model_annualized_terminal_exits"]) / 1e6 for row in rows],
        color="#b22222",
        linewidth=2.0,
        label="Terminal-age exits",
    )
    axes[2].plot(
        years,
        [float(row["model_annualized_net_head_flow"]) / 1e6 for row in rows],
        color="#777777",
        linewidth=1.4,
        linestyle="--",
        label="Net",
    )
    axes[2].axhline(0.0, color="black", linewidth=0.7)
    axes[2].set(
        title="Model head-flow accounting",
        xlabel="Calendar year",
        ylabel="Millions per year",
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
    for name, expected in headship_certification["satellite_source_hashes"].items():
        if sha256(source_dir / name) != expected:
            raise RuntimeError(f"Source fingerprint differs for {name}")
    anchor = float(
        headship_certification[
            "vintage_2025_implied_heads_age_18_85_at_2023_rates"
        ]
    )
    years = np.arange(int(args.start_year), int(args.maximum_year) + 1)
    baseline = read_csv(source_dir / "rebated_tax1_transition_path.csv")
    model_mass = interpolate(baseline, "adult_population", years)
    scale = anchor / model_mass[0]
    model_heads = model_mass * scale
    model_births = (
        interpolate(baseline, "birth_children_topcode_adjusted", years)
        * scale
        / 4.0
    )
    model_entrants = (
        interpolate(baseline, "effective_mature_entrant_flow_B", years)
        * scale
        / 4.0
    )
    model_exits = interpolate(baseline, "adult_deaths", years) * scale / 4.0

    zero_birth_rows = read_csv(source_dir / "births_zero.csv")
    mid_birth_rows = read_csv(source_dir / "births_mid.csv")
    zero_births_by_sex = census_births_by_sex(zero_birth_rows)
    zero_births_total = national_total_births(zero_birth_rows)
    mid_births_total = national_total_births(mid_birth_rows)
    survival = native_survival_by_sex(read_csv(source_dir / "survival.csv"))
    queue_raw = queue_mapped_heads(read_csv(source_dir / "population_zero.csv"))
    queue_scale = anchor / queue_raw[int(years[0])]
    census_zero_queue = np.asarray(
        [queue_raw[int(year)] * queue_scale for year in years]
    )

    signed_cohorts: dict[int, dict[int, float]] = {1: {}, 2: {}}
    signed_queue_heads: list[float] = []
    for index, year_value in enumerate(years):
        year = int(year_value)
        total_census_births = sum(zero_births_by_sex[(sex, year)] for sex in SEXES)
        birth_difference = model_births[index] - total_census_births
        for sex in SEXES:
            schedule = survival[(sex, year)]
            signed_cohorts[sex] = {
                age + 1: mass * schedule[age + 1]
                for age, mass in signed_cohorts[sex].items()
                if age < 100
            }
            sex_share = zero_births_by_sex[(sex, year)] / total_census_births
            signed_cohorts[sex][0] = birth_difference * sex_share * schedule[0]
        signed_queue_heads.append(
            sum(
                mass / 2.1
                for sex in SEXES
                for age, mass in signed_cohorts[sex].items()
                if 20 <= age <= 85
            )
            * queue_scale
        )
    census_zero_with_model_births = census_zero_queue + np.asarray(signed_queue_heads)

    rows: list[dict[str, Any]] = []
    for index, year_value in enumerate(years):
        year = int(year_value)
        birth_wedge = census_zero_queue[index] - census_zero_with_model_births[index]
        residual = census_zero_with_model_births[index] - model_heads[index]
        total = census_zero_queue[index] - model_heads[index]
        rows.append(
            {
                "calendar_year": year,
                "model_heads_anchored": model_heads[index],
                "census_zero_queue_heads_anchored": census_zero_queue[index],
                "census_zero_queue_heads_with_model_births": census_zero_with_model_births[
                    index
                ],
                "model_topcode_adjusted_annual_births": model_births[index],
                "census_zero_annual_births": zero_births_total[year],
                "census_mid_annual_births": mid_births_total[year],
                "birth_flow_wedge": birth_wedge,
                "initial_state_aging_exit_residual": residual,
                "total_population_law_wedge": total,
                "decomposition_residual": total - birth_wedge - residual,
                "model_annualized_mature_entrants": model_entrants[index],
                "model_annualized_terminal_exits": model_exits[index],
                "model_annualized_net_head_flow": model_entrants[index]
                - model_exits[index],
            }
        )

    write_csv(outdir / "population_law_residual_audit.csv", rows)
    make_figure(rows, outdir / "population_law_residual_audit.png")
    gates = {
        "decomposition_identity": bool(
            max(abs(float(row["decomposition_residual"])) for row in rows) <= 1e-6
        ),
        "signed_cohort_linearity": bool(
            np.all(np.isfinite(census_zero_with_model_births))
        ),
        "positive_level_paths": bool(
            np.min(census_zero_with_model_births) > 0.0 and np.min(model_heads) > 0.0
        ),
    }
    if not all(gates.values()):
        raise RuntimeError(f"Population-law audit gates failed: {gates}")

    selected = {
        int(row["calendar_year"]): row
        for row in rows
        if int(row["calendar_year"]) in SUMMARY_YEARS
    }
    baseline_by_year = {
        int(row["calendar_year"]): row for row in baseline
    }
    initial_row = baseline_by_year[2023]
    current_youngest_mass = float(initial_row["entry_flow_E"])
    scheduled_next_entrant_mass = float(
        initial_row["effective_mature_entrant_flow_B"]
    )
    late_exit_rows = [
        row
        for row in baseline
        if 2075 <= int(row["calendar_year"]) <= 2095
    ]
    exit_trough = min(late_exit_rows, key=lambda row: float(row["adult_deaths"]))
    exit_rebound = max(
        (
            row
            for row in late_exit_rows
            if int(row["calendar_year"]) > int(exit_trough["calendar_year"])
        ),
        key=lambda row: float(row["adult_deaths"]),
    )
    certification = {
        "status": "complete_exploratory_population_law_residual_audit_not_promoted",
        "anchor_heads_age_18_85": anchor,
        "model_scale_per_adult_mass_unit": scale,
        "initial_entry_discontinuity": {
            "current_2023_youngest_head_mass": current_youngest_mass,
            "scheduled_next_entrant_mass": scheduled_next_entrant_mass,
            "scheduled_over_current_youngest_mass": scheduled_next_entrant_mass
            / current_youngest_mass,
            "late_terminal_exit_trough_year": int(exit_trough["calendar_year"]),
            "late_terminal_exit_trough_mass": float(exit_trough["adult_deaths"]),
            "next_terminal_exit_rebound_year": int(exit_rebound["calendar_year"]),
            "next_terminal_exit_rebound_mass": float(exit_rebound["adult_deaths"]),
        },
        "accounting_gates": gates,
        "source_hashes": {
            name: sha256(source_dir / name)
            for name in (
                "rebated_tax1_transition_path.csv",
                "births_zero.csv",
                "births_mid.csv",
                "population_zero.csv",
                "survival.csv",
            )
        },
        "limitations": [
            "The birth-flow exercise is a linear accounting substitution, not a re-solved equilibrium.",
            "Signed birth-cohort differences use Census native-born survival and the Census zero-immigration sex ratio.",
            "The remaining initial-state/aging/exit term is a residual and must not be labeled mortality.",
            "Model four-year flows are linearly interpolated and divided by four only for annual display.",
            "The underlying perfect-foresight path remains preliminary and unpromoted.",
        ],
    }
    (outdir / "certification.json").write_text(
        json.dumps(certification, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    readme = f"""# Exploratory E5F population-law residual audit

Status: `complete_exploratory_population_law_residual_audit_not_promoted`.

This diagnostic decomposes the previously reported model-to-Census-zero queue
mapping wedge. It replaces Census zero-immigration births with the model's
top-code-adjusted births, follows the signed cohort difference with official
native-born survival, and leaves the initial Census population unchanged.

The active transition also has a direct initialization discontinuity. Its 2023
youngest household-head cell is {current_youngest_mass:.6f}, while the queue
schedules {scheduled_next_entrant_mass:.6f} heads for the next cohort, a ratio
of {scheduled_next_entrant_mass / current_youngest_mass:.2f}. When these cohorts
reach the terminal age, exits fall to {float(exit_trough['adult_deaths']):.6f}
in {int(exit_trough['calendar_year'])} and rebound to
{float(exit_rebound['adult_deaths']):.6f} in
{int(exit_rebound['calendar_year'])}.

| Year | Total population-law wedge | Birth-flow wedge | Initial-state / aging / exit residual | Model entrants per year | Model terminal exits per year |
|---|---:|---:|---:|---:|---:|
"""
    for year in SUMMARY_YEARS:
        row = selected[year]
        readme += (
            f"| {year} | {float(row['total_population_law_wedge']) / 1e6:+.3f}m | "
            f"{float(row['birth_flow_wedge']) / 1e6:+.3f}m | "
            f"{float(row['initial_state_aging_exit_residual']) / 1e6:+.3f}m | "
            f"{float(row['model_annualized_mature_entrants']) / 1e6:.3f}m | "
            f"{float(row['model_annualized_terminal_exits']) / 1e6:.3f}m |\n"
        )
    readme += """

## Interpretation

The near-term model contraction is not primarily caused by its birth-flow
shortfall. On the zero-immigration comparison, replacing Census births with the
model birth path explains only a small share of the 2050 gap. By 2080 the model
birth path is above the Census zero-immigration birth path and therefore closes,
rather than creates, part of the head-count gap.

Most of the remaining near-term wedge is the way the initial household-head
state ages and exits. This residual combines the inherited age distribution,
the absence of nonhead-to-head household formation, deterministic terminal-age
exit, and four-year timing. The 2023 youngest-cell/next-entrant discontinuity
shows that the observed headship age profile and the birth queue are not a
single coherent population law. It cannot be assigned to mortality without an
explicit person-level cohort law. No active model or slide file is changed.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
