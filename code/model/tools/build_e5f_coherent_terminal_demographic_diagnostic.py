#!/usr/bin/env python3
"""Test E5F terminal fertility rates against the coherent person law.

This is a demographic fixed-point diagnostic, not a household-equilibrium
solve.  It asks whether the current impact and previously imposed terminal
birth rates admit a finite population under fixed late-Census net migration,
age-specific survival, and the ACS headship mapping.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from demographic_transition import (
    solve_demographic_steady_state,
    stationary_demographic_components,
)
from build_e5f_coherent_person_cohort_path import (
    births_path,
    fixed_2023_headship,
    read_csv,
    survival_path,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE_DIR = (
    ROOT
    / "output/model/e5f_persons_demographic_satellite_20260826b/source_data"
)
DEFAULT_HEADSHIP_DIR = (
    ROOT / "output/model/e5f_headship_demographic_bridge_20260826a"
)
DEFAULT_COHORT_DIR = (
    ROOT / "output/model/e5f_coherent_person_cohort_path_20260826a"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--headship-dir", type=Path, default=DEFAULT_HEADSHIP_DIR)
    parser.add_argument("--cohort-dir", type=Path, default=DEFAULT_COHORT_DIR)
    parser.add_argument("--tax1-summary", type=Path, required=True)
    parser.add_argument("--tax2-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--terminal-average-years", type=int, default=5)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_terminal_case(path: Path) -> tuple[str, dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    states = dict(payload.get("terminal_steady_states") or {})
    if len(states) != 1:
        raise RuntimeError(f"Expected one terminal state in {path}; found {states}")
    name, row = next(iter(states.items()))
    gates = dict(row.get("fixed_point_gates") or {})
    if gates.get("status") != "passed" or not all(
        bool(value) for value in dict(gates.get("checks") or {}).values()
    ):
        raise RuntimeError(f"Source terminal household fixed point did not pass: {path}")
    return str(name), dict(row)


def migration_path(rows: list[dict[str, str]]) -> dict[int, np.ndarray]:
    result: dict[int, np.ndarray] = {}
    for row in rows:
        year = int(row["calendar_year"])
        result.setdefault(year, np.zeros((2, 101), dtype=float))
        result[year][int(row["sex"]) - 1, int(row["age"])] = float(
            row["net_migration_residual"]
        )
    return result


def impact_birth_rate(path: Path) -> float:
    rows = read_csv(path)
    first = rows[0]
    if int(first["period"]) != 0:
        raise RuntimeError(f"Policy path does not start at period zero: {path}")
    four_year_rate = float(first["birth_children_topcode_adjusted"]) / float(
        first["adult_population"]
    )
    return four_year_rate / 4.0


def make_figure(rows: list[dict[str, Any]], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
    labels = [str(row["case"]) for row in rows]
    ratios = [float(row["renewal_ratio"]) for row in rows]
    colors = ["#26547c" if value < 1.0 else "#b22222" for value in ratios]
    axes[0].bar(np.arange(len(rows)), ratios, color=colors)
    axes[0].axhline(1.0, color="black", linewidth=1.0, linestyle="--")
    axes[0].set_xticks(np.arange(len(rows)), labels, rotation=25, ha="right")
    axes[0].set(title="Demographic renewal ratio", ylabel="Birth rate × heads per birth")

    valid = [row for row in rows if bool(row["finite_fixed_migration_ss"])]
    x = np.arange(len(valid), dtype=float)
    axes[1].bar(
        x - 0.18,
        [float(row["stationary_resident_persons"]) / 1e6 for row in valid],
        width=0.36,
        color="#26547c",
        label="Resident persons",
    )
    axes[1].bar(
        x + 0.18,
        [float(row["stationary_heads_age_18_85"]) / 1e6 for row in valid],
        width=0.36,
        color="#7aa6c2",
        label="Heads 18–85",
    )
    axes[1].set_xticks(x, [str(row["case"]) for row in valid], rotation=20, ha="right")
    axes[1].set(title="Fixed-migration steady states", ylabel="Millions")
    axes[1].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.18, axis="y")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    source_dir = args.source_dir.resolve()
    headship_dir = args.headship_dir.resolve()
    cohort_dir = args.cohort_dir.resolve()
    tax1_summary = args.tax1_summary.resolve()
    tax2_summary = args.tax2_summary.resolve()
    outdir = args.output_dir.resolve()
    required = (
        source_dir / "survival.csv",
        source_dir / "births_mid.csv",
        source_dir / "rebated_tax1_transition_path.csv",
        source_dir / "rebated_tax2_transition_path.csv",
        headship_dir / "acs_headship_profiles.csv",
        cohort_dir / "net_migration_residual_by_age_sex.csv",
        cohort_dir / "certification.json",
        tax1_summary,
        tax2_summary,
    )
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    copied = outdir / "source_data"
    copied.mkdir()
    for source, name in (
        (tax1_summary, "rebated_tax1_terminal_summary.json"),
        (tax2_summary, "rebated_tax2_terminal_summary.json"),
    ):
        with source.open("rb") as src, (copied / name).open("xb") as dst:
            shutil.copyfileobj(src, dst)

    terminal_years = int(args.terminal_average_years)
    if terminal_years < 1 or terminal_years > 20:
        raise ValueError("terminal-average-years must lie in [1, 20]")
    final_years = list(range(2101 - terminal_years, 2101))
    survival = survival_path(read_csv(source_dir / "survival.csv"), nativity=0)
    births = births_path(read_csv(source_dir / "births_mid.csv"))
    migrations = migration_path(
        read_csv(cohort_dir / "net_migration_residual_by_age_sex.csv")
    )
    average_survival = np.mean([survival[year] for year in final_years], axis=0)
    average_births = np.mean([births[year] for year in final_years], axis=0)
    birth_shares = average_births / float(np.sum(average_births))
    average_migration = np.mean([migrations[year] for year in final_years], axis=0)
    headship = fixed_2023_headship(
        read_csv(headship_dir / "acs_headship_profiles.csv")
    )

    components = stationary_demographic_components(
        birth_sex_shares=birth_shares,
        survival=average_survival,
        net_migration=average_migration,
        headship_rates=headship,
    )
    critical_annual_rate = 1.0 / components.birth_head_multiplier
    critical_four_year_rate = 4.0 * critical_annual_rate

    tax1_name, tax1_terminal = load_terminal_case(tax1_summary)
    tax2_name, tax2_terminal = load_terminal_case(tax2_summary)
    candidates = [
        (
            "2023 rebated 1%",
            impact_birth_rate(source_dir / "rebated_tax1_transition_path.csv"),
            "common reconstructed 2023 state",
        ),
        (
            "2023 rebated 2%",
            impact_birth_rate(source_dir / "rebated_tax2_transition_path.csv"),
            "common reconstructed 2023 state",
        ),
        (
            f"old terminal {tax1_name}",
            float(tax1_terminal["births_per_adult"]) / 4.0,
            "previous queue-based terminal household fixed point",
        ),
        (
            f"old terminal {tax2_name}",
            float(tax2_terminal["births_per_adult"]) / 4.0,
            "previous queue-based terminal household fixed point",
        ),
    ]
    rows: list[dict[str, Any]] = []
    for label, rate, source_label in candidates:
        renewal = rate * components.birth_head_multiplier
        finite = renewal < 1.0
        solved = None
        if finite:
            solved = solve_demographic_steady_state(
                annual_births_per_head=rate,
                birth_sex_shares=birth_shares,
                survival=average_survival,
                net_migration=average_migration,
                headship_rates=headship,
            )
        rows.append(
            {
                "case": label,
                "source": source_label,
                "annual_births_per_head": rate,
                "four_year_births_per_head": 4.0 * rate,
                "birth_head_multiplier": components.birth_head_multiplier,
                "renewal_ratio": renewal,
                "critical_annual_births_per_head": critical_annual_rate,
                "critical_four_year_births_per_head": critical_four_year_rate,
                "finite_fixed_migration_ss": finite,
                "stationary_resident_persons": ""
                if solved is None
                else float(np.sum(solved.persons)),
                "stationary_heads_age_18_85": ""
                if solved is None
                else float(np.sum(solved.heads[:, 18:86])),
                "stationary_annual_births": ""
                if solved is None
                else solved.annual_births,
                "population_fixed_point_residual": ""
                if solved is None
                else solved.population_fixed_point_residual,
            }
        )

    write_csv(outdir / "terminal_demographic_fixed_point.csv", rows)
    make_figure(rows, outdir / "terminal_demographic_fixed_point.png")
    gates = {
        "source_terminal_household_fixed_points_passed": True,
        "current_2023_rates_have_finite_fixed_migration_ss": all(
            bool(row["finite_fixed_migration_ss"]) for row in rows[:2]
        ),
        "old_queue_terminal_rates_rejected_by_person_law": all(
            not bool(row["finite_fixed_migration_ss"]) for row in rows[2:]
        ),
        "valid_fixed_points_close": all(
            row["population_fixed_point_residual"] == ""
            or float(row["population_fixed_point_residual"]) <= 1e-5
            for row in rows
        ),
    }
    if not all(gates.values()):
        raise RuntimeError(f"Terminal demographic diagnostic gates failed: {gates}")

    source_hashes = {str(path): sha256(path) for path in required}
    certification = {
        "status": "complete_isolated_terminal_demographic_diagnostic_not_integrated",
        "source_hashes": source_hashes,
        "terminal_average_years": final_years,
        "average_annual_net_migration": float(np.sum(average_migration)),
        "birth_head_multiplier": components.birth_head_multiplier,
        "fixed_migration_head_base": components.fixed_migration_head_base,
        "critical_annual_births_per_head": critical_annual_rate,
        "critical_four_year_births_per_head": critical_four_year_rate,
        "accounting_gates": gates,
        "limitations": [
            "This holds household birth rates fixed and does not yet re-solve prices, transfers, policies, or the stationary household distribution.",
            "Late-Census net migration is held fixed in absolute age-sex units; alternative migration closures will imply different population levels.",
            "Total-population survival is used for the terminal stock because the person state does not yet split native and foreign-born residents.",
            "Fixed 2023 ACS headship rates are an accounting mapping, not an estimated long-run formation model.",
        ],
    }
    (outdir / "certification.json").write_text(
        json.dumps(certification, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    current1, current2, old1, old2 = rows
    readme = f"""# Coherent terminal demographic fixed-point diagnostic

Status: `complete_isolated_terminal_demographic_diagnostic_not_integrated`.

One annual birth generates `{components.birth_head_multiplier:.3f}` future
head-years under late-Census survival and fixed 2023 ACS headship. The critical
birth rate is therefore `{critical_annual_rate:.6f}` per head per year, or
`{critical_four_year_rate:.6f}` per four-year model period.

The reconstructed 2023 rebated-1% and rebated-2% rates imply renewal ratios
`{float(current1['renewal_ratio']):.4f}` and
`{float(current2['renewal_ratio']):.4f}`. With the average 2096–2100 net
migration residual held fixed, both admit finite positive demographic steady
states.

The old queue-based terminal rates imply renewal ratios
`{float(old1['renewal_ratio']):.4f}` and
`{float(old2['renewal_ratio']):.4f}`. Both exceed one, so they do **not** admit
a finite steady state under the coherent person law with positive fixed
migration. The old terminal preference value was chosen to close the
`births / 2.1` household-entry queue; it cannot be retained as the terminal
condition after that queue is removed.

This identifies the next algorithmic step: jointly re-solve the terminal
household policy, housing price, fiscal transfer, and demographic scale under
the new person/head fixed point. The levels reported for the valid 2023-rate
diagnostics are closure sensitivities, not forecasts.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
