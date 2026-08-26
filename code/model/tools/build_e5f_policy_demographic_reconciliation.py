#!/usr/bin/env python3
"""Reconcile the E5F policy effect across persons and household-head units."""

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
DEFAULT_HEADSHIP = ROOT / "output/model/e5f_headship_demographic_bridge_20260826a"
DEFAULT_BASELINE_GAP = (
    ROOT / "output/model/e5f_baseline_demographic_gap_decomposition_20260826c"
)
SUMMARY_YEARS = (2050, 2080, 2100)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--headship-dir", type=Path, default=DEFAULT_HEADSHIP)
    parser.add_argument("--baseline-gap-dir", type=Path, default=DEFAULT_BASELINE_GAP)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def make_figure(rows: list[dict[str, Any]], path: Path) -> None:
    years = [int(row["calendar_year"]) for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.4))
    axes[0].plot(
        years,
        [float(row["extra_resident_persons"]) / 1e6 for row in rows],
        color="#7aa6c2",
        linewidth=2.1,
        label="Resident persons",
    )
    axes[0].plot(
        years,
        [float(row["extra_children_age_0_17"]) / 1e6 for row in rows],
        color="#c7dce8",
        linewidth=1.7,
        label="Children 0--17",
    )
    axes[0].plot(
        years,
        [float(row["extra_adult_persons_age_18_85"]) / 1e6 for row in rows],
        color="#26547c",
        linewidth=1.8,
        label="Adult persons 18--85",
    )
    axes[0].set(
        title="Policy-born resident persons",
        xlabel="Calendar year",
        ylabel="Millions",
    )
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].plot(
        years,
        [float(row["model_transition_extra_heads"]) / 1e3 for row in rows],
        color="#b22222",
        linewidth=2.1,
        label="Active transition",
    )
    axes[1].plot(
        years,
        [float(row["queue_proxy_extra_heads"]) / 1e3 for row in rows],
        color="#777777",
        linewidth=1.5,
        linestyle="--",
        label="Birth queue only",
    )
    axes[1].plot(
        years,
        [float(row["acs_mapped_extra_heads"]) / 1e3 for row in rows],
        color="#26547c",
        linewidth=2.1,
        label="ACS headship mapping",
    )
    axes[1].set(
        title="Policy effect in household-head units",
        xlabel="Calendar year",
        ylabel="Thousands",
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
    headship = args.headship_dir.resolve()
    baseline_gap = args.baseline_gap_dir.resolve()
    outdir = args.output_dir.resolve()
    for path in (headship, baseline_gap):
        if not path.exists():
            raise FileNotFoundError(path)
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    policy_rows = {
        int(row["calendar_year"]): row
        for row in read_csv(headship / "policy_person_to_head_wedge.csv")
        if row["immigration_scenario"] == "mid"
    }
    model_rows = {
        int(row["calendar_year"]): row
        for row in read_csv(
            baseline_gap / "baseline_demographic_gap_decomposition.csv"
        )
    }
    years = sorted(set(policy_rows) & set(model_rows))
    rows: list[dict[str, Any]] = []
    for year in years:
        policy = policy_rows[year]
        model = model_rows[year]
        children = float(policy["extra_children_age_0_17"])
        adults = float(policy["extra_adult_persons_age_18_85"])
        mapped_heads = float(policy["extra_heads_2023_acs_headship"])
        queue_heads = float(
            policy[
                "queue_proxy_extra_heads_births_divided_by_2_1_at_age_20"
            ]
        )
        model_heads = float(model["model_policy_wedge_heads"])
        rows.append(
            {
                "calendar_year": year,
                "extra_children_age_0_17": children,
                "extra_adult_persons_age_18_85": adults,
                "extra_resident_persons": children + adults,
                "acs_mapped_extra_heads": mapped_heads,
                "queue_proxy_extra_heads": queue_heads,
                "model_transition_extra_heads": model_heads,
                "model_minus_acs_mapped_heads": model_heads - mapped_heads,
                "model_over_acs_mapped_heads": (
                    model_heads / mapped_heads if mapped_heads > 0.0 else None
                ),
            }
        )
    write_csv(outdir / "policy_demographic_reconciliation.csv", rows)
    make_figure(rows, outdir / "policy_demographic_reconciliation.png")

    gates = {
        "persons_equal_children_plus_adults": all(
            abs(
                float(row["extra_resident_persons"])
                - float(row["extra_children_age_0_17"])
                - float(row["extra_adult_persons_age_18_85"])
            )
            <= 1e-6
            for row in rows
        ),
        "mapped_heads_do_not_exceed_adult_persons": all(
            float(row["acs_mapped_extra_heads"])
            <= float(row["extra_adult_persons_age_18_85"]) + 1e-6
            for row in rows
        ),
        "common_year_support": bool(years and years[0] == 2025 and years[-1] == 2100),
    }
    if not all(gates.values()):
        raise RuntimeError(f"Policy demographic reconciliation gates failed: {gates}")
    certification = {
        "status": "complete_exploratory_policy_demographic_reconciliation_not_promoted",
        "accounting_gates": gates,
        "source_hashes": {
            "policy_person_to_head_wedge.csv": sha256(
                headship / "policy_person_to_head_wedge.csv"
            ),
            "baseline_demographic_gap_decomposition.csv": sha256(
                baseline_gap / "baseline_demographic_gap_decomposition.csv"
            ),
        },
        "limitations": [
            "The persons and ACS-headship paths are satellites and do not feed back into equilibrium.",
            "The active transition head wedge includes its endogenous queue and equilibrium distribution response.",
            "The reconciliation diagnoses timing and unit differences; it is not a welfare or causal decomposition.",
            "The underlying perfect-foresight paths remain preliminary and unpromoted.",
        ],
    }
    (outdir / "certification.json").write_text(
        json.dumps(certification, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    selected = {int(row["calendar_year"]): row for row in rows}
    readme = """# Exploratory E5F policy demographic reconciliation

Status: `complete_exploratory_policy_demographic_reconciliation_not_promoted`.

The table keeps resident persons, adult persons, and household heads in
separate units. The ACS-mapped head count applies frozen 2023 age- and
sex-specific headship to the policy-born resident cohorts. The active
transition count is the scaled difference in its endogenous adult-head mass.

| Year | Extra resident persons | Extra adult persons | ACS-mapped extra heads | Queue proxy heads | Active-transition extra heads | Active / ACS heads |
|---|---:|---:|---:|---:|---:|---:|
"""
    for year in SUMMARY_YEARS:
        row = selected[year]
        readme += (
            f"| {year} | {float(row['extra_resident_persons']):,.0f} | "
            f"{float(row['extra_adult_persons_age_18_85']):,.0f} | "
            f"{float(row['acs_mapped_extra_heads']):,.0f} | "
            f"{float(row['queue_proxy_extra_heads']):,.0f} | "
            f"{float(row['model_transition_extra_heads']):,.0f} | "
            f"{float(row['model_over_acs_mapped_heads']):.2f} |\n"
        )
    readme += """

The active model puts the policy-born cohorts into household-head status too
quickly. The discrepancy is large in 2050, smaller in 2080, and nearly gone by
2100. This pattern supports a timing diagnosis: the long-run number of heads
per birth is similar, but the single age-20 queue conversion is not a credible
short- and medium-run household-formation schedule. No active solver or slide
file is changed.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
