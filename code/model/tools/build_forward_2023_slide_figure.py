#!/usr/bin/env python3
"""Build the fitted-economy forward path used in the transition slides."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--path-csv", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--scenario", default="closed_national_benchmark")
    parser.add_argument("--last-year", type=int, default=2063)
    return parser.parse_args()


def load_path(path: Path, scenario: str, last_year: int) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        raw_rows = list(csv.DictReader(stream))
    rows = [
        row
        for row in raw_rows
        if row["scenario"] == scenario
        and 2023 <= int(float(row["calendar_year"])) <= last_year
    ]
    if not rows or int(float(rows[0]["calendar_year"])) != 2023:
        raise ValueError("The path must begin in 2023.")
    return [
        {
            "year": float(row["calendar_year"]),
            "price": float(row["asset_price"]),
            "births": float(row["topcode_adjusted_births_per_adult"]),
            "population": float(row["model_population_units"]),
        }
        for row in rows
    ]


def main() -> None:
    args = parse_args()
    rows = load_path(args.path_csv, args.scenario, args.last_year)
    years = [row["year"] for row in rows]
    fields = (
        ("House price", "price"),
        ("Births per adult household", "births"),
        ("Population: adults and children", "population"),
    )

    plt.rcParams.update(
        {
            "font.size": 10.5,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.15), constrained_layout=True)
    for axis, (title, field) in zip(axes, fields, strict=True):
        values = [row[field] / rows[0][field] for row in rows]
        axis.plot(
            years,
            values,
            color="#26527a",
            marker="o",
            linewidth=2.2,
            markersize=5,
        )
        axis.axhline(1.0, color="#999999", linewidth=0.9, linestyle="--")
        axis.set_title(title)
        axis.set_xlabel("Year")
        axis.set_ylabel("Index (2023 = 1)")
        axis.set_xticks(years[::2])
        axis.grid(alpha=0.22)
        axis.spines[["top", "right"]].set_visible(False)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=220, bbox_inches="tight")


if __name__ == "__main__":
    main()
