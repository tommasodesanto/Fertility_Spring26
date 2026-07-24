#!/usr/bin/env python3
"""Build the two deterministic calibration-note figures for the eqscale note."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REFERENCE_MEAN_AGE = 25.310560799362
DARK_COLOR = "#1f2937"
PDF_METADATA = {"CreationDate": None, "ModDate": None}


def repository_root() -> Path:
    """Return the repository root from this file's code/model/tools location."""
    return Path(__file__).resolve().parents[2].parent


def save_fecundity_figure(outdir: Path) -> Path:
    """Plot the model fecundity schedule and its four-year attempt ages."""
    ages = np.linspace(18.0, 46.0, 561)
    probabilities = np.where(
        ages < 45.0,
        np.clip(1.0 - 0.02 * np.exp(0.134 * (ages - 18.0)), 0.0, 1.0),
        0.0,
    )
    attempt_ages = np.array([18, 22, 26, 30, 34, 38, 42, 46], dtype=float)
    attempt_probabilities = np.where(
        attempt_ages < 45.0,
        np.clip(1.0 - 0.02 * np.exp(0.134 * (attempt_ages - 18.0)), 0.0, 1.0),
        0.0,
    )

    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.plot(ages, probabilities, color=DARK_COLOR, linewidth=1.8)
    ax.plot(
        attempt_ages,
        attempt_probabilities,
        linestyle="none",
        marker="o",
        markersize=4.2,
        color=DARK_COLOR,
    )
    ax.axvline(45, color=DARK_COLOR, linestyle="--", linewidth=1.0)
    ax.text(44.7, 0.98, "end of fertile window", fontsize=8, ha="right", va="top")
    ax.set_xlabel("Age")
    ax.set_ylabel("Per-period conception probability")
    ax.set_xlim(18, 46)
    ax.set_ylim(0, 1.02)
    fig.tight_layout()

    output_path = outdir / "eqscale_note_fecundity.pdf"
    fig.savefig(output_path, format="pdf", metadata=PDF_METADATA)
    plt.close(fig)
    return output_path


def first_birth_distribution(csv_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return single-year first-birth shares for the 1979--1984 birth cohorts."""
    counts_by_age: defaultdict[int, int] = defaultdict(int)
    with csv_path.open(newline="") as csv_file:
        for row in csv.DictReader(csv_file):
            year = int(row["year"])
            age = int(row["age"])
            if 1979 <= year - age <= 1984 and 12 <= age <= 49:
                counts_by_age[age] += int(row["n_first_births"])

    ages = np.array(sorted(counts_by_age), dtype=float)
    counts = np.array([counts_by_age[int(age)] for age in ages], dtype=float)
    return ages, counts / counts.sum()


def save_first_birth_age_figure(outdir: Path, csv_path: Path) -> tuple[Path, float]:
    """Plot the first-birth age distribution and return its count-weighted mean."""
    ages, shares = first_birth_distribution(csv_path)
    mean_age = float(np.dot(ages, shares))

    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.plot(ages, shares, color=DARK_COLOR, linewidth=1.8)
    ax.axvspan(30, 49, color="#9ca3af", alpha=0.18, linewidth=0)
    ax.axvline(REFERENCE_MEAN_AGE, color=DARK_COLOR, linestyle="--", linewidth=1.0)
    ax.text(REFERENCE_MEAN_AGE + 0.25, shares.max() * 0.95, "mean 25.31", fontsize=8, va="top")
    ax.text(37.2, shares.max() * 0.55, "share 30+ = 0.270", fontsize=8, ha="center")
    ax.set_xlabel("Mother's age at first birth")
    ax.set_ylabel("Share of first births")
    ax.set_xlim(12, 49)
    ax.set_ylim(bottom=0)
    fig.tight_layout()

    output_path = outdir / "eqscale_note_first_birth_age.pdf"
    fig.savefig(output_path, format="pdf", metadata=PDF_METADATA)
    plt.close(fig)
    return output_path, mean_age


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=repository_root() / "latex" / "figures",
        help="Directory for generated PDFs (default: latex/figures).",
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    fecundity_path = save_fecundity_figure(args.outdir)
    first_birth_path, mean_age = save_first_birth_age_figure(
        args.outdir,
        repository_root() / "code" / "data" / "nchs_natality_timing" / "first_birth_counts_year_age.csv",
    )
    print(
        "Consistency check: plotted mean age = "
        f"{mean_age:.12f}; reference = {REFERENCE_MEAN_AGE:.12f}"
    )
    if abs(mean_age - REFERENCE_MEAN_AGE) > 0.02:
        raise RuntimeError("Plotted mean age differs from the reference by more than 0.02.")
    print(f"Wrote {fecundity_path}")
    print(f"Wrote {first_birth_path}")


if __name__ == "__main__":
    main()
