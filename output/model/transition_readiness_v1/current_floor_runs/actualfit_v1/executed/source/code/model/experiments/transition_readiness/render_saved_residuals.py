#!/usr/bin/env python3
"""Plot saved six-date housing and pension residuals; performs no model solves."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CASES = ("baseline", "trial", "fresh_replay")
COLORS = {"baseline": "#4C78A8", "trial": "#F58518", "fresh_replay": "#54A24B"}
LABELS = {"baseline": "Baseline", "trial": "Trial", "fresh_replay": "Fresh replay"}


def load_case(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = json.loads(path.read_text())
    rows = data.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"{path}: expected nonempty rows")
    years = np.asarray([row["calendar_year"] for row in rows], dtype=float)
    housing = np.asarray([row["relative_market_residual"] for row in rows], dtype=float)
    pension = np.asarray([row["scaled_pension_budget_residual"] for row in rows], dtype=float)
    if not (len(years) == len(housing) == len(pension)) or not np.isfinite(
        np.concatenate((years, housing, pension))
    ).all():
        raise ValueError(f"{path}: residual rows have inconsistent lengths or nonfinite values")
    return years, housing, pension


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root",
        type=Path,
        default=Path("output/model/transition_readiness_v1/prior_evidence/smoke"),
        help="Folder containing baseline/, trial/, and fresh_replay/ mapping.json files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/model/transition_readiness_v1/prior_evidence"),
    )
    args = parser.parse_args()

    cases = {
        name: load_case(args.input_root / name / "mapping.json") for name in CASES
    }
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.2), sharex=True, constrained_layout=True)
    specifications = (
        (axes[0], 1, "Relative housing market residual", 2e-4),
        (axes[1], 2, "Scaled pension budget residual", 1e-6),
    )
    for ax, column, ylabel, gate in specifications:
        for name, (years, housing, pension) in cases.items():
            residual = (housing, pension)[column - 1]
            ax.plot(
                years,
                residual,
                marker="o",
                markersize=4,
                linewidth=1.5,
                color=COLORS[name],
                label=LABELS[name],
            )
        ax.axhline(0.0, color="black", linewidth=0.8)
        ax.axhline(gate, color="#777777", linestyle="--", linewidth=0.9,
                   label=f"Gate ±{gate:g}")
        ax.axhline(-gate, color="#777777", linestyle="--", linewidth=0.9)
        ax.set_ylabel(ylabel)
        ax.set_yscale("symlog", linthresh=gate / 10)
        ax.grid(True, which="both", alpha=0.22)
        ax.legend(frameon=False, ncol=2, fontsize=8, loc="best")
    axes[1].set_xlabel("Calendar year")
    fig.suptitle("Supplemental historical fixed-reference residuals", fontsize=12)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(args.output_dir / f"historical_fixed_reference_residuals.{suffix}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
