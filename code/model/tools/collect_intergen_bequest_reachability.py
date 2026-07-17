#!/usr/bin/env python3
"""Collect the wide theta1 bequest reachability frontier."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from tools.profile_intergen_bequest_reachability import DISPERSION, FAMILY_GAP, MEDIAN, THETA1_GRID


MEDIAN_SE = 0.2319582116443
FAMILY_GAP_SE = 0.563010076238571


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    paths = sorted(args.results_dir.glob("cell_*/summary.json"))
    if len(paths) != len(THETA1_GRID):
        raise RuntimeError(f"expected {len(THETA1_GRID)} cell summaries, found {len(paths)}")

    rows: list[dict[str, Any]] = []
    eligible = 0
    for path in paths:
        summary = json.loads(path.read_text())
        meta = summary["metadata"]
        tight = summary.get("best_tight")
        repeat = summary.get("tight_repeat_check") or {}
        valid = bool(
            isinstance(tight, dict)
            and tight.get("strict_converged")
            and repeat.get("both_strict")
            and float(repeat.get("loss_abs_difference", 1.0)) <= 1e-10
            and float(repeat.get("max_abs_moment_difference", 1.0)) <= 1e-10
        )
        if valid:
            eligible += 1
        record = tight if valid else summary.get("best_search") or {}
        moments = record.get("moments") or {}
        gaps = record.get("gaps") or {}
        theta = record.get("theta") or {}
        median_gap = gaps.get(MEDIAN)
        family_gap = gaps.get(FAMILY_GAP)
        admissible_1se = bool(
            valid
            and median_gap is not None
            and family_gap is not None
            and abs(float(median_gap)) <= MEDIAN_SE
            and abs(float(family_gap)) <= FAMILY_GAP_SE
        )
        rows.append(
            {
                "cell": int(meta["cell"]),
                "theta1": float(meta["theta1"]),
                "eligible_tight": valid,
                "admissible_1se": admissible_1se,
                "theta0": theta.get("theta0"),
                "theta_n": theta.get("theta_n"),
                "profile_loss": record.get("profile_loss"),
                "market_residual": record.get("market_residual"),
                "median_model": moments.get(MEDIAN),
                "median_target": meta["targets"][MEDIAN],
                "median_gap": median_gap,
                "dispersion_model": moments.get(DISPERSION),
                "dispersion_target": meta["targets"][DISPERSION],
                "dispersion_gap": gaps.get(DISPERSION),
                "family_gap_model": moments.get(FAMILY_GAP),
                "family_gap_target": meta["targets"][FAMILY_GAP],
                "family_gap_gap": family_gap,
                "n_cases": summary["n_cases_completed"],
                "elapsed_sec": summary["elapsed_sec"],
            }
        )
    rows.sort(key=lambda row: float(row["theta1"]))
    non_strict_cells = [int(row["cell"]) for row in rows if not row["eligible_tight"]]
    admissible = [row for row in rows if row["admissible_1se"]]
    max_row = max(admissible, key=lambda row: float(row["dispersion_model"])) if admissible else None
    closest_row = min(rows, key=lambda row: float(row["profile_loss"]))
    target = float(rows[0]["dispersion_target"])
    target_reached = bool(max_row is not None and float(max_row["dispersion_model"]) >= target)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "reachability_frontier.csv", rows)
    result = {
        "status": "complete",
        "eligible_cells": eligible,
        "non_strict_cells": non_strict_cells,
        "admissible_1se_cells": len(admissible),
        "dispersion_target": target,
        "maximum_admissible_dispersion": max_row,
        "closest_two_target_fit": closest_row,
        "dispersion_target_reached": target_reached,
        "frontier": rows,
    }
    (args.outdir / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True))

    lines = [
        "# Bequest reachability frontier",
        "",
        f"Strict, exactly repeated cells: {eligible}/{len(THETA1_GRID)}.",
        f"Non-strict cells (excluded from admissibility): {non_strict_cells}.",
        f"Cells matching the median and family gap within one bootstrap SE: {len(admissible)}.",
        f"Dispersion target reached: {target_reached}.",
        "",
    ]
    if max_row is not None:
        lines.extend(
            [
                "| theta1 | theta0 | theta_n | Median | Family gap | p90/median | Target |",
                "|---:|---:|---:|---:|---:|---:|---:|",
                f"| {float(max_row['theta1']):.4g} | {float(max_row['theta0']):.4g} | "
                f"{float(max_row['theta_n']):.4g} | {float(max_row['median_model']):.4f} | "
                f"{float(max_row['family_gap_model']):.4f} | {float(max_row['dispersion_model']):.4f} | {target:.4f} |",
                "",
            ]
        )
    lines.append("The complete strict frontier is in `reachability_frontier.csv`.")
    (args.outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
