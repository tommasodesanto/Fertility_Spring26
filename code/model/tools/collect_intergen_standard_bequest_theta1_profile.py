#!/usr/bin/env python3
"""Collect the strict five-cell conditional theta1 profile for M4."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


TARGET_SET = "candidate_replacement_bequest_median_composition_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--winner-json", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def exact_tight(summary: dict[str, Any]) -> dict[str, Any]:
    repeats = summary.get("tight_repeats") or []
    check = summary.get("tight_repeat_check") or {}
    if not (
        len(repeats) == 2
        and all(bool(row.get("strict_converged", False)) for row in repeats)
        and bool(check.get("both_strict", False))
        and float(check.get("loss_abs_difference", math.inf)) == 0.0
        and float(check.get("max_abs_moment_difference", math.inf)) == 0.0
    ):
        raise RuntimeError("profile cell lacks two bit-identical strict tight repeats")
    return repeats[0]


def main() -> None:
    args = parse_args()
    payload = json.loads(args.winner_json.read_text())
    winner = (payload.get("winners") or {}).get("M4")
    if not isinstance(winner, dict):
        raise RuntimeError("winner JSON has no collected M4 winner")
    winner_loss = float(winner["rank_loss"])
    winner_theta1 = float(winner["theta"]["theta1"])

    summary_paths = sorted(args.results_dir.glob("task_*/summary.json"))
    if len(summary_paths) != 5:
        raise RuntimeError(f"expected exactly five theta1 profile cells, found {len(summary_paths)}")
    rows: list[dict[str, Any]] = []
    for path in summary_paths:
        summary = json.loads(path.read_text())
        metadata = summary.get("metadata") or {}
        cell = json.loads((path.parent / "profile_cell.json").read_text())
        fixed = metadata.get("fixed_parameters") or {}
        if not (
            metadata.get("arm") == "M4_PROFILE"
            and metadata.get("target_set") == TARGET_SET
            and int(metadata.get("free_parameter_count", -1)) == 12
            and int(metadata.get("target_count", -1)) == 14
            and float(fixed.get("theta1", math.nan)) == float(cell["theta1"])
            and float(fixed.get("theta_n", math.nan)) == 0.0
            and float(fixed.get("tenure_choice_kappa", math.nan)) == 0.0
        ):
            raise RuntimeError(f"profile contract mismatch in {path.parent}")
        tight = exact_tight(summary)
        loss = float(tight["rank_loss"])
        rows.append(
            {
                "task_id": int(cell["task_id"]),
                "theta1_unit_offset": float(cell["theta1_unit_offset"]),
                "theta1": float(cell["theta1"]),
                "theta0_profile_estimate": float(tight["theta"]["theta0"]),
                "tight_loss": loss,
                "loss_gap_from_joint_winner": loss - winner_loss,
                "market_residual": float(tight["market_residual"]),
                "search_cases": int(summary["n_cases_completed"]),
                "search_strict": int(summary["n_strict"]),
                "exact_tight_repeat": True,
            }
        )
    rows.sort(key=lambda row: float(row["theta1_unit_offset"]))
    if len({float(row["theta1"]) for row in rows}) != 5:
        raise RuntimeError("theta1 profile cells are not distinct")
    center = [row for row in rows if float(row["theta1_unit_offset"]) == 0.0]
    if len(center) != 1:
        raise RuntimeError("theta1 profile lacks exactly one central winner cell")
    if float(center[0]["theta1"]) != winner_theta1:
        raise RuntimeError("central theta1 profile cell does not reproduce the joint winner value")
    if abs(float(center[0]["tight_loss"]) - winner_loss) > 1e-10:
        raise RuntimeError("central theta1 profile cell does not reproduce the joint winner loss")

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "theta1_profile.csv", rows)
    summary = {
        "status": "complete",
        "winner_loss": winner_loss,
        "winner_theta1": winner_theta1,
        "profile_rows": rows,
        "joint_search_quality_pass": min(
            float(row["loss_gap_from_joint_winner"]) for row in rows
        ) >= -1e-10,
        "max_loss_gap_within_profile": max(float(row["loss_gap_from_joint_winner"]) for row in rows),
        "min_noncentral_loss_gap": min(
            float(row["loss_gap_from_joint_winner"])
            for row in rows
            if float(row["theta1_unit_offset"]) != 0.0
        ),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    lines = [
        "# M4 conditional theta1 profile",
        "",
        "Each noncentral cell fixes `theta1` and re-estimates the other 12 M4 parameters. ",
        "Every reported loss comes from two bit-identical strict tight solves.",
        f"Joint-search quality gate (no profile improvement): "
        f"{summary['joint_search_quality_pass']}.",
        "",
        "| Unit offset | theta1 | theta0 | Tight loss | Delta loss | Residual |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {float(row['theta1_unit_offset']):.3f} | {float(row['theta1']):.6f} | "
            f"{float(row['theta0_profile_estimate']):.6f} | {float(row['tight_loss']):.6f} | "
            f"{float(row['loss_gap_from_joint_winner']):.6f} | "
            f"{float(row['market_residual']):.2e} |"
        )
    (args.outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
