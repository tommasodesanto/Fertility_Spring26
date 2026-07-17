#!/usr/bin/env python3
"""Collect the diagnostic retirement-income-dispersion frontier."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from tools.profile_intergen_bequest_reachability import (
    DISPERSION,
    FAMILY_GAP,
    MEDIAN,
    RETIREMENT_INCOME_SCALE_GRID,
)
from tools.run_intergen_bequest_exit_chain import (
    BASE_DOMAIN,
    PERIOD_YEARS,
    THETA0_DOMAIN,
    THETA1_DOMAIN,
    THETA_N_DOMAIN,
    target_system,
)


BOOTSTRAP_SE = {
    MEDIAN: 0.2319582116443,
    DISPERSION: 0.132477154936326,
    FAMILY_GAP: 0.563010076238571,
}
BEQUEST_NAMES = (MEDIAN, DISPERSION, FAMILY_GAP)


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


def strict_record(summary: dict[str, Any]) -> tuple[dict[str, Any], bool]:
    tight = summary.get("best_tight")
    repeat = summary.get("tight_repeat_check") or {}
    valid = bool(
        isinstance(tight, dict)
        and tight.get("strict_converged")
        and repeat.get("both_strict")
        and float(repeat.get("loss_abs_difference", 1.0)) <= 1e-10
        and float(repeat.get("max_abs_moment_difference", 1.0)) <= 1e-10
    )
    return (tight if valid else summary.get("best_search") or {}), valid


def parameter_rows(record: dict[str, Any]) -> list[dict[str, Any]]:
    theta = record["theta"]
    rows: list[dict[str, Any]] = []
    for name, lo, hi, transform in (*BASE_DOMAIN, THETA0_DOMAIN, THETA1_DOMAIN, THETA_N_DOMAIN):
        key = "beta" if name == "beta_annual" else name
        estimate = float(theta[key]) ** (1.0 / PERIOD_YEARS) if name == "beta_annual" else float(theta[key])
        rows.append(
            {
                "parameter": name,
                "estimate": estimate,
                "lower": lo,
                "upper": hi,
                "transform": transform,
                "role": "reoptimized" if name in {"theta0", "theta1", "theta_n"} else "fixed_at_M3_winner",
                "near_bound": min(estimate - lo, hi - estimate) <= 0.02 * (hi - lo),
            }
        )
    rows.append(
        {
            "parameter": "retirement_income_z_scale",
            "estimate": float(record["retirement_income_z_scale"]),
            "lower": min(RETIREMENT_INCOME_SCALE_GRID),
            "upper": max(RETIREMENT_INCOME_SCALE_GRID),
            "transform": "external diagnostic grid",
            "role": "profiled_not_calibrated",
            "near_bound": float(record["retirement_income_z_scale"]) in {
                min(RETIREMENT_INCOME_SCALE_GRID),
                max(RETIREMENT_INCOME_SCALE_GRID),
            },
        }
    )
    return rows


def main() -> None:
    args = parse_args()
    paths = sorted(args.results_dir.glob("cell_*/summary.json"))
    if len(paths) != len(RETIREMENT_INCOME_SCALE_GRID):
        raise RuntimeError(
            f"expected {len(RETIREMENT_INCOME_SCALE_GRID)} cell summaries, found {len(paths)}"
        )

    targets, weights = target_system("candidate_replacement_bequest_internal_v1")
    rows: list[dict[str, Any]] = []
    strict_candidates: list[dict[str, Any]] = []
    for path in paths:
        summary = json.loads(path.read_text())
        meta = summary["metadata"]
        if meta.get("experiment") != "retirement_income_dispersion":
            raise RuntimeError(f"wrong experiment in {path}")
        record, valid = strict_record(summary)
        moments = record.get("moments") or {}
        gaps = record.get("gaps") or {}
        within_1se = bool(
            valid
            and all(
                name in gaps and abs(float(gaps[name])) <= BOOTSTRAP_SE[name]
                for name in BEQUEST_NAMES
            )
        )
        row = {
            "cell": int(meta["cell"]),
            "retirement_income_z_scale": float(meta["retirement_income_z_scale"]),
            "eligible_tight": valid,
            "all_three_within_1se": within_1se,
            "theta0": (record.get("theta") or {}).get("theta0"),
            "theta1": (record.get("theta") or {}).get("theta1"),
            "theta_n": (record.get("theta") or {}).get("theta_n"),
            "three_target_loss": record.get("profile_loss"),
            "market_residual": record.get("market_residual"),
            "median_model": moments.get(MEDIAN),
            "median_target": targets[MEDIAN],
            "median_gap": gaps.get(MEDIAN),
            "dispersion_model": moments.get(DISPERSION),
            "dispersion_target": targets[DISPERSION],
            "dispersion_gap": gaps.get(DISPERSION),
            "family_gap_model": moments.get(FAMILY_GAP),
            "family_gap_target": targets[FAMILY_GAP],
            "family_gap_gap": gaps.get(FAMILY_GAP),
            "n_cases": summary["n_cases_completed"],
            "elapsed_sec": summary["elapsed_sec"],
        }
        rows.append(row)
        if valid:
            strict_candidates.append({"row": row, "record": record})

    rows.sort(key=lambda row: float(row["retirement_income_z_scale"]))
    if not strict_candidates:
        raise RuntimeError("no strict, exactly repeated diagnostic cell")
    selected = min(strict_candidates, key=lambda item: float(item["row"]["three_target_loss"]))
    selected_row, selected_record = selected["row"], selected["record"]
    all_moments = selected_record.get("all_moments") or {}
    missing = set(targets) - set(all_moments)
    if missing:
        raise RuntimeError(f"selected record is missing full target-system moments: {sorted(missing)}")

    fit_rows: list[dict[str, Any]] = []
    for name, target in targets.items():
        model = float(all_moments[name])
        gap = model - float(target)
        weight = float(weights[name])
        fit_rows.append(
            {
                "moment": name,
                "target": float(target),
                "model": model,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap * gap,
            }
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "retirement_income_dispersion_frontier.csv", rows)
    write_csv(args.outdir / "target_fit_full.csv", fit_rows)
    write_csv(args.outdir / "parameter_table_full.csv", parameter_rows(selected_record))
    result = {
        "status": "complete",
        "experiment_role": "diagnostic_only_no_automatic_promotion",
        "economic_object": (
            "persistent retirement-income heterogeneity linked to the existing Markov income state; "
            "mean scheduled pension unchanged"
        ),
        "strict_cells": len(strict_candidates),
        "non_strict_cells": [int(row["cell"]) for row in rows if not row["eligible_tight"]],
        "cells_matching_all_three_within_1se": sum(bool(row["all_three_within_1se"]) for row in rows),
        "selected_by_three_target_loss": selected_row,
        "frontier": rows,
    }
    (args.outdir / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True))
    lines = [
        "# Retirement-income dispersion diagnostic",
        "",
        "This is a diagnostic-only profile; it does not identify or promote a new production model.",
        f"Strict, exactly repeated cells: {len(strict_candidates)}/{len(rows)}.",
        f"Cells matching all three bequest targets within one bootstrap SE: {result['cells_matching_all_three_within_1se']}.",
        "",
        "The complete frontier, 15-moment fit, and 15-row parameter/restriction table are in this folder.",
    ]
    (args.outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
