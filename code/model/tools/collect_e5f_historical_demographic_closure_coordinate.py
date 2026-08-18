#!/usr/bin/env python3
"""Collect and verify the four-group historical coordinate panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from intergen_eqscale_seq_optimized.e5_profile import e5_target_system


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    plan_path = args.plan.resolve()
    results_dir = args.results_dir.resolve()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("status") != "complete" or int(plan.get("task_count", -1)) != 23:
        raise RuntimeError("Coordinate plan is incomplete")
    target = e5_target_system()
    expected_names = list(target.moment_names)
    summaries: list[dict[str, Any]] = []
    all_fit_rows: list[dict[str, Any]] = []
    common: dict[str, str] = {}
    for task in plan["tasks"]:
        task_id = int(task["task_id"])
        task_dir = results_dir / f"task_{task_id:03d}"
        summary_path = task_dir / "summary.json"
        fit_path = task_dir / "target_fit.csv"
        if not summary_path.is_file() or not fit_path.is_file():
            raise FileNotFoundError(f"Task {task_id} is incomplete: {task_dir}")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("status") != "complete_historical_four_group_candidate":
            raise RuntimeError(f"Task {task_id} status is invalid")
        candidate = dict(summary["candidate"])
        if (
            candidate["candidate_id"] != task["candidate_id"]
            or summary["candidate_sha256"] != task["candidate_sha256"]
        ):
            raise RuntimeError(f"Task {task_id} candidate identity failed")
        for field in (
            "selected_report_sha256",
            "source_sha256",
            "age_path_sha256",
            "scientific_bundle_sha256",
            "target_fingerprint",
        ):
            value = str(summary[field])
            if field in common and common[field] != value:
                raise RuntimeError(f"Mixed {field} across coordinate tasks")
            common[field] = value
        if summary["target_fingerprint"] != target.fingerprint:
            raise RuntimeError(f"Task {task_id} target fingerprint is stale")
        rows = read_csv(fit_path)
        if [row["moment"] for row in rows] != expected_names:
            raise RuntimeError(f"Task {task_id} target rows are incomplete or reordered")
        loss = 0.0
        for row, target_value, weight in zip(
            rows, target.target_values, target.weights, strict=True
        ):
            model_value = float(row["model"])
            gap = model_value - float(target_value)
            contribution = float(weight) * gap * gap
            if not (
                math.isclose(float(row["target"]), float(target_value), rel_tol=0.0, abs_tol=1e-13)
                and math.isclose(float(row["weight"]), float(weight), rel_tol=0.0, abs_tol=1e-10)
                and math.isclose(float(row["gap"]), gap, rel_tol=0.0, abs_tol=2e-12)
                and math.isclose(
                    float(row["loss_contribution"]),
                    contribution,
                    rel_tol=0.0,
                    abs_tol=2e-9,
                )
            ):
                raise RuntimeError(f"Task {task_id} target algebra failed for {row['moment']}")
            loss += contribution
            all_fit_rows.append({"task_id": task_id, **row})
        reported = float(summary["transition_loss_at_selected_parameters"])
        if not math.isclose(loss, reported, rel_tol=0.0, abs_tol=2e-8):
            raise RuntimeError(f"Task {task_id} loss identity failed")
        summaries.append(
            {
                "task_id": task_id,
                "candidate_id": task["candidate_id"],
                "changed_parameter": task["changed_parameter"],
                "direction": int(task["direction"]),
                "loss": reported,
                "old_psi_child": float(summary["old_psi_child"]),
                "new_psi_child": float(summary["new_psi_child"]),
                "model_solve_count": int(summary["model_solve_count"]),
                "old_stationary_solve_count": int(
                    summary["renewal_old_state"]["old_fertility_normalization"][
                        "stationary_solves"
                    ]
                ),
                "total_solve_count": int(summary["model_solve_count"])
                + int(
                    summary["renewal_old_state"]["old_fertility_normalization"][
                        "stationary_solves"
                    ]
                ),
                "elapsed_seconds": float(summary["elapsed_seconds"]),
                "summary_sha256": sha256(summary_path),
            }
        )
    anchor = next(row for row in summaries if int(row["task_id"]) == 1)
    best = min(summaries, key=lambda row: float(row["loss"]))
    for row in summaries:
        row["loss_change_from_anchor"] = float(row["loss"]) - float(anchor["loss"])
    write_csv(outdir / "coordinate_results.csv", summaries)
    write_csv(outdir / "all_target_fit.csv", all_fit_rows)
    report = {
        "status": "complete_historical_four_group_coordinate_collection",
        "plan": str(plan_path),
        "plan_sha256": sha256(plan_path),
        "results_dir": str(results_dir),
        "task_count": len(summaries),
        "anchor": anchor,
        "best": best,
        "common_contract": common,
        "total_model_solves": sum(int(row["total_solve_count"]) for row in summaries),
        "total_task_seconds": sum(float(row["elapsed_seconds"]) for row in summaries),
    }
    (outdir / "summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "HISTORICAL_FOUR_GROUP_COORDINATE_COLLECTION_COMPLETE "
        f"best={best['candidate_id']} loss={best['loss']:.9g} "
        f"anchor={anchor['loss']:.9g}",
        flush=True,
    )


if __name__ == "__main__":
    main()
