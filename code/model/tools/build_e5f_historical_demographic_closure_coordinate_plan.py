#!/usr/bin/env python3
"""Build a local coordinate plan around the selected transition estimate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import run_e5f_transition_calibration as calibration


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REPORT = (
    ROOT
    / "output/model/e5f_transition_ridge_refinement_"
    "jump11_polish_r2_20260818/report/summary.json"
)
DEFAULT_OUTDIR = ROOT / "output/model/e5f_historical_four_group_coordinate_plan"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--unit-radius", type=float, default=0.015)
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def candidate_from_units(
    *,
    center_theta: dict[str, float],
    estimated: list[dict[str, Any]],
    units: list[float],
) -> tuple[dict[str, float], float]:
    theta = {str(key): float(value) for key, value in center_theta.items()}
    delta = math.nan
    for row, unit in zip(estimated, units, strict=True):
        name = str(row["parameter"])
        value = calibration.transform_unit(
            float(unit),
            float(row["lower"]),
            float(row["upper"]),
            str(row["transform"]),
        )
        if name == "beta_annual":
            theta["beta"] = value**4
        elif name == "psi_child_change_2023":
            delta = value
        else:
            theta[name] = value
    if not math.isfinite(delta):
        raise RuntimeError("Coordinate plan omitted the terminal preference change")
    return theta, delta


def main() -> None:
    args = parse_args()
    radius = float(args.unit_radius)
    if not 0.0 < radius <= 0.05:
        raise ValueError("--unit-radius must lie in (0, 0.05]")
    report_path = args.selected_report.resolve()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    if not report_path.is_file():
        raise FileNotFoundError(report_path)
    outdir.mkdir(parents=True, exist_ok=True)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("status") != "complete_refinement_with_two_independent_identity_repeats":
        raise RuntimeError("Selected report is not complete")
    estimated = sorted(
        [dict(row) for row in report["estimated_parameters"]],
        key=lambda row: int(row["dimension"]),
    )
    if len(estimated) != 11:
        raise RuntimeError(f"Expected 11 estimated coordinates, found {len(estimated)}")
    center_units = [float(row["unit_coordinate"]) for row in estimated]
    if any(not radius < unit < 1.0 - radius for unit in center_units):
        raise RuntimeError("The declared radius leaves the parameter domain")
    center_theta = dict(report["best_candidate"]["theta"])
    tasks: list[dict[str, Any]] = []
    unit_vectors: list[list[float]] = [list(center_units)]
    changes: list[tuple[str, int]] = [("anchor", 0)]
    for dimension, row in enumerate(estimated):
        for direction in (-1, 1):
            values = list(center_units)
            values[dimension] += direction * radius
            unit_vectors.append(values)
            changes.append((str(row["parameter"]), direction))

    for task_id, (units, (changed, direction)) in enumerate(
        zip(unit_vectors, changes, strict=True), start=1
    ):
        theta, delta = candidate_from_units(
            center_theta=center_theta,
            estimated=estimated,
            units=units,
        )
        candidate_id = f"candidate_{task_id:03d}"
        payload = {
            "schema": "e5f_historical_four_group_candidate_v1",
            "candidate_id": candidate_id,
            "task_id": task_id,
            "changed_parameter": changed,
            "direction": direction,
            "unit_radius": radius,
            "unit_coordinates": units,
            "psi_child_change_2023": delta,
            "theta": theta,
            "selected_report_sha256": sha256(report_path),
            "scientific_bundle_sha256": report["scientific_contract"][
                "code_fingerprints"
            ]["bundle_sha256"],
            "target_fingerprint": report["scientific_contract"][
                "target_fingerprint"
            ],
        }
        path = outdir / f"{candidate_id}.json"
        write_json(path, payload)
        tasks.append(
            {
                "task_id": task_id,
                "candidate_id": candidate_id,
                "candidate_file": str(path),
                "candidate_sha256": sha256(path),
                "changed_parameter": changed,
                "direction": direction,
            }
        )
    plan = {
        "schema": "e5f_historical_four_group_coordinate_plan_v1",
        "status": "complete",
        "selected_report": str(report_path),
        "selected_report_sha256": sha256(report_path),
        "unit_radius": radius,
        "parameter_count": len(estimated),
        "task_count": len(tasks),
        "estimated_parameters": estimated,
        "tasks": tasks,
    }
    write_json(outdir / "coordinate_plan.json", plan)
    print(
        f"HISTORICAL_FOUR_GROUP_COORDINATE_PLAN_COMPLETE tasks={len(tasks)} "
        f"radius={radius:g}",
        flush=True,
    )


if __name__ == "__main__":
    main()
