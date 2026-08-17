#!/usr/bin/env python3
"""Collect and gate a Torch E5F dated-transition calibration panel."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--expected-tasks", type=int, required=True)
    parser.add_argument("--require-complete", action="store_true")
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def target_measurement_provenance(recorded: Any, target_set: str) -> dict[str, Any]:
    if not isinstance(recorded, dict) or not recorded:
        raise RuntimeError("task summary omits target-measurement definitions")
    return {
        "contract": str(target_set),
        "recorded_model_timing": recorded,
        "authoritative_source_ledgers": {
            "fertility_housing_tenure": "docs/model/e5_target_review_20260724.md",
            "model_statistic_definitions": "docs/model/intergen_target_object_audit.md",
            "wealth_level": "docs/model/intergen_wealth_target_beta_audit_20260723.md",
        },
        "empirical_hold": (
            "housing_increment_0to1 remains under the documented August 2026 "
            "rooms-code review; its household-and-year-FE estimate and clustered "
            "standard error remain in the active contract until an author-approved "
            "replacement is refit"
        ),
    }


def main() -> None:
    args = parse_args()
    results_dir = args.results_dir.resolve()
    report_dir = results_dir / "report"
    expected = int(args.expected_tasks)
    summaries: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    contracts: set[tuple[Any, ...]] = set()
    provenance_contracts: set[str] = set()

    for task in range(1, expected + 1):
        task_dir = results_dir / f"task_{task:03d}"
        summary_path = task_dir / "summary.json"
        if not summary_path.exists():
            failure_path = task_dir / "failure.json"
            if failure_path.exists():
                failure = read_json(failure_path)
                failures.append(
                    {
                        "task_id": task,
                        "reason": "classified_invalid_candidate",
                        "error_type": failure.get("error_type"),
                        "error": failure.get("error"),
                    }
                )
            else:
                failures.append({"task_id": task, "reason": "missing_summary"})
            continue
        try:
            summary = read_json(summary_path)
            candidate = dict(summary["best_candidate"])
            panel = dict(summary["panel_design"])
            profile = dict(summary.get("model_profile") or {})
            provenance = {
                field: summary[field]
                for field in (
                    "population_bridge",
                    "population_validation_status",
                )
            }
            provenance["target_measurements"] = target_measurement_provenance(
                summary["target_measurements"], str(summary["target_set"])
            )
            provenance_contracts.add(
                json.dumps(provenance, sort_keys=True, separators=(",", ":"))
            )
            contract = (
                str(summary["source_sha256"]),
                str(summary["target_set"]),
                str(summary["target_fingerprint"]),
                int(summary["target_count"]),
                int(panel["panel_seed"]),
                str(profile.get("name", "e5f-floor")),
                str(profile.get("profile_id", "none")),
                str(profile.get("permanent_income_log_variance", "none")),
                int(profile.get("income_state_count", 5)),
                str(profile.get("first_birth_fixed_cost_semantics", "none")),
                str(panel.get("center_sha256", "legacy_unspecified")),
                str(
                    (summary.get("code_fingerprints") or {}).get(
                        "bundle_sha256", "legacy_unspecified"
                    )
                ),
                float(summary["outside_origin_entry_share"]),
            )
            contracts.add(contract)
            if int(panel["task_id"]) != task or int(panel["panel_size"]) != expected:
                raise RuntimeError("panel task metadata mismatch")
            fit = read_csv(task_dir / "target_fit_long.csv")
            if len(fit) != int(summary["target_count"]):
                raise RuntimeError("target-fit row count mismatch")
            for row in fit:
                row["task_id"] = task
                target_rows.append(row)
            candidate["task_id"] = task
            candidate["design"] = panel["design"]
            candidate["model_profile"] = str(
                (summary.get("model_profile") or {}).get("name", "e5f-floor")
            )
            candidate["stationary_measurement_max_abs_gap"] = float(
                summary["stationary_measurement_max_abs_gap"]
            )
            candidate["valid"] = bool(
                math.isfinite(float(candidate["transition_loss"]))
                and float(candidate["max_market_residual"]) <= 2e-4
                and abs(float(candidate["population_target_gap"])) <= 2e-10
                and float(candidate["max_mass_residual"]) <= 2e-10
                and float(candidate["stationary_measurement_max_abs_gap"]) <= 2e-8
            )
            candidate_rows.append(candidate)
            summaries.append(summary)
        except Exception as error:
            failures.append(
                {"task_id": task, "reason": f"invalid_artifact: {type(error).__name__}: {error}"}
            )

    if len(contracts) > 1:
        raise RuntimeError(f"Mixed source/target/panel contracts: {sorted(contracts)}")
    if len(provenance_contracts) > 1:
        raise RuntimeError("Mixed population or target-measurement provenance")
    valid = [row for row in candidate_rows if bool(row["valid"])]
    best = (
        dict(min(valid, key=lambda row: float(row["transition_loss"])))
        if valid
        else None
    )
    if best is not None:
        selected_summary = next(
            summary
            for summary in summaries
            if int(summary["panel_design"]["task_id"]) == int(best["task_id"])
        )
        best.update(
            source=selected_summary["source"],
            source_sha256=selected_summary["source_sha256"],
            target_set=selected_summary["target_set"],
            target_fingerprint=selected_summary["target_fingerprint"],
            target_count=selected_summary["target_count"],
            model_profile=selected_summary.get("model_profile"),
            income_profile_gates=selected_summary.get("income_profile_gates"),
            code_fingerprints=selected_summary.get("code_fingerprints"),
            population_bridge=selected_summary["population_bridge"],
            population_validation_status=selected_summary[
                "population_validation_status"
            ],
            target_measurements=target_measurement_provenance(
                selected_summary["target_measurements"],
                str(selected_summary["target_set"]),
            ),
            outside_origin_entry_share=float(
                selected_summary["outside_origin_entry_share"]
            ),
            old_psi_child=selected_summary["old_psi_child"],
            old_completed_fertility=selected_summary["old_model_completed_fertility"],
        )
    candidate_rows.sort(key=lambda row: float(row["transition_loss"]))
    write_csv(report_dir / "all_candidates.csv", candidate_rows)
    write_csv(report_dir / "all_target_fit.csv", target_rows)
    if best is not None:
        write_json(report_dir / "best_candidate.json", best)
        write_json(
            report_dir / "calibration_provenance.json",
            {
                "population_bridge": best["population_bridge"],
                "population_validation_status": best[
                    "population_validation_status"
                ],
                "target_measurements": best["target_measurements"],
            },
        )
        best_fit = [
            row for row in target_rows if int(row["task_id"]) == int(best["task_id"])
        ]
        write_csv(report_dir / "best_target_fit.csv", best_fit)
    report = {
        "status": "complete" if not failures else "partial",
        "results_dir": str(results_dir),
        "expected_tasks": expected,
        "completed_tasks": len(candidate_rows),
        "valid_tasks": len(valid),
        "failed_or_missing_tasks": failures,
        "contract": list(next(iter(contracts))) if contracts else None,
        "outside_origin_entry_share": (
            float(best["outside_origin_entry_share"]) if best is not None else None
        ),
        "code_fingerprints": (
            best.get("code_fingerprints") if best is not None else None
        ),
        "model_profile": best.get("model_profile") if best is not None else None,
        "population_bridge": (
            best.get("population_bridge") if best is not None else None
        ),
        "population_validation_status": (
            best.get("population_validation_status") if best is not None else None
        ),
        "target_measurements": (
            best.get("target_measurements") if best is not None else None
        ),
        "best_candidate": best,
    }
    write_json(report_dir / "summary.json", report)
    print(
        "TRANSITION_PANEL_COLLECTED "
        f"completed={len(candidate_rows)}/{expected} valid={len(valid)} "
        f"best_loss={float(best['transition_loss']) if best else math.inf:.6f}",
        flush=True,
    )
    if args.require_complete and failures:
        raise SystemExit(2)
    if best is None:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
