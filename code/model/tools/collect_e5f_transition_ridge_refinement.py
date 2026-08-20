#!/usr/bin/env python3
"""Collect, select, and independently repeat an E5F ridge refinement.

The first pass evaluates all seven predeclared centers.  ``--validation-only``
writes an immutable selection receipt used by the repeat launcher.  A final
pass must see exactly two fresh repeat tasks and proves parameter, target,
loss, contract, and execution identity before declaring the result complete.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from build_e5f_transition_ridge_refinement import (
    CENTRAL_STEP,
    ContractError,
    ELIGIBILITY_TOLERANCES,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    canonical_bytes,
    file_sha256,
    parse_target_rows,
    read_csv,
    read_json,
    scientific_contract,
    validate_candidate_eligibility,
    validate_preference_metadata,
    validate_renewal_accounting,
)


EXPECTED_CANDIDATES = 7
EXPECTED_REPEATS = 2
LOSS_REPEAT_TOL = 1.0e-9
MOMENT_REPEAT_TOL = 1.0e-10
THETA_REPEAT_TOL = 5.0e-12
DIAGNOSTIC_REPEAT_TOL = 1.0e-10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--expected-plan-sha256", required=True)
    parser.add_argument("--report-dir", type=Path, default=None)
    parser.add_argument("--validation-only", action="store_true")
    parser.add_argument("--repeat-results-dir", type=Path, default=None)
    parser.add_argument("--expected-selection-sha256", default=None)
    return parser.parse_args()


def load_and_validate_plan(plan_path: Path, expected_sha256: str) -> tuple[dict[str, Any], str]:
    plan_path = plan_path.resolve()
    actual_sha = file_sha256(plan_path)
    if actual_sha != expected_sha256:
        raise ContractError(
            f"Plan hash mismatch: actual={actual_sha}, expected={expected_sha256}"
        )
    plan = read_json(plan_path)
    if plan.get("schema") != "e5f_transition_ridge_refinement_plan_v1":
        raise ContractError("Unsupported refinement-plan schema")
    candidates = list(plan.get("candidates") or [])
    if len(candidates) != EXPECTED_CANDIDATES:
        raise ContractError("Plan must contain exactly seven candidates")
    if [int(row.get("candidate_id", -1)) for row in candidates] != list(
        range(1, EXPECTED_CANDIDATES + 1)
    ):
        raise ContractError("Plan candidate ids must be exactly 1..7")
    if int(plan["proposal_contract"].get("formula_candidate_count", -1)) != 6:
        raise ContractError("Plan must contain six formula-generated candidates")
    if int(plan["proposal_contract"].get("observed_incumbent_count", -1)) != 1:
        raise ContractError("Plan must contain one observed incumbent")
    coordinate = dict(plan.get("coordinate_panel_contract") or {})
    dimensions = int(coordinate.get("dimensions", -1))
    if dimensions not in (10, 11):
        raise ContractError("Plan dimension must be one of the supported 10 or 11")
    if int(coordinate.get("task_count", -1)) != 1 + 2 * dimensions:
        raise ContractError("Plan task count is inconsistent with its dimension")
    central_step = float(coordinate.get("central_step", math.nan))
    if not math.isfinite(central_step) or not 0.0 < central_step <= 0.1:
        raise ContractError("Plan central-difference step must lie in (0,0.1]")
    rank = dict(plan["jacobian_diagnostics"]["full"]["relative_ranks"])
    if int(rank.get("relative_0.001", -1)) != dimensions:
        raise ContractError("Plan failed its declared full-rank gate")
    if not bool(plan["jacobian_diagnostics"]["identification_gate"].get("passed")):
        raise ContractError("Plan identification gate is not marked passed")

    plan_dir = plan_path.parent
    for row in candidates:
        candidate_path = (plan_dir / str(row["candidate_file"])).resolve()
        if not candidate_path.is_file():
            raise ContractError(f"Candidate file is missing: {candidate_path}")
        candidate_sha = file_sha256(candidate_path)
        if candidate_sha != str(row["candidate_sha256"]):
            raise ContractError(
                f"Candidate {row['candidate_id']} hash mismatch: {candidate_sha}"
            )
        payload = read_json(candidate_path)
        if payload.get("schema") != "e5f_transition_ridge_candidate_v1":
            raise ContractError(f"Candidate {row['candidate_id']} has the wrong schema")
        if int(payload.get("candidate_id", -1)) != int(row["candidate_id"]):
            raise ContractError(f"Candidate {row['candidate_id']} id mismatch")
        _compare_numeric_vector(
            payload.get("unit_vector"),
            row.get("unit_vector"),
            f"candidate {row['candidate_id']} unit vector",
            0.0,
        )
        if payload.get("scientific_contract_sha256") != plan.get(
            "scientific_contract_sha256"
        ):
            raise ContractError(f"Candidate {row['candidate_id']} contract hash mismatch")
    return plan, actual_sha


def _compare_numeric_vector(
    actual: Any,
    expected: Any,
    label: str,
    tolerance: float,
) -> None:
    try:
        left = np.asarray(actual, dtype=float)
        right = np.asarray(expected, dtype=float)
    except (TypeError, ValueError) as error:
        raise ContractError(f"{label} is not numeric") from error
    if left.shape != right.shape or not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise ContractError(f"{label} has incompatible shape or nonfinite values")
    if not np.allclose(left, right, rtol=0.0, atol=tolerance):
        raise ContractError(
            f"{label} differs by {float(np.max(np.abs(left - right)))}; "
            f"tolerance={tolerance}"
        )


def _compare_theta(
    actual: dict[str, Any], expected: dict[str, Any], label: str, tolerance: float
) -> None:
    if set(actual) != set(expected):
        raise ContractError(
            f"{label} parameter keys differ: actual={sorted(actual)}, expected={sorted(expected)}"
        )
    for name in sorted(actual):
        left = float(actual[name])
        right = float(expected[name])
        if not math.isfinite(left) or not math.isfinite(right):
            raise ContractError(f"{label}.{name} is nonfinite")
        if not math.isclose(left, right, rel_tol=0.0, abs_tol=tolerance):
            raise ContractError(
                f"{label}.{name} differs: {left} versus {right}; tolerance={tolerance}"
            )


def _candidate_path(plan_path: Path, record: dict[str, Any]) -> Path:
    return (plan_path.resolve().parent / str(record["candidate_file"])).resolve()


def require_exact_task_directories(root: Path, prefix: str, count: int) -> None:
    expected = {f"{prefix}_{task:03d}" for task in range(1, count + 1)}
    observed = {
        path.name
        for path in root.glob(f"{prefix}_*")
        if path.is_dir()
    }
    missing = sorted(expected.difference(observed))
    unexpected = sorted(observed.difference(expected))
    if missing or unexpected:
        raise ContractError(
            f"Expected exactly {count} {prefix} directories; "
            f"missing={missing}, unexpected={unexpected}"
        )


def _expected_target_contract(plan: dict[str, Any]) -> tuple[list[str], np.ndarray, np.ndarray]:
    rows = list(plan.get("target_contract") or [])
    if len(rows) != 12:
        raise ContractError("Plan target contract must contain twelve rows")
    names = [str(row["moment"]) for row in rows]
    targets = np.asarray([float(row["target"]) for row in rows], dtype=float)
    weights = np.asarray([float(row["weight"]) for row in rows], dtype=float)
    return names, targets, weights


def validate_scientific_task(
    task_dir: Path,
    plan: dict[str, Any],
    plan_path: Path,
    plan_sha256: str,
    candidate_record: dict[str, Any],
    *,
    mode: str,
    outer_task_id: int,
    replicate_id: int | None,
    selection_sha256: str | None,
) -> dict[str, Any]:
    task_dir = task_dir.resolve()
    if (task_dir / "failure.json").exists():
        raise ContractError(f"Scientific task has a failure artifact: {task_dir}")
    required_paths = (
        task_dir / "summary.json",
        task_dir / "target_fit_long.csv",
        task_dir / "refinement_task_contract.json",
    )
    missing = [str(path) for path in required_paths if not path.is_file()]
    if missing:
        raise ContractError(f"Scientific task omits required artifacts: {missing}")
    summary = read_json(task_dir / "summary.json")
    wrapper = read_json(task_dir / "refinement_task_contract.json")
    candidate_path = _candidate_path(plan_path, candidate_record)
    candidate_payload = read_json(candidate_path)
    candidate_sha = str(candidate_record["candidate_sha256"])

    wrapper_expected = {
        "schema": "e5f_transition_ridge_task_contract_v1",
        "mode": mode,
        "outer_task_id": outer_task_id,
        "candidate_id": int(candidate_record["candidate_id"]),
        "candidate_sha256": candidate_sha,
        "plan_sha256": plan_sha256,
        "scientific_contract_sha256": plan["scientific_contract_sha256"],
        "source_sha256": plan["scientific_contract"]["source_sha256"],
        "target_set": plan["scientific_contract"]["target_set"],
        "target_fingerprint": plan["scientific_contract"]["target_fingerprint"],
        "code_bundle_sha256": plan["scientific_contract"]["code_fingerprints"][
            "bundle_sha256"
        ],
        "renewal_contract_sha256": plan["renewal_contract_sha256"],
        "dated_contract_sha256": plan["dated_contract_sha256"],
    }
    for name, expected in wrapper_expected.items():
        if wrapper.get(name) != expected:
            raise ContractError(
                f"{task_dir.name}: wrapper {name}={wrapper.get(name)!r}; expected {expected!r}"
            )
    if mode == "repeat":
        if int(wrapper.get("replicate_id", -1)) != int(replicate_id or -1):
            raise ContractError(f"{task_dir.name}: repeat id mismatch")
        if wrapper.get("selection_sha256") != selection_sha256:
            raise ContractError(f"{task_dir.name}: selection hash mismatch")
    elif wrapper.get("replicate_id") is not None:
        raise ContractError(f"{task_dir.name}: validation task has a repeat id")
    execution_identity = str(wrapper.get("execution_identity") or "")
    if not execution_identity:
        raise ContractError(f"{task_dir.name}: execution identity is missing")

    if summary.get("status") != "complete_transition_calibration_panel_task":
        raise ContractError(f"{task_dir.name}: scientific summary is incomplete")
    panel = dict(summary.get("panel_design") or {})
    if int(panel.get("task_id", -1)) != 1 or int(panel.get("panel_size", -1)) != 1:
        raise ContractError(f"{task_dir.name}: driver did not run one-candidate mode")
    if panel.get("design") != "anchor" or panel.get("panel_design") != "mixed":
        raise ContractError(f"{task_dir.name}: driver did not run anchor mode")
    if int(panel.get("panel_seed", -1)) != int(
        plan["coordinate_panel_contract"]["panel_seed"]
    ):
        raise ContractError(f"{task_dir.name}: panel seed differs from plan")
    if not math.isclose(
        float(panel.get("local_radius", math.nan)),
        CENTRAL_STEP,
        rel_tol=0.0,
        abs_tol=1.0e-14,
    ):
        raise ContractError(f"{task_dir.name}: local radius is not .02")
    if panel.get("center_sha256") != candidate_sha:
        raise ContractError(f"{task_dir.name}: driver center hash is not candidate hash")
    if canonical_bytes(panel.get("domain")) != canonical_bytes(plan.get("domain")):
        raise ContractError(f"{task_dir.name}: parameter domain differs from plan")
    _compare_numeric_vector(
        panel.get("unit_vector"),
        candidate_record.get("unit_vector"),
        f"{task_dir.name} transformed unit vector",
        THETA_REPEAT_TOL,
    )

    validate_renewal_accounting(summary)
    validate_candidate_eligibility(summary)
    actual_contract = scientific_contract(summary)
    if canonical_bytes(actual_contract) != canonical_bytes(plan["scientific_contract"]):
        raise ContractError(f"{task_dir.name}: scientific contract differs from plan")
    candidate = dict(summary.get("best_candidate") or {})
    validate_preference_metadata(summary, candidate, panel)
    expected_theta = dict(candidate_payload["best_candidate"]["theta"])
    _compare_theta(
        dict(candidate.get("theta") or {}),
        expected_theta,
        f"{task_dir.name} theta",
        THETA_REPEAT_TOL,
    )
    terminal_delta = float(candidate_payload["best_candidate"]["new_psi_child"]) - float(
        candidate_payload["best_candidate"]["old_psi_child"]
    )
    realized_delta = float(panel["new_psi_child"]) - float(panel["old_psi_child"])
    if not math.isclose(
        realized_delta, terminal_delta, rel_tol=0.0, abs_tol=THETA_REPEAT_TOL
    ):
        raise ContractError(f"{task_dir.name}: terminal preference change differs from plan")

    parsed_rows, residual, table_loss = parse_target_rows(
        read_csv(task_dir / "target_fit_long.csv"),
        12,
        outer_task_id,
        str(candidate.get("candidate")),
    )
    expected_names, expected_targets, expected_weights = _expected_target_contract(plan)
    if [row["moment"] for row in parsed_rows] != expected_names:
        raise ContractError(f"{task_dir.name}: target order differs from plan")
    _compare_numeric_vector(
        [row["target"] for row in parsed_rows],
        expected_targets,
        f"{task_dir.name} target values",
        0.0,
    )
    _compare_numeric_vector(
        [row["weight"] for row in parsed_rows],
        expected_weights,
        f"{task_dir.name} target weights",
        0.0,
    )
    loss = float(candidate["transition_loss"])
    if not math.isclose(loss, table_loss, rel_tol=1.0e-9, abs_tol=1.0e-8):
        raise ContractError(f"{task_dir.name}: loss differs from its complete target table")
    return {
        "mode": mode,
        "outer_task_id": outer_task_id,
        "replicate_id": replicate_id,
        "task_dir": str(task_dir),
        "candidate_id": int(candidate_record["candidate_id"]),
        "candidate_sha256": candidate_sha,
        "execution_identity": execution_identity,
        "summary": summary,
        "candidate": candidate,
        "target_rows": parsed_rows,
        "residual": residual,
        "loss": loss,
        "wrapper": wrapper,
    }


def compare_repeat(reference: dict[str, Any], repeat: dict[str, Any]) -> None:
    label = f"repeat {repeat['replicate_id']}"
    if repeat["candidate_id"] != reference["candidate_id"]:
        raise ContractError(f"{label}: candidate id differs from selected validation")
    if repeat["candidate_sha256"] != reference["candidate_sha256"]:
        raise ContractError(f"{label}: candidate hash differs from selected validation")
    if not math.isclose(
        float(repeat["loss"]),
        float(reference["loss"]),
        rel_tol=0.0,
        abs_tol=LOSS_REPEAT_TOL,
    ):
        raise ContractError(f"{label}: loss is not reproducible")
    _compare_theta(
        dict(repeat["candidate"]["theta"]),
        dict(reference["candidate"]["theta"]),
        f"{label} theta identity",
        THETA_REPEAT_TOL,
    )
    for path in (
        ("old_psi_child",),
        ("old_model_completed_fertility",),
        ("panel_design", "old_psi_child"),
        ("panel_design", "new_psi_child"),
    ):
        left: Any = reference["summary"]
        right: Any = repeat["summary"]
        for key in path:
            left = left[key]
            right = right[key]
        if not math.isclose(
            float(left), float(right), rel_tol=0.0, abs_tol=THETA_REPEAT_TOL
        ):
            raise ContractError(f"{label}: {'.'.join(path)} differs")
    reference_renewal = dict(reference["summary"]["renewal_accounting_old_state"])
    repeat_renewal = dict(repeat["summary"]["renewal_accounting_old_state"])
    if set(reference_renewal) != set(repeat_renewal):
        raise ContractError(f"{label}: old renewal-state keys differ")
    for name in sorted(reference_renewal):
        left = reference_renewal[name]
        right = repeat_renewal[name]
        if isinstance(left, (int, float)) and isinstance(right, (int, float)):
            if not math.isclose(
                float(left), float(right), rel_tol=0.0, abs_tol=DIAGNOSTIC_REPEAT_TOL
            ):
                raise ContractError(f"{label}: renewal state {name} differs")
        elif left != right:
            raise ContractError(f"{label}: renewal state {name} differs")
    reference_rows = list(reference["target_rows"])
    repeat_rows = list(repeat["target_rows"])
    if len(reference_rows) != len(repeat_rows):
        raise ContractError(f"{label}: target table length differs")
    for left, right in zip(reference_rows, repeat_rows):
        if left["moment"] != right["moment"]:
            raise ContractError(f"{label}: target moment order differs")
        for name in ("target", "weight"):
            if float(left[name]) != float(right[name]):
                raise ContractError(f"{label}: {left['moment']} {name} differs")
        for name in ("model", "gap", "standardized_gap"):
            if not math.isclose(
                float(left[name]),
                float(right[name]),
                rel_tol=0.0,
                abs_tol=MOMENT_REPEAT_TOL,
            ):
                raise ContractError(f"{label}: {left['moment']} {name} differs")
        if not math.isclose(
            float(left["loss_contribution"]),
            float(right["loss_contribution"]),
            rel_tol=0.0,
            abs_tol=LOSS_REPEAT_TOL,
        ):
            raise ContractError(f"{label}: {left['moment']} contribution differs")
    reference_summary = reference["summary"]
    repeat_summary = repeat["summary"]
    diagnostic_paths = (
        ("stationary_measurement_max_abs_gap",),
        ("best_candidate", "max_market_residual"),
        ("best_candidate", "max_mass_residual"),
        ("best_candidate", "population_target_gap"),
        ("best_candidate", "terminal_synthetic_childless_consistency_gap"),
    )
    for path in diagnostic_paths:
        left: Any = reference_summary
        right: Any = repeat_summary
        for key in path:
            left = left[key]
            right = right[key]
        if not math.isclose(
            float(left), float(right), rel_tol=0.0, abs_tol=DIAGNOSTIC_REPEAT_TOL
        ):
            raise ContractError(f"{label}: diagnostic {'.'.join(path)} differs")


def selection_payload(
    plan_sha256: str,
    selected: dict[str, Any],
) -> dict[str, Any]:
    return {
        "schema": "e5f_transition_ridge_selection_v1",
        "status": "validation_complete_repeats_pending",
        "plan_sha256": plan_sha256,
        "selected_candidate_id": int(selected["candidate_id"]),
        "selected_candidate_sha256": str(selected["candidate_sha256"]),
        "validation_task_id": int(selected["outer_task_id"]),
        "validation_task_dir": str(selected["task_dir"]),
        "validation_execution_identity": str(selected["execution_identity"]),
        "transition_loss": float(selected["loss"]),
        "theta": selected["candidate"]["theta"],
        "target_models": [
            {"moment": row["moment"], "model": float(row["model"])}
            for row in selected["target_rows"]
        ],
    }


def write_or_validate_selection(
    path: Path,
    payload: dict[str, Any],
    expected_sha256: str | None,
) -> str:
    if path.exists():
        existing = read_json(path)
        if canonical_bytes(existing) != canonical_bytes(payload):
            raise ContractError(
                "Existing immutable selection receipt differs from current validation results"
            )
    else:
        atomic_write_json(path, payload)
    actual_sha = file_sha256(path)
    if expected_sha256 is not None and actual_sha != expected_sha256:
        raise ContractError(
            f"Selection hash mismatch: actual={actual_sha}, expected={expected_sha256}"
        )
    atomic_write_text(path.with_suffix(".sha256"), f"{actual_sha}  {path.name}\n")
    return actual_sha


def candidate_report_rows(
    validations: Sequence[dict[str, Any]], plan: dict[str, Any]
) -> list[dict[str, Any]]:
    plan_by_id = {int(row["candidate_id"]): row for row in plan["candidates"]}
    anchor_loss = float(plan["anchor"]["loss"])
    rows: list[dict[str, Any]] = []
    for result in validations:
        proposal = plan_by_id[int(result["candidate_id"])]
        predicted_loss = proposal.get("predicted_loss")
        predicted_improvement = (
            None if predicted_loss is None else anchor_loss - float(predicted_loss)
        )
        actual_improvement = anchor_loss - float(result["loss"])
        ratio = (
            actual_improvement / predicted_improvement
            if predicted_improvement is not None and predicted_improvement > 0.0
            else None
        )
        candidate = result["candidate"]
        rows.append(
            {
                "candidate_id": int(result["candidate_id"]),
                "label": proposal["label"],
                "kind": proposal["kind"],
                "candidate_sha256": result["candidate_sha256"],
                "transition_loss": float(result["loss"]),
                "anchor_loss": anchor_loss,
                "actual_improvement": actual_improvement,
                "predicted_loss": predicted_loss,
                "predicted_improvement": predicted_improvement,
                "actual_to_predicted_improvement": ratio,
                "ridge_relative": proposal.get("ridge_relative"),
                "trust_radius": proposal.get("trust_radius"),
                "step_norm": proposal.get("step_norm"),
                "max_market_residual": candidate["max_market_residual"],
                "population_target_gap": candidate["population_target_gap"],
                "max_mass_residual": candidate["max_mass_residual"],
                "stationary_measurement_max_abs_gap": result["summary"][
                    "stationary_measurement_max_abs_gap"
                ],
                "terminal_synthetic_childless_consistency_gap": candidate[
                    "terminal_synthetic_childless_consistency_gap"
                ],
                "execution_identity": result["execution_identity"],
            }
        )
    return sorted(rows, key=lambda row: float(row["transition_loss"]))


def parameter_rows(selected: dict[str, Any], plan: dict[str, Any]) -> list[dict[str, Any]]:
    theta = dict(selected["candidate"]["theta"])
    unit = np.asarray(
        next(
            row["unit_vector"]
            for row in plan["candidates"]
            if int(row["candidate_id"]) == int(selected["candidate_id"])
        ),
        dtype=float,
    )
    rows: list[dict[str, Any]] = []
    for index, domain in enumerate(plan["domain"]):
        name = str(domain["name"])
        if name == "beta_annual":
            value = float(theta["beta"]) ** 0.25
        elif name == "psi_child_change_2023":
            panel = selected["summary"]["panel_design"]
            value = float(panel["new_psi_child"]) - float(panel["old_psi_child"])
        else:
            value = float(theta[name])
        rows.append(
            {
                "dimension": index,
                "parameter": name,
                "estimate": value,
                "lower": float(domain["lower"]),
                "upper": float(domain["upper"]),
                "transform": str(domain["transform"]),
                "unit_coordinate": float(unit[index]),
                "distance_to_nearest_unit_bound": float(min(unit[index], 1.0 - unit[index])),
                "near_bound_within_0.02": bool(min(unit[index], 1.0 - unit[index]) <= 0.02),
                "status": "estimated_in_dated_transition_refinement",
            }
        )
    return rows


def copy_selected_artifacts(selected: dict[str, Any], report_dir: Path) -> None:
    task_dir = Path(selected["task_dir"])
    case_name = str(selected["candidate"]["candidate"])
    case_dir = task_dir / "cases" / case_name
    artifacts = {
        "transition_path.csv": "selected_transition_path.csv",
        "dated_model_moments.csv": "selected_dated_model_moments.csv",
        "dated_period_fertility.csv": "selected_dated_period_fertility.csv",
        "dated_period_fertility_by_age.csv": "selected_dated_period_fertility_by_age.csv",
        "dated_first_birth_housing_did.csv": "selected_dated_first_birth_housing_did.csv",
        "dated_cohort_timing_ledgers.json": "selected_dated_cohort_timing_ledgers.json",
    }
    for source_name, destination_name in artifacts.items():
        source = case_dir / source_name
        if not source.is_file():
            raise ContractError(f"Selected validation omits required artifact: {source}")
        shutil.copy2(source, report_dir / destination_name)
    parameter_table = task_dir / "parameter_table.csv"
    if not parameter_table.is_file():
        raise ContractError("Selected validation omits parameter_table.csv")
    shutil.copy2(parameter_table, report_dir / "selected_driver_parameter_table.csv")


def main() -> None:
    args = parse_args()
    if args.validation_only:
        if args.repeat_results_dir is not None or args.expected_selection_sha256 is not None:
            raise ContractError("Validation-only collection cannot consume repeat arguments")
    elif args.repeat_results_dir is None or args.expected_selection_sha256 is None:
        raise ContractError(
            "Final collection requires --repeat-results-dir and --expected-selection-sha256"
        )

    plan_path = args.plan.resolve()
    plan, plan_sha = load_and_validate_plan(plan_path, str(args.expected_plan_sha256))
    results_dir = args.results_dir.resolve()
    require_exact_task_directories(results_dir, "task", EXPECTED_CANDIDATES)
    report_dir = (
        args.report_dir.resolve()
        if args.report_dir is not None
        else results_dir / "report"
    )
    report_dir.mkdir(parents=True, exist_ok=True)
    candidates = list(plan["candidates"])

    validations: list[dict[str, Any]] = []
    for task_id, candidate_record in enumerate(candidates, start=1):
        task_dir = results_dir / f"task_{task_id:03d}"
        validations.append(
            validate_scientific_task(
                task_dir,
                plan,
                plan_path,
                plan_sha,
                candidate_record,
                mode="validation",
                outer_task_id=task_id,
                replicate_id=None,
                selection_sha256=None,
            )
        )
    selected = min(validations, key=lambda row: float(row["loss"]))
    selection = selection_payload(plan_sha, selected)
    selection_path = report_dir / "selection.json"
    selection_sha = write_or_validate_selection(
        selection_path,
        selection,
        None if args.validation_only else str(args.expected_selection_sha256),
    )

    all_candidate_rows = candidate_report_rows(validations, plan)
    all_target_rows: list[dict[str, Any]] = []
    for validation in validations:
        for row in validation["target_rows"]:
            all_target_rows.append(
                {
                    "candidate_id": validation["candidate_id"],
                    **{key: value for key, value in row.items() if key != "index"},
                }
            )
    atomic_write_csv(report_dir / "all_candidates.csv", all_candidate_rows)
    atomic_write_csv(report_dir / "all_target_fit.csv", all_target_rows)
    atomic_write_csv(
        report_dir / "selected_target_fit.csv",
        [{key: value for key, value in row.items() if key != "index"} for row in selected["target_rows"]],
    )
    atomic_write_csv(report_dir / "selected_parameters_and_bounds.csv", parameter_rows(selected, plan))
    copy_selected_artifacts(selected, report_dir)

    repeats: list[dict[str, Any]] = []
    repeat_status: dict[str, Any]
    if args.validation_only:
        repeat_status = {
            "required": EXPECTED_REPEATS,
            "completed": 0,
            "identity_gate_passed": False,
            "numerical_identity_gate_passed": False,
        }
        status = "validation_complete_repeats_pending"
        promotion_eligible = False
    else:
        repeat_root = args.repeat_results_dir.resolve()
        require_exact_task_directories(repeat_root, "repeat", EXPECTED_REPEATS)
        selected_record = next(
            row
            for row in candidates
            if int(row["candidate_id"]) == int(selected["candidate_id"])
        )
        for replicate_id in range(1, EXPECTED_REPEATS + 1):
            repeats.append(
                validate_scientific_task(
                    repeat_root / f"repeat_{replicate_id:03d}",
                    plan,
                    plan_path,
                    plan_sha,
                    selected_record,
                    mode="repeat",
                    outer_task_id=replicate_id,
                    replicate_id=replicate_id,
                    selection_sha256=selection_sha,
                )
            )
        identities = [str(row["execution_identity"]) for row in repeats]
        identities.append(str(selected["execution_identity"]))
        if len(set(identities)) != len(identities):
            raise ContractError(
                "The selected validation and two repeats do not have three distinct execution identities"
            )
        for repeat in repeats:
            compare_repeat(selected, repeat)
        compare_repeat(repeats[0], repeats[1])
        repeat_status = {
            "required": EXPECTED_REPEATS,
            "completed": EXPECTED_REPEATS,
            "execution_identities": [row["execution_identity"] for row in repeats],
            "identity_gate_passed": True,
            "numerical_identity_gate_passed": True,
            "loss_abs_tolerance": LOSS_REPEAT_TOL,
            "model_moment_abs_tolerance": MOMENT_REPEAT_TOL,
            "theta_abs_tolerance": THETA_REPEAT_TOL,
        }
        status = "complete_refinement_with_two_independent_identity_repeats"
        promotion_eligible = True

    summary = {
        "schema": "e5f_transition_ridge_refinement_report_v1",
        "status": status,
        "promotion_eligible": promotion_eligible,
        "collector": str(Path(__file__).resolve()),
        "collector_sha256": file_sha256(Path(__file__).resolve()),
        "results_dir": str(results_dir),
        "plan": str(plan_path),
        "plan_sha256": plan_sha,
        "selection": str(selection_path),
        "selection_sha256": selection_sha,
        "scientific_contract": plan["scientific_contract"],
        "scientific_contract_sha256": plan["scientific_contract_sha256"],
        "renewal_contract_sha256": plan["renewal_contract_sha256"],
        "dated_contract_sha256": plan["dated_contract_sha256"],
        "coordinate_jacobian_diagnostics": plan["jacobian_diagnostics"],
        "validation_candidate_count": len(validations),
        "all_validation_candidates_eligible": True,
        "selected_candidate_id": int(selected["candidate_id"]),
        "selected_candidate_sha256": selected["candidate_sha256"],
        "selected_transition_loss": float(selected["loss"]),
        "selected_validation_execution_identity": selected["execution_identity"],
        "repeat_gate": repeat_status,
        "eligibility_tolerances": ELIGIBILITY_TOLERANCES,
        "target_fit_table": [
            {key: value for key, value in row.items() if key != "index"}
            for row in selected["target_rows"]
        ],
        "estimated_parameters": parameter_rows(selected, plan),
        "best_candidate": {
            **selected["candidate"],
            "candidate_id": int(selected["candidate_id"]),
            "candidate_sha256": selected["candidate_sha256"],
            "old_psi_child": selected["summary"]["old_psi_child"],
            "old_completed_fertility": selected["summary"][
                "old_model_completed_fertility"
            ],
            "renewal_accounting_contract": selected["summary"][
                "renewal_accounting_contract"
            ],
            "renewal_accounting_old_state": selected["summary"][
                "renewal_accounting_old_state"
            ],
            "dated_measurement_contract": selected["summary"][
                "dated_measurement_contract"
            ],
        },
    }
    summary_path = report_dir / "summary.json"
    atomic_write_json(summary_path, summary)
    summary_sha = file_sha256(summary_path)
    atomic_write_text(report_dir / "summary.sha256", f"{summary_sha}  summary.json\n")
    print(
        "E5F_RIDGE_REFINEMENT_COLLECTED "
        f"status={status} selected={selected['candidate_id']} "
        f"loss={selected['loss']:.9f} selection_sha256={selection_sha}",
        flush=True,
    )


if __name__ == "__main__":
    try:
        main()
    except ContractError as error:
        print(f"E5F_RIDGE_REFINEMENT_COLLECTION_REJECTED {error}", file=sys.stderr, flush=True)
        raise SystemExit(2) from error
