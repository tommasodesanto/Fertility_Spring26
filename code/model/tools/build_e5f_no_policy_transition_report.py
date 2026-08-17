#!/usr/bin/env python3
"""Build the canonical no-policy E5F transition report.

This is a deliberately narrow presentation layer.  It accepts only a fully
collected ridge-refinement result with two independent repeats, the separately
validated 2007--2023 historical packet, and the paired closed/open post-2023
continuation packet.  It does not solve, refit, rescale, or run policy cases.
Every input contract is checked before a new output directory is created.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


EXPECTED_TARGETS = 12
EXPECTED_PARAMETERS = 10
EXPECTED_YEARS = (2007, 2011, 2015, 2019, 2023)
LOSS_TOL = 1.0e-9
MOMENT_TOL = 1.0e-10
PARAMETER_TOL = 5.0e-12
DIAGNOSTIC_TOL = 1.0e-10
TARGET_PROVENANCE_STATUSES = {
    "measured",
    "provisional",
    "external_normalization",
    "empirical_hold",
}

HISTORICAL_CLASSIFICATIONS = {
    "imposed_bridge_matched_by_construction_not_fit": "imposed bridge",
    "imposed_bridge_numerical_audit_not_fit": "imposed bridge",
    "terminal_maintained_target": "calibrated 2023 moment",
    "descriptive_noncomparable_holdout": "untargeted holdout",
    "untargeted_holdout": "untargeted holdout",
    "untargeted_holdout_data_only": "untargeted holdout",
    "untargeted_holdout_cross_section_only_coding_break": "untargeted holdout",
    "untargeted_holdout_longitudinal_2011_2023": "untargeted holdout",
    "descriptive_noncomparable_diagnostic": "untargeted holdout",
    "numerical_residual_audit_not_fit": "numerical residual",
}

MOMENT_LABELS = {
    "tfr": "Completed fertility, age-42 cohort",
    "childless_rate": "Childless share, age-42 cohort",
    "mean_age_first_birth": "Mean age at first birth",
    "share_first_births_age30plus": "Share of first births at age 30+",
    "housing_increment_0to1": "Rooms response around first birth",
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": "Rooms gap: 3+ versus 1--2 children",
    "own_family_gap": "Ownership gap: family versus no family",
    "own_rate": "Ownership rate",
    "aggregate_mean_occupied_rooms_18_85": "Mean occupied rooms, ages 18--85",
    "aggregate_wealth_to_annual_gross_labor_earnings": "Wealth / annual gross labor earnings",
    "annual_bequest_flow_to_aggregate_wealth": "Annual bequests / aggregate wealth",
    "old_total_wealth_to_annual_income_p90_p50_7684": "Old-age wealth p90/p50",
}

PARAMETER_LABELS = {
    "beta_annual": "Annual discount factor",
    "kappa_fert": "First-birth logit scale",
    "kappa_fert_continuation": "Continuation-birth logit scale",
    "chi": "Housing preference weight",
    "H0": "Housing-supply scale",
    "theta0": "Bequest intercept",
    "theta1": "Bequest slope",
    "hbar_child_rooms": "Per-child room floor",
    "first_birth_fixed_cost": "One-time first-birth utility cost",
    "psi_child_change_2023": "2007--2023 fertility-preference change",
}


class ContractError(RuntimeError):
    """Raised when an input is incomplete, mixed, or outside report scope."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--selected-report",
        type=Path,
        required=True,
        help="Final ridge report summary.json or its directory.",
    )
    parser.add_argument(
        "--repeat-task",
        type=Path,
        action="append",
        default=[],
        help="One exact-repeat task directory (must be supplied exactly twice).",
    )
    parser.add_argument(
        "--refinement-plan",
        type=Path,
        required=True,
        help="The candidate_plan.json used by the final ridge collector.",
    )
    parser.add_argument(
        "--target-provenance-csv",
        type=Path,
        required=True,
        help="Verified twelve-row empirical target/provenance ledger.",
    )
    parser.add_argument(
        "--historical-packet",
        type=Path,
        required=True,
        help="Completed output of build_e5f_transition_historical_validation.py.",
    )
    parser.add_argument(
        "--continuation-packet",
        type=Path,
        required=True,
        help="Completed paired no-policy closed/open continuation packet.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode()


def object_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ContractError(f"Cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ContractError(f"JSON root must be an object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
    except OSError as error:
        raise ContractError(f"Cannot read CSV {path}: {error}") from error
    if not rows:
        raise ContractError(f"CSV is empty: {path}")
    return rows


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ContractError(f"Refusing to write empty table: {path.name}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def finite_float(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ContractError(f"{label} is not numeric: {value!r}") from error
    if not math.isfinite(number):
        raise ContractError(f"{label} is not finite")
    return number


def integer(value: Any, label: str) -> int:
    number = finite_float(value, label)
    rounded = int(round(number))
    if not math.isclose(number, rounded, rel_tol=0.0, abs_tol=1.0e-12):
        raise ContractError(f"{label} is not an integer: {number}")
    return rounded


def require_hash(value: Any, label: str) -> str:
    text = str(value or "")
    if len(text) != 64 or any(character not in "0123456789abcdef" for character in text):
        raise ContractError(f"{label} is not a lowercase SHA-256 digest")
    return text


def require_columns(
    rows: Sequence[Mapping[str, Any]], fields: Iterable[str], label: str
) -> None:
    missing = sorted(set(fields).difference(rows[0]))
    if missing:
        raise ContractError(f"{label} omits columns {missing}")


def assert_close(
    actual: Any, expected: Any, label: str, tolerance: float
) -> None:
    left = finite_float(actual, f"{label} actual")
    right = finite_float(expected, f"{label} expected")
    if not math.isclose(left, right, rel_tol=0.0, abs_tol=tolerance):
        raise ContractError(
            f"{label} differs: actual={left:.17g}, expected={right:.17g}, "
            f"tolerance={tolerance:g}"
        )


def assert_same_object(actual: Any, expected: Any, label: str) -> None:
    if canonical_bytes(actual) != canonical_bytes(expected):
        raise ContractError(f"{label} differs across inputs")


def resolve_summary(path: Path, label: str) -> tuple[Path, Path]:
    path = path.resolve()
    summary = path / "summary.json" if path.is_dir() else path
    if not summary.is_file() or summary.name != "summary.json":
        raise ContractError(f"{label} must resolve to summary.json: {path}")
    return summary, summary.parent


def resolve_repeat_task(path: Path) -> tuple[Path, Path]:
    path = path.resolve()
    if path.is_file():
        if path.name != "summary.json":
            raise ContractError(f"Repeat file is not summary.json: {path}")
        return path, path.parent
    if not path.is_dir():
        raise ContractError(f"Repeat task does not exist: {path}")
    direct = path / "summary.json"
    if direct.is_file():
        return direct, path
    matches = sorted(path.glob("repeat_*/summary.json"))
    if len(matches) != 1:
        raise ContractError(
            f"Repeat input must identify exactly one task; found {len(matches)} in {path}"
        )
    return matches[0].resolve(), matches[0].resolve().parent


def receipt(path: Path) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file():
        raise ContractError(f"Missing input artifact: {path}")
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }


def validate_sha_sidecar(path: Path) -> Path:
    sidecar = path.with_suffix(".sha256")
    if not sidecar.is_file():
        raise ContractError(f"Missing SHA-256 sidecar: {sidecar}")
    expected = f"{file_sha256(path)}  {path.name}\n"
    if sidecar.read_text() != expected:
        raise ContractError(f"SHA-256 sidecar does not reproduce: {sidecar}")
    return sidecar


def task_scientific_contract(summary: Mapping[str, Any]) -> dict[str, Any]:
    fields = (
        "source",
        "source_sha256",
        "source_metadata",
        "target_set",
        "target_fingerprint",
        "target_count",
        "code_fingerprints",
        "model_profile",
        "income_profile_gates",
        "target_measurements",
        "renewal_accounting_contract",
        "dated_measurement_contract",
        "population_validation_status",
        "outside_origin_entry_share",
        "old_completed_fertility_reference",
        "policy_case",
        "post_2023_periods",
        "stationary_free_parameter_count",
        "transition_free_parameter_count",
    )
    missing = [field for field in fields if field not in summary]
    if missing:
        raise ContractError(f"Scientific task summary omits contract fields {missing}")
    bridge = summary.get("population_bridge")
    if not isinstance(bridge, dict) or not isinstance(
        bridge.get("initial_age_reweight"), dict
    ):
        raise ContractError("Scientific task omits population-bridge provenance")
    initial = bridge["initial_age_reweight"]
    bridge_fields = ("status", "target_indices", "total_households", "age_distribution")
    initial_fields = (
        "mapping",
        "model_age_grid",
        "sample",
        "source",
        "source_path",
        "source_receipt",
        "source_sha256",
        "source_year",
        "units",
        "period_years",
        "horizon_periods",
        "horizon_calendar_year",
        "future_entry_rule",
    )
    missing_bridge = [field for field in bridge_fields if field not in bridge]
    missing_bridge += [
        f"initial_age_reweight.{field}" for field in initial_fields if field not in initial
    ]
    if missing_bridge:
        raise ContractError(f"Population bridge omits fields {missing_bridge}")
    projected_bridge = {
        **{field: bridge[field] for field in bridge_fields},
        "initial_age_reweight": {field: initial[field] for field in initial_fields},
    }
    return {
        "source": str(summary["source"]),
        "source_sha256": str(summary["source_sha256"]),
        "source_metadata": summary["source_metadata"],
        "target_set": str(summary["target_set"]),
        "target_fingerprint": str(summary["target_fingerprint"]),
        "target_count": integer(summary["target_count"], "task target_count"),
        "code_fingerprints": summary["code_fingerprints"],
        "model_profile": summary["model_profile"],
        "income_profile_gates": summary["income_profile_gates"],
        "target_measurements": summary["target_measurements"],
        "renewal_accounting_contract": summary["renewal_accounting_contract"],
        "dated_measurement_contract": summary["dated_measurement_contract"],
        "population_bridge": projected_bridge,
        "population_validation_status": str(summary["population_validation_status"]),
        "outside_origin_entry_share": finite_float(
            summary["outside_origin_entry_share"], "outside-origin entry share"
        ),
        "old_completed_fertility_reference": finite_float(
            summary["old_completed_fertility_reference"], "old fertility reference"
        ),
        "policy_case": str(summary["policy_case"]),
        "post_2023_periods": integer(summary["post_2023_periods"], "post_2023_periods"),
        "stationary_free_parameter_count": integer(
            summary["stationary_free_parameter_count"], "stationary free parameters"
        ),
        "transition_free_parameter_count": integer(
            summary["transition_free_parameter_count"], "transition free parameters"
        ),
    }


def parse_target_rows(rows: Sequence[Mapping[str, Any]], label: str) -> list[dict[str, Any]]:
    if len(rows) != EXPECTED_TARGETS:
        raise ContractError(f"{label} must contain exactly twelve rows")
    require_columns(
        rows,
        ("moment", "target", "model", "gap", "weight", "loss_contribution"),
        label,
    )
    parsed: list[dict[str, Any]] = []
    names: set[str] = set()
    for index, row in enumerate(rows):
        name = str(row["moment"])
        if not name or name in names:
            raise ContractError(f"{label} has an empty or duplicate moment {name!r}")
        names.add(name)
        target = finite_float(row["target"], f"{label}.{name}.target")
        model = finite_float(row["model"], f"{label}.{name}.model")
        gap = finite_float(row["gap"], f"{label}.{name}.gap")
        weight = finite_float(row["weight"], f"{label}.{name}.weight")
        contribution = finite_float(
            row["loss_contribution"], f"{label}.{name}.loss_contribution"
        )
        if weight <= 0.0:
            raise ContractError(f"{label}.{name}.weight must be positive")
        assert_close(gap, model - target, f"{label}.{name}.gap identity", 1.0e-11)
        assert_close(
            contribution,
            weight * gap * gap,
            f"{label}.{name}.loss identity",
            1.0e-8,
        )
        parsed.append(
            {
                "order": index + 1,
                "moment": name,
                "label": MOMENT_LABELS.get(name, name.replace("_", " ")),
                "target": target,
                "model": model,
                "gap": gap,
                "weight": weight,
                "loss_contribution": contribution,
                "standardized_gap": math.copysign(math.sqrt(contribution), gap),
            }
        )
    return parsed


def parse_parameter_rows(
    rows: Sequence[Mapping[str, Any]], plan: Mapping[str, Any]
) -> list[dict[str, Any]]:
    if len(rows) != EXPECTED_PARAMETERS:
        raise ContractError("Selected parameter table must contain ten estimated rows")
    require_columns(
        rows,
        (
            "dimension",
            "parameter",
            "estimate",
            "lower",
            "upper",
            "transform",
            "unit_coordinate",
            "distance_to_nearest_unit_bound",
            "near_bound_within_0.02",
            "status",
        ),
        "selected parameter table",
    )
    domain = list(plan.get("domain") or [])
    if len(domain) != EXPECTED_PARAMETERS:
        raise ContractError("Refinement plan must have a ten-dimensional domain")
    parsed: list[dict[str, Any]] = []
    for index, (row, restriction) in enumerate(zip(rows, domain)):
        if integer(row["dimension"], "parameter dimension") != index:
            raise ContractError("Parameter dimensions are not ordered 0..9")
        name = str(row["parameter"])
        if name != str(restriction.get("name")):
            raise ContractError(f"Parameter {index} differs from the plan domain")
        lower = finite_float(row["lower"], f"{name}.lower")
        upper = finite_float(row["upper"], f"{name}.upper")
        estimate = finite_float(row["estimate"], f"{name}.estimate")
        unit = finite_float(row["unit_coordinate"], f"{name}.unit_coordinate")
        distance = finite_float(
            row["distance_to_nearest_unit_bound"], f"{name}.bound distance"
        )
        assert_close(lower, restriction["lower"], f"{name}.lower", 0.0)
        assert_close(upper, restriction["upper"], f"{name}.upper", 0.0)
        if str(row["transform"]) != str(restriction["transform"]):
            raise ContractError(f"{name}.transform differs from the plan")
        if not lower <= estimate <= upper or not 0.0 <= unit <= 1.0:
            raise ContractError(f"{name} is outside its reported bounds")
        assert_close(distance, min(unit, 1.0 - unit), f"{name}.bound distance", 1e-12)
        near = distance <= 0.02
        recorded_near = str(row["near_bound_within_0.02"]).strip().lower() in {
            "true",
            "1",
        }
        if near != recorded_near:
            raise ContractError(f"{name}.near-bound flag is inconsistent")
        if row["status"] != "estimated_in_dated_transition_refinement":
            raise ContractError(f"{name} is not classified as estimated")
        parsed.append(
            {
                "dimension": index,
                "parameter": name,
                "label": PARAMETER_LABELS.get(name, name.replace("_", " ")),
                "estimate": estimate,
                "lower": lower,
                "upper": upper,
                "transform": str(row["transform"]),
                "unit_coordinate": unit,
                "distance_to_nearest_unit_bound": distance,
                "near_bound_within_0.02": near,
                "status": "estimated/free in dated transition",
            }
        )
    return parsed


def validate_target_provenance(
    path: Path, selected: Mapping[str, Any]
) -> dict[str, Any]:
    path = path.resolve()
    rows = read_csv(path)
    if len(rows) != EXPECTED_TARGETS:
        raise ContractError("Target-provenance ledger must contain exactly twelve rows")
    fields = (
        "moment",
        "target",
        "authoritative_builder",
        "sample_geography_vintage",
        "estimator_measurement",
        "fixed_effects",
        "clustering",
        "uncertainty_type",
        "standard_error_or_scale",
        "weight",
        "status",
        "caveat",
    )
    require_columns(rows, fields, "target-provenance ledger")
    parsed: list[dict[str, Any]] = []
    for index, (row, target_row) in enumerate(zip(rows, selected["target_rows"])):
        moment = str(row["moment"]).strip()
        if moment != target_row["moment"]:
            raise ContractError(
                f"Target-provenance row {index + 1} is {moment!r}; "
                f"expected {target_row['moment']!r}"
            )
        assert_close(
            row["target"],
            target_row["target"],
            f"target provenance {moment} target",
            0.0,
        )
        assert_close(
            row["weight"],
            target_row["weight"],
            f"target provenance {moment} weight",
            0.0,
        )
        text = {field: str(row[field]).strip() for field in fields if field not in {"target", "weight", "standard_error_or_scale"}}
        missing = [field for field, value in text.items() if not value]
        if missing:
            raise ContractError(f"Target provenance {moment} has empty fields {missing}")
        status = text["status"]
        if status not in TARGET_PROVENANCE_STATUSES:
            raise ContractError(
                f"Target provenance {moment} has unsupported status {status!r}"
            )
        scale = finite_float(
            row["standard_error_or_scale"], f"target provenance {moment} scale"
        )
        if scale <= 0.0:
            raise ContractError(f"Target provenance {moment} scale must be positive")
        weight = finite_float(row["weight"], f"target provenance {moment} weight")
        if not math.isclose(weight, 1.0 / (scale * scale), rel_tol=2.0e-8, abs_tol=1.0e-9):
            raise ContractError(
                f"Target provenance {moment} weight is not inverse scale squared"
            )
        parsed.append(
            {
                "order": index + 1,
                "moment": moment,
                "label": target_row["label"],
                "target": target_row["target"],
                "weight": target_row["weight"],
                "authoritative_builder": text["authoritative_builder"],
                "sample_geography_vintage": text["sample_geography_vintage"],
                "estimator_measurement": text["estimator_measurement"],
                "fixed_effects": text["fixed_effects"],
                "clustering": text["clustering"],
                "uncertainty_type": text["uncertainty_type"],
                "standard_error_or_scale": scale,
                "status": status,
                "caveat": text["caveat"],
            }
        )
    counts = Counter(row["status"] for row in parsed)
    return {
        "path": path,
        "rows": parsed,
        "status_counts": dict(sorted(counts.items())),
    }


def validate_jacobian(
    plan_path: Path, plan: Mapping[str, Any]
) -> dict[str, Any]:
    plan_dir = plan_path.parent
    paths = {
        name: plan_dir / name
        for name in (
            "weighted_jacobian.csv",
            "singular_values.csv",
            "side_consistency.csv",
            "weak_directions.csv",
        )
    }
    for path in paths.values():
        if not path.is_file():
            raise ContractError(f"Refinement diagnostic is missing: {path}")

    target_contract = list(plan.get("target_contract") or [])
    domain = list(plan.get("domain") or [])
    if len(target_contract) != EXPECTED_TARGETS or len(domain) != EXPECTED_PARAMETERS:
        raise ContractError("Refinement plan has the wrong target/domain dimensions")
    rows = read_csv(paths["weighted_jacobian.csv"])
    if len(rows) != EXPECTED_TARGETS:
        raise ContractError("Weighted Jacobian must have twelve rows")
    parameter_names = [str(row["name"]) for row in domain]
    require_columns(rows, ("moment_index", "moment", *parameter_names), "weighted Jacobian")
    matrix = np.empty((EXPECTED_TARGETS, EXPECTED_PARAMETERS), dtype=float)
    for index, (row, target) in enumerate(zip(rows, target_contract)):
        if integer(row["moment_index"], "Jacobian moment_index") != index:
            raise ContractError("Weighted Jacobian moment indices are not ordered")
        if row["moment"] != target["moment"]:
            raise ContractError("Weighted Jacobian moment order differs from target contract")
        for dimension, name in enumerate(parameter_names):
            matrix[index, dimension] = finite_float(row[name], f"Jacobian[{index},{name}]")

    singular = np.linalg.svd(matrix, compute_uv=False)
    diagnostics = dict(plan.get("jacobian_diagnostics") or {})
    full = dict(diagnostics.get("full") or {})
    recorded = np.asarray(full.get("singular_values") or [], dtype=float)
    if recorded.shape != singular.shape or not np.allclose(
        singular, recorded, rtol=1.0e-10, atol=1.0e-10
    ):
        raise ContractError("Recomputed Jacobian singular values differ from the plan")
    leading = float(singular[0])
    rank_1e3 = int(np.sum(singular >= 1.0e-3 * leading)) if leading > 0 else 0
    reported_ranks = dict(full.get("relative_ranks") or {})
    if rank_1e3 != EXPECTED_PARAMETERS or integer(
        reported_ranks.get("relative_0.001"), "reported rank at 1e-3"
    ) != EXPECTED_PARAMETERS:
        raise ContractError("Weighted Jacobian fails the ten-dimensional rank gate")
    gate = dict(diagnostics.get("identification_gate") or {})
    if gate.get("passed") is not True:
        raise ContractError("Refinement identification gate is not marked passed")

    side_rows = read_csv(paths["side_consistency.csv"])
    if len(side_rows) != EXPECTED_PARAMETERS:
        raise ContractError("Side-consistency table must have ten rows")
    require_columns(
        side_rows,
        ("dimension", "parameter", "side_consistency_cosine", "frozen_for_step"),
        "side consistency",
    )
    frozen: list[int] = []
    for index, (row, name) in enumerate(zip(side_rows, parameter_names)):
        if integer(row["dimension"], "side dimension") != index or row["parameter"] != name:
            raise ContractError("Side-consistency dimensions differ from the plan domain")
        finite_float(row["side_consistency_cosine"], f"{name}.side cosine")
        if str(row["frozen_for_step"]).lower() in {"true", "1"}:
            frozen.append(index)
    if frozen != [int(value) for value in diagnostics.get("frozen_dimensions", [])]:
        raise ContractError("Side-consistency frozen dimensions differ from the plan")

    singular_rows = read_csv(paths["singular_values.csv"])
    if len(singular_rows) != EXPECTED_PARAMETERS:
        raise ContractError("Singular-value table must have ten rows")
    for index, row in enumerate(singular_rows):
        assert_close(
            row["singular_value"], singular[index], f"singular value {index}", 1.0e-10
        )
    weak_rows = read_csv(paths["weak_directions.csv"])
    require_columns(
        weak_rows,
        ("weak_order", "singular_value_index", "parameter", "loading", "squared_loading_share"),
        "weak directions",
    )
    for row in weak_rows:
        finite_float(row["loading"], "weak-direction loading")
        finite_float(row["squared_loading_share"], "weak-direction loading share")

    return {
        "paths": paths,
        "matrix": matrix,
        "singular_values": singular,
        "rank_at_relative_1e3": rank_1e3,
        "condition_number": finite_float(full.get("condition_number"), "condition number"),
        "frozen_dimensions": frozen,
        "parameter_names": parameter_names,
        "moment_names": [str(row["moment"]) for row in target_contract],
    }


def validate_selected_report(
    selected_path: Path, plan_path: Path
) -> dict[str, Any]:
    summary_path, report_dir = resolve_summary(selected_path, "Selected report")
    summary = read_json(summary_path)
    if summary.get("schema") != "e5f_transition_ridge_refinement_report_v1":
        raise ContractError("Selected report has an unsupported schema")
    if summary.get("status") != "complete_refinement_with_two_independent_identity_repeats":
        raise ContractError("Selected ridge refinement is not complete with two repeats")
    if summary.get("promotion_eligible") is not True:
        raise ContractError("Selected ridge refinement is not promotion eligible")
    repeat_gate = dict(summary.get("repeat_gate") or {})
    if any(
        (
            integer(repeat_gate.get("required"), "required repeats") != 2,
            integer(repeat_gate.get("completed"), "completed repeats") != 2,
            repeat_gate.get("identity_gate_passed") is not True,
            repeat_gate.get("numerical_identity_gate_passed") is not True,
        )
    ):
        raise ContractError("Selected report's two-repeat gate did not pass")

    plan_path = plan_path.resolve()
    if not plan_path.is_file() or plan_path.name != "candidate_plan.json":
        raise ContractError("--refinement-plan must be candidate_plan.json")
    plan_sha = file_sha256(plan_path)
    if require_hash(summary.get("plan_sha256"), "selected plan_sha256") != plan_sha:
        raise ContractError("Selected report and refinement plan hashes differ")
    plan_sha_path = validate_sha_sidecar(plan_path)
    plan = read_json(plan_path)
    if plan.get("schema") != "e5f_transition_ridge_refinement_plan_v1":
        raise ContractError("Unsupported refinement-plan schema")
    assert_same_object(summary.get("scientific_contract"), plan.get("scientific_contract"), "scientific contract")
    if require_hash(summary.get("scientific_contract_sha256"), "scientific contract hash") != object_sha256(plan["scientific_contract"]):
        raise ContractError("Scientific-contract hash does not reproduce")
    for summary_key, plan_key, payload_key in (
        ("renewal_contract_sha256", "renewal_contract_sha256", "renewal_accounting_contract"),
        ("dated_contract_sha256", "dated_contract_sha256", "dated_measurement_contract"),
    ):
        expected = object_sha256(plan["scientific_contract"][payload_key])
        if require_hash(summary.get(summary_key), summary_key) != expected or require_hash(plan.get(plan_key), plan_key) != expected:
            raise ContractError(f"{summary_key} does not reproduce")
    scientific = dict(summary["scientific_contract"])
    if scientific.get("policy_case") != "none" or integer(
        scientific.get("post_2023_periods"), "scientific post_2023_periods"
    ) != 0:
        raise ContractError("Calibration report is not the five-date no-policy case")
    if integer(scientific.get("target_count"), "scientific target_count") != EXPECTED_TARGETS:
        raise ContractError("Scientific contract does not contain twelve targets")
    if integer(
        scientific.get("transition_free_parameter_count"), "transition free parameters"
    ) != EXPECTED_PARAMETERS:
        raise ContractError("Scientific contract does not contain ten free parameters")

    target_rows = parse_target_rows(
        list(summary.get("target_fit_table") or []), "selected target table"
    )
    loss = math.fsum(row["loss_contribution"] for row in target_rows)
    assert_close(loss, summary.get("selected_transition_loss"), "selected loss", 1.0e-8)
    csv_target_rows = parse_target_rows(
        read_csv(report_dir / "selected_target_fit.csv"), "selected_target_fit.csv"
    )
    assert_same_object(
        [{key: row[key] for key in ("moment", "target", "model", "gap", "weight", "loss_contribution")} for row in target_rows],
        [{key: row[key] for key in ("moment", "target", "model", "gap", "weight", "loss_contribution")} for row in csv_target_rows],
        "embedded and CSV target tables",
    )
    parameter_rows = parse_parameter_rows(
        list(summary.get("estimated_parameters") or []), plan
    )
    csv_parameter_rows = parse_parameter_rows(
        read_csv(report_dir / "selected_parameters_and_bounds.csv"), plan
    )
    assert_same_object(parameter_rows, csv_parameter_rows, "embedded and CSV parameter tables")

    selection_path = report_dir / "selection.json"
    selection = read_json(selection_path)
    if require_hash(summary.get("selection_sha256"), "selection hash") != file_sha256(selection_path):
        raise ContractError("Selection receipt hash differs from selected report")
    selection_sha_path = validate_sha_sidecar(selection_path)
    summary_sha_path = validate_sha_sidecar(summary_path)
    for field in ("selected_candidate_id", "selected_candidate_sha256"):
        if selection.get(field) != summary.get(field):
            raise ContractError(f"Selection receipt disagrees on {field}")
    assert_close(
        selection.get("transition_loss"), loss, "selection transition loss", 1.0e-8
    )
    best = dict(summary.get("best_candidate") or {})
    theta = dict(best.get("theta") or {})
    if not theta:
        raise ContractError("Selected report omits the selected theta")
    assert_same_object(selection.get("theta"), theta, "selection theta")

    jacobian = validate_jacobian(plan_path, plan)
    required_selected_files = (
        "selected_transition_path.csv",
        "selected_dated_model_moments.csv",
        "selected_dated_period_fertility.csv",
        "selected_dated_period_fertility_by_age.csv",
        "selected_dated_cohort_timing_ledgers.json",
    )
    for name in required_selected_files:
        if not (report_dir / name).is_file():
            raise ContractError(f"Selected report omits {name}")
    return {
        "summary_path": summary_path,
        "report_dir": report_dir,
        "summary": summary,
        "plan_path": plan_path,
        "plan": plan,
        "selection_path": selection_path,
        "selection_sha_path": selection_sha_path,
        "summary_sha_path": summary_sha_path,
        "plan_sha_path": plan_sha_path,
        "selection": selection,
        "scientific_contract": scientific,
        "target_rows": target_rows,
        "parameter_rows": parameter_rows,
        "loss": loss,
        "best": best,
        "theta": theta,
        "jacobian": jacobian,
        "required_selected_files": required_selected_files,
    }


def validate_repeat(
    repeat_path: Path,
    selected: Mapping[str, Any],
    expected_replicate_id: int,
) -> dict[str, Any]:
    summary_path, task_dir = resolve_repeat_task(repeat_path)
    summary = read_json(summary_path)
    wrapper_path = task_dir / "refinement_task_contract.json"
    wrapper = read_json(wrapper_path)
    if wrapper.get("schema") != "e5f_transition_ridge_task_contract_v1" or wrapper.get("mode") != "repeat":
        raise ContractError(f"Repeat {expected_replicate_id} lacks the repeat-task contract")
    if integer(wrapper.get("replicate_id"), "repeat replicate_id") != expected_replicate_id:
        raise ContractError("Repeat inputs must be supplied in replicate-id order 1, 2")
    for field, expected in (
        ("candidate_id", selected["summary"]["selected_candidate_id"]),
        ("candidate_sha256", selected["summary"]["selected_candidate_sha256"]),
        ("plan_sha256", selected["summary"]["plan_sha256"]),
        ("scientific_contract_sha256", selected["summary"]["scientific_contract_sha256"]),
        ("renewal_contract_sha256", selected["summary"]["renewal_contract_sha256"]),
        ("dated_contract_sha256", selected["summary"]["dated_contract_sha256"]),
        ("selection_sha256", selected["summary"]["selection_sha256"]),
    ):
        if wrapper.get(field) != expected:
            raise ContractError(f"Repeat {expected_replicate_id} disagrees on {field}")
    execution_identity = str(wrapper.get("execution_identity") or "")
    expected_identities = list(selected["summary"]["repeat_gate"].get("execution_identities") or [])
    if len(expected_identities) != 2 or execution_identity != expected_identities[expected_replicate_id - 1]:
        raise ContractError("Repeat execution identity differs from the final collector")
    if summary.get("status") != "complete_transition_calibration_panel_task":
        raise ContractError(f"Repeat {expected_replicate_id} scientific solve is incomplete")
    assert_same_object(
        task_scientific_contract(summary),
        selected["scientific_contract"],
        f"repeat {expected_replicate_id} scientific contract",
    )
    if summary.get("policy_case") != "none" or integer(
        summary.get("post_2023_periods"), "repeat post_2023_periods"
    ) != 0:
        raise ContractError("A repeat includes policy or a post-2023 path")
    candidate = dict(summary.get("best_candidate") or {})
    if str(candidate.get("candidate")) == "":
        raise ContractError("Repeat candidate identifier is empty")
    assert_same_object(candidate.get("theta"), selected["theta"], "repeat theta")
    assert_close(
        summary.get("old_psi_child"),
        selected["best"].get("old_psi_child"),
        "repeat old_psi_child",
        PARAMETER_TOL,
    )
    assert_close(
        candidate.get("new_psi_child"),
        selected["best"].get("new_psi_child"),
        "repeat new_psi_child",
        PARAMETER_TOL,
    )
    rows = parse_target_rows(read_csv(task_dir / "target_fit_long.csv"), f"repeat {expected_replicate_id} target table")
    for reference, observed in zip(selected["target_rows"], rows):
        if observed["moment"] != reference["moment"]:
            raise ContractError("Repeat target order differs from selected report")
        for field in ("target", "weight"):
            assert_close(observed[field], reference[field], f"repeat {expected_replicate_id} {reference['moment']} {field}", 0.0)
        for field in ("model", "gap", "standardized_gap"):
            assert_close(observed[field], reference[field], f"repeat {expected_replicate_id} {reference['moment']} {field}", MOMENT_TOL)
        assert_close(observed["loss_contribution"], reference["loss_contribution"], f"repeat {expected_replicate_id} {reference['moment']} contribution", LOSS_TOL)
    loss = math.fsum(row["loss_contribution"] for row in rows)
    assert_close(loss, selected["loss"], f"repeat {expected_replicate_id} loss", LOSS_TOL)
    diagnostics = {
        "max_market_residual": candidate.get("max_market_residual"),
        "max_mass_residual": candidate.get("max_mass_residual"),
        "population_target_gap": candidate.get("population_target_gap"),
        "terminal_synthetic_childless_consistency_gap": candidate.get(
            "terminal_synthetic_childless_consistency_gap"
        ),
        "stationary_measurement_max_abs_gap": summary.get(
            "stationary_measurement_max_abs_gap"
        ),
    }
    for name, value in diagnostics.items():
        finite_float(value, f"repeat {expected_replicate_id} {name}")
    return {
        "summary_path": summary_path,
        "task_dir": task_dir,
        "summary": summary,
        "wrapper_path": wrapper_path,
        "wrapper": wrapper,
        "execution_identity": execution_identity,
        "candidate": candidate,
        "target_rows": rows,
        "loss": loss,
        "diagnostics": diagnostics,
    }


def local_artifact(packet_dir: Path, recorded: Mapping[str, Any], label: str) -> Path:
    expected_sha = require_hash(recorded.get("sha256"), f"{label}.sha256")
    expected_size = recorded.get("size_bytes")
    recorded_path = Path(str(recorded.get("path") or ""))
    candidates = [recorded_path, packet_dir / recorded_path.name]
    valid = []
    for candidate in candidates:
        candidate = candidate.resolve()
        if (
            candidate.is_file()
            and (expected_size is None or candidate.stat().st_size == integer(expected_size, f"{label}.size_bytes"))
            and file_sha256(candidate) == expected_sha
            and candidate not in valid
        ):
            valid.append(candidate)
    if len(valid) != 1:
        raise ContractError(
            f"Cannot resolve exactly one hash-matching local artifact for {label}: {candidates}"
        )
    return valid[0]


def validate_historical_packet(
    packet_path: Path, selected: Mapping[str, Any]
) -> dict[str, Any]:
    packet_dir = packet_path.resolve()
    provenance_path = packet_dir / "classification_and_provenance.json"
    provenance = read_json(provenance_path)
    if provenance.get("status") != "complete_e5f_2007_2023_historical_validation_packet":
        raise ContractError("Historical-validation packet is not complete")
    if tuple(provenance.get("calendar_years") or []) != EXPECTED_YEARS:
        raise ContractError("Historical packet does not contain the canonical five dates")
    gates = dict(provenance.get("validation_gates") or {})
    if not gates or any(value is not True for value in gates.values()):
        failed = sorted(name for name, value in gates.items() if value is not True)
        raise ContractError(f"Historical packet has failed validation gates: {failed}")
    input_records = dict(provenance.get("inputs") or {})
    if not input_records:
        raise ContractError("Historical packet has no input receipts")
    input_paths = {
        name: local_artifact(packet_dir, dict(record), f"historical input {name}")
        for name, record in input_records.items()
    }
    bundle_payload = "\n".join(
        f"{name}:{input_records[name]['sha256']}" for name in sorted(input_records)
    )
    if hashlib.sha256(bundle_payload.encode()).hexdigest() != require_hash(
        provenance.get("input_bundle_sha256"), "historical input bundle"
    ):
        raise ContractError("Historical input-bundle hash does not reproduce")
    model = dict(provenance.get("model_contract") or {})
    if model.get("no_policy_historical_path") is not True:
        raise ContractError("Historical packet is not explicitly no policy")
    scientific = selected["scientific_contract"]
    for historical_field, scientific_field in (
        ("source_sha256", "source_sha256"),
        ("target_fingerprint", "target_fingerprint"),
        ("code_bundle_sha256", None),
        ("target_set", "target_set"),
    ):
        expected = (
            scientific["code_fingerprints"]["bundle_sha256"]
            if scientific_field is None
            else scientific[scientific_field]
        )
        if model.get(historical_field) != expected:
            raise ContractError(f"Historical packet disagrees on {historical_field}")
    assert_same_object(
        model.get("renewal_accounting_contract"),
        scientific["renewal_accounting_contract"],
        "historical renewal contract",
    )
    reconciliation = dict(model.get("target_and_loss_reconciliation") or {})
    if integer(reconciliation.get("target_count"), "historical target_count") != EXPECTED_TARGETS:
        raise ContractError("Historical packet does not reconcile twelve terminal targets")
    if reconciliation.get("target_fingerprint") != scientific["target_fingerprint"]:
        raise ContractError("Historical target fingerprint differs")
    assert_close(
        reconciliation.get("recomputed_loss"),
        selected["loss"],
        "historical selected loss",
        LOSS_TOL,
    )

    outputs = dict(provenance.get("outputs") or {})
    comparison_path = local_artifact(
        packet_dir, dict(outputs.get("comparison_csv") or {}), "historical comparison"
    )
    fit_path = local_artifact(
        packet_dir, dict(outputs.get("fit_statistics_csv") or {}), "historical fit statistics"
    )
    comparison = read_csv(comparison_path)
    require_columns(
        comparison,
        ("calendar_year", "series_id", "series_label", "classification", "data_value", "model_value"),
        "historical comparison",
    )
    classifications = {row["classification"] for row in comparison}
    unknown = sorted(classifications.difference(HISTORICAL_CLASSIFICATIONS))
    if unknown:
        raise ContractError(f"Historical comparison has unknown classifications {unknown}")
    if not {"imposed bridge", "calibrated 2023 moment", "untargeted holdout", "numerical residual"}.issubset(
        {HISTORICAL_CLASSIFICATIONS[name] for name in classifications}
    ):
        raise ContractError("Historical packet omits one of the four required classifications")
    terminal = [row for row in comparison if row["classification"] == "terminal_maintained_target"]
    if len(terminal) != EXPECTED_TARGETS or {integer(row["calendar_year"], "terminal year") for row in terminal} != {2023}:
        raise ContractError("Historical packet must contain twelve terminal 2023 targets")
    terminal_by_name = {row["series_id"].removeprefix("terminal_target__"): row for row in terminal}
    if set(terminal_by_name) != {row["moment"] for row in selected["target_rows"]}:
        raise ContractError("Historical terminal target names differ from calibration")
    for target in selected["target_rows"]:
        row = terminal_by_name[target["moment"]]
        assert_close(row["data_value"], target["target"], f"historical {target['moment']} target", 0.0)
        assert_close(row["model_value"], target["model"], f"historical {target['moment']} model", MOMENT_TOL)

    fit_rows = read_csv(fit_path)
    if len(fit_rows) != 14:
        raise ContractError("Historical fit-statistics table must contain fourteen rows")
    require_columns(
        fit_rows,
        (
            "series_id",
            "classification",
            "evaluation_years",
            "scalar_statistics_status",
            "rmse_evaluation_window",
            "notes",
        ),
        "historical fit statistics",
    )
    if len({row["series_id"] for row in fit_rows}) != len(fit_rows):
        raise ContractError("Historical fit-statistics table has duplicate series")
    fit_classes = {
        "imposed_bridge_matched_by_construction_not_fit",
        "descriptive_noncomparable_holdout",
        "descriptive_noncomparable_holdout_sensitivity",
        "untargeted_holdout",
        "untargeted_holdout_data_only",
        "descriptive_noncomparable_diagnostic",
        "numerical_residual_audit_not_fit",
    }
    unknown_fit = sorted({row["classification"] for row in fit_rows}.difference(fit_classes))
    if unknown_fit:
        raise ContractError(f"Historical fit table has unknown classifications {unknown_fit}")
    for row in fit_rows:
        classification = row["classification"]
        status = row["scalar_statistics_status"]
        if classification in {
            "descriptive_noncomparable_holdout",
            "descriptive_noncomparable_holdout_sensitivity",
            "descriptive_noncomparable_diagnostic",
        } and status != "withheld_noncomparable_measurement":
            raise ContractError(
                f"Historical scalar claims were not withheld for {row['series_id']}"
            )
        if classification == "untargeted_holdout_data_only" and status != "not_available_data_only":
            raise ContractError(f"Data-only row {row['series_id']} has a model fit statistic")

    figures_record = dict(outputs.get("figures") or {})
    expected_figures = {
        "imposed_inputs_matched_by_construction",
        "untargeted_fertility_housing_validation",
        "untargeted_price_rent_validation",
        "numerical_residual_audit",
    }
    if set(figures_record) != expected_figures:
        raise ContractError("Historical packet does not contain the stable four-figure set")
    figure_paths: dict[str, list[Path]] = {}
    for name in sorted(figures_record):
        records = list(figures_record[name] or [])
        paths = [local_artifact(packet_dir, dict(record), f"historical figure {name}") for record in records]
        if {path.suffix for path in paths} != {".png", ".pdf"}:
            raise ContractError(f"Historical figure {name} must have PNG and PDF versions")
        figure_paths[name] = paths
    return {
        "packet_dir": packet_dir,
        "provenance_path": provenance_path,
        "provenance": provenance,
        "input_paths": input_paths,
        "comparison_path": comparison_path,
        "comparison": comparison,
        "fit_path": fit_path,
        "fit_rows": fit_rows,
        "figure_paths": figure_paths,
    }


def validate_continuation_packet(
    packet_path: Path,
    selected: Mapping[str, Any],
    historical: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate the exact packet written by the no-policy continuation driver."""
    packet_dir = packet_path.resolve()
    summary_path = packet_dir / "summary.json"
    manifest_path = packet_dir / "manifest.json"
    summary = read_json(summary_path)
    manifest = read_json(manifest_path)
    if summary.get("status") != "complete_no_policy_post2023_continuation_pair":
        raise ContractError("Paired continuation packet is not complete")
    if summary.get("policy_case") != "none" or summary.get("fiscal_change") != "none":
        raise ContractError("Continuation packet activates policy or a fiscal change")
    if manifest.get("status") != "complete_no_policy_post2023_continuation_manifest":
        raise ContractError("Continuation manifest is not complete")

    artifacts = dict(manifest.get("artifacts") or {})
    required_artifacts = {
        "README.md",
        "summary.json",
        "paired_continuation_path.csv",
        "shared_2007_2023_history.csv",
        "history_reproduction_audit.csv",
        "outside_share_invariance_audit.csv",
        "closed_stationary_schedule.csv",
        "closed_stationary_schedule_progress.csv",
        "closed_stationary_endpoint.json",
        "continuation_progress.csv",
        "latest_completed_period.json",
        "latest_endpoint_search.json",
        "open_stationary_endpoint.json",
        "open_endpoint/stationary_endpoint_search.csv",
        "open_endpoint/stationary_open_endpoint.csv",
        "paired_continuation_levels.png",
        "paired_continuation_levels.pdf",
        "paired_renewal_diagnostics.png",
        "paired_renewal_diagnostics.pdf",
        "closed_stationary_renewal_schedule.png",
        "closed_stationary_renewal_schedule.pdf",
    }
    if set(artifacts) != required_artifacts:
        raise ContractError(
            "Continuation manifest's required artifact set differs: "
            f"missing={sorted(required_artifacts.difference(artifacts))}, "
            f"extra={sorted(set(artifacts).difference(required_artifacts))}"
        )
    artifact_paths: dict[str, Path] = {}
    for name, expected_sha in artifacts.items():
        path = (packet_dir / name).resolve()
        if not path.is_file():
            raise ContractError(f"Continuation artifact is missing: {name}")
        if file_sha256(path) != require_hash(expected_sha, f"continuation {name}"):
            raise ContractError(f"Continuation artifact hash differs: {name}")
        artifact_paths[name] = path
    driver_path = Path(__file__).resolve().with_name(
        "run_e5f_post2023_no_policy_continuations.py"
    )
    if not driver_path.is_file():
        raise ContractError("Local no-policy continuation driver is missing")
    if file_sha256(driver_path) != require_hash(
        manifest.get("driver_sha256"), "continuation driver hash"
    ):
        raise ContractError("Continuation packet was built by a different driver")

    provenance = dict(summary.get("provenance") or {})
    assert_same_object(
        manifest.get("input_provenance"), provenance, "continuation input provenance"
    )
    scientific = selected["scientific_contract"]
    expected_provenance = {
        "selected_report_sha256": file_sha256(selected["summary_path"]),
        "selected_case_transition_sha256": file_sha256(
            selected["report_dir"] / "selected_transition_path.csv"
        ),
        "source_sha256": scientific["source_sha256"],
        "renewal_contract_sha256": selected["summary"]["renewal_contract_sha256"],
        "target_fingerprint": scientific["target_fingerprint"],
        "code_bundle_sha256": scientific["code_fingerprints"]["bundle_sha256"],
        "scientific_contract_sha256": selected["summary"][
            "scientific_contract_sha256"
        ],
        "selection_sha256": selected["summary"]["selection_sha256"],
    }
    for field, expected in expected_provenance.items():
        if provenance.get(field) != expected:
            raise ContractError(f"Continuation packet disagrees on {field}")
    if "selected_task_summary_sha256" in provenance:
        historical_task_sha = historical["provenance"].get("inputs", {}).get(
            "case_summary", {}
        ).get("sha256")
        if provenance.get("selected_task_summary_sha256") != historical_task_sha:
            raise ContractError(
                "Continuation and historical packets use different selected task summaries"
            )
    for field in provenance:
        require_hash(provenance[field], f"continuation provenance {field}")
    assert_same_object(
        summary.get("renewal_accounting_contract"),
        scientific["renewal_accounting_contract"],
        "continuation renewal contract",
    )
    if summary.get("target_set") != scientific["target_set"] or summary.get(
        "target_fingerprint"
    ) != scientific["target_fingerprint"]:
        raise ContractError("Continuation target contract differs from calibration")
    assert_close(
        summary.get("outside_origin_entry_share"),
        scientific["outside_origin_entry_share"],
        "continuation outside-origin share",
        0.0,
    )
    if "not identified" not in str(summary.get("outside_origin_entry_status", "")):
        raise ContractError("Continuation does not disclose that the outside share is external")

    history_audit = dict(summary.get("history_reproduction_audit") or {})
    if history_audit.get("status") != "passed":
        raise ContractError("Continuation failed the 2007--2023 reproduction audit")
    maximum_history_gap = finite_float(
        history_audit.get("maximum_absolute_gap"), "maximum history gap"
    )
    tolerance = finite_float(history_audit.get("tolerance"), "history tolerance")
    if tolerance > 5.0e-10 or maximum_history_gap > tolerance:
        raise ContractError("Continuation history reproduction exceeds 5e-10")
    paired = dict(summary.get("paired_initial_state_audit") or {})
    if paired.get("status") != "passed" or finite_float(
        paired.get("shared_2007_2023_history_maximum_absolute_gap"),
        "shared 2007--2023 history gap",
    ) != 0.0:
        raise ContractError("Closed and open paths do not share the exact 2023 state")
    if integer(paired.get("common_history_periods"), "common history periods") != 5:
        raise ContractError("Continuation does not preserve the five dated history states")
    post_2023_periods = integer(paired.get("post_2023_periods"), "post-2023 periods")
    if post_2023_periods != 40:
        raise ContractError("Canonical continuation must contain exactly 40 post-2023 dates")
    require_hash(
        paired.get("matched_2023_pre_fertility_distribution_sha256"),
        "matched 2023 distribution hash",
    )
    require_hash(
        paired.get("matched_2023_birth_vintage_queue_sha256"),
        "matched 2023 birth-vintage queue hash",
    )

    paths_summary = dict(summary.get("paths") or {})
    expected_paths = {"closed_national_benchmark", "open_cbsa_sensitivity"}
    if set(paths_summary) != expected_paths:
        raise ContractError("Continuation packet must contain exactly closed and open paths")
    for name, record in paths_summary.items():
        gates = dict(record.get("path_gates") or {})
        if finite_float(gates.get("maximum_market_residual"), f"{name} market residual") > 2.0e-4:
            raise ContractError(f"{name} fails the market-residual gate")
        if abs(finite_float(gates.get("maximum_mass_residual"), f"{name} mass residual")) > 2.0e-8:
            raise ContractError(f"{name} fails the mass-accounting gate")
        if finite_float(gates.get("maximum_feasibility_projection_mass"), f"{name} projection mass") > 1.0e-6:
            raise ContractError(f"{name} fails the feasibility-projection gate")
        if finite_float(gates.get("minimum_distribution_mass"), f"{name} minimum mass") < -1.0e-14:
            raise ContractError(f"{name} has negative distribution mass")
        if integer(gates.get("maximum_nonfinite_distribution_count"), f"{name} nonfinite count") != 0:
            raise ContractError(f"{name} has nonfinite distribution entries")
        if record.get("2023", {}).get("policy_case") != "none" or record.get(
            "last_date", {}
        ).get("policy_case") != "none":
            raise ContractError(f"{name} activates policy")
    closed = dict(summary.get("closed_stationary_endpoint") or {})
    open_endpoint = dict(summary.get("open_stationary_endpoint") or {})
    if closed.get("policy_case") != "none":
        raise ContractError("Closed stationary endpoint activates policy")
    label_allowed = closed.get("between_steady_states_label_allowed") is True
    statement = str(summary.get("between_steady_states_statement") or "")
    if label_allowed:
        if closed.get("usable_closed_root") is not True or "verified positive" not in statement:
            raise ContractError("Between-steady-states label lacks a verified closed root")
    elif "not described as a transition between steady states" not in statement:
        raise ContractError("No-root closed case is not described conservatively")
    if open_endpoint.get("usable_open_endpoint") is not True:
        raise ContractError("Open stationary endpoint is not usable")
    open_outside_flow = finite_float(
        open_endpoint.get("outside_flow_M"), "open endpoint outside flow"
    )
    open_population = finite_float(
        open_endpoint.get("stationary_population_scale"),
        "open endpoint population scale",
    )
    open_entry_per_adult = finite_float(
        open_endpoint.get("entry_households_per_adult"),
        "open endpoint entry per adult",
    )
    open_realized_outside_share = open_outside_flow / (
        open_population * open_entry_per_adult
    )
    if not 0.0 < open_realized_outside_share < 1.0:
        raise ContractError("Open endpoint realized outside-origin share is invalid")
    invariance = dict(summary.get("outside_share_invariance_audit") or {})
    if invariance.get("status") != "passed" or finite_float(
        invariance.get("maximum_absolute_identity_residual"), "outside-share identity"
    ) > 2.0e-10:
        raise ContractError("Outside-share invariance audit failed")

    paired_path = artifact_paths["paired_continuation_path.csv"]
    paired_rows = read_csv(paired_path)
    require_columns(
        paired_rows,
        (
            "scenario",
            "phase",
            "policy_case",
            "calendar_year",
            "population_index_2023",
            "asset_price_index",
            "topcode_adjusted_births_per_adult",
            "owner_rate",
        ),
        "paired continuation path",
    )
    if {row["scenario"] for row in paired_rows} != expected_paths or any(
        row["policy_case"] != "none" for row in paired_rows
    ):
        raise ContractError("Paired continuation CSV has the wrong scenarios or policy scope")
    expected_years = [2007, 2011, 2015, 2019, 2023] + [
        2027 + 4 * index for index in range(post_2023_periods)
    ]
    expected_phases = (
        ["shared_calibrated_history"] * 4
        + ["shared_matched_2023_state"]
        + ["post_2023_continuation"] * post_2023_periods
    )
    for scenario in sorted(expected_paths):
        scenario_rows = [row for row in paired_rows if row["scenario"] == scenario]
        years = [integer(row["calendar_year"], f"{scenario} calendar year") for row in scenario_rows]
        phases = [row["phase"] for row in scenario_rows]
        if years != expected_years or phases != expected_phases:
            raise ContractError(
                f"{scenario} must contain the exact 2007--2183 dated path"
            )

    comparison_rows: list[dict[str, Any]] = []
    metrics = (
        ("adult-household population index", "population_index_2023"),
        ("housing-cost index", "asset_price_index"),
        ("births per adult household", "topcode_adjusted_births_per_adult"),
        ("ownership rate", "owner_rate"),
    )
    for closure in sorted(paths_summary):
        record = paths_summary[closure]
        start = record["2023"]
        terminal = record["last_date"]
        endpoint = closed if closure == "closed_national_benchmark" else open_endpoint
        endpoint_usable = (
            endpoint.get("usable_closed_root") is True
            if closure == "closed_national_benchmark"
            else endpoint.get("usable_open_endpoint") is True
        )
        for label, key in metrics:
            stationary_value: Any = ""
            if endpoint_usable:
                if key == "population_index_2023" and "stationary_population_scale" in endpoint:
                    stationary_value = finite_float(
                        endpoint["stationary_population_scale"], "stationary scale"
                    ) / finite_float(start["adult_population"], "2023 adult population")
                elif key == "asset_price_index" and "price_ratio" in endpoint:
                    stationary_value = finite_float(endpoint["price_ratio"], "price ratio")
            comparison_rows.append(
                {
                    "closure": closure,
                    "metric": label,
                    "start_2023": finite_float(start[key], f"{closure} 2023 {key}"),
                    "finite_horizon_terminal": finite_float(
                        terminal[key], f"{closure} terminal {key}"
                    ),
                    "terminal_year": integer(terminal["calendar_year"], "terminal year"),
                    "stationary_endpoint": stationary_value,
                    "endpoint_status": endpoint.get("status"),
                    "between_steady_states_label_allowed": (
                        label_allowed if closure == "closed_national_benchmark" else False
                    ),
                    "interpretation": (
                        "M=0 and rho=1; national closed benchmark"
                        if closure == "closed_national_benchmark"
                        else (
                            "old-state normalization of "
                            f"{finite_float(summary['outside_origin_entry_share'], 'outside share'):.6g} "
                            "with fixed M and rho thereafter; endpoint realized share "
                            f"{open_realized_outside_share:.6g}; open/CBSA sensitivity"
                        )
                    ),
                }
            )

    figure_names = (
        "paired_continuation_levels",
        "paired_renewal_diagnostics",
        "closed_stationary_renewal_schedule",
    )
    figure_paths = {
        name: [
            artifact_paths[f"{name}.png"],
            artifact_paths[f"{name}.pdf"],
        ]
        for name in figure_names
    }
    return {
        "packet_dir": packet_dir,
        "summary_path": summary_path,
        "manifest_path": manifest_path,
        "summary": summary,
        "manifest": manifest,
        "artifact_paths": artifact_paths,
        "paired_path": paired_path,
        "comparison_rows": comparison_rows,
        "open_realized_outside_share": open_realized_outside_share,
        "figure_paths": figure_paths,
        "driver_path": driver_path,
    }


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def make_target_figure(rows: Sequence[Mapping[str, Any]], output_dir: Path) -> None:
    labels = [str(row["label"]) for row in rows]
    values = [float(row["standardized_gap"]) for row in rows]
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    positions = np.arange(len(rows))
    colors = ["#b4493d" if value < 0 else "#2f6f9f" for value in values]
    ax.barh(positions, values, color=colors)
    ax.axvline(0.0, color="#333333", linewidth=0.8)
    ax.set_yticks(positions, labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Standardized gap: sqrt(weight) x (model - target)")
    ax.set_title("2023 calibration fit: all twelve maintained moments")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / f"calibrated_2023_target_fit.{suffix}", dpi=220 if suffix == "png" else None)
    plt.close(fig)


def make_jacobian_figure(jacobian: Mapping[str, Any], output_dir: Path) -> None:
    singular = np.asarray(jacobian["singular_values"], dtype=float)
    relative = singular / singular[0]
    fig, ax = plt.subplots(figsize=(6.5, 3.9))
    ax.semilogy(np.arange(1, len(relative) + 1), relative, "o-", color="#2f6f9f")
    ax.axhline(1.0e-3, color="#b4493d", linestyle="--", label="declared rank threshold")
    ax.set_xticks(np.arange(1, len(relative) + 1))
    ax.set_xlabel("Singular-value order")
    ax.set_ylabel("Relative to largest")
    ax.set_title("Local weighted-Jacobian spectrum")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / f"local_identification_spectrum.{suffix}", dpi=220 if suffix == "png" else None)
    plt.close(fig)


def numerical_residual_rows(
    selected: Mapping[str, Any],
    repeats: Sequence[Mapping[str, Any]],
    historical: Mapping[str, Any],
    continuation: Mapping[str, Any],
) -> list[dict[str, Any]]:
    tolerances = dict(selected["summary"].get("eligibility_tolerances") or {})
    required = {
        "max_market_residual",
        "population_target_gap",
        "max_mass_residual",
        "stationary_measurement_max_abs_gap",
        "terminal_synthetic_childless_consistency_gap",
    }
    if set(tolerances) != required:
        raise ContractError("Selected report omits the exact numerical eligibility contract")
    rows: list[dict[str, Any]] = []

    def append(scope: str, check: str, value: Any, tolerance: Any) -> None:
        absolute_value = abs(finite_float(value, f"{scope} {check}"))
        bound = finite_float(tolerance, f"{scope} {check} tolerance")
        if bound < 0.0 or absolute_value > bound:
            raise ContractError(
                f"Numerical gate failed: {scope} {check}={absolute_value:g} > {bound:g}"
            )
        rows.append(
            {
                "scope": scope,
                "check": check,
                "absolute_value": absolute_value,
                "tolerance": bound,
                "status": "passed_numerical_gate_not_empirical_moment",
            }
        )

    for check in (
        "max_market_residual",
        "population_target_gap",
        "max_mass_residual",
        "terminal_synthetic_childless_consistency_gap",
    ):
        append("selected validation", check, selected["best"].get(check), tolerances[check])
    for index, repeat in enumerate(repeats, start=1):
        for check, value in repeat["diagnostics"].items():
            append(f"repeat {index}", check, value, tolerances[check])

    historical_tolerances = {
        "market_clearing_residual": 2.0e-4,
        "mass_accounting_residual": 2.0e-8,
        "population_bridge_gap": 2.0e-10,
        "age_bridge_max_gap": 2.0e-10,
    }
    for series, tolerance in historical_tolerances.items():
        values = [
            finite_float(row["model_value"], f"historical {series}")
            for row in historical["comparison"]
            if row["series_id"] == series
        ]
        if values:
            append("historical 2007--2023", series, max(abs(value) for value in values), tolerance)

    continuation_tolerances = {
        "maximum_market_residual": 2.0e-4,
        "maximum_mass_residual": 2.0e-8,
        "maximum_feasibility_projection_mass": 1.0e-6,
        "maximum_nonfinite_distribution_count": 0.0,
    }
    for closure, record in continuation["summary"]["paths"].items():
        gates = record["path_gates"]
        for check, tolerance in continuation_tolerances.items():
            append(f"post-2023 {closure}", check, gates[check], tolerance)
    append(
        "post-2023 shared history",
        "maximum_history_reproduction_gap",
        continuation["summary"]["history_reproduction_audit"]["maximum_absolute_gap"],
        5.0e-10,
    )
    append(
        "post-2023 paired origin",
        "shared_2023_state_maximum_absolute_gap",
        continuation["summary"]["paired_initial_state_audit"][
            "shared_2007_2023_history_maximum_absolute_gap"
        ],
        0.0,
    )
    return rows


def markdown_number(value: Any, digits: int = 6) -> str:
    number = finite_float(value, "markdown value")
    if number == 0:
        return "0"
    if abs(number) >= 10000 or abs(number) < 1.0e-4:
        return f"{number:.3e}"
    return f"{number:.{digits}g}"


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> list[str]:
    output = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join("---" for _ in headers) + "|",
    ]
    output.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return output


def build_readme(
    selected: Mapping[str, Any],
    target_provenance: Mapping[str, Any],
    repeats: Sequence[Mapping[str, Any]],
    historical: Mapping[str, Any],
    continuation: Mapping[str, Any],
    numerical: Sequence[Mapping[str, Any]],
) -> str:
    target_rows = selected["target_rows"]
    parameter_rows = selected["parameter_rows"]
    classification_counts = Counter(
        HISTORICAL_CLASSIFICATIONS[row["classification"]]
        for row in historical["comparison"]
    )
    outside_share = finite_float(
        selected["scientific_contract"]["outside_origin_entry_share"],
        "outside-origin share",
    )
    provenance_counts = ", ".join(
        f"{name}={count}"
        for name, count in target_provenance["status_counts"].items()
    )
    lines = [
        "# E5F transition: canonical no-policy quantitative report",
        "",
        (
            f"The selected transition calibration has loss **{selected['loss']:.6f}**. "
            "It is a 2023 calibration performed along the 2007--2023 transition, "
            "not a stationary calibration relabeled as 2023. Two fresh executions "
            "reproduce the selected parameters, every target moment, and the loss "
            "within the declared tolerances."
        ),
        "",
        "This directory contains no funded policy exercise. The historical population "
        "path is an imposed Census/ACS bridge; the data series not in the 2023 objective "
        "remain holdouts; and the post-2023 closed and open paths are alternative closure "
        "assumptions, not two fitted forecasts.",
        "",
        "## What is imposed, calibrated, and checked",
        "",
    ]
    lines += markdown_table(
        ("Object", "Role", "Where to inspect"),
        (
            ("Census/ACS household totals and head-age shares, 2007--2023", "Imposed bridge; matched by construction", "historical_model_data.csv"),
            ("Twelve maintained moments measured on the 2023 model state", "Calibration objective", "calibrated_2023_target_fit.csv"),
            ("Historical fertility, timing, tenure, rooms, prices, and rents", "Untargeted holdouts or noncomparable context", "historical_model_data.csv"),
            ("Market, mass, and bridge gaps", "Numerical residuals, not empirical moments", "historical_model_data.csv"),
            ("Closed national continuation", "No outside inflow; alternative post-2023 closure", "continuation_comparison.csv"),
            (
                "Open continuation",
                f"Old-state outside-origin normalization {outside_share:.6g}; fixed M and rho thereafter",
                "continuation_comparison.csv",
            ),
        ),
    )
    lines += [
        "",
        "Historical row counts by presentation class: "
        + ", ".join(f"{name}={classification_counts[name]}" for name in sorted(classification_counts))
        + ".",
        "",
        "## 2023 target fit",
        "",
    ]
    lines += markdown_table(
        ("Moment", "Target", "Model", "Gap", "Weight", "Contribution"),
        [
            (
                row["label"],
                markdown_number(row["target"]),
                markdown_number(row["model"]),
                markdown_number(row["gap"]),
                markdown_number(row["weight"]),
                markdown_number(row["loss_contribution"]),
            )
            for row in target_rows
        ],
    )
    lines += [
        "",
        f"The contributions sum to {selected['loss']:.9f}. The table is complete: no "
        "target row is suppressed or reclassified after estimation.",
        "",
        "## Target provenance and measurement",
        "",
        f"All twelve calibration rows have an explicit source/measurement receipt. Status counts: {provenance_counts}.",
        "",
    ]
    lines += markdown_table(
        ("Moment", "Status", "Authoritative builder", "Sample", "Estimator / measurement", "Uncertainty", "Caveat"),
        [
            (
                row["label"],
                row["status"],
                row["authoritative_builder"],
                row["sample_geography_vintage"],
                row["estimator_measurement"],
                f"{row['uncertainty_type']}; scale={markdown_number(row['standard_error_or_scale'])}",
                row["caveat"],
            )
            for row in target_provenance["rows"]
        ],
    )
    lines += [
        "",
        "## Estimated parameters and bounds",
        "",
    ]
    lines += markdown_table(
        ("Parameter", "Estimate", "Bounds", "Transform", "Bound status", "Status"),
        [
            (
                row["label"],
                markdown_number(row["estimate"]),
                f"[{markdown_number(row['lower'])}, {markdown_number(row['upper'])}]",
                row["transform"],
                "within 2% of unit bound" if row["near_bound_within_0.02"] else "interior",
                row["status"],
            )
            for row in parameter_rows
        ],
    )
    jacobian = selected["jacobian"]
    lines += [
        "",
        "## Local identification and numerical reproducibility",
        "",
        (
            f"The 12-by-10 weighted Jacobian has rank {jacobian['rank_at_relative_1e3']} "
            f"at the predeclared relative singular-value threshold 1e-3 and condition "
            f"number {jacobian['condition_number']:.4g}. "
            f"Side-consistency froze {len(jacobian['frozen_dimensions'])} coordinate(s) "
            "for the ridge proposal. This is local numerical rank for the maintained "
            "moment system; it is not a claim that the causal housing-cost-to-fertility "
            "slope is empirically identified."
        ),
        "",
    ]
    lines += markdown_table(
        ("Execution", "Identity", "Loss", "Status"),
        [
            ("selected validation", selected["summary"]["selected_validation_execution_identity"], f"{selected['loss']:.12f}", "passed"),
            *[
                (f"repeat {index}", repeat["execution_identity"], f"{repeat['loss']:.12f}", "passed")
                for index, repeat in enumerate(repeats, start=1)
            ],
        ],
    )
    lines += [
        "",
        "## Numerical residuals",
        "",
        "These are solver, accounting, and reproduction gates. They are not empirical moments and never enter the calibration loss.",
        "",
    ]
    lines += markdown_table(
        ("Scope", "Numerical check", "Absolute value", "Tolerance"),
        [
            (
                row["scope"],
                row["check"],
                markdown_number(row["absolute_value"]),
                markdown_number(row["tolerance"]),
            )
            for row in numerical
        ],
    )
    lines += [
        "",
        "## 2007--2023 historical validation",
        "",
        "The five-date table keeps the source builder's classification on every row. "
        "Population and age composition are inputs, not evidence of fit. Historical "
        "fertility, first-birth timing, ownership, rooms, and price/rent series were "
        "not added to the objective after seeing the transition.",
        "",
        "The observed price-to-rent denominator is rent CPI, while the model denominator "
        "is its implicit user-cost flow. The packet therefore withholds a scalar "
        "price-to-rent fit claim.",
        "",
        "## Post-2023 no-policy continuations",
        "",
    ]
    continuation_rows = []
    for row in continuation["comparison_rows"]:
        continuation_rows.append(
            (
                row["closure"],
                row["metric"],
                markdown_number(row["start_2023"]),
                markdown_number(row["finite_horizon_terminal"]),
                str(integer(row["terminal_year"], "terminal year")),
                row["interpretation"],
            )
        )
    lines += markdown_table(
        ("Closure", "Metric", "2023", "Terminal", "Year", "Interpretation"),
        continuation_rows,
    )
    closed_endpoint = continuation["summary"]["closed_stationary_endpoint"]
    open_endpoint = continuation["summary"]["open_stationary_endpoint"]
    open_realized_outside_share = finite_float(
        continuation["open_realized_outside_share"],
        "open endpoint realized outside-origin share",
    )
    lines += [
        "",
        (
            f"Closed endpoint status: `{closed_endpoint['status']}`. "
            f"{continuation['summary']['between_steady_states_statement']}"
        ),
        "",
        f"Open endpoint status: `{open_endpoint['status']}`. This endpoint is an open/CBSA sensitivity, not a national demographic forecast.",
        f"Its realized outside-origin share is {open_realized_outside_share:.6f}; "
        f"the value {outside_share:.6g} applies only to the old-state normalization.",
        "",
        "The closed path sets outside entry to zero. In the open path, "
        f"{outside_share:.6g} normalizes the old-state decomposition; M and rho are then fixed, "
        "so the realized outside-origin share can change. Their difference is a closure sensitivity, not "
        "a policy effect. A finite terminal date is called a steady state only when the "
        "continuation packet's explicit endpoint/root and drift gates pass.",
        "",
        "## Files and provenance",
        "",
        "- `calibrated_2023_target_fit.csv`: all twelve target rows.",
        "- `target_provenance.csv`: authoritative builder, sample, estimator, "
        "  uncertainty, status, and caveat for every target row.",
        "- `estimated_parameters_and_bounds.csv`: all ten free parameters.",
        "- `weighted_jacobian.csv`, `singular_values.csv`, `side_consistency.csv`, and "
        "  `weak_directions.csv`: local numerical-identification diagnostics.",
        "- `numerical_residuals.csv`: solver, accounting, repeat, and path gates kept "
        "  separate from empirical moments.",
        "- `historical_model_data.csv` and `historical_fit_statistics.csv`: classified "
        "  five-date evidence.",
        "- `continuation_*.csv`: the two no-policy paths and their comparison.",
        "- `report_manifest.json`: exact input/output paths, sizes, and SHA-256 hashes.",
        "",
        "The report builder is presentation-only: it does not alter estimates, solve the "
        "model, or rescale any series.",
        "",
    ]
    return "\n".join(lines)


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    if len(args.repeat_task) != 2:
        raise ContractError("Exactly two --repeat-task inputs are required")
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ContractError(f"Refusing to overwrite existing output directory: {output_dir}")

    selected = validate_selected_report(args.selected_report, args.refinement_plan)
    target_provenance = validate_target_provenance(
        args.target_provenance_csv, selected
    )
    repeats = [
        validate_repeat(path, selected, index)
        for index, path in enumerate(args.repeat_task, start=1)
    ]
    identities = [selected["summary"]["selected_validation_execution_identity"]] + [
        repeat["execution_identity"] for repeat in repeats
    ]
    if len(set(identities)) != 3:
        raise ContractError("Selected validation and repeats lack three distinct identities")
    historical = validate_historical_packet(args.historical_packet, selected)
    continuation = validate_continuation_packet(
        args.continuation_packet, selected, historical
    )
    numerical = numerical_residual_rows(
        selected, repeats, historical, continuation
    )

    # No output is touched until all cross-packet gates above pass.
    output_dir.mkdir(parents=True)
    write_csv(output_dir / "calibrated_2023_target_fit.csv", selected["target_rows"])
    write_csv(
        output_dir / "target_provenance.csv", target_provenance["rows"]
    )
    write_csv(output_dir / "estimated_parameters_and_bounds.csv", selected["parameter_rows"])
    write_csv(output_dir / "numerical_residuals.csv", numerical)
    write_csv(output_dir / "repeat_checks.csv", [
        {
            "replicate_id": index,
            "execution_identity": repeat["execution_identity"],
            "loss": repeat["loss"],
            "loss_gap": repeat["loss"] - selected["loss"],
            "max_model_moment_gap": max(
                abs(row["model"] - reference["model"])
                for row, reference in zip(repeat["target_rows"], selected["target_rows"])
            ),
            "status": "passed_exact_repeat_gate",
        }
        for index, repeat in enumerate(repeats, start=1)
    ])
    copy_file(selected["plan_path"], output_dir / "refinement_plan.json")
    copy_file(selected["selection_path"], output_dir / "selection.json")
    for name, path in selected["jacobian"]["paths"].items():
        copy_file(path, output_dir / name)
    for name in selected["required_selected_files"]:
        copy_file(selected["report_dir"] / name, output_dir / name)
    copy_file(historical["comparison_path"], output_dir / "historical_model_data.csv")
    copy_file(historical["fit_path"], output_dir / "historical_fit_statistics.csv")
    copy_file(
        historical["provenance_path"], output_dir / "historical_classification_and_provenance.json"
    )
    classification_counts = Counter(
        HISTORICAL_CLASSIFICATIONS[row["classification"]]
        for row in historical["comparison"]
    )
    write_csv(output_dir / "historical_classification_summary.csv", [
        {"presentation_class": name, "row_count": count}
        for name, count in sorted(classification_counts.items())
    ])
    for figure_name, paths in historical["figure_paths"].items():
        for path in paths:
            copy_file(path, output_dir / f"historical_{figure_name}{path.suffix}")

    copy_file(continuation["summary_path"], output_dir / "continuation_summary.json")
    copy_file(continuation["manifest_path"], output_dir / "continuation_manifest.json")
    write_csv(output_dir / "continuation_comparison.csv", continuation["comparison_rows"])
    for name, path in continuation["artifact_paths"].items():
        if name in {
            "README.md",
            "summary.json",
            "paired_continuation_levels.png",
            "paired_continuation_levels.pdf",
            "paired_renewal_diagnostics.png",
            "paired_renewal_diagnostics.pdf",
            "closed_stationary_renewal_schedule.png",
            "closed_stationary_renewal_schedule.pdf",
        }:
            continue
        destination = "continuation_" + name.replace("/", "_")
        copy_file(path, output_dir / destination)
    for figure_name, paths in continuation["figure_paths"].items():
        for path in paths:
            copy_file(path, output_dir / f"continuation_{figure_name}{path.suffix}")

    make_target_figure(selected["target_rows"], output_dir)
    make_jacobian_figure(selected["jacobian"], output_dir)
    (output_dir / "README.md").write_text(
        build_readme(
            selected,
            target_provenance,
            repeats,
            historical,
            continuation,
            numerical,
        )
    )

    input_paths = {
        "selected_report_summary": selected["summary_path"],
        "selected_report_summary_sha256_sidecar": selected["summary_sha_path"],
        "selection": selected["selection_path"],
        "selection_sha256_sidecar": selected["selection_sha_path"],
        "refinement_plan": selected["plan_path"],
        "refinement_plan_sha256_sidecar": selected["plan_sha_path"],
        "target_provenance": target_provenance["path"],
        **{f"refinement_{name}": path for name, path in selected["jacobian"]["paths"].items()},
        **{f"selected_{name}": selected["report_dir"] / name for name in selected["required_selected_files"]},
        **{f"repeat_{index}_summary": repeat["summary_path"] for index, repeat in enumerate(repeats, start=1)},
        **{f"repeat_{index}_contract": repeat["wrapper_path"] for index, repeat in enumerate(repeats, start=1)},
        **{f"repeat_{index}_target_fit": repeat["task_dir"] / "target_fit_long.csv" for index, repeat in enumerate(repeats, start=1)},
        "historical_provenance": historical["provenance_path"],
        "historical_comparison": historical["comparison_path"],
        "historical_fit_statistics": historical["fit_path"],
        **{
            f"historical_source_{name}": path
            for name, path in historical["input_paths"].items()
        },
        "continuation_summary": continuation["summary_path"],
        "continuation_manifest": continuation["manifest_path"],
        "continuation_driver": continuation["driver_path"],
        **{
            f"continuation_artifact_{name.replace('/', '__')}": path
            for name, path in continuation["artifact_paths"].items()
        },
    }
    for group, figure_paths in historical["figure_paths"].items():
        for path in figure_paths:
            input_paths[f"historical_figure_{group}_{path.suffix.lstrip('.')}"] = path
    for group, figure_paths in continuation["figure_paths"].items():
        for path in figure_paths:
            input_paths[f"continuation_figure_{group}_{path.suffix.lstrip('.')}"] = path
    input_receipts = {name: receipt(path) for name, path in sorted(input_paths.items())}
    input_bundle = object_sha256(
        [[name, row["sha256"]] for name, row in sorted(input_receipts.items())]
    )
    output_receipts = {
        path.name: receipt(path)
        for path in sorted(output_dir.iterdir())
        if path.is_file() and path.name not in {"report_manifest.json", "report_manifest.sha256"}
    }
    manifest = {
        "schema": "e5f_canonical_no_policy_transition_report_v1",
        "status": "complete_canonical_no_policy_transition_report",
        "scope": "theory-facing quantitative evidence; no funded policy material",
        "selected_transition_loss": selected["loss"],
        "target_count": EXPECTED_TARGETS,
        "free_parameter_count": EXPECTED_PARAMETERS,
        "target_set": selected["scientific_contract"]["target_set"],
        "target_fingerprint": selected["scientific_contract"]["target_fingerprint"],
        "target_provenance_status_counts": target_provenance["status_counts"],
        "scientific_contract_sha256": selected["summary"]["scientific_contract_sha256"],
        "renewal_contract_sha256": selected["summary"]["renewal_contract_sha256"],
        "dated_contract_sha256": selected["summary"]["dated_contract_sha256"],
        "selection_sha256": selected["summary"]["selection_sha256"],
        "selected_candidate_sha256": selected["summary"]["selected_candidate_sha256"],
        "input_bundle_sha256": input_bundle,
        "input_receipts": input_receipts,
        "output_receipts": output_receipts,
        "validation_gates": {
            "selected_refinement_complete": True,
            "two_independent_exact_repeats": True,
            "twelve_target_loss_reconciled": True,
            "twelve_row_target_provenance_complete": True,
            "all_free_parameters_and_bounds_present": True,
            "weighted_jacobian_recomputed": True,
            "local_rank_gate_passed": True,
            "historical_packet_complete_and_classified": True,
            "imposed_bridge_not_relabelled_as_fit": True,
            "untargeted_holdouts_not_added_to_objective": True,
            "paired_closed_open_continuations_complete": True,
            "all_inputs_no_policy": True,
            "no_posthoc_rescaling_or_refit": True,
        },
    }
    manifest_path = output_dir / "report_manifest.json"
    write_json(manifest_path, manifest)
    manifest_sha = file_sha256(manifest_path)
    (output_dir / "report_manifest.sha256").write_text(
        f"{manifest_sha}  report_manifest.json\n"
    )
    return manifest


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    result = build_report(args)
    print(
        json.dumps(
            {
                "status": result["status"],
                "loss": result["selected_transition_loss"],
                "target_fingerprint": result["target_fingerprint"],
                "input_bundle_sha256": result["input_bundle_sha256"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    try:
        main()
    except ContractError as error:
        print(f"E5F_CANONICAL_NO_POLICY_REPORT_REJECTED {error}", file=sys.stderr)
        raise SystemExit(2) from error
