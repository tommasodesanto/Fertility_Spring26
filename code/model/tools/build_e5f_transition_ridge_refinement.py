#!/usr/bin/env python3
"""Build a fail-closed ridge/trust refinement from an E5F coordinate panel.

The input is an already evaluated central-coordinate panel.  This script does
no model solving: it reconstructs the weighted residual Jacobian, reports its
rank and side consistency, and writes seven fully hashed centers for a separate
scientific validation run (one observed incumbent and six predeclared
ridge/trust proposals).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np


EXPECTED_TASKS = 21
EXPECTED_DIMENSIONS = 10
CENTRAL_STEP = 0.02
EXPECTED_TARGET_COUNT = 12
RANK_THRESHOLDS = (1.0e-2, 1.0e-3, 1.0e-4)
IDENTIFICATION_THRESHOLD = 1.0e-3
RIDGE_RELATIVE_GRID = (1.0e-3, 1.0e-2, 1.0e-1)
TRUST_RADIUS_GRID = (0.01, 0.02)
ELIGIBILITY_TOLERANCES = {
    "max_market_residual": 2.0e-4,
    "population_target_gap": 2.0e-10,
    "max_mass_residual": 2.0e-10,
    "stationary_measurement_max_abs_gap": 2.0e-8,
    "terminal_synthetic_childless_consistency_gap": 2.0e-10,
}


class ContractError(RuntimeError):
    """Raised when an input or output violates the declared contract."""


@dataclass(frozen=True)
class Expectations:
    source_sha256: str
    target_set: str
    target_fingerprint: str
    code_bundle_sha256: str
    model_profile: str
    panel_seed: int
    dimensions: int = EXPECTED_DIMENSIONS
    source: str | None = None
    renewal_contract_sha256: str | None = None
    dated_contract_sha256: str | None = None


@dataclass
class PanelTask:
    task_id: int
    task_dir: Path
    summary: dict[str, Any]
    candidate: dict[str, Any]
    panel: dict[str, Any]
    target_rows: list[dict[str, Any]]
    residual: np.ndarray
    unit: np.ndarray
    loss: float


@dataclass
class CoordinatePanel:
    results_dir: Path
    tasks: list[PanelTask]
    contract: dict[str, Any]
    contract_sha256: str
    domain: list[dict[str, Any]]
    moment_names: list[str]
    targets: np.ndarray
    weights: np.ndarray
    input_artifact_sha256: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coordinate-results-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-source", default=None)
    parser.add_argument("--expected-target-set", required=True)
    parser.add_argument("--expected-target-fingerprint", required=True)
    parser.add_argument("--expected-code-bundle-sha256", required=True)
    parser.add_argument("--expected-model-profile", required=True)
    parser.add_argument("--expected-panel-seed", type=int, required=True)
    parser.add_argument(
        "--expected-dimensions",
        type=int,
        default=EXPECTED_DIMENSIONS,
        help="Joint transition-parameter count; coordinate tasks must equal 1+2d.",
    )
    parser.add_argument("--expected-renewal-contract-sha256", required=True)
    parser.add_argument("--expected-dated-contract-sha256", required=True)
    return parser.parse_args()


def canonical_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")


def object_sha256(payload: Any) -> str:
    return hashlib.sha256(canonical_bytes(payload)).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContractError(f"Cannot read JSON artifact {path}: {error}") from error
    if not isinstance(payload, dict):
        raise ContractError(f"JSON artifact must contain an object: {path}")
    return payload


def read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except OSError as error:
        raise ContractError(f"Cannot read CSV artifact {path}: {error}") from error


def atomic_write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def atomic_write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ContractError(f"Refusing to write an empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _finite_float(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ContractError(f"{label} is not numeric: {value!r}") from error
    if not math.isfinite(result):
        raise ContractError(f"{label} is nonfinite: {value!r}")
    return result


def _exact_float(value: Any, expected: float, label: str, tol: float = 1.0e-14) -> None:
    actual = _finite_float(value, label)
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tol):
        raise ContractError(f"{label}={actual!r}; expected {expected!r}")


def population_bridge_contract(recorded: Any) -> dict[str, Any]:
    """Retain only exogenous population-bridge provenance.

    Candidate-dependent stationary masses are diagnostics and cannot define a
    cross-candidate contract.
    """
    if not isinstance(recorded, dict):
        raise ContractError("Summary omits the population-bridge contract")
    initial = recorded.get("initial_age_reweight")
    if not isinstance(initial, dict):
        raise ContractError("Population bridge omits initial_age_reweight")
    bridge_fields = (
        "status",
        "target_indices",
        "total_households",
        "age_distribution",
    )
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
    missing = [name for name in bridge_fields if name not in recorded]
    missing += [f"initial_age_reweight.{name}" for name in initial_fields if name not in initial]
    if missing:
        raise ContractError(f"Population bridge omits required fields: {missing}")
    return {
        **{name: recorded[name] for name in bridge_fields},
        "initial_age_reweight": {name: initial[name] for name in initial_fields},
    }


def validate_renewal_accounting(summary: dict[str, Any]) -> None:
    contract = dict(summary.get("renewal_accounting_contract") or {})
    old = dict(summary.get("renewal_accounting_old_state") or {})
    replacement = _finite_float(contract.get("replacement_fertility"), "replacement fertility")
    conversion = _finite_float(
        contract.get("effective_birth_to_household_conversion"),
        "birth-to-household conversion",
    )
    _exact_float(replacement, 2.1, "replacement fertility")
    _exact_float(conversion, 1.0 / 2.1, "birth-to-household conversion", 1.0e-15)
    if int(contract.get("birth_vintage_queue_waiting_slots", -1)) != 4:
        raise ContractError("Renewal contract must use four queue waiting slots")
    if int(contract.get("birth_to_entry_effect_lag_dates", -1)) != 5:
        raise ContractError("Renewal contract must report a five-date effect lag")
    _exact_float(
        contract.get("birth_to_entry_effect_lag_years"),
        20.0,
        "birth-to-entry lag years",
    )
    names = (
        "old_entry_flow_E",
        "old_queue_mature_flow_B",
        "old_queue_B_over_E",
        "old_renewal_residual",
        "outside_flow_M",
        "outside_origin_entry_share",
        "retention_rho",
    )
    values = {name: _finite_float(old.get(name), f"renewal old state {name}") for name in names}
    entry = values["old_entry_flow_E"]
    births = values["old_queue_mature_flow_B"]
    outside = values["outside_flow_M"]
    share = values["outside_origin_entry_share"]
    retention = values["retention_rho"]
    if min(entry, births, outside) <= 0.0 or not 0.0 < share < 1.0:
        raise ContractError("Old renewal accounting has invalid flows or outside share")
    if not 0.0 <= retention <= 1.0:
        raise ContractError("Old renewal retention must lie in [0,1]")
    completed = _finite_float(
        summary.get("old_model_completed_fertility"), "old completed fertility"
    )
    identities = (
        (values["old_queue_B_over_E"], births / entry, "reported B/E"),
        (values["old_queue_B_over_E"], completed / replacement, "fertility/replacement"),
        (outside, share * entry, "outside flow"),
        (retention, (1.0 - share) * entry / births, "retention"),
        (
            values["old_renewal_residual"],
            entry - outside - retention * births,
            "renewal residual",
        ),
    )
    for reported, reconstructed, label in identities:
        if not math.isclose(reported, reconstructed, rel_tol=0.0, abs_tol=2.0e-10):
            raise ContractError(
                f"Old renewal {label} failed: {reported} versus {reconstructed}"
            )
    if abs(values["old_queue_B_over_E"] - 1.0) > 5.0e-4:
        raise ContractError("Old state is not within 5e-4 of replacement")


def scientific_contract(summary: dict[str, Any]) -> dict[str, Any]:
    """Return the exact cross-task contract used by builder and collector."""
    required = (
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
        "population_bridge",
        "population_validation_status",
        "outside_origin_entry_share",
        "old_completed_fertility_reference",
        "policy_case",
        "post_2023_periods",
        "stationary_free_parameter_count",
        "transition_free_parameter_count",
    )
    missing = [name for name in required if name not in summary]
    if missing:
        raise ContractError(f"Summary omits exact contract fields: {missing}")
    return {
        "source": str(summary["source"]),
        "source_sha256": str(summary["source_sha256"]),
        "source_metadata": summary["source_metadata"],
        "target_set": str(summary["target_set"]),
        "target_fingerprint": str(summary["target_fingerprint"]),
        "target_count": int(summary["target_count"]),
        "code_fingerprints": summary["code_fingerprints"],
        "model_profile": summary["model_profile"],
        "income_profile_gates": summary["income_profile_gates"],
        "target_measurements": summary["target_measurements"],
        "renewal_accounting_contract": summary["renewal_accounting_contract"],
        "dated_measurement_contract": summary["dated_measurement_contract"],
        "population_bridge": population_bridge_contract(summary["population_bridge"]),
        "population_validation_status": str(summary["population_validation_status"]),
        "outside_origin_entry_share": _finite_float(
            summary["outside_origin_entry_share"], "outside-origin entry share"
        ),
        "old_completed_fertility_reference": _finite_float(
            summary["old_completed_fertility_reference"],
            "old completed fertility reference",
        ),
        "policy_case": str(summary["policy_case"]),
        "post_2023_periods": int(summary["post_2023_periods"]),
        "stationary_free_parameter_count": int(summary["stationary_free_parameter_count"]),
        "transition_free_parameter_count": int(summary["transition_free_parameter_count"]),
    }


def check_expectations(contract: dict[str, Any], expectations: Expectations) -> None:
    checks = {
        "source_sha256": expectations.source_sha256,
        "target_set": expectations.target_set,
        "target_fingerprint": expectations.target_fingerprint,
        "code_fingerprints.bundle_sha256": expectations.code_bundle_sha256,
        "model_profile.name": expectations.model_profile,
    }
    actual = {
        "source_sha256": contract["source_sha256"],
        "target_set": contract["target_set"],
        "target_fingerprint": contract["target_fingerprint"],
        "code_fingerprints.bundle_sha256": str(
            dict(contract["code_fingerprints"]).get("bundle_sha256")
        ),
        "model_profile.name": str(dict(contract["model_profile"]).get("name")),
    }
    for label, expected in checks.items():
        if actual[label] != str(expected):
            raise ContractError(f"{label}={actual[label]!r}; expected {expected!r}")
    if expectations.source is not None and contract["source"] != expectations.source:
        raise ContractError(
            f"source={contract['source']!r}; expected {expectations.source!r}"
        )
    if contract["target_count"] != EXPECTED_TARGET_COUNT:
        raise ContractError(
            f"target_count={contract['target_count']}; expected {EXPECTED_TARGET_COUNT}"
        )
    if contract["transition_free_parameter_count"] != expectations.dimensions:
        raise ContractError("Transition parameter count differs from the declared dimension")
    if contract["stationary_free_parameter_count"] != expectations.dimensions:
        raise ContractError("Stationary parameter count differs from the declared dimension")
    if contract["policy_case"] != "none" or contract["post_2023_periods"] != 0:
        raise ContractError("Coordinate panel must be a no-policy 2007--2023 calibration")
    _exact_float(
        contract["old_completed_fertility_reference"],
        2.1,
        "old completed-fertility reference",
    )
    if expectations.renewal_contract_sha256 is None:
        raise ContractError("An expected renewal-contract hash is mandatory")
    actual_hash = object_sha256(contract["renewal_accounting_contract"])
    if actual_hash != expectations.renewal_contract_sha256:
        raise ContractError(
            f"renewal contract hash={actual_hash}; expected "
            f"{expectations.renewal_contract_sha256}"
        )
    if expectations.dated_contract_sha256 is None:
        raise ContractError("An expected dated-contract hash is mandatory")
    actual_hash = object_sha256(contract["dated_measurement_contract"])
    if actual_hash != expectations.dated_contract_sha256:
        raise ContractError(
            f"dated contract hash={actual_hash}; expected "
            f"{expectations.dated_contract_sha256}"
        )


def validate_candidate_eligibility(summary: dict[str, Any]) -> None:
    candidate = dict(summary.get("best_candidate") or {})
    loss = _finite_float(candidate.get("transition_loss"), "transition loss")
    if loss < 0.0:
        raise ContractError("Transition loss cannot be negative")
    for name, tolerance in ELIGIBILITY_TOLERANCES.items():
        if name == "stationary_measurement_max_abs_gap":
            value = abs(_finite_float(summary.get(name), name))
        else:
            if name not in candidate:
                raise ContractError(f"Best candidate omits gate {name}")
            value = abs(_finite_float(candidate[name], name))
        if value > tolerance:
            raise ContractError(f"Candidate fails {name}: {value} > {tolerance}")


def parse_target_rows(
    raw_rows: list[dict[str, str]],
    target_count: int,
    task_id: int,
    candidate_label: str | None = None,
) -> tuple[list[dict[str, Any]], np.ndarray, float]:
    required = {
        "candidate",
        "moment",
        "target",
        "model",
        "gap",
        "weight",
        "loss_contribution",
        "standardized_gap",
    }
    if len(raw_rows) != target_count:
        raise ContractError(
            f"task {task_id}: target table has {len(raw_rows)} rows, expected {target_count}"
        )
    parsed: list[dict[str, Any]] = []
    residual: list[float] = []
    seen: set[str] = set()
    for index, raw in enumerate(raw_rows):
        missing = sorted(required.difference(raw))
        if missing:
            raise ContractError(f"task {task_id}: target row omits {missing}")
        moment = str(raw["moment"])
        if candidate_label is not None and str(raw["candidate"]) != candidate_label:
            raise ContractError(
                f"task {task_id}: target row candidate={raw['candidate']!r}; "
                f"expected {candidate_label!r}"
            )
        if not moment or moment in seen:
            raise ContractError(f"task {task_id}: duplicate or empty target {moment!r}")
        seen.add(moment)
        target = _finite_float(raw["target"], f"task {task_id} {moment} target")
        model = _finite_float(raw["model"], f"task {task_id} {moment} model")
        gap = _finite_float(raw["gap"], f"task {task_id} {moment} gap")
        weight = _finite_float(raw["weight"], f"task {task_id} {moment} weight")
        contribution = _finite_float(
            raw["loss_contribution"], f"task {task_id} {moment} loss contribution"
        )
        standardized = _finite_float(
            raw["standardized_gap"], f"task {task_id} {moment} standardized gap"
        )
        if weight <= 0.0:
            raise ContractError(f"task {task_id} {moment}: weight must be positive")
        reconstructed_gap = model - target
        reconstructed_standardized = math.sqrt(weight) * reconstructed_gap
        reconstructed_contribution = reconstructed_standardized**2
        if not math.isclose(gap, reconstructed_gap, rel_tol=1.0e-11, abs_tol=1.0e-12):
            raise ContractError(f"task {task_id} {moment}: gap is inconsistent")
        if not math.isclose(
            standardized,
            reconstructed_standardized,
            rel_tol=1.0e-10,
            abs_tol=1.0e-11,
        ):
            raise ContractError(f"task {task_id} {moment}: standardized gap is inconsistent")
        if not math.isclose(
            contribution,
            reconstructed_contribution,
            rel_tol=1.0e-9,
            abs_tol=1.0e-10,
        ):
            raise ContractError(f"task {task_id} {moment}: loss contribution is inconsistent")
        parsed.append(
            {
                "index": index,
                "moment": moment,
                "target": target,
                "model": model,
                "gap": gap,
                "weight": weight,
                "loss_contribution": contribution,
                "standardized_gap": standardized,
            }
        )
        residual.append(reconstructed_standardized)
    return parsed, np.asarray(residual, dtype=float), float(sum(row["loss_contribution"] for row in parsed))


def _assert_close_vector(
    actual: np.ndarray, expected: np.ndarray, label: str, tol: float = 2.0e-12
) -> None:
    if actual.shape != expected.shape or not np.allclose(
        actual, expected, rtol=0.0, atol=tol
    ):
        max_gap = math.inf if actual.shape != expected.shape else float(np.max(np.abs(actual - expected)))
        raise ContractError(f"{label} mismatch; maximum absolute gap={max_gap}")


def inverse_unit(value: float, lower: float, upper: float, kind: str) -> float:
    if not math.isfinite(value) or value < lower - 1.0e-12 or value > upper + 1.0e-12:
        raise ContractError(
            f"Parameter value {value!r} lies outside declared bounds [{lower},{upper}]"
        )
    value = min(upper, max(lower, value))
    if kind == "log":
        return math.log(value / lower) / math.log(upper / lower)
    if kind == "discount":
        return 1.0 - math.sqrt(max(0.0, 1.0 - (value - lower) / (upper - lower)))
    if kind == "asinh":
        return (math.asinh(value) - math.asinh(lower)) / (
            math.asinh(upper) - math.asinh(lower)
        )
    if kind == "softzero":
        return math.sqrt(max(0.0, (value - lower) / (upper - lower)))
    raise ContractError(f"Unsupported unit transform: {kind}")


def reconstruct_candidate_unit(
    candidate: dict[str, Any],
    panel: dict[str, Any],
    domain: list[dict[str, Any]],
) -> np.ndarray:
    theta = dict(candidate.get("theta") or {})
    if not theta:
        raise ContractError("Best candidate omits theta")
    if "old_psi_child" not in panel or "new_psi_child" not in panel:
        raise ContractError("Panel metadata omits realized old/new child preferences")
    terminal_change = _finite_float(panel["new_psi_child"], "new psi_child") - _finite_float(
        panel["old_psi_child"], "old psi_child"
    )
    values: list[float] = []
    for row in domain:
        name = str(row["name"])
        if name == "beta_annual":
            quarterly = _finite_float(theta.get("beta"), "theta.beta")
            if quarterly <= 0.0:
                raise ContractError("theta.beta must be positive")
            parameter_value = quarterly**0.25
        elif name == "psi_child_change_2023":
            parameter_value = terminal_change
        else:
            parameter_value = _finite_float(theta.get(name), f"theta.{name}")
        values.append(
            inverse_unit(
                parameter_value,
                float(row["lower"]),
                float(row["upper"]),
                str(row["transform"]),
            )
        )
    return np.asarray(values, dtype=float)


def validate_preference_metadata(
    summary: dict[str, Any], candidate: dict[str, Any], panel: dict[str, Any]
) -> None:
    summary_old = _finite_float(summary.get("old_psi_child"), "summary old_psi_child")
    panel_old = _finite_float(panel.get("old_psi_child"), "panel old_psi_child")
    panel_new = _finite_float(panel.get("new_psi_child"), "panel new_psi_child")
    candidate_new = _finite_float(
        candidate.get("new_psi_child"), "candidate new_psi_child"
    )
    panel_change = _finite_float(
        panel.get("psi_child_change_2023"), "panel psi_child_change_2023"
    )
    identities = (
        (summary_old, panel_old, "summary old psi versus panel old psi"),
        (candidate_new, panel_new, "candidate new psi versus panel new psi"),
        (panel_change, panel_new - panel_old, "panel preference change versus new-old"),
    )
    for left, right, label in identities:
        if not math.isclose(left, right, rel_tol=0.0, abs_tol=5.0e-12):
            raise ContractError(f"{label} failed: {left} versus {right}")


def load_coordinate_panel(results_dir: Path, expectations: Expectations) -> CoordinatePanel:
    results_dir = results_dir.resolve()
    if not results_dir.is_dir():
        raise ContractError(f"Coordinate results directory does not exist: {results_dir}")
    dimensions = int(expectations.dimensions)
    if dimensions < 1 or dimensions >= EXPECTED_TARGET_COUNT:
        raise ContractError(
            "The declared parameter dimension must be positive and smaller than "
            "the target count"
        )
    expected_tasks = 1 + 2 * dimensions
    expected_names = {f"task_{task:03d}" for task in range(1, expected_tasks + 1)}
    observed_names = {
        path.name
        for path in results_dir.glob("task_*")
        if path.is_dir()
    }
    missing = sorted(expected_names.difference(observed_names))
    unexpected = sorted(observed_names.difference(expected_names))
    if missing or unexpected:
        raise ContractError(
            f"Coordinate panel must contain exactly task_001..task_{expected_tasks:03d}; "
            f"missing={missing}, unexpected={unexpected}"
        )

    tasks: list[PanelTask] = []
    contracts: list[dict[str, Any]] = []
    artifact_hashes: list[dict[str, Any]] = []
    moment_names: list[str] | None = None
    target_values: np.ndarray | None = None
    weights: np.ndarray | None = None
    common_domain: list[dict[str, Any]] | None = None
    common_center_hash: str | None = None

    for task_id in range(1, expected_tasks + 1):
        task_dir = results_dir / f"task_{task_id:03d}"
        failure = task_dir / "failure.json"
        if failure.exists():
            raise ContractError(f"task {task_id} has a failure artifact: {failure}")
        summary_path = task_dir / "summary.json"
        fit_path = task_dir / "target_fit_long.csv"
        if not summary_path.is_file() or not fit_path.is_file():
            raise ContractError(f"task {task_id} lacks summary.json or target_fit_long.csv")
        summary = read_json(summary_path)
        if summary.get("status") != "complete_transition_calibration_panel_task":
            raise ContractError(f"task {task_id} has non-complete status")
        candidate = dict(summary.get("best_candidate") or {})
        panel = dict(summary.get("panel_design") or {})
        if int(panel.get("task_id", -1)) != task_id:
            raise ContractError(f"task {task_id}: internal task id mismatch")
        if int(panel.get("panel_size", -1)) != expected_tasks:
            raise ContractError(f"task {task_id}: panel size is not {expected_tasks}")
        if int(panel.get("panel_seed", -1)) != expectations.panel_seed:
            raise ContractError(f"task {task_id}: panel seed mismatch")
        if panel.get("panel_design") != "coordinate":
            raise ContractError(f"task {task_id}: panel design is not coordinate")
        _exact_float(panel.get("local_radius"), CENTRAL_STEP, f"task {task_id} local radius")
        if panel.get("terminal_preference_coordinate") != "psi_child_change_2023":
            raise ContractError(f"task {task_id}: wrong terminal preference coordinate")
        domain = list(panel.get("domain") or [])
        if len(domain) != dimensions:
            raise ContractError(f"task {task_id}: domain does not have {dimensions} rows")
        domain_names = [str(row.get("name")) for row in domain]
        if len(set(domain_names)) != dimensions:
            raise ContractError(f"task {task_id}: duplicate domain names")
        for row in domain:
            if set(row) != {"name", "lower", "upper", "transform"}:
                raise ContractError(f"task {task_id}: malformed domain row {row}")
            lower = _finite_float(row["lower"], "domain lower")
            upper = _finite_float(row["upper"], "domain upper")
            if not lower < upper:
                raise ContractError(f"task {task_id}: non-increasing parameter bounds")
        if common_domain is None:
            common_domain = domain
        elif canonical_bytes(domain) != canonical_bytes(common_domain):
            raise ContractError(f"task {task_id}: domain differs across panel")
        center_hash = str(panel.get("center_sha256") or "")
        if len(center_hash) != 64:
            raise ContractError(f"task {task_id}: center_sha256 is missing")
        if common_center_hash is None:
            common_center_hash = center_hash
        elif center_hash != common_center_hash:
            raise ContractError(f"task {task_id}: center hash differs across panel")
        unit = np.asarray(panel.get("unit_vector"), dtype=float)
        if unit.shape != (dimensions,) or not np.all(np.isfinite(unit)):
            raise ContractError(f"task {task_id}: invalid unit vector")
        if np.any(unit < 0.0) or np.any(unit > 1.0):
            raise ContractError(f"task {task_id}: unit vector lies outside [0,1]")
        validate_preference_metadata(summary, candidate, panel)
        reconstructed_unit = reconstruct_candidate_unit(candidate, panel, domain)
        _assert_close_vector(
            unit,
            reconstructed_unit,
            f"task {task_id} parameter-to-unit reconstruction",
            5.0e-12,
        )

        validate_renewal_accounting(summary)
        validate_candidate_eligibility(summary)
        contract = scientific_contract(summary)
        check_expectations(contract, expectations)
        contracts.append(contract)

        parsed_rows, residual, table_loss = parse_target_rows(
            read_csv(fit_path),
            int(summary["target_count"]),
            task_id,
            str(candidate.get("candidate")),
        )
        candidate_loss = _finite_float(candidate.get("transition_loss"), "candidate loss")
        if not math.isclose(candidate_loss, table_loss, rel_tol=1.0e-9, abs_tol=1.0e-8):
            raise ContractError(
                f"task {task_id}: candidate loss {candidate_loss} != target-table loss {table_loss}"
            )
        names = [row["moment"] for row in parsed_rows]
        target = np.asarray([row["target"] for row in parsed_rows], dtype=float)
        weight = np.asarray([row["weight"] for row in parsed_rows], dtype=float)
        if moment_names is None:
            moment_names = names
            target_values = target
            weights = weight
        else:
            if names != moment_names:
                raise ContractError(f"task {task_id}: target order or names differ")
            _assert_close_vector(target, target_values, f"task {task_id} target values", 0.0)
            _assert_close_vector(weight, weights, f"task {task_id} target weights", 0.0)
        tasks.append(
            PanelTask(
                task_id=task_id,
                task_dir=task_dir,
                summary=summary,
                candidate=candidate,
                panel=panel,
                target_rows=parsed_rows,
                residual=residual,
                unit=unit,
                loss=candidate_loss,
            )
        )
        artifact_hashes.append(
            {
                "task_id": task_id,
                "summary_sha256": file_sha256(summary_path),
                "target_fit_sha256": file_sha256(fit_path),
            }
        )

    reference_contract = contracts[0]
    reference_bytes = canonical_bytes(reference_contract)
    for task_id, contract in enumerate(contracts[1:], start=2):
        if canonical_bytes(contract) != reference_bytes:
            raise ContractError(f"task {task_id}: mixed scientific contract")

    anchor = tasks[0].unit
    if np.any(anchor <= CENTRAL_STEP) or np.any(anchor >= 1.0 - CENTRAL_STEP):
        raise ContractError(
            "Anchor is within h=.02 of a unit bound; central differences would be clipped"
        )
    if tasks[0].panel.get("design") != "anchor":
        raise ContractError("task 1 is not the anchor")
    _assert_close_vector(tasks[0].unit, anchor, "anchor")
    for dimension in range(dimensions):
        minus_task = tasks[1 + 2 * dimension]
        plus_task = tasks[2 + 2 * dimension]
        expected_minus = anchor.copy()
        expected_plus = anchor.copy()
        expected_minus[dimension] -= CENTRAL_STEP
        expected_plus[dimension] += CENTRAL_STEP
        expected_minus_design = f"coordinate_{dimension}_minus"
        expected_plus_design = f"coordinate_{dimension}_plus"
        if minus_task.panel.get("design") != expected_minus_design:
            raise ContractError(
                f"task {minus_task.task_id}: expected {expected_minus_design}"
            )
        if plus_task.panel.get("design") != expected_plus_design:
            raise ContractError(
                f"task {plus_task.task_id}: expected {expected_plus_design}"
            )
        _assert_close_vector(
            minus_task.unit, expected_minus, f"coordinate {dimension} minus unit"
        )
        _assert_close_vector(
            plus_task.unit, expected_plus, f"coordinate {dimension} plus unit"
        )

    assert common_domain is not None
    assert moment_names is not None and target_values is not None and weights is not None
    return CoordinatePanel(
        results_dir=results_dir,
        tasks=tasks,
        contract=reference_contract,
        contract_sha256=object_sha256(reference_contract),
        domain=common_domain,
        moment_names=moment_names,
        targets=target_values,
        weights=weights,
        input_artifact_sha256=object_sha256(artifact_hashes),
    )


def matrix_rank_report(matrix: np.ndarray) -> dict[str, Any]:
    if matrix.ndim != 2 or min(matrix.shape) < 1:
        raise ContractError("Rank report requires a nonempty matrix")
    singular = np.linalg.svd(matrix, compute_uv=False)
    if singular.size == 0 or not np.all(np.isfinite(singular)):
        raise ContractError("Jacobian has invalid singular values")
    leading = float(singular[0])
    ranks = {
        f"relative_{threshold:g}": int(np.sum(singular >= threshold * leading))
        if leading > 0.0
        else 0
        for threshold in RANK_THRESHOLDS
    }
    smallest = float(singular[-1])
    return {
        "shape": [int(value) for value in matrix.shape],
        "singular_values": singular.tolist(),
        "relative_ranks": ranks,
        "condition_number": leading / smallest if smallest > 0.0 else None,
    }


def central_jacobian(panel: CoordinatePanel) -> tuple[np.ndarray, list[dict[str, Any]], np.ndarray]:
    anchor_residual = panel.tasks[0].residual
    columns: list[np.ndarray] = []
    side_rows: list[dict[str, Any]] = []
    frozen = np.zeros(len(panel.domain), dtype=bool)
    for dimension, domain_row in enumerate(panel.domain):
        minus = panel.tasks[1 + 2 * dimension].residual
        plus = panel.tasks[2 + 2 * dimension].residual
        forward = (plus - anchor_residual) / CENTRAL_STEP
        backward = (anchor_residual - minus) / CENTRAL_STEP
        central = (plus - minus) / (2.0 * CENTRAL_STEP)
        forward_norm = float(np.linalg.norm(forward))
        backward_norm = float(np.linalg.norm(backward))
        dot = float(forward @ backward)
        if forward_norm == 0.0 or backward_norm == 0.0:
            consistency = 0.0
        else:
            consistency = dot / (forward_norm * backward_norm)
        freeze = bool(dot <= 0.0 or forward_norm == 0.0 or backward_norm == 0.0)
        frozen[dimension] = freeze
        columns.append(central)
        side_rows.append(
            {
                "dimension": dimension,
                "parameter": str(domain_row["name"]),
                "minus_task": 2 + 2 * dimension,
                "plus_task": 3 + 2 * dimension,
                "forward_norm": forward_norm,
                "backward_norm": backward_norm,
                "forward_backward_dot": dot,
                "side_consistency_cosine": consistency,
                "relative_side_difference": float(
                    np.linalg.norm(forward - backward)
                    / max(forward_norm + backward_norm, np.finfo(float).tiny)
                ),
                "frozen_for_step": freeze,
                "freeze_rule": "forward_backward_dot_nonpositive",
            }
        )
    return np.column_stack(columns), side_rows, frozen


def trust_and_box_scale(
    raw_step: np.ndarray,
    anchor: np.ndarray,
    trust_radius: float,
) -> tuple[np.ndarray, float, str]:
    raw_norm = float(np.linalg.norm(raw_step))
    factors = [(1.0, "unscaled")]
    if raw_norm > trust_radius:
        factors.append((trust_radius / raw_norm, "trust_radius"))
    for index, value in enumerate(raw_step):
        if value > 0.0:
            factors.append(((1.0 - anchor[index]) / value, f"upper_bound_{index}"))
        elif value < 0.0:
            factors.append((anchor[index] / -value, f"lower_bound_{index}"))
    scale, reason = min(factors, key=lambda item: item[0])
    scale = max(0.0, min(1.0, float(scale)))
    step = raw_step * scale
    candidate = anchor + step
    if np.any(candidate < -1.0e-13) or np.any(candidate > 1.0 + 1.0e-13):
        raise ContractError("Trust/box scaling failed to keep a candidate inside [0,1]")
    return step, scale, reason


def ridge_trust_proposals(
    anchor_residual: np.ndarray,
    jacobian: np.ndarray,
    anchor: np.ndarray,
    frozen: np.ndarray,
) -> list[dict[str, Any]]:
    active_indices = np.flatnonzero(~frozen)
    if active_indices.size < 1:
        raise ContractError("All coordinates fail the predeclared side-consistency gate")
    active = jacobian[:, active_indices]
    active_rank = matrix_rank_report(active)
    active_key = f"relative_{IDENTIFICATION_THRESHOLD:g}"
    if int(active_rank["relative_ranks"][active_key]) != int(active_indices.size):
        raise ContractError(
            "Active Jacobian is rank deficient at the 1e-3 relative threshold; "
            "ridge regularization cannot substitute for identification"
        )
    leading = float(active_rank["singular_values"][0])
    gram = active.T @ active
    score = active.T @ anchor_residual
    proposals: list[dict[str, Any]] = []
    for ridge_relative in RIDGE_RELATIVE_GRID:
        ridge_absolute = ridge_relative * leading**2
        try:
            active_raw = np.linalg.solve(
                gram + ridge_absolute * np.eye(active_indices.size), -score
            )
        except np.linalg.LinAlgError as error:
            raise ContractError(f"Ridge solve failed: {error}") from error
        raw = np.zeros(anchor.size, dtype=float)
        raw[active_indices] = active_raw
        for trust_radius in TRUST_RADIUS_GRID:
            step, scale, scale_reason = trust_and_box_scale(raw, anchor, trust_radius)
            candidate = anchor + step
            predicted_residual = anchor_residual + jacobian @ step
            proposals.append(
                {
                    "kind": "ridge_trust",
                    "ridge_relative": ridge_relative,
                    "ridge_absolute": ridge_absolute,
                    "trust_radius": trust_radius,
                    "raw_step": raw.tolist(),
                    "raw_step_norm": float(np.linalg.norm(raw)),
                    "applied_scale": scale,
                    "applied_scale_reason": scale_reason,
                    "step_vector": step.tolist(),
                    "step_norm": float(np.linalg.norm(step)),
                    "unit_vector": candidate.tolist(),
                    "predicted_residual": predicted_residual.tolist(),
                    "predicted_loss": float(predicted_residual @ predicted_residual),
                    "active_dimensions": active_indices.tolist(),
                    "frozen_dimensions": np.flatnonzero(frozen).tolist(),
                }
            )
    if len(proposals) != 6:
        raise AssertionError("The predeclared ridge/trust grid must generate six proposals")
    return proposals


def transform_unit(u: float, lower: float, upper: float, kind: str) -> float:
    if not 0.0 <= u <= 1.0:
        raise ContractError(f"Unit coordinate outside [0,1]: {u}")
    if kind == "log":
        return lower * (upper / lower) ** u
    if kind == "discount":
        return lower + (upper - lower) * (1.0 - (1.0 - u) ** 2)
    if kind == "asinh":
        return math.sinh((1.0 - u) * math.asinh(lower) + u * math.asinh(upper))
    if kind == "softzero":
        return lower + (upper - lower) * u**2
    raise ContractError(f"Unsupported unit transform: {kind}")


def center_payload(
    panel: CoordinatePanel,
    unit: np.ndarray,
    candidate_id: int,
    construction: dict[str, Any],
) -> dict[str, Any]:
    anchor_theta = dict(panel.tasks[0].candidate.get("theta") or {})
    if not anchor_theta:
        raise ContractError("Anchor best candidate omits theta")
    theta = dict(anchor_theta)
    terminal_delta: float | None = None
    if len(unit) != len(panel.domain):
        raise ContractError("Candidate unit vector and domain have different lengths")
    for coordinate, row in zip(unit, panel.domain):
        name = str(row["name"])
        value = transform_unit(
            float(coordinate),
            float(row["lower"]),
            float(row["upper"]),
            str(row["transform"]),
        )
        if name == "beta_annual":
            theta["beta"] = value**4
        elif name == "psi_child_change_2023":
            terminal_delta = value
        else:
            theta[name] = value
    if terminal_delta is None:
        raise ContractError("Domain omits psi_child_change_2023")
    return {
        "schema": "e5f_transition_ridge_candidate_v1",
        "candidate_id": candidate_id,
        "scientific_contract_sha256": panel.contract_sha256,
        "coordinate_panel_input_sha256": panel.input_artifact_sha256,
        "unit_vector": unit.tolist(),
        "construction": construction,
        "best_candidate": {
            "theta": theta,
            # Only the difference is consumed by the repaired scientific driver.
            "old_psi_child": 0.0,
            "new_psi_child": terminal_delta,
        },
    }


def weak_direction_rows(jacobian: np.ndarray, domain: list[dict[str, Any]]) -> list[dict[str, Any]]:
    _, singular, right = np.linalg.svd(jacobian, full_matrices=False)
    rows: list[dict[str, Any]] = []
    for order, singular_index in enumerate(range(len(singular) - 1, max(-1, len(singular) - 4), -1), start=1):
        vector = right[singular_index]
        for dimension, loading in enumerate(vector):
            rows.append(
                {
                    "weak_order": order,
                    "singular_value_index": singular_index,
                    "singular_value": float(singular[singular_index]),
                    "dimension": dimension,
                    "parameter": str(domain[dimension]["name"]),
                    "loading": float(loading),
                    "squared_loading_share": float(loading**2),
                }
            )
    return rows


def build_refinement(
    coordinate_results_dir: Path,
    outdir: Path,
    expectations: Expectations,
) -> dict[str, Any]:
    outdir = outdir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise ContractError(f"Refusing to overwrite nonempty output directory: {outdir}")
    panel = load_coordinate_panel(coordinate_results_dir, expectations)
    dimensions = len(panel.domain)
    expected_tasks = 1 + 2 * dimensions
    jacobian, side_rows, frozen = central_jacobian(panel)
    full_rank = matrix_rank_report(jacobian)
    rank_key = f"relative_{IDENTIFICATION_THRESHOLD:g}"
    if int(full_rank["relative_ranks"][rank_key]) != dimensions:
        raise ContractError(
            "The weighted Jacobian is rank deficient at the 1e-3 "
            "relative threshold; no refinement plan was created"
        )
    active_rank = matrix_rank_report(jacobian[:, ~frozen])

    anchor_task = panel.tasks[0]
    incumbent_task = min(panel.tasks, key=lambda task: task.loss)
    formula_proposals = ridge_trust_proposals(
        anchor_task.residual, jacobian, anchor_task.unit, frozen
    )
    proposal_specs: list[dict[str, Any]] = [
        {
            "kind": "observed_incumbent",
            "observed_source_task": incumbent_task.task_id,
            "observed_loss": incumbent_task.loss,
            "unit_vector": incumbent_task.unit.tolist(),
            "predicted_loss": None,
            "step_vector": (incumbent_task.unit - anchor_task.unit).tolist(),
            "step_norm": float(np.linalg.norm(incumbent_task.unit - anchor_task.unit)),
        }
    ] + formula_proposals

    # No output is touched until the complete input, rank, and proposal gates pass.
    outdir.mkdir(parents=True, exist_ok=True)
    candidates_dir = outdir / "candidates"
    candidates_dir.mkdir()

    candidate_records: list[dict[str, Any]] = []
    for candidate_id, specification in enumerate(proposal_specs, start=1):
        unit = np.asarray(specification["unit_vector"], dtype=float)
        construction = {key: value for key, value in specification.items() if key != "unit_vector"}
        payload = center_payload(panel, unit, candidate_id, construction)
        relative_path = Path("candidates") / f"candidate_{candidate_id:03d}.json"
        candidate_path = outdir / relative_path
        atomic_write_json(candidate_path, payload)
        candidate_records.append(
            {
                "candidate_id": candidate_id,
                "label": (
                    "observed_incumbent"
                    if candidate_id == 1
                    else f"ridge_{specification['ridge_relative']:g}_trust_{specification['trust_radius']:g}"
                ),
                "kind": specification["kind"],
                "candidate_file": str(relative_path),
                "candidate_sha256": file_sha256(candidate_path),
                "unit_vector": unit.tolist(),
                "step_vector": specification["step_vector"],
                "step_norm": specification["step_norm"],
                "predicted_loss": specification.get("predicted_loss"),
                "observed_source_task": specification.get("observed_source_task"),
                "observed_loss": specification.get("observed_loss"),
                "ridge_relative": specification.get("ridge_relative"),
                "ridge_absolute": specification.get("ridge_absolute"),
                "trust_radius": specification.get("trust_radius"),
                "applied_scale": specification.get("applied_scale"),
                "applied_scale_reason": specification.get("applied_scale_reason"),
            }
        )

    jacobian_rows = []
    for moment_index, moment in enumerate(panel.moment_names):
        row: dict[str, Any] = {
            "moment_index": moment_index,
            "moment": moment,
            "anchor_weighted_residual": float(anchor_task.residual[moment_index]),
        }
        for dimension, domain_row in enumerate(panel.domain):
            row[str(domain_row["name"])] = float(jacobian[moment_index, dimension])
        jacobian_rows.append(row)
    singular_rows = [
        {
            "singular_value_index": index,
            "singular_value": value,
            "relative_to_largest": value / full_rank["singular_values"][0],
        }
        for index, value in enumerate(full_rank["singular_values"])
    ]
    atomic_write_csv(outdir / "weighted_jacobian.csv", jacobian_rows)
    atomic_write_csv(outdir / "singular_values.csv", singular_rows)
    atomic_write_csv(outdir / "side_consistency.csv", side_rows)
    atomic_write_csv(outdir / "weak_directions.csv", weak_direction_rows(jacobian, panel.domain))

    target_contract = [
        {
            "moment_index": index,
            "moment": moment,
            "target": float(panel.targets[index]),
            "weight": float(panel.weights[index]),
        }
        for index, moment in enumerate(panel.moment_names)
    ]
    plan = {
        "schema": "e5f_transition_ridge_refinement_plan_v1",
        "status": "ready_for_independent_scientific_validation",
        "refinement_builder": str(Path(__file__).resolve()),
        "refinement_builder_sha256": file_sha256(Path(__file__).resolve()),
        "coordinate_results_dir": str(panel.results_dir),
        "coordinate_panel_contract": {
            "task_count": expected_tasks,
            "dimensions": dimensions,
            "central_step": CENTRAL_STEP,
            "panel_seed": expectations.panel_seed,
            "center_sha256": panel.tasks[0].panel["center_sha256"],
            "center_role": (
                "prior_screen_warm_start_only; never a prior estimate; task_001 "
                "is freshly solved and rescored under this plan's exact contract"
            ),
            "input_artifact_sha256": panel.input_artifact_sha256,
        },
        "scientific_contract": panel.contract,
        "scientific_contract_sha256": panel.contract_sha256,
        "renewal_contract_sha256": object_sha256(
            panel.contract["renewal_accounting_contract"]
        ),
        "dated_contract_sha256": object_sha256(
            panel.contract["dated_measurement_contract"]
        ),
        "target_contract": target_contract,
        "target_contract_sha256": object_sha256(target_contract),
        "domain": panel.domain,
        "anchor": {
            "task_id": anchor_task.task_id,
            "status": "fresh_current_contract_solve_and_score",
            "loss": anchor_task.loss,
            "unit_vector": anchor_task.unit.tolist(),
            "weighted_residual": anchor_task.residual.tolist(),
        },
        "observed_incumbent": {
            "task_id": incumbent_task.task_id,
            "loss": incumbent_task.loss,
            "unit_vector": incumbent_task.unit.tolist(),
        },
        "jacobian_diagnostics": {
            "definition": "J[:,j]=(r(u+h e_j)-r(u-h e_j))/(2h), r=sqrt(weight)*(model-target)",
            "full": full_rank,
            "active_after_side_gate": active_rank,
            "identification_gate": {
                "relative_singular_value_threshold": IDENTIFICATION_THRESHOLD,
                "required_full_rank": dimensions,
                "passed": True,
            },
            "side_consistency_rule": "freeze a step column iff forward dot backward is nonpositive or either side has zero norm",
            "frozen_dimensions": np.flatnonzero(frozen).tolist(),
            "active_dimensions": np.flatnonzero(~frozen).tolist(),
            "side_consistency": side_rows,
        },
        "proposal_contract": {
            "ridge_relative_grid": list(RIDGE_RELATIVE_GRID),
            "ridge_absolute_formula": "lambda_abs=lambda_relative*sigma_max(J_active)^2",
            "trust_radius_grid": list(TRUST_RADIUS_GRID),
            "trust_norm": "Euclidean norm in transformed unit coordinates",
            "unit_box_rule": "scale the complete step along its ray before adding it to the anchor",
            "formula_candidate_count": 6,
            "observed_incumbent_count": 1,
            "total_candidate_count": 7,
            "no_automatic_second_round": True,
        },
        "validation_contract": {
            "scientific_driver": "code/model/tools/run_e5f_transition_calibration.py",
            "driver_mode": "one-candidate anchor: panel-task-id=1, panel-size=1, panel-design=mixed",
            "fresh_validation_tasks": 7,
            "fresh_identity_repeats_of_selected_candidate": 2,
            "eligibility_tolerances": ELIGIBILITY_TOLERANCES,
            "repeat_tolerances": {
                "loss_abs": 1.0e-9,
                "model_moment_abs": 1.0e-10,
                "theta_abs": 5.0e-12,
                "diagnostic_abs": 1.0e-10,
            },
        },
        "candidates": candidate_records,
    }
    plan_path = outdir / "candidate_plan.json"
    atomic_write_json(plan_path, plan)
    plan_sha = file_sha256(plan_path)
    atomic_write_text(outdir / "candidate_plan.sha256", f"{plan_sha}  candidate_plan.json\n")
    candidate_manifest = [
        f"{row['candidate_sha256']}  {row['candidate_file']}" for row in candidate_records
    ]
    atomic_write_text(outdir / "candidates.sha256", "\n".join(candidate_manifest) + "\n")
    return {
        "plan": plan,
        "plan_path": str(plan_path),
        "plan_sha256": plan_sha,
    }


def main() -> None:
    args = parse_args()
    result = build_refinement(
        args.coordinate_results_dir,
        args.outdir,
        Expectations(
            source_sha256=str(args.expected_source_sha256),
            source=str(args.expected_source) if args.expected_source is not None else None,
            target_set=str(args.expected_target_set),
            target_fingerprint=str(args.expected_target_fingerprint),
            code_bundle_sha256=str(args.expected_code_bundle_sha256),
            model_profile=str(args.expected_model_profile),
            panel_seed=int(args.expected_panel_seed),
            dimensions=int(args.expected_dimensions),
            renewal_contract_sha256=(
                str(args.expected_renewal_contract_sha256)
                if args.expected_renewal_contract_sha256 is not None
                else None
            ),
            dated_contract_sha256=(
                str(args.expected_dated_contract_sha256)
                if args.expected_dated_contract_sha256 is not None
                else None
            ),
        ),
    )
    plan = result["plan"]
    rank = plan["jacobian_diagnostics"]["full"]["relative_ranks"]
    frozen = plan["jacobian_diagnostics"]["frozen_dimensions"]
    print(
        "E5F_RIDGE_REFINEMENT_BUILT "
        f"plan_sha256={result['plan_sha256']} candidates=7 "
        f"rank_1e-3={rank['relative_0.001']}/"
        f"{plan['coordinate_panel_contract']['dimensions']} frozen={frozen}",
        flush=True,
    )


if __name__ == "__main__":
    try:
        main()
    except ContractError as error:
        print(f"E5F_RIDGE_REFINEMENT_REJECTED {error}", file=sys.stderr, flush=True)
        raise SystemExit(2) from error
