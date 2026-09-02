#!/usr/bin/env python3
"""Paired no-policy continuations from a certified 2023 transition state.

The driver deliberately does one narrow job.  It reconstructs the selected
2007--2023 transition calibration, verifies that reconstruction against the
selected case, and then branches from the same 2023 pre-fertility household
distribution and birth-vintage queue.  The two post-2023 population laws are

* a national closed benchmark, ``M=0`` and ``rho=1``; and
* an open/CBSA sensitivity using the selected report's externally fixed
  outside-origin share.

There are no policy or fiscal changes.  Both population laws retain the
sequential household solver, the date-by-date housing-market clearing rule,
and the top-code-adjusted births divided by 2.1 renewal accounting used in the
transition calibration.  A closed finite-horizon path is always a legitimate
benchmark.  It is called a transition between steady states only if the
post-shock stationary schedule contains a verified positive-price root of
``births / (2.1 * entrants) = 1``.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

import audit_closed_reproductive_closure as closure
import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_transition_calibration as calibration

from intergen_eqscale_seq_optimized.e5_profile import e5_target_system


START_YEAR = 2007
TERMINAL_YEAR = 2023
TRANSITION_PERIODS = 4
QUEUE_WAITING_SLOTS = 4
DEFAULT_OUTPUT = ROOT / "output/model/e5f_post2023_no_policy_continuations"
RECONSTRUCTION_TOLERANCE = 5e-10
ACCOUNTING_TOLERANCE = 2e-10
ROOT_RESIDUAL_TOLERANCE = 2.5e-5

HISTORY_STATE_COLUMNS = (
    "psi_child",
    "asset_price",
    "asset_price_index",
    "adult_population",
    "population_index",
    "entry_flow_E",
    "birth_children",
    "birth_children_topcode_adjusted",
    "effective_mature_entrant_flow_B",
    "housing_demand",
    "housing_supply",
    "owner_rate",
    "dependent_child_owner_rate",
    "parent_owner_rate",
    "childless_owner_rate",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-report", type=Path, required=True)
    parser.add_argument("--selected-task-summary", type=Path)
    parser.add_argument("--selected-case-dir", type=Path)
    parser.add_argument(
        "--selected-case-transition",
        type=Path,
        help=(
            "Selected transition_path.csv. For a ridge report this defaults to "
            "the report's copied selected_transition_path.csv."
        ),
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument("--expected-task-summary-sha256")
    parser.add_argument("--expected-case-transition-sha256", required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-target-fingerprint", required=True)
    parser.add_argument("--expected-code-bundle-sha256", required=True)
    parser.add_argument("--expected-renewal-contract-sha256", required=True)
    parser.add_argument("--expected-scientific-contract-sha256")
    parser.add_argument("--expected-selection-sha256")
    parser.add_argument("--post-2023-periods", type=int, default=40)
    parser.add_argument("--market-tol", type=float, default=2e-4)
    parser.add_argument("--market-max-iter", type=int, default=30)
    parser.add_argument("--closed-price-min-ratio", type=float, default=1e-4)
    parser.add_argument("--closed-price-max-ratio", type=float, default=3.0)
    parser.add_argument("--closed-grid-points", type=int, default=25)
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(jsonable(rows))
    temporary.replace(path)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        jsonable(payload), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def array_sha256(array: np.ndarray) -> str:
    values = np.ascontiguousarray(np.asarray(array, dtype="<f8"))
    digest = hashlib.sha256()
    digest.update(json.dumps(list(values.shape), separators=(",", ":")).encode("ascii"))
    digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def require_hash(path: Path, expected: str, label: str) -> str:
    actual = file_sha256(path)
    if actual != str(expected):
        raise RuntimeError(
            f"{label} hash gate failed: actual={actual}, expected={expected}, path={path}"
        )
    return actual


def _selected(report: dict[str, Any], key: str) -> Any:
    if key in report:
        return report[key]
    best = report.get("best_candidate", {})
    if key in best:
        return best[key]
    raise KeyError(f"Selected report does not contain {key!r}")


def required_renewal_contract() -> dict[str, Any]:
    return {
        "replacement_fertility": transition.REPLACEMENT_FERTILITY,
        "effective_birth_to_household_conversion": (
            transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
        ),
        "birth_measure": "topcode_adjusted_birth_children",
        "top_bin_mean_children": 3.602359422009,
        "birth_vintage_queue_waiting_slots": QUEUE_WAITING_SLOTS,
        "birth_to_entry_effect_lag_dates": QUEUE_WAITING_SLOTS + 1,
        "birth_to_entry_effect_lag_years": 20.0,
        "outside_flow_formula": "M = outside_origin_entry_share * E_old",
        "retention_formula": (
            "rho = (1 - outside_origin_entry_share) * E_old / B_old"
        ),
        "status": (
            "external replacement normalization; household KFE maturation "
            "is diagnostic and does not determine population renewal"
        ),
    }


def validate_renewal_contract(contract: dict[str, Any]) -> None:
    required = required_renewal_contract()
    if set(contract) != set(required):
        raise RuntimeError(
            "Renewal contract field gate failed: "
            f"actual={sorted(contract)}, required={sorted(required)}"
        )
    for key, expected in required.items():
        actual = contract[key]
        if isinstance(expected, float):
            if not math.isclose(
                float(actual), expected, rel_tol=0.0, abs_tol=1e-14
            ):
                raise RuntimeError(
                    f"Renewal contract gate failed for {key}: "
                    f"actual={actual}, required={expected}"
                )
        elif actual != expected:
            raise RuntimeError(
                f"Renewal contract gate failed for {key}: "
                f"actual={actual!r}, required={expected!r}"
            )


def compare_selected_candidate(
    report_candidate: dict[str, Any], task_candidate: dict[str, Any]
) -> None:
    required = (
        "candidate",
        "theta",
        "new_psi_child",
        "transition_loss",
        "policy_case",
        "post_2023_periods",
    )
    for key in required:
        if key not in report_candidate:
            raise RuntimeError(f"Collected best candidate is missing {key!r}")
        task_value = task_candidate.get(key)
        report_value = report_candidate[key]
        if isinstance(report_value, (float, int)) and not isinstance(
            report_value, bool
        ):
            if not math.isclose(
                float(report_value), float(task_value), rel_tol=0.0, abs_tol=1e-13
            ):
                raise RuntimeError(
                    f"Report/task candidate mismatch for {key}: "
                    f"report={report_value}, task={task_value}"
                )
        elif report_value != task_value:
            raise RuntimeError(
                f"Report/task candidate mismatch for {key}: "
                f"report={report_value!r}, task={task_value!r}"
            )


def validate_input_contracts(
    *,
    report_path: Path,
    task_summary_path: Path | None,
    case_dir: Path | None,
    case_transition_path: Path | None,
    source_path: Path,
    expected_report_sha256: str,
    expected_task_summary_sha256: str | None,
    expected_case_transition_sha256: str,
    expected_source_sha256: str,
    expected_target_fingerprint: str,
    expected_code_bundle_sha256: str,
    expected_renewal_contract_sha256: str,
    expected_scientific_contract_sha256: str | None,
    expected_selection_sha256: str | None,
    chain: Any,
    model: Any,
) -> dict[str, Any]:
    report_path = report_path.resolve()
    source_path = source_path.resolve()
    raw_report = json.loads(report_path.read_text(encoding="utf-8"))
    ridge_schema = raw_report.get("schema") == "e5f_transition_ridge_refinement_report_v1"
    if ridge_schema:
        if task_summary_path is not None or case_dir is not None:
            raise RuntimeError(
                "A ridge report must use its copied selected transition, not a "
                "separate task summary or case directory"
            )
        case_transition_path = (
            case_transition_path.resolve()
            if case_transition_path is not None
            else report_path.parent / "selected_transition_path.csv"
        )
    else:
        if task_summary_path is None or case_dir is None:
            raise RuntimeError(
                "A collector-style report requires --selected-task-summary and "
                "--selected-case-dir"
            )
        task_summary_path = task_summary_path.resolve()
        case_dir = case_dir.resolve()
        inferred = case_dir / "transition_path.csv"
        if case_transition_path is not None and case_transition_path.resolve() != inferred:
            raise RuntimeError("Collector-style case-transition path differs from case directory")
        case_transition_path = inferred
    required_paths = [report_path, case_transition_path, source_path]
    if task_summary_path is not None:
        required_paths.append(task_summary_path)
    for path in required_paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    hashes = {
        "selected_report_sha256": require_hash(
            report_path, expected_report_sha256, "Selected report"
        ),
        "selected_case_transition_sha256": require_hash(
            case_transition_path,
            expected_case_transition_sha256,
            "Selected case transition",
        ),
        "source_sha256": require_hash(
            source_path, expected_source_sha256, "Stationary source"
        ),
    }
    if ridge_schema:
        if expected_task_summary_sha256 is not None:
            raise RuntimeError("A ridge report must not declare a task-summary hash")
        if (
            expected_scientific_contract_sha256 is None
            or expected_selection_sha256 is None
        ):
            raise RuntimeError(
                "A ridge report requires expected scientific-contract and selection hashes"
            )
        if raw_report.get("status") != "complete_refinement_with_two_independent_identity_repeats":
            raise RuntimeError(
                f"Ridge report is not complete: status={raw_report.get('status')!r}"
            )
        if not bool(raw_report.get("promotion_eligible", False)):
            raise RuntimeError("Ridge report is not promotion eligible")
        repeat_gate = dict(raw_report.get("repeat_gate") or {})
        if not (
            int(repeat_gate.get("required", -1)) == 2
            and int(repeat_gate.get("completed", -1)) == 2
            and bool(repeat_gate.get("identity_gate_passed", False))
            and bool(repeat_gate.get("numerical_identity_gate_passed", False))
        ):
            raise RuntimeError(f"Ridge repeat gate failed: {repeat_gate}")
        scientific = dict(raw_report.get("scientific_contract") or {})
        scientific_sha = canonical_json_sha256(scientific)
        if not (
            scientific_sha
            == str(raw_report.get("scientific_contract_sha256"))
            == str(expected_scientific_contract_sha256)
        ):
            raise RuntimeError(
                "Scientific-contract hash gate failed: "
                f"actual={scientific_sha}, report={raw_report.get('scientific_contract_sha256')}, "
                f"expected={expected_scientific_contract_sha256}"
            )
        if str(raw_report.get("selection_sha256")) != str(expected_selection_sha256):
            raise RuntimeError(
                "Selection hash gate failed: "
                f"report={raw_report.get('selection_sha256')}, "
                f"expected={expected_selection_sha256}"
            )
        hashes["scientific_contract_sha256"] = scientific_sha
        hashes["selection_sha256"] = str(expected_selection_sha256)
        hashes["selected_candidate_sha256"] = str(
            raw_report.get("selected_candidate_sha256")
        )
        hashes["plan_sha256"] = str(raw_report.get("plan_sha256"))
        hashes["dated_contract_sha256"] = str(
            raw_report.get("dated_contract_sha256")
        )
        report_best = dict(raw_report.get("best_candidate") or {})
        if str(report_best.get("candidate_sha256")) != hashes[
            "selected_candidate_sha256"
        ]:
            raise RuntimeError(
                "Ridge selected-candidate hash differs between report fields"
            )
        report = dict(scientific)
        report.update(
            best_candidate=report_best,
            old_psi_child=report_best.get("old_psi_child"),
            renewal_accounting_old_state=report_best.get(
                "renewal_accounting_old_state"
            ),
        )
        task = copy.deepcopy(report)
        task_best = copy.deepcopy(report_best)
    else:
        if expected_task_summary_sha256 is None:
            raise RuntimeError("Collector-style report requires a task-summary hash")
        hashes["selected_task_summary_sha256"] = require_hash(
            task_summary_path,
            expected_task_summary_sha256,
            "Selected task summary",
        )
        report = raw_report
        task = json.loads(task_summary_path.read_text(encoding="utf-8"))
        if report.get("status") != "complete":
            raise RuntimeError(
                f"Selected report is not complete: status={report.get('status')!r}"
            )
        if not str(task.get("status", "")).startswith("complete_transition_calibration"):
            raise RuntimeError(
                f"Selected task is not complete: status={task.get('status')!r}"
            )
        report_best = report.get("best_candidate")
        task_best = task.get("best_candidate")
        if not isinstance(report_best, dict) or not isinstance(task_best, dict):
            raise RuntimeError("Selected artifacts do not contain best_candidate records")
        compare_selected_candidate(report_best, task_best)
        if str(report_best["candidate"]) != case_dir.name:
            raise RuntimeError(
                "Selected case directory does not match the collected winner: "
                f"case={case_dir.name}, winner={report_best['candidate']}"
            )
    if not isinstance(report_best, dict) or not isinstance(task_best, dict):
        raise RuntimeError("Selected artifacts do not contain best_candidate records")
    if not bool(report_best.get("valid", True)):
        raise RuntimeError("Collected best candidate is marked invalid")

    for artifact_name, artifact in (("report", report), ("task", task)):
        policy_case = str(_selected(artifact, "policy_case"))
        post_periods = int(_selected(artifact, "post_2023_periods"))
        if policy_case != "none" or post_periods != 0:
            raise RuntimeError(
                f"{artifact_name} scope gate failed: policy={policy_case}, "
                f"post_2023_periods={post_periods}"
            )

    report_target = str(_selected(report, "target_fingerprint"))
    task_target = str(_selected(task, "target_fingerprint"))
    live_target = e5_target_system()
    if not (
        report_target
        == task_target
        == live_target.fingerprint
        == str(expected_target_fingerprint)
    ):
        raise RuntimeError(
            "Target fingerprint gate failed: "
            f"report={report_target}, task={task_target}, "
            f"live={live_target.fingerprint}, expected={expected_target_fingerprint}"
        )

    report_old_psi = float(_selected(report, "old_psi_child"))
    task_old_psi = float(_selected(task, "old_psi_child"))
    if not math.isclose(report_old_psi, task_old_psi, rel_tol=0.0, abs_tol=1e-13):
        raise RuntimeError(
            "Report/task old preference intercept mismatch: "
            f"report={report_old_psi}, task={task_old_psi}"
        )

    report_code = str(_selected(report, "code_fingerprints")["bundle_sha256"])
    task_code = str(_selected(task, "code_fingerprints")["bundle_sha256"])
    live_code_contract = calibration.code_fingerprint_contract(model)
    live_code = str(live_code_contract["bundle_sha256"])
    if not (
        report_code == task_code == live_code == str(expected_code_bundle_sha256)
    ):
        raise RuntimeError(
            "Code-bundle fingerprint gate failed: "
            f"report={report_code}, task={task_code}, live={live_code}, "
            f"expected={expected_code_bundle_sha256}"
        )

    report_renewal = dict(_selected(report, "renewal_accounting_contract"))
    task_renewal = dict(_selected(task, "renewal_accounting_contract"))
    validate_renewal_contract(report_renewal)
    validate_renewal_contract(task_renewal)
    if canonical_json_sha256(report_renewal) != canonical_json_sha256(task_renewal):
        raise RuntimeError("Report and selected task use different renewal contracts")
    renewal_sha = canonical_json_sha256(report_renewal)
    reported_renewal_sha = (
        str(raw_report.get("renewal_contract_sha256")) if ridge_schema else renewal_sha
    )
    if not (
        renewal_sha
        == reported_renewal_sha
        == str(expected_renewal_contract_sha256)
    ):
        raise RuntimeError(
            "Renewal-contract hash gate failed: "
            f"actual={renewal_sha}, report={reported_renewal_sha}, "
            f"expected={expected_renewal_contract_sha256}"
        )
    hashes["renewal_contract_sha256"] = renewal_sha
    hashes["target_fingerprint"] = live_target.fingerprint
    hashes["code_bundle_sha256"] = live_code

    report_source = str(_selected(report, "source_sha256"))
    task_source = str(_selected(task, "source_sha256"))
    if not (
        report_source
        == task_source
        == hashes["source_sha256"]
        == str(expected_source_sha256)
    ):
        raise RuntimeError(
            "Stationary-source provenance gate failed: "
            f"report={report_source}, task={task_source}, "
            f"actual={hashes['source_sha256']}, expected={expected_source_sha256}"
        )

    report_old = dict(_selected(report, "renewal_accounting_old_state"))
    task_old = dict(_selected(task, "renewal_accounting_old_state"))
    for key in (
        "old_entry_flow_E",
        "old_queue_mature_flow_B",
        "outside_flow_M",
        "outside_origin_entry_share",
        "retention_rho",
    ):
        if not math.isclose(
            float(report_old[key]),
            float(task_old[key]),
            rel_tol=0.0,
            abs_tol=1e-13,
        ):
            raise RuntimeError(f"Report/task old-renewal mismatch for {key}")
    outside_share = float(_selected(report, "outside_origin_entry_share"))
    if not 0.0 < outside_share < 1.0:
        raise RuntimeError(f"Outside-origin share is invalid: {outside_share}")
    if not math.isclose(
        outside_share,
        float(report_old["outside_origin_entry_share"]),
        rel_tol=0.0,
        abs_tol=1e-14,
    ):
        raise RuntimeError("Outside-origin share is inconsistent within the report")
    E_old = float(report_old["old_entry_flow_E"])
    B_old = float(report_old["old_queue_mature_flow_B"])
    M_old = float(report_old["outside_flow_M"])
    rho_old = float(report_old["retention_rho"])
    if E_old <= 0.0 or B_old <= 0.0 or not 0.0 <= rho_old <= 1.0:
        raise RuntimeError("Selected old renewal flows or retention are invalid")
    formula_gaps = {
        "outside_flow": M_old - outside_share * E_old,
        "retention": rho_old - (1.0 - outside_share) * E_old / B_old,
        "renewal_identity": E_old - M_old - rho_old * B_old,
    }
    if max(abs(value) for value in formula_gaps.values()) > 2e-12:
        raise RuntimeError(
            f"Selected old renewal accounting fails its formulas: {formula_gaps}"
        )

    return {
        "report": report,
        "raw_report": raw_report,
        "task": task,
        "report_best": report_best,
        "task_best": task_best,
        "renewal_contract": report_renewal,
        "renewal_old_state": report_old,
        "outside_origin_entry_share": outside_share,
        "hashes": hashes,
        "paths": {
            "selected_report": str(report_path),
            "selected_task_summary": (
                str(task_summary_path) if task_summary_path is not None else None
            ),
            "selected_case_dir": str(case_dir) if case_dir is not None else None,
            "selected_case_transition": str(case_transition_path),
            "source": str(source_path),
        },
    }


def outside_share_invariance_audit(
    E_old: float, B_old: float, selected_share: float = 0.169
) -> dict[str, Any]:
    if E_old <= 0.0 or B_old <= 0.0:
        raise ValueError("Old-state entry and mature flow must be positive")
    rows = []
    shares = sorted(set((0.0, float(selected_share), 0.5, 0.9)))
    for share in shares:
        if not 0.0 <= share < 1.0:
            raise ValueError(f"Outside-origin share is invalid: {share}")
        outside_flow = share * E_old
        retention = (1.0 - share) * E_old / B_old
        implied_entry = outside_flow + retention * B_old
        rows.append(
            {
                "outside_origin_entry_share": share,
                "outside_flow_M": outside_flow,
                "retention_rho": retention,
                "implied_entry": implied_entry,
                "identity_residual": implied_entry - E_old,
            }
        )
    maximum = max(abs(float(row["identity_residual"])) for row in rows)
    if maximum > ACCOUNTING_TOLERANCE:
        raise RuntimeError(
            f"Outside-share invariance accounting failed: max residual={maximum:.3e}"
        )
    return {
        "status": "passed",
        "maximum_absolute_identity_residual": maximum,
        "rows": rows,
        "interpretation": (
            "The outside-origin share only decomposes old steady-state entry "
            "into M and rho*B. It is externally fixed and is not identified by "
            "the 2007--2023 transition, whose pre-2023 entrant cohorts come from "
            "the initialized birth-vintage queue."
        ),
    }


@dataclass
class DynamicState:
    g_pre: np.ndarray
    scheduled_entries: list[float]
    scheduled_raw_entries: list[float]
    price_guess: float
    initial_policy: Any | None = None


@dataclass
class PreparedModel:
    chain: Any
    model: Any
    base_overrides: dict[str, Any]
    parameters: SimpleNamespace
    b_grid: np.ndarray
    old_price: float
    old_policy: Any
    initial_g_pre_2007: np.ndarray
    initial_mass_2007: float
    births_old: float
    births_old_raw: float
    E_old: float
    B_old: float
    supply_rule: calendar.HousingSupplyRule
    supply_normalization: dict[str, Any]
    operator_gates: dict[str, Any]
    stationary_gates: dict[str, Any]


def prepare_model(
    contracts: dict[str, Any], source_path: Path, chain: Any, model: Any
) -> PreparedModel:
    best = contracts["report_best"]
    theta = {str(key): float(value) for key, value in best["theta"].items()}
    source_winner, _ = closure.load_winner(source_path)
    source_theta = closure.theta_from_winner(source_winner)
    missing_source_coordinates = sorted(set(source_theta) - set(theta))
    if missing_source_coordinates:
        raise RuntimeError(
            "Selected transition estimate omits source coordinates: "
            f"{missing_source_coordinates}"
        )
    profile_name = str(_selected(contracts["report"], "model_profile")["name"])
    _, profile_overrides, _ = calibration.activate_model_profile(profile_name, theta)
    base = closure.make_overrides(chain, theta, nb=120, profile="e5f-floor")
    base.update(profile_overrides)
    base.update(theta)
    positive_tenure_dispersion = float(theta.get("tenure_choice_kappa", 0.0)) > 0.0
    if positive_tenure_dispersion:
        # Match the audited transition-calibration numerical contract. Positive
        # tenure smoothing creates only solver-scale redistribution roundoff;
        # normalize it after the same strict mass gate used in calibration.
        base["normalize_transition_mass_roundoff"] = True
    base.update(
        parent_dp_waiver=False,
        parent_dp_waiver_phi=1.0,
        parent_dp_waiver_locations=np.array([], dtype=int),
        parent_dp_waiver_owner_rungs=np.array([], dtype=int),
        parent_dp_waiver_birth_state_only=False,
    )
    old_psi = float(best["old_psi_child"])
    old_base = dict(base)
    old_base["psi_child"] = old_psi
    old_solution, old_parameters, old_price_array, _ = closure.solve_ge(chain, old_base)
    old_price = float(old_price_array[0])
    old_moments = chain.extract_moments(old_solution, old_parameters)
    old_completed = float(old_moments["tfr"])
    reported_completed = float(best["old_completed_fertility"])
    if not math.isclose(
        old_completed, reported_completed, rel_tol=0.0, abs_tol=2e-8
    ):
        raise RuntimeError(
            "Old steady-state reproduction failed to reproduce the selected case: "
            f"solved={old_completed}, selected={reported_completed}"
        )
    if abs(old_completed - transition.REPLACEMENT_FERTILITY) > 5e-4:
        raise RuntimeError(
            f"Selected old steady state is not at replacement: {old_completed}"
        )
    if int(old_parameters.J) != 17 or int(old_parameters.Nb) != 120:
        raise RuntimeError("Production dimension gate failed")

    b_grid = np.asarray(old_solution.b_grid, dtype=float)
    shared = model.precompute_shared(old_parameters, b_grid)
    old_parameters._fert2_probs = np.asarray(
        old_solution.fert2_probs, dtype=float
    ).copy()
    old_policy = calendar.policy_from_solution(
        old_solution, old_price_array, old_parameters, b_grid, shared
    )
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = (
        transition.advance_sequential_calendar_distribution
    )
    calendar.distribution_rows = transition.independent_child_distribution_rows
    stationary_g_pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(
        old_solution, old_policy, old_parameters, b_grid, shared
    )
    operator_gates = transition.operator_gates(
        old_solution,
        old_policy,
        stationary_g_pre,
        old_parameters,
        b_grid,
        shared,
    )
    operator_gates.update(reconstruction)
    gate_names = (
        "stationary_post_fertility_nesting_l1",
        "one_step_constant_path_nesting_l1",
        "mature_flow_abs_error",
        "birth_flow_abs_error",
        "topcode_adjusted_birth_flow_abs_error",
    )
    if max(abs(float(operator_gates[name])) for name in gate_names) > 5e-9:
        raise RuntimeError(f"Stationary transition-operator gates failed: {operator_gates}")

    renewal = closure.topcode_consistent_renewal_accounting(
        old_solution, old_parameters
    )
    conversion = transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
    births_old = float(renewal["topcode_adjusted_birth_children"])
    births_old_raw = float(renewal["raw_birth_children"])
    E_old = float(old_solution.entry_rate)
    B_old = conversion * births_old
    reported_old = contracts["renewal_old_state"]
    stationary_gates = {
        "old_completed_fertility": old_completed,
        "old_entry_flow_E": E_old,
        "old_queue_mature_flow_B": B_old,
        "old_queue_B_over_E": B_old / E_old,
        "old_price": old_price,
        "reported_old_entry_flow_E": float(reported_old["old_entry_flow_E"]),
        "reported_old_queue_mature_flow_B": float(
            reported_old["old_queue_mature_flow_B"]
        ),
    }
    if (
        abs(E_old - stationary_gates["reported_old_entry_flow_E"]) > 2e-10
        or abs(B_old - stationary_gates["reported_old_queue_mature_flow_B"]) > 2e-10
    ):
        raise RuntimeError(
            "Re-solved old renewal state does not match the selected report: "
            f"{stationary_gates}"
        )

    ages = (
        float(old_parameters.age_start)
        + np.arange(int(old_parameters.J), dtype=float) * float(old_parameters.da)
    )
    stationary_age_mass = np.sum(
        stationary_g_pre, axis=(0, 1, 2, 4, 5, 6)
    )
    structural_survival = None
    if positive_tenure_dispersion:
        structural_survival = (
            np.asarray(old_parameters.survival_probs, dtype=float)[
                : int(old_parameters.J) - 1
            ]
            if bool(old_parameters.use_age_survival)
            else np.ones(int(old_parameters.J) - 1, dtype=float)
        )
        calibration.validated_structural_stationary_age_mass(
            stationary_age_mass,
            entry_flow=E_old,
            structural_survival=structural_survival,
        )
    age_reweight = transition.acs_2007_age_reweight_diagnostic(
        stationary_age_mass,
        ages,
        E_old,
        periods=TRANSITION_PERIODS,
        period_years=float(old_parameters.period_years),
        structural_survival=structural_survival,
    )
    initial_g_pre = transition.reweight_distribution_to_acs_2007_ages(
        stationary_g_pre, age_reweight
    )
    supply_rule, supply_normalization = calendar.normalize_date0_housing_supply(
        initial_g_pre,
        old_policy,
        old_parameters,
        b_grid,
        shared,
        "static-elastic",
    )
    try:
        external_closure = _selected(
            contracts["report"], "external_closure_contract"
        )
    except KeyError:
        external_closure = None
    if external_closure is not None:
        external_closure = dict(external_closure)
        selected_elasticity = float(
            external_closure["housing_supply_elasticity"]
        )
        if (
            not math.isfinite(selected_elasticity)
            or selected_elasticity <= 0.0
            or str(external_closure.get("housing_supply_elasticity_status"))
            != "externally_fixed_profile_not_estimated"
        ):
            raise RuntimeError("Selected external supply elasticity is invalid")
        retained_elasticity = float(supply_rule.elasticity)
        supply_rule = calendar.HousingSupplyRule(
            mode=supply_rule.mode,
            initial_price=supply_rule.initial_price,
            initial_stock=supply_rule.initial_stock,
            elasticity=selected_elasticity,
        )
        supply_normalization["retained_housing_supply_elasticity"] = (
            retained_elasticity
        )
        supply_normalization["transition_housing_supply_elasticity"] = (
            selected_elasticity
        )
        supply_normalization["elasticity_status"] = (
            "externally_fixed_profile_not_estimated"
        )
    return PreparedModel(
        chain=chain,
        model=model,
        base_overrides=base,
        parameters=old_parameters,
        b_grid=b_grid,
        old_price=old_price,
        old_policy=old_policy,
        initial_g_pre_2007=initial_g_pre,
        initial_mass_2007=float(np.sum(initial_g_pre)),
        births_old=births_old,
        births_old_raw=births_old_raw,
        E_old=E_old,
        B_old=B_old,
        supply_rule=supply_rule,
        supply_normalization=supply_normalization,
        operator_gates=operator_gates,
        stationary_gates=stationary_gates,
    )


def clear_market(
    state: DynamicState,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    counter: calendar.SolveCounter,
    supply_rule: calendar.HousingSupplyRule,
    market_tol: float,
    market_max_iter: int,
) -> tuple[calendar.PeriodEvaluation, bool]:
    target_tolerance = min(float(market_tol), 5e-5)
    try:
        evaluation = calendar.clear_scalar_housing_market(
            state.g_pre,
            state.price_guess,
            P,
            b_grid,
            shared,
            counter,
            target_tolerance,
            min(int(market_max_iter), 18),
            supply_rule,
            initial_policy=state.initial_policy,
        )
        return evaluation, False
    except RuntimeError as error:
        if (
            float(market_tol) <= target_tolerance
            or "Housing market did not clear" not in str(error)
        ):
            raise
        evaluation = calendar.clear_scalar_housing_market(
            state.g_pre,
            state.price_guess,
            P,
            b_grid,
            shared,
            counter,
            float(market_tol),
            int(market_max_iter),
            supply_rule,
            initial_policy=state.initial_policy,
        )
        return evaluation, True


def evaluate_state(
    state: DynamicState,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    supply_rule: calendar.HousingSupplyRule,
    counter: calendar.SolveCounter,
    market_tol: float,
    market_max_iter: int,
) -> tuple[calendar.PeriodEvaluation, SimpleNamespace, bool]:
    shared = calendar.model.precompute_shared(P, b_grid)
    evaluation, fallback = clear_market(
        state,
        P,
        b_grid,
        shared,
        counter,
        supply_rule,
        market_tol,
        market_max_iter,
    )
    return evaluation, shared, fallback


def advance_from_evaluation(
    *,
    label: str,
    period_from_2007: int,
    evaluation: calendar.PeriodEvaluation,
    state: DynamicState,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    supply_rule: calendar.HousingSupplyRule,
    outside_flow: float,
    retention: float,
    initial_mass_2007: float,
    mass_2023: float,
    next_bridge_year: int | None,
    grid_fallback: bool,
) -> tuple[dict[str, Any], DynamicState]:
    conversion = transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
    empty_next, model_mature_raw_by_loc, deaths, _ = (
        transition.advance_sequential_calendar_distribution(
            evaluation, np.zeros(int(P.I)), P, b_grid, shared
        )
    )
    accounting = transition.calendar_topcode_birth_accounting(
        evaluation.g_pre,
        evaluation.g_post_fertility,
        float(evaluation.births),
        P,
    )
    adjusted_births = float(accounting["topcode_adjusted_birth_children"])
    scheduled_B, next_queue = transition.advance_birth_vintage_queue(
        state.scheduled_entries, adjusted_births, conversion
    )
    scheduled_raw_B, next_raw_queue = transition.advance_birth_vintage_queue(
        state.scheduled_raw_entries, float(evaluation.births), conversion
    )
    entrants_next_by_loc = np.asarray(P.entry_shares, dtype=float).reshape(-1)
    entrants_next_by_loc /= np.sum(entrants_next_by_loc)
    entrants_next_by_loc *= float(outside_flow)
    retained_B = float(retention) * scheduled_B
    entrants_next_by_loc[0] += retained_B
    empty_next[:, :, :, 0, :, :, :] = calendar.entrant_cohort(
        entrants_next_by_loc, P, b_grid
    )
    bridge_audit = None
    if next_bridge_year is not None:
        ages = float(P.age_start) + np.arange(int(P.J), dtype=float) * float(P.da)
        next_pre, bridge_audit = transition.reweight_distribution_to_observed_age_path(
            empty_next,
            ages,
            year=int(next_bridge_year),
            initial_mass=initial_mass_2007,
        )
    else:
        next_pre = empty_next

    health = calendar.distribution_health(
        {
            "pre_fertility": evaluation.g_pre,
            "post_fertility": evaluation.g_post_fertility,
            "current_choice": evaluation.g_current,
            "next_pre_fertility": next_pre,
        }
    )
    minimum = health["min_distribution_mass"]
    nonfinite = int(health["nonfinite_distribution_count"])
    if minimum is None or nonfinite > 0 or float(minimum) < -1e-14:
        raise RuntimeError(f"{label} distribution-health gate failed: {health}")
    expected_next_mass = (
        float(np.sum(evaluation.g_post_fertility))
        - deaths
        + float(np.sum(entrants_next_by_loc))
        + (float(bridge_audit["net_residual"]) if bridge_audit is not None else 0.0)
    )
    mass_residual = float(np.sum(next_pre)) - expected_next_mass
    if abs(mass_residual) > 2e-8:
        raise RuntimeError(
            f"{label} mass-accounting gate failed: residual={mass_residual:.3e}"
        )
    if float(evaluation.relative_market_residual) > 2e-4:
        raise RuntimeError(
            f"{label} market gate failed: residual={evaluation.relative_market_residual:.3e}"
        )
    if float(evaluation.feasibility_projection_mass) > 1e-6:
        raise RuntimeError(
            f"{label} feasibility-projection gate failed: "
            f"mass={evaluation.feasibility_projection_mass:.3e}"
        )

    adult_population = float(np.sum(evaluation.g_post_fertility))
    current_mass = float(np.sum(evaluation.g_current))
    dependent_mass = float(np.sum(evaluation.g_current[:, :, :, :, :, :, 1:]))
    parent_mass = float(np.sum(evaluation.g_current[:, :, :, :, :, 1:, :]))
    childless_mass = float(np.sum(evaluation.g_current[:, :, :, :, :, 0, :]))
    owner_mass = float(np.sum(evaluation.g_current[:, 1:, :, :, :, :, :]))
    dependent_owner_mass = float(
        np.sum(evaluation.g_current[:, 1:, :, :, :, :, 1:])
    )
    parent_owner_mass = float(
        np.sum(evaluation.g_current[:, 1:, :, :, :, 1:, :])
    )
    childless_owner_mass = float(
        np.sum(evaluation.g_current[:, 1:, :, :, :, 0, :])
    )
    entry_flow = float(np.sum(evaluation.g_pre[:, :, :, 0, :, :, :]))
    demand = float(np.sum(evaluation.demand_by_loc))
    supply = float(np.sum(evaluation.supply_by_loc))
    year = START_YEAR + int(period_from_2007) * int(P.period_years)
    row = {
        "scenario": label,
        "policy_case": "none",
        "period_from_2007": int(period_from_2007),
        "calendar_year": year,
        "years_from_2023": year - TERMINAL_YEAR,
        "psi_child": float(P.psi_child),
        "asset_price": float(evaluation.policy.price[0]),
        "asset_price_index": float(evaluation.policy.price[0])
        / float(supply_rule.initial_price),
        "adult_population": adult_population,
        "population_index": adult_population / max(initial_mass_2007, 1e-15),
        "population_index_2023": adult_population / max(mass_2023, 1e-15),
        "entry_flow_E": entry_flow,
        "birth_children": float(evaluation.births),
        "birth_children_topcode_adjusted": adjusted_births,
        "topcode_adjusted_births_per_adult": adjusted_births
        / max(adult_population, 1e-15),
        "effective_mature_entrant_flow_B": scheduled_B,
        "raw_state_scheduled_mature_entrant_flow_B": scheduled_raw_B,
        "queue_B_over_current_E": scheduled_B / max(entry_flow, 1e-15),
        "outside_entry_flow_M": float(outside_flow),
        "renewal_retention_rho": float(retention),
        "retained_mature_entrants": retained_B,
        "entrant_flow_next": float(np.sum(entrants_next_by_loc)),
        "adult_deaths": deaths,
        "housing_demand": demand,
        "housing_demand_per_adult": demand / max(adult_population, 1e-15),
        "housing_supply": supply,
        "relative_market_residual": float(evaluation.relative_market_residual),
        "mass_accounting_residual": mass_residual,
        "owner_rate": owner_mass / max(current_mass, 1e-15),
        "dependent_child_owner_rate": dependent_owner_mass
        / max(dependent_mass, 1e-15),
        "parent_owner_rate": parent_owner_mass / max(parent_mass, 1e-15),
        "childless_owner_rate": childless_owner_mass
        / max(childless_mass, 1e-15),
        "birth_queue_scheduled_flows": json.dumps(next_queue),
        "birth_queue_raw_state_scheduled_flows": json.dumps(next_raw_queue),
        "model_state_same_period_mature_flow_B": float(P.entrant_conversion_factor)
        * float(np.sum(model_mature_raw_by_loc)),
        "census_age_bridge_net_residual": (
            float(bridge_audit["net_residual"]) if bridge_audit is not None else 0.0
        ),
        "census_age_bridge_maximum_target_gap": (
            float(bridge_audit["maximum_absolute_target_gap"])
            if bridge_audit is not None
            else 0.0
        ),
        "grid_resolution_fallback": bool(grid_fallback),
        "feasibility_frontier_projection_mass": float(
            evaluation.feasibility_projection_mass
        ),
        "min_distribution_mass": float(minimum),
        "nonfinite_distribution_count": nonfinite,
    }
    next_state = DynamicState(
        g_pre=next_pre,
        scheduled_entries=list(next_queue),
        scheduled_raw_entries=list(next_raw_queue),
        price_guess=float(evaluation.policy.price[0]),
        initial_policy=evaluation.policy,
    )
    return row, next_state


def reconstruct_and_branch(
    prepared: PreparedModel,
    contracts: dict[str, Any],
    *,
    post_2023_periods: int,
    market_tol: float,
    market_max_iter: int,
    progress_dir: Path | None = None,
) -> tuple[
    list[dict[str, Any]],
    dict[str, list[dict[str, Any]]],
    dict[str, Any],
]:
    if post_2023_periods < 1:
        raise ValueError("At least one post-2023 period is required")
    best = contracts["report_best"]
    old_psi = float(best["old_psi_child"])
    new_psi = float(best["new_psi_child"])
    report_old = contracts["renewal_old_state"]
    outside_flow = float(report_old["outside_flow_M"])
    retention = float(report_old["retention_rho"])
    initial_queue = [prepared.B_old] * QUEUE_WAITING_SLOTS
    initial_raw_B = (
        transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
        * prepared.births_old_raw
    )
    state = DynamicState(
        g_pre=prepared.initial_g_pre_2007.copy(),
        scheduled_entries=initial_queue,
        scheduled_raw_entries=[initial_raw_B] * QUEUE_WAITING_SLOTS,
        price_guess=prepared.old_price,
        initial_policy=None,
    )
    P_history = copy.deepcopy(prepared.parameters)
    counter = calendar.SolveCounter()
    target_age_years = transition.census_household_age_target_years()
    common_history: list[dict[str, Any]] = []
    previous_psi: float | None = None
    for period in range(TRANSITION_PERIODS):
        dated_psi = transition.preference_shifter_at_date(
            period, old_psi, new_psi, TRANSITION_PERIODS
        )
        P_history.psi_child = dated_psi
        P_history.parent_dp_waiver = False
        P_history.parent_dp_waiver_phi = 1.0
        P_history.parent_dp_waiver_locations = np.array([], dtype=int)
        P_history.parent_dp_waiver_owner_rungs = np.array([], dtype=int)
        P_history.parent_dp_waiver_birth_state_only = False
        if previous_psi is not None and not math.isclose(
            dated_psi, previous_psi, rel_tol=0.0, abs_tol=1e-14
        ):
            state.initial_policy = None
        previous_psi = dated_psi
        evaluation, shared, fallback = evaluate_state(
            state,
            P_history,
            prepared.b_grid,
            prepared.supply_rule,
            counter,
            market_tol,
            market_max_iter,
        )
        row, state = advance_from_evaluation(
            label="shared_calibrated_history",
            period_from_2007=period,
            evaluation=evaluation,
            state=state,
            P=P_history,
            b_grid=prepared.b_grid,
            shared=shared,
            supply_rule=prepared.supply_rule,
            outside_flow=outside_flow,
            retention=retention,
            initial_mass_2007=prepared.initial_mass_2007,
            mass_2023=1.0,
            next_bridge_year=target_age_years[period + 1],
            grid_fallback=fallback,
        )
        common_history.append(row)
        if progress_dir is not None:
            write_csv(progress_dir / "continuation_progress.csv", common_history)
            write_json(
                progress_dir / "latest_completed_period.json",
                {
                    "status": "running",
                    "phase": "reconstructing_selected_2007_2023_history",
                    "calendar_year": row["calendar_year"],
                    "completed_model_dates": len(common_history),
                    "latest": row,
                },
            )

    state_2023 = DynamicState(
        g_pre=state.g_pre.copy(),
        scheduled_entries=list(state.scheduled_entries),
        scheduled_raw_entries=list(state.scheduled_raw_entries),
        price_guess=state.price_guess,
        initial_policy=None,
    )
    mass_2023 = float(np.sum(state_2023.g_pre))
    matched_state_sha256 = array_sha256(state_2023.g_pre)
    P_2023 = copy.deepcopy(P_history)
    P_2023.psi_child = new_psi
    shared_evaluation, shared_2023, fallback_2023 = evaluate_state(
        state_2023,
        P_2023,
        prepared.b_grid,
        prepared.supply_rule,
        counter,
        market_tol,
        market_max_iter,
    )

    closure_specs = {
        "closed_national_benchmark": {"outside_flow": 0.0, "retention": 1.0},
        "open_cbsa_sensitivity": {
            "outside_flow": outside_flow,
            "retention": retention,
        },
    }
    paths: dict[str, list[dict[str, Any]]] = {}
    branch_states: dict[str, DynamicState] = {}
    branch_parameters: dict[str, SimpleNamespace] = {}
    for label, spec in closure_specs.items():
        branch_state = DynamicState(
            g_pre=state_2023.g_pre.copy(),
            scheduled_entries=list(state_2023.scheduled_entries),
            scheduled_raw_entries=list(state_2023.scheduled_raw_entries),
            price_guess=state_2023.price_guess,
            initial_policy=None,
        )
        branch_P = copy.deepcopy(P_2023)
        if array_sha256(branch_state.g_pre) != matched_state_sha256:
            raise RuntimeError(f"{label} did not receive the exact matched 2023 state")
        row_2023, next_state = advance_from_evaluation(
            label=label,
            period_from_2007=TRANSITION_PERIODS,
            evaluation=shared_evaluation,
            state=branch_state,
            P=branch_P,
            b_grid=prepared.b_grid,
            shared=shared_2023,
            supply_rule=prepared.supply_rule,
            outside_flow=float(spec["outside_flow"]),
            retention=float(spec["retention"]),
            initial_mass_2007=prepared.initial_mass_2007,
            mass_2023=mass_2023,
            next_bridge_year=None,
            grid_fallback=fallback_2023,
        )
        paths[label] = [row_2023]
        branch_states[label] = next_state
        branch_parameters[label] = branch_P
    if progress_dir is not None:
        write_csv(
            progress_dir / "continuation_progress.csv",
            build_paired_rows(common_history, paths),
        )
        write_json(
            progress_dir / "latest_completed_period.json",
            {
                "status": "running",
                "phase": "branched_from_exact_matched_2023_state",
                "calendar_year": TERMINAL_YEAR,
                "completed_post_2023_dates_per_path": 0,
            },
        )

    for post_index in range(1, post_2023_periods + 1):
        period = TRANSITION_PERIODS + post_index
        for label, spec in closure_specs.items():
            branch_state = branch_states[label]
            branch_P = branch_parameters[label]
            branch_P.psi_child = new_psi
            evaluation, shared, fallback = evaluate_state(
                branch_state,
                branch_P,
                prepared.b_grid,
                prepared.supply_rule,
                counter,
                market_tol,
                market_max_iter,
            )
            row, next_state = advance_from_evaluation(
                label=label,
                period_from_2007=period,
                evaluation=evaluation,
                state=branch_state,
                P=branch_P,
                b_grid=prepared.b_grid,
                shared=shared,
                supply_rule=prepared.supply_rule,
                outside_flow=float(spec["outside_flow"]),
                retention=float(spec["retention"]),
                initial_mass_2007=prepared.initial_mass_2007,
                mass_2023=mass_2023,
                next_bridge_year=None,
                grid_fallback=fallback,
            )
            paths[label].append(row)
            branch_states[label] = next_state
            print(
                f"NO_POLICY_CONTINUATION {label} year={row['calendar_year']} "
                f"p={row['asset_price']:.6f} pop={row['population_index_2023']:.6f} "
                f"births={row['birth_children_topcode_adjusted']:.6f}",
                flush=True,
            )
        if progress_dir is not None:
            write_csv(
                progress_dir / "continuation_progress.csv",
                build_paired_rows(common_history, paths),
            )
            write_json(
                progress_dir / "latest_completed_period.json",
                {
                    "status": "running",
                    "phase": "paired_no_policy_post_2023_continuations",
                    "calendar_year": TERMINAL_YEAR + 4 * post_index,
                    "completed_post_2023_dates_per_path": post_index,
                    "requested_post_2023_dates_per_path": post_2023_periods,
                    "latest": {label: paths[label][-1] for label in sorted(paths)},
                },
            )

    for row in common_history:
        row["population_index_2023"] = float(row["adult_population"]) / mass_2023
    paired_history = {
        label: [dict(row) for row in common_history] + [rows[0]]
        for label, rows in paths.items()
    }
    state_gap = paired_history_state_gap(paired_history, HISTORY_STATE_COLUMNS)
    if state_gap["maximum_absolute_gap"] != 0.0:
        raise RuntimeError(
            "The two continuations do not begin from an identical 2023 state: "
            f"{state_gap}"
        )
    run_diagnostics = {
        "status": "passed",
        "shared_2007_2023_history_maximum_absolute_gap": state_gap[
            "maximum_absolute_gap"
        ],
        "shared_2007_2023_history_columns": list(HISTORY_STATE_COLUMNS),
        "common_history_periods": len(common_history) + 1,
        "post_2023_periods": post_2023_periods,
        "model_solve_count": counter.total,
        "mass_2023": mass_2023,
        "matched_2023_pre_fertility_distribution_sha256": matched_state_sha256,
        "matched_2023_birth_vintage_queue_sha256": canonical_json_sha256(
            state_2023.scheduled_entries
        ),
    }
    if progress_dir is not None:
        write_json(
            progress_dir / "latest_completed_period.json",
            {
                "status": "complete",
                "phase": "paired_no_policy_post_2023_continuations",
                "terminal_calendar_year": TERMINAL_YEAR + 4 * post_2023_periods,
                "diagnostics": run_diagnostics,
            },
        )
    return common_history, paths, run_diagnostics


def paired_history_state_gap(
    paths: dict[str, list[dict[str, Any]]], columns: tuple[str, ...]
) -> dict[str, Any]:
    labels = sorted(paths)
    if len(labels) != 2 or not paths[labels[0]] or not paths[labels[1]]:
        raise ValueError("Paired history audit requires two nonempty paths")
    left_rows = paths[labels[0]]
    right_rows = paths[labels[1]]
    if len(left_rows) != len(right_rows):
        raise ValueError("Paired histories have different lengths")
    detailed = []
    gaps = {key: 0.0 for key in columns}
    for row_index, (left, right) in enumerate(
        zip(left_rows, right_rows, strict=True)
    ):
        if int(left.get("calendar_year", row_index)) != int(
            right.get("calendar_year", row_index)
        ):
            raise ValueError("Paired histories have different calendar dates")
        for key in columns:
            gap = abs(float(left[key]) - float(right[key]))
            gaps[key] = max(gaps[key], gap)
            detailed.append(
                {
                    "row_index": row_index,
                    "calendar_year": left.get("calendar_year"),
                    "field": key,
                    "absolute_gap": gap,
                }
            )
    return {
        "labels": labels,
        "gaps": gaps,
        "maximum_absolute_gap": max(gaps.values(), default=0.0),
        "compared_rows": len(left_rows),
        "details": detailed,
    }


def validate_reconstructed_history(
    common_history: list[dict[str, Any]],
    open_path: list[dict[str, Any]],
    case_transition_path: Path,
) -> dict[str, Any]:
    reconstructed = list(common_history) + [open_path[0]]
    saved = read_csv(case_transition_path)
    if len(saved) != TRANSITION_PERIODS + 1:
        raise RuntimeError(
            "Selected calibration case must end at 2023 with exactly five dates: "
            f"found={len(saved)}"
        )
    if [int(row["period"]) for row in saved] != list(range(TRANSITION_PERIODS + 1)):
        raise RuntimeError("Selected case transition periods are not 0,...,4")
    comparisons: list[dict[str, Any]] = []
    maximum = 0.0
    for period, (actual, expected) in enumerate(zip(reconstructed, saved, strict=True)):
        for key in HISTORY_STATE_COLUMNS:
            gap = abs(float(actual[key]) - float(expected[key]))
            comparisons.append({"period": period, "field": key, "absolute_gap": gap})
            maximum = max(maximum, gap)
    if maximum > RECONSTRUCTION_TOLERANCE:
        worst = max(comparisons, key=lambda row: float(row["absolute_gap"]))
        raise RuntimeError(
            "Selected 2007--2023 case did not reproduce under the pinned code: "
            f"maximum_gap={maximum:.3e}, worst={worst}"
        )
    return {
        "status": "passed",
        "periods": TRANSITION_PERIODS + 1,
        "fields": list(HISTORY_STATE_COLUMNS),
        "maximum_absolute_gap": maximum,
        "tolerance": RECONSTRUCTION_TOLERANCE,
        "comparisons": comparisons,
    }


def closed_schedule_brackets(rows: list[dict[str, Any]]) -> list[tuple[float, float]]:
    ordered = sorted(rows, key=lambda row: float(row["asset_price"]))
    brackets: list[tuple[float, float]] = []
    for left, right in zip(ordered[:-1], ordered[1:]):
        a = float(left["renewal_residual_ratio"])
        b = float(right["renewal_residual_ratio"])
        if not math.isfinite(a) or not math.isfinite(b):
            continue
        if a == 0.0:
            brackets.append((float(left["asset_price"]), float(left["asset_price"])))
        elif a * b < 0.0:
            brackets.append((float(left["asset_price"]), float(right["asset_price"])))
    if ordered and float(ordered[-1]["renewal_residual_ratio"]) == 0.0:
        value = float(ordered[-1]["asset_price"])
        brackets.append((value, value))
    deduplicated: list[tuple[float, float]] = []
    for bracket in brackets:
        if bracket not in deduplicated:
            deduplicated.append(bracket)
    return deduplicated


def solve_closed_stationary_endpoint(
    *,
    chain: Any,
    base_overrides: dict[str, Any],
    old_price: float,
    new_psi_child: float,
    supply_rule: calendar.HousingSupplyRule,
    price_min_ratio: float,
    price_max_ratio: float,
    grid_points: int,
    outdir: Path,
    root_residual_tolerance: float = ROOT_RESIDUAL_TOLERANCE,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not 0.0 < price_min_ratio < price_max_ratio or grid_points < 5:
        raise ValueError("Closed-endpoint price grid is invalid")
    if not 0.0 < float(root_residual_tolerance) <= ROOT_RESIDUAL_TOLERANCE:
        raise ValueError(
            "The endpoint root tolerance must be positive and no looser than "
            f"{ROOT_RESIDUAL_TOLERANCE:g}."
        )
    cache: dict[float, dict[str, Any]] = {}

    def evaluate(asset_price: float) -> dict[str, Any]:
        key = round(float(asset_price), 14)
        if key not in cache:
            solution, parameters, price, elapsed = closure.solve_pe(
                chain,
                base_overrides,
                asset_price=float(asset_price),
                psi_child=float(new_psi_child),
            )
            readout = closure.readout(
                chain,
                solution,
                parameters,
                price,
                label="postshock_closed_endpoint_search",
                price_ratio=float(asset_price / old_price),
                psi_child=float(new_psi_child),
                elapsed=elapsed,
            )
            renewal = closure.topcode_consistent_renewal_accounting(
                solution, parameters
            )
            E = float(solution.entry_rate) / max(float(solution.total_mass), 1e-15)
            adjusted_births = float(
                renewal["topcode_adjusted_birth_children"]
            ) / max(float(solution.total_mass), 1e-15)
            B = (
                transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
                * adjusted_births
            )
            ratio = B / max(E, 1e-15)
            demand_per_adult = float(readout["housing_demand_per_adult"])
            housing_supply = float(supply_rule.quantity(price)[0])
            scale = housing_supply / max(demand_per_adult, 1e-15)
            row = dict(readout)
            row.update(
                policy_case="none",
                closure="closed_national_benchmark",
                outside_flow_M=0.0,
                retention_rho=1.0,
                replacement_fertility=transition.REPLACEMENT_FERTILITY,
                effective_birth_to_household_conversion=(
                    transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
                ),
                stationary_entry_flow_E=E,
                topcode_adjusted_birth_children_per_adult=adjusted_births,
                queue_mature_flow_B=B,
                queue_B_over_E=ratio,
                renewal_residual_ratio=ratio - 1.0,
                stationary_housing_supply=housing_supply,
                stationary_population_scale=scale,
                scaled_housing_demand=scale * demand_per_adult,
                housing_clearing_residual=(
                    scale * demand_per_adult - housing_supply
                ),
            )
            cache[key] = row
            progress_rows = sorted(
                cache.values(), key=lambda item: float(item["asset_price"])
            )
            write_csv(
                outdir / "closed_stationary_schedule_progress.csv", progress_rows
            )
            write_json(
                outdir / "latest_endpoint_search.json",
                {
                    "status": "running",
                    "closure": "closed_national_benchmark",
                    "evaluated_prices": len(progress_rows),
                    "latest": row,
                },
            )
        return cache[key]

    prices = old_price * np.geomspace(
        float(price_min_ratio), float(price_max_ratio), int(grid_points)
    )
    schedule = [evaluate(float(price)) for price in prices]
    brackets = closed_schedule_brackets(schedule)
    if len(brackets) == 0:
        endpoint = {
            "status": "no_positive_renewal_root_on_audited_grid",
            "usable_closed_root": False,
            "between_steady_states_label_allowed": False,
            "policy_case": "none",
            "closure": "closed_national_benchmark",
            "price_min_ratio": price_min_ratio,
            "price_max_ratio": price_max_ratio,
            "audit_grid_points": grid_points,
            "maximum_queue_B_over_E": max(
                float(row["queue_B_over_E"]) for row in schedule
            ),
            "minimum_queue_B_over_E": min(
                float(row["queue_B_over_E"]) for row in schedule
            ),
            "interpretation": (
                "The finite-horizon national closed path is reported, but the "
                "post-shock stationary schedule has no sign-changing positive-price "
                "renewal root on the declared grid. It is not labeled a transition "
                "between steady states."
            ),
        }
        write_csv(outdir / "closed_stationary_schedule.csv", schedule)
        write_json(outdir / "closed_stationary_endpoint.json", endpoint)
        write_json(
            outdir / "latest_endpoint_search.json",
            {"status": "complete", "closed_stationary_endpoint": endpoint},
        )
        return endpoint, schedule
    if len(brackets) > 1:
        endpoint = {
            "status": "multiple_positive_renewal_roots_on_audited_grid",
            "usable_closed_root": False,
            "between_steady_states_label_allowed": False,
            "policy_case": "none",
            "closure": "closed_national_benchmark",
            "brackets": brackets,
            "interpretation": (
                "More than one sign change appears on the audited stationary "
                "schedule. No unique between-steady-states interpretation is assigned."
            ),
        }
        write_csv(outdir / "closed_stationary_schedule.csv", schedule)
        write_json(outdir / "closed_stationary_endpoint.json", endpoint)
        write_json(
            outdir / "latest_endpoint_search.json",
            {"status": "complete", "closed_stationary_endpoint": endpoint},
        )
        return endpoint, schedule

    lower, upper = brackets[0]
    if lower == upper:
        root = evaluate(lower)
    else:
        low_row = evaluate(lower)
        low_gap = float(low_row["renewal_residual_ratio"])
        root = low_row
        for _ in range(40):
            midpoint = 0.5 * (lower + upper)
            root = evaluate(midpoint)
            gap = float(root["renewal_residual_ratio"])
            if abs(gap) <= float(root_residual_tolerance):
                break
            if low_gap * gap <= 0.0:
                upper = midpoint
            else:
                lower = midpoint
                low_gap = gap
    scale = float(root["stationary_population_scale"])
    renewal_gap = abs(float(root["renewal_residual_ratio"]))
    housing_gap = abs(float(root["housing_clearing_residual"]))
    usable = (
        float(root["asset_price"]) > 0.0
        and math.isfinite(scale)
        and scale > 0.0
        and renewal_gap <= float(root_residual_tolerance)
        and housing_gap <= 1e-10
    )
    endpoint = dict(root)
    endpoint.update(
        status="complete_usable_positive_root" if usable else "root_failed_acceptance_gates",
        usable_closed_root=usable,
        between_steady_states_label_allowed=usable,
        renewal_root_absolute_residual=renewal_gap,
        renewal_root_declared_tolerance=float(root_residual_tolerance),
        housing_clearing_absolute_residual=housing_gap,
        audited_root_bracket=brackets[0],
        interpretation=(
            "A positive closed stationary endpoint is verified; the closed path "
            "may be described as movement toward a new steady state."
            if usable
            else "A sign change was found, but the numerical root failed a declared "
            "acceptance gate and is not used as a steady-state endpoint."
        ),
    )
    schedule = sorted(cache.values(), key=lambda row: float(row["asset_price"]))
    write_csv(outdir / "closed_stationary_schedule.csv", schedule)
    write_json(outdir / "closed_stationary_endpoint.json", endpoint)
    write_json(
        outdir / "latest_endpoint_search.json",
        {"status": "complete", "closed_stationary_endpoint": endpoint},
    )
    return endpoint, schedule


def verify_open_endpoint(
    endpoint: dict[str, Any], outside_flow: float, retention: float
) -> dict[str, Any]:
    verified = dict(endpoint)
    if endpoint.get("status") != "complete":
        verified.update(
            usable_open_endpoint=False,
            renewal_identity_residual=None,
            interpretation=(
                "The open stationary search did not return a usable endpoint; "
                "no stationary level is attached to the open continuation."
            ),
        )
        return verified
    if not math.isclose(
        float(endpoint.get("total_mass", math.nan)),
        1.0,
        rel_tol=0.0,
        abs_tol=2e-10,
    ):
        raise RuntimeError(
            "Open endpoint stationary distribution is not unit normalized: "
            f"total_mass={endpoint.get('total_mass')}"
        )
    scale = float(endpoint["stationary_population_scale"])
    B = float(endpoint["queue_mature_entrant_flow_B"])
    denominator = float(endpoint["renewal_denominator"])
    E = denominator + float(retention) * B
    residual = scale * E - float(outside_flow) - float(retention) * scale * B
    usable = (
        scale > 0.0
        and math.isfinite(scale)
        and abs(residual) <= ACCOUNTING_TOLERANCE
        and abs(float(endpoint["fixed_stock_relative_market_gap"])) <= 2.5e-5
    )
    if not usable:
        raise RuntimeError(
            "Open stationary endpoint failed renewal or market gates: "
            f"renewal={residual:.3e}, market={endpoint['fixed_stock_relative_market_gap']}"
        )
    verified.update(
        usable_open_endpoint=True,
        renewal_scale_formula="S = M / (E - rho * B)",
        renewal_identity_residual=residual,
        outside_flow_M=float(outside_flow),
        retention_rho=float(retention),
        interpretation=(
            "Open/CBSA sensitivity endpoint under an externally fixed outside-origin "
            "share; the share and resulting level are not identified by the calibration."
        ),
    )
    return verified


def build_paired_rows(
    common_history: list[dict[str, Any]],
    paths: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for label in sorted(paths):
        for row in common_history:
            copied = dict(row)
            copied["scenario"] = label
            copied["phase"] = "shared_calibrated_history"
            rows.append(copied)
        for index, row in enumerate(paths[label]):
            copied = dict(row)
            copied["phase"] = (
                "shared_matched_2023_state" if index == 0 else "post_2023_continuation"
            )
            rows.append(copied)
    return sorted(rows, key=lambda row: (str(row["scenario"]), int(row["calendar_year"])))


def make_figures(
    *,
    common_history: list[dict[str, Any]],
    paths: dict[str, list[dict[str, Any]]],
    closed_endpoint: dict[str, Any],
    open_endpoint: dict[str, Any],
    closed_schedule: list[dict[str, Any]],
    outdir: Path,
) -> None:
    colors = {
        "closed_national_benchmark": "#1f4e79",
        "open_cbsa_sensitivity": "#b45f06",
    }
    labels = {
        "closed_national_benchmark": "National closed benchmark",
        "open_cbsa_sensitivity": "Open/CBSA sensitivity",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)
    history_years = [int(row["calendar_year"]) for row in common_history]
    history_population = [float(row["population_index_2023"]) for row in common_history]
    history_price = [float(row["asset_price_index"]) for row in common_history]
    for axis, history, key, title in (
        (axes[0], history_population, "population_index_2023", "Adult-household population"),
        (axes[1], history_price, "asset_price_index", "Housing-cost index"),
    ):
        axis.plot(history_years, history, color="0.45", linewidth=2.0)
        for scenario, rows in paths.items():
            axis.plot(
                [int(row["calendar_year"]) for row in rows],
                [float(row[key]) for row in rows],
                color=colors[scenario],
                linewidth=2.0,
                label=labels[scenario],
            )
        axis.axvline(TERMINAL_YEAR, color="0.2", linestyle=":", linewidth=1.0)
        axis.axhline(1.0, color="0.75", linewidth=0.8)
        axis.set_title(title)
        axis.set_xlabel("Calendar year")
        axis.grid(alpha=0.2)
    if bool(closed_endpoint.get("usable_closed_root", False)):
        axes[0].axhline(
            float(closed_endpoint["stationary_population_scale"])
            / max(float(paths["closed_national_benchmark"][0]["adult_population"]), 1e-15),
            color=colors["closed_national_benchmark"],
            linestyle="--",
            linewidth=1.0,
        )
        axes[1].axhline(
            float(closed_endpoint["price_ratio"]),
            color=colors["closed_national_benchmark"],
            linestyle="--",
            linewidth=1.0,
        )
    if bool(open_endpoint.get("usable_open_endpoint", False)):
        axes[0].axhline(
            float(open_endpoint["stationary_population_scale"])
            / max(float(paths["open_cbsa_sensitivity"][0]["adult_population"]), 1e-15),
            color=colors["open_cbsa_sensitivity"],
            linestyle="--",
            linewidth=1.0,
        )
        axes[1].axhline(
            float(open_endpoint["price_ratio"]),
            color=colors["open_cbsa_sensitivity"],
            linestyle="--",
            linewidth=1.0,
        )
    axes[0].legend(frameon=False, fontsize=8)
    fig.savefig(outdir / "paired_continuation_levels.png", dpi=220)
    fig.savefig(outdir / "paired_continuation_levels.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), constrained_layout=True)
    for scenario, rows in paths.items():
        years = [int(row["calendar_year"]) for row in rows]
        axes[0].plot(
            years,
            [float(row["topcode_adjusted_births_per_adult"]) for row in rows],
            color=colors[scenario],
            linewidth=2.0,
            label=labels[scenario],
        )
        axes[1].plot(
            years,
            [float(row["queue_B_over_current_E"]) for row in rows],
            color=colors[scenario],
            linewidth=2.0,
        )
    axes[0].set_title("Top-code-adjusted births per adult household")
    axes[1].set_title("Maturing birth-vintage flow / current entry")
    for axis in axes:
        axis.axvline(TERMINAL_YEAR, color="0.2", linestyle=":", linewidth=1.0)
        axis.set_xlabel("Calendar year")
        axis.grid(alpha=0.2)
    axes[1].axhline(1.0, color="0.45", linestyle="--", linewidth=1.0)
    axes[0].legend(frameon=False, fontsize=8)
    fig.savefig(outdir / "paired_renewal_diagnostics.png", dpi=220)
    fig.savefig(outdir / "paired_renewal_diagnostics.pdf")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(6.2, 4.2), constrained_layout=True)
    axis.plot(
        [float(row["price_ratio"]) for row in closed_schedule],
        [float(row["queue_B_over_E"]) for row in closed_schedule],
        color=colors["closed_national_benchmark"],
        linewidth=2.0,
    )
    axis.axhline(1.0, color="0.2", linestyle="--", linewidth=1.0)
    axis.set_xscale("log")
    axis.set_xlabel("Housing-cost ratio to 2007")
    axis.set_ylabel(r"Top-code-adjusted births / (2.1 $\times$ entry)")
    axis.set_title("Closed stationary renewal schedule")
    axis.grid(alpha=0.2)
    fig.savefig(outdir / "closed_stationary_renewal_schedule.png", dpi=220)
    fig.savefig(outdir / "closed_stationary_renewal_schedule.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if int(args.post_2023_periods) < 1:
        raise ValueError("--post-2023-periods must be positive")
    if not 0.0 < float(args.market_tol) <= 2e-4:
        raise ValueError("--market-tol must lie in (0, 2e-4]")
    started = time.perf_counter()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    chain, model = transition.configure_sequential_model()
    contracts = validate_input_contracts(
        report_path=args.selected_report,
        task_summary_path=args.selected_task_summary,
        case_dir=args.selected_case_dir,
        case_transition_path=args.selected_case_transition,
        source_path=args.source,
        expected_report_sha256=args.expected_report_sha256,
        expected_task_summary_sha256=args.expected_task_summary_sha256,
        expected_case_transition_sha256=args.expected_case_transition_sha256,
        expected_source_sha256=args.expected_source_sha256,
        expected_target_fingerprint=args.expected_target_fingerprint,
        expected_code_bundle_sha256=args.expected_code_bundle_sha256,
        expected_renewal_contract_sha256=args.expected_renewal_contract_sha256,
        expected_scientific_contract_sha256=(
            args.expected_scientific_contract_sha256
        ),
        expected_selection_sha256=args.expected_selection_sha256,
        chain=chain,
        model=model,
    )
    prepared = prepare_model(contracts, args.source.resolve(), chain, model)
    outside_audit = outside_share_invariance_audit(
        prepared.E_old,
        prepared.B_old,
        contracts["outside_origin_entry_share"],
    )
    common_history, paths, transition_gates = reconstruct_and_branch(
        prepared,
        contracts,
        post_2023_periods=int(args.post_2023_periods),
        market_tol=float(args.market_tol),
        market_max_iter=int(args.market_max_iter),
        progress_dir=outdir,
    )
    history_reproduction = validate_reconstructed_history(
        common_history,
        paths["open_cbsa_sensitivity"],
        Path(contracts["paths"]["selected_case_transition"]),
    )
    write_csv(outdir / "shared_2007_2023_history.csv", common_history + [paths["open_cbsa_sensitivity"][0]])
    write_csv(
        outdir / "paired_continuation_path.csv",
        build_paired_rows(common_history, paths),
    )
    write_csv(outdir / "history_reproduction_audit.csv", history_reproduction["comparisons"])
    write_csv(outdir / "outside_share_invariance_audit.csv", outside_audit["rows"])

    best = contracts["report_best"]
    closed_endpoint, closed_schedule = solve_closed_stationary_endpoint(
        chain=prepared.chain,
        base_overrides=prepared.base_overrides,
        old_price=prepared.old_price,
        new_psi_child=float(best["new_psi_child"]),
        supply_rule=prepared.supply_rule,
        price_min_ratio=float(args.closed_price_min_ratio),
        price_max_ratio=float(args.closed_price_max_ratio),
        grid_points=int(args.closed_grid_points),
        outdir=outdir,
    )
    report_old = contracts["renewal_old_state"]
    open_dir = outdir / "open_endpoint"
    open_dir.mkdir(parents=True, exist_ok=True)
    open_raw = transition.solve_stationary_open_endpoint(
        chain=prepared.chain,
        base_overrides=prepared.base_overrides,
        old_price=prepared.old_price,
        new_psi_child=float(best["new_psi_child"]),
        outside_flow=float(report_old["outside_flow_M"]),
        retention=float(report_old["retention_rho"]),
        effective_birth_to_household_conversion=(
            transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
        ),
        supply_rule=prepared.supply_rule,
        policy_case="none",
        outdir=open_dir,
    )
    open_endpoint = verify_open_endpoint(
        open_raw,
        float(report_old["outside_flow_M"]),
        float(report_old["retention_rho"]),
    )
    write_json(outdir / "open_stationary_endpoint.json", open_endpoint)
    make_figures(
        common_history=common_history,
        paths=paths,
        closed_endpoint=closed_endpoint,
        open_endpoint=open_endpoint,
        closed_schedule=closed_schedule,
        outdir=outdir,
    )

    def path_gates(rows: list[dict[str, Any]]) -> dict[str, Any]:
        return {
            "maximum_market_residual": max(
                abs(float(row["relative_market_residual"])) for row in rows
            ),
            "maximum_mass_residual": max(
                abs(float(row["mass_accounting_residual"])) for row in rows
            ),
            "maximum_feasibility_projection_mass": max(
                abs(float(row["feasibility_frontier_projection_mass"])) for row in rows
            ),
            "minimum_distribution_mass": min(
                float(row["min_distribution_mass"]) for row in rows
            ),
            "maximum_nonfinite_distribution_count": max(
                int(row["nonfinite_distribution_count"]) for row in rows
            ),
        }

    summary = {
        "status": "complete_no_policy_post2023_continuation_pair",
        "paper_scope": "new equilibrium and transition validation; no policy experiment",
        "policy_case": "none",
        "fiscal_change": "none",
        "selected_inputs": contracts["paths"],
        "provenance": contracts["hashes"],
        "model_profile": _selected(contracts["report"], "model_profile"),
        "target_set": _selected(contracts["report"], "target_set"),
        "target_fingerprint": contracts["hashes"]["target_fingerprint"],
        "renewal_accounting_contract": contracts["renewal_contract"],
        "replacement_interpretation": (
            "The externally fixed 2.1 replacement statistic bundles sex composition, "
            "survival to household entry, and household formation. The model does not "
            "add an explicit child-death state."
        ),
        "outside_origin_entry_share": contracts["outside_origin_entry_share"],
        "outside_origin_entry_status": (
            "external open/CBSA sensitivity; not identified by the 2007--2023 "
            "transition calibration"
        ),
        "outside_share_invariance_audit": outside_audit,
        "history_reproduction_audit": {
            key: value
            for key, value in history_reproduction.items()
            if key != "comparisons"
        },
        "paired_initial_state_audit": transition_gates,
        "supply_normalization": prepared.supply_normalization,
        "stationary_operator_gates": prepared.operator_gates,
        "old_stationary_gates": prepared.stationary_gates,
        "post_2023_periods": int(args.post_2023_periods),
        "terminal_calendar_year": TERMINAL_YEAR + 4 * int(args.post_2023_periods),
        "paths": {
            label: {
                "closure": (
                    "M=0, rho=1"
                    if label == "closed_national_benchmark"
                    else "M=s_out*E_old, rho=(1-s_out)*E_old/B_old"
                ),
                "path_gates": path_gates(rows),
                "2023": rows[0],
                "last_date": rows[-1],
            }
            for label, rows in paths.items()
        },
        "closed_stationary_endpoint": closed_endpoint,
        "open_stationary_endpoint": open_endpoint,
        "between_steady_states_statement": (
            "The national closed path has a verified positive stationary endpoint "
            "and can be described as movement between steady states."
            if bool(closed_endpoint.get("between_steady_states_label_allowed", False))
            else "The national closed finite-horizon path is a benchmark only. No "
            "usable positive closed stationary root was verified, so it is not "
            "described as a transition between steady states."
        ),
        "timing_interpretation": (
            "Both paths use the same calibrated 2007--2023 history and exact matched "
            "2023 pre-fertility distribution. Population-law differences affect only "
            "entrant cohorts after 2023."
        ),
        "asset_price_interpretation": (
            "A sequence of temporary housing-market equilibria under the maintained "
            "stationary user-cost identity; not perfect-foresight asset pricing or welfare."
        ),
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json(outdir / "summary.json", summary)
    readme = f"""# Paired no-policy post-2023 continuations

This packet reconstructs the selected transition calibration through 2023 and
branches from the exact same matched 2023 state.  The national benchmark sets
`M=0` and `rho=1`.  The open/CBSA sensitivity uses the selected report's
outside-origin share `{contracts['outside_origin_entry_share']:.6g}`.  That
share is external and is **not identified** by the 2007--2023 calibration.

There is no housing policy, tax change, transfer, grant, or fiscal experiment.
Renewal uses top-code-adjusted births divided by 2.1 and a four-slot queue,
which makes a birth affect household entry five model dates (20 years) later.

Closed endpoint status: `{closed_endpoint['status']}`.  {summary['between_steady_states_statement']}

The path is a sequence of temporary equilibria under the maintained user-cost
identity, not a perfect-foresight asset-price or welfare calculation.

Core artifacts:

- `summary.json`: contracts, gates, endpoint classifications, and terminal rows;
- `paired_continuation_path.csv`: the common 2007--2023 history and both continuations;
- `closed_stationary_schedule.csv`: the directly solved closed renewal schedule;
- `closed_stationary_endpoint.json` and `open_stationary_endpoint.json`;
- three stable diagnostic figures in PDF and PNG.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")
    artifact_names = [
        "README.md",
        "summary.json",
        "paired_continuation_path.csv",
        "shared_2007_2023_history.csv",
        "history_reproduction_audit.csv",
        "outside_share_invariance_audit.csv",
        "continuation_progress.csv",
        "latest_completed_period.json",
        "closed_stationary_schedule.csv",
        "closed_stationary_schedule_progress.csv",
        "closed_stationary_endpoint.json",
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
    ]
    missing_artifacts = [name for name in artifact_names if not (outdir / name).is_file()]
    if missing_artifacts:
        raise RuntimeError(f"Continuation packet is incomplete: {missing_artifacts}")
    manifest = {
        "status": "complete_no_policy_post2023_continuation_manifest",
        "driver": str(Path(__file__).resolve()),
        "driver_sha256": file_sha256(Path(__file__).resolve()),
        "input_provenance": contracts["hashes"],
        "artifacts": {
            name: file_sha256(outdir / name) for name in artifact_names
        },
    }
    write_json(outdir / "manifest.json", manifest)
    print(
        "NO_POLICY_CONTINUATIONS_COMPLETE "
        f"closed_endpoint={closed_endpoint['status']} "
        f"open_endpoint={open_endpoint['status']} output={outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
