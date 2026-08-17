#!/usr/bin/env python3
"""Fixed-parameter fertility-slope diagnostic around the final E5F transition.

This is orchestration, not a calibration routine.  For exactly one predeclared
common multiplier on the first-birth and continuation-birth logit scales, it

1. calls the existing dated-transition driver, which independently derives the
   2007 preference intercept so completed fertility is 2.1 and then evaluates
   all twelve target moments in 2023; and
2. calls the existing closed-stationary endpoint helper on the resulting 2023
   preference regime over the exact 25-point price grid.

Every other parameter and the final ridge estimate's 2007--2023 preference
change are held fixed.  The output is therefore a slope/fit frontier diagnostic
and must never be classified as an estimated or promotion-eligible calibration.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))


SCHEMA = "e5f_final_closed_root_slope_frontier_factor_v1"
COLLECTION_SCHEMA = "e5f_final_closed_root_slope_frontier_report_v1"
FACTORS = (1.0, 0.5, 0.25, 0.1, 0.05)
PRICE_MIN_RATIO = 1.0e-4
PRICE_MAX_RATIO = 3.0
PRICE_GRID_POINTS = 25
REPLACEMENT_FERTILITY = 2.1
MODEL_PROFILE = "e5f-income-entry"
TARGET_COUNT = 12
SELECTED_COMPARISON_TOLERANCE = 5.0e-12


class ContractError(RuntimeError):
    """Raised when an input or generated artifact violates the fixed contract."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ridge-report", type=Path, required=True)
    parser.add_argument("--selected-task-summary", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--factor-index", type=int)
    parser.add_argument("--collect", action="store_true")
    parser.add_argument(
        "--run-stage", choices=("smoke", "factor", "collect"), required=True
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        help="Root containing factor_001,...,factor_005; required with --collect.",
    )
    parser.add_argument("--expected-ridge-report-sha256", required=True)
    parser.add_argument("--expected-task-summary-sha256", required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-target-set", required=True)
    parser.add_argument("--expected-target-fingerprint", required=True)
    parser.add_argument("--expected-code-bundle-sha256", required=True)
    parser.add_argument("--expected-scientific-contract-sha256", required=True)
    parser.add_argument("--expected-selection-sha256", required=True)
    parser.add_argument("--expected-selected-candidate-sha256", required=True)
    parser.add_argument("--expected-selected-candidate-id", type=int, required=True)
    parser.add_argument("--expected-renewal-contract-sha256", required=True)
    parser.add_argument("--expected-dated-contract-sha256", required=True)
    parser.add_argument("--expected-frontier-driver-sha256", required=True)
    parser.add_argument("--expected-continuation-driver-sha256", required=True)
    parser.add_argument("--expected-launcher-sha256", required=True)
    parser.add_argument("--expected-diagnostic-bundle-sha256", required=True)
    parser.add_argument("--smoke-receipt", type=Path)
    parser.add_argument("--expected-smoke-receipt-sha256")
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (np.integer, np.floating)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def atomic_write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def atomic_write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ContractError(f"Refusing to write an empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(jsonable(rows))
    temporary.replace(path)


def read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContractError(f"Cannot read JSON object {path}: {error}") from error
    if not isinstance(payload, dict):
        raise ContractError(f"JSON artifact must contain an object: {path}")
    return payload


def read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except OSError as error:
        raise ContractError(f"Cannot read CSV artifact {path}: {error}") from error


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        jsonable(payload),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def require_file_hash(path: Path, expected: str, label: str) -> str:
    if not path.is_file():
        raise ContractError(f"{label} is unavailable: {path}")
    actual = file_sha256(path)
    if actual != str(expected):
        raise ContractError(
            f"{label} hash mismatch: actual={actual}, expected={expected}, path={path}"
        )
    return actual


def diagnostic_code_contract(
    frontier_driver: Path, continuation_driver: Path, launcher: Path
) -> dict[str, str]:
    """Hash every executable that defines this diagnostic's numerical loop."""
    return {
        "frontier_driver_sha256": file_sha256(frontier_driver),
        "continuation_driver_sha256": file_sha256(continuation_driver),
        "launcher_sha256": file_sha256(launcher),
    }


def validate_smoke_receipt(
    receipt_path: Path,
    expected_receipt_sha256: str,
    expected_hashes: dict[str, str],
) -> dict[str, str]:
    """Require a complete factor-0.5 smoke under the identical input contract."""
    receipt_path = receipt_path.resolve()
    receipt_sha = require_file_hash(
        receipt_path, expected_receipt_sha256, "Exact-loop smoke receipt"
    )
    receipt = read_json(receipt_path)
    require_equal(receipt.get("schema"), SCHEMA, "smoke receipt schema")
    require_equal(
        receipt.get("status"),
        "complete_fixed_parameter_slope_diagnostic",
        "smoke receipt status",
    )
    require_equal(receipt.get("run_stage"), "smoke", "smoke receipt run stage")
    if bool(receipt.get("promotion_eligible", True)):
        raise ContractError("Smoke receipt is mislabeled promotion eligible")
    require_equal(receipt.get("factor_index"), 2, "smoke receipt factor index")
    require_close(receipt.get("kappa_factor"), 0.5, "smoke receipt factor", 1.0e-14)
    if int(receipt.get("target_count", -1)) != TARGET_COUNT:
        raise ContractError("Smoke receipt does not retain all twelve targets")
    require_close(
        receipt.get("old_completed_fertility"),
        REPLACEMENT_FERTILITY,
        "smoke receipt old fertility",
        5.0e-4,
    )
    require_close(
        receipt.get("old_queue_B_over_E"),
        1.0,
        "smoke receipt old B/E",
        5.0e-4,
    )
    price_contract = dict(receipt.get("price_grid_contract") or {})
    require_close(
        price_contract.get("minimum_ratio"),
        PRICE_MIN_RATIO,
        "smoke receipt minimum price ratio",
        1.0e-15,
    )
    require_close(
        price_contract.get("maximum_ratio"),
        PRICE_MAX_RATIO,
        "smoke receipt maximum price ratio",
        1.0e-14,
    )
    require_equal(
        price_contract.get("declared_grid_points"),
        PRICE_GRID_POINTS,
        "smoke receipt grid points",
    )
    if not math.isfinite(float(receipt.get("transition_loss", math.nan))):
        raise ContractError("Smoke receipt transition loss is nonfinite")
    if not math.isfinite(float(receipt.get("maximum_B_over_E", math.nan))):
        raise ContractError("Smoke receipt maximum B/E is nonfinite")
    provenance = dict(receipt.get("input_provenance") or {})
    for name, expected in expected_hashes.items():
        require_equal(provenance.get(name), expected, f"smoke provenance {name}")

    manifest_path = receipt_path.with_name("manifest.json")
    manifest = read_json(manifest_path)
    require_equal(
        manifest.get("schema"),
        "e5f_final_closed_root_slope_frontier_factor_manifest_v1",
        "smoke manifest schema",
    )
    require_equal(
        manifest.get("status"),
        "complete_fixed_parameter_slope_diagnostic_manifest",
        "smoke manifest status",
    )
    require_equal(manifest.get("run_stage"), "smoke", "smoke manifest run stage")
    require_equal(manifest.get("factor_index"), 2, "smoke manifest factor index")
    require_close(manifest.get("kappa_factor"), 0.5, "smoke manifest factor", 1.0e-14)
    diagnostic_contract = dict(manifest.get("diagnostic_code_contract") or {})
    for name in (
        "frontier_driver_sha256",
        "continuation_driver_sha256",
        "launcher_sha256",
    ):
        require_equal(
            diagnostic_contract.get(name),
            expected_hashes[name],
            f"smoke manifest {name}",
        )
    require_equal(
        manifest.get("diagnostic_bundle_sha256"),
        expected_hashes["diagnostic_bundle_sha256"],
        "smoke manifest diagnostic bundle",
    )
    artifacts = dict(manifest.get("artifacts") or {})
    require_equal(
        artifacts.get(receipt_path.name), receipt_sha, "smoke manifest receipt hash"
    )
    for name, expected_hash in artifacts.items():
        require_equal(
            file_sha256(receipt_path.parent / name),
            expected_hash,
            f"smoke artifact hash {name}",
        )
    return {
        "smoke_receipt_sha256": receipt_sha,
        "smoke_manifest_sha256": file_sha256(manifest_path),
    }


def require_equal(actual: Any, expected: Any, label: str) -> None:
    if str(actual) != str(expected):
        raise ContractError(f"{label} mismatch: {actual!r} versus {expected!r}")


def require_close(
    actual: Any, expected: Any, label: str, tolerance: float = 1.0e-13
) -> None:
    if not math.isclose(
        float(actual), float(expected), rel_tol=0.0, abs_tol=float(tolerance)
    ):
        raise ContractError(f"{label} mismatch: {actual!r} versus {expected!r}")


def factor_for_index(index: int) -> float:
    if not 1 <= int(index) <= len(FACTORS):
        raise ContractError(f"factor index must lie in [1,{len(FACTORS)}]: {index}")
    return float(FACTORS[int(index) - 1])


def compare_theta(
    left: dict[str, Any], right: dict[str, Any], label: str
) -> None:
    if set(left) != set(right):
        raise ContractError(
            f"{label} coordinate mismatch: left={sorted(left)}, right={sorted(right)}"
        )
    for name in sorted(left):
        require_close(
            left[name],
            right[name],
            f"{label}.{name}",
            SELECTED_COMPARISON_TOLERANCE,
        )


def required_renewal_contract() -> dict[str, Any]:
    return {
        "replacement_fertility": 2.1,
        "effective_birth_to_household_conversion": 1.0 / 2.1,
        "birth_measure": "topcode_adjusted_birth_children",
        "top_bin_mean_children": 3.602359422009,
        "birth_vintage_queue_waiting_slots": 4,
        "birth_to_entry_effect_lag_dates": 5,
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


def validate_exact_renewal_contract(contract: dict[str, Any], label: str) -> None:
    required = required_renewal_contract()
    if set(contract) != set(required):
        raise ContractError(
            f"{label} fields mismatch: actual={sorted(contract)}, "
            f"required={sorted(required)}"
        )
    for name, expected in required.items():
        actual = contract[name]
        if isinstance(expected, float):
            require_close(actual, expected, f"{label}.{name}", 1.0e-14)
        else:
            require_equal(actual, expected, f"{label}.{name}")


def build_factor_center(
    selected_candidate: dict[str, Any], factor: float
) -> tuple[dict[str, Any], dict[str, float], float]:
    selected_theta = {
        str(name): float(value)
        for name, value in dict(selected_candidate["theta"]).items()
    }
    scaled_theta = dict(selected_theta)
    for name in ("kappa_fert", "kappa_fert_continuation"):
        scaled_theta[name] = float(selected_theta[name]) * float(factor)
        if scaled_theta[name] < 0.02:
            raise ContractError(
                f"factor {factor:g} sends {name} below the active lower bound: "
                f"{scaled_theta[name]:.12g}"
            )
    old_psi = float(selected_candidate["old_psi_child"])
    new_psi = float(selected_candidate["new_psi_child"])
    preference_change = new_psi - old_psi
    center = {
        "schema": "e5f_final_closed_root_slope_frontier_center_v1",
        "classification": "fixed_parameter_diagnostic_not_calibration",
        "kappa_factor": float(factor),
        "best_candidate": {
            "theta": scaled_theta,
            # The production driver reads the difference between these values.
            # It then derives a fresh old intercept and holds this difference fixed.
            "old_psi_child": old_psi,
            "new_psi_child": new_psi,
        },
    }
    return center, scaled_theta, preference_change


def validate_inputs(args: argparse.Namespace, *, require_live_code: bool) -> dict[str, Any]:
    report_path = args.ridge_report.resolve()
    task_path = args.selected_task_summary.resolve()
    source_path = args.source.resolve()
    frontier_driver_path = Path(__file__).resolve()
    continuation_driver_path = (
        TOOLS_ROOT / "run_e5f_post2023_no_policy_continuations.py"
    ).resolve()
    launcher_path = (
        ROOT / "code/cluster/submit_e5f_final_closed_root_slope_frontier.sh"
    ).resolve()
    diagnostic_hashes = diagnostic_code_contract(
        frontier_driver_path, continuation_driver_path, launcher_path
    )
    diagnostic_expectations = {
        "frontier_driver_sha256": args.expected_frontier_driver_sha256,
        "continuation_driver_sha256": args.expected_continuation_driver_sha256,
        "launcher_sha256": args.expected_launcher_sha256,
    }
    for name, expected in diagnostic_expectations.items():
        require_equal(diagnostic_hashes[name], expected, f"live diagnostic {name}")
    diagnostic_bundle_sha = canonical_json_sha256(diagnostic_hashes)
    require_equal(
        diagnostic_bundle_sha,
        args.expected_diagnostic_bundle_sha256,
        "live diagnostic bundle",
    )
    hashes = {
        "ridge_report_sha256": require_file_hash(
            report_path, args.expected_ridge_report_sha256, "Final ridge report"
        ),
        "selected_task_summary_sha256": require_file_hash(
            task_path, args.expected_task_summary_sha256, "Selected task summary"
        ),
        "source_sha256": require_file_hash(
            source_path, args.expected_source_sha256, "Stationary source"
        ),
        **diagnostic_hashes,
        "diagnostic_bundle_sha256": diagnostic_bundle_sha,
    }
    report = read_json(report_path)
    task = read_json(task_path)
    require_equal(
        report.get("schema"),
        "e5f_transition_ridge_refinement_report_v1",
        "ridge report schema",
    )
    require_equal(
        report.get("status"),
        "complete_refinement_with_two_independent_identity_repeats",
        "ridge report status",
    )
    if not bool(report.get("promotion_eligible", False)):
        raise ContractError("Final ridge report is not promotion eligible")
    repeat = dict(report.get("repeat_gate") or {})
    if not (
        int(repeat.get("required", -1)) == 2
        and int(repeat.get("completed", -1)) == 2
        and bool(repeat.get("identity_gate_passed", False))
        and bool(repeat.get("numerical_identity_gate_passed", False))
    ):
        raise ContractError(f"Final ridge repeat gate failed: {repeat}")
    require_equal(
        report.get("selection_sha256"),
        args.expected_selection_sha256,
        "selection hash",
    )
    require_equal(
        report.get("selected_candidate_sha256"),
        args.expected_selected_candidate_sha256,
        "selected-candidate hash",
    )
    require_equal(
        report.get("selected_candidate_id"),
        args.expected_selected_candidate_id,
        "selected-candidate id",
    )
    if task_path.parent.name != f"task_{int(args.expected_selected_candidate_id):03d}":
        raise ContractError(
            "Selected task path does not identify the expected ridge task: "
            f"{task_path.parent}"
        )
    require_equal(
        task.get("status"),
        "complete_transition_calibration_panel_task",
        "selected task status",
    )

    scientific = dict(report.get("scientific_contract") or {})
    scientific_sha = canonical_json_sha256(scientific)
    if not (
        scientific_sha
        == str(report.get("scientific_contract_sha256"))
        == str(args.expected_scientific_contract_sha256)
    ):
        raise ContractError(
            "Scientific-contract hash mismatch: "
            f"computed={scientific_sha}, report={report.get('scientific_contract_sha256')}, "
            f"expected={args.expected_scientific_contract_sha256}"
        )
    hashes["scientific_contract_sha256"] = scientific_sha

    best = dict(report.get("best_candidate") or {})
    task_best = dict(task.get("best_candidate") or {})
    if not best or not task_best:
        raise ContractError("Ridge report or selected task omits best_candidate")
    require_equal(
        best.get("candidate_id"),
        args.expected_selected_candidate_id,
        "best-candidate id",
    )
    require_equal(
        best.get("candidate_sha256"),
        args.expected_selected_candidate_sha256,
        "best-candidate hash",
    )
    panel = dict(task.get("panel_design") or {})
    require_equal(
        panel.get("center_sha256"),
        args.expected_selected_candidate_sha256,
        "selected task center hash",
    )
    compare_theta(dict(best["theta"]), dict(task_best["theta"]), "report/task theta")
    for name in (
        "new_psi_child",
        "transition_loss",
        "terminal_tfr",
        "terminal_childless_rate",
        "terminal_first_birth_housing_response",
    ):
        require_close(best[name], task_best[name], f"report/task {name}")
    require_close(
        best["old_psi_child"], task["old_psi_child"], "report/task old psi"
    )

    report_target_set = str(scientific.get("target_set"))
    report_target_fingerprint = str(scientific.get("target_fingerprint"))
    report_code_bundle = str(
        dict(scientific.get("code_fingerprints") or {}).get("bundle_sha256")
    )
    report_profile = str(dict(scientific.get("model_profile") or {}).get("name"))
    checks = (
        (report_target_set, args.expected_target_set, "report target set"),
        (
            report_target_fingerprint,
            args.expected_target_fingerprint,
            "report target fingerprint",
        ),
        (report_code_bundle, args.expected_code_bundle_sha256, "report code bundle"),
        (report_profile, MODEL_PROFILE, "report model profile"),
        (task.get("target_set"), args.expected_target_set, "task target set"),
        (
            task.get("target_fingerprint"),
            args.expected_target_fingerprint,
            "task target fingerprint",
        ),
        (
            dict(task.get("code_fingerprints") or {}).get("bundle_sha256"),
            args.expected_code_bundle_sha256,
            "task code bundle",
        ),
        (
            dict(task.get("model_profile") or {}).get("name"),
            MODEL_PROFILE,
            "task model profile",
        ),
        (scientific.get("source"), str(source_path), "report source path"),
        (scientific.get("source_sha256"), hashes["source_sha256"], "report source hash"),
        (task.get("source"), str(source_path), "task source path"),
        (task.get("source_sha256"), hashes["source_sha256"], "task source hash"),
    )
    for actual, expected, label in checks:
        require_equal(actual, expected, label)
    if int(task.get("target_count", -1)) != TARGET_COUNT:
        raise ContractError(f"Selected task target count is not {TARGET_COUNT}")
    if str(task.get("policy_case")) != "none" or int(task.get("post_2023_periods", -1)) != 0:
        raise ContractError("Selected task is not the no-policy 2007--2023 calibration")

    renewal_contracts = (
        (dict(scientific.get("renewal_accounting_contract") or {}), "report scientific renewal"),
        (dict(task.get("renewal_accounting_contract") or {}), "task renewal"),
        (dict(best.get("renewal_accounting_contract") or {}), "selected renewal"),
    )
    for contract, label in renewal_contracts:
        validate_exact_renewal_contract(contract, label)
        require_equal(
            canonical_json_sha256(contract),
            args.expected_renewal_contract_sha256,
            f"{label} hash",
        )
    require_equal(
        report.get("renewal_contract_sha256"),
        args.expected_renewal_contract_sha256,
        "report renewal hash",
    )
    require_close(
        scientific.get("old_completed_fertility_reference"),
        REPLACEMENT_FERTILITY,
        "report replacement fertility",
        1.0e-14,
    )
    require_close(
        task.get("old_completed_fertility_reference"),
        REPLACEMENT_FERTILITY,
        "task replacement fertility",
        1.0e-14,
    )
    require_close(
        best.get("old_completed_fertility"),
        REPLACEMENT_FERTILITY,
        "selected old completed fertility",
        5.0e-4,
    )

    dated_contracts = (
        dict(scientific.get("dated_measurement_contract") or {}),
        dict(task.get("dated_measurement_contract") or {}),
        dict(best.get("dated_measurement_contract") or {}),
    )
    for index, contract in enumerate(dated_contracts, start=1):
        require_equal(
            canonical_json_sha256(contract),
            args.expected_dated_contract_sha256,
            f"dated contract {index} hash",
        )
    require_equal(
        report.get("dated_contract_sha256"),
        args.expected_dated_contract_sha256,
        "report dated-contract hash",
    )

    chain = model = None
    if require_live_code:
        import run_e5f_open_population_transition as transition
        import run_e5f_transition_calibration as calibration
        from intergen_eqscale_seq_optimized.e5_profile import e5_target_system

        chain, model = transition.configure_sequential_model()
        live_code = calibration.code_fingerprint_contract(model)["bundle_sha256"]
        require_equal(live_code, args.expected_code_bundle_sha256, "live code bundle")
        target_system = e5_target_system()
        require_equal(target_system.name, args.expected_target_set, "live target set")
        require_equal(
            target_system.fingerprint,
            args.expected_target_fingerprint,
            "live target fingerprint",
        )
        if int(target_system.count) != TARGET_COUNT:
            raise ContractError(f"Live target count is not {TARGET_COUNT}")

    contract_hashes = {
        **hashes,
        "target_fingerprint": args.expected_target_fingerprint,
        "code_bundle_sha256": args.expected_code_bundle_sha256,
        "selection_sha256": args.expected_selection_sha256,
        "selected_candidate_sha256": args.expected_selected_candidate_sha256,
        "renewal_contract_sha256": args.expected_renewal_contract_sha256,
        "dated_contract_sha256": args.expected_dated_contract_sha256,
    }
    smoke_receipt_path: Path | None = None
    if args.run_stage == "smoke":
        if (
            args.smoke_receipt is not None
            or args.expected_smoke_receipt_sha256 is not None
        ):
            raise ContractError("Smoke stage forbids a recursive smoke receipt")
    else:
        if args.smoke_receipt is None or args.expected_smoke_receipt_sha256 is None:
            raise ContractError(
                f"{args.run_stage} stage requires a completed exact-loop smoke receipt"
            )
        smoke_receipt_path = args.smoke_receipt.resolve()
        contract_hashes.update(
            validate_smoke_receipt(
                smoke_receipt_path,
                args.expected_smoke_receipt_sha256,
                contract_hashes,
            )
        )

    return {
        "report_path": report_path,
        "task_path": task_path,
        "source_path": source_path,
        "report": report,
        "task": task,
        "scientific_contract": scientific,
        "selected_candidate": best,
        "smoke_receipt_path": smoke_receipt_path,
        "hashes": contract_hashes,
        "chain": chain,
        "model": model,
    }


class Heartbeat:
    def __init__(self, outdir: Path, factor_index: int | None, factor: float | None):
        self.outdir = outdir
        self.factor_index = factor_index
        self.factor = factor
        self.started = time.perf_counter()
        self.phase = "initializing"
        self._stop = threading.Event()
        self._lock = threading.RLock()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def start(self) -> None:
        self.write()
        self._thread.start()

    def set_phase(self, phase: str) -> None:
        with self._lock:
            self.phase = str(phase)
        self.write()

    def write(self) -> None:
        with self._lock:
            phase = self.phase
            status = (
                "complete"
                if phase == "complete"
                else "failed"
                if phase == "failed"
                else "running"
            )
            atomic_write_json(
                self.outdir / "heartbeat.json",
                {
                    "status": status,
                    "phase": phase,
                    "factor_index": self.factor_index,
                    "kappa_factor": self.factor,
                    "elapsed_seconds": time.perf_counter() - self.started,
                    "pid": os.getpid(),
                    "hostname": os.environ.get("HOSTNAME", "unknown"),
                    "updated_unix_seconds": time.time(),
                },
            )

    def _run(self) -> None:
        while not self._stop.wait(60.0):
            self.write()
            print(
                "FRONTIER_HEARTBEAT "
                f"factor={self.factor} phase={self.phase} "
                f"elapsed={time.perf_counter() - self.started:.1f}",
                flush=True,
            )

    def stop(self) -> None:
        self._stop.set()
        self._thread.join(timeout=2.0)


def run_subprocess_with_heartbeat(
    command: list[str], outdir: Path, heartbeat: Heartbeat
) -> float:
    stdout_path = outdir / "transition_stdout.log"
    stderr_path = outdir / "transition_stderr.log"
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH", "")
    prefixes = [str(MODEL_ROOT), str(TOOLS_ROOT)]
    if existing_pythonpath:
        prefixes.append(existing_pythonpath)
    environment["PYTHONPATH"] = os.pathsep.join(prefixes)
    for name in (
        "NUMBA_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
    ):
        environment[name] = "1"
    started = time.perf_counter()
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr:
        process = subprocess.Popen(
            command,
            cwd=ROOT,
            env=environment,
            stdout=stdout,
            stderr=stderr,
            text=True,
        )
        while process.poll() is None:
            time.sleep(5.0)
        return_code = int(process.returncode)
    elapsed = time.perf_counter() - started
    if return_code != 0:
        tail = ""
        try:
            tail = "\n".join(stderr_path.read_text(encoding="utf-8").splitlines()[-30:])
        except OSError:
            pass
        raise RuntimeError(
            f"Dated-transition driver failed with exit code {return_code}:\n{tail}"
        )
    heartbeat.write()
    return elapsed


def validate_transition_output(
    transition_dir: Path,
    contracts: dict[str, Any],
    scaled_theta: dict[str, float],
    preference_change: float,
    factor: float,
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    required_files = (
        "summary.json",
        "target_fit_long.csv",
        "parameter_table.csv",
        "candidate_summary.csv",
    )
    missing = [name for name in required_files if not (transition_dir / name).is_file()]
    if missing:
        raise ContractError(f"Transition driver omitted artifacts: {missing}")
    summary = read_json(transition_dir / "summary.json")
    require_equal(
        summary.get("status"),
        "complete_transition_calibration_panel_task",
        "frontier transition status",
    )
    require_equal(summary.get("target_set"), contracts["scientific_contract"]["target_set"], "frontier target set")
    require_equal(
        summary.get("target_fingerprint"),
        contracts["hashes"]["target_fingerprint"],
        "frontier target fingerprint",
    )
    require_equal(
        dict(summary.get("code_fingerprints") or {}).get("bundle_sha256"),
        contracts["hashes"]["code_bundle_sha256"],
        "frontier code bundle",
    )
    require_equal(
        dict(summary.get("model_profile") or {}).get("name"),
        MODEL_PROFILE,
        "frontier model profile",
    )
    require_equal(summary.get("policy_case"), "none", "frontier policy case")
    if int(summary.get("post_2023_periods", -1)) != 0:
        raise ContractError("Frontier transition extends beyond 2023")
    if int(summary.get("target_count", -1)) != TARGET_COUNT:
        raise ContractError("Frontier transition does not retain all twelve targets")
    validate_exact_renewal_contract(
        dict(summary.get("renewal_accounting_contract") or {}),
        "frontier renewal contract",
    )
    require_close(
        summary.get("old_model_completed_fertility"),
        REPLACEMENT_FERTILITY,
        "frontier old completed fertility",
        5.0e-4,
    )
    require_close(
        dict(summary.get("renewal_accounting_old_state") or {}).get(
            "old_queue_B_over_E"
        ),
        1.0,
        "frontier old B/E",
        5.0e-4,
    )
    best = dict(summary.get("best_candidate") or {})
    compare_theta(dict(best.get("theta") or {}), scaled_theta, "frontier scaled theta")
    require_close(
        float(best["new_psi_child"]) - float(summary["old_psi_child"]),
        preference_change,
        "frontier preference change",
        2.0e-12,
    )
    require_close(
        float(best["theta"]["kappa_fert"])
        / float(contracts["selected_candidate"]["theta"]["kappa_fert"]),
        factor,
        "frontier first-birth kappa factor",
        2.0e-12,
    )
    require_close(
        float(best["theta"]["kappa_fert_continuation"])
        / float(
            contracts["selected_candidate"]["theta"]["kappa_fert_continuation"]
        ),
        factor,
        "frontier continuation kappa factor",
        2.0e-12,
    )
    fit = read_csv(transition_dir / "target_fit_long.csv")
    if len(fit) != TARGET_COUNT or len({row["moment"] for row in fit}) != TARGET_COUNT:
        raise ContractError(
            f"Frontier target-fit table is not complete: rows={len(fit)}"
        )
    loss_from_rows = sum(float(row["loss_contribution"]) for row in fit)
    require_close(
        loss_from_rows,
        best["transition_loss"],
        "frontier target-loss identity",
        2.0e-8,
    )
    return summary, fit


def build_parameter_table(
    transition_dir: Path,
    selected_theta: dict[str, Any],
    scaled_theta: dict[str, Any],
    factor: float,
) -> list[dict[str, Any]]:
    source_rows = read_csv(transition_dir / "parameter_table.csv")
    rows: list[dict[str, Any]] = []
    for row in source_rows:
        name = str(row["parameter"])
        diagnostic_status = "held_fixed_at_final_ridge_selection"
        selected_value: float | None = None
        if name in selected_theta:
            selected_value = float(selected_theta[name])
        elif name == "beta_annual" and "beta" in selected_theta:
            selected_value = float(selected_theta["beta"]) ** 0.25
        if name in {"kappa_fert", "kappa_fert_continuation"}:
            diagnostic_status = "fixed_common_multiplier_diagnostic_not_estimated"
        elif name == "psi_child_2007":
            diagnostic_status = "renormalized_to_external_replacement_2p1"
        elif name in {"psi_child_2023", "psi_child_change_2023"}:
            diagnostic_status = "derived_holding_selected_preference_change_fixed"
        elif str(row.get("is_free_parameter", "")).strip().lower() != "true":
            diagnostic_status = "externally_fixed_or_derived_unchanged"
        rows.append(
            {
                "parameter": name,
                "selected_ridge_value": selected_value,
                "diagnostic_value": float(row["value"]),
                "kappa_factor": float(factor),
                "lower_bound": row.get("lower_bound"),
                "upper_bound": row.get("upper_bound"),
                "transform": row.get("transform"),
                "diagnostic_status": diagnostic_status,
                "near_bound": row.get("near_bound"),
            }
        )
    expected_scaled = {
        name: float(value)
        for name, value in scaled_theta.items()
        if name in {"kappa_fert", "kappa_fert_continuation"}
    }
    observed = {
        row["parameter"]: float(row["diagnostic_value"])
        for row in rows
        if row["parameter"] in expected_scaled
    }
    compare_theta(observed, expected_scaled, "diagnostic parameter table")
    return rows


def solve_closed_endpoint(
    transition_summary: dict[str, Any],
    contracts: dict[str, Any],
    outdir: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]], float]:
    import audit_closed_reproductive_closure as closure
    import run_dynamic_population_transition as calendar
    import run_e5f_open_population_transition as transition
    import run_e5f_post2023_no_policy_continuations as continuation
    import run_e5f_transition_calibration as calibration
    from intergen_eqscale_seq_optimized import parameters as parameter_module

    chain = contracts["chain"]
    model = contracts["model"]
    if chain is None or model is None:
        chain, model = transition.configure_sequential_model()
    best = dict(transition_summary["best_candidate"])
    theta = {str(name): float(value) for name, value in best["theta"].items()}
    _, profile_overrides, _ = calibration.activate_model_profile(MODEL_PROFILE, theta)
    base = closure.make_overrides(chain, theta, nb=120, profile="e5f-floor")
    base.update(profile_overrides)
    base.update(theta)
    base.update(
        parent_dp_waiver=False,
        parent_dp_waiver_phi=1.0,
        parent_dp_waiver_locations=np.array([], dtype=int),
        parent_dp_waiver_owner_rungs=np.array([], dtype=int),
        parent_dp_waiver_birth_state_only=False,
    )
    prototype = parameter_module.setup_parameters()
    prototype = parameter_module.apply_overrides(prototype, profile_overrides)
    prototype = parameter_module.apply_overrides(prototype, theta)
    elasticity = float(np.asarray(prototype.xi_supply, dtype=float).reshape(-1)[0])
    supply_meta = dict(transition_summary["housing_supply"])
    require_equal(
        supply_meta.get("housing_supply_mode"),
        "static-elastic",
        "frontier housing-supply mode",
    )
    old_price = float(supply_meta["retained_date0_asset_price"])
    supply_rule = calendar.HousingSupplyRule(
        str(supply_meta["housing_supply_mode"]),
        old_price,
        float(supply_meta["date0_normalized_housing_stock"]),
        elasticity,
    )
    started = time.perf_counter()
    endpoint, schedule = continuation.solve_closed_stationary_endpoint(
        chain=chain,
        base_overrides=base,
        old_price=old_price,
        new_psi_child=float(best["new_psi_child"]),
        supply_rule=supply_rule,
        price_min_ratio=PRICE_MIN_RATIO,
        price_max_ratio=PRICE_MAX_RATIO,
        grid_points=PRICE_GRID_POINTS,
        outdir=outdir,
    )
    elapsed = time.perf_counter() - started
    return endpoint, schedule, elapsed


def artifact_hashes(outdir: Path, names: Iterable[str]) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for name in names:
        path = outdir / name
        if not path.is_file():
            raise ContractError(f"Required frontier artifact is missing: {path}")
        hashes[name] = file_sha256(path)
    return hashes


def run_factor(args: argparse.Namespace) -> None:
    if args.factor_index is None or args.collect:
        raise ContractError("Factor mode requires --factor-index and forbids --collect")
    if args.run_stage not in {"smoke", "factor"}:
        raise ContractError("Factor execution requires smoke or factor run stage")
    if args.results_root is not None:
        raise ContractError("Factor mode does not use --results-root")
    factor_index = int(args.factor_index)
    factor = factor_for_index(factor_index)
    if args.run_stage == "smoke" and (factor_index != 2 or factor != 0.5):
        raise ContractError("Exact-loop smoke is restricted to factor index 2 (0.5)")
    outdir = args.output_dir.resolve()
    if outdir.exists():
        raise ContractError(f"Refusing to overwrite an existing output path: {outdir}")
    outdir.mkdir(parents=True)
    heartbeat = Heartbeat(outdir, factor_index, factor)
    heartbeat.start()
    started = time.perf_counter()
    try:
        heartbeat.set_phase("validating_final_ridge_contract")
        contracts = validate_inputs(args, require_live_code=True)
        selected = contracts["selected_candidate"]
        center, scaled_theta, preference_change = build_factor_center(selected, factor)
        center_path = outdir / "factor_center.json"
        atomic_write_json(center_path, center)
        center_sha = file_sha256(center_path)
        atomic_write_json(
            outdir / "latest_completed_case.json",
            {
                "status": "running",
                "completed_phase": "input_contracts",
                "factor_index": factor_index,
                "kappa_factor": factor,
                "center_sha256": center_sha,
                "elapsed_seconds": time.perf_counter() - started,
            },
        )

        transition_dir = outdir / "transition"
        heartbeat.set_phase("dated_transition_2007_2023")
        command = [
            sys.executable,
            str(TOOLS_ROOT / "run_e5f_transition_calibration.py"),
            "--source",
            str(contracts["source_path"]),
            "--expected-source-sha256",
            str(args.expected_source_sha256),
            "--expected-target-set",
            str(args.expected_target_set),
            "--expected-target-fingerprint",
            str(args.expected_target_fingerprint),
            "--expected-code-bundle-sha256",
            str(args.expected_code_bundle_sha256),
            "--outdir",
            str(transition_dir),
            "--panel-task-id",
            "1",
            "--panel-size",
            "1",
            "--panel-seed",
            "2026081701",
            "--panel-local-radius",
            "0.02",
            "--panel-design",
            "mixed",
            "--panel-center-json",
            str(center_path),
            "--old-psi-child",
            str(float(selected["old_psi_child"])),
            "--replacement-fertility",
            "2.1",
            "--old-completed-fertility-target",
            "2.1",
            "--post-2023-periods",
            "0",
            "--policy-case",
            "none",
            "--model-profile",
            MODEL_PROFILE,
            "--outside-origin-entry-share",
            str(float(contracts["scientific_contract"]["outside_origin_entry_share"])),
        ]
        transition_seconds = run_subprocess_with_heartbeat(command, outdir, heartbeat)
        transition_summary, fit_rows = validate_transition_output(
            transition_dir,
            contracts,
            scaled_theta,
            preference_change,
            factor,
        )
        shutil.copy2(transition_dir / "target_fit_long.csv", outdir / "target_fit.csv")
        parameter_rows = build_parameter_table(
            transition_dir,
            dict(selected["theta"]),
            scaled_theta,
            factor,
        )
        atomic_write_csv(outdir / "parameters.csv", parameter_rows)
        atomic_write_json(
            outdir / "latest_completed_case.json",
            {
                "status": "running",
                "completed_phase": "dated_transition_2007_2023",
                "factor_index": factor_index,
                "kappa_factor": factor,
                "transition_loss": float(
                    transition_summary["best_candidate"]["transition_loss"]
                ),
                "elapsed_seconds": time.perf_counter() - started,
            },
        )

        heartbeat.set_phase("closed_stationary_25_point_schedule")
        endpoint_dir = outdir / "closed_endpoint"
        endpoint_dir.mkdir()
        endpoint, schedule, endpoint_seconds = solve_closed_endpoint(
            transition_summary, contracts, endpoint_dir
        )
        shutil.copy2(
            endpoint_dir / "closed_stationary_schedule.csv",
            outdir / "closed_stationary_schedule.csv",
        )
        shutil.copy2(
            endpoint_dir / "closed_stationary_endpoint.json",
            outdir / "closed_stationary_endpoint.json",
        )
        if int(endpoint.get("audit_grid_points", PRICE_GRID_POINTS)) != PRICE_GRID_POINTS:
            raise ContractError("Closed endpoint did not use the exact 25-point grid")
        ratios = np.asarray(
            [float(row["queue_B_over_E"]) for row in schedule], dtype=float
        )
        prices = np.asarray([float(row["price_ratio"]) for row in schedule], dtype=float)
        if ratios.size < PRICE_GRID_POINTS or np.any(~np.isfinite(ratios)):
            raise ContractError("Closed endpoint schedule is incomplete or nonfinite")
        maximum_index = int(np.argmax(ratios))
        best_transition = dict(transition_summary["best_candidate"])
        total_seconds = time.perf_counter() - started
        summary = {
            "schema": SCHEMA,
            "status": "complete_fixed_parameter_slope_diagnostic",
            "classification": "fixed_parameter_diagnostic_not_calibration",
            "promotion_eligible": False,
            "policy_case": "none",
            "run_stage": args.run_stage,
            "factor_index": factor_index,
            "kappa_factor": factor,
            "factors_contract": list(FACTORS),
            "selected_ridge_theta": selected["theta"],
            "diagnostic_theta": best_transition["theta"],
            "selected_preference_change_2023": preference_change,
            "diagnostic_old_psi_child": float(transition_summary["old_psi_child"]),
            "diagnostic_new_psi_child": float(best_transition["new_psi_child"]),
            "transition_loss": float(best_transition["transition_loss"]),
            "target_count": TARGET_COUNT,
            "old_completed_fertility": float(
                transition_summary["old_model_completed_fertility"]
            ),
            "old_queue_B_over_E": float(
                transition_summary["renewal_accounting_old_state"][
                    "old_queue_B_over_E"
                ]
            ),
            "closed_root_status": endpoint["status"],
            "usable_closed_root": bool(endpoint.get("usable_closed_root", False)),
            "between_steady_states_label_allowed": bool(
                endpoint.get("between_steady_states_label_allowed", False)
            ),
            "closed_root_price_ratio": (
                float(endpoint["price_ratio"])
                if bool(endpoint.get("usable_closed_root", False))
                else None
            ),
            "maximum_B_over_E": float(ratios[maximum_index]),
            "maximum_B_over_E_price_ratio": float(prices[maximum_index]),
            "minimum_B_over_E": float(np.min(ratios)),
            "renewal_schedule_monotone_nonincreasing": bool(
                np.all(np.diff(ratios[np.argsort(prices)]) <= 1.0e-10)
            ),
            "price_grid_contract": {
                "minimum_ratio": PRICE_MIN_RATIO,
                "maximum_ratio": PRICE_MAX_RATIO,
                "declared_grid_points": PRICE_GRID_POINTS,
                "saved_schedule_rows_including_root_refinement": int(len(schedule)),
            },
            "renewal_accounting_contract": required_renewal_contract(),
            "runtime": {
                "dated_transition_seconds": transition_seconds,
                "closed_endpoint_seconds": endpoint_seconds,
                "total_seconds": total_seconds,
                "old_stationary_solves": int(
                    transition_summary["old_fertility_normalization"][
                        "stationary_solves"
                    ]
                ),
                "dated_transition_model_solves": int(
                    best_transition["model_solve_count"]
                ),
                "closed_endpoint_model_solves": int(len(schedule)),
            },
            "input_provenance": {
                **contracts["hashes"],
                "ridge_report": str(contracts["report_path"]),
                "selected_task_summary": str(contracts["task_path"]),
                "source": str(contracts["source_path"]),
                "factor_center_sha256": center_sha,
            },
            "target_fit": [dict(row) for row in fit_rows],
            "parameters": parameter_rows,
            "closed_stationary_endpoint": endpoint,
        }
        atomic_write_json(outdir / "summary.json", summary)
        artifact_names = (
            "summary.json",
            "factor_center.json",
            "target_fit.csv",
            "parameters.csv",
            "closed_stationary_schedule.csv",
            "closed_stationary_endpoint.json",
            "transition_stdout.log",
            "transition_stderr.log",
        )
        manifest = {
            "schema": "e5f_final_closed_root_slope_frontier_factor_manifest_v1",
            "status": "complete_fixed_parameter_slope_diagnostic_manifest",
            "driver": str(Path(__file__).resolve()),
            "driver_sha256": file_sha256(Path(__file__).resolve()),
            "run_stage": args.run_stage,
            "factor_index": factor_index,
            "kappa_factor": factor,
            "diagnostic_code_contract": {
                name: contracts["hashes"][name]
                for name in (
                    "frontier_driver_sha256",
                    "continuation_driver_sha256",
                    "launcher_sha256",
                )
            },
            "diagnostic_bundle_sha256": contracts["hashes"][
                "diagnostic_bundle_sha256"
            ],
            "artifacts": artifact_hashes(outdir, artifact_names),
        }
        atomic_write_json(outdir / "manifest.json", manifest)
        atomic_write_json(
            outdir / "latest_completed_case.json",
            {
                "status": "complete",
                "completed_phase": "closed_stationary_25_point_schedule",
                "factor_index": factor_index,
                "kappa_factor": factor,
                "transition_loss": summary["transition_loss"],
                "maximum_B_over_E": summary["maximum_B_over_E"],
                "usable_closed_root": summary["usable_closed_root"],
                "elapsed_seconds": total_seconds,
            },
        )
        heartbeat.set_phase("complete")
        print(
            "FINAL_CLOSED_ROOT_SLOPE_FACTOR_COMPLETE "
            f"index={factor_index} factor={factor:g} "
            f"loss={summary['transition_loss']:.9f} "
            f"max_B_over_E={summary['maximum_B_over_E']:.9f} "
            f"root={summary['usable_closed_root']} output={outdir}",
            flush=True,
        )
    except Exception as error:
        atomic_write_json(
            outdir / "failure.json",
            {
                "status": "failed_fixed_parameter_slope_diagnostic",
                "factor_index": factor_index,
                "kappa_factor": factor,
                "error_type": type(error).__name__,
                "error": str(error),
                "traceback": traceback.format_exc(),
                "elapsed_seconds": time.perf_counter() - started,
            },
        )
        heartbeat.set_phase("failed")
        raise
    finally:
        heartbeat.stop()


def collect_frontier(args: argparse.Namespace) -> None:
    if not args.collect or args.factor_index is not None:
        raise ContractError("Collection mode requires --collect and forbids --factor-index")
    if args.run_stage != "collect":
        raise ContractError("Collection mode requires collect run stage")
    if args.results_root is None:
        raise ContractError("Collection mode requires --results-root")
    results_root = args.results_root.resolve()
    outdir = args.output_dir.resolve()
    if outdir.exists():
        raise ContractError(f"Refusing to overwrite an existing report path: {outdir}")
    outdir.mkdir(parents=True)
    contracts = validate_inputs(args, require_live_code=True)
    summaries: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    parameter_rows: list[dict[str, Any]] = []
    schedule_rows: list[dict[str, Any]] = []
    factor_artifacts: dict[str, dict[str, str]] = {}
    for factor_index, factor in enumerate(FACTORS, start=1):
        factor_dir = results_root / f"factor_{factor_index:03d}"
        summary_path = factor_dir / "summary.json"
        manifest_path = factor_dir / "manifest.json"
        summary = read_json(summary_path)
        manifest = read_json(manifest_path)
        require_equal(summary.get("schema"), SCHEMA, f"factor {factor_index} schema")
        require_equal(
            summary.get("status"),
            "complete_fixed_parameter_slope_diagnostic",
            f"factor {factor_index} status",
        )
        if bool(summary.get("promotion_eligible", True)):
            raise ContractError(f"Factor {factor_index} is mislabeled promotion eligible")
        require_equal(summary.get("run_stage"), "factor", "factor run stage")
        require_equal(summary.get("factor_index"), factor_index, "factor index")
        require_close(summary.get("kappa_factor"), factor, "kappa factor", 1.0e-14)
        if tuple(float(value) for value in summary.get("factors_contract", [])) != FACTORS:
            raise ContractError(f"Factor {factor_index} carries the wrong factor contract")
        if int(summary.get("target_count", -1)) != TARGET_COUNT:
            raise ContractError(f"Factor {factor_index} target count is incomplete")
        provenance = dict(summary.get("input_provenance") or {})
        for name, expected in contracts["hashes"].items():
            require_equal(
                provenance.get(name), expected, f"factor {factor_index} provenance {name}"
            )
        require_equal(
            manifest.get("status"),
            "complete_fixed_parameter_slope_diagnostic_manifest",
            f"factor {factor_index} manifest",
        )
        require_equal(manifest.get("run_stage"), "factor", "factor manifest run stage")
        diagnostic_contract = dict(manifest.get("diagnostic_code_contract") or {})
        for name in (
            "frontier_driver_sha256",
            "continuation_driver_sha256",
            "launcher_sha256",
        ):
            require_equal(
                diagnostic_contract.get(name),
                contracts["hashes"][name],
                f"factor {factor_index} manifest {name}",
            )
        require_equal(
            manifest.get("diagnostic_bundle_sha256"),
            contracts["hashes"]["diagnostic_bundle_sha256"],
            f"factor {factor_index} manifest diagnostic bundle",
        )
        for name, expected_hash in dict(manifest.get("artifacts") or {}).items():
            require_equal(
                file_sha256(factor_dir / name),
                expected_hash,
                f"factor {factor_index} artifact hash {name}",
            )
        summaries.append(summary)
        factor_artifacts[f"factor_{factor_index:03d}/summary.json"] = {
            "sha256": file_sha256(summary_path)
        }
        factor_artifacts[f"factor_{factor_index:03d}/manifest.json"] = {
            "sha256": file_sha256(manifest_path)
        }
        for row in read_csv(factor_dir / "target_fit.csv"):
            target_rows.append(
                {"factor_index": factor_index, "kappa_factor": factor, **row}
            )
        for row in read_csv(factor_dir / "parameters.csv"):
            parameter_rows.append(
                {"factor_index": factor_index, "kappa_factor": factor, **row}
            )
        for row in read_csv(factor_dir / "closed_stationary_schedule.csv"):
            schedule_rows.append(
                {"factor_index": factor_index, "kappa_factor": factor, **row}
            )
    if len(target_rows) != len(FACTORS) * TARGET_COUNT:
        raise ContractError(
            f"Collected target table has {len(target_rows)} rows; expected "
            f"{len(FACTORS) * TARGET_COUNT}"
        )
    frontier_rows: list[dict[str, Any]] = []
    for summary in summaries:
        runtime = dict(summary["runtime"])
        frontier_rows.append(
            {
                "factor_index": int(summary["factor_index"]),
                "kappa_factor": float(summary["kappa_factor"]),
                "kappa_fert": float(summary["diagnostic_theta"]["kappa_fert"]),
                "kappa_fert_continuation": float(
                    summary["diagnostic_theta"]["kappa_fert_continuation"]
                ),
                "old_psi_child": float(summary["diagnostic_old_psi_child"]),
                "new_psi_child": float(summary["diagnostic_new_psi_child"]),
                "transition_loss": float(summary["transition_loss"]),
                "maximum_B_over_E": float(summary["maximum_B_over_E"]),
                "maximum_B_over_E_price_ratio": float(
                    summary["maximum_B_over_E_price_ratio"]
                ),
                "closed_root_status": str(summary["closed_root_status"]),
                "usable_closed_root": bool(summary["usable_closed_root"]),
                "closed_root_price_ratio": summary["closed_root_price_ratio"],
                "total_model_solves": int(runtime["old_stationary_solves"])
                + int(runtime["dated_transition_model_solves"])
                + int(runtime["closed_endpoint_model_solves"]),
                "runtime_seconds": float(runtime["total_seconds"]),
            }
        )
    root_rows = [row for row in frontier_rows if bool(row["usable_closed_root"])]
    first_root = max(root_rows, key=lambda row: float(row["kappa_factor"])) if root_rows else None
    atomic_write_csv(outdir / "frontier_summary.csv", frontier_rows)
    atomic_write_csv(outdir / "target_fit_full.csv", target_rows)
    atomic_write_csv(outdir / "parameters_full.csv", parameter_rows)
    atomic_write_csv(outdir / "closed_stationary_schedules_full.csv", schedule_rows)
    report = {
        "schema": COLLECTION_SCHEMA,
        "status": "complete_fixed_parameter_slope_frontier_diagnostic",
        "classification": "fixed_parameter_diagnostic_not_calibration",
        "promotion_eligible": False,
        "policy_case": "none",
        "factors": list(FACTORS),
        "factor_count": len(frontier_rows),
        "target_count_per_factor": TARGET_COUNT,
        "selected_preference_change_2023": float(
            contracts["selected_candidate"]["new_psi_child"]
        )
        - float(contracts["selected_candidate"]["old_psi_child"]),
        "renewal_accounting_contract": required_renewal_contract(),
        "price_grid_contract": {
            "minimum_ratio": PRICE_MIN_RATIO,
            "maximum_ratio": PRICE_MAX_RATIO,
            "declared_grid_points": PRICE_GRID_POINTS,
        },
        "first_root_factor": (
            float(first_root["kappa_factor"]) if first_root is not None else None
        ),
        "first_root_price_ratio": (
            float(first_root["closed_root_price_ratio"])
            if first_root is not None
            else None
        ),
        "frontier": frontier_rows,
        "input_provenance": {
            **contracts["hashes"],
            "ridge_report": str(contracts["report_path"]),
            "selected_task_summary": str(contracts["task_path"]),
            "source": str(contracts["source_path"]),
        },
        "interpretation": (
            "Every row changes only the two fertility-logit scales by a common "
            "predeclared factor, independently re-normalizes the old preference "
            "intercept to 2.1, and holds the selected 2007--2023 preference change "
            "fixed. The rows are fixed-parameter diagnostics, not re-estimated "
            "calibrations and not promotion candidates."
        ),
    }
    atomic_write_json(outdir / "summary.json", report)
    artifact_names = (
        "summary.json",
        "frontier_summary.csv",
        "target_fit_full.csv",
        "parameters_full.csv",
        "closed_stationary_schedules_full.csv",
    )
    manifest = {
        "schema": "e5f_final_closed_root_slope_frontier_report_manifest_v1",
        "status": "complete_fixed_parameter_slope_frontier_manifest",
        "driver": str(Path(__file__).resolve()),
        "driver_sha256": file_sha256(Path(__file__).resolve()),
        "diagnostic_code_contract": {
            name: contracts["hashes"][name]
            for name in (
                "frontier_driver_sha256",
                "continuation_driver_sha256",
                "launcher_sha256",
            )
        },
        "diagnostic_bundle_sha256": contracts["hashes"][
            "diagnostic_bundle_sha256"
        ],
        "factor_artifacts": factor_artifacts,
        "artifacts": artifact_hashes(outdir, artifact_names),
    }
    atomic_write_json(outdir / "manifest.json", manifest)
    print(
        "FINAL_CLOSED_ROOT_SLOPE_FRONTIER_COLLECTED "
        f"factors={len(frontier_rows)} first_root_factor="
        f"{report['first_root_factor']} output={outdir}",
        flush=True,
    )


def main() -> None:
    args = parse_args()
    if args.collect:
        collect_frontier(args)
    else:
        run_factor(args)


if __name__ == "__main__":
    main()
