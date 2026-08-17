from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest
import numpy as np

import run_e5f_post2023_no_policy_continuations as driver


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")


def _contract_fixture(tmp_path: Path) -> tuple[dict, dict, Path, Path, Path, Path]:
    renewal = driver.required_renewal_contract()
    old = {
        "old_entry_flow_E": 0.06,
        "old_queue_mature_flow_B": 0.06,
        "outside_flow_M": 0.01014,
        "outside_origin_entry_share": 0.169,
        "retention_rho": 0.831,
    }
    task_best = {
        "candidate": "task_001",
        "theta": {"beta": 0.99},
        "new_psi_child": -0.1,
        "transition_loss": 12.0,
        "policy_case": "none",
        "post_2023_periods": 0,
    }
    report_best = {
        **task_best,
        "old_psi_child": 0.2,
        "target_fingerprint": "target-hash",
        "code_fingerprints": {"bundle_sha256": "code-hash"},
        "source_sha256": "placeholder",
        "renewal_accounting_contract": renewal,
        "renewal_accounting_old_state": old,
        "outside_origin_entry_share": 0.169,
        "valid": True,
    }
    source_path = tmp_path / "source.json"
    source_path.write_text("source\n", encoding="utf-8")
    source_hash = driver.file_sha256(source_path)
    report_best["source_sha256"] = source_hash
    report = {
        "status": "complete",
        "best_candidate": report_best,
        "code_fingerprints": {"bundle_sha256": "code-hash"},
        "renewal_accounting_contract": renewal,
        "renewal_accounting_old_state": old,
        "outside_origin_entry_share": 0.169,
    }
    task = {
        "status": "complete_transition_calibration_panel_task",
        "best_candidate": task_best,
        "old_psi_child": 0.2,
        "target_fingerprint": "target-hash",
        "code_fingerprints": {"bundle_sha256": "code-hash"},
        "source_sha256": source_hash,
        "renewal_accounting_contract": renewal,
        "renewal_accounting_old_state": old,
        "outside_origin_entry_share": 0.169,
        "policy_case": "none",
        "post_2023_periods": 0,
    }
    report_path = tmp_path / "report.json"
    task_path = tmp_path / "task.json"
    case_dir = tmp_path / "task_001"
    case_path = case_dir / "transition_path.csv"
    _write_json(report_path, report)
    _write_json(task_path, task)
    case_dir.mkdir(parents=True)
    case_path.write_text("period,asset_price\n0,1.0\n", encoding="utf-8")
    return report, task, report_path, task_path, case_dir, source_path


def _validate_fixture(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutate=None,
) -> dict:
    report, task, report_path, task_path, case_dir, source_path = _contract_fixture(tmp_path)
    if mutate is not None:
        mutate(report, task)
        _write_json(report_path, report)
        _write_json(task_path, task)
    monkeypatch.setattr(
        driver, "e5_target_system", lambda: SimpleNamespace(fingerprint="target-hash")
    )
    monkeypatch.setattr(
        driver.calibration,
        "code_fingerprint_contract",
        lambda model: {"bundle_sha256": "code-hash", "files": {}},
    )
    return driver.validate_input_contracts(
        report_path=report_path,
        task_summary_path=task_path,
        case_dir=case_dir,
        case_transition_path=None,
        source_path=source_path,
        expected_report_sha256=driver.file_sha256(report_path),
        expected_task_summary_sha256=driver.file_sha256(task_path),
        expected_case_transition_sha256=driver.file_sha256(
            case_dir / "transition_path.csv"
        ),
        expected_source_sha256=driver.file_sha256(source_path),
        expected_target_fingerprint="target-hash",
        expected_code_bundle_sha256="code-hash",
        expected_renewal_contract_sha256=driver.canonical_json_sha256(
            driver.required_renewal_contract()
        ),
        expected_scientific_contract_sha256=None,
        expected_selection_sha256=None,
        chain=object(),
        model=object(),
    )


def test_contract_fixture_passes_all_fail_closed_gates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    validated = _validate_fixture(tmp_path, monkeypatch)
    assert validated["report_best"]["candidate"] == "task_001"
    assert validated["hashes"]["target_fingerprint"] == "target-hash"
    assert validated["hashes"]["code_bundle_sha256"] == "code-hash"


def test_report_byte_mutation_is_rejected_before_use(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, _, report_path, task_path, case_dir, source_path = _contract_fixture(tmp_path)
    expected_report_hash = driver.file_sha256(report_path)
    report_path.write_text(report_path.read_text() + " ", encoding="utf-8")
    monkeypatch.setattr(
        driver, "e5_target_system", lambda: SimpleNamespace(fingerprint="target-hash")
    )
    monkeypatch.setattr(
        driver.calibration,
        "code_fingerprint_contract",
        lambda model: {"bundle_sha256": "code-hash", "files": {}},
    )
    with pytest.raises(RuntimeError, match="Selected report hash gate failed"):
        driver.validate_input_contracts(
            report_path=report_path,
            task_summary_path=task_path,
            case_dir=case_dir,
            case_transition_path=None,
            source_path=source_path,
            expected_report_sha256=expected_report_hash,
            expected_task_summary_sha256=driver.file_sha256(task_path),
            expected_case_transition_sha256=driver.file_sha256(
                case_dir / "transition_path.csv"
            ),
            expected_source_sha256=driver.file_sha256(source_path),
            expected_target_fingerprint="target-hash",
            expected_code_bundle_sha256="code-hash",
            expected_renewal_contract_sha256=driver.canonical_json_sha256(
                driver.required_renewal_contract()
            ),
            expected_scientific_contract_sha256=None,
            expected_selection_sha256=None,
            chain=object(),
            model=object(),
        )


def test_final_ridge_report_schema_is_consumed_without_a_task_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "source.json"
    source.write_text("source\n", encoding="utf-8")
    source_hash = driver.file_sha256(source)
    renewal = driver.required_renewal_contract()
    old = {
        "old_entry_flow_E": 0.06,
        "old_queue_mature_flow_B": 0.06,
        "outside_flow_M": 0.01014,
        "outside_origin_entry_share": 0.169,
        "retention_rho": 0.831,
    }
    scientific = {
        "source": str(source),
        "source_sha256": source_hash,
        "target_set": "target-set",
        "target_fingerprint": "target-hash",
        "code_fingerprints": {"bundle_sha256": "code-hash"},
        "model_profile": {"name": "e5f-income-entry"},
        "renewal_accounting_contract": renewal,
        "outside_origin_entry_share": 0.169,
        "policy_case": "none",
        "post_2023_periods": 0,
    }
    scientific_hash = driver.canonical_json_sha256(scientific)
    report = {
        "schema": "e5f_transition_ridge_refinement_report_v1",
        "status": "complete_refinement_with_two_independent_identity_repeats",
        "promotion_eligible": True,
        "repeat_gate": {
            "required": 2,
            "completed": 2,
            "identity_gate_passed": True,
            "numerical_identity_gate_passed": True,
        },
        "scientific_contract": scientific,
        "scientific_contract_sha256": scientific_hash,
        "selection_sha256": "selection-hash",
        "selected_candidate_sha256": "candidate-hash",
        "plan_sha256": "plan-hash",
        "dated_contract_sha256": "dated-hash",
        "renewal_contract_sha256": driver.canonical_json_sha256(renewal),
        "best_candidate": {
            "candidate": "task_001",
            "theta": {"beta": 0.99},
            "new_psi_child": -0.1,
            "old_psi_child": 0.2,
            "transition_loss": 12.0,
            "candidate_sha256": "candidate-hash",
            "renewal_accounting_contract": renewal,
            "renewal_accounting_old_state": old,
            "valid": True,
        },
    }
    report_path = tmp_path / "summary.json"
    _write_json(report_path, report)
    case_path = tmp_path / "selected_transition_path.csv"
    case_path.write_text("period,asset_price\n0,1.0\n", encoding="utf-8")
    monkeypatch.setattr(
        driver, "e5_target_system", lambda: SimpleNamespace(fingerprint="target-hash")
    )
    monkeypatch.setattr(
        driver.calibration,
        "code_fingerprint_contract",
        lambda model: {"bundle_sha256": "code-hash", "files": {}},
    )
    validated = driver.validate_input_contracts(
        report_path=report_path,
        task_summary_path=None,
        case_dir=None,
        case_transition_path=None,
        source_path=source,
        expected_report_sha256=driver.file_sha256(report_path),
        expected_task_summary_sha256=None,
        expected_case_transition_sha256=driver.file_sha256(case_path),
        expected_source_sha256=source_hash,
        expected_target_fingerprint="target-hash",
        expected_code_bundle_sha256="code-hash",
        expected_renewal_contract_sha256=driver.canonical_json_sha256(renewal),
        expected_scientific_contract_sha256=scientific_hash,
        expected_selection_sha256="selection-hash",
        chain=object(),
        model=object(),
    )
    assert validated["raw_report"]["schema"] == "e5f_transition_ridge_refinement_report_v1"
    assert validated["paths"]["selected_task_summary"] is None
    assert validated["hashes"]["selection_sha256"] == "selection-hash"


def test_policy_or_post2023_calibration_artifact_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def mutate(report, task) -> None:
        report["best_candidate"]["policy_case"] = "dependent-child-ltv95"
        task["best_candidate"]["policy_case"] = "dependent-child-ltv95"
        task["policy_case"] = "dependent-child-ltv95"

    with pytest.raises(RuntimeError, match="scope gate"):
        _validate_fixture(tmp_path, monkeypatch, mutate=mutate)


def test_renewal_contract_is_exact_and_uses_births_over_2p1() -> None:
    contract = driver.required_renewal_contract()
    driver.validate_renewal_contract(contract)
    assert contract["replacement_fertility"] == 2.1
    assert math.isclose(
        contract["effective_birth_to_household_conversion"],
        1.0 / 2.1,
        rel_tol=0.0,
        abs_tol=1e-15,
    )
    corrupted = copy.deepcopy(contract)
    corrupted["effective_birth_to_household_conversion"] = 0.5
    with pytest.raises(RuntimeError, match="effective_birth_to_household_conversion"):
        driver.validate_renewal_contract(corrupted)


def test_outside_share_only_decomposes_old_replacement_entry() -> None:
    audit = driver.outside_share_invariance_audit(0.06, 0.06, 0.247)
    assert audit["status"] == "passed"
    assert audit["maximum_absolute_identity_residual"] <= 1e-15
    assert any(
        row["outside_origin_entry_share"] == 0.247 for row in audit["rows"]
    )
    assert "not identified" in audit["interpretation"]


def test_paired_2023_audit_requires_exact_state_identity() -> None:
    left = {key: float(index) for index, key in enumerate(driver.HISTORY_STATE_COLUMNS)}
    right = dict(left)
    paths = {"closed": [left], "open": [right]}
    assert (
        driver.paired_history_state_gap(paths, driver.HISTORY_STATE_COLUMNS)[
            "maximum_absolute_gap"
        ]
        == 0.0
    )
    right["asset_price"] += 1e-12
    assert (
        driver.paired_history_state_gap(paths, driver.HISTORY_STATE_COLUMNS)[
            "maximum_absolute_gap"
        ]
        > 0.0
    )


def test_closed_schedule_does_not_invent_a_root() -> None:
    no_root = [
        {"asset_price": 0.1, "renewal_residual_ratio": -0.2},
        {"asset_price": 1.0, "renewal_residual_ratio": -0.1},
        {"asset_price": 2.0, "renewal_residual_ratio": -0.05},
    ]
    assert driver.closed_schedule_brackets(no_root) == []
    one_root = [
        {"asset_price": 0.1, "renewal_residual_ratio": 0.1},
        {"asset_price": 1.0, "renewal_residual_ratio": -0.1},
        {"asset_price": 2.0, "renewal_residual_ratio": -0.2},
    ]
    assert driver.closed_schedule_brackets(one_root) == [(0.1, 1.0)]


def test_open_endpoint_verifies_scale_identity() -> None:
    M = 0.02
    rho = 0.8
    B = 0.05
    E = 0.06
    denominator = E - rho * B
    scale = M / denominator
    endpoint = {
        "status": "complete",
        "total_mass": 1.0,
        "stationary_population_scale": scale,
        "queue_mature_entrant_flow_B": B,
        "renewal_denominator": denominator,
        "fixed_stock_relative_market_gap": 1e-7,
    }
    verified = driver.verify_open_endpoint(endpoint, M, rho)
    assert verified["usable_open_endpoint"]
    assert abs(verified["renewal_identity_residual"]) <= 1e-15
    assert "not identified" in verified["interpretation"]


def test_json_output_replaces_nonfinite_values_with_null() -> None:
    assert driver.jsonable({"x": math.nan, "y": math.inf}) == {"x": None, "y": None}


def test_matched_state_hash_pins_shape_and_values() -> None:
    state = np.arange(12.0).reshape(3, 4)
    assert driver.array_sha256(state) == driver.array_sha256(state.copy())
    changed = state.copy()
    changed[0, 0] = 1e-12
    assert driver.array_sha256(state) != driver.array_sha256(changed)
    assert driver.array_sha256(state) != driver.array_sha256(state.reshape(2, 6))
