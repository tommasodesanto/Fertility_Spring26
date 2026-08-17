#!/usr/bin/env python3
"""Build the canonical E5F transition-calibration and funded-policy report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.e5f_floor_profile import E5F_DOMAIN
from intergen_eqscale_seq_optimized.e5f_income_entry_profile import (
    E5F_INCOME_ENTRY_DOMAIN,
)
from intergen_eqscale_seq_optimized.e5_profile import e5_target_system


BASELINE_MODEL_PROFILE = "e5f-floor"
REPAIRED_MODEL_PROFILE = "e5f-income-entry"
EXPECTED_TARGET_COUNT = 12
EXPECTED_FUNDED_PERIODS = 45
EXPECTED_FUNDED_ENDPOINT_YEAR = 2183
EXPECTED_FUNDED_CASES = (
    "unrebated_tax1_status_quo",
    "rebated_tax1_baseline",
    "rebated_tax2",
    "rebated_tax2_grant0p4_Hge6",
)
EXPECTED_FUNDED_TOOL_CODE_FILES = (
    "code/model/tools/run_e5f_transition_funded_policy.py",
    "code/model/tools/run_e5f_transition_calibration.py",
    "code/model/tools/run_e5f_open_population_transition.py",
    "code/model/tools/run_dynamic_population_transition.py",
    "code/model/tools/audit_closed_reproductive_closure.py",
)
EXPECTED_FUNDED_CODE_BUNDLE_FILES = EXPECTED_FUNDED_TOOL_CODE_FILES + tuple(
    path.relative_to(ROOT).as_posix()
    for path in sorted(
        (MODEL_ROOT / "intergen_eqscale_seq_optimized").glob("*.py"),
        key=lambda item: item.name,
    )
)
FUNDED_CASE_CONVENTIONS = {
    "unrebated_tax1_status_quo": (0.01, False, False),
    "rebated_tax1_baseline": (0.01, False, True),
    "rebated_tax2": (0.02, False, True),
    "rebated_tax2_grant0p4_Hge6": (0.02, True, True),
}
FUNDED_CASE_LABELS = {
    "unrebated_tax1_status_quo": "Status quo: 1% tax, no rebate",
    "rebated_tax1_baseline": "Funded 1% tax",
    "rebated_tax2": "Funded 2% tax",
    "rebated_tax2_grant0p4_Hge6": "Funded 2% tax + purchase grant",
}
FUNDED_CASE_COLORS = {
    "unrebated_tax1_status_quo": "#1f1f1f",
    "rebated_tax1_baseline": "#377eb8",
    "rebated_tax2": "#e41a1c",
    "rebated_tax2_grant0p4_Hge6": "#4daf4a",
}
FUNDED_CASE_LINESTYLES = {
    "unrebated_tax1_status_quo": "-",
    "rebated_tax1_baseline": "--",
    "rebated_tax2": "-.",
    "rebated_tax2_grant0p4_Hge6": ":",
}
FUNDED_PATH_FINITE_FIELDS = (
    "period",
    "year",
    "years_from_start",
    "psi_child",
    "asset_price",
    "asset_price_index",
    "housing_user_cost",
    "adult_population",
    "population_index",
    "mean_adult_age",
    "birth_children",
    "birth_children_raw_explicit_states",
    "top_bin_entry_birth_flow",
    "birth_children_topcode_adjusted",
    "births_per_adult",
    "topcode_adjusted_births_per_adult",
    "births_over_entry",
    "topcode_adjusted_births_over_entry",
    "entry_flow_E",
    "mature_children_raw",
    "effective_mature_entrant_flow_B",
    "raw_state_scheduled_mature_entrant_flow_B",
    "model_state_same_period_mature_flow_B",
    "birth_queue_new_potential_flow_B",
    "birth_queue_new_raw_state_potential_flow_B",
    "dual_clock_raw_flow_gap_percent",
    "outside_entry_flow_M",
    "outside_entry_flow_share_of_initial_mass",
    "retained_mature_entrants",
    "entrant_flow_next",
    "youngest_age_mass_next",
    "census_age_bridge_net_residual",
    "census_age_bridge_absolute_group_residual",
    "census_age_bridge_maximum_target_gap",
    "adult_deaths",
    "housing_demand",
    "housing_demand_per_adult",
    "housing_supply",
    "relative_market_residual",
    "owner_rate",
    "dependent_child_owner_rate",
    "parent_owner_rate",
    "childless_owner_rate",
    "mass_accounting_residual",
    "feasibility_frontier_projection_mass",
    "min_distribution_mass",
    "nonfinite_distribution_count",
    "lump_sum_transfer_period_units",
    "annual_property_tax_rate",
    "property_tax_revenue",
    "purchase_grant_recipient_mass",
    "purchase_grant_outlays",
    "lump_sum_transfer_outlays",
    "government_budget_residual",
    "scaled_government_budget_residual",
    "balanced_budget_implied_transfer",
    "balanced_budget_transfer_gap",
    "fiscal_root_evaluations",
)
FUNDED_PATH_JSON_FINITE_FIELDS = (
    "birth_queue_scheduled_flows",
    "birth_queue_raw_state_scheduled_flows",
)
OPTIONAL_PATH_NONFINITE_AS_MISSING = ("population_target_gap",)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--selected-report",
        type=Path,
        required=True,
        help="Completed collector summary used to launch the funded-policy run.",
    )
    parser.add_argument(
        "--round-dir",
        type=Path,
        action="append",
        default=[],
        help="Optional calibration-panel result directory for the search trace.",
    )
    parser.add_argument("--repeat-dir", type=Path, action="append", default=[])
    parser.add_argument(
        "--funded-policy-dir",
        type=Path,
        required=True,
        help="Completed four-case funded-policy packet.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def live_calibration_code_fingerprints() -> dict[str, Any]:
    """Rebuild the calibration driver's deterministic source contract."""
    tools_root = MODEL_ROOT / "tools"
    package_root = MODEL_ROOT / "intergen_eqscale_seq_optimized"
    paths = {
        "transition_calibration_driver": tools_root
        / "run_e5f_transition_calibration.py",
        "calendar_operator": tools_root / "run_dynamic_population_transition.py",
        "e5f_transition_operator": tools_root
        / "run_e5f_open_population_transition.py",
        "closure_audit": tools_root / "audit_closed_reproductive_closure.py",
    }
    paths.update(
        {
            f"sequential_package/{path.name}": path
            for path in sorted(package_root.glob("*.py"))
        }
    )
    if not paths or any(not path.is_file() for path in paths.values()):
        missing = [str(path) for path in paths.values() if not path.is_file()]
        raise FileNotFoundError(f"Live calibration code contract is incomplete: {missing}")
    files = {name: sha256(path) for name, path in paths.items()}
    encoded = json.dumps(files, sort_keys=True, separators=(",", ":")).encode()
    return {
        "files": files,
        "bundle_sha256": hashlib.sha256(encoded).hexdigest(),
    }


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def require_sha256(value: Any, label: str) -> str:
    result = str(value)
    if len(result) != 64 or any(character not in "0123456789abcdef" for character in result):
        raise RuntimeError(f"Malformed {label} SHA-256: {result!r}")
    return result


def require_finite(value: Any, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"Nonfinite {label}: {value!r}")
    return result


def assert_finite_tree(value: Any, label: str) -> None:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return
    if isinstance(value, (int, float)):
        require_finite(value, label)
        return
    if isinstance(value, dict):
        for key, item in value.items():
            assert_finite_tree(item, f"{label}.{key}")
        return
    if isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            assert_finite_tree(item, f"{label}[{index}]")
        return
    raise RuntimeError(f"Unsupported value in {label}: {type(value).__name__}")


def close(
    observed: Any,
    expected: Any,
    *,
    absolute: float = 1.0e-12,
    relative: float = 1.0e-12,
) -> bool:
    return math.isclose(
        require_finite(observed, "observed identity value"),
        require_finite(expected, "expected identity value"),
        rel_tol=relative,
        abs_tol=absolute,
    )


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
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
        writer.writerows(rows)
    temporary.replace(path)


def consistent_value(
    payload: dict[str, Any], best: dict[str, Any], name: str
) -> Any:
    values = [item[name] for item in (best, payload) if item.get(name) is not None]
    if not values:
        raise RuntimeError(f"Selected report omits required contract field {name}")
    encoded = [canonical_json(value) for value in values]
    if len(set(encoded)) != 1:
        raise RuntimeError(f"Selected report has conflicting values for {name}")
    return values[0]


def contract_from_payload(
    payload: dict[str, Any], best: dict[str, Any]
) -> dict[str, Any]:
    profile = consistent_value(payload, best, "model_profile")
    if not isinstance(profile, dict):
        raise RuntimeError("Production reports must contain full model_profile metadata")
    profile = dict(profile)
    profile_name = str(profile.get("name", ""))
    if profile_name not in {BASELINE_MODEL_PROFILE, REPAIRED_MODEL_PROFILE}:
        raise RuntimeError(f"Unsupported model profile: {profile_name!r}")
    income_gates = consistent_value(payload, best, "income_profile_gates")
    if not isinstance(income_gates, dict):
        raise RuntimeError("Production reports must contain income_profile_gates")
    income_gates = dict(income_gates)
    code = consistent_value(payload, best, "code_fingerprints")
    if not isinstance(code, dict):
        raise RuntimeError("Production reports must contain calibration code_fingerprints")
    code = dict(code)
    require_sha256(code.get("bundle_sha256"), "calibration code-bundle")
    files = code.get("files")
    if not isinstance(files, dict) or not files:
        raise RuntimeError("Calibration code contract must contain per-file fingerprints")
    for filename, fingerprint in files.items():
        require_sha256(fingerprint, f"calibration file {filename}")
    target_count = int(consistent_value(payload, best, "target_count"))
    if target_count != EXPECTED_TARGET_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_TARGET_COUNT} target moments, found {target_count}"
        )
    outside_share = float(
        consistent_value(payload, best, "outside_origin_entry_share")
    )
    if not math.isfinite(outside_share) or not 0.0 < outside_share < 1.0:
        raise RuntimeError(f"Invalid outside-origin entry share: {outside_share}")
    source = str(consistent_value(payload, best, "source"))
    source_sha = require_sha256(
        consistent_value(payload, best, "source_sha256"), "source"
    )
    target_set = str(consistent_value(payload, best, "target_set"))
    target_fingerprint = require_sha256(
        consistent_value(payload, best, "target_fingerprint"), "target-system"
    )
    if not target_set or not target_fingerprint:
        raise RuntimeError("Target-set contract is empty")
    if profile_name == REPAIRED_MODEL_PROFILE:
        required_profile = (
            "profile_id",
            "permanent_income_log_variance",
            "income_state_count",
            "first_birth_fixed_cost",
            "first_birth_fixed_cost_semantics",
        )
        missing = [name for name in required_profile if name not in profile]
        if missing:
            raise RuntimeError(
                "Repaired model-profile contract omits " + ", ".join(missing)
            )
        if int(profile["income_state_count"]) != 15:
            raise RuntimeError("Repaired profile must use the measured 15-state income grid")
        if not bool(income_gates.get("permanent_income_levels_enabled", False)):
            raise RuntimeError("Repaired profile did not activate permanent income levels")
        if int(income_gates.get("income_state_count", -1)) != 15:
            raise RuntimeError("Repaired income-state gate is not 15")
        stationary_gap = float(
            income_gates.get("stationary_weight_max_abs_gap", math.inf)
        )
        if not math.isfinite(stationary_gap) or stationary_gap > 1.0e-12:
            raise RuntimeError("Repaired stationary income-weight gate failed")
    return {
        "source": source,
        "source_sha256": source_sha,
        "target_set": target_set,
        "target_fingerprint": target_fingerprint,
        "target_count": target_count,
        "model_profile": profile,
        "income_profile_gates": income_gates,
        "code_fingerprints": code,
        "outside_origin_entry_share": outside_share,
    }


def contract_identity(contract: dict[str, Any]) -> str:
    fields = {
        name: contract[name]
        for name in (
            "source_sha256",
            "target_set",
            "target_fingerprint",
            "target_count",
            "model_profile",
            "income_profile_gates",
            "code_fingerprints",
            "outside_origin_entry_share",
        )
    }
    return canonical_json(fields)


def assert_same_contract(
    expected: dict[str, Any], observed: dict[str, Any], label: str
) -> None:
    if contract_identity(expected) != contract_identity(observed):
        raise RuntimeError(f"{label} does not match the selected calibration contract")


def candidate_records(round_dirs: list[Path]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for round_index, directory in enumerate(round_dirs, start=1):
        for summary_path in sorted(directory.glob("task_*/summary.json")):
            summary = read_json(summary_path)
            candidate = dict(summary["best_candidate"])
            panel_raw = summary.get("panel_design")
            if not isinstance(panel_raw, dict):
                raise RuntimeError(f"Missing panel_design in {summary_path}")
            panel = dict(panel_raw)
            contract = contract_from_payload(summary, candidate)
            candidate.update(
                {
                    **contract,
                    "contract": contract,
                    "round": round_index,
                    "round_name": directory.name,
                    "task_dir": str(summary_path.parent),
                    "old_psi_child": summary["old_psi_child"],
                    "old_completed_fertility": summary[
                        "old_model_completed_fertility"
                    ],
                    "theta_json": canonical_json(candidate["theta"]),
                    "task_id": int(panel["task_id"]),
                    "design": str(panel["design"]),
                }
            )
            records.append(candidate)
    return records


def portable_provenance(
    report_path: Path,
    payload: dict[str, Any],
    best: dict[str, Any],
) -> dict[str, Any]:
    sources = [best, payload]
    sidecar_path = report_path.parent / "calibration_provenance.json"
    if sidecar_path.is_file():
        sources.append(read_json(sidecar_path))
    result: dict[str, Any] = {}
    for field in (
        "population_bridge",
        "population_validation_status",
        "target_measurements",
    ):
        values = [source[field] for source in sources if source.get(field) is not None]
        if not values:
            raise RuntimeError(
                "Portable calibration provenance omits "
                f"{field}; record it in the collector report or "
                f"{sidecar_path.name}"
            )
        encoded = [canonical_json(value) for value in values]
        if len(set(encoded)) != 1:
            raise RuntimeError(f"Portable calibration provenance conflicts at {field}")
        result[field] = values[0]
    if not isinstance(result["target_measurements"], (dict, list)) or not result[
        "target_measurements"
    ]:
        raise RuntimeError("Portable target-measurement provenance is empty or malformed")
    for field in ("population_bridge", "population_validation_status"):
        value = result[field]
        if not isinstance(value, (str, dict)) or not value:
            raise RuntimeError(f"Portable calibration provenance is malformed at {field}")
    assert_finite_tree(result, "portable calibration provenance")
    return result


def load_selected_report(
    path: Path,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
]:
    payload = read_json(path)
    collector_status = str(payload.get("status", ""))
    if collector_status not in {"complete", "partial"}:
        raise RuntimeError(f"Selected collector report is not complete: {path}")
    expected = int(payload.get("expected_tasks", -1))
    completed = int(payload.get("completed_tasks", -1))
    failures = payload.get("failed_or_missing_tasks")
    if expected < 1 or completed < 1 or not isinstance(failures, list):
        raise RuntimeError("Selected collector report has an incomplete task ledger")
    if collector_status == "complete" and (
        completed != expected or failures
    ):
        raise RuntimeError("Complete collector report does not certify every task")
    if collector_status == "partial":
        if completed + len(failures) != expected:
            raise RuntimeError("Partial collector report does not account for every task")
        bad_reasons = [
            row.get("reason")
            for row in failures
            if row.get("reason") != "classified_invalid_candidate"
        ]
        if bad_reasons:
            raise RuntimeError(
                "Collector report has missing or malformed tasks rather than "
                f"classified invalid candidates: {bad_reasons}"
            )
    best_raw = payload.get("best_candidate")
    if not isinstance(best_raw, dict) or not isinstance(best_raw.get("theta"), dict):
        raise RuntimeError(f"Selected collector report has no best candidate: {path}")
    best = dict(best_raw)
    if not bool(best.get("valid", True)):
        raise RuntimeError("Selected collector winner is not marked valid")
    assert_finite_tree(best, "selected calibration winner")
    contract = contract_from_payload(payload, best)
    if not (path.parent / "best_target_fit.csv").is_file():
        raise FileNotFoundError(
            "Collector report omits the required portable target-fit sidecar: "
            f"{path.parent / 'best_target_fit.csv'}"
        )
    old_psi = best.get("old_psi_child")
    old_fertility = best.get("old_completed_fertility")
    if old_psi is None or old_fertility is None:
        raise RuntimeError("Collector report omits the old-fertility normalization")
    best.update(
        **contract,
        contract=contract,
        old_psi_child=float(old_psi),
        old_completed_fertility=float(old_fertility),
        theta_json=canonical_json(best["theta"]),
    )
    provenance = portable_provenance(path, payload, best)
    return payload, best, contract, provenance


def target_fit_rows(path: Path) -> list[dict[str, Any]]:
    csv_path = path / "target_fit_long.csv" if path.is_dir() else path
    raw = read_csv(csv_path)
    if len(raw) != EXPECTED_TARGET_COUNT:
        raise RuntimeError(
            f"{csv_path}: expected {EXPECTED_TARGET_COUNT} target rows, found {len(raw)}"
        )
    required = (
        "moment",
        "target",
        "model",
        "gap",
        "weight",
        "loss_contribution",
        "standardized_gap",
    )
    rows: list[dict[str, Any]] = []
    for row in raw:
        missing = [field for field in required if field not in row]
        if missing:
            raise RuntimeError(f"{csv_path}: target row omits {missing}")
        converted = {"moment": str(row["moment"])}
        converted.update({field: float(row[field]) for field in required[1:]})
        if not converted["moment"] or not all(
            math.isfinite(float(converted[field])) for field in required[1:]
        ):
            raise RuntimeError(f"{csv_path}: malformed target row {row}")
        target = float(converted["target"])
        model = float(converted["model"])
        weight = float(converted["weight"])
        gap = model - target
        contribution = weight * gap**2
        standardized = math.copysign(math.sqrt(contribution), gap)
        if weight <= 0.0:
            raise RuntimeError(f"{csv_path}: target weight must be positive")
        algebra = (
            ("gap", converted["gap"], gap),
            ("loss_contribution", converted["loss_contribution"], contribution),
            ("standardized_gap", converted["standardized_gap"], standardized),
        )
        for field, observed, expected in algebra:
            if not close(observed, expected, absolute=2.0e-12, relative=2.0e-12):
                raise RuntimeError(
                    f"{csv_path}: {field} identity fails for {converted['moment']}"
                )
        rows.append(converted)
    moments = [str(row["moment"]) for row in rows]
    if len(set(moments)) != EXPECTED_TARGET_COUNT:
        raise RuntimeError(f"{csv_path}: target moment names are duplicated")
    return rows


def recomputed_target_fingerprint(rows: list[dict[str, Any]]) -> str:
    payload = [
        [str(row["moment"]), float(row["target"]), float(row["weight"])]
        for row in rows
    ]
    encoded = json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def validate_selected_target_fit(
    rows: list[dict[str, Any]], best: dict[str, Any], contract: dict[str, Any]
) -> float:
    fingerprint = recomputed_target_fingerprint(rows)
    if fingerprint != contract["target_fingerprint"]:
        raise RuntimeError(
            "Target rows do not reproduce the selected target fingerprint: "
            f"rows={fingerprint}, report={contract['target_fingerprint']}"
        )
    loss = math.fsum(float(row["loss_contribution"]) for row in rows)
    if not close(
        loss,
        best["transition_loss"],
        absolute=2.0e-10,
        relative=2.0e-12,
    ):
        raise RuntimeError(
            "Target loss contributions do not reproduce the selected transition loss"
        )
    live_targets = e5_target_system()
    if (
        live_targets.name != contract["target_set"]
        or live_targets.count != EXPECTED_TARGET_COUNT
        or live_targets.fingerprint != fingerprint
    ):
        raise RuntimeError("Selected target rows differ from the live target-system contract")
    return loss


def target_fit_contract(rows: list[dict[str, Any]]) -> list[tuple[str, float, float]]:
    return [
        (str(row["moment"]), float(row["target"]), float(row["weight"]))
        for row in rows
    ]


def assert_same_target_fit_contract(
    expected: list[dict[str, Any]], observed: list[dict[str, Any]], label: str
) -> None:
    left = target_fit_contract(expected)
    right = target_fit_contract(observed)
    if len(left) != len(right):
        raise RuntimeError(f"{label} target row count differs")
    for expected_row, observed_row in zip(left, right, strict=True):
        if expected_row[0] != observed_row[0] or not all(
            math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-14)
            for a, b in zip(expected_row[1:], observed_row[1:], strict=True)
        ):
            raise RuntimeError(f"{label} target contract differs at {expected_row[0]}")


def fit_lookup(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    return {
        str(row["moment"]): {
            key: float(row[key])
            for key in (
                "target",
                "model",
                "gap",
                "weight",
                "loss_contribution",
                "standardized_gap",
            )
        }
        for row in rows
    }


def transition_domain(best: dict[str, Any]) -> tuple[tuple[str, float, float, str], ...]:
    profile = dict(best["model_profile"])
    name = str(profile["name"])
    stationary_domain = (
        E5F_INCOME_ENTRY_DOMAIN if name == REPAIRED_MODEL_PROFILE else E5F_DOMAIN
    )
    terminal = (
        ("psi_child_change_2023", -1.50, 0.20, "asinh")
        if name == REPAIRED_MODEL_PROFILE
        else ("psi_child_2023", -1.25, 0.20, "asinh")
    )
    domain = tuple(row for row in stationary_domain if row[0] != "psi_child") + (
        terminal,
    )
    expected_size = 10 if name == REPAIRED_MODEL_PROFILE else 9
    if len(domain) != expected_size:
        raise RuntimeError(
            f"{name} transition domain has {len(domain)} coordinates, expected {expected_size}"
        )
    if name == REPAIRED_MODEL_PROFILE and "first_birth_fixed_cost" not in best["theta"]:
        raise RuntimeError("Repaired candidate omits first_birth_fixed_cost")
    return domain


def validate_solved_model_profile(
    metadata: dict[str, Any], best: dict[str, Any]
) -> None:
    profile = dict(best["model_profile"])
    active_profile = metadata.get("active_model_profile")
    if not isinstance(active_profile, dict) or canonical_json(active_profile) != canonical_json(
        profile
    ):
        raise RuntimeError("Funded solve used different model-profile metadata")
    expected_stationary_domain = (
        E5F_INCOME_ENTRY_DOMAIN
        if profile["name"] == REPAIRED_MODEL_PROFILE
        else E5F_DOMAIN
    )
    recorded_domain = metadata.get("active_profile_domain")
    if canonical_json(recorded_domain) != canonical_json(
        [list(row) for row in expected_stationary_domain]
    ):
        raise RuntimeError("Funded solve used a different stationary parameter domain")
    expected_income_states = int(profile.get("income_state_count", 5))
    if int(metadata.get("income_state_count", -1)) != expected_income_states:
        raise RuntimeError("Funded solve used a different income-state count")
    expected_cost = float(best["theta"].get("first_birth_fixed_cost", 0.0))
    solved_cost = require_finite(
        metadata.get("first_birth_fixed_cost"), "solved first-birth cost"
    )
    if not close(solved_cost, expected_cost, absolute=1.0e-12, relative=0.0):
        raise RuntimeError("Funded solve did not preserve the selected first-birth cost")
    if profile["name"] == REPAIRED_MODEL_PROFILE:
        if expected_income_states != 15:
            raise RuntimeError("Repaired profile does not declare 15 income states")
        for field in (
            "profile_id",
            "permanent_income_log_variance",
            "income_state_count",
            "first_birth_fixed_cost",
            "first_birth_fixed_cost_semantics",
        ):
            if field not in active_profile or canonical_json(
                active_profile[field]
            ) != canonical_json(profile[field]):
                raise RuntimeError(f"Funded repaired-profile field differs at {field}")


def parameter_table(
    best: dict[str, Any], domain: tuple[tuple[str, float, float, str], ...]
) -> list[dict[str, Any]]:
    theta = dict(best["theta"])
    rows: list[dict[str, Any]] = []
    for name, lower, upper, _ in domain:
        if name == "psi_child_2023":
            value = float(best["new_psi_child"])
        elif name == "psi_child_change_2023":
            value = float(best["new_psi_child"]) - float(best["old_psi_child"])
        elif name == "beta_annual":
            value = float(theta["beta"]) ** 0.25
        else:
            if name not in theta:
                raise RuntimeError(f"Selected theta omits free parameter {name}")
            value = float(theta[name])
        if not math.isfinite(value) or value < lower - 1.0e-12 or value > upper + 1.0e-12:
            raise RuntimeError(
                f"Selected {name}={value} lies outside its [{lower}, {upper}] domain"
            )
        span = upper - lower
        rows.append(
            {
                "parameter": name,
                "estimate": value,
                "lower_bound": lower,
                "upper_bound": upper,
                "near_bound": min(value - lower, upper - value) <= 0.02 * span,
                "status": "free_in_bounded_transition_search",
            }
        )
    old_psi = float(best["old_psi_child"])
    rows.append(
        {
            "parameter": "psi_child_2007",
            "estimate": old_psi,
            "lower_bound": -3.0,
            "upper_bound": 3.0,
            "near_bound": min(old_psi + 3.0, 3.0 - old_psi) <= 0.12,
            "status": "derived_to_match_old_completed_fertility",
        }
    )
    if str(best["model_profile"]["name"]) == REPAIRED_MODEL_PROFILE:
        rows.extend(
            [
                {
                    "parameter": "psi_child_2023",
                    "estimate": float(best["new_psi_child"]),
                    "lower_bound": None,
                    "upper_bound": None,
                    "near_bound": False,
                    "status": "derived_as_psi_child_2007_plus_preference_change",
                },
                {
                    "parameter": "permanent_income_log_variance",
                    "estimate": float(
                        best["model_profile"]["permanent_income_log_variance"]
                    ),
                    "lower_bound": None,
                    "upper_bound": None,
                    "near_bound": False,
                    "status": "externally_measured_not_estimated",
                },
            ]
        )
    rows.append(
        {
            "parameter": "outside_origin_entry_share",
            "estimate": float(best["outside_origin_entry_share"]),
            "lower_bound": None,
            "upper_bound": None,
            "near_bound": False,
            "status": "externally_fixed_provisional_anchor",
        }
    )
    return rows


def validate_selected_transition_path(
    rows: list[dict[str, str]], best: dict[str, Any]
) -> None:
    if len(rows) < 5:
        raise RuntimeError("Selected transition path does not reach the matched 2023 date")
    first_five = rows[:5]
    if [int(row["period"]) for row in first_five] != list(range(5)) or [
        int(float(row["years_from_start"])) for row in first_five
    ] != [0, 4, 8, 12, 16]:
        raise RuntimeError("Selected transition path has an invalid 2007--2023 calendar")
    required = (
        "population_index",
        "asset_price_index",
        "psi_child",
        "birth_children_topcode_adjusted",
    )
    if not all(
        math.isfinite(float(row[field])) for row in rows for field in required
    ):
        raise RuntimeError("Selected transition path contains nonfinite report objects")
    gates = (
        (first_five[0]["psi_child"], best["old_psi_child"]),
        (first_five[4]["psi_child"], best["new_psi_child"]),
        (first_five[4]["population_index"], best["terminal_population_index"]),
        (first_five[4]["asset_price_index"], best["terminal_housing_cost_index"]),
    )
    if any(
        not math.isclose(float(observed), float(expected), rel_tol=0.0, abs_tol=2.0e-8)
        for observed, expected in gates
    ):
        raise RuntimeError("Selected transition path disagrees with its reported endpoints")


def make_figures(
    output_dir: Path,
    comparison: list[dict[str, Any]],
    transition_rows: list[dict[str, str]],
    round_rows: list[dict[str, Any]],
    *,
    show_anchor: bool,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    short_labels = {
        "tfr": "Completed fertility",
        "childless_rate": "Childlessness",
        "mean_age_first_birth": "Mean age at first birth",
        "share_first_births_age30plus": "First births at 30+",
        "housing_increment_0to1": "Rooms after first birth",
        "prime30_55_parent_3plus_minus_1to2_mean_rooms": "Room gap: 3+ vs. 1–2 children",
        "own_family_gap": "Family ownership gap",
        "own_rate": "Ownership rate",
        "aggregate_mean_occupied_rooms_18_85": "Mean occupied rooms",
        "aggregate_wealth_to_annual_gross_labor_earnings": "Wealth / earnings",
        "annual_bequest_flow_to_aggregate_wealth": "Bequests / wealth",
        "old_total_wealth_to_annual_income_p90_p50_7684": "Old wealth p90/p50",
    }
    moments = [row["moment"] for row in comparison]
    positions = np.arange(len(moments))
    fig, ax = plt.subplots(figsize=(8.6, 6.2), constrained_layout=True)
    if show_anchor:
        ax.scatter(
            [row["anchor_standardized_gap"] for row in comparison],
            positions - 0.13,
            label="Transition anchor",
        )
    ax.scatter(
        [row["best_standardized_gap"] for row in comparison],
        positions + (0.13 if show_anchor else 0.0),
        label="Selected calibration",
    )
    ax.axvline(0.0, color="black", lw=0.8)
    ax.set_yticks(positions)
    ax.set_yticklabels([short_labels.get(moment, moment) for moment in moments])
    ax.invert_yaxis()
    ax.set_xlabel("Signed square root of loss contribution")
    ax.set_title("Dated 2023 target fit")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(output_dir / "target_fit_comparison.png", dpi=220)
    fig.savefig(output_dir / "target_fit_comparison.pdf")
    plt.close(fig)

    years = np.asarray([float(row["years_from_start"]) for row in transition_rows])
    population = np.asarray([float(row["population_index"]) for row in transition_rows])
    housing = np.asarray([float(row["asset_price_index"]) for row in transition_rows])
    psi = np.asarray([float(row["psi_child"]) for row in transition_rows])
    births = np.asarray(
        [float(row["birth_children_topcode_adjusted"]) for row in transition_rows]
    )
    births /= births[0]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    axes[0].plot(years, population, marker="o", label="Households")
    axes[0].plot(years, housing, marker="o", label="Housing cost")
    axes[0].axhline(1.0, color="black", lw=0.7)
    axes[0].set(xlabel="Years from 2007", ylabel="Index, 2007 = 1")
    axes[0].set_title("Households and housing costs both rise")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)
    birth_line = axes[1].plot(years, births, marker="o", label="Birth flow index")
    preference_axis = axes[1].twinx()
    preference_line = preference_axis.plot(
        years, psi, marker="o", color="tab:orange", label=r"Preference shifter $\psi_t$"
    )
    axes[1].set(xlabel="Years from 2007", ylabel="Birth flow, 2007 = 1")
    preference_axis.set_ylabel(r"Preference shifter $\psi_t$")
    axes[1].set_title("Preference decline lowers births")
    lines = birth_line + preference_line
    axes[1].legend(lines, [line.get_label() for line in lines], frameon=False)
    axes[1].grid(alpha=0.25)
    fig.savefig(output_dir / "best_transition_path.png", dpi=220)
    fig.savefig(output_dir / "best_transition_path.pdf")
    plt.close(fig)

    if round_rows:
        fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
        ax.plot(
            [row["round"] for row in round_rows],
            [row["best_loss"] for row in round_rows],
            marker="o",
        )
        ax.set(
            xlabel="Search round",
            ylabel="Best transition loss",
            xticks=[row["round"] for row in round_rows],
        )
        ax.grid(alpha=0.25)
        fig.savefig(output_dir / "search_progress.png", dpi=220)
        fig.savefig(output_dir / "search_progress.pdf")
        plt.close(fig)


def make_funded_policy_figure(
    output_dir: Path, rows: list[dict[str, Any]]
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    grouped: dict[str, list[dict[str, Any]]] = {case: [] for case in EXPECTED_FUNDED_CASES}
    for row in rows:
        case = str(row.get("case") or row.get("scenario") or "")
        if case in grouped:
            grouped[case].append(row)
    panels = (
        ("asset_price_index", "Housing-cost index", 1.0),
        ("population_index", "Adult-household population index", 1.0),
        ("topcode_adjusted_births_per_adult", "Births per adult", 1.0),
        ("dependent_child_owner_rate", "Dependent-child ownership (%)", 100.0),
    )
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 7.0), sharex=True)
    for axis, (field, title, scale) in zip(axes.flat, panels, strict=True):
        for case in EXPECTED_FUNDED_CASES:
            case_rows = sorted(grouped[case], key=lambda row: int(float(row["year"])))
            axis.plot(
                [int(float(row["year"])) for row in case_rows],
                [scale * float(row[field]) for row in case_rows],
                label=FUNDED_CASE_LABELS[case],
                color=FUNDED_CASE_COLORS[case],
                linestyle=FUNDED_CASE_LINESTYLES[case],
                linewidth=2.0,
            )
        axis.axvline(2027, color="#777777", linewidth=1.0, alpha=0.8)
        axis.set_title(title)
        axis.grid(alpha=0.22, linewidth=0.6)
        axis.tick_params(labelsize=9)
    for axis in axes[-1, :]:
        axis.set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
    )
    figure.suptitle("Transition paths; reforms begin in 2027", fontsize=12)
    figure.tight_layout(rect=(0.0, 0.08, 1.0, 0.96))
    figure.savefig(output_dir / "funded_policy_transition_paths.png", dpi=220)
    figure.savefig(output_dir / "funded_policy_transition_paths.pdf")
    plt.close(figure)


def parse_bool(value: Any, label: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise RuntimeError(f"Malformed Boolean {label}: {value!r}")


def assert_endpoint_row_matches(
    metadata_endpoint: dict[str, Any], disk_row: dict[str, str]
) -> None:
    for key, expected in metadata_endpoint.items():
        if key not in disk_row:
            raise RuntimeError(f"Stationary endpoint CSV omits metadata field {key}")
        observed = disk_row[key]
        if expected is None:
            if str(observed).strip():
                raise RuntimeError(f"Stationary endpoint differs at {key}")
        elif isinstance(expected, bool):
            if parse_bool(observed, key) != expected:
                raise RuntimeError(f"Stationary endpoint differs at {key}")
        elif isinstance(expected, (int, float)):
            observed_number = float(observed)
            expected_number = float(expected)
            if math.isnan(expected_number):
                equal = math.isnan(observed_number)
            else:
                equal = math.isclose(
                    observed_number,
                    expected_number,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-12,
                )
            if not equal:
                raise RuntimeError(f"Stationary endpoint differs at {key}")
        elif str(observed) != str(expected):
            raise RuntimeError(f"Stationary endpoint differs at {key}")


def assert_csv_row_is_payload_subset(
    payload: dict[str, Any], disk_row: dict[str, str], label: str
) -> None:
    for key, observed in disk_row.items():
        if key not in payload:
            raise RuntimeError(f"{label} contains unbound field {key}")
        expected = payload[key]
        if expected is None:
            equal = not str(observed).strip()
        elif isinstance(expected, bool):
            equal = parse_bool(observed, f"{label}.{key}") == expected
        elif isinstance(expected, (int, float)):
            observed_number = require_finite(observed, f"{label}.{key}")
            expected_number = require_finite(expected, f"{label}.{key}")
            equal = math.isclose(
                observed_number,
                expected_number,
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            )
        else:
            equal = str(observed) == str(expected)
        if not equal:
            raise RuntimeError(f"{label} differs at {key}")


def finite_csv_fields(
    row: dict[str, Any], fields: tuple[str, ...], label: str
) -> dict[str, float]:
    missing = [field for field in fields if field not in row]
    if missing:
        raise RuntimeError(f"{label} omits numeric fields: {missing}")
    return {field: require_finite(row[field], f"{label}.{field}") for field in fields}


def validate_stationary_endpoint_accounting(
    endpoint: dict[str, Any], metadata: dict[str, Any], best: dict[str, Any]
) -> None:
    common_fields = (
        "asset_price",
        "price_ratio",
        "housing_user_cost",
        "psi_child",
        "total_mass",
        "entry_households_per_adult",
        "birth_children_per_adult",
        "top_bin_entry_birth_flow_per_adult",
        "topcode_adjusted_birth_children_per_adult",
        "mature_entrant_households_per_adult",
        "topcode_adjusted_mature_entrant_households_per_adult",
        "reproduction_ratio_B_over_E",
        "raw_state_B_over_E",
        "reproduction_residual_B_minus_E",
        "topcode_adjusted_B_over_E_diagnostic",
        "topcode_consistent_B_over_E",
        "extra_children_per_top_bin_family",
        "threeplus_share_of_mature_flow",
        "maturation_exit_yield",
        "entrant_conversion_factor",
        "queue_birth_children_raw_explicit_states",
        "queue_birth_children_topcode_adjusted",
        "queue_mature_entrant_flow_B",
        "renewal_denominator",
        "stationary_population_scale",
        "housing_demand_per_adult",
        "housing_supply",
        "stationary_housing_supply",
        "fixed_stock_market_gap",
        "fixed_stock_relative_market_gap",
        "tfr_topcoded",
        "mean_completed_fertility_raw",
        "childless_rate",
        "own_rate",
    )
    values = {
        field: require_finite(endpoint.get(field), f"stationary endpoint.{field}")
        for field in common_fields
    }
    positive = (
        "asset_price",
        "price_ratio",
        "housing_user_cost",
        "total_mass",
        "entry_households_per_adult",
        "birth_children_per_adult",
        "topcode_adjusted_birth_children_per_adult",
        "mature_entrant_households_per_adult",
        "topcode_adjusted_mature_entrant_households_per_adult",
        "maturation_exit_yield",
        "entrant_conversion_factor",
        "renewal_denominator",
        "stationary_population_scale",
        "housing_demand_per_adult",
        "housing_supply",
        "stationary_housing_supply",
        "tfr_topcoded",
        "mean_completed_fertility_raw",
    )
    if any(values[field] <= 0.0 for field in positive):
        raise RuntimeError("Stationary endpoint has a nonpositive primitive or flow")
    if endpoint.get("status") != "complete" or not str(endpoint.get("label", "")):
        raise RuntimeError("Stationary endpoint has an incomplete status or label")
    if not str(endpoint.get("interpretation", "")):
        raise RuntimeError("Stationary endpoint omits its interpretation")
    if endpoint.get("policy_case") != "none":
        raise RuntimeError("Stationary endpoint is not the status-quo policy case")
    if endpoint.get("housing_supply_mode") != "static-elastic":
        raise RuntimeError("Stationary endpoint does not use static-elastic supply")
    if endpoint.get("renewal_child_accounting_method") != (
        "sequential_top_bin_birth_flow_common_maturation_yield"
    ):
        raise RuntimeError("Stationary endpoint uses the wrong child-accounting method")
    if not (
        0.0 <= values["childless_rate"] <= 1.0
        and 0.0 <= values["own_rate"] <= 1.0
    ):
        raise RuntimeError("Stationary endpoint has an invalid household rate")
    endpoint_contract = metadata.get("stationary_open_endpoint_contract")
    renewal = metadata.get("renewal")
    if not isinstance(endpoint_contract, dict) or not isinstance(renewal, dict):
        raise RuntimeError("Stationary endpoint omits renewal contracts")
    outside_flow = require_finite(
        endpoint_contract.get("outside_entry_flow"), "stationary outside flow"
    )
    retention = require_finite(
        endpoint_contract.get("retention"), "stationary retention"
    )
    conversion = require_finite(
        endpoint_contract.get("entrant_conversion_factor"),
        "stationary conversion",
    )
    maturation_yield = require_finite(
        endpoint_contract.get("topcode_consistent_maturation_survival_yield"),
        "stationary maturation yield",
    )
    if not (outside_flow > 0.0 and 0.0 <= retention <= 1.0 and conversion > 0.0):
        raise RuntimeError("Stationary renewal contract has invalid primitives")
    for observed, expected, label in (
        (renewal.get("outside_flow"), outside_flow, "renewal outside flow"),
        (renewal.get("retention"), retention, "renewal retention"),
        (renewal.get("entrant_conversion_factor"), conversion, "renewal conversion"),
        (values["entrant_conversion_factor"], conversion, "endpoint conversion"),
        (values["maturation_exit_yield"], maturation_yield, "endpoint maturation yield"),
        (values["psi_child"], best["new_psi_child"], "terminal preference"),
        (
            endpoint_contract.get("terminal_psi_child"),
            best["new_psi_child"],
            "endpoint-contract terminal preference",
        ),
    ):
        if not close(observed, expected, absolute=2.0e-12, relative=2.0e-12):
            raise RuntimeError(f"Stationary {label} identity failed")
    if endpoint_contract.get("housing_supply_mode") != "static-elastic":
        raise RuntimeError("Stationary endpoint contract uses the wrong supply mode")
    old_adjusted_mature = require_finite(
        renewal.get("topcode_adjusted_old_mature_flow"),
        "old adjusted mature flow",
    )
    old_raw_mature = require_finite(
        renewal.get("raw_old_mature_flow"), "old raw mature flow"
    )
    if min(old_adjusted_mature, old_raw_mature) <= 0.0:
        raise RuntimeError("Stationary renewal metadata has nonpositive old flows")
    old_normalization = metadata.get("old_fertility_normalization")
    if not isinstance(old_normalization, dict):
        raise RuntimeError("Funded metadata omits its old-fertility normalization")
    old_completed = require_finite(
        old_normalization.get("completed_fertility"),
        "old-normalization completed fertility",
    )
    old_target = require_finite(
        old_normalization.get("target"), "old-normalization target"
    )
    old_gap = require_finite(
        old_normalization.get("absolute_gap"), "old-normalization absolute gap"
    )
    old_solves = require_finite(
        old_normalization.get("stationary_solves"), "old-normalization solve count"
    )
    old_seconds = require_finite(
        old_normalization.get("stationary_solve_seconds"),
        "old-normalization solve seconds",
    )
    if (
        old_normalization.get("status") != "fixed_intercept"
        or old_solves < 1.0
        or old_solves != int(old_solves)
        or old_seconds < 0.0
    ):
        raise RuntimeError("Old-fertility normalization has an invalid solve contract")
    for observed, expected, label in (
        (old_normalization.get("psi_child"), best["old_psi_child"], "old preference"),
        (old_completed, best["old_completed_fertility"], "old fertility"),
        (old_target, best["old_completed_fertility"], "old fertility target"),
        (old_gap, abs(old_completed - old_target), "old fertility gap"),
    ):
        if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
            raise RuntimeError(f"Stationary {label} identity failed")

    raw_births = values["birth_children_per_adult"]
    adjusted_births = values["topcode_adjusted_birth_children_per_adult"]
    top_bin_births = values["top_bin_entry_birth_flow_per_adult"]
    extra = values["extra_children_per_top_bin_family"]
    raw_mature = values["mature_entrant_households_per_adult"]
    adjusted_mature = values[
        "topcode_adjusted_mature_entrant_households_per_adult"
    ]
    entry = values["entry_households_per_adult"]
    identities = (
        (adjusted_births, raw_births + extra * top_bin_births, "top-code births"),
        (values["queue_birth_children_raw_explicit_states"], raw_births, "raw queue births"),
        (
            values["queue_birth_children_topcode_adjusted"],
            adjusted_births,
            "adjusted queue births",
        ),
        (raw_mature, conversion * maturation_yield * raw_births, "raw mature flow"),
        (
            adjusted_mature,
            conversion * maturation_yield * adjusted_births,
            "adjusted mature flow",
        ),
        (values["queue_mature_entrant_flow_B"], adjusted_mature, "queue mature flow"),
        (values["reproduction_ratio_B_over_E"], raw_mature / entry, "raw B/E"),
        (values["raw_state_B_over_E"], raw_mature / entry, "raw-state B/E"),
        (
            values["reproduction_residual_B_minus_E"],
            raw_mature - entry,
            "raw renewal residual",
        ),
        (
            values["topcode_adjusted_B_over_E_diagnostic"],
            adjusted_mature / entry,
            "adjusted B/E diagnostic",
        ),
        (
            values["topcode_consistent_B_over_E"],
            adjusted_mature / entry,
            "top-code-consistent B/E",
        ),
        (
            values["threeplus_share_of_mature_flow"],
            top_bin_births / raw_births,
            "three-plus flow share",
        ),
    )
    for observed, expected, label in identities:
        if not close(observed, expected, absolute=3.0e-11, relative=3.0e-11):
            raise RuntimeError(f"Stationary {label} identity failed")
    denominator = entry - retention * adjusted_mature
    population = outside_flow / denominator
    supply = values["stationary_housing_supply"]
    demand = values["housing_demand_per_adult"]
    gap = population * demand - supply
    for observed, expected, label in (
        (values["renewal_denominator"], denominator, "renewal denominator"),
        (values["stationary_population_scale"], population, "population scale"),
        (values["fixed_stock_market_gap"], gap, "housing gap"),
        (values["fixed_stock_relative_market_gap"], gap / supply, "relative housing gap"),
    ):
        if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
            raise RuntimeError(f"Stationary {label} identity failed")
    if abs(values["fixed_stock_relative_market_gap"]) > 2.5e-5:
        raise RuntimeError("Stationary housing market does not clear")
    method = str(endpoint.get("market_clearing_method", ""))
    if method == "pure-price root":
        return
    if method != "persistent entry-type convexification at a discrete household-policy jump":
        raise RuntimeError(f"Unknown stationary market-clearing method: {method!r}")
    if not close(
        endpoint.get("scale_mapping_Hs_over_D"),
        supply / demand,
        absolute=3.0e-10,
    ) or not close(
        endpoint.get("stationary_total_housing_demand"),
        population * demand,
        absolute=3.0e-10,
    ):
        raise RuntimeError("Persistent endpoint housing-scale identity failed")
    persistent_fields = (
        "persistent_type_population_share_positive_branch",
        "persistent_type_population_share_negative_branch",
        "positive_branch_renewal_denominator",
        "negative_branch_renewal_denominator",
        "positive_branch_population",
        "negative_branch_population",
        "positive_branch_outside_entry_flow",
        "negative_branch_outside_entry_flow",
        "positive_branch_renewal_identity_residual",
        "negative_branch_renewal_identity_residual",
        "implied_outside_entrant_share_positive_branch",
        "implied_outside_entrant_share_negative_branch",
        "outside_entry_flow_identity_residual",
        "indifference_price_bracket_lower",
        "indifference_price_bracket_upper",
        "indifference_price_bracket_width",
        "pure_branch_positive_relative_market_gap",
        "pure_branch_negative_relative_market_gap",
        "existing_solver_pure_branch_asset_price",
        "existing_solver_pure_branch_relative_market_gap",
    )
    persistent = {
        field: require_finite(endpoint.get(field), f"persistent endpoint.{field}")
        for field in persistent_fields
    }
    w_pos = persistent["persistent_type_population_share_positive_branch"]
    w_neg = persistent["persistent_type_population_share_negative_branch"]
    pop_pos = persistent["positive_branch_population"]
    pop_neg = persistent["negative_branch_population"]
    den_pos = persistent["positive_branch_renewal_denominator"]
    den_neg = persistent["negative_branch_renewal_denominator"]
    flow_pos = persistent["positive_branch_outside_entry_flow"]
    flow_neg = persistent["negative_branch_outside_entry_flow"]
    share_pos = persistent["implied_outside_entrant_share_positive_branch"]
    share_neg = persistent["implied_outside_entrant_share_negative_branch"]
    if min(w_pos, w_neg, pop_pos, pop_neg, den_pos, den_neg, flow_pos, flow_neg) <= 0.0:
        raise RuntimeError("Persistent endpoint has a nonpositive type object")
    if max(w_pos, w_neg, share_pos, share_neg) >= 1.0:
        raise RuntimeError("Persistent endpoint has an invalid type or entrant share")
    persistent_identities = (
        (w_pos + w_neg, 1.0, "population shares"),
        (pop_pos + pop_neg, population, "branch populations"),
        (pop_pos / population, w_pos, "positive population share"),
        (pop_neg / population, w_neg, "negative population share"),
        (flow_pos, pop_pos * den_pos, "positive branch outside flow"),
        (flow_neg, pop_neg * den_neg, "negative branch outside flow"),
        (flow_pos + flow_neg, outside_flow, "aggregate outside flow"),
        (share_pos, flow_pos / outside_flow, "positive entrant share"),
        (share_neg, flow_neg / outside_flow, "negative entrant share"),
        (share_pos + share_neg, 1.0, "entrant shares"),
        (
            persistent["outside_entry_flow_identity_residual"],
            0.0,
            "outside-flow residual",
        ),
        (
            persistent["positive_branch_renewal_identity_residual"],
            0.0,
            "positive renewal residual",
        ),
        (
            persistent["negative_branch_renewal_identity_residual"],
            0.0,
            "negative renewal residual",
        ),
    )
    for observed, expected, label in persistent_identities:
        if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
            raise RuntimeError(f"Persistent endpoint {label} identity failed")
    lower = persistent["indifference_price_bracket_lower"]
    upper = persistent["indifference_price_bracket_upper"]
    width = persistent["indifference_price_bracket_width"]
    if not (
        lower > 0.0
        and upper > lower
        and width > 0.0
        and width <= 1.0e-6 * max(1.0, values["asset_price"])
        and close(width, upper - lower, absolute=1.0e-12)
        and close(values["asset_price"], 0.5 * (lower + upper), absolute=1.0e-12)
    ):
        raise RuntimeError("Persistent endpoint has an invalid price bracket")
    if persistent["pure_branch_positive_relative_market_gap"] <= 0.0 or persistent[
        "pure_branch_negative_relative_market_gap"
    ] >= 0.0:
        raise RuntimeError("Persistent endpoint pure branches do not bracket zero")
    expected_condition = (
        "type assigned at entry and never changes; locally renewed entrants retain "
        "the source type; outside entrants are split by the reported implied shares"
    )
    if endpoint.get("persistent_type_stationarity_condition") != expected_condition:
        raise RuntimeError("Persistent endpoint omits its type-stationarity condition")


def validate_stationary_endpoint_diagnostics(
    directory: Path, endpoint: dict[str, Any]
) -> None:
    endpoint_dir = directory / "stationary_endpoint_status_quo"
    pure_rows = read_csv(endpoint_dir / "stationary_pure_branch_endpoint.csv")
    search_rows = read_csv(endpoint_dir / "stationary_endpoint_search.csv")
    if len(pure_rows) != 1 or not search_rows:
        raise RuntimeError("Stationary endpoint diagnostics are empty or duplicated")
    pure = pure_rows[0]
    pure_values = finite_csv_fields(
        pure,
        (
            "asset_price",
            "psi_child",
            "entry_households_per_adult",
            "topcode_adjusted_birth_children_per_adult",
            "renewal_denominator",
            "stationary_population_scale",
            "housing_demand_per_adult",
            "stationary_housing_supply",
            "fixed_stock_market_gap",
            "fixed_stock_relative_market_gap",
        ),
        "stationary pure-branch diagnostic",
    )
    if pure.get("status") != "complete" or pure.get("policy_case") != "none":
        raise RuntimeError("Stationary pure-branch diagnostic has the wrong contract")
    finite_search: list[tuple[float, float]] = []
    for row in search_rows:
        try:
            price = float(row["asset_price"])
            gap = float(row["fixed_stock_relative_market_gap"])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(price) and math.isfinite(gap):
            finite_search.append((price, gap))
    if not finite_search:
        raise RuntimeError("Stationary endpoint search has no finite audited row")
    method = str(endpoint["market_clearing_method"])
    if method == "pure-price root":
        assert_csv_row_is_payload_subset(
            endpoint, pure, "stationary pure-price endpoint"
        )
        if not any(
            close(price, endpoint["asset_price"], absolute=1.0e-12)
            and close(
                gap,
                endpoint["fixed_stock_relative_market_gap"],
                absolute=1.0e-12,
            )
            for price, gap in finite_search
        ):
            raise RuntimeError("Pure stationary endpoint is absent from its search trace")
        return
    if method != "persistent entry-type convexification at a discrete household-policy jump":
        raise RuntimeError("Stationary endpoint diagnostic has an unknown clearing method")
    for observed, expected, label in (
        (
            endpoint["existing_solver_pure_branch_asset_price"],
            pure_values["asset_price"],
            "existing pure-branch price",
        ),
        (
            endpoint["existing_solver_pure_branch_relative_market_gap"],
            pure_values["fixed_stock_relative_market_gap"],
            "existing pure-branch gap",
        ),
    ):
        if not close(observed, expected, absolute=2.0e-12, relative=2.0e-12):
            raise RuntimeError(f"Persistent endpoint {label} differs from its audit")
    lower = require_finite(
        endpoint["indifference_price_bracket_lower"], "persistent lower price"
    )
    upper = require_finite(
        endpoint["indifference_price_bracket_upper"], "persistent upper price"
    )
    lower_rows = [gap for price, gap in finite_search if close(price, lower)]
    upper_rows = [gap for price, gap in finite_search if close(price, upper)]
    if len(lower_rows) != 1 or len(upper_rows) != 1:
        raise RuntimeError("Persistent endpoint bracket is absent from its search trace")
    bracket_gaps = (lower_rows[0], upper_rows[0])
    positive = max(bracket_gaps)
    negative = min(bracket_gaps)
    if positive <= 0.0 or negative >= 0.0:
        raise RuntimeError("Persistent endpoint search rows do not bracket zero")
    if not close(
        positive,
        endpoint["pure_branch_positive_relative_market_gap"],
        absolute=2.0e-12,
        relative=2.0e-12,
    ) or not close(
        negative,
        endpoint["pure_branch_negative_relative_market_gap"],
        absolute=2.0e-12,
        relative=2.0e-12,
    ):
        raise RuntimeError("Persistent endpoint branch gaps differ from its search trace")


def validate_funded_paths_and_summaries(
    raw_paths: list[dict[str, str]],
    summaries: list[dict[str, str]],
    metadata: dict[str, Any],
    directory: Path,
    best: dict[str, Any],
) -> list[dict[str, Any]]:
    if [row.get("case") for row in summaries] != list(EXPECTED_FUNDED_CASES):
        raise RuntimeError("Funded summary does not preserve the ordered four-case contract")
    expected_long_order = [
        case
        for case in EXPECTED_FUNDED_CASES
        for _ in range(EXPECTED_FUNDED_PERIODS)
    ]
    observed_long_order = [
        str(row.get("case") or row.get("scenario") or "") for row in raw_paths
    ]
    if observed_long_order != expected_long_order:
        raise RuntimeError("Funded long path is not the ordered four-by-45 packet")
    numerical = metadata["numerical_gates"]
    renewal = metadata.get("renewal")
    if not isinstance(renewal, dict):
        raise RuntimeError("Funded metadata omits its renewal accounting")
    conversion = require_finite(
        renewal.get("entrant_conversion_factor"), "funded entrant conversion"
    )
    retention = require_finite(renewal.get("retention"), "funded retention")
    outside_flow = require_finite(renewal.get("outside_flow"), "funded outside flow")
    endpoint_contract = metadata.get("stationary_open_endpoint_contract")
    if not isinstance(endpoint_contract, dict):
        raise RuntimeError("Funded metadata omits its stationary endpoint contract")
    maturation_yield = require_finite(
        endpoint_contract.get("topcode_consistent_maturation_survival_yield"),
        "funded maturation yield",
    )
    if (
        conversion <= 0.0
        or not 0.0 <= retention <= 1.0
        or outside_flow <= 0.0
        or not 0.0 < maturation_yield <= 1.0 + 1.0e-10
    ):
        raise RuntimeError("Funded renewal accounting has invalid primitives")
    market_tolerance = require_finite(
        numerical["market_tolerance"], "market tolerance"
    )
    fiscal_tolerance = require_finite(
        numerical["fiscal_absolute_tolerance"], "fiscal tolerance"
    )
    mass_tolerance = require_finite(numerical["mass_tolerance"], "mass tolerance")
    if not (
        close(market_tolerance, 2.0e-4, absolute=1.0e-15, relative=0.0)
        and close(fiscal_tolerance, 2.5e-5, absolute=1.0e-15, relative=0.0)
        and close(mass_tolerance, 2.0e-8, absolute=1.0e-15, relative=0.0)
    ):
        raise RuntimeError("Funded packet uses noncanonical numerical tolerances")

    paths: dict[str, list[dict[str, Any]]] = {
        case: [] for case in EXPECTED_FUNDED_CASES
    }
    canonical_paths: list[dict[str, Any]] = []
    for raw in raw_paths:
        case = str(raw.get("case") or raw.get("scenario") or "")
        if case not in EXPECTED_FUNDED_CASES:
            raise RuntimeError(f"Funded long path has an unknown case: {case!r}")
        if raw.get("scenario") != case or raw.get("policy_case") != case:
            raise RuntimeError(f"{case}: scenario or policy-case label differs")
        if raw.get("closure") != "open_birth_vintage":
            raise RuntimeError(f"{case}: path uses the wrong population closure")
        if raw.get("housing_supply_mode") != "static-elastic":
            raise RuntimeError(f"{case}: path uses the wrong housing-supply convention")
        row: dict[str, Any] = {**raw, "case": case}
        for field, item in tuple(row.items()):
            token = str(item).strip().lower()
            if token not in {"nan", "+nan", "-nan", "inf", "+inf", "-inf"}:
                continue
            target_missing = not str(row.get("population_target_index", "")).strip()
            if field in OPTIONAL_PATH_NONFINITE_AS_MISSING and target_missing:
                row[field] = ""
                continue
            raise RuntimeError(f"{case}: nonfinite path field {field}")
        for field in FUNDED_PATH_JSON_FINITE_FIELDS:
            try:
                sequence = json.loads(str(row[field]))
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
                raise RuntimeError(f"{case}: malformed path sequence {field}") from error
            if not isinstance(sequence, list) or len(sequence) != 4:
                raise RuntimeError(f"{case}: path sequence {field} is not length four")
            assert_finite_tree(sequence, f"{case} path.{field}")
        values = finite_csv_fields(row, FUNDED_PATH_FINITE_FIELDS, f"{case} path")
        period = int(values["period"])
        if values["period"] != period or period != len(paths[case]):
            raise RuntimeError(f"{case}: periods are not the ordered integers 0--44")
        if values["year"] != 2007 + 4 * period or values[
            "years_from_start"
        ] != 4 * period:
            raise RuntimeError(f"{case}: invalid calendar at period {period}")
        expected_psi = float(best["old_psi_child"]) + min(period / 4.0, 1.0) * (
            float(best["new_psi_child"]) - float(best["old_psi_child"])
        )
        if not close(values["psi_child"], expected_psi, absolute=2.0e-12):
            raise RuntimeError(f"{case}: incorrect preference path at period {period}")
        active = parse_bool(row["policy_active"], f"{case} policy_active")
        expected_tax, expected_grant, funded_reform = FUNDED_CASE_CONVENTIONS[case]
        expected_active = bool(funded_reform and period >= 5)
        if active != expected_active:
            raise RuntimeError(f"{case}: incorrect policy activation at period {period}")
        expected_convention = (
            "balanced_budget_equal_rebate"
            if expected_active
            else "status_quo_unrebated_property_tax"
        )
        if row.get("fiscal_convention") != expected_convention:
            raise RuntimeError(f"{case}: incorrect fiscal convention at period {period}")
        active_tax = expected_tax if expected_active else 0.01
        if not close(values["annual_property_tax_rate"], active_tax, absolute=1.0e-14):
            raise RuntimeError(f"{case}: incorrect property-tax rate at period {period}")
        grant = parse_bool(row["purchase_grant"], f"{case} purchase_grant")
        if grant != bool(expected_active and expected_grant):
            raise RuntimeError(f"{case}: incorrect grant status at period {period}")
        fiscal_evaluations = int(values["fiscal_root_evaluations"])
        if values["fiscal_root_evaluations"] != fiscal_evaluations or (
            expected_active and fiscal_evaluations <= 0
        ) or (not expected_active and fiscal_evaluations != 0):
            raise RuntimeError(f"{case}: invalid fiscal solve count at period {period}")
        if expected_active:
            iterations = require_finite(
                row.get("joint_root_iterations"),
                f"{case} period {period} joint iterations",
            )
            if iterations < 1.0 or iterations != int(iterations):
                raise RuntimeError(f"{case}: invalid joint-root iterations")
        first = row if not paths[case] else paths[case][0]
        initial_population = float(first["adult_population"])
        initial_price = float(first["asset_price"])
        positive_fields = (
            "asset_price",
            "housing_user_cost",
            "adult_population",
            "entry_flow_E",
            "housing_demand",
            "housing_demand_per_adult",
            "housing_supply",
        )
        if any(values[field] <= 0.0 for field in positive_fields):
            raise RuntimeError(f"{case}: path has a nonpositive price, flow, or quantity")
        nonnegative_fields = (
            "birth_children",
            "birth_children_raw_explicit_states",
            "top_bin_entry_birth_flow",
            "birth_children_topcode_adjusted",
            "mature_children_raw",
            "effective_mature_entrant_flow_B",
            "raw_state_scheduled_mature_entrant_flow_B",
            "model_state_same_period_mature_flow_B",
            "birth_queue_new_potential_flow_B",
            "birth_queue_new_raw_state_potential_flow_B",
            "outside_entry_flow_M",
            "retained_mature_entrants",
            "entrant_flow_next",
            "youngest_age_mass_next",
            "adult_deaths",
            "relative_market_residual",
            "feasibility_frontier_projection_mass",
        )
        if any(values[field] < 0.0 for field in nonnegative_fields):
            raise RuntimeError(f"{case}: path has a negative demographic flow")
        rate_fields = (
            "owner_rate",
            "dependent_child_owner_rate",
            "parent_owner_rate",
            "childless_owner_rate",
        )
        if any(not 0.0 <= values[field] <= 1.0 for field in rate_fields):
            raise RuntimeError(f"{case}: path has an invalid ownership rate")
        path_identities = (
            (values["asset_price_index"], values["asset_price"] / initial_price, "price index"),
            (
                values["population_index"],
                values["adult_population"] / initial_population,
                "population index",
            ),
            (
                values["birth_children_raw_explicit_states"],
                values["birth_children"],
                "raw birth flow",
            ),
            (
                values["births_per_adult"],
                values["birth_children"] / values["adult_population"],
                "raw birth rate",
            ),
            (
                values["topcode_adjusted_births_per_adult"],
                values["birth_children_topcode_adjusted"] / values["adult_population"],
                "adjusted birth rate",
            ),
            (
                values["births_over_entry"],
                values["birth_children"] / values["entry_flow_E"],
                "raw births over entry",
            ),
            (
                values["topcode_adjusted_births_over_entry"],
                values["birth_children_topcode_adjusted"] / values["entry_flow_E"],
                "adjusted births over entry",
            ),
            (
                values["housing_demand_per_adult"],
                values["housing_demand"] / values["adult_population"],
                "housing demand per adult",
            ),
            (
                values["outside_entry_flow_share_of_initial_mass"],
                values["outside_entry_flow_M"] / initial_population,
                "outside-entry share",
            ),
            (
                values["entrant_flow_next"],
                values["retained_mature_entrants"] + values["outside_entry_flow_M"],
                "next entrant flow",
            ),
            (
                values["outside_entry_flow_M"],
                outside_flow,
                "outside-entry flow",
            ),
            (
                values["birth_queue_new_potential_flow_B"],
                conversion
                * maturation_yield
                * values["birth_children_topcode_adjusted"],
                "adjusted new birth-vintage flow",
            ),
            (
                values["birth_queue_new_raw_state_potential_flow_B"],
                conversion * maturation_yield * values["birth_children"],
                "raw new birth-vintage flow",
            ),
            (
                values["raw_state_scheduled_mature_entrant_flow_B"],
                conversion * values["mature_children_raw"],
                "scheduled raw mature flow",
            ),
            (
                values["retained_mature_entrants"],
                retention * values["effective_mature_entrant_flow_B"],
                "retained mature entrants",
            ),
            (
                values["dual_clock_raw_flow_gap_percent"],
                100.0
                * (
                    values["raw_state_scheduled_mature_entrant_flow_B"]
                    - values["model_state_same_period_mature_flow_B"]
                )
                / max(values["model_state_same_period_mature_flow_B"], 1.0e-15),
                "dual-clock raw flow gap",
            ),
        )
        for observed, expected, label in path_identities:
            if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
                raise RuntimeError(f"{case}: {label} identity failed at period {period}")
        if (
            values["birth_children_topcode_adjusted"]
            < values["birth_children"] - 3.0e-12
            or values["top_bin_entry_birth_flow"]
            > values["birth_children"] + 3.0e-12
        ):
            raise RuntimeError(f"{case}: top-code birth accounting is infeasible")
        reconstructed_market = abs(
            values["housing_demand"] - values["housing_supply"]
        ) / values["housing_supply"]
        if not close(
            values["relative_market_residual"],
            reconstructed_market,
            absolute=3.0e-10,
            relative=3.0e-10,
        ) or values["relative_market_residual"] > market_tolerance:
            raise RuntimeError(f"{case}: per-date housing-market gate failed")
        if abs(values["mass_accounting_residual"]) > mass_tolerance:
            raise RuntimeError(f"{case}: per-date mass gate failed")
        if values["nonfinite_distribution_count"] != 0.0:
            raise RuntimeError(f"{case}: distribution contains nonfinite cells")
        if values["min_distribution_mass"] < -1.0e-14:
            raise RuntimeError(f"{case}: negative distribution mass")
        if values["feasibility_frontier_projection_mass"] > 1.0e-6:
            raise RuntimeError(f"{case}: feasibility-projection gate failed")
        revenue = values["property_tax_revenue"]
        grant_outlays = values["purchase_grant_outlays"]
        transfer_outlays = values["lump_sum_transfer_outlays"]
        fiscal_residual = revenue - grant_outlays - transfer_outlays
        fiscal_scale = max(abs(revenue), abs(grant_outlays) + abs(transfer_outlays), 1.0e-12)
        implied_transfer = (revenue - grant_outlays) / values["adult_population"]
        fiscal_identities = (
            (
                transfer_outlays,
                values["lump_sum_transfer_period_units"] * values["adult_population"],
                "transfer outlays",
            ),
            (values["government_budget_residual"], fiscal_residual, "budget residual"),
            (
                values["scaled_government_budget_residual"],
                fiscal_residual / fiscal_scale,
                "scaled budget residual",
            ),
            (
                values["balanced_budget_implied_transfer"],
                implied_transfer,
                "implied transfer",
            ),
            (
                values["balanced_budget_transfer_gap"],
                values["lump_sum_transfer_period_units"] - implied_transfer,
                "transfer gap",
            ),
        )
        for observed, expected, label in fiscal_identities:
            if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
                raise RuntimeError(f"{case}: {label} identity failed at period {period}")
        fiscal_nonnegative = (
            revenue,
            grant_outlays,
            transfer_outlays,
            values["purchase_grant_recipient_mass"],
            values["lump_sum_transfer_period_units"],
        )
        if min(fiscal_nonnegative) < 0.0:
            raise RuntimeError(f"{case}: negative fiscal flow at period {period}")
        expected_grant_outlays = (
            0.4 * values["purchase_grant_recipient_mass"] if grant else 0.0
        )
        if not close(
            grant_outlays,
            expected_grant_outlays,
            absolute=3.0e-10,
            relative=3.0e-10,
        ):
            raise RuntimeError(f"{case}: purchase-grant accounting failed")
        if expected_active and (
            abs(fiscal_residual) > fiscal_tolerance
            or abs(values["balanced_budget_transfer_gap"])
            * values["adult_population"]
            > fiscal_tolerance
        ):
            raise RuntimeError(f"{case}: per-date funded fiscal gate failed")
        paths[case].append(row)
        canonical_paths.append(row)

    summary_by_case = {str(row["case"]): row for row in summaries}
    baseline_endpoint = paths[EXPECTED_FUNDED_CASES[0]][-1]
    summary_numeric_fields = (
        "policy_start_period",
        "policy_start_year",
        "annual_property_tax_rate",
        "periods",
        "endpoint_year",
        "endpoint_asset_price_index",
        "endpoint_population_index",
        "endpoint_births_per_adult",
        "endpoint_owner_rate",
        "endpoint_dependent_child_owner_rate",
        "maximum_market_residual",
        "maximum_mass_residual",
        "maximum_fiscal_residual",
        "maximum_scaled_fiscal_residual",
        "nonfinite_reported_objects",
        "maximum_nonfinite_distribution_count",
        "minimum_distribution_mass",
        "maximum_feasibility_frontier_projection_mass",
        "model_solve_count",
        "endpoint_asset_price_change_percent",
        "endpoint_population_change_percent",
        "endpoint_birth_rate_change_percent",
        "endpoint_owner_rate_change_pp",
        "endpoint_dependent_child_owner_rate_change_pp",
        "baseline_endpoint_year",
    )
    for case in EXPECTED_FUNDED_CASES:
        rows = paths[case]
        if len(rows) != EXPECTED_FUNDED_PERIODS or [int(row["period"]) for row in rows] != list(
            range(EXPECTED_FUNDED_PERIODS)
        ):
            raise RuntimeError(f"{case}: path is not exactly 45 ordered dates")
        case_file_rows = read_csv(directory / "cases" / case / "transition_path.csv")
        long_case_rows = [
            {key: value for key, value in row.items() if key != "case"}
            for row in rows
        ]
        rows_match = len(case_file_rows) == len(long_case_rows)
        if rows_match:
            for case_row, long_row in zip(case_file_rows, long_case_rows):
                fields = set(case_row) | set(long_row)
                if any(
                    (
                        ""
                        if field in OPTIONAL_PATH_NONFINITE_AS_MISSING
                        and str(case_row.get(field, "")).strip().lower()
                        in {"nan", "+nan", "-nan"}
                        else str(case_row.get(field, ""))
                    )
                    != str(long_row.get(field, ""))
                    for field in fields
                ):
                    rows_match = False
                    break
        if not rows_match:
            raise RuntimeError(f"{case}: long and case-specific path files differ")
        summary = summary_by_case[case]
        if summary.get("status") != "complete":
            raise RuntimeError(f"{case}: summary is incomplete")
        sv = finite_csv_fields(summary, summary_numeric_fields, f"{case} summary")
        expected_tax, expected_grant, expected_funded = FUNDED_CASE_CONVENTIONS[case]
        if parse_bool(summary["purchase_grant"], f"{case} summary grant") != expected_grant:
            raise RuntimeError(f"{case}: summary grant convention differs")
        if parse_bool(summary["funded_reform"], f"{case} funded flag") != expected_funded:
            raise RuntimeError(f"{case}: summary funded convention differs")
        if not close(sv["annual_property_tax_rate"], expected_tax, absolute=1.0e-14):
            raise RuntimeError(f"{case}: summary property-tax convention differs")
        if not (
            sv["policy_start_period"] == 5.0
            and sv["policy_start_year"] == 2027.0
            and sv["periods"] == EXPECTED_FUNDED_PERIODS
            and sv["endpoint_year"] == EXPECTED_FUNDED_ENDPOINT_YEAR
            and sv["baseline_endpoint_year"] == EXPECTED_FUNDED_ENDPOINT_YEAR
            and sv["model_solve_count"] > 0.0
            and sv["model_solve_count"] == int(sv["model_solve_count"])
            and sv["nonfinite_reported_objects"] == 0.0
            and sv["maximum_nonfinite_distribution_count"] == 0.0
        ):
            raise RuntimeError(f"{case}: summary timing or solve-count contract differs")
        endpoint = rows[-1]
        endpoint_identities = (
            (sv["endpoint_asset_price_index"], endpoint["asset_price_index"], "price"),
            (sv["endpoint_population_index"], endpoint["population_index"], "population"),
            (
                sv["endpoint_births_per_adult"],
                endpoint["topcode_adjusted_births_per_adult"],
                "birth rate",
            ),
            (sv["endpoint_owner_rate"], endpoint["owner_rate"], "ownership"),
            (
                sv["endpoint_dependent_child_owner_rate"],
                endpoint["dependent_child_owner_rate"],
                "dependent-child ownership",
            ),
            (
                sv["maximum_market_residual"],
                max(float(row["relative_market_residual"]) for row in rows),
                "maximum market residual",
            ),
            (
                sv["maximum_mass_residual"],
                max(abs(float(row["mass_accounting_residual"])) for row in rows),
                "maximum mass residual",
            ),
            (
                sv["maximum_fiscal_residual"],
                max(
                    (
                        abs(float(row["government_budget_residual"]))
                        for row in rows
                        if parse_bool(row["policy_active"], f"{case} active")
                    ),
                    default=0.0,
                ),
                "maximum fiscal residual",
            ),
            (
                sv["maximum_scaled_fiscal_residual"],
                max(
                    (
                        abs(float(row["scaled_government_budget_residual"]))
                        for row in rows
                        if parse_bool(row["policy_active"], f"{case} active")
                    ),
                    default=0.0,
                ),
                "maximum scaled fiscal residual",
            ),
            (
                sv["minimum_distribution_mass"],
                min(float(row["min_distribution_mass"]) for row in rows),
                "minimum mass",
            ),
            (
                sv["maximum_feasibility_frontier_projection_mass"],
                max(float(row["feasibility_frontier_projection_mass"]) for row in rows),
                "maximum feasibility projection",
            ),
        )
        for observed, expected, label in endpoint_identities:
            if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
                raise RuntimeError(f"{case}: summary {label} identity failed")
        baseline_values = {
            "price": float(baseline_endpoint["asset_price_index"]),
            "population": float(baseline_endpoint["population_index"]),
            "birth_rate": float(baseline_endpoint["topcode_adjusted_births_per_adult"]),
            "owner": float(baseline_endpoint["owner_rate"]),
            "dependent_owner": float(baseline_endpoint["dependent_child_owner_rate"]),
        }
        effects = (
            (
                sv["endpoint_asset_price_change_percent"],
                100.0 * (float(endpoint["asset_price_index"]) / baseline_values["price"] - 1.0),
                "price effect",
            ),
            (
                sv["endpoint_population_change_percent"],
                100.0
                * (float(endpoint["population_index"]) / baseline_values["population"] - 1.0),
                "population effect",
            ),
            (
                sv["endpoint_birth_rate_change_percent"],
                100.0
                * (
                    float(endpoint["topcode_adjusted_births_per_adult"])
                    / baseline_values["birth_rate"]
                    - 1.0
                ),
                "birth-rate effect",
            ),
            (
                sv["endpoint_owner_rate_change_pp"],
                100.0 * (float(endpoint["owner_rate"]) - baseline_values["owner"]),
                "ownership effect",
            ),
            (
                sv["endpoint_dependent_child_owner_rate_change_pp"],
                100.0
                * (
                    float(endpoint["dependent_child_owner_rate"])
                    - baseline_values["dependent_owner"]
                ),
                "dependent-child ownership effect",
            ),
        )
        for observed, expected, label in effects:
            if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
                raise RuntimeError(f"{case}: summary {label} identity failed")
    aggregate_identities = (
        (
            numerical["maximum_market_residual"],
            max(float(row["maximum_market_residual"]) for row in summaries),
            "metadata maximum market residual",
        ),
        (
            numerical["maximum_fiscal_residual"],
            max(float(row["maximum_fiscal_residual"]) for row in summaries),
            "metadata maximum fiscal residual",
        ),
        (
            numerical["maximum_mass_residual"],
            max(float(row["maximum_mass_residual"]) for row in summaries),
            "metadata maximum mass residual",
        ),
        (
            numerical["minimum_distribution_mass"],
            min(float(row["minimum_distribution_mass"]) for row in summaries),
            "metadata minimum mass",
        ),
        (
            numerical["maximum_feasibility_frontier_projection_mass"],
            max(
                float(row["maximum_feasibility_frontier_projection_mass"])
                for row in summaries
            ),
            "metadata feasibility projection",
        ),
    )
    for observed, expected, label in aggregate_identities:
        if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
            raise RuntimeError(f"Funded {label} identity failed")
    if int(numerical["maximum_nonfinite_reported_objects"]) != 0 or int(
        numerical["maximum_nonfinite_distribution_count"]
    ) != 0:
        raise RuntimeError("Funded metadata reports nonfinite model objects")
    return canonical_paths


def live_funded_code_bundle(relative_files: list[str]) -> str:
    if not relative_files or len(set(relative_files)) != len(relative_files):
        raise RuntimeError("Funded code-bundle file list is empty or duplicated")
    if tuple(relative_files) != EXPECTED_FUNDED_CODE_BUNDLE_FILES:
        raise RuntimeError("Funded code bundle differs from the canonical ordered file list")
    digest = hashlib.sha256()
    for relative in relative_files:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(f"Funded code-bundle file is unavailable: {path}")
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def load_funded_policy_packet(
    directory: Path,
    *,
    selected_report: Path,
    selected_contract: dict[str, Any],
    best: dict[str, Any],
) -> dict[str, Any]:
    metadata = read_json(directory / "metadata.json")
    launch = read_json(directory / "launch_contract.json")
    launcher_status = read_json(directory / "launcher_status.json")
    latest_case = read_json(directory / "latest_completed_case.json")
    summaries = read_csv(directory / "summary.csv")
    raw_paths = read_csv(directory / "transition_path_long.csv")
    endpoint_rows = read_csv(
        directory
        / "stationary_endpoint_status_quo"
        / "stationary_open_endpoint.csv"
    )
    assert_finite_tree(metadata, "funded metadata")
    assert_finite_tree(launch, "funded launch contract")
    assert_finite_tree(launcher_status, "funded launcher status")
    assert_finite_tree(latest_case, "funded completed-case ledger")
    if metadata.get("status") != "complete" or bool(metadata.get("smoke", True)):
        raise RuntimeError("Funded-policy packet is incomplete or marked as a smoke run")
    if launcher_status.get("status") != "complete" or not str(
        launcher_status.get("completed_utc", "")
    ):
        raise RuntimeError("Funded launcher did not certify a complete final packet")
    if launch.get("status") != "launched" or not str(launch.get("launched_utc", "")):
        raise RuntimeError("Funded launch contract is missing its launch certification")
    if (
        latest_case.get("status") != "complete"
        or int(latest_case.get("completed_cases", -1)) != len(EXPECTED_FUNDED_CASES)
        or int(latest_case.get("total_cases", -1)) != len(EXPECTED_FUNDED_CASES)
    ):
        raise RuntimeError("Funded completed-case ledger is not exactly four of four")
    launch_job = str(launch.get("slurm_job_id", ""))
    if not launch_job or str(launcher_status.get("slurm_job_id")) != launch_job:
        raise RuntimeError("Launcher completion and launch contracts use different jobs")
    if (
        tuple(str(case) for case in launch.get("policy_cases", []))
        != EXPECTED_FUNDED_CASES
        or int(launch.get("post_2023_periods", -1)) != 40
        or int(launch.get("total_calendar_points", -1))
        != len(EXPECTED_FUNDED_CASES) * EXPECTED_FUNDED_PERIODS
    ):
        raise RuntimeError("Launch contract is not the exact four-by-45 production job")
    if int(metadata.get("periods", -1)) != EXPECTED_FUNDED_PERIODS or int(
        metadata.get("endpoint_year", -1)
    ) != EXPECTED_FUNDED_ENDPOINT_YEAR:
        raise RuntimeError("Funded metadata is not the exact 45-date production packet")
    cases = tuple(str(case) for case in metadata.get("cases", []))
    if cases != EXPECTED_FUNDED_CASES:
        raise RuntimeError(f"Funded packet does not contain the four canonical cases: {cases}")
    exact_metadata_strings = {
        "policy_start": "first model date after the matched 2023 economy (2027)",
        "relative_effect_baseline": EXPECTED_FUNDED_CASES[0],
        "prepolicy_fiscal_convention": (
            "calibrated unrebated 1 percent property tax through 2023"
        ),
        "policy_fiscal_rule": (
            "property-tax revenue equals equal per-household transfers plus "
            "targeted purchase-grant outlays at every policy date"
        ),
        "purchase_grant": (
            "0.4 model units for eligible dependent-child renters purchasing "
            "an owner rung with at least six rooms"
        ),
        "housing_price_expectations": "temporary equilibrium/static expectations",
    }
    for field, expected in exact_metadata_strings.items():
        if metadata.get(field) != expected:
            raise RuntimeError(f"Funded metadata convention differs at {field}")
    report_hash = sha256(selected_report)
    funded_contract = metadata.get("contract")
    if not isinstance(funded_contract, dict):
        raise RuntimeError("Funded metadata omits its transition contract")
    direct_contract_checks = {
        "report_sha256": report_hash,
        "source_sha256": selected_contract["source_sha256"],
        "target_set": selected_contract["target_set"],
        "target_fingerprint": selected_contract["target_fingerprint"],
    }
    for field, expected in direct_contract_checks.items():
        if str(funded_contract.get(field)) != str(expected):
            raise RuntimeError(f"Funded transition contract differs at {field}")
    if canonical_json(funded_contract.get("theta")) != canonical_json(best["theta"]):
        raise RuntimeError("Funded transition contract uses a different theta")
    for field, expected in (
        ("transition_loss", best["transition_loss"]),
        ("old_psi_child", best["old_psi_child"]),
        ("new_psi_child", best["new_psi_child"]),
        ("old_completed_fertility", best["old_completed_fertility"]),
        ("matched_2023_population_index", best["terminal_population_index"]),
        ("matched_2023_housing_cost_index", best["terminal_housing_cost_index"]),
    ):
        if not close(
            funded_contract.get(field), expected, absolute=1.0e-12, relative=0.0
        ):
            raise RuntimeError(f"Funded transition contract differs at {field}")
    funded_profile = funded_contract.get("model_profile")
    if not isinstance(funded_profile, dict):
        raise RuntimeError("Funded transition contract omits model_profile")
    if str(funded_profile.get("name")) != str(selected_contract["model_profile"]["name"]):
        raise RuntimeError("Funded transition contract uses a different model profile")
    if bool(funded_profile.get("legacy_default", True)):
        raise RuntimeError("Funded transition contract used a legacy profile exception")
    if canonical_json(funded_profile.get("report_metadata")) != canonical_json(
        selected_contract["model_profile"]
    ):
        raise RuntimeError("Funded model-profile metadata differs from calibration")
    if canonical_json(funded_contract.get("income_profile_gates")) != canonical_json(
        selected_contract["income_profile_gates"]
    ):
        raise RuntimeError("Funded income-profile gates differ from calibration")
    provenance = funded_contract.get("provenance")
    if not isinstance(provenance, dict):
        raise RuntimeError("Funded transition contract omits provenance")
    if bool(provenance.get("legacy_provenance_exception", True)) or provenance.get(
        "legacy_missing_fields"
    ) != []:
        raise RuntimeError("Funded transition contract used a legacy provenance exception")
    outside_share = float(selected_contract["outside_origin_entry_share"])
    for observed in (
        metadata.get("outside_origin_entry_share"),
        provenance.get("outside_origin_entry_share"),
        provenance.get("outside_origin_entry_share_reported"),
    ):
        if not close(observed, outside_share, relative=0.0, absolute=1.0e-15):
            raise RuntimeError("Funded outside-origin entry share differs from calibration")
    if metadata.get("outside_origin_entry_share_status") != (
        "selected transition-report contract"
    ):
        raise RuntimeError("Funded packet used a legacy outside-share exception")
    selected_code = selected_contract["code_fingerprints"]
    for field in (
        "calibration_code_fingerprints_reported",
        "calibration_code_fingerprints_live",
    ):
        if canonical_json(provenance.get(field)) != canonical_json(selected_code):
            raise RuntimeError(f"Funded {field} differs from calibration")
    if canonical_json(metadata.get("transition_report_provenance")) != canonical_json(
        provenance
    ):
        raise RuntimeError("Funded metadata duplicates conflicting report provenance")

    funded_code_hash = require_sha256(
        metadata.get("code_bundle_sha256"), "funded-policy code-bundle"
    )
    funded_code_files = metadata.get("code_bundle_files")
    if not isinstance(funded_code_files, list) or not funded_code_files:
        raise RuntimeError("Funded-policy metadata omits its code-bundle file list")
    funded_code_files = [str(path) for path in funded_code_files]
    if live_funded_code_bundle(funded_code_files) != funded_code_hash:
        raise RuntimeError("Live funded-policy code differs from the completed packet")
    launch_checks = {
        "report_sha256": report_hash,
        "source_sha256": selected_contract["source_sha256"],
        "target_set": selected_contract["target_set"],
        "target_fingerprint": selected_contract["target_fingerprint"],
        "model_profile": selected_contract["model_profile"]["name"],
        "code_bundle_sha256": funded_code_hash,
    }
    for field, expected in launch_checks.items():
        if str(launch.get(field)) != str(expected):
            raise RuntimeError(f"Funded launch contract differs at {field}")
    if canonical_json(launch.get("code_bundle_files")) != canonical_json(funded_code_files):
        raise RuntimeError("Funded launch and runtime code-bundle files differ")
    if not close(
        launch.get("outside_origin_entry_share"),
        outside_share,
        relative=0.0,
        absolute=1.0e-15,
    ):
        raise RuntimeError("Funded launch uses a different outside-origin share")
    numerical = metadata.get("numerical_gates")
    if not isinstance(numerical, dict):
        raise RuntimeError("Funded metadata omits numerical gates")
    supply_normalization = metadata.get("supply_normalization")
    if not isinstance(supply_normalization, dict):
        raise RuntimeError("Funded metadata omits its dated supply normalization")
    if (
        supply_normalization.get("status")
        != "empirical_initialization_normalization_not_recalibration"
        or supply_normalization.get("housing_supply_mode") != "static-elastic"
        or supply_normalization.get("rule")
        != "rescale the static supply intercept so supply at p0 equals date-0 demand"
    ):
        raise RuntimeError("Funded supply-normalization convention differs")
    supply_values = {
        field: require_finite(
            supply_normalization.get(field), f"supply normalization.{field}"
        )
        for field in (
            "retained_date0_asset_price",
            "date0_housing_demand",
            "date0_normalized_housing_stock",
            "original_stationary_supply_at_date0_price",
            "static_supply_intercept_scale",
            "date0_market_residual_after_normalization",
        )
    }
    if min(
        supply_values["retained_date0_asset_price"],
        supply_values["date0_housing_demand"],
        supply_values["date0_normalized_housing_stock"],
        supply_values["original_stationary_supply_at_date0_price"],
        supply_values["static_supply_intercept_scale"],
    ) <= 0.0:
        raise RuntimeError("Funded supply normalization has a nonpositive primitive")
    for observed, expected, label in (
        (
            supply_values["date0_normalized_housing_stock"],
            supply_values["date0_housing_demand"],
            "normalized date-0 stock",
        ),
        (
            supply_values["static_supply_intercept_scale"],
            supply_values["date0_housing_demand"]
            / supply_values["original_stationary_supply_at_date0_price"],
            "supply-intercept scale",
        ),
    ):
        if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
            raise RuntimeError(f"Funded {label} identity failed")
    if abs(supply_values["date0_market_residual_after_normalization"]) > 1.0e-12:
        raise RuntimeError("Funded date-0 supply normalization fails its market gate")
    for maximum, tolerance in (
        ("maximum_market_residual", "market_tolerance"),
        ("maximum_fiscal_residual", "fiscal_absolute_tolerance"),
        ("maximum_mass_residual", "mass_tolerance"),
    ):
        maximum_value = require_finite(numerical.get(maximum), maximum)
        tolerance_value = require_finite(numerical.get(tolerance), tolerance)
        if maximum_value > tolerance_value:
            raise RuntimeError(f"Funded metadata fails {maximum}")
    if int(numerical.get("maximum_nonfinite_reported_objects", -1)) != 0 or int(
        numerical.get("maximum_nonfinite_distribution_count", -1)
    ) != 0:
        raise RuntimeError("Funded metadata reports nonfinite model objects")
    if require_finite(
        numerical.get("minimum_distribution_mass"), "minimum distribution mass"
    ) < -1.0e-14:
        raise RuntimeError("Funded metadata reports negative distribution mass")
    if require_finite(
        numerical.get("maximum_feasibility_frontier_projection_mass"),
        "maximum feasibility projection",
    ) > 1.0e-6:
        raise RuntimeError("Funded metadata fails its feasibility-projection gate")
    canonical_paths = validate_funded_paths_and_summaries(
        raw_paths, summaries, metadata, directory, best
    )
    path_by_case = {
        case: [row for row in canonical_paths if row["case"] == case]
        for case in EXPECTED_FUNDED_CASES
    }
    baseline = path_by_case[EXPECTED_FUNDED_CASES[0]]
    renewal = metadata["renewal"]
    date0_identities = (
        (
            baseline[0]["asset_price"],
            supply_values["retained_date0_asset_price"],
            "date-0 retained price",
        ),
        (
            baseline[0]["housing_demand"],
            supply_values["date0_housing_demand"],
            "date-0 housing demand",
        ),
        (
            baseline[0]["effective_mature_entrant_flow_B"],
            renewal["topcode_adjusted_old_mature_flow"],
            "date-0 adjusted mature flow",
        ),
        (
            baseline[0]["raw_state_scheduled_mature_entrant_flow_B"],
            renewal["raw_old_mature_flow"],
            "date-0 raw mature flow",
        ),
    )
    for observed, expected, label in date0_identities:
        if not close(observed, expected, absolute=3.0e-10, relative=3.0e-10):
            raise RuntimeError(f"Funded {label} identity failed")
    for case in EXPECTED_FUNDED_CASES[1:]:
        for period in range(5):
            stripped_baseline = {
                key: value
                for key, value in baseline[period].items()
                if key not in {"case", "scenario", "policy_case"}
            }
            stripped_case = {
                key: value
                for key, value in path_by_case[case][period].items()
                if key not in {"case", "scenario", "policy_case"}
            }
            if canonical_json(stripped_case) != canonical_json(stripped_baseline):
                raise RuntimeError(
                    f"{case}: prepolicy row {period} is not exactly the shared history"
                )
    identity_fields = (
        "asset_price",
        "adult_population",
        "birth_children_topcode_adjusted",
        "entry_flow_E",
        "owner_rate",
        "dependent_child_owner_rate",
    )
    identity_gap = max(
        abs(float(path_by_case[case][period][field]) - float(baseline[period][field]))
        for case in EXPECTED_FUNDED_CASES[1:]
        for period in range(5)
        for field in identity_fields
    )
    recorded_identity_gap = require_finite(
        metadata.get("prepolicy_path_identity_max_abs_gap"),
        "prepolicy identity gap",
    )
    if identity_gap > 1.0e-12 or not close(
        recorded_identity_gap,
        identity_gap,
        absolute=1.0e-14,
        relative=0.0,
    ):
        raise RuntimeError("Funded paths do not share the same history through 2023")
    timing_gates = metadata.get("matched_2023_and_policy_timing_gates")
    expected_timing_gates = {
        "matched_2023_population_gap": float(baseline[4]["population_index"])
        - float(best["terminal_population_index"]),
        "matched_2023_housing_cost_gap": float(baseline[4]["asset_price_index"])
        - float(best["terminal_housing_cost_index"]),
    }
    if not isinstance(timing_gates, dict) or set(timing_gates) != {
        "matched_2023_population_gap",
        "matched_2023_housing_cost_gap",
    }:
        raise RuntimeError("Funded packet omits its matched-2023 timing gates")
    for field, expected in expected_timing_gates.items():
        observed = require_finite(timing_gates[field], f"matched-2023 gate {field}")
        if abs(observed) > 2.0e-8 or not close(
            observed, expected, absolute=3.0e-12, relative=3.0e-12
        ):
            raise RuntimeError("Funded packet fails its matched-2023 timing gates")
    if max(abs(value) for value in expected_timing_gates.values()) > 2.0e-8:
        raise RuntimeError("Funded packet fails its matched-2023 timing gates")
    if not close(
        baseline[4]["population_index"],
        best["terminal_population_index"],
        absolute=2.0e-8,
    ) or not close(
        baseline[4]["asset_price_index"],
        best["terminal_housing_cost_index"],
        absolute=2.0e-8,
    ):
        raise RuntimeError("Funded 2023 state differs from the selected calibration")

    if len(endpoint_rows) != 1:
        raise RuntimeError("Stationary endpoint CSV must contain exactly one row")
    endpoint = metadata.get("stationary_open_endpoint")
    if not isinstance(endpoint, dict) or endpoint.get("status") != "complete":
        raise RuntimeError("Funded metadata omits a complete stationary endpoint")
    assert_endpoint_row_matches(endpoint, endpoint_rows[0])
    validate_stationary_endpoint_accounting(endpoint, metadata, best)
    validate_stationary_endpoint_diagnostics(directory, endpoint)
    endpoint_contract = metadata["stationary_open_endpoint_contract"]
    if endpoint_contract.get("policy") != "unrebated 1 percent status quo (no reform)":
        raise RuntimeError("Stationary endpoint uses the wrong policy convention")
    for required_file in (
        directory
        / "stationary_endpoint_status_quo"
        / "stationary_pure_branch_endpoint.csv",
        directory
        / "stationary_endpoint_status_quo"
        / "stationary_endpoint_search.csv",
        directory / "funded_policy_transition_paths.png",
        directory / "funded_policy_transition_paths.pdf",
    ):
        if not required_file.is_file() or required_file.stat().st_size <= 0:
            raise RuntimeError(f"Funded final packet omits diagnostic {required_file}")
    validate_solved_model_profile(metadata, best)
    operator_gates = metadata.get("operator_gates")
    if not isinstance(operator_gates, dict):
        raise RuntimeError("Funded metadata omits stationary operator gates")
    operator_names = (
        "stationary_post_fertility_nesting_l1",
        "one_step_constant_path_nesting_l1",
        "mature_flow_abs_error",
        "birth_flow_abs_error",
        "topcode_adjusted_birth_flow_abs_error",
    )
    if max(
        require_finite(operator_gates.get(name), f"operator gate {name}")
        for name in operator_names
    ) > 5.0e-9:
        raise RuntimeError("Funded stationary-operator gate failed")
    return {
        "metadata": metadata,
        "launch_contract": launch,
        "launcher_status": launcher_status,
        "summaries": summaries,
        "paths": canonical_paths,
        "stationary_endpoint": endpoint,
        "contract": {
            "transition_report": str(selected_report),
            "transition_report_sha256": report_hash,
            "funded_code_bundle_sha256": funded_code_hash,
            "funded_code_bundle_files": funded_code_files,
            "policy_start_year": 2027,
            "prepolicy_fiscal_convention": metadata["prepolicy_fiscal_convention"],
            "policy_fiscal_rule": metadata["policy_fiscal_rule"],
            "housing_price_expectations": metadata["housing_price_expectations"],
            "outside_origin_entry_share": outside_share,
        },
    }


def main() -> None:
    args = parse_args()
    selected_report = args.selected_report.resolve()
    funded_dir = args.funded_policy_dir.resolve()
    round_dirs = [path.resolve() for path in args.round_dir]
    repeat_dirs = [path.resolve() for path in args.repeat_dir]
    output_dir = args.output_dir.resolve()

    collector, best, selected_contract, calibration_provenance = (
        load_selected_report(selected_report)
    )
    live_code = live_calibration_code_fingerprints()
    if canonical_json(live_code) != canonical_json(
        selected_contract["code_fingerprints"]
    ):
        raise RuntimeError(
            "Live calibration code differs from the selected recorded contract; "
            "the parameter domain cannot be interpreted safely"
        )
    domain = transition_domain(best)
    parameters = parameter_table(best, domain)
    best_target_rows = target_fit_rows(selected_report.parent / "best_target_fit.csv")
    recomputed_loss = validate_selected_target_fit(
        best_target_rows, best, selected_contract
    )
    best_fit = fit_lookup(best_target_rows)

    records = candidate_records(round_dirs)
    for record in records:
        assert_same_contract(selected_contract, record["contract"], "Search trace")
    if records and not any(
        canonical_json(record["theta"]) == canonical_json(best["theta"])
        and math.isclose(
            float(record["transition_loss"]),
            float(best["transition_loss"]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        for record in records
    ):
        raise RuntimeError("Selected collector winner is absent from the supplied search trace")

    anchor = best
    anchor_target_rows = best_target_rows
    show_anchor = False
    if records:
        anchors = [
            row for row in records if row["round"] == 1 and int(row["task_id"]) == 1
        ]
        if len(anchors) != 1:
            raise RuntimeError("The optional search trace lacks one round-one anchor")
        anchor = anchors[0]
        anchor_target_rows = target_fit_rows(Path(anchor["task_dir"]))
        assert_same_target_fit_contract(
            best_target_rows, anchor_target_rows, "Search anchor"
        )
        show_anchor = True
    anchor_fit = fit_lookup(anchor_target_rows)
    comparison: list[dict[str, Any]] = []
    for row in best_target_rows:
        moment = str(row["moment"])
        a, b = anchor_fit[moment], best_fit[moment]
        comparison.append(
            {
                "moment": moment,
                "target": b["target"],
                "anchor_model": a["model"],
                "anchor_gap": a["gap"],
                "anchor_loss_contribution": a["loss_contribution"],
                "anchor_standardized_gap": a["standardized_gap"],
                "best_model": b["model"],
                "best_gap": b["gap"],
                "weight": b["weight"],
                "best_loss_contribution": b["loss_contribution"],
                "best_standardized_gap": b["standardized_gap"],
            }
        )

    round_rows: list[dict[str, Any]] = []
    candidate_table: list[dict[str, Any]] = []
    for index, directory in enumerate(round_dirs, start=1):
        subset = [row for row in records if row["round"] == index]
        if not subset:
            raise RuntimeError(f"Search round contains no completed candidates: {directory}")
        round_rows.append(
            {
                "round": index,
                "round_name": directory.name,
                "completed_candidates": len(subset),
                "best_loss": min(float(row["transition_loss"]) for row in subset),
            }
        )
    for row in sorted(records, key=lambda item: float(item["transition_loss"])):
        candidate_table.append(
            {
                "round": row["round"],
                "round_name": row["round_name"],
                "task_id": row["task_id"],
                "design": row["design"],
                "transition_loss": row["transition_loss"],
                "new_psi_child": row["new_psi_child"],
                "old_psi_child": row["old_psi_child"],
                "terminal_tfr": row["terminal_tfr"],
                "terminal_childless_rate": row["terminal_childless_rate"],
                "terminal_mean_age_first_birth": row["terminal_mean_age_first_birth"],
                "terminal_share_first_births_age30plus": row[
                    "terminal_share_first_births_age30plus"
                ],
                "terminal_housing_cost_index": row["terminal_housing_cost_index"],
                "theta_json": row["theta_json"],
            }
        )

    repeat_rows: list[dict[str, Any]] = []
    for directory in repeat_dirs:
        repeat_summary_path = directory / "task_001" / "summary.json"
        repeat_summary = read_json(repeat_summary_path)
        repeat = dict(repeat_summary["best_candidate"])
        repeat_contract = contract_from_payload(repeat_summary, repeat)
        assert_same_contract(selected_contract, repeat_contract, f"Repeat {directory.name}")
        repeat_target_rows = target_fit_rows(repeat_summary_path.parent)
        assert_same_target_fit_contract(
            best_target_rows, repeat_target_rows, f"Repeat {directory.name}"
        )
        if canonical_json(repeat["theta"]) != canonical_json(best["theta"]):
            raise RuntimeError(f"Repeat {directory.name} uses a different theta")
        repeat_old = float(repeat_summary["old_psi_child"])
        for label, observed, expected in (
            ("old_psi_child", repeat_old, best["old_psi_child"]),
            ("new_psi_child", repeat["new_psi_child"], best["new_psi_child"]),
        ):
            if not math.isclose(
                float(observed), float(expected), rel_tol=0.0, abs_tol=1.0e-14
            ):
                raise RuntimeError(f"Repeat {directory.name} differs in {label}")
        repeat_fit = fit_lookup(repeat_target_rows)
        max_moment_gap = max(
            abs(repeat_fit[name]["model"] - best_fit[name]["model"])
            for name in best_fit
        )
        loss_gap = float(repeat["transition_loss"]) - float(best["transition_loss"])
        if abs(loss_gap) > 1.0e-10 or max_moment_gap > 1.0e-10:
            raise RuntimeError(
                f"Repeat {directory.name} fails deterministic reproduction: "
                f"loss_gap={loss_gap:.3e}, moment_gap={max_moment_gap:.3e}"
            )
        repeat_rows.append(
            {
                "repeat": directory.name,
                "loss": float(repeat["transition_loss"]),
                "loss_gap_from_selected": loss_gap,
                "max_market_residual": float(repeat["max_market_residual"]),
                "max_mass_residual": float(repeat["max_mass_residual"]),
                "max_abs_target_moment_gap_from_selected": max_moment_gap,
                "theta_matches_selected": True,
                "contract_matches_selected": True,
                "status": "passed_exact_repeat_gate_1e-10",
            }
        )

    funded = load_funded_policy_packet(
        funded_dir,
        selected_report=selected_report,
        selected_contract=selected_contract,
        best=best,
    )
    transition_rows = sorted(
        [
            row
            for row in funded["paths"]
            if row["case"] == EXPECTED_FUNDED_CASES[0] and int(row["period"]) <= 4
        ],
        key=lambda row: int(row["period"]),
    )
    validate_selected_transition_path(transition_rows, best)

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "best_target_fit.csv", best_target_rows)
    write_csv(output_dir / "target_fit_anchor_vs_best.csv", comparison)
    write_csv(output_dir / "best_parameter_table.csv", parameters)
    write_csv(output_dir / "best_transition_path.csv", transition_rows)
    write_csv(output_dir / "funded_policy_summary.csv", funded["summaries"])
    write_csv(output_dir / "funded_policy_paths.csv", funded["paths"])
    write_csv(
        output_dir / "stationary_open_endpoint.csv",
        [funded["stationary_endpoint"]],
    )
    if candidate_table:
        write_csv(output_dir / "all_completed_candidates.csv", candidate_table)
        write_csv(output_dir / "search_rounds.csv", round_rows)
    if repeat_rows:
        write_csv(output_dir / "repeat_checks.csv", repeat_rows)
    make_figures(
        output_dir,
        comparison,
        transition_rows,
        round_rows,
        show_anchor=show_anchor,
    )
    make_funded_policy_figure(output_dir, funded["paths"])

    best_loss = float(best["transition_loss"])
    anchor_loss = float(anchor["transition_loss"])
    result = {
        "status": "complete_canonical_transition_quantitative_packet",
        "claim_status": "calibrated_along_transition_with_funded_policy_paths",
        "selected_report": str(selected_report),
        "selected_report_sha256": sha256(selected_report),
        "selected_collector_status": collector["status"],
        "selected_collector_classified_invalid_candidates": len(
            collector.get("failed_or_missing_tasks", [])
        ),
        "selected_collector_failure_ledger": collector.get(
            "failed_or_missing_tasks", []
        ),
        "contract": selected_contract,
        "source_recorded_path": selected_contract["source"],
        "source_sha256": selected_contract["source_sha256"],
        "source_verification": "selected_report_funded_runtime_and_launcher_sha256_agree",
        "target_set": selected_contract["target_set"],
        "target_fingerprint": selected_contract["target_fingerprint"],
        "target_count": selected_contract["target_count"],
        "model_profile": selected_contract["model_profile"],
        "income_profile_gates": selected_contract["income_profile_gates"],
        "code_fingerprints": selected_contract["code_fingerprints"],
        "live_code_fingerprints": live_code,
        "outside_origin_entry_share": selected_contract[
            "outside_origin_entry_share"
        ],
        "outside_origin_entry_status": (
            "externally_fixed_provisional_across_market_anchor"
        ),
        "transition_free_parameter_count": len(domain),
        "transition_domain": [list(row) for row in domain],
        "transition_domain_code_bundle_sha256": live_code["bundle_sha256"],
        "rounds": round_rows,
        "completed_candidate_count": len(records) if records else None,
        "anchor_loss": anchor_loss if show_anchor else None,
        "best_loss": best_loss,
        "recomputed_best_loss": recomputed_loss,
        "loss_improvement": anchor_loss - best_loss if show_anchor else None,
        "loss_improvement_percent": (
            100.0 * (anchor_loss - best_loss) / anchor_loss
            if show_anchor and anchor_loss != 0.0
            else None
        ),
        "best_candidate": best,
        "repeat_checks": repeat_rows,
        "funded_policy_contract": funded["contract"],
        "funded_policy_summary": funded["summaries"],
        "stationary_open_endpoint": funded["stationary_endpoint"],
        "stationary_open_endpoint_contract": funded["metadata"][
            "stationary_open_endpoint_contract"
        ],
        "population_bridge": calibration_provenance["population_bridge"],
        "population_validation_status": calibration_provenance[
            "population_validation_status"
        ],
        "target_measurements": calibration_provenance["target_measurements"],
        "largest_best_loss_contributions": sorted(
            [
                {
                    "moment": row["moment"],
                    "loss_contribution": row["best_loss_contribution"],
                }
                for row in comparison
            ],
            key=lambda row: float(row["loss_contribution"]),
            reverse=True,
        ),
    }
    write_json(output_dir / "summary.json", result)
    (output_dir / "README.md").write_text(
        "\n".join(
            [
                "# Canonical transition calibration and funded-policy packet",
                "",
                "This is the single paper-facing quantitative report. The model starts",
                "from the old steady state, reweights the 2007 household-age distribution,",
                "introduces the reduced-form fertility-preference change, and measures all",
                "twelve calibration moments on the simulated 2023 economy. Census household",
                "totals and ACS householder-age shares form an imposed population bridge;",
                "post-2007 fertility alone does not generate those historical population data.",
                "",
                f"The selected {selected_contract['model_profile']['name']} calibration has",
                f"loss {best_loss:.6f} and {len(domain)} free transition coordinates. The",
                "complete target table reports every target, model moment, gap, weight, and",
                "loss contribution. The parameter table reports every free coordinate with",
                "its bounds and boundary status, followed by derived and external objects.",
                "",
                "The four funded-policy paths share the matched 2023 history and introduce",
                "reforms at the first subsequent model date, 2027. Property-tax revenue pays",
                "equal household rebates and, where applicable, the targeted purchase grant",
                "at every policy date. These are temporary-equilibrium paths with static",
                "housing-price expectations, not perfect-foresight or welfare calculations.",
                "",
                "The finite transition endpoint is not labeled a steady state. The separate",
                "`stationary_open_endpoint.csv` is the exact status-quo open-population steady",
                "state under the terminal fertility preference, dated static-elastic housing",
                "supply, and fixed migration closure. Its market-clearing method is recorded",
                "explicitly because a discrete policy jump can require persistent entry types.",
                "",
                "Primary files are `summary.json`, `best_target_fit.csv`,",
                "`best_parameter_table.csv`, `best_transition_path.csv`,",
                "`funded_policy_summary.csv`, `funded_policy_paths.csv`,",
                "`stationary_open_endpoint.csv`, and the dated-fit and funded-path figures.",
                "Search-trace and repeat files appear only when those optional inputs are",
                "supplied. The outside-origin entry share remains an externally fixed,",
                "provisional closure normalization rather than a national migration estimate.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(
        "TRANSITION_CALIBRATION_REPORT_COMPLETE "
        f"candidates={len(records)} best={best_loss:.6f} "
        f"repeats={len(repeat_rows)} funded_cases={len(funded['summaries'])}",
        flush=True,
    )


if __name__ == "__main__":
    main()
