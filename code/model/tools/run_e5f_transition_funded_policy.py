#!/usr/bin/env python3
"""Run funded property-tax policies along the matched E5F transition.

The household model and calendar operator are unchanged.  This wrapper loads
the selected dated-transition parameter vector, reconstructs the same 2007
initial state and Census/ACS bridge, preserves the selected matched 2023
economy, and introduces one of the established funded policy cases at the
first subsequent model date, 2027.  At every funded-policy date it jointly
clears the contemporaneous housing market and the government budget by solving
for the house price and equal lump-sum transfer.

This remains the paper's temporary-equilibrium transition: households treat
current prices and policy primitives as permanent when making each date's
choices.  It is not a perfect-foresight asset-pricing or welfare calculation.
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
from typing import Any, Callable

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for path in (MODEL_ROOT, TOOLS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import audit_closed_reproductive_closure as closure
import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_transition_calibration as calibration

from intergen_eqscale_seq_optimized.e5_profile import e5_target_system


DEFAULT_REPORT = ROOT / "output/model/e5f_transition_calibration_report/summary.json"
DEFAULT_OUTDIR = ROOT / "output/model/e5f_transition_funded_policy"
POLICY_CASES: dict[str, tuple[float, bool, bool]] = {
    "unrebated_tax1_status_quo": (0.01, False, False),
    "rebated_tax1_baseline": (0.01, False, True),
    "rebated_tax2": (0.02, False, True),
    "rebated_tax2_grant0p4_Hge6": (0.02, True, True),
}
BASELINE_CASE = "unrebated_tax1_status_quo"
CASE_LABELS = {
    "unrebated_tax1_status_quo": "Status quo: 1% tax, no rebate",
    "rebated_tax1_baseline": "Funded 1% tax",
    "rebated_tax2": "Funded 2% tax",
    "rebated_tax2_grant0p4_Hge6": "Funded 2% tax + purchase grant",
}
CASE_COLORS = {
    "unrebated_tax1_status_quo": "#1f1f1f",
    "rebated_tax1_baseline": "#377eb8",
    "rebated_tax2": "#e41a1c",
    "rebated_tax2_grant0p4_Hge6": "#4daf4a",
}
CASE_LINESTYLES = {
    "unrebated_tax1_status_quo": "-",
    "rebated_tax1_baseline": "--",
    "rebated_tax2": "-.",
    "rebated_tax2_grant0p4_Hge6": ":",
}
POLICY_START_PERIOD = calibration.TRANSITION_PERIODS + 1
POLICY_START_YEAR = calibration.START_YEAR + 4 * POLICY_START_PERIOD
FISCAL_ABSOLUTE_TOLERANCE = 2.5e-5
MASS_TOLERANCE = 2.0e-8
TOOL_CODE_BUNDLE_FILES = (
    "code/model/tools/run_e5f_transition_funded_policy.py",
    "code/model/tools/run_e5f_transition_calibration.py",
    "code/model/tools/run_e5f_open_population_transition.py",
    "code/model/tools/run_dynamic_population_transition.py",
    "code/model/tools/audit_closed_reproductive_closure.py",
)
PACKAGE_CODE_BUNDLE_FILES = tuple(
    path.relative_to(ROOT).as_posix()
    for path in sorted(
        (MODEL_ROOT / "intergen_eqscale_seq_optimized").glob("*.py"),
        key=lambda item: item.name,
    )
)
CODE_BUNDLE_FILES = TOOL_CODE_BUNDLE_FILES + PACKAGE_CODE_BUNDLE_FILES


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=closure.E5F_FLOOR_SOURCE)
    parser.add_argument("--transition-report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--expected-report-sha256", default=None)
    parser.add_argument("--expected-source-sha256", default=None)
    parser.add_argument("--expected-target-set", default=None)
    parser.add_argument("--expected-target-fingerprint", default=None)
    parser.add_argument("--expected-code-bundle-sha256", default=None)
    parser.add_argument(
        "--expected-model-profile",
        choices=(calibration.BASELINE_MODEL_PROFILE, calibration.REPAIRED_MODEL_PROFILE),
        default=None,
    )
    parser.add_argument(
        "--case",
        action="append",
        choices=tuple(POLICY_CASES),
        help="Case(s) to run. The unrebated 1 percent status quo is always included.",
    )
    parser.add_argument("--post-2023-periods", type=int, default=40)
    parser.add_argument(
        "--outside-origin-entry-share",
        type=float,
        required=True,
        help="Externally fixed share; must equal the selected transition report.",
    )
    parser.add_argument("--market-tol", type=float, default=2.0e-4)
    parser.add_argument("--market-max-iter", type=int, default=30)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help=(
            "Keep exact J=17 and Nb=120 but always stop at the first reform date."
        ),
    )
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    return value


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def code_bundle_sha256(root: Path = ROOT) -> str:
    """Hash ordered relative paths and bytes for all model-critical transition code."""
    digest = hashlib.sha256()
    for relative in CODE_BUNDLE_FILES:
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"Code-bundle file is missing: {path}")
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def selected_cases(raw: list[str] | None) -> list[str]:
    requested = list(raw or POLICY_CASES)
    ordered = [BASELINE_CASE]
    ordered.extend(case for case in POLICY_CASES if case in requested and case != BASELINE_CASE)
    return ordered


def effective_post_2023_periods(requested: int, smoke: bool) -> int:
    """A smoke always includes t=5, the first date on which reform is active."""
    return 1 if smoke else int(requested)


def report_model_profile(
    payload: dict[str, Any],
    best: dict[str, Any],
) -> dict[str, Any]:
    """Read the profile contract, defaulting only legacy reports to E5F floor."""
    raw_profiles = [
        value
        for value in (best.get("model_profile"), payload.get("model_profile"))
        if value is not None
    ]
    if not raw_profiles:
        return {
            "name": calibration.BASELINE_MODEL_PROFILE,
            "legacy_default": True,
            "report_metadata": {},
        }
    names: list[str] = []
    metadata: dict[str, Any] = {}
    for raw in raw_profiles:
        if isinstance(raw, str):
            names.append(raw)
        elif isinstance(raw, dict):
            name = raw.get("name")
            if not isinstance(name, str) or not name:
                raise ValueError(f"Malformed report model_profile: {raw}")
            names.append(name)
            for key, value in raw.items():
                if key in metadata and json.dumps(
                    jsonable(metadata[key]), sort_keys=True
                ) != json.dumps(jsonable(value), sort_keys=True):
                    raise RuntimeError(
                        f"Conflicting model-profile field {key}: "
                        f"{metadata[key]!r} versus {value!r}"
                    )
                metadata[key] = value
        else:
            raise ValueError(f"Malformed report model_profile: {raw!r}")
    if len(set(names)) != 1:
        raise RuntimeError(f"Conflicting model profiles in transition report: {names}")
    name = names[0]
    if name not in {
        calibration.BASELINE_MODEL_PROFILE,
        calibration.REPAIRED_MODEL_PROFILE,
    }:
        raise ValueError(f"Unsupported transition-report model profile: {name}")
    metadata["name"] = name
    return {"name": name, "legacy_default": False, "report_metadata": metadata}


def _report_contract_dicts(
    payload: dict[str, Any],
    best: dict[str, Any],
) -> list[dict[str, Any]]:
    dictionaries = [best, payload]
    for parent in (best, payload):
        nested = parent.get("contract")
        if isinstance(nested, dict):
            dictionaries.append(nested)
    return dictionaries


def report_outside_origin_share(
    payload: dict[str, Any],
    best: dict[str, Any],
) -> float | None:
    values = [
        float(item["outside_origin_entry_share"])
        for item in _report_contract_dicts(payload, best)
        if item.get("outside_origin_entry_share") is not None
    ]
    if not values:
        return None
    if any(not math.isfinite(value) for value in values):
        raise RuntimeError("Transition report has a nonfinite outside-origin entry share")
    if max(values) - min(values) > 1.0e-15:
        raise RuntimeError(f"Conflicting outside-origin entry shares in report: {values}")
    return values[0]


def report_calibration_code_fingerprints(
    payload: dict[str, Any],
    best: dict[str, Any],
) -> dict[str, Any] | None:
    candidates: list[dict[str, Any]] = []
    for item in _report_contract_dicts(payload, best):
        for key in (
            "code_fingerprints",
            "calibration_code_fingerprints",
            "code_fingerprint_contract",
        ):
            value = item.get(key)
            if isinstance(value, dict) and value.get("bundle_sha256") is not None:
                candidates.append(value)
    if not candidates:
        return None
    bundles = [str(value["bundle_sha256"]) for value in candidates]
    if len(set(bundles)) != 1:
        raise RuntimeError(f"Conflicting calibration code bundles in report: {bundles}")
    if len(bundles[0]) != 64 or any(character not in "0123456789abcdef" for character in bundles[0]):
        raise RuntimeError(f"Malformed calibration code-bundle SHA-256: {bundles[0]}")
    file_maps = [value["files"] for value in candidates if isinstance(value.get("files"), dict)]
    canonical_files = [
        json.dumps(value, sort_keys=True, separators=(",", ":")) for value in file_maps
    ]
    if len(set(canonical_files)) > 1:
        raise RuntimeError("Conflicting calibration file fingerprints in report")
    richest = max(candidates, key=lambda value: int(isinstance(value.get("files"), dict)))
    return dict(richest)


def report_income_profile_gates(
    payload: dict[str, Any],
    best: dict[str, Any],
) -> dict[str, Any] | None:
    candidates = [
        item["income_profile_gates"]
        for item in _report_contract_dicts(payload, best)
        if item.get("income_profile_gates") is not None
    ]
    if not candidates:
        return None
    if any(not isinstance(value, dict) for value in candidates):
        raise RuntimeError("Malformed transition-report income_profile_gates")
    canonical = [
        json.dumps(jsonable(value), sort_keys=True, separators=(",", ":"))
        for value in candidates
    ]
    if len(set(canonical)) != 1:
        raise RuntimeError("Conflicting income_profile_gates in transition report")
    return dict(candidates[0])


def validate_transition_report_provenance(
    payload: dict[str, Any],
    best: dict[str, Any],
    model_profile: dict[str, Any],
    *,
    cli_outside_origin_share: float,
    live_calibration_code_fingerprints: dict[str, Any],
    expected_model_profile: str | None,
) -> dict[str, Any]:
    """Gate transition-only closure objects before any household solve."""
    explicit_legacy = bool(
        model_profile.get("legacy_default")
        and expected_model_profile == calibration.BASELINE_MODEL_PROFILE
    )
    if model_profile.get("legacy_default") and not explicit_legacy:
        raise RuntimeError(
            "A legacy transition report without model_profile must be routed "
            "explicitly with --expected-model-profile e5f-floor"
        )
    recorded_outside = report_outside_origin_share(payload, best)
    recorded_code = report_calibration_code_fingerprints(payload, best)
    missing = []
    if recorded_outside is None:
        missing.append("outside_origin_entry_share")
    if recorded_code is None:
        missing.append("code_fingerprints.bundle_sha256")
    if missing and not explicit_legacy:
        raise RuntimeError(
            "Selected transition report omits production provenance fields: "
            + ", ".join(missing)
        )
    cli_share = float(cli_outside_origin_share)
    if recorded_outside is not None and not math.isclose(
        recorded_outside, cli_share, rel_tol=0.0, abs_tol=1.0e-15
    ):
        raise RuntimeError(
            "Outside-origin entry share gate failed: "
            f"report={recorded_outside}, CLI={cli_share}"
        )
    live_bundle = str(live_calibration_code_fingerprints.get("bundle_sha256", ""))
    if len(live_bundle) != 64:
        raise RuntimeError("Live calibration code fingerprint is malformed")
    recorded_bundle = None if recorded_code is None else str(recorded_code["bundle_sha256"])
    if recorded_bundle is not None and recorded_bundle != live_bundle:
        raise RuntimeError(
            "Calibration code-bundle gate failed: "
            f"report={recorded_bundle}, live={live_bundle}"
        )
    if recorded_code is not None and model_profile["name"] == calibration.REPAIRED_MODEL_PROFILE:
        recorded_files = recorded_code.get("files")
        live_files = live_calibration_code_fingerprints.get("files")
        if not isinstance(recorded_files, dict) or not isinstance(live_files, dict):
            raise RuntimeError("Repaired report must record calibration code-fingerprint files")
        if recorded_files != live_files:
            raise RuntimeError("Recorded calibration file fingerprints disagree with live code")
    return {
        "outside_origin_entry_share": (
            cli_share if recorded_outside is None else recorded_outside
        ),
        "outside_origin_entry_share_reported": recorded_outside,
        "calibration_code_fingerprints_reported": recorded_code,
        "calibration_code_fingerprints_live": live_calibration_code_fingerprints,
        "legacy_provenance_exception": bool(missing),
        "legacy_missing_fields": missing,
    }


def activate_contract_profile(
    contract: dict[str, Any],
    theta: dict[str, float],
) -> tuple[tuple[tuple[str, float, float, str], ...], dict[str, Any], dict[str, Any]]:
    """Activate exactly the profile declared by the selected transition report."""
    name = str(contract["model_profile"]["name"])
    if name == calibration.REPAIRED_MODEL_PROFILE and "first_birth_fixed_cost" not in theta:
        raise RuntimeError(
            "An e5f-income-entry report must estimate and save first_birth_fixed_cost"
        )
    if name == calibration.REPAIRED_MODEL_PROFILE:
        cost = float(theta["first_birth_fixed_cost"])
        if not math.isfinite(cost) or not 0.0 <= cost <= 8.0:
            raise RuntimeError(
                "The estimated first_birth_fixed_cost must lie in its declared [0,8] domain"
            )
    domain, overrides, metadata = calibration.activate_model_profile(name, theta)
    declared = dict(contract["model_profile"].get("report_metadata") or {})
    if declared.get("profile_id") not in (None, metadata.get("profile_id")):
        raise RuntimeError(
            "Report profile identifier disagrees with the live calibration contract: "
            f"report={declared.get('profile_id')}, live={metadata.get('profile_id')}"
        )
    if name == calibration.REPAIRED_MODEL_PROFILE:
        required = (
            "profile_id",
            "permanent_income_log_variance",
            "income_state_count",
            "first_birth_fixed_cost",
            "first_birth_fixed_cost_semantics",
        )
        missing = [field for field in required if field not in declared]
        if missing:
            raise RuntimeError(
                "Income-entry report omits profile contract fields: "
                + ", ".join(missing)
            )
        if not math.isclose(
            float(declared["permanent_income_log_variance"]),
            float(metadata["permanent_income_log_variance"]),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ):
            raise RuntimeError(
                "Report permanent-income variance disagrees with the live profile"
            )
        if int(declared["income_state_count"]) != int(metadata["income_state_count"]):
            raise RuntimeError("Report income-state count disagrees with the live profile")
        for field in ("first_birth_fixed_cost", "first_birth_fixed_cost_semantics"):
            if str(declared[field]) != str(metadata[field]):
                raise RuntimeError(
                    f"Report {field} disagrees with the live profile: "
                    f"report={declared[field]!r}, live={metadata[field]!r}"
                )
        if not bool(overrides.get("permanent_income_levels_enabled", False)):
            raise RuntimeError("Income-entry routing did not activate permanent income levels")
        if len(np.asarray(overrides.get("z_grid", []), dtype=float)) != 15:
            raise RuntimeError("Income-entry routing did not recover the measured 15-state income grid")
    return domain, overrides, metadata


def load_transition_contract(
    report: Path,
    source: Path,
    *,
    expected_report_sha256: str | None = None,
    expected_source_sha256: str | None = None,
    expected_target_set: str | None = None,
    expected_target_fingerprint: str | None = None,
    expected_model_profile: str | None = None,
    outside_origin_entry_share: float,
    live_calibration_code_fingerprints: dict[str, Any],
) -> dict[str, Any]:
    payload = json.loads(report.read_text(encoding="utf-8"))
    best = payload.get("best_candidate")
    if not isinstance(best, dict) or not isinstance(best.get("theta"), dict):
        raise ValueError(f"Transition report has no best_candidate.theta: {report}")
    theta = {str(key): float(value) for key, value in best["theta"].items()}
    if not theta or any(not math.isfinite(value) for value in theta.values()):
        raise ValueError("Selected transition theta is empty or nonfinite")
    actual_report_hash = sha256(report)
    if expected_report_sha256 is not None and actual_report_hash != expected_report_sha256:
        raise RuntimeError(
            "Report fingerprint gate failed: "
            f"actual={actual_report_hash}, expected={expected_report_sha256}"
        )
    actual_source_hash = sha256(source)
    expected_source_hash = str(best.get("source_sha256", ""))
    if actual_source_hash != expected_source_hash:
        raise RuntimeError(
            "Source fingerprint gate failed: "
            f"actual={actual_source_hash}, expected={expected_source_hash}"
        )
    if expected_source_sha256 is not None and actual_source_hash != expected_source_sha256:
        raise RuntimeError(
            "Explicit source fingerprint gate failed: "
            f"actual={actual_source_hash}, expected={expected_source_sha256}"
        )
    targets = e5_target_system()
    if str(best.get("target_set")) != str(targets.name):
        raise RuntimeError(
            f"Target-set gate failed: report={best.get('target_set')}, live={targets.name}"
        )
    if str(best.get("target_fingerprint")) != str(targets.fingerprint):
        raise RuntimeError(
            "Target fingerprint gate failed: "
            f"report={best.get('target_fingerprint')}, live={targets.fingerprint}"
        )
    if expected_target_set is not None and str(targets.name) != expected_target_set:
        raise RuntimeError(
            f"Explicit target-set gate failed: live={targets.name}, expected={expected_target_set}"
        )
    if (
        expected_target_fingerprint is not None
        and str(targets.fingerprint) != expected_target_fingerprint
    ):
        raise RuntimeError(
            "Explicit target fingerprint gate failed: "
            f"live={targets.fingerprint}, expected={expected_target_fingerprint}"
        )
    model_profile = report_model_profile(payload, best)
    if expected_model_profile is not None and model_profile["name"] != expected_model_profile:
        raise RuntimeError(
            "Model-profile gate failed: "
            f"report={model_profile['name']}, expected={expected_model_profile}"
        )
    # Validate the declared repaired-profile variance, state count, and cost
    # semantics against the live profile before reading dependent gates.
    activate_contract_profile({"model_profile": model_profile}, dict(theta))
    provenance = validate_transition_report_provenance(
        payload,
        best,
        model_profile,
        cli_outside_origin_share=float(outside_origin_entry_share),
        live_calibration_code_fingerprints=live_calibration_code_fingerprints,
        expected_model_profile=expected_model_profile,
    )
    income_profile_gates = report_income_profile_gates(payload, best)
    if model_profile["name"] == calibration.REPAIRED_MODEL_PROFILE:
        if income_profile_gates is None:
            raise RuntimeError("Repaired transition report omits income_profile_gates")
        if not bool(income_profile_gates.get("permanent_income_levels_enabled", False)):
            raise RuntimeError("Repaired report does not activate permanent income levels")
        if int(income_profile_gates.get("income_state_count", -1)) != int(
            model_profile["report_metadata"]["income_state_count"]
        ):
            raise RuntimeError("Repaired report has inconsistent income-state gates")
        stationary_gap = float(
            income_profile_gates.get("stationary_weight_max_abs_gap", math.inf)
        )
        if not math.isfinite(stationary_gap) or stationary_gap > 1.0e-12:
            raise RuntimeError("Repaired report stationary-income-weight gate failed")
    return {
        "report": str(report),
        "report_sha256": actual_report_hash,
        "source": str(source),
        "source_sha256": actual_source_hash,
        "target_set": str(targets.name),
        "target_fingerprint": str(targets.fingerprint),
        "transition_loss": float(best["transition_loss"]),
        "theta": theta,
        "old_psi_child": float(best["old_psi_child"]),
        "new_psi_child": float(best["new_psi_child"]),
        "old_completed_fertility": float(best["old_completed_fertility"]),
        "matched_2023_population_index": float(best["terminal_population_index"]),
        "matched_2023_housing_cost_index": float(best["terminal_housing_cost_index"]),
        "model_profile": model_profile,
        "income_profile_gates": income_profile_gates,
        "provenance": provenance,
    }


def set_fiscal_policy(
    P: SimpleNamespace,
    *,
    annual_tax: float,
    grant: bool,
    transfer: float,
) -> None:
    """Set the exact fiscal primitives used by the stationary funded driver."""
    P.tau_H = float(annual_tax) * float(P.period_years)
    P.user_cost_rate = float(P.q) + float(P.delta) + float(P.tau_H)
    P.property_tax_lump_sum_transfer = float(transfer)
    P.birth_entry_grant = bool(grant)
    P.birth_entry_grant_amount = 0.4 if grant else 0.0
    P.birth_entry_grant_locations = np.array([], dtype=int)
    eligible = np.flatnonzero(np.asarray(P.H_own, dtype=float) >= 6.0) + 1
    P.birth_entry_grant_owner_rungs = eligible.astype(int) if grant else np.array([], dtype=int)
    P.propagate_birth_entry_grant = True


def fiscal_ledger(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
    shared: SimpleNamespace,
) -> dict[str, float]:
    model = calendar.model
    revenue = float(
        model.property_tax_revenue_from_distribution(
            evaluation.g_current,
            evaluation.policy.hR_pol,
            evaluation.policy.price,
            P,
        )
    )
    recipient_mass, grant_outlays = model.markov_grant_outlays(
        evaluation.g_post_fertility,
        evaluation.policy.tenure_choice,
        evaluation.policy.tenure_probs,
        P,
        shared,
    )
    transfer_outlays = float(P.property_tax_lump_sum_transfer) * float(
        np.sum(evaluation.g_current)
    )
    current_mass = float(np.sum(evaluation.g_current))
    residual = revenue - float(grant_outlays) - transfer_outlays
    scale = max(abs(revenue), abs(float(grant_outlays)) + abs(transfer_outlays), 1.0e-12)
    implied_transfer = (revenue - float(grant_outlays)) / max(current_mass, 1.0e-12)
    return {
        "property_tax_revenue": revenue,
        "purchase_grant_recipient_mass": float(recipient_mass),
        "purchase_grant_outlays": float(grant_outlays),
        "lump_sum_transfer_outlays": transfer_outlays,
        "government_budget_residual": residual,
        "scaled_government_budget_residual": residual / scale,
        "balanced_budget_implied_transfer": implied_transfer,
        "balanced_budget_transfer_gap": float(P.property_tax_lump_sum_transfer)
        - implied_transfer,
    }


def update_namespace(target: SimpleNamespace, source: SimpleNamespace) -> None:
    """Make the caller's shared-data object match the policy used in the solve."""
    target.__dict__.clear()
    target.__dict__.update(source.__dict__)


def solve_joint_funded_period(
    *,
    g_pre: np.ndarray,
    price_guess: float,
    transfer_guess: float,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    counter: calendar.SolveCounter,
    supply_rule: calendar.HousingSupplyRule,
    annual_tax: float,
    grant: bool,
    tolerance: float,
) -> SimpleNamespace:
    """Jointly solve log house price and the balanced-budget transfer.

    Each residual evaluation is one Bellman solve at a fixed candidate price.
    This is the transient-distribution analogue of the established stationary
    funded-policy root and avoids nesting a complete housing-market bisection
    inside a fiscal bisection.
    """
    cache: dict[tuple[float, float], SimpleNamespace] = {}

    def evaluate(point: np.ndarray) -> SimpleNamespace:
        log_price = float(point[0])
        transfer = float(point[1])
        price = float(math.exp(log_price))
        key = (float(round(log_price, 10)), float(round(transfer, 10)))
        if key not in cache:
            set_fiscal_policy(P, annual_tax=annual_tax, grant=grant, transfer=transfer)
            candidate_shared = calendar.model.precompute_shared(P, b_grid)
            candidate = calendar.evaluate_period(
                np.array([price], dtype=float),
                g_pre,
                P,
                b_grid,
                candidate_shared,
                counter,
                supply_rule=supply_rule,
            )
            ledger = fiscal_ledger(candidate, P, candidate_shared)
            supply = float(candidate.supply_by_loc[0])
            signed_market = float(candidate.demand_by_loc[0] - supply) / max(
                abs(supply), 1.0e-14
            )
            fiscal_scale = max(
                abs(float(ledger["property_tax_revenue"])),
                abs(float(ledger["purchase_grant_outlays"]))
                + abs(float(ledger["lump_sum_transfer_outlays"])),
                0.1,
            )
            cache[key] = SimpleNamespace(
                point=np.array([log_price, transfer], dtype=float),
                price=price,
                transfer=transfer,
                evaluation=candidate,
                shared=candidate_shared,
                ledger=ledger,
                residual=np.array(
                    [
                        signed_market,
                        float(ledger["government_budget_residual"]) / fiscal_scale,
                    ],
                    dtype=float,
                ),
            )
        return cache[key]

    point = np.array(
        [math.log(float(price_guess)), max(float(transfer_guess), 0.0)], dtype=float
    )
    result = evaluate(point)
    steps = np.array([0.01, 0.01], dtype=float)
    jacobian = np.column_stack(
        [
            (evaluate(point + np.array([steps[0], 0.0])).residual - result.residual)
            / steps[0],
            (evaluate(point + np.array([0.0, steps[1]])).residual - result.residual)
            / steps[1],
        ]
    )
    iterations = 0
    for iteration in range(16):
        if not np.all(np.isfinite(result.residual)):
            raise RuntimeError(
                "Joint funded-policy root produced nonfinite residuals: "
                f"{result.residual.tolist()}"
            )
        if float(np.max(np.abs(result.residual))) <= tolerance:
            break
        try:
            direction = np.linalg.solve(jacobian, -result.residual)
        except np.linalg.LinAlgError:
            direction = np.linalg.lstsq(jacobian, -result.residual, rcond=None)[0]
        direction = np.clip(direction, [-0.25, -0.15], [0.25, 0.15])
        current_norm = float(np.linalg.norm(result.residual))
        accepted = None
        accepted_point = None
        for damping in (1.0, 0.5, 0.25, 0.125):
            trial_point = point + damping * direction
            trial_point[0] = float(
                np.clip(
                    trial_point[0],
                    math.log(max(float(getattr(P, "p_min", 1.0e-4)), 1.0e-8)),
                    math.log(float(getattr(P, "p_max", 100.0))),
                )
            )
            trial_point[1] = float(np.clip(trial_point[1], 0.0, 8.0))
            trial = evaluate(trial_point)
            if float(np.linalg.norm(trial.residual)) < current_norm:
                accepted = trial
                accepted_point = trial_point
                break
        if accepted is None or accepted_point is None:
            steps *= 0.5
            jacobian = np.column_stack(
                [
                    (
                        evaluate(point + np.array([steps[0], 0.0])).residual
                        - result.residual
                    )
                    / steps[0],
                    (
                        evaluate(point + np.array([0.0, steps[1]])).residual
                        - result.residual
                    )
                    / steps[1],
                ]
            )
            continue
        displacement = accepted_point - point
        residual_change = accepted.residual - result.residual
        denominator = float(displacement @ displacement)
        if denominator > 1.0e-14:
            jacobian += np.outer(
                residual_change - jacobian @ displacement,
                displacement,
            ) / denominator
        point = accepted_point
        result = accepted
        iterations = iteration + 1
    if not np.all(np.isfinite(result.residual)):
        raise RuntimeError(
            "Joint funded-policy root ended with nonfinite residuals: "
            f"{result.residual.tolist()}"
        )
    if float(np.max(np.abs(result.residual))) > tolerance:
        raise RuntimeError(
            "Housing market did not clear jointly with the fiscal budget: "
            f"residuals={result.residual.tolist()}, evaluations={len(cache)}"
        )
    if not all(math.isfinite(float(value)) for value in result.ledger.values()):
        raise RuntimeError(f"Joint funded-policy root produced a nonfinite ledger: {result.ledger}")
    if abs(float(result.ledger["government_budget_residual"])) > FISCAL_ABSOLUTE_TOLERANCE:
        raise RuntimeError(
            "Joint root passed scaled residuals but failed the absolute fiscal gate: "
            f"{result.ledger}"
        )
    result.joint_iterations = iterations
    result.joint_model_evaluations = len(cache)
    return result


@dataclass
class FundedMarketClearer:
    """Drop-in calendar market clearer with an outer fiscal-transfer root."""

    original_clear: Callable[..., calendar.PeriodEvaluation]
    case: str
    policy_start_period: int
    prepolicy_annual_tax: float
    fiscal_tolerance: float = FISCAL_ABSOLUTE_TOLERANCE

    def __post_init__(self) -> None:
        if self.case not in POLICY_CASES:
            raise ValueError(f"Unknown funded policy case: {self.case}")
        self.period = 0
        self.records: list[dict[str, Any]] = []
        self.last_transfer = 0.0

    def __call__(
        self,
        g_pre: np.ndarray,
        price_guess: float,
        P: SimpleNamespace,
        b_grid: np.ndarray,
        shared: SimpleNamespace,
        counter: calendar.SolveCounter,
        market_tol: float,
        max_iter: int,
        supply_rule: calendar.HousingSupplyRule,
        initial_policy: calendar.PolicyBundle | None = None,
    ) -> calendar.PeriodEvaluation:
        period = int(self.period)
        case_tax, case_grant, case_is_funded = POLICY_CASES[self.case]
        active = bool(case_is_funded and period >= int(self.policy_start_period))
        annual_tax, grant = (
            (case_tax, case_grant)
            if active
            else (float(self.prepolicy_annual_tax), False)
        )
        if not active:
            set_fiscal_policy(P, annual_tax=annual_tax, grant=False, transfer=0.0)
            fresh = calendar.model.precompute_shared(P, b_grid)
            evaluation = self.original_clear(
                g_pre,
                price_guess,
                P,
                b_grid,
                fresh,
                counter,
                market_tol,
                max_iter,
                supply_rule,
                initial_policy=initial_policy,
            )
            ledger = fiscal_ledger(evaluation, P, fresh)
            update_namespace(shared, fresh)
            record = {
                "period": period,
                "policy_case": self.case,
                "policy_active": False,
                "fiscal_convention": "status_quo_unrebated_property_tax",
                "annual_property_tax_rate": annual_tax,
                "lump_sum_transfer_period_units": 0.0,
                "purchase_grant": False,
                "fiscal_root_evaluations": 0,
                **ledger,
            }
        else:
            # In the no-grant cases, market clearing implies the exact budget
            # identity T=tau_H*p*H^S(p)/N.  Use it as the first transfer guess;
            # the joint root below imposes the identity exactly and also covers
            # the grant case, whose outlays depend on endogenous purchases.
            if self.last_transfer > 0.0:
                transfer_guess = self.last_transfer
            else:
                approximate_revenue = (
                    float(annual_tax)
                    * float(P.period_years)
                    * float(price_guess)
                    * float(supply_rule.quantity(np.array([price_guess]))[0])
                )
                transfer_guess = approximate_revenue / max(float(np.sum(g_pre)), 1.0e-12)
            solved = solve_joint_funded_period(
                g_pre=g_pre,
                price_guess=price_guess,
                transfer_guess=transfer_guess,
                P=P,
                b_grid=b_grid,
                counter=counter,
                supply_rule=supply_rule,
                annual_tax=annual_tax,
                grant=grant,
                tolerance=float(market_tol),
            )
            transfer = float(solved.transfer)
            evaluation = solved.evaluation
            ledger = solved.ledger
            set_fiscal_policy(P, annual_tax=annual_tax, grant=grant, transfer=transfer)
            update_namespace(shared, solved.shared)
            if abs(float(ledger["government_budget_residual"])) > self.fiscal_tolerance:
                raise RuntimeError(
                    f"{self.case} t={period}: fiscal gate failed: {ledger}"
                )
            if (
                abs(float(ledger["balanced_budget_transfer_gap"]))
                * float(np.sum(evaluation.g_current))
                > self.fiscal_tolerance
            ):
                raise RuntimeError(
                    f"{self.case} t={period}: transfer-identity gate failed: {ledger}"
                )
            self.last_transfer = transfer
            record = {
                "period": period,
                "policy_case": self.case,
                "policy_active": True,
                "fiscal_convention": "balanced_budget_equal_rebate",
                "annual_property_tax_rate": annual_tax,
                "lump_sum_transfer_period_units": transfer,
                "purchase_grant": bool(grant),
                "fiscal_root_evaluations": int(solved.joint_model_evaluations),
                "joint_root_iterations": int(solved.joint_iterations),
                **ledger,
            }
        if not math.isfinite(float(evaluation.relative_market_residual)):
            raise RuntimeError(f"{self.case} t={period}: nonfinite market residual")
        if evaluation.relative_market_residual > float(market_tol):
            raise RuntimeError(
                f"{self.case} t={period}: market gate failed: "
                f"{evaluation.relative_market_residual:.3e}"
            )
        self.records.append(record)
        self.period += 1
        return evaluation


def reconstruct_transition_inputs(
    contract: dict[str, Any],
    *,
    outside_share: float,
) -> dict[str, Any]:
    chain, model = transition.configure_sequential_model()
    theta = dict(contract["theta"])
    active_domain, profile_overrides, active_profile = activate_contract_profile(
        contract, theta
    )
    base = closure.make_overrides(chain, theta, nb=120, profile="e5f-floor")
    base.update(profile_overrides)
    base.update(theta)
    base.update(
        property_tax_lump_sum_transfer=0.0,
        birth_entry_grant=False,
        birth_entry_grant_amount=0.0,
    )
    old_solution, old_parameters, old_price, _, old_normalization = (
        calibration.solve_old_steady_state(
            chain,
            base,
            initial_psi=float(contract["old_psi_child"]),
            completed_fertility_target=float(contract["old_completed_fertility"]),
            completed_fertility_tolerance=5.0e-4,
            normalize=False,
        )
    )
    if int(old_parameters.J) != 17 or int(old_parameters.Nb) != 120:
        raise RuntimeError("Production dimension gate failed")
    if str(active_profile["name"]) == calibration.REPAIRED_MODEL_PROFILE:
        declared_profile = dict(contract["model_profile"]["report_metadata"])
        if not bool(getattr(old_parameters, "permanent_income_levels_enabled", False)):
            raise RuntimeError("Solved parameters dropped the measured permanent-income process")
        if int(old_parameters.Nz) != int(declared_profile["income_state_count"]):
            raise RuntimeError(
                "Solved income-state count does not match the selected report profile"
            )
        if int(np.asarray(old_parameters.z_grid).size) != int(old_parameters.Nz):
            raise RuntimeError("Solved income grid and Nz disagree")
        if not math.isclose(
            float(old_parameters.permanent_income_log_variance),
            float(declared_profile["permanent_income_log_variance"]),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ):
            raise RuntimeError(
                "Solved permanent-income variance does not match the selected report profile"
            )
        if not math.isclose(
            float(getattr(old_parameters, "first_birth_fixed_cost", math.nan)),
            float(theta["first_birth_fixed_cost"]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise RuntimeError("Solved first-birth cost does not match the report estimate")
        declared_gates = contract.get("income_profile_gates")
        if declared_gates is not None:
            if not bool(declared_gates.get("permanent_income_levels_enabled", False)):
                raise RuntimeError("Report income-profile gate does not activate permanent levels")
            if int(declared_gates.get("income_state_count", -1)) != int(old_parameters.Nz):
                raise RuntimeError("Report income-profile Nz gate disagrees with the solve")
            if float(declared_gates.get("stationary_weight_max_abs_gap", math.inf)) > 1.0e-12:
                raise RuntimeError("Report income-profile stationary-weight gate failed")
    old_moments = chain.extract_moments(old_solution, old_parameters)
    old_fertility_gap = abs(
        float(old_moments["tfr"]) - float(contract["old_completed_fertility"])
    )
    if old_fertility_gap > 5.0e-4:
        raise RuntimeError(
            f"Old steady-state fertility gate failed: gap={old_fertility_gap:.3e}"
        )
    b_grid = np.asarray(old_solution.b_grid, dtype=float)
    shared = model.precompute_shared(old_parameters, b_grid)
    old_parameters._fert2_probs = np.asarray(old_solution.fert2_probs, dtype=float).copy()
    old_policy = calendar.policy_from_solution(
        old_solution, old_price, old_parameters, b_grid, shared
    )
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    calendar.distribution_rows = transition.independent_child_distribution_rows
    stationary_g_pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(
        old_solution, old_policy, old_parameters, b_grid, shared
    )
    gates = transition.operator_gates(
        old_solution,
        old_policy,
        stationary_g_pre,
        old_parameters,
        b_grid,
        shared,
    )
    gates.update(reconstruction)
    nesting_names = (
        "stationary_post_fertility_nesting_l1",
        "one_step_constant_path_nesting_l1",
        "mature_flow_abs_error",
        "birth_flow_abs_error",
        "topcode_adjusted_birth_flow_abs_error",
    )
    if max(float(gates[name]) for name in nesting_names) > 5.0e-9:
        raise RuntimeError(f"Stationary operator gates failed: {gates}")

    renewal = closure.topcode_consistent_renewal_accounting(old_solution, old_parameters)
    entry = float(old_solution.entry_rate)
    raw_mature = float(old_solution.entrants_mature_total)
    adjusted_mature = float(renewal["topcode_adjusted_mature_entrant_households"])
    raw_births = float(old_solution.total_births_kfe)
    adjusted_births = float(renewal["topcode_adjusted_birth_children"])
    retention = (1.0 - float(outside_share)) * entry / adjusted_mature
    if not 0.0 <= retention <= 1.0:
        raise RuntimeError(f"Outside share implies invalid retention: {retention}")

    model_ages = float(old_parameters.age_start) + np.arange(int(old_parameters.J)) * float(
        old_parameters.da
    )
    stationary_age_mass = np.sum(stationary_g_pre, axis=(0, 1, 2, 4, 5, 6))
    age_reweight = transition.acs_2007_age_reweight_diagnostic(
        stationary_age_mass,
        model_ages,
        entry,
        periods=calibration.TRANSITION_PERIODS,
        period_years=float(old_parameters.period_years),
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
    prepolicy_annual_tax = float(old_parameters.tau_H) / float(old_parameters.period_years)
    if not math.isclose(prepolicy_annual_tax, 0.01, rel_tol=0.0, abs_tol=1.0e-12):
        raise RuntimeError(
            f"Expected the calibrated pre-policy 1 percent tax; found {prepolicy_annual_tax}"
        )
    return {
        "chain": chain,
        "base_overrides": base,
        "parameters": old_parameters,
        "b_grid": b_grid,
        "initial_g_pre": initial_g_pre,
        "baseline_price": float(np.asarray(old_price).reshape(-1)[0]),
        "baseline_B": adjusted_mature,
        "baseline_raw_B": raw_mature,
        "baseline_births": adjusted_births,
        "baseline_raw_births": raw_births,
        "entry": entry,
        "outside_flow": float(outside_share) * entry,
        "retention": retention,
        "conversion": float(old_parameters.entrant_conversion_factor),
        "maturation_survival_yield": adjusted_mature
        / (float(old_parameters.entrant_conversion_factor) * adjusted_births),
        "raw_maturation_survival_yield": raw_mature
        / (float(old_parameters.entrant_conversion_factor) * raw_births),
        "supply_rule": supply_rule,
        "supply_normalization": supply_normalization,
        "operator_gates": gates,
        "old_fertility_normalization": old_normalization,
        "prepolicy_annual_tax": prepolicy_annual_tax,
        "active_model_profile": active_profile,
        "active_profile_domain": [list(row) for row in active_domain],
        "income_state_count": int(np.asarray(old_parameters.z_grid).size),
        "first_birth_fixed_cost": float(
            getattr(old_parameters, "first_birth_fixed_cost", 0.0)
        ),
    }


def run_case(
    case: str,
    inputs: dict[str, Any],
    contract: dict[str, Any],
    *,
    periods: int,
    market_tol: float,
    market_max_iter: int,
    outdir: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    case_dir = outdir / "cases" / case
    case_dir.mkdir(parents=True, exist_ok=True)
    counter = calendar.SolveCounter()
    original_clear = calendar.clear_scalar_housing_market
    funded_clearer = FundedMarketClearer(
        original_clear=original_clear,
        case=case,
        policy_start_period=POLICY_START_PERIOD,
        prepolicy_annual_tax=float(inputs["prepolicy_annual_tax"]),
    )
    try:
        calendar.clear_scalar_housing_market = funded_clearer
        path, ages, children, scenario = transition.run_birth_vintage_scenario(
            label=case,
            initial_g_pre=np.asarray(inputs["initial_g_pre"], dtype=float),
            baseline_price=float(inputs["baseline_price"]),
            baseline_B=float(inputs["baseline_B"]),
            baseline_births=float(inputs["baseline_births"]),
            baseline_raw_B=float(inputs["baseline_raw_B"]),
            baseline_raw_births=float(inputs["baseline_raw_births"]),
            parameters=copy.deepcopy(inputs["parameters"]),
            b_grid=np.asarray(inputs["b_grid"], dtype=float),
            periods=int(periods),
            outside_flow=float(inputs["outside_flow"]),
            retention=float(inputs["retention"]),
            conversion=float(inputs["conversion"]),
            delay_periods=4,
            initial_birth_pipeline_multiplier=1.0,
            population_target_indices=transition.census_household_target_indices(),
            population_age_target_years=transition.census_household_age_target_years(),
            policy_case="none",
            policy_start_period=POLICY_START_PERIOD,
            old_psi_child=float(contract["old_psi_child"]),
            new_psi_child=float(contract["new_psi_child"]),
            preference_transition_periods=calibration.TRANSITION_PERIODS,
            supply_rule=inputs["supply_rule"],
            market_tol=float(market_tol),
            market_max_iter=int(market_max_iter),
            outdir=case_dir,
            counter=counter,
        )
    finally:
        calendar.clear_scalar_housing_market = original_clear
    if len(funded_clearer.records) != len(path):
        raise RuntimeError(
            f"Fiscal ledger is incomplete for {case}: "
            f"records={len(funded_clearer.records)}, periods={len(path)}"
        )
    for row, fiscal in zip(path, funded_clearer.records):
        if int(row["period"]) != int(fiscal["period"]):
            raise RuntimeError(f"Fiscal/path period mismatch for {case}")
        row.update(fiscal)
        row["year"] = calibration.START_YEAR + int(row["years_from_start"])
    active = [row for row in path if bool(row["policy_active"])]
    max_fiscal = max(
        (abs(float(row["government_budget_residual"])) for row in active),
        default=0.0,
    )
    max_scaled_fiscal = max(
        (abs(float(row["scaled_government_budget_residual"])) for row in active),
        default=0.0,
    )
    max_mass = float(scenario["max_abs_mass_accounting_residual"])
    max_market = float(scenario["max_market_residual"])
    finite_fields = (
        "asset_price",
        "asset_price_index",
        "housing_user_cost",
        "adult_population",
        "population_index",
        "birth_children_topcode_adjusted",
        "topcode_adjusted_births_per_adult",
        "entry_flow_E",
        "housing_demand",
        "housing_supply",
        "relative_market_residual",
        "owner_rate",
        "dependent_child_owner_rate",
        "mass_accounting_residual",
        "min_distribution_mass",
        "lump_sum_transfer_period_units",
        "property_tax_revenue",
        "purchase_grant_outlays",
        "government_budget_residual",
        "scaled_government_budget_residual",
    )
    nonfinite = sum(
        int(not math.isfinite(float(row[key])))
        for row in path
        for key in finite_fields
    )
    maximum_distribution_nonfinite = max(
        int(row["nonfinite_distribution_count"]) for row in path
    )
    minimum_distribution_mass = min(float(row["min_distribution_mass"]) for row in path)
    maximum_feasibility_projection = max(
        float(row["feasibility_frontier_projection_mass"]) for row in path
    )
    aggregate_gates = (max_fiscal, max_scaled_fiscal, max_mass, max_market)
    if not all(math.isfinite(value) for value in aggregate_gates):
        raise RuntimeError(f"{case}: nonfinite aggregate numerical gate {aggregate_gates}")
    if max_fiscal > FISCAL_ABSOLUTE_TOLERANCE:
        raise RuntimeError(f"{case}: maximum fiscal residual {max_fiscal:.3e}")
    if max_mass > MASS_TOLERANCE:
        raise RuntimeError(f"{case}: maximum mass residual {max_mass:.3e}")
    if max_market > float(market_tol):
        raise RuntimeError(f"{case}: maximum market residual {max_market:.3e}")
    if nonfinite:
        raise RuntimeError(f"{case}: found {nonfinite} nonfinite reported objects")
    if maximum_distribution_nonfinite:
        raise RuntimeError(
            f"{case}: distribution contains {maximum_distribution_nonfinite} nonfinite cells"
        )
    if minimum_distribution_mass < -1.0e-14:
        raise RuntimeError(
            f"{case}: negative distribution mass {minimum_distribution_mass:.3e}"
        )
    if maximum_feasibility_projection > 1.0e-6:
        raise RuntimeError(
            f"{case}: excessive feasibility projection {maximum_feasibility_projection:.3e}"
        )
    summary = {
        "case": case,
        "status": "complete",
        "policy_start_period": POLICY_START_PERIOD,
        "policy_start_year": POLICY_START_YEAR,
        "annual_property_tax_rate": POLICY_CASES[case][0],
        "purchase_grant": POLICY_CASES[case][1],
        "funded_reform": POLICY_CASES[case][2],
        "periods": len(path),
        "endpoint_year": int(path[-1]["year"]),
        "endpoint_asset_price_index": float(path[-1]["asset_price_index"]),
        "endpoint_population_index": float(path[-1]["population_index"]),
        "endpoint_births_per_adult": float(path[-1]["topcode_adjusted_births_per_adult"]),
        "endpoint_owner_rate": float(path[-1]["owner_rate"]),
        "endpoint_dependent_child_owner_rate": float(path[-1]["dependent_child_owner_rate"]),
        "maximum_market_residual": max_market,
        "maximum_mass_residual": max_mass,
        "maximum_fiscal_residual": max_fiscal,
        "maximum_scaled_fiscal_residual": max_scaled_fiscal,
        "nonfinite_reported_objects": nonfinite,
        "maximum_nonfinite_distribution_count": maximum_distribution_nonfinite,
        "minimum_distribution_mass": minimum_distribution_mass,
        "maximum_feasibility_frontier_projection_mass": maximum_feasibility_projection,
        "model_solve_count": int(counter.total),
    }
    calendar.write_csv(case_dir / "transition_path.csv", path)
    calendar.write_csv(case_dir / "adult_age_distribution.csv", ages)
    calendar.write_csv(case_dir / "child_pipeline.csv", children)
    write_json(case_dir / "scenario_summary.json", {**scenario, **summary})
    return path, summary


def prepolicy_identity_gap(paths: dict[str, list[dict[str, Any]]]) -> float:
    baseline = paths[BASELINE_CASE]
    fields = (
        "asset_price",
        "adult_population",
        "birth_children_topcode_adjusted",
        "entry_flow_E",
        "owner_rate",
        "dependent_child_owner_rate",
    )
    maximum = 0.0
    for case, path in paths.items():
        if case == BASELINE_CASE:
            continue
        for period in range(POLICY_START_PERIOD):
            for field in fields:
                maximum = max(
                    maximum,
                    abs(float(path[period][field]) - float(baseline[period][field])),
                )
    return maximum


def matched_2023_and_policy_timing_gates(
    paths: dict[str, list[dict[str, Any]]],
    contract: dict[str, Any],
) -> dict[str, float]:
    baseline = paths[BASELINE_CASE]
    matched = baseline[calibration.TRANSITION_PERIODS]
    population_gap = abs(
        float(matched["population_index"])
        - float(contract["matched_2023_population_index"])
    )
    housing_gap = abs(
        float(matched["asset_price_index"])
        - float(contract["matched_2023_housing_cost_index"])
    )
    if max(population_gap, housing_gap) > 2.0e-8:
        raise RuntimeError(
            "The status-quo path does not preserve the selected matched 2023 economy: "
            f"population_gap={population_gap:.3e}, housing_gap={housing_gap:.3e}"
        )
    for case, path in paths.items():
        if int(path[calibration.TRANSITION_PERIODS]["year"]) != calibration.TERMINAL_YEAR:
            raise RuntimeError(f"{case}: matched period is not labeled 2023")
        if any(bool(row["policy_active"]) for row in path[:POLICY_START_PERIOD]):
            raise RuntimeError(f"{case}: policy activates on or before the matched 2023 date")
        expected_active = bool(POLICY_CASES[case][2])
        if len(path) > POLICY_START_PERIOD:
            first_reform = path[POLICY_START_PERIOD]
            if int(first_reform["year"]) != POLICY_START_YEAR:
                raise RuntimeError(f"{case}: first post-2023 date is not 2027")
            if bool(first_reform["policy_active"]) != expected_active:
                raise RuntimeError(f"{case}: first-reform-date activation is incorrect")
    return {
        "matched_2023_population_gap": population_gap,
        "matched_2023_housing_cost_gap": housing_gap,
    }


def add_relative_effects(
    summaries: list[dict[str, Any]],
    paths: dict[str, list[dict[str, Any]]],
) -> None:
    baseline_summary = next(row for row in summaries if row["case"] == BASELINE_CASE)
    baseline_path = paths[BASELINE_CASE]
    for row in summaries:
        case_path = paths[str(row["case"])]
        row["endpoint_asset_price_change_percent"] = 100.0 * (
            float(row["endpoint_asset_price_index"])
            / float(baseline_summary["endpoint_asset_price_index"])
            - 1.0
        )
        row["endpoint_population_change_percent"] = 100.0 * (
            float(row["endpoint_population_index"])
            / float(baseline_summary["endpoint_population_index"])
            - 1.0
        )
        row["endpoint_birth_rate_change_percent"] = 100.0 * (
            float(row["endpoint_births_per_adult"])
            / float(baseline_summary["endpoint_births_per_adult"])
            - 1.0
        )
        row["endpoint_owner_rate_change_pp"] = 100.0 * (
            float(row["endpoint_owner_rate"])
            - float(baseline_summary["endpoint_owner_rate"])
        )
        row["endpoint_dependent_child_owner_rate_change_pp"] = 100.0 * (
            float(row["endpoint_dependent_child_owner_rate"])
            - float(baseline_summary["endpoint_dependent_child_owner_rate"])
        )
        maximum_path_fiscal = max(
            (
                abs(float(item["government_budget_residual"]))
                for item in case_path
                if bool(item["policy_active"])
            ),
            default=0.0,
        )
        row["maximum_fiscal_residual"] = maximum_path_fiscal
        row["baseline_endpoint_year"] = int(baseline_path[-1]["year"])


def solve_status_quo_stationary_endpoint(
    inputs: dict[str, Any],
    contract: dict[str, Any],
    outdir: Path,
) -> dict[str, Any]:
    """Solve and gate the no-reform open-population steady state."""
    if str(inputs["supply_rule"].mode) != "static-elastic":
        raise RuntimeError("The stationary endpoint requires the dated static-elastic supply rule")
    maturation_yield = float(inputs["maturation_survival_yield"])
    raw_maturation_yield = float(inputs["raw_maturation_survival_yield"])
    if not math.isfinite(maturation_yield) or not 0.0 < maturation_yield <= 1.0 + 1.0e-10:
        raise RuntimeError(f"Invalid top-code-consistent maturation yield: {maturation_yield}")
    if not math.isclose(
        maturation_yield,
        raw_maturation_yield,
        rel_tol=0.0,
        abs_tol=2.0e-12,
    ):
        raise RuntimeError(
            "Top-code-adjusted and raw maturation yields disagree: "
            f"{maturation_yield} versus {raw_maturation_yield}"
        )
    endpoint = transition.solve_stationary_open_endpoint(
        chain=inputs["chain"],
        base_overrides=inputs["base_overrides"],
        old_price=float(inputs["baseline_price"]),
        new_psi_child=float(contract["new_psi_child"]),
        outside_flow=float(inputs["outside_flow"]),
        retention=float(inputs["retention"]),
        conversion=float(inputs["conversion"]),
        maturation_survival_yield=maturation_yield,
        supply_rule=inputs["supply_rule"],
        policy_case="none",
        outdir=outdir / "stationary_endpoint_status_quo",
    )
    if endpoint.get("status") != "complete":
        raise RuntimeError(f"No stationary status-quo endpoint: {endpoint}")
    calendar.write_csv(
        outdir
        / "stationary_endpoint_status_quo"
        / "stationary_pure_branch_endpoint.csv",
        [dict(endpoint)],
    )
    if abs(float(endpoint.get("fixed_stock_relative_market_gap", math.inf))) > 2.5e-5:
        endpoint = refine_stationary_jump_persistent_types(
            endpoint,
            outdir / "stationary_endpoint_status_quo" / "stationary_endpoint_search.csv",
            inputs,
        )
    else:
        endpoint["market_clearing_method"] = "pure-price root"
    finite_names = (
        "asset_price",
        "stationary_population_scale",
        "renewal_denominator",
        "fixed_stock_relative_market_gap",
        "psi_child",
    )
    if not all(math.isfinite(float(endpoint[name])) for name in finite_names):
        raise RuntimeError(f"Nonfinite stationary status-quo endpoint: {endpoint}")
    if float(endpoint["renewal_denominator"]) <= 0.0:
        raise RuntimeError(f"Nonpositive stationary renewal denominator: {endpoint}")
    if float(endpoint["asset_price"]) <= 0.0 or float(endpoint["stationary_population_scale"]) <= 0.0:
        raise RuntimeError(f"Nonpositive stationary price or population: {endpoint}")
    if abs(float(endpoint["fixed_stock_relative_market_gap"])) > 2.5e-5:
        raise RuntimeError(f"Stationary housing-market gate failed: {endpoint}")
    if not math.isclose(
        float(endpoint["psi_child"]),
        float(contract["new_psi_child"]),
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        raise RuntimeError(f"Stationary terminal-preference gate failed: {endpoint}")
    endpoint["maturation_survival_yield"] = maturation_yield
    endpoint["interpretation"] = (
        "No-reform open-population steady state under the terminal fertility "
        "preference, dated static-elastic housing supply, and fixed migration closure. "
        + (
            "At the discrete policy jump, stationarity is conditional on persistent "
            "tie-breaking types, type inheritance along local renewal lines, and the "
            "reported split of outside entrants."
            if endpoint["market_clearing_method"].startswith("persistent entry-type")
            else "The housing market clears at a pure household-policy branch."
        )
    )
    calendar.write_csv(
        outdir / "stationary_endpoint_status_quo" / "stationary_open_endpoint.csv",
        [endpoint],
    )
    return endpoint


def refine_stationary_jump_persistent_types(
    endpoint: dict[str, Any],
    search_path: Path,
    inputs: dict[str, Any],
) -> dict[str, Any]:
    """Clear a policy jump using persistent entry types, not kernel randomization.

    At the indifference price, each type permanently follows one limiting policy
    kernel. Locally renewed entrants retain their type, and outside entrants are
    split across types so each branch's renewal deficit is funded separately.
    Under precisely those conditions, a convex combination of the two invariant
    branch distributions remains stationary. The population-type share is chosen
    jointly with the aggregate renewal scale so the housing market clears.
    """
    with search_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    finite_rows = [
        row
        for row in rows
        if math.isfinite(float(row["asset_price"]))
        and math.isfinite(float(row["fixed_stock_market_gap"]))
    ]
    finite_rows.sort(key=lambda row: float(row["asset_price"]))
    brackets = [
        (left, right)
        for left, right in zip(finite_rows[:-1], finite_rows[1:])
        if float(left["fixed_stock_market_gap"])
        * float(right["fixed_stock_market_gap"])
        < 0.0
    ]
    if not brackets:
        raise RuntimeError(
            "The stationary endpoint failed the market gate and has no audited "
            "discrete-policy jump that brackets zero"
        )
    left, right = min(
        brackets,
        key=lambda pair: float(pair[1]["asset_price"]) - float(pair[0]["asset_price"]),
    )
    price_left = float(left["asset_price"])
    price_right = float(right["asset_price"])
    bracket_width = price_right - price_left
    if bracket_width > 1.0e-6 * max(1.0, abs(price_left), abs(price_right)):
        raise RuntimeError(
            "Stationary sign change is too wide to identify a discrete-policy "
            f"indifference price: width={bracket_width:.3e}"
        )
    positive = left if float(left["fixed_stock_market_gap"]) > 0.0 else right
    negative = right if positive is left else left
    price = 0.5 * (price_left + price_right)
    supply = float(inputs["supply_rule"].quantity(np.array([price], dtype=float))[0])
    outside_flow = float(inputs["outside_flow"])
    retention = float(inputs["retention"])
    conversion = float(inputs["conversion"])
    maturation_yield = float(inputs["maturation_survival_yield"])

    def branch_renewal_denominator(row: dict[str, Any]) -> float:
        entry = float(row["entry_households_per_adult"])
        births = float(row["topcode_adjusted_birth_children_per_adult"])
        mature = conversion * maturation_yield * births
        return entry - retention * mature

    positive_denominator = branch_renewal_denominator(positive)
    negative_denominator = branch_renewal_denominator(negative)
    if not (
        math.isfinite(positive_denominator)
        and math.isfinite(negative_denominator)
        and positive_denominator > 0.0
        and negative_denominator > 0.0
    ):
        raise RuntimeError(
            "Persistent endpoint types require positive branch renewal deficits: "
            f"d_pos={positive_denominator}, d_neg={negative_denominator}"
        )

    def blend(field: str, positive_weight: float) -> float:
        return positive_weight * float(positive[field]) + (1.0 - positive_weight) * float(
            negative[field]
        )

    def mixture(positive_weight: float) -> dict[str, float]:
        entry = blend("entry_households_per_adult", positive_weight)
        adjusted_births = blend(
            "topcode_adjusted_birth_children_per_adult", positive_weight
        )
        queue_mature = conversion * maturation_yield * adjusted_births
        denominator = entry - retention * queue_mature
        if denominator <= 0.0:
            raise RuntimeError(
                "Persistent-type mixture implies a nonpositive renewal denominator"
            )
        population = outside_flow / denominator
        demand_per_adult = blend("housing_demand_per_adult", positive_weight)
        gap = population * demand_per_adult - supply
        return {
            "entry_households_per_adult": entry,
            "topcode_adjusted_birth_children_per_adult": adjusted_births,
            "queue_mature_entrant_flow_B": queue_mature,
            "renewal_denominator": denominator,
            "stationary_population_scale": population,
            "housing_demand_per_adult": demand_per_adult,
            "fixed_stock_market_gap": gap,
            "fixed_stock_relative_market_gap": gap / max(abs(supply), 1.0e-15),
        }

    low_weight = 0.0
    high_weight = 1.0
    low = mixture(low_weight)
    high = mixture(high_weight)
    if low["fixed_stock_market_gap"] >= 0.0 or high["fixed_stock_market_gap"] <= 0.0:
        raise RuntimeError(
            "Recomputed persistent-type branches do not bracket housing-market clearing: "
            f"negative_branch={low}, positive_branch={high}"
        )
    mixed = low
    weight = 0.5
    for _ in range(80):
        weight = 0.5 * (low_weight + high_weight)
        mixed = mixture(weight)
        if abs(mixed["fixed_stock_relative_market_gap"]) <= 1.0e-12:
            break
        if mixed["fixed_stock_market_gap"] > 0.0:
            high_weight = weight
        else:
            low_weight = weight
    if abs(mixed["fixed_stock_relative_market_gap"]) > 1.0e-10:
        raise RuntimeError(f"Persistent-type endpoint did not clear housing: {mixed}")

    population = mixed["stationary_population_scale"]
    positive_population = population * weight
    negative_population = population * (1.0 - weight)
    positive_outside_flow = positive_population * positive_denominator
    negative_outside_flow = negative_population * negative_denominator
    positive_entry_per_adult = float(positive["entry_households_per_adult"])
    negative_entry_per_adult = float(negative["entry_households_per_adult"])
    positive_mature_per_adult = (
        conversion
        * maturation_yield
        * float(positive["topcode_adjusted_birth_children_per_adult"])
    )
    negative_mature_per_adult = (
        conversion
        * maturation_yield
        * float(negative["topcode_adjusted_birth_children_per_adult"])
    )
    positive_branch_flow_residual = (
        positive_population * positive_entry_per_adult
        - retention * positive_population * positive_mature_per_adult
        - positive_outside_flow
    )
    negative_branch_flow_residual = (
        negative_population * negative_entry_per_adult
        - retention * negative_population * negative_mature_per_adult
        - negative_outside_flow
    )
    outside_flow_residual = positive_outside_flow + negative_outside_flow - outside_flow
    if not all(
        math.isfinite(value)
        for value in (
            positive_population,
            negative_population,
            positive_outside_flow,
            negative_outside_flow,
            positive_branch_flow_residual,
            negative_branch_flow_residual,
            outside_flow_residual,
        )
    ):
        raise RuntimeError("Persistent-type endpoint accounting is nonfinite")
    if min(
        positive_population,
        negative_population,
        positive_outside_flow,
        negative_outside_flow,
    ) <= 0.0:
        raise RuntimeError("Persistent endpoint types must both have positive mass and inflow")
    if abs(outside_flow_residual) > 1.0e-12 * max(1.0, abs(outside_flow)):
        raise RuntimeError(
            "Persistent-type outside flows do not reproduce the aggregate closure: "
            f"residual={outside_flow_residual:.3e}"
        )
    if max(
        abs(positive_branch_flow_residual),
        abs(negative_branch_flow_residual),
    ) > 1.0e-12 * max(1.0, abs(outside_flow)):
        raise RuntimeError(
            "Persistent-type branch renewal identities failed: "
            f"positive={positive_branch_flow_residual:.3e}, "
            f"negative={negative_branch_flow_residual:.3e}"
        )
    positive_outside_share = positive_outside_flow / outside_flow
    negative_outside_share = negative_outside_flow / outside_flow
    if not (
        0.0 < positive_outside_share < 1.0
        and 0.0 < negative_outside_share < 1.0
        and abs(positive_outside_share + negative_outside_share - 1.0) <= 1.0e-12
    ):
        raise RuntimeError("Persistent-type implied outside-entrant shares are invalid")

    raw_births = blend("birth_children_per_adult", weight)
    adjusted_births = mixed["topcode_adjusted_birth_children_per_adult"]
    top_bin_birth_flow = blend("top_bin_entry_birth_flow_per_adult", weight)
    extra_positive = float(positive["extra_children_per_top_bin_family"])
    extra_negative = float(negative["extra_children_per_top_bin_family"])
    if not math.isclose(extra_positive, extra_negative, rel_tol=0.0, abs_tol=1.0e-12):
        raise RuntimeError("Persistent branches disagree on the top-bin fertility weight")
    extra_per_top_family = extra_positive
    expected_adjusted_births = raw_births + extra_per_top_family * top_bin_birth_flow
    if not math.isclose(
        adjusted_births,
        expected_adjusted_births,
        rel_tol=0.0,
        abs_tol=2.0e-12,
    ):
        raise RuntimeError(
            "Blended top-code birth accounting is inconsistent: "
            f"reported={adjusted_births}, reconstructed={expected_adjusted_births}"
        )
    accounting_positive = str(positive["renewal_child_accounting_method"])
    accounting_negative = str(negative["renewal_child_accounting_method"])
    if accounting_positive != accounting_negative:
        raise RuntimeError("Persistent branches use different child-accounting methods")
    raw_mature = conversion * maturation_yield * raw_births
    adjusted_mature = conversion * maturation_yield * adjusted_births
    if not math.isclose(
        adjusted_mature,
        mixed["queue_mature_entrant_flow_B"],
        rel_tol=0.0,
        abs_tol=2.0e-12,
    ):
        raise RuntimeError("Blended top-code mature flow is internally inconsistent")
    entry = mixed["entry_households_per_adult"]
    if min(raw_births, adjusted_births, entry, mixed["housing_demand_per_adult"]) <= 0.0:
        raise RuntimeError("Persistent endpoint has a nonpositive displayed flow or demand")
    raw_B_over_E = raw_mature / entry
    adjusted_B_over_E = adjusted_mature / entry
    threeplus_share = top_bin_birth_flow / raw_births
    scale_mapping = supply / mixed["housing_demand_per_adult"]
    if abs(scale_mapping - mixed["stationary_population_scale"]) > 1.0e-10:
        raise RuntimeError(
            "Persistent endpoint Hs/D scale disagrees with the renewal population scale"
        )
    user_cost_rate = float(inputs["parameters"].user_cost_rate)
    if not math.isfinite(user_cost_rate) or user_cost_rate <= 0.0:
        raise RuntimeError("Persistent endpoint has an invalid common user-cost rate")

    result = {
        "status": "complete",
        "label": "postshock_open_endpoint_persistent_types",
        "policy_case": "none",
        "psi_child": float(endpoint["psi_child"]),
        "total_mass": 1.0,
        "asset_price": price,
        "price_ratio": price / float(inputs["baseline_price"]),
        "housing_user_cost": user_cost_rate * price,
        "entry_households_per_adult": entry,
        "birth_children_per_adult": raw_births,
        "top_bin_entry_birth_flow_per_adult": top_bin_birth_flow,
        "topcode_adjusted_birth_children_per_adult": adjusted_births,
        "mature_entrant_households_per_adult": raw_mature,
        "topcode_adjusted_mature_entrant_households_per_adult": adjusted_mature,
        "reproduction_ratio_B_over_E": raw_B_over_E,
        "raw_state_B_over_E": raw_B_over_E,
        "reproduction_residual_B_minus_E": raw_mature - entry,
        "topcode_adjusted_B_over_E_diagnostic": adjusted_B_over_E,
        "topcode_consistent_B_over_E": adjusted_B_over_E,
        "renewal_child_accounting_method": accounting_positive,
        "extra_children_per_top_bin_family": extra_per_top_family,
        "threeplus_share_of_mature_flow": threeplus_share,
        "maturation_exit_yield": maturation_yield,
        "entrant_conversion_factor": conversion,
        "queue_birth_children_raw_explicit_states": raw_births,
        "queue_birth_children_topcode_adjusted": adjusted_births,
        "queue_mature_entrant_flow_B": adjusted_mature,
        "renewal_denominator": mixed["renewal_denominator"],
        "stationary_population_scale": mixed["stationary_population_scale"],
        "housing_supply_mode": str(inputs["supply_rule"].mode),
        "housing_demand_per_adult": mixed["housing_demand_per_adult"],
        "housing_supply": supply,
        "scale_mapping_Hs_over_D": scale_mapping,
        "stationary_housing_supply": supply,
        "stationary_total_housing_demand": (
            mixed["stationary_population_scale"] * mixed["housing_demand_per_adult"]
        ),
        "fixed_stock_market_gap": mixed["fixed_stock_market_gap"],
        "fixed_stock_relative_market_gap": mixed["fixed_stock_relative_market_gap"],
        "market_clearing_method": (
            "persistent entry-type convexification at a discrete household-policy jump"
        ),
        "persistent_type_population_share_positive_branch": weight,
        "persistent_type_population_share_negative_branch": 1.0 - weight,
        "positive_branch_renewal_denominator": positive_denominator,
        "negative_branch_renewal_denominator": negative_denominator,
        "positive_branch_population": positive_population,
        "negative_branch_population": negative_population,
        "positive_branch_outside_entry_flow": positive_outside_flow,
        "negative_branch_outside_entry_flow": negative_outside_flow,
        "positive_branch_renewal_identity_residual": positive_branch_flow_residual,
        "negative_branch_renewal_identity_residual": negative_branch_flow_residual,
        "implied_outside_entrant_share_positive_branch": positive_outside_share,
        "implied_outside_entrant_share_negative_branch": negative_outside_share,
        "outside_entry_flow_identity_residual": outside_flow_residual,
        "persistent_type_stationarity_condition": (
            "type assigned at entry and never changes; locally renewed entrants retain "
            "the source type; outside entrants are split by the reported implied shares"
        ),
        "indifference_price_bracket_lower": price_left,
        "indifference_price_bracket_upper": price_right,
        "indifference_price_bracket_width": bracket_width,
        "pure_branch_positive_relative_market_gap": float(
            positive["fixed_stock_relative_market_gap"]
        ),
        "pure_branch_negative_relative_market_gap": float(
            negative["fixed_stock_relative_market_gap"]
        ),
        "existing_solver_pure_branch_asset_price": float(endpoint["asset_price"]),
        "existing_solver_pure_branch_relative_market_gap": float(
            endpoint["fixed_stock_relative_market_gap"]
        ),
    }
    for field in (
        "own_rate",
        "tfr_topcoded",
        "mean_completed_fertility_raw",
        "childless_rate",
    ):
        if field in positive and field in negative:
            result[field] = blend(field, weight)
    return result


def write_transition_diagnostic(
    paths: dict[str, list[dict[str, Any]]],
    outdir: Path,
) -> dict[str, str]:
    """Write the stable paper-diagnostic 2x2 path figure in PNG and PDF."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = (
        ("asset_price_index", "Housing-cost index", 1.0),
        ("population_index", "Adult-household population index", 1.0),
        ("topcode_adjusted_births_per_adult", "Births per adult", 1.0),
        ("dependent_child_owner_rate", "Dependent-child ownership (%)", 100.0),
    )
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 7.0), sharex=True)
    for axis, (field, title, scale) in zip(axes.flat, panels):
        for case in POLICY_CASES:
            if case not in paths:
                continue
            rows = paths[case]
            axis.plot(
                [int(row["year"]) for row in rows],
                [scale * float(row[field]) for row in rows],
                label=CASE_LABELS[case],
                color=CASE_COLORS[case],
                linestyle=CASE_LINESTYLES[case],
                linewidth=2.0,
            )
        axis.axvline(POLICY_START_YEAR, color="#777777", linewidth=1.0, alpha=0.8)
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
    png = outdir / "funded_policy_transition_paths.png"
    pdf = outdir / "funded_policy_transition_paths.pdf"
    figure.savefig(png, dpi=220, bbox_inches="tight")
    figure.savefig(pdf, bbox_inches="tight")
    plt.close(figure)
    return {"png": str(png), "pdf": str(pdf)}


def main() -> None:
    args = parse_args()
    if int(args.post_2023_periods) < 0:
        raise ValueError("--post-2023-periods must be nonnegative")
    if not 0.0 < float(args.outside_origin_entry_share) < 1.0:
        raise ValueError("--outside-origin-entry-share must lie in (0,1)")
    args.post_2023_periods = effective_post_2023_periods(
        int(args.post_2023_periods), bool(args.smoke)
    )
    source = args.source.resolve()
    report = args.transition_report.resolve()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    cases = selected_cases(args.case)
    started = time.perf_counter()
    actual_code_bundle_sha256 = code_bundle_sha256()
    if (
        args.expected_code_bundle_sha256 is not None
        and actual_code_bundle_sha256 != str(args.expected_code_bundle_sha256)
    ):
        raise RuntimeError(
            "Code-bundle fingerprint gate failed: "
            f"actual={actual_code_bundle_sha256}, "
            f"expected={args.expected_code_bundle_sha256}"
        )
    _, fingerprint_model = transition.configure_sequential_model()
    live_calibration_code_fingerprints = calibration.code_fingerprint_contract(
        fingerprint_model
    )
    contract = load_transition_contract(
        report,
        source,
        expected_report_sha256=args.expected_report_sha256,
        expected_source_sha256=args.expected_source_sha256,
        expected_target_set=args.expected_target_set,
        expected_target_fingerprint=args.expected_target_fingerprint,
        expected_model_profile=args.expected_model_profile,
        outside_origin_entry_share=float(args.outside_origin_entry_share),
        live_calibration_code_fingerprints=live_calibration_code_fingerprints,
    )
    metadata = {
        "status": "running",
        "contract": contract,
        "cases": cases,
        "policy_start": "first model date after the matched 2023 economy (2027)",
        "relative_effect_baseline": BASELINE_CASE,
        "prepolicy_fiscal_convention": "calibrated unrebated 1 percent property tax through 2023",
        "policy_fiscal_rule": (
            "property-tax revenue equals equal per-household transfers plus "
            "targeted purchase-grant outlays at every policy date"
        ),
        "purchase_grant": (
            "0.4 model units for eligible dependent-child renters purchasing "
            "an owner rung with at least six rooms"
        ),
        "housing_price_expectations": "temporary equilibrium/static expectations",
        "outside_origin_entry_share": float(
            contract["provenance"]["outside_origin_entry_share"]
        ),
        "outside_origin_entry_share_status": (
            "selected transition-report contract"
            if contract["provenance"]["outside_origin_entry_share_reported"] is not None
            else "explicit legacy-report routing exception"
        ),
        "transition_report_provenance": contract["provenance"],
        "code_bundle_sha256": actual_code_bundle_sha256,
        "code_bundle_files": list(CODE_BUNDLE_FILES),
        "smoke": bool(args.smoke),
    }
    write_json(outdir / "metadata.json", metadata)
    inputs = reconstruct_transition_inputs(
        contract,
        outside_share=float(contract["provenance"]["outside_origin_entry_share"]),
    )
    total_periods = calibration.TRANSITION_PERIODS + 1 + int(args.post_2023_periods)
    paths: dict[str, list[dict[str, Any]]] = {}
    summaries: list[dict[str, Any]] = []
    for index, case in enumerate(cases, start=1):
        print(f"FUNDED_TRANSITION_CASE {index}/{len(cases)} {case}", flush=True)
        path, summary = run_case(
            case,
            inputs,
            contract,
            periods=total_periods,
            market_tol=float(args.market_tol),
            market_max_iter=int(args.market_max_iter),
            outdir=outdir,
        )
        paths[case] = path
        summaries.append(summary)
        calendar.write_csv(outdir / "summary.csv", summaries)
        write_json(
            outdir / "latest_completed_case.json",
            {
                "status": "running",
                "completed_cases": index,
                "total_cases": len(cases),
                "latest": summary,
            },
        )
    identity_gap = prepolicy_identity_gap(paths)
    if identity_gap > 1.0e-12:
        raise RuntimeError(
            "Policy paths do not share an identical history through the matched "
            f"2023 economy: {identity_gap:.3e}"
        )
    timing_gates = matched_2023_and_policy_timing_gates(paths, contract)
    add_relative_effects(summaries, paths)
    stationary_endpoint = solve_status_quo_stationary_endpoint(inputs, contract, outdir)
    diagnostic_figure = write_transition_diagnostic(paths, outdir)
    combined = [row for case in cases for row in paths[case]]
    calendar.write_csv(outdir / "summary.csv", summaries)
    calendar.write_csv(outdir / "transition_path_long.csv", combined)
    metadata.update(
        status="complete",
        elapsed_seconds=time.perf_counter() - started,
        periods=total_periods,
        endpoint_year=calibration.START_YEAR + 4 * (total_periods - 1),
        prepolicy_path_identity_max_abs_gap=identity_gap,
        matched_2023_and_policy_timing_gates=timing_gates,
        stationary_open_endpoint=stationary_endpoint,
        stationary_open_endpoint_contract={
            "policy": "unrebated 1 percent status quo (no reform)",
            "terminal_psi_child": float(contract["new_psi_child"]),
            "housing_supply_mode": str(inputs["supply_rule"].mode),
            "outside_entry_flow": float(inputs["outside_flow"]),
            "retention": float(inputs["retention"]),
            "entrant_conversion_factor": float(inputs["conversion"]),
            "topcode_consistent_maturation_survival_yield": float(
                inputs["maturation_survival_yield"]
            ),
        },
        diagnostic_figure=diagnostic_figure,
        supply_normalization=inputs["supply_normalization"],
        operator_gates=inputs["operator_gates"],
        old_fertility_normalization=inputs["old_fertility_normalization"],
        active_model_profile=inputs["active_model_profile"],
        active_profile_domain=inputs["active_profile_domain"],
        income_state_count=inputs["income_state_count"],
        first_birth_fixed_cost=inputs["first_birth_fixed_cost"],
        renewal={
            "topcode_adjusted_old_mature_flow": inputs["baseline_B"],
            "raw_old_mature_flow": inputs["baseline_raw_B"],
            "outside_flow": inputs["outside_flow"],
            "retention": inputs["retention"],
            "entrant_conversion_factor": inputs["conversion"],
        },
        numerical_gates={
            "market_tolerance": float(args.market_tol),
            "fiscal_absolute_tolerance": FISCAL_ABSOLUTE_TOLERANCE,
            "mass_tolerance": MASS_TOLERANCE,
            "maximum_market_residual": max(
                float(row["maximum_market_residual"]) for row in summaries
            ),
            "maximum_fiscal_residual": max(
                float(row["maximum_fiscal_residual"]) for row in summaries
            ),
            "maximum_mass_residual": max(
                float(row["maximum_mass_residual"]) for row in summaries
            ),
            "maximum_nonfinite_reported_objects": max(
                int(row["nonfinite_reported_objects"]) for row in summaries
            ),
            "maximum_nonfinite_distribution_count": max(
                int(row["maximum_nonfinite_distribution_count"]) for row in summaries
            ),
            "minimum_distribution_mass": min(
                float(row["minimum_distribution_mass"]) for row in summaries
            ),
            "maximum_feasibility_frontier_projection_mass": max(
                float(row["maximum_feasibility_frontier_projection_mass"])
                for row in summaries
            ),
        },
    )
    write_json(outdir / "metadata.json", metadata)
    write_json(
        outdir / "latest_completed_case.json",
        {"status": "complete", "completed_cases": len(cases), "total_cases": len(cases)},
    )
    print(
        f"FUNDED_TRANSITION_COMPLETE cases={len(cases)} "
        f"periods={total_periods} seconds={metadata['elapsed_seconds']:.2f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
