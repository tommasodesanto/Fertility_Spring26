#!/usr/bin/env python3
"""Build the fixed 2007--2023 E5F model-versus-data validation packet.

The builder performs no model solve and never refits or rescales a model path.
Every input is supplied on the command line.  The only allowed indexing is the
declared 2007=1 normalization already present in the source data, plus the same
transparent normalization of the model's implicit rental-cost flow used solely
to construct a price-to-implicit-cost diagnostic.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.e5_profile import e5_target_system


CALENDAR_YEARS = (2007, 2011, 2015, 2019, 2023)
HASH_PATTERN = re.compile(r"^[0-9a-f]{64}$")
EXPECTED_MODEL_PROFILE = "e5f-income-entry"
EXPECTED_MODEL_PROFILE_ID = "eqscale_seq_e5f_income_entry_20260816"
ACTIVE_TARGET_SYSTEM = e5_target_system()
if ACTIVE_TARGET_SYSTEM.count != 12:
    raise RuntimeError("The live E5 profile must contain exactly twelve targets")
EXPECTED_TARGET_SET = ACTIVE_TARGET_SYSTEM.name
EXPECTED_TARGET_COUNT = 12
EXPECTED_TARGET_FINGERPRINT = ACTIVE_TARGET_SYSTEM.fingerprint
EXPECTED_NCHS_METADATA_SHA256 = (
    "631454d5901399e9ee26ec914d13a5c3cb657b77b428308d04cf65aa644ed60d"
)
EXPECTED_NCHS_CONTRACT_ID = (
    "nchs_first_birth_timing_v4_boundary_collapsed_midpoint_1987_2023_"
    "lbo_unknown_fixed"
)
EXPECTED_NCHS_CACHE_CONTRACT = (
    "nchs_first_birth_counts_v4_lbo_unknown_code_fixed_raw_sha256"
)
EXPECTED_NCHS_SOURCE_BUNDLE_SHA256 = (
    "98bf7bb5339f5ff599a3b0fe135c20784321df3573cefe9d2a890a0a2251cb54"
)
EXPECTED_NCHS_CONTRACT_BUNDLE_SHA256 = (
    "1f52a6280e23e090099f31a07393b9d4954c413f79a8c047d084c2a55911e938"
)
EXPECTED_NCHS_PRIMARY_TARGET = {
    "mean_age_first_birth": 26.0446272574833,
    "share_30plus": 0.260327401666964,
    "n_first_births": 10_578_871,
    "mean_standard_error": 0.15,
    "share_30plus_standard_error": 0.01,
    "mean_inverse_variance_weight": 44.4444444444444,
    "share_30plus_inverse_variance_weight": 10_000.0,
}
EXPECTED_NCHS_DOCUMENT_SHA256 = {
    "official_1987_nchs_documentation_pdf": (
        "145d8496d764ee66583ef1520ae56804e1c3373bea35d3a570ecbc1426bd5c79"
    ),
    "official_1989_nchs_documentation_pdf": (
        "92dab8115baec71eec3633239cbd042b2079ad6b80bd1b3a3a43c3276ac3a7cb"
    ),
    "official_2003_nchs_documentation_pdf": (
        "82246197e30d54c56a69314bfdcb8f553e6ca0a0d509f8ce8112c2e996b5b2f5"
    ),
}
REPLACEMENT_FERTILITY = 2.1
RENEWAL_CONVERSION = 1.0 / REPLACEMENT_FERTILITY
ENTRY_DELAY_YEARS = 20.0
NCHS_AGE_BINS = tuple((lower, lower + 3) for lower in range(18, 43, 4))
SCOPE_WARNING = (
    "The imposed bridge and historical holdouts are national series. The twelve "
    "maintained 2023 targets combine national aggregates with mixed-scope PSID, "
    "NCHS, and other calibration objects; they are not a national time-series fit."
)
ROOM_MOMENT = "aggregate_mean_occupied_rooms_18_85"
TERMINAL_ROLE = "calibration_moment_at_2023_model_state"
DIAGNOSTIC_ROLE = "untargeted_transition_diagnostic"
PERIOD_FERTILITY_ROLE = "untargeted_current-rate_diagnostic"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case-dir",
        type=Path,
        required=True,
        help=(
            "Selected transition case containing transition_path.csv, "
            "dated_model_moments.csv, and dated_period_fertility.csv."
        ),
    )
    parser.add_argument("--acs-housing-csv", type=Path, required=True)
    parser.add_argument("--acs-female-exposure-csv", type=Path, required=True)
    parser.add_argument("--historical-comparison-csv", type=Path, required=True)
    parser.add_argument("--historical-metadata-json", type=Path, required=True)
    parser.add_argument("--nchs-counts-csv", type=Path, required=True)
    parser.add_argument("--nchs-manifest-csv", type=Path, required=True)
    parser.add_argument("--nchs-metadata-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    for row in rows:
        if list(row) != fieldnames:
            raise ValueError(f"Inconsistent output schema for {path}")
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"CSV has no observations: {path}")
    return rows


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def require_columns(
    rows: Sequence[Mapping[str, str]], required: Iterable[str], label: str
) -> None:
    columns = set(rows[0])
    missing = sorted(set(required) - columns)
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def finite_float(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} is not numeric: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"{label} is not finite: {value!r}")
    return result


def optional_float(value: Any, label: str) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    return finite_float(value, label)


def integer(value: Any, label: str) -> int:
    number = finite_float(value, label)
    rounded = int(round(number))
    if abs(number - rounded) > 1e-9:
        raise ValueError(f"{label} is not an integer: {value!r}")
    return rounded


def require_hash(value: Any, label: str) -> str:
    result = str(value or "").lower()
    if not HASH_PATTERN.fullmatch(result):
        raise ValueError(f"{label} is not a SHA-256 fingerprint")
    return result


def canonical_rows_by_year(
    rows: Sequence[Mapping[str, str]],
    *,
    label: str,
    year_column: str = "calendar_year",
    allow_other_years: bool = False,
) -> dict[int, Mapping[str, str]]:
    require_columns(rows, [year_column], label)
    by_year: dict[int, Mapping[str, str]] = {}
    for row in rows:
        year = integer(row[year_column], f"{label}.{year_column}")
        if year in by_year:
            raise ValueError(f"{label} has duplicate calendar year {year}")
        by_year[year] = row
    missing = sorted(set(CALENDAR_YEARS) - set(by_year))
    if missing:
        raise ValueError(f"{label} is missing calendar years {missing}")
    if not allow_other_years and set(by_year) != set(CALENDAR_YEARS):
        extra = sorted(set(by_year) - set(CALENDAR_YEARS))
        raise ValueError(f"{label} has unexpected calendar years {extra}")
    return {year: by_year[year] for year in CALENDAR_YEARS}


def find_case_summary(case_dir: Path) -> Path:
    candidates = [
        case_dir / "summary.json",
        case_dir.parent / "summary.json",
        case_dir.parent.parent / "summary.json",
    ]
    existing = []
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.is_file() and resolved not in existing:
            existing.append(resolved)
    if not existing:
        raise FileNotFoundError(
            f"No summary.json found in or up to two levels above {case_dir}"
        )
    if len(existing) > 1:
        raise ValueError(f"Ambiguous case summary paths: {existing}")
    return existing[0]


def validate_case_summary(summary: Mapping[str, Any], summary_path: Path) -> dict[str, str]:
    required = {
        "status",
        "source_sha256",
        "code_fingerprints",
        "target_fingerprint",
        "target_set",
        "target_measurements",
        "population_bridge",
        "population_validation_status",
        "renewal_accounting_contract",
        "renewal_accounting_old_state",
        "renewal_retention",
        "outside_origin_entry_share",
        "old_model_completed_fertility",
        "old_completed_fertility_reference",
        "old_fertility_normalization",
        "model_profile",
        "target_count",
        "transition_free_parameter_count",
        "policy_case",
        "post_2023_periods",
    }
    missing = sorted(required - set(summary))
    if missing:
        raise ValueError(f"{summary_path} lacks provenance fields: {missing}")
    if not str(summary["status"]).startswith("complete_"):
        raise ValueError(f"Selected model task is not complete: {summary['status']!r}")
    code = summary["code_fingerprints"]
    if not isinstance(code, dict) or not isinstance(code.get("files"), dict):
        raise ValueError("code_fingerprints must contain the per-file fingerprint map")
    fingerprints = {
        "source_sha256": require_hash(summary["source_sha256"], "source_sha256"),
        "target_fingerprint": require_hash(
            summary["target_fingerprint"], "target_fingerprint"
        ),
        "code_bundle_sha256": require_hash(
            code.get("bundle_sha256"), "code_fingerprints.bundle_sha256"
        ),
    }
    for name, value in code["files"].items():
        require_hash(value, f"code_fingerprints.files[{name!r}]")
    if summary["target_set"] != EXPECTED_TARGET_SET:
        raise ValueError(f"Unexpected target set: {summary['target_set']!r}")
    if fingerprints["target_fingerprint"] != EXPECTED_TARGET_FINGERPRINT:
        raise ValueError("Selected task does not use the pinned target fingerprint")
    if integer(summary["target_count"], "target_count") != EXPECTED_TARGET_COUNT:
        raise ValueError("Selected task does not have the pinned twelve-target contract")
    free_count = integer(
        summary["transition_free_parameter_count"], "transition_free_parameter_count"
    )
    if free_count != 10:
        raise ValueError("Selected transition calibration does not have ten free parameters")
    profile = summary["model_profile"]
    if not isinstance(profile, dict):
        raise ValueError("model_profile must be an object")
    if profile.get("name") != EXPECTED_MODEL_PROFILE:
        raise ValueError(f"Historical packet requires {EXPECTED_MODEL_PROFILE}")
    if profile.get("profile_id") != EXPECTED_MODEL_PROFILE_ID:
        raise ValueError("Unexpected e5f-income-entry profile identifier")
    if integer(profile.get("income_state_count"), "model_profile.income_state_count") != 15:
        raise ValueError("e5f-income-entry must have fifteen income states")
    if profile.get("first_birth_fixed_cost_semantics") != (
        "one-time utility cost paid only when the first child arrives"
    ):
        raise ValueError("Unexpected first-birth fixed-cost semantics")
    if summary["policy_case"] != "none" or integer(
        summary["post_2023_periods"], "post_2023_periods"
    ) != 0:
        raise ValueError("Historical validation requires the five-date no-policy path")
    bridge = summary["population_bridge"]
    if not isinstance(bridge, dict) or bridge.get("status") != "externally_normalized_not_estimated":
        raise ValueError("Population bridge is not classified as externally normalized")
    measured = summary["target_measurements"]
    expected_measurements = {
        "tfr",
        "childless_rate",
        "mean_age_first_birth",
        "share_first_births_age30plus",
        "housing_increment_0to1",
        "remaining_targets",
    }
    if not isinstance(measured, dict) or set(measured) != expected_measurements or not all(
        isinstance(value, str) and value.strip() for value in measured.values()
    ):
        raise ValueError("target_measurements does not contain the pinned measurement ledger")
    renewal = summary["renewal_accounting_contract"]
    if not isinstance(renewal, dict):
        raise ValueError("renewal_accounting_contract is incomplete")
    for key in (
        "replacement_fertility",
        "effective_birth_to_household_conversion",
        "birth_to_entry_effect_lag_years",
        "birth_to_entry_effect_lag_dates",
        "birth_vintage_queue_waiting_slots",
        "birth_measure",
    ):
        if key not in renewal:
            raise ValueError(f"renewal_accounting_contract lacks {key}")
    replacement = finite_float(
        renewal["replacement_fertility"], "renewal replacement fertility"
    )
    conversion = finite_float(
        renewal["effective_birth_to_household_conversion"], "renewal conversion"
    )
    delay_years = finite_float(
        renewal["birth_to_entry_effect_lag_years"], "birth-to-entry effect lag years"
    )
    if not math.isclose(replacement, REPLACEMENT_FERTILITY, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("Replacement fertility is not exactly the pinned 2.1 benchmark")
    if not math.isclose(conversion, 1.0 / replacement, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("Birth-to-household conversion is not the inverse of 2.1")
    if not math.isclose(conversion, RENEWAL_CONVERSION, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("Birth-to-household conversion changed from 1/2.1")
    if not math.isclose(delay_years, ENTRY_DELAY_YEARS, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("Birth-to-entry delay is not twenty years")
    if integer(renewal["birth_to_entry_effect_lag_dates"], "effect lag dates") != 5:
        raise ValueError("Twenty-year birth-to-entry lag must span five four-year dates")
    if integer(renewal["birth_vintage_queue_waiting_slots"], "queue waiting slots") != 4:
        raise ValueError("Birth-vintage queue does not use four waiting slots")
    if renewal["birth_measure"] != "topcode_adjusted_birth_children":
        raise ValueError("Population renewal is not using topcode-adjusted births")

    old_reference = finite_float(
        summary["old_completed_fertility_reference"], "old fertility reference"
    )
    old_fertility = finite_float(
        summary["old_model_completed_fertility"], "old model completed fertility"
    )
    normalization = summary["old_fertility_normalization"]
    if not isinstance(normalization, dict) or normalization.get("status") != "derived_intercept":
        raise ValueError("Old fertility intercept was not transparently derived")
    if not math.isclose(old_reference, replacement, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("Old steady-state fertility reference does not equal 2.1")
    if not math.isclose(
        finite_float(normalization.get("target"), "old normalization target"),
        replacement,
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise ValueError("Old fertility normalization target is not 2.1")
    if not math.isclose(
        finite_float(normalization.get("completed_fertility"), "old normalized fertility"),
        old_fertility,
        rel_tol=0.0,
        abs_tol=1e-10,
    ):
        raise ValueError("Old completed-fertility receipts disagree")

    old_state = summary["renewal_accounting_old_state"]
    if not isinstance(old_state, dict):
        raise ValueError("renewal_accounting_old_state must be an object")
    entry = finite_float(old_state.get("old_entry_flow_E"), "old entry flow")
    mature = finite_float(old_state.get("old_queue_mature_flow_B"), "old queue flow")
    ratio = finite_float(old_state.get("old_queue_B_over_E"), "old queue B/E")
    outside = finite_float(old_state.get("outside_flow_M"), "outside flow")
    outside_share = finite_float(
        summary["outside_origin_entry_share"], "outside-origin entry share"
    )
    retention = finite_float(summary["renewal_retention"], "renewal retention")
    if min(entry, mature) <= 0 or not 0 <= outside_share < 1 or not 0 <= retention <= 1:
        raise ValueError("Renewal accounting has an invalid flow/share")
    if not math.isclose(mature / entry, ratio, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("Saved old queue B/E does not equal B divided by E")
    if not math.isclose(
        ratio, old_fertility / replacement, rel_tol=0.0, abs_tol=2e-12
    ):
        raise ValueError("Old queue renewal ratio does not implement fertility divided by 2.1")
    if not math.isclose(outside, outside_share * entry, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("Outside entry flow does not equal its declared share times E")
    if not math.isclose(
        retention, (1.0 - outside_share) * entry / mature, rel_tol=0.0, abs_tol=1e-12
    ):
        raise ValueError("Retention does not implement the declared renewal algebra")
    renewal_residual = entry - outside - retention * mature
    if abs(renewal_residual) > 1e-12 or abs(
        finite_float(old_state.get("old_renewal_residual"), "saved renewal residual")
    ) > 1e-12:
        raise ValueError("Old steady-state renewal identity does not close")
    fingerprints["renewal_contract"] = (
        "replacement=2.1; conversion=1/2.1; effect lag=20 years/5 dates; "
        "queue waiting slots=4; "
        "E=M+retention*B"
    )
    return fingerprints


def validate_selected_case(summary: Mapping[str, Any], case_dir: Path) -> str:
    best = summary.get("best_candidate")
    if not isinstance(best, dict) or not best.get("candidate"):
        raise ValueError("Task summary lacks the selected best_candidate identifier")
    candidate = str(best["candidate"])
    if case_dir.name != candidate:
        raise ValueError(
            f"Case directory {case_dir.name!r} is not the task's selected candidate {candidate!r}"
        )
    return candidate


def load_transition(case_dir: Path) -> tuple[dict[int, dict[str, float]], Path]:
    path = (case_dir / "transition_path.csv").resolve()
    rows = read_csv(path)
    require_columns(
        rows,
        [
            "period",
            "years_from_start",
            "asset_price_index",
            "housing_user_cost",
            "population_index",
            "population_target_index",
            "population_target_gap",
            "owner_rate",
            "relative_market_residual",
            "mass_accounting_residual",
            "census_age_bridge_maximum_target_gap",
            "policy_case",
            "policy_active",
            "closure",
        ],
        "transition_path",
    )
    if len(rows) != len(CALENDAR_YEARS):
        raise ValueError("No-policy historical transition must contain exactly five dates")
    by_year: dict[int, dict[str, float]] = {}
    for row in rows:
        years_from_start = finite_float(
            row["years_from_start"], "transition_path.years_from_start"
        )
        year = int(round(CALENDAR_YEARS[0] + years_from_start))
        if abs(CALENDAR_YEARS[0] + years_from_start - year) > 1e-9:
            raise ValueError("transition_path has a noninteger calendar year")
        if year not in CALENDAR_YEARS:
            raise ValueError(f"transition_path has unexpected calendar year {year}")
        if year in by_year:
            raise ValueError(f"transition_path has duplicate calendar year {year}")
        if row["policy_case"] != "none" or row["policy_active"].strip().lower() not in {
            "false",
            "0",
        }:
            raise ValueError(f"Historical transition activates policy in {year}")
        if row["closure"] != "open_birth_vintage":
            raise ValueError(f"Unexpected population closure in {year}")
        by_year[year] = {
            "period": finite_float(row["period"], "transition_path.period"),
            "asset_price_index": finite_float(
                row["asset_price_index"], "transition_path.asset_price_index"
            ),
            "housing_user_cost": finite_float(
                row["housing_user_cost"], "transition_path.housing_user_cost"
            ),
            "population_index": finite_float(
                row["population_index"], "transition_path.population_index"
            ),
            "population_target_index": finite_float(
                row["population_target_index"],
                "transition_path.population_target_index",
            ),
            "population_target_gap": finite_float(
                row["population_target_gap"], "transition_path.population_target_gap"
            ),
            "owner_rate": finite_float(row["owner_rate"], "transition_path.owner_rate"),
            "relative_market_residual": finite_float(
                row["relative_market_residual"],
                "transition_path.relative_market_residual",
            ),
            "mass_accounting_residual": finite_float(
                row["mass_accounting_residual"],
                "transition_path.mass_accounting_residual",
            ),
            "census_age_bridge_maximum_target_gap": finite_float(
                row["census_age_bridge_maximum_target_gap"],
                "transition_path.census_age_bridge_maximum_target_gap",
            ),
        }
    missing = sorted(set(CALENDAR_YEARS) - set(by_year))
    if missing:
        raise ValueError(f"transition_path is missing calendar years {missing}")
    for index, year in enumerate(CALENDAR_YEARS):
        if abs(by_year[year]["period"] - index) > 1e-9:
            raise ValueError(f"transition_path period/year mismatch at {year}")
    for key in ("asset_price_index", "population_index", "population_target_index"):
        if not math.isclose(by_year[2007][key], 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(f"{key} is not declared on a 2007=1 scale")
    for year in CALENDAR_YEARS:
        row = by_year[year]
        implied_gap = row["population_index"] - row["population_target_index"]
        if not math.isclose(
            implied_gap, row["population_target_gap"], rel_tol=0.0, abs_tol=1e-10
        ):
            raise ValueError(f"Population bridge gap is internally inconsistent in {year}")
    return by_year, path


def load_dated_moments(
    case_dir: Path,
    expected_candidate: str,
) -> tuple[dict[int, dict[str, dict[str, Any]]], Path]:
    path = (case_dir / "dated_model_moments.csv").resolve()
    rows = read_csv(path)
    require_columns(
        rows,
        [
            "period",
            "calendar_year",
            "moment",
            "model",
            "role",
            "target",
            "weight",
            "gap",
            "loss_contribution",
            "candidate",
        ],
        "dated_model_moments",
    )
    by_year: dict[int, dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in rows:
        if row["candidate"] != expected_candidate:
            raise ValueError("dated_model_moments belongs to a different candidate")
        year = integer(row["calendar_year"], "dated_model_moments.calendar_year")
        if year not in CALENDAR_YEARS:
            raise ValueError(f"dated_model_moments has unexpected calendar year {year}")
        period = integer(row["period"], "dated_model_moments.period")
        if period != CALENDAR_YEARS.index(year):
            raise ValueError(f"dated_model_moments period/year mismatch in {year}")
        name = str(row["moment"])
        if not name or name in by_year[year]:
            raise ValueError(f"Duplicate/empty dated moment {name!r} in {year}")
        model_value = finite_float(row["model"], f"dated moment {name}.model")
        role = str(row["role"])
        target = optional_float(row["target"], f"dated moment {name}.target")
        weight = optional_float(row["weight"], f"dated moment {name}.weight")
        gap = optional_float(row["gap"], f"dated moment {name}.gap")
        loss = optional_float(
            row["loss_contribution"], f"dated moment {name}.loss_contribution"
        )
        if year == 2023:
            if role != TERMINAL_ROLE or None in (target, weight, gap, loss):
                raise ValueError(f"Terminal maintained moment {name} lacks its target contract")
            assert target is not None and weight is not None and gap is not None and loss is not None
            if not math.isclose(model_value - target, gap, rel_tol=0.0, abs_tol=1e-10):
                raise ValueError(f"Terminal gap mismatch for {name}")
            if not math.isclose(weight * gap * gap, loss, rel_tol=1e-9, abs_tol=1e-10):
                raise ValueError(f"Terminal loss-contribution mismatch for {name}")
        else:
            if role != DIAGNOSTIC_ROLE or any(
                value is not None for value in (target, weight, gap, loss)
            ):
                raise ValueError(f"Nonterminal moment {name} is incorrectly classified")
        by_year[year][name] = {
            "model": model_value,
            "role": role,
            "target": target,
            "weight": weight,
            "gap": gap,
            "loss_contribution": loss,
        }
    if set(by_year) != set(CALENDAR_YEARS):
        raise ValueError("dated_model_moments does not cover the five canonical dates")
    moment_sets = {tuple(sorted(values)) for values in by_year.values()}
    if len(moment_sets) != 1:
        raise ValueError("dated_model_moments changes its moment set across dates")
    if ROOM_MOMENT not in by_year[2007]:
        raise ValueError(f"dated_model_moments lacks {ROOM_MOMENT}")
    return dict(by_year), path


def validate_target_and_loss_contract(
    summary: Mapping[str, Any],
    moments: Mapping[int, Mapping[str, Mapping[str, Any]]],
) -> dict[str, Any]:
    terminal = moments[2023]
    if len(terminal) != EXPECTED_TARGET_COUNT:
        raise ValueError("Terminal ledger does not contain exactly twelve target rows")
    payload = [
        [name, float(record["target"]), float(record["weight"])]
        for name, record in terminal.items()
    ]
    encoded = json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode()
    fingerprint = hashlib.sha256(encoded).hexdigest()
    if fingerprint != EXPECTED_TARGET_FINGERPRINT:
        raise ValueError("Terminal target rows do not reproduce the pinned fingerprint")
    if fingerprint != summary["target_fingerprint"]:
        raise ValueError("Terminal target rows disagree with the task fingerprint")
    loss = math.fsum(
        float(record["loss_contribution"]) for record in terminal.values()
    )
    best = summary["best_candidate"]
    saved_loss = finite_float(best.get("transition_loss"), "best transition loss")
    if not math.isclose(loss, saved_loss, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError("Terminal target rows do not reproduce the selected loss")
    if best.get("policy_case") != "none" or integer(
        best.get("post_2023_periods"), "best post_2023_periods"
    ) != 0:
        raise ValueError("Selected best candidate is not the no-policy five-date case")
    terminal_links = {
        "tfr": "terminal_tfr",
        "childless_rate": "terminal_childless_rate",
        "mean_age_first_birth": "terminal_mean_age_first_birth",
        "share_first_births_age30plus": "terminal_share_first_births_age30plus",
    }
    for moment_name, summary_name in terminal_links.items():
        if not math.isclose(
            float(terminal[moment_name]["model"]),
            finite_float(best.get(summary_name), f"best_candidate.{summary_name}"),
            rel_tol=0.0,
            abs_tol=1e-10,
        ):
            raise ValueError(f"Selected summary does not reconcile {moment_name}")
    return {
        "target_count": len(terminal),
        "target_fingerprint": fingerprint,
        "recomputed_loss": loss,
        "saved_loss": saved_loss,
    }


def validate_population_bridge_contract(
    summary: Mapping[str, Any],
    transition: Mapping[int, Mapping[str, float]],
    historical: Mapping[int, Mapping[str, float]],
) -> dict[str, Any]:
    """Require the observed household series to be the exact imposed bridge.

    This is an accounting/provenance gate, not a model-fit test.  A historical
    population series that differs from the bridge cannot enter this packet
    under the label "matched by construction."
    """

    bridge = summary["population_bridge"]
    targets = bridge.get("target_indices")
    if not isinstance(targets, dict):
        raise ValueError("Population bridge lacks its five target indices")
    receipts: dict[str, Any] = {}
    for period, year in enumerate(CALENDAR_YEARS):
        declared = finite_float(
            targets.get(str(period)), f"population bridge target period {period}"
        )
        target = float(transition[year]["population_target_index"])
        model = float(transition[year]["population_index"])
        observed = float(historical[year]["observed_population_index"])
        for label, value in (
            ("transition target", target),
            ("transition population", model),
            ("observed household population", observed),
        ):
            if not math.isclose(value, declared, rel_tol=0.0, abs_tol=2e-12):
                raise ValueError(
                    f"{label} in {year} does not equal the imposed population target"
                )
        receipts[str(year)] = {
            "period": period,
            "declared_target_index": declared,
            "observed_population_index": observed,
            "model_population_index": model,
            "maximum_absolute_gap": max(
                abs(target - declared),
                abs(model - declared),
                abs(observed - declared),
            ),
        }
    return {
        "status": "observed household path equals imposed target at every date",
        "classification": "accounting_identity_not_fit",
        "receipts": receipts,
    }


def load_dated_fertility(case_dir: Path) -> tuple[dict[int, dict[str, float]], Path]:
    path = (case_dir / "dated_period_fertility.csv").resolve()
    rows = read_csv(path)
    require_columns(
        rows,
        [
            "period",
            "calendar_year",
            "period_tfr_explicit_states",
            "period_tfr_topcode_adjusted",
            "period_first_birth_mean_age",
            "period_first_birth_share_age30plus",
            "role",
            "denominator",
        ],
        "dated_period_fertility",
    )
    by_year = canonical_rows_by_year(rows, label="dated_period_fertility")
    output: dict[int, dict[str, float]] = {}
    for index, year in enumerate(CALENDAR_YEARS):
        row = by_year[year]
        if integer(row["period"], "dated_period_fertility.period") != index:
            raise ValueError(f"dated_period_fertility period/year mismatch in {year}")
        if row["role"] != PERIOD_FERTILITY_ROLE:
            raise ValueError("Period fertility is not classified as an untargeted diagnostic")
        if "adult-household" not in row["denominator"]:
            raise ValueError("Period fertility denominator is not explicitly documented")
        output[year] = {
            key: finite_float(row[key], f"dated_period_fertility.{key}")
            for key in (
                "period_tfr_explicit_states",
                "period_tfr_topcode_adjusted",
                "period_first_birth_mean_age",
                "period_first_birth_share_age30plus",
            )
        }
    return output, path


def load_acs_housing(
    path: Path,
) -> tuple[dict[int, dict[str, float]], dict[str, Any], Path, Path]:
    path = path.resolve()
    rows = read_csv(path)
    require_columns(
        rows,
        [
            "calendar_year",
            "sample_code",
            "base_head_records",
            "base_head_weight",
            "tenure_valid_records",
            "tenure_valid_weight",
            "owner_records",
            "owner_weight",
            "ownership_rate",
            "rooms_valid_records",
            "rooms_valid_weight",
            "rooms_weighted_sum",
            "mean_rooms_literal",
            "ownership_rate_index_2007",
            "mean_rooms_literal_index_2007",
        ],
        "ACS housing path",
    )
    by_year = canonical_rows_by_year(rows, label="ACS housing path")
    metadata_path = path.parent / "metadata.json"
    verification_path = path.parent / "verification_report.json"
    metadata = read_json(metadata_path)
    verification = read_json(verification_path)
    for key in ("contract_sha256", "source", "builder", "output_files"):
        if key not in metadata:
            raise ValueError(f"ACS metadata lacks {key}")
    require_hash(metadata["contract_sha256"], "ACS contract_sha256")
    if not isinstance(metadata["source"], dict) or not isinstance(metadata["builder"], dict):
        raise ValueError("ACS source and builder receipts must be objects")
    require_hash(metadata["source"].get("sha256"), "ACS source.sha256")
    require_hash(metadata["builder"].get("sha256"), "ACS builder.sha256")
    output_files = metadata["output_files"]
    if not isinstance(output_files, dict):
        raise ValueError("ACS metadata.output_files must be an object")
    receipt = output_files.get(path.name)
    if not isinstance(receipt, dict) or "sha256" not in receipt:
        raise ValueError(f"ACS metadata has no output fingerprint for {path.name}")
    recorded_hash = require_hash(receipt["sha256"], "ACS output CSV sha256")
    if recorded_hash != sha256(path):
        raise ValueError("ACS path fingerprint does not match metadata")
    if verification.get("status") not in {"pass", "passed", "complete", "verified"}:
        raise ValueError("ACS verification report does not have a passing status")
    for key in (
        "canonical_csv_byte_exact",
        "canonical_female_exposure_csv_byte_exact",
        "canonical_metadata_byte_exact",
    ):
        if verification.get(key) is not True:
            raise ValueError(f"ACS verification report failed {key}")
    checks = verification.get("self_checks")
    if not isinstance(checks, list) or not checks or not all(
        isinstance(check, dict) and check.get("passed") is True for check in checks
    ):
        raise ValueError("ACS verification report has a failed or missing self-check")
    output: dict[int, dict[str, float]] = {}
    sample_codes: set[str] = set()
    for year in CALENDAR_YEARS:
        row = by_year[year]
        sample_codes.add(row["sample_code"])
        values = {
            name: finite_float(row[name], f"ACS.{name}")
            for name in (
                "base_head_records",
                "base_head_weight",
                "tenure_valid_records",
                "tenure_valid_weight",
                "owner_records",
                "owner_weight",
                "ownership_rate",
                "rooms_valid_records",
                "rooms_valid_weight",
                "rooms_weighted_sum",
                "mean_rooms_literal",
                "ownership_rate_index_2007",
                "mean_rooms_literal_index_2007",
            )
        }
        if min(
            values["base_head_records"],
            values["base_head_weight"],
            values["tenure_valid_records"],
            values["tenure_valid_weight"],
            values["rooms_valid_records"],
            values["rooms_valid_weight"],
        ) <= 0:
            raise ValueError(f"ACS has a nonpositive sample receipt in {year}")
        if not math.isclose(
            values["owner_weight"] / values["tenure_valid_weight"],
            values["ownership_rate"],
            rel_tol=0.0,
            abs_tol=1e-10,
        ):
            raise ValueError(f"ACS ownership receipt does not reproduce the rate in {year}")
        if not math.isclose(
            values["rooms_weighted_sum"] / values["rooms_valid_weight"],
            values["mean_rooms_literal"],
            rel_tol=0.0,
            abs_tol=1e-10,
        ):
            raise ValueError(f"ACS rooms receipt does not reproduce the mean in {year}")
        output[year] = values
    if len(sample_codes) != 5:
        raise ValueError("ACS sample_code must identify each annual ACS sample uniquely")
    for key in ("ownership_rate_index_2007", "mean_rooms_literal_index_2007"):
        if not math.isclose(output[2007][key], 1.0, rel_tol=0.0, abs_tol=1e-10):
            raise ValueError(f"ACS {key} is not indexed to 2007=1")
    return output, metadata, metadata_path, verification_path


def load_standardized_first_birth_timing(
    case_dir: Path,
    female_exposure_path: Path,
    fertility: Mapping[int, Mapping[str, float]],
    acs_metadata: Mapping[str, Any],
) -> tuple[dict[int, dict[str, float]], Path, dict[str, Any]]:
    """Standardize model age-specific first-birth rates to ACS female exposure.

    Raw model first-birth counts inherit the imposed distribution of household
    heads by age and are therefore not a valid timing comparison.  This routine
    divides the model's first-birth flow by its own age-cell mass, then applies
    national ACS female exposure shares over the matching four-year age bins.
    """
    model_path = (case_dir / "dated_period_fertility_by_age.csv").resolve()
    model_rows = read_csv(model_path)
    require_columns(
        model_rows,
        [
            "period",
            "calendar_year",
            "age_index",
            "age",
            "age_cell_start",
            "age_cell_midpoint",
            "age_mass",
            "birth_flow_first",
            "birth_flow_explicit",
            "birth_flow_topcode_adjusted",
            "age_specific_birth_rate_explicit",
            "age_specific_birth_rate_topcode_adjusted",
        ],
        "dated_period_fertility_by_age",
    )
    exposure_path = female_exposure_path.resolve()
    exposure_rows = read_csv(exposure_path)
    require_columns(
        exposure_rows,
        [
            "calendar_year",
            "sample_code",
            "age_lower",
            "age_upper",
            "female_records",
            "female_perwt_mass",
            "female_exposure_share_18_45",
        ],
        "ACS female exposure path",
    )
    output_files = acs_metadata.get("output_files")
    if not isinstance(output_files, dict):
        raise ValueError("ACS metadata lacks output_files")
    receipt = output_files.get(exposure_path.name)
    if not isinstance(receipt, dict):
        raise ValueError(f"ACS metadata has no receipt for {exposure_path.name}")
    recorded_hash = require_hash(receipt.get("sha256"), "ACS female exposure sha256")
    if recorded_hash != sha256(exposure_path):
        raise ValueError("ACS female exposure fingerprint does not match metadata")

    exposure_by_year: dict[int, dict[tuple[int, int], dict[str, float]]] = defaultdict(dict)
    for raw in exposure_rows:
        year = integer(raw["calendar_year"], "ACS exposure.calendar_year")
        if year not in CALENDAR_YEARS:
            raise ValueError(f"ACS exposure has unexpected calendar year {year}")
        lower = integer(raw["age_lower"], "ACS exposure.age_lower")
        upper = integer(raw["age_upper"], "ACS exposure.age_upper")
        if upper - lower != 3 or lower not in range(18, 43, 4):
            raise ValueError(f"ACS exposure has an unexpected age bin {lower}--{upper}")
        key = (lower, upper)
        if key in exposure_by_year[year]:
            raise ValueError(f"ACS exposure duplicates {year}, ages {lower}--{upper}")
        records = finite_float(raw["female_records"], "ACS exposure.female_records")
        mass = finite_float(raw["female_perwt_mass"], "ACS exposure.female_perwt_mass")
        share = finite_float(
            raw["female_exposure_share_18_45"], "ACS exposure.share"
        )
        if records <= 0 or mass <= 0 or not 0 < share < 1:
            raise ValueError(f"ACS exposure has an invalid receipt in {year}, {key}")
        exposure_by_year[year][key] = {
            "female_records": records,
            "female_perwt_mass": mass,
            "share": share,
        }
    expected_bins = {(lower, lower + 3) for lower in range(18, 43, 4)}
    if set(exposure_by_year) != set(CALENDAR_YEARS):
        raise ValueError("ACS female exposure does not cover all canonical dates")
    for year in CALENDAR_YEARS:
        if set(exposure_by_year[year]) != expected_bins:
            raise ValueError(f"ACS female exposure has incomplete age bins in {year}")
        share_sum = sum(cell["share"] for cell in exposure_by_year[year].values())
        if not math.isclose(share_sum, 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"ACS female exposure shares do not sum to one in {year}")

    model_by_year: dict[int, dict[int, dict[str, float]]] = defaultdict(dict)
    for raw in model_rows:
        year = integer(raw["calendar_year"], "dated fertility age.calendar_year")
        if year not in CALENDAR_YEARS:
            raise ValueError(
                f"dated_period_fertility_by_age has unexpected calendar year {year}"
            )
        period = integer(raw["period"], "dated fertility age.period")
        if period != CALENDAR_YEARS.index(year):
            raise ValueError(f"Dated fertility age period/year mismatch in {year}")
        age = integer(raw["age"], "dated fertility age.age")
        age_start = integer(
            raw["age_cell_start"], "dated fertility age.age_cell_start"
        )
        age_midpoint = integer(
            raw["age_cell_midpoint"], "dated fertility age.age_cell_midpoint"
        )
        if age != age_midpoint or age_midpoint != age_start + 2:
            raise ValueError(
                f"Dated fertility age labels are not explicit four-year midpoints in {year}"
            )
        age_index = integer(raw["age_index"], "dated fertility age.age_index")
        if age_start != 18 + 4 * age_index:
            raise ValueError(
                f"Dated fertility age index does not map to its four-year cell in {year}"
            )
        if age_start in model_by_year[year]:
            raise ValueError(
                f"Dated fertility age table duplicates {year}, cell start {age_start}"
            )
        age_mass = finite_float(raw["age_mass"], "dated fertility age.age_mass")
        if age_mass <= 0:
            raise ValueError(f"Nonpositive model age mass in {year}, age {age}")
        first_flow = finite_float(
            raw["birth_flow_first"], "dated fertility age.birth_flow_first"
        )
        explicit_flow = finite_float(
            raw["birth_flow_explicit"], "dated fertility age.birth_flow_explicit"
        )
        adjusted_flow = finite_float(
            raw["birth_flow_topcode_adjusted"],
            "dated fertility age.birth_flow_topcode_adjusted",
        )
        explicit_rate = finite_float(
            raw["age_specific_birth_rate_explicit"],
            "dated fertility age.explicit_rate",
        )
        adjusted_rate = finite_float(
            raw["age_specific_birth_rate_topcode_adjusted"],
            "dated fertility age.adjusted_rate",
        )
        if min(first_flow, explicit_flow, adjusted_flow) < -1e-14:
            raise ValueError(f"Negative birth flow in {year}, age {age}")
        if not math.isclose(explicit_flow / age_mass, explicit_rate, abs_tol=1e-12):
            raise ValueError(f"Explicit birth rate receipt fails in {year}, age {age}")
        if not math.isclose(adjusted_flow / age_mass, adjusted_rate, abs_tol=1e-12):
            raise ValueError(f"Adjusted birth rate receipt fails in {year}, age {age}")
        model_by_year[year][age_start] = {
            "age_index": float(age_index),
            "age_cell_midpoint": float(age_midpoint),
            "age_mass": age_mass,
            "first_birth_rate": first_flow / age_mass,
            "explicit_rate": explicit_rate,
            "adjusted_rate": adjusted_rate,
        }
    if set(model_by_year) != set(CALENDAR_YEARS):
        raise ValueError("Dated age-specific fertility does not cover all canonical dates")

    output: dict[int, dict[str, float]] = {}
    detailed_receipts: dict[str, Any] = {}
    for year in CALENDAR_YEARS:
        ages = sorted(model_by_year[year])
        if ages != list(range(18, 83, 4)):
            raise ValueError(f"Model fertility ages do not match the e5f grid in {year}")
        if any(
            abs(model_by_year[year][age]["first_birth_rate"]) > 1e-12
            for age in ages
            if age > 42
        ):
            raise ValueError(
                f"Model has first births outside the NCHS-comparable age cells in {year}"
            )
        explicit_sum = sum(cell["explicit_rate"] for cell in model_by_year[year].values())
        adjusted_sum = sum(cell["adjusted_rate"] for cell in model_by_year[year].values())
        if not math.isclose(
            explicit_sum,
            fertility[year]["period_tfr_explicit_states"],
            rel_tol=0.0,
            abs_tol=1e-10,
        ):
            raise ValueError(f"Age rates do not reconstruct explicit period TFR in {year}")
        if not math.isclose(
            adjusted_sum,
            fertility[year]["period_tfr_topcode_adjusted"],
            rel_tol=0.0,
            abs_tol=1e-10,
        ):
            raise ValueError(f"Age rates do not reconstruct adjusted period TFR in {year}")
        standardized_cells = []
        for age in range(18, 43, 4):
            exposure = exposure_by_year[year][(age, age + 3)]["share"]
            rate = model_by_year[year][age]["first_birth_rate"]
            midpoint = model_by_year[year][age]["age_cell_midpoint"]
            standardized_cells.append(
                {
                    "decision_age": age,
                    "cell_midpoint": midpoint,
                    "female_exposure_share": exposure,
                    "rate": rate,
                    "flow": exposure * rate,
                }
            )
        total = sum(cell["flow"] for cell in standardized_cells)
        if total <= 0:
            raise ValueError(f"Standardized model first-birth flow is nonpositive in {year}")
        mean_age = sum(
            cell["cell_midpoint"] * cell["flow"] for cell in standardized_cells
        ) / total
        share_30plus = sum(
            cell["flow"]
            for cell in standardized_cells
            if cell["decision_age"] >= 30
        ) / total
        output[year] = {
            "period_first_birth_mean_age": mean_age,
            "period_first_birth_share_age30plus": share_30plus,
        }
        detailed_receipts[str(year)] = {
            "standardized_first_birth_intensity": total,
            "age_cells": standardized_cells,
            "explicit_period_tfr_reconstructed": explicit_sum,
            "topcode_adjusted_period_tfr_reconstructed": adjusted_sum,
        }
    provenance = {
        "method": (
            "model first-birth flow divided by model age-cell mass, multiplied by "
            "ACS national female exposure share over matching four-year ages 18--45"
        ),
        "timing_operator": "boundary_collapsed_four_year_cells_midpoint_labels",
        "age_label": (
            "dated_period_fertility_by_age.csv must contain explicit age_cell_start "
            "and age_cell_midpoint columns, with age equal to the midpoint; the "
            "reported midpoint is used directly with no post-processing age shift"
        ),
        "boundary_treatment": (
            "the model has seven cells starting 18,22,...,42; NCHS observations "
            "below age 22 enter the first cell and observations age 42 or older "
            "enter the terminal cell"
        ),
        "female_exposure_csv_sha256": recorded_hash,
        "dated_period_fertility_by_age_sha256": sha256(model_path),
        "receipts": detailed_receipts,
    }
    return output, model_path, provenance


def load_historical_observations(
    comparison_path: Path, metadata_path: Path
) -> tuple[dict[int, dict[str, float]], dict[str, Any], dict[str, Any]]:
    comparison_path = comparison_path.resolve()
    metadata_path = metadata_path.resolve()
    rows = read_csv(comparison_path)
    require_columns(
        rows,
        [
            "calendar_year",
            "observed_population_index",
            "observed_nominal_house_price_index",
            "observed_nominal_rent_index",
            "observed_real_house_price_index",
            "observed_real_rent_index",
            "observed_price_to_rent_index",
            "observed_tfr",
        ],
        "historical comparison",
    )
    by_year = canonical_rows_by_year(rows, label="historical comparison")
    metadata = read_json(metadata_path)
    if metadata.get("aggregation") != (
        "arithmetic mean of nonmissing observations within calendar year"
    ):
        raise ValueError("Historical aggregation contract is missing or changed")
    comparison_summary_path = comparison_path.parent / "historical_comparison_summary.json"
    comparison_summary = read_json(comparison_summary_path)
    if comparison_summary.get("normalization") != "annual averages indexed to 2007=1":
        raise ValueError("Historical comparison normalization contract changed")
    if integer(comparison_summary.get("number_of_model_dates"), "historical date count") != 5:
        raise ValueError("Historical comparison summary does not cover five dates")
    if integer(comparison_summary.get("start_year"), "historical start year") != 2007 or integer(
        comparison_summary.get("last_overlapping_year"), "historical end year"
    ) != 2023:
        raise ValueError("Historical comparison summary has a different date window")
    if comparison_summary.get("population_units") != (
        "Census HH-3 households with heads ages 18--85 versus model adult-household mass"
    ):
        raise ValueError("Historical household-population measurement contract changed")
    series = metadata.get("series")
    if not isinstance(series, dict):
        raise ValueError("Historical metadata lacks the FRED/World Bank series ledger")
    expected_ids = {
        "consumer_price": "CPIAUCSL",
        "fertility": "SPDYNTFRTINUSA",
        "house_price": "CSUSHPINSA",
        "rent": "CUUR0000SEHA",
    }
    annual_values: dict[str, dict[int, float]] = {}
    raw_receipts: dict[str, Any] = {}
    summary_series = comparison_summary.get("sources")
    if not isinstance(summary_series, dict):
        raise ValueError("Historical comparison summary lacks its source ledger")
    for name, expected_id in expected_ids.items():
        receipt = series.get(name)
        if not isinstance(receipt, dict):
            raise ValueError(f"Historical metadata lacks series {name}")
        if receipt.get("id") != expected_id:
            raise ValueError(f"Historical series {name} is not {expected_id}")
        expected_hash = require_hash(
            receipt.get("raw_sha256"), f"historical {name} raw_sha256"
        )
        if not receipt.get("download_url") or not receipt.get("fred_page"):
            raise ValueError(f"Historical series {name} lacks source URLs")
        summary_receipt = summary_series.get(name)
        if not isinstance(summary_receipt, dict):
            raise ValueError(f"Historical comparison summary lacks {name}")
        for field in ("id", "raw_sha256", "download_url", "fred_page"):
            if summary_receipt.get(field) != receipt.get(field):
                raise ValueError(f"Historical {name} metadata and summary disagree on {field}")
        recorded_raw = Path(str(receipt.get("raw_file", "")))
        sibling_raw = metadata_path.parent / "observed_fred_raw" / f"{expected_id}.csv"
        raw_path = recorded_raw if recorded_raw.is_file() else sibling_raw
        if not raw_path.is_file():
            raise FileNotFoundError(f"Missing raw historical series: {raw_path}")
        if sha256(raw_path) != expected_hash:
            raise ValueError(f"Raw historical series fingerprint changed for {name}")
        raw_rows = read_csv(raw_path)
        require_columns(raw_rows, ["observation_date", expected_id], f"raw {expected_id}")
        values_by_year: dict[int, list[float]] = defaultdict(list)
        for raw in raw_rows:
            date = raw["observation_date"]
            if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", date):
                raise ValueError(f"Malformed observation date in {expected_id}: {date!r}")
            value = optional_float(raw[expected_id], f"raw {expected_id}")
            if value is not None:
                values_by_year[int(date[:4])].append(value)
        missing_raw_years = [year for year in CALENDAR_YEARS if not values_by_year[year]]
        if missing_raw_years:
            raise ValueError(f"Raw {expected_id} lacks years {missing_raw_years}")
        annual_values[name] = {
            year: sum(values_by_year[year]) / len(values_by_year[year])
            for year in CALENDAR_YEARS
        }
        raw_receipts[name] = {
            "path": str(raw_path.resolve()),
            "sha256": expected_hash,
            "series_id": expected_id,
            "annual_nonmissing_counts": {
                str(year): len(values_by_year[year]) for year in CALENDAR_YEARS
            },
        }

    households = summary_series.get("households")
    if not isinstance(households, dict) or not isinstance(
        households.get("model_age_totals_thousands"), dict
    ):
        raise ValueError("Historical comparison lacks household population receipts")
    if not isinstance(households.get("published_hh3_totals_thousands"), dict):
        raise ValueError("Historical comparison lacks published HH-3 receipts")
    if not str(households.get("url", "")).startswith("https://www2.census.gov/"):
        raise ValueError("Historical household receipt lacks its Census source URL")
    if households.get("revision_note") != "The revised 2011 row is used.":
        raise ValueError("Historical household receipt does not declare the revised 2011 row")
    household_totals = {
        year: finite_float(
            households["model_age_totals_thousands"].get(str(year)),
            f"historical household total {year}",
        )
        for year in CALENDAR_YEARS
    }
    published_households = {
        year: finite_float(
            households["published_hh3_totals_thousands"].get(str(year)),
            f"published HH-3 total {year}",
        )
        for year in CALENDAR_YEARS
    }
    if any(
        household_totals[year] <= 0
        or published_households[year] <= 0
        or household_totals[year] > published_households[year]
        for year in CALENDAR_YEARS
    ):
        raise ValueError("Historical household receipts contain invalid age-adjusted totals")
    base_price = annual_values["house_price"][2007]
    base_rent = annual_values["rent"][2007]
    base_cpi = annual_values["consumer_price"][2007]
    base_households = household_totals[2007]
    output: dict[int, dict[str, float]] = {}
    for year in CALENDAR_YEARS:
        row = by_year[year]
        output[year] = {
            key: finite_float(row[key], f"historical.{key}")
            for key in (
                "observed_population_index",
                "observed_nominal_house_price_index",
                "observed_nominal_rent_index",
                "observed_real_house_price_index",
                "observed_real_rent_index",
                "observed_price_to_rent_index",
                "observed_tfr",
            )
        }
        price_index = annual_values["house_price"][year] / base_price
        rent_index = annual_values["rent"][year] / base_rent
        cpi_index = annual_values["consumer_price"][year] / base_cpi
        reconstructed = {
            "observed_population_index": household_totals[year] / base_households,
            "observed_nominal_house_price_index": price_index,
            "observed_nominal_rent_index": rent_index,
            "observed_real_house_price_index": price_index / cpi_index,
            "observed_real_rent_index": rent_index / cpi_index,
            "observed_price_to_rent_index": price_index / rent_index,
            "observed_tfr": annual_values["fertility"][year],
        }
        for key, expected in reconstructed.items():
            if not math.isclose(
                output[year][key], expected, rel_tol=0.0, abs_tol=2e-12
            ):
                raise ValueError(
                    f"Historical comparison {key} in {year} does not reproduce raw receipts"
                )
    for key in (
        "observed_population_index",
        "observed_real_house_price_index",
        "observed_real_rent_index",
        "observed_price_to_rent_index",
    ):
        if not math.isclose(output[2007][key], 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(f"Historical {key} is not declared on a 2007=1 scale")
    for endpoint_name, year in (("first_row", 2007), ("last_row", 2023)):
        endpoint = comparison_summary.get(endpoint_name)
        if not isinstance(endpoint, dict):
            raise ValueError(f"Historical comparison summary lacks {endpoint_name}")
        for key, value in output[year].items():
            if not math.isclose(
                finite_float(endpoint.get(key), f"{endpoint_name}.{key}"),
                value,
                rel_tol=0.0,
                abs_tol=2e-12,
            ):
                raise ValueError(f"Historical summary endpoint disagrees on {key}")
    contract = {
        "comparison_summary_path": str(comparison_summary_path.resolve()),
        "comparison_summary_sha256": sha256(comparison_summary_path),
        "raw_series": raw_receipts,
        "household_receipt": {
            "source": households.get("title"),
            "url": households.get("url"),
            "model_age_totals_thousands": {
                str(year): household_totals[year] for year in CALENDAR_YEARS
            },
            "published_hh3_totals_thousands": {
                str(year): published_households[year] for year in CALENDAR_YEARS
            },
        },
        "reconstruction_status": "all five dates reproduced from raw annual means",
    }
    return output, metadata, contract


def load_nchs_first_birth_timing(
    counts_path: Path, manifest_path: Path, metadata_path: Path
) -> tuple[dict[int, dict[str, float]], dict[str, Any]]:
    counts_path = counts_path.resolve()
    manifest_path = manifest_path.resolve()
    metadata_path = metadata_path.resolve()
    counts = read_csv(counts_path)
    manifest = read_csv(manifest_path)
    metadata = read_json(metadata_path)
    if sha256(metadata_path) != EXPECTED_NCHS_METADATA_SHA256:
        raise ValueError("Supplied NCHS metadata fingerprint is not the pinned contract")
    require_columns(counts, ["year", "age", "n_first_births"], "NCHS counts")
    require_columns(
        manifest,
        ["year", "file", "size_bytes", "sha256", "cache_contract"],
        "NCHS manifest",
    )

    expected_contract_id = EXPECTED_NCHS_CONTRACT_ID
    expected_cache_contract = EXPECTED_NCHS_CACHE_CONTRACT
    if integer(metadata.get("schema_version"), "NCHS metadata schema_version") != 1:
        raise ValueError("Unexpected NCHS metadata schema")
    if metadata.get("contract_id") != expected_contract_id:
        raise ValueError("Unexpected NCHS timing contract identifier")
    if metadata.get("cache_contract") != expected_cache_contract:
        raise ValueError("Unexpected NCHS cache contract")
    vintage = metadata.get("data_vintage")
    if not isinstance(vintage, dict) or integer(
        vintage.get("year_min"), "NCHS year_min"
    ) != 1987 or integer(vintage.get("year_max"), "NCHS year_max") != 2023:
        raise ValueError("NCHS cache does not use the frozen 1987--2023 vintage")
    if metadata.get("contains_five_dated_period_targets") is not False:
        raise ValueError("NCHS metadata must distinguish pooled targets from dated periods")
    if metadata.get("target_interpretation") != (
        "pooled cohort timing target; annual period comparisons must apply the "
        "declared operator directly to the year-age cache"
    ):
        raise ValueError("NCHS annual-period interpretation is missing or changed")
    raw_receipt = metadata.get("raw_microdata_receipt")
    if not isinstance(raw_receipt, dict) or raw_receipt != {
        "status": "every 1987--2023 national raw file is content-hashed before cache reuse",
        "manifest": "first_birth_counts_manifest.csv",
        "manifest_columns": [
            "year",
            "file",
            "size_bytes",
            "sha256",
            "cache_contract",
        ],
    }:
        raise ValueError("NCHS raw-microdata content-hash receipt changed")
    live_birth_order = metadata.get("live_birth_order_codes")
    if not isinstance(live_birth_order, dict) or live_birth_order != {
        "first_birth": 1,
        "legacy_dlivord_unknown": 99,
        "lbo_rec_unknown": 9,
        "missing_values_also_excluded": True,
    }:
        raise ValueError("NCHS live-birth-order code contract changed")

    file_hashes = metadata.get("file_sha256")
    if not isinstance(file_hashes, dict):
        raise ValueError("NCHS metadata lacks its output hash ledger")
    expected_files = {
        "first_birth_counts_year_age.csv",
        "first_birth_counts_manifest.csv",
        "order_exclusion_shares_by_year.csv",
        "cohort_truncation.csv",
        "timing_targets.csv",
        "timing_targets_exact_age.csv",
        "timing_targets_model_comparable.csv",
        "timing_target_contract.csv",
        "timing_age_bin_counts.csv",
        "build_first_birth_timing.R",
        "test_first_birth_timing_targets.R",
    }
    if set(file_hashes) != expected_files:
        raise ValueError("NCHS metadata output hash ledger is incomplete or expanded")
    metadata_dir = metadata_path.parent
    for name in sorted(expected_files):
        recorded = require_hash(file_hashes[name], f"NCHS file_sha256.{name}")
        path = metadata_dir / name
        if not path.is_file() or sha256(path) != recorded:
            raise ValueError(f"NCHS output fingerprint changed for {name}")
    if counts_path.name != "first_birth_counts_year_age.csv" or sha256(
        counts_path
    ) != file_hashes["first_birth_counts_year_age.csv"]:
        raise ValueError("Supplied NCHS counts do not match the metadata contract")
    if manifest_path.name != "first_birth_counts_manifest.csv" or sha256(
        manifest_path
    ) != file_hashes["first_birth_counts_manifest.csv"]:
        raise ValueError("Supplied NCHS manifest does not match the metadata contract")

    source_members = metadata.get("source_bundle_members")
    expected_source_keys = (
        "first_birth_counts_year_age.csv",
        "first_birth_counts_manifest.csv",
        "official_1987_nchs_documentation_pdf",
        "official_1989_nchs_documentation_pdf",
        "official_2003_nchs_documentation_pdf",
        "cache_contract",
    )
    if not isinstance(source_members, dict) or set(source_members) != set(
        expected_source_keys
    ):
        raise ValueError("NCHS source bundle member ledger changed")
    if (
        source_members["first_birth_counts_year_age.csv"]
        != file_hashes["first_birth_counts_year_age.csv"]
        or source_members["first_birth_counts_manifest.csv"]
        != file_hashes["first_birth_counts_manifest.csv"]
        or source_members["cache_contract"] != expected_cache_contract
        or any(
            source_members[name] != expected_hash
            for name, expected_hash in EXPECTED_NCHS_DOCUMENT_SHA256.items()
        )
    ):
        raise ValueError("NCHS source bundle members do not reconcile")
    pre_2003 = metadata.get("pre_2003_measurement")
    if not isinstance(pre_2003, dict) or any(
        pre_2003.get(key) != value
        for key, value in {
            "maternal_age_field": "dmage",
            "maternal_age_field_label": "Age of Mother",
            "live_birth_order_field": "dlivord",
            "live_birth_order_field_label": "Detail Live Birth Order",
            "years": "1987--2002",
            "official_1987_document_url": (
                "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/"
                "Dataset_Documentation/DVS/natality/Nat1987doc.pdf"
            ),
            "official_1987_document_sha256": EXPECTED_NCHS_DOCUMENT_SHA256[
                "official_1987_nchs_documentation_pdf"
            ],
            "official_1989_document_url": (
                "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/"
                "Dataset_Documentation/DVS/natality/Nat1989doc.pdf"
            ),
            "official_1989_document_sha256": EXPECTED_NCHS_DOCUMENT_SHA256[
                "official_1989_nchs_documentation_pdf"
            ],
        }.items()
    ):
        raise ValueError("NCHS 1987--2002 field/documentation contract changed")
    source_2003 = metadata.get("maternal_age_source_2003")
    if not isinstance(source_2003, dict) or any(
        source_2003.get(key) != value
        for key, value in {
            "field": "MAGER41",
            "transformation": (
                "code 01 assigned 14 only to preserve its under-15 boundary status; "
                "codes 02:41 decoded as code+13"
            ),
            "official_document_url": (
                "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/"
                "Dataset_Documentation/DVS/natality/Nat2003doc.pdf"
            ),
            "official_document_sha256": EXPECTED_NCHS_DOCUMENT_SHA256[
                "official_2003_nchs_documentation_pdf"
            ],
            "official_document_size_bytes": 16_116_001,
        }.items()
    ):
        raise ValueError("NCHS 2003 MAGER41 documentation contract changed")
    source_payload = "\n".join(
        f"{name}={source_members[name]}" for name in expected_source_keys
    )
    source_bundle_sha256 = require_hash(
        metadata.get("source_bundle_sha256"), "NCHS source_bundle_sha256"
    )
    if source_bundle_sha256 != EXPECTED_NCHS_SOURCE_BUNDLE_SHA256 or (
        hashlib.sha256(source_payload.encode()).hexdigest() != source_bundle_sha256
    ):
        raise ValueError("NCHS source bundle fingerprint does not reproduce")

    expected_contract_files = (
        "timing_target_contract.csv",
        "timing_targets_model_comparable.csv",
        "timing_age_bin_counts.csv",
        "timing_targets_exact_age.csv",
        "build_first_birth_timing.R",
        "test_first_birth_timing_targets.R",
    )
    contract_members = metadata.get("contract_bundle_members")
    if not isinstance(contract_members, dict) or set(contract_members) != set(
        expected_contract_files
    ):
        raise ValueError("NCHS contract bundle member ledger changed")
    if any(contract_members[name] != file_hashes[name] for name in expected_contract_files):
        raise ValueError("NCHS contract bundle members do not match output hashes")
    contract_payload = "\n".join(
        f"{name}={contract_members[name]}" for name in expected_contract_files
    )
    contract_bundle_sha256 = require_hash(
        metadata.get("contract_bundle_sha256"), "NCHS contract_bundle_sha256"
    )
    if contract_bundle_sha256 != EXPECTED_NCHS_CONTRACT_BUNDLE_SHA256 or (
        hashlib.sha256(contract_payload.encode()).hexdigest() != contract_bundle_sha256
    ):
        raise ValueError("NCHS contract bundle fingerprint does not reproduce")

    operator = metadata.get("recommended_operator")
    if not isinstance(operator, dict) or operator.get("id") != (
        "boundary_collapsed_midpoint_v1"
    ):
        raise ValueError("NCHS recommended timing operator changed")
    expected_boundary = (
        "exact/decoded ages <=21 map to cell start 18; 22-25 to 22; 26-29 to 26; "
        "30-33 to 30; 34-37 to 34; 38-41 to 38; ages >=42 map to cell start 42"
    )
    if operator.get("boundary_treatment") != expected_boundary:
        raise ValueError("NCHS boundary treatment changed")
    if [integer(value, "NCHS bin start") for value in operator.get("bin_starts", [])] != [
        18,
        22,
        26,
        30,
        34,
        38,
        42,
    ] or [
        integer(value, "NCHS midpoint")
        for value in operator.get("reported_midpoint_labels", [])
    ] != [20, 24, 28, 32, 36, 40, 44]:
        raise ValueError("NCHS four-year cells or midpoint labels changed")
    if operator.get("mean_formula") != (
        "first-birth-count-weighted mean of midpoint labels 20,24,...,44"
    ) or operator.get("share_30plus_formula") != (
        "sum first-birth counts with bin_start >= 30 divided by all used first-birth counts"
    ):
        raise ValueError("NCHS timing formulas changed")
    if operator.get("boundary_exclusions") != "none":
        raise ValueError("NCHS operator silently excludes boundary ages")

    target_contract = read_csv(metadata_dir / "timing_target_contract.csv")
    comparable_targets = read_csv(metadata_dir / "timing_targets_model_comparable.csv")
    require_columns(
        target_contract,
        [
            "target_status",
            "timing_operator",
            "window",
            "moment",
            "target_value",
            "declared_standard_error",
            "inverse_variance_weight",
            "support_definition",
            "age_label_convention",
            "age_30plus_rule",
        ],
        "NCHS timing target contract",
    )
    require_columns(
        comparable_targets,
        [
            "window",
            "timing_operator",
            "bin_label_convention",
            "support_definition",
            "mean_age_first_birth",
            "share_30plus",
            "n_first_births_used",
        ],
        "NCHS comparable targets",
    )
    recommended_rows = [
        row
        for row in target_contract
        if row["target_status"] == "recommended_model_comparable"
    ]
    if len(recommended_rows) != 4 or {
        (row["window"], row["moment"]) for row in recommended_rows
    } != {
        (window, moment)
        for window in ("primary_1979_1984", "sensitivity_1975_1980")
        for moment in ("mean_age_first_birth", "share_30plus")
    } or any(row["timing_operator"] != "boundary_collapsed" for row in recommended_rows):
        raise ValueError("NCHS recommended target rows do not implement the operator")
    comparable_boundary = [
        row for row in comparable_targets if row["timing_operator"] == "boundary_collapsed"
    ]
    if len(comparable_boundary) != 2 or {
        row["window"] for row in comparable_boundary
    } != {"primary_1979_1984", "sensitivity_1975_1980"} or any(
        row["bin_label_convention"] != "midpoint" for row in comparable_boundary
    ):
        raise ValueError("NCHS comparable target receipt changed")
    primary_rows = [
        row for row in comparable_boundary if row["window"] == "primary_1979_1984"
    ]
    primary = primary_rows[0]
    primary_metadata = operator.get("primary_target")
    if not isinstance(primary_metadata, dict):
        raise ValueError("NCHS metadata lacks its primary target receipt")
    if set(primary_metadata) != set(EXPECTED_NCHS_PRIMARY_TARGET):
        raise ValueError("NCHS primary target receipt schema changed")
    for key, expected in EXPECTED_NCHS_PRIMARY_TARGET.items():
        if not math.isclose(
            finite_float(primary_metadata.get(key), f"NCHS primary {key}"),
            float(expected),
            rel_tol=0.0,
            abs_tol=2e-12,
        ):
            raise ValueError(f"NCHS primary target changed for {key}")
    for metadata_key, csv_key in (
        ("mean_age_first_birth", "mean_age_first_birth"),
        ("share_30plus", "share_30plus"),
        ("n_first_births", "n_first_births_used"),
    ):
        if not math.isclose(
            finite_float(primary_metadata.get(metadata_key), f"NCHS primary {metadata_key}"),
            finite_float(primary[csv_key], f"NCHS comparable {csv_key}"),
            rel_tol=0.0,
            abs_tol=2e-12,
        ):
            raise ValueError("NCHS primary target metadata and CSV disagree")
    recommended_by_moment = {
        row["moment"]: row
        for row in recommended_rows
        if row["window"] == "primary_1979_1984"
    }
    primary_row_contract = {
        "mean_age_first_birth": (
            EXPECTED_NCHS_PRIMARY_TARGET["mean_age_first_birth"],
            EXPECTED_NCHS_PRIMARY_TARGET["mean_standard_error"],
            EXPECTED_NCHS_PRIMARY_TARGET["mean_inverse_variance_weight"],
        ),
        "share_30plus": (
            EXPECTED_NCHS_PRIMARY_TARGET["share_30plus"],
            EXPECTED_NCHS_PRIMARY_TARGET["share_30plus_standard_error"],
            EXPECTED_NCHS_PRIMARY_TARGET["share_30plus_inverse_variance_weight"],
        ),
    }
    if set(recommended_by_moment) != set(primary_row_contract):
        raise ValueError("NCHS primary target-contract rows changed")
    for moment, expected_values in primary_row_contract.items():
        receipt = recommended_by_moment[moment]
        actual_values = (
            finite_float(receipt["target_value"], f"NCHS {moment} target"),
            finite_float(
                receipt["declared_standard_error"], f"NCHS {moment} standard error"
            ),
            finite_float(receipt["inverse_variance_weight"], f"NCHS {moment} weight"),
        )
        if any(
            not math.isclose(actual, float(expected), rel_tol=0.0, abs_tol=2e-12)
            for actual, expected in zip(actual_values, expected_values)
        ):
            raise ValueError(f"NCHS primary target contract changed for {moment}")

    active_targets = ACTIVE_TARGET_SYSTEM.targets_dict()
    active_weights = ACTIVE_TARGET_SYSTEM.weights_dict()
    for profile_name, target_key, weight_key in (
        (
            "mean_age_first_birth",
            "mean_age_first_birth",
            "mean_inverse_variance_weight",
        ),
        (
            "share_first_births_age30plus",
            "share_30plus",
            "share_30plus_inverse_variance_weight",
        ),
    ):
        if not math.isclose(
            active_targets[profile_name],
            float(EXPECTED_NCHS_PRIMARY_TARGET[target_key]),
            rel_tol=0.0,
            abs_tol=2e-12,
        ) or not math.isclose(
            active_weights[profile_name],
            float(EXPECTED_NCHS_PRIMARY_TARGET[weight_key]),
            rel_tol=0.0,
            abs_tol=2e-12,
        ):
            raise ValueError(f"Active E5 target profile disagrees with NCHS for {profile_name}")

    manifest_by_year: dict[int, Mapping[str, str]] = {}
    for row in manifest:
        year = integer(row["year"], "NCHS manifest.year")
        if year in manifest_by_year:
            raise ValueError(f"NCHS manifest has duplicate year {year}")
        if not row["file"]:
            raise ValueError(f"NCHS manifest lacks the source file for {year}")
        if row["cache_contract"] != expected_cache_contract:
            raise ValueError(f"NCHS manifest has the wrong cache contract in {year}")
        raw_sha256 = require_hash(row["sha256"], f"NCHS manifest.sha256[{year}]")
        size = integer(row["size_bytes"], "NCHS manifest.size_bytes")
        if size <= 0:
            raise ValueError(f"NCHS manifest has nonpositive source size in {year}")
        raw_path = Path(row["file"])
        if raw_path.exists() and raw_path.stat().st_size != size:
            raise ValueError(f"NCHS source size changed for {year}: {raw_path}")
        manifest_by_year[year] = {**row, "sha256": raw_sha256}
    grouped: dict[int, list[tuple[int, float]]] = defaultdict(list)
    seen: set[tuple[int, int]] = set()
    for row in counts:
        year = integer(row["year"], "NCHS counts.year")
        age = integer(row["age"], "NCHS counts.age")
        count = finite_float(row["n_first_births"], "NCHS counts.n_first_births")
        if not 12 <= age <= 49 or count < 0:
            raise ValueError(f"Invalid NCHS age/count in {year}: age={age}, count={count}")
        key = (year, age)
        if key in seen:
            raise ValueError(f"NCHS cache has duplicate year-age cell {key}")
        seen.add(key)
        grouped[year].append((age, count))
    expected_vintage = set(range(1987, 2024))
    if set(grouped) != expected_vintage or set(manifest_by_year) != expected_vintage:
        raise ValueError("NCHS cache/manifest does not contain exactly 1987--2023")
    expected_ages = set(range(12, 50))
    for year in sorted(expected_vintage):
        if {age for age, _ in grouped[year]} != expected_ages:
            raise ValueError(f"NCHS cache has incomplete exact/decoded ages in {year}")
    output: dict[int, dict[str, float]] = {}
    dated_receipts: dict[str, Any] = {}
    for year in CALENDAR_YEARS:
        cells = grouped[year]
        total = sum(count for _, count in cells)
        if total <= 0:
            raise ValueError(f"NCHS first-birth count is nonpositive in {year}")
        binned = {start: 0.0 for start, _ in NCHS_AGE_BINS}
        for age, count in cells:
            if age <= 21:
                start = 18
            elif age >= 42:
                start = 42
            else:
                start = 22 + 4 * ((age - 22) // 4)
            binned[start] += count
        if not math.isclose(math.fsum(binned.values()), total, abs_tol=1e-9):
            raise ValueError(f"NCHS boundary-collapse accounting fails in {year}")
        output[year] = {
            "period_first_birth_mean_age": math.fsum(
                (start + 2) * count for start, count in binned.items()
            )
            / total,
            "period_first_birth_share_age30plus": math.fsum(
                count for start, count in binned.items() if start >= 30
            )
            / total,
            "n_first_births": total,
        }
        dated_receipts[str(year)] = {
            "total_first_births": total,
            "boundary_exclusions": 0,
            "bin_counts": {str(start): count for start, count in binned.items()},
        }
    provenance = {
        "counts_sha256": sha256(counts_path),
        "manifest_sha256": sha256(manifest_path),
        "metadata_sha256": sha256(metadata_path),
        "source_bundle_sha256": metadata["source_bundle_sha256"],
        "contract_bundle_sha256": metadata["contract_bundle_sha256"],
        "contract_id": expected_contract_id,
        "manifest_years": [1987, 2023],
        "annual_period_operator": {
            "id": operator["id"],
            "boundary_treatment": expected_boundary,
            "bin_starts": [18, 22, 26, 30, 34, 38, 42],
            "reported_midpoint_labels": [20, 24, 28, 32, 36, 40, 44],
            "boundary_exclusions": "none",
            "share_30plus_rule": "bin_start >= 30",
        },
        "dated_receipts": dated_receipts,
        "calendar_source_files": {
            str(year): {
                "path": manifest_by_year[year]["file"],
                "size_bytes": integer(
                    manifest_by_year[year]["size_bytes"], "NCHS source size"
                ),
                "sha256": manifest_by_year[year]["sha256"],
            }
            for year in CALENDAR_YEARS
        },
    }
    return output, provenance


def row(
    *,
    section_order: int,
    series_order: int,
    year: int,
    series_id: str,
    series_label: str,
    classification: str,
    data_value: float | None,
    model_value: float | None,
    model_sensitivity_value: float | None = None,
    units: str,
    normalization: str,
    data_source: str,
    model_source: str,
    comparability: str,
    notes: str,
) -> dict[str, Any]:
    return {
        "section_order": section_order,
        "series_order": series_order,
        "calendar_year": year,
        "series_id": series_id,
        "series_label": series_label,
        "classification": classification,
        "data_value": "" if data_value is None else data_value,
        "model_value": "" if model_value is None else model_value,
        "model_sensitivity_value": (
            "" if model_sensitivity_value is None else model_sensitivity_value
        ),
        "units": units,
        "normalization": normalization,
        "data_source": data_source,
        "model_source": model_source,
        "comparability": comparability,
        "notes": notes,
    }


def build_comparison_rows(
    transition: Mapping[int, Mapping[str, float]],
    moments: Mapping[int, Mapping[str, Mapping[str, Any]]],
    fertility: Mapping[int, Mapping[str, float]],
    standardized_timing: Mapping[int, Mapping[str, float]],
    acs: Mapping[int, Mapping[str, float]],
    historical: Mapping[int, Mapping[str, float]],
    nchs: Mapping[int, Mapping[str, float]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for year in CALENDAR_YEARS:
        rows.append(
            row(
                section_order=1,
                series_order=1,
                year=year,
                series_id="household_population_index",
                series_label="Adult-household population",
                classification="imposed_bridge_matched_by_construction_not_fit",
                data_value=historical[year]["observed_population_index"],
                model_value=transition[year]["population_index"],
                units="index",
                normalization="source-declared 2007=1",
                data_source="Census HH-3 adjusted to heads ages 18--85",
                model_source="transition_path.population_index",
                comparability="exact bridge object",
                notes="Imposed household-formation/migration path; never count as model fit.",
            )
        )
        rows.append(
            row(
                section_order=1,
                series_order=2,
                year=year,
                series_id="householder_age_bridge_max_gap",
                series_label="Householder-age bridge maximum gap",
                classification="imposed_bridge_numerical_audit_not_fit",
                data_value=0.0,
                model_value=transition[year]["census_age_bridge_maximum_target_gap"],
                units="adult-household mass",
                normalization="none",
                data_source="bridge identity benchmark",
                model_source="transition_path.census_age_bridge_maximum_target_gap",
                comparability="numerical residual",
                notes="Zero is the imposed age-distribution accounting identity.",
            )
        )

    for year in CALENDAR_YEARS:
        rows.append(
            row(
                section_order=2,
                series_order=1,
                year=year,
                series_id="period_tfr_explicit_states",
                series_label="Period fertility rate",
                classification="descriptive_noncomparable_holdout",
                data_value=historical[year]["observed_tfr"],
                model_value=fertility[year]["period_tfr_explicit_states"],
                model_sensitivity_value=fertility[year]["period_tfr_topcode_adjusted"],
                units="births per exposure",
                normalization="levels; no indexing",
                data_source="World Bank period TFR (FRED SPDYNTFRTINUSA)",
                model_source="dated_period_fertility.csv",
                comparability="descriptive only: denominator mismatch",
                notes=(
                    "Main model line uses explicit parity states; sensitivity assigns the "
                    "observed mean to the 3+ bin. Model exposure is adult-household mass, "
                    "not literal female exposure."
                ),
            )
        )
        rows.append(
            row(
                section_order=2,
                series_order=2,
                year=year,
                series_id="period_first_birth_mean_age",
                series_label="Mean age among period first births",
                classification="untargeted_holdout",
                data_value=nchs[year]["period_first_birth_mean_age"],
                model_value=standardized_timing[year]["period_first_birth_mean_age"],
                units="years",
                normalization="levels; no indexing",
                data_source="NCHS natality first-birth records",
                model_source=(
                    "dated_period_fertility_by_age.csv standardized with ACS female exposure"
                ),
                comparability="period first-birth timing after female-exposure standardization",
                notes=(
                    "Both lines use seven four-year cells labeled by midpoints "
                    "20,24,...,44. NCHS ages below 22 collapse into the first cell "
                    "and ages 42+ into the terminal cell."
                ),
            )
        )
        rows.append(
            row(
                section_order=2,
                series_order=3,
                year=year,
                series_id="period_first_birth_share_age30plus",
                series_label="Period first births at age 30+",
                classification="untargeted_holdout",
                data_value=nchs[year]["period_first_birth_share_age30plus"],
                model_value=standardized_timing[year][
                    "period_first_birth_share_age30plus"
                ],
                units="share",
                normalization="levels; no indexing",
                data_source="NCHS natality first-birth records",
                model_source=(
                    "dated_period_fertility_by_age.csv standardized with ACS female exposure"
                ),
                comparability="period first-birth timing after female-exposure standardization",
                notes=(
                    "Both lines use the boundary-collapsed four-year operator; age "
                    "30+ means cell start at least 30. This is not the cohort target."
                ),
            )
        )
        rows.append(
            row(
                section_order=2,
                series_order=4,
                year=year,
                series_id="aggregate_ownership_rate_18_85",
                series_label="Homeownership, household heads ages 18--85",
                classification="untargeted_holdout",
                data_value=acs[year]["ownership_rate"],
                model_value=transition[year]["owner_rate"],
                units="share",
                normalization="levels; no indexing",
                data_source="national ACS household heads, valid tenure",
                model_source="transition_path.owner_rate",
                comparability="aggregate household-head ownership",
                notes="Distinct from the maintained prime-age ownership target.",
            )
        )
        rows.append(
            row(
                section_order=2,
                series_order=5,
                year=year,
                series_id="mean_rooms_literal_18_85",
                series_label="Mean occupied rooms, household heads ages 18--85",
                classification=(
                    "untargeted_holdout_cross_section_only_coding_break"
                    if year == 2007
                    else "untargeted_holdout_longitudinal_2011_2023"
                ),
                data_value=acs[year]["mean_rooms_literal"],
                model_value=moments[year][ROOM_MOMENT]["model"],
                units="literal ACS rooms",
                normalization="levels; no indexing",
                data_source="national ACS household heads with ROOMS>0",
                model_source=f"dated_model_moments.{ROOM_MOMENT}",
                comparability=(
                    "cross-section only; excluded from longitudinal statistics"
                    if year == 2007
                    else "same age range and literal room unit; 2011--2023 window"
                ),
                notes=(
                    "This independently built national ACS path is not the maintained "
                    "target series. The 2007 code 9 means 9+, creating a coding break "
                    "relative to later samples."
                ),
            )
        )

    first_cost = transition[2007]["housing_user_cost"]
    if first_cost <= 0:
        raise ValueError("The model's initial implicit rental cost is nonpositive")
    implicit_cost_index = {
        year: transition[year]["housing_user_cost"] / first_cost
        for year in CALENDAR_YEARS
    }
    model_price_to_implicit_cost = {
        year: transition[year]["asset_price_index"] / implicit_cost_index[year]
        for year in CALENDAR_YEARS
    }
    for year in CALENDAR_YEARS:
        rows.append(
            row(
                section_order=3,
                series_order=1,
                year=year,
                series_id="real_house_price_index",
                series_label="Real national house price",
                classification="untargeted_holdout",
                data_value=historical[year]["observed_real_house_price_index"],
                model_value=transition[year]["asset_price_index"],
                units="index",
                normalization="source-declared 2007=1",
                data_source="Case-Shiller CSUSHPINSA deflated by CPIAUCSL",
                model_source="transition_path.asset_price_index",
                comparability="real asset-price index",
                notes="No endpoint rescaling or price target is used.",
            )
        )
        rows.append(
            row(
                section_order=3,
                series_order=2,
                year=year,
                series_id="real_rent_cpi_proxy_index",
                series_label="Real rent CPI proxy",
                classification="untargeted_holdout_data_only",
                data_value=historical[year]["observed_real_rent_index"],
                model_value=None,
                units="index",
                normalization="source-declared 2007=1",
                data_source="rent-of-primary-residence CPI deflated by CPIAUCSL",
                model_source="none: model has no separate market-rent series",
                comparability="data context only",
                notes="Never substitute the model user-cost flow for observed rent CPI.",
            )
        )
        rows.append(
            row(
                section_order=3,
                series_order=3,
                year=year,
                series_id="price_to_rent_index",
                series_label="Price-to-rent index",
                classification="descriptive_noncomparable_diagnostic",
                data_value=historical[year]["observed_price_to_rent_index"],
                model_value=model_price_to_implicit_cost[year],
                units="index",
                normalization="2007=1; model denominator is implicit rental cost",
                data_source="CSUSHPINSA divided by rent CPI",
                model_source="asset price divided by fixed-rate implicit rental cost",
                comparability="diagnostic: model denominator is not observed rent CPI",
                notes=(
                    "The maintained user-cost identity mechanically limits this model "
                    "series. No scalar fit or RMSE claim is computed."
                ),
            )
        )

    residual_specs = (
        ("population_bridge_gap", "Population bridge gap", "population_target_gap"),
        (
            "age_bridge_max_gap",
            "Householder-age bridge maximum gap",
            "census_age_bridge_maximum_target_gap",
        ),
        ("market_clearing_residual", "Housing-market residual", "relative_market_residual"),
        ("mass_accounting_residual", "Mass-accounting residual", "mass_accounting_residual"),
    )
    for series_order, (series_id, label, key) in enumerate(residual_specs, start=1):
        for year in CALENDAR_YEARS:
            rows.append(
                row(
                    section_order=4,
                    series_order=series_order,
                    year=year,
                    series_id=series_id,
                    series_label=label,
                    classification="numerical_residual_audit_not_fit",
                    data_value=0.0,
                    model_value=transition[year][key],
                    units="relative residual or model mass",
                    normalization="none",
                    data_source="numerical identity benchmark",
                    model_source=f"transition_path.{key}",
                    comparability="numerical residual",
                    notes="Reported for numerical validation; never an empirical fit moment.",
                )
            )

    for series_order, name in enumerate(sorted(moments[2023]), start=1):
        record = moments[2023][name]
        rows.append(
            row(
                section_order=5,
                series_order=series_order,
                year=2023,
                series_id=f"terminal_target__{name}",
                series_label=name,
                classification="terminal_maintained_target",
                data_value=record["target"],
                model_value=record["model"],
                units="target-specific level",
                normalization="maintained target definition",
                data_source="active target contract in selected task summary",
                model_source="dated_model_moments.csv, 2023 model state",
                comparability="maintained target object",
                notes=(
                    f"weight={record['weight']}; "
                    f"loss_contribution={record['loss_contribution']}. {SCOPE_WARNING}"
                ),
            )
        )
    return sorted(
        rows,
        key=lambda item: (
            int(item["section_order"]),
            int(item["series_order"]),
            int(item["calendar_year"]),
        ),
    )


def sign(value: float, tolerance: float = 1e-12) -> int:
    if abs(value) <= tolerance:
        return 0
    return 1 if value > 0 else -1


def fit_statistic_row(
    series_id: str,
    label: str,
    classification: str,
    data: Sequence[float] | None,
    model: Sequence[float] | None,
    *,
    evaluation_years: Sequence[int] = CALENDAR_YEARS,
    withhold_scalar_statistics: bool = False,
    notes: str = "",
) -> dict[str, Any]:
    blank = ""
    base = {
        "series_id": series_id,
        "series_label": label,
        "classification": classification,
        "evaluation_years": ";".join(str(year) for year in evaluation_years),
        "n_common_dates": 0 if data is None or model is None else len(evaluation_years),
        "scalar_statistics_status": "computed",
        "gap_2007_model_minus_data": blank,
        "gap_2011_model_minus_data": blank,
        "gap_2023_model_minus_data": blank,
        "change_start_year": blank,
        "change_end_year": blank,
        "change_gap_model_minus_data": blank,
        "rmse_five_dates": blank,
        "rmse_evaluation_window": blank,
        "interval_sign_agreement_fraction": blank,
        "interval_sign_agreement_count": blank,
        "interval_count": blank,
        "notes": notes,
    }
    if data is None or model is None:
        base["scalar_statistics_status"] = "not_available_data_only"
        return base
    if len(data) != len(evaluation_years) or len(model) != len(evaluation_years):
        raise ValueError(f"Fit statistic {series_id} has a date/value length mismatch")
    if len(evaluation_years) < 2 or tuple(sorted(evaluation_years)) != tuple(
        evaluation_years
    ):
        raise ValueError(f"Fit statistic {series_id} has an invalid evaluation window")
    data_array = np.asarray(data, dtype=float)
    model_array = np.asarray(model, dtype=float)
    if not np.all(np.isfinite(data_array)) or not np.all(np.isfinite(model_array)):
        raise ValueError(f"Fit statistic {series_id} contains nonfinite values")
    if withhold_scalar_statistics:
        base["scalar_statistics_status"] = "withheld_noncomparable_measurement"
        return base
    gaps = model_array - data_array
    matches = [
        sign(model_array[index + 1] - model_array[index])
        == sign(data_array[index + 1] - data_array[index])
        for index in range(len(evaluation_years) - 1)
    ]
    gap_by_year = dict(zip(evaluation_years, gaps))
    base.update(
        {
            "gap_2007_model_minus_data": gap_by_year.get(2007, blank),
            "gap_2011_model_minus_data": gap_by_year.get(2011, blank),
            "gap_2023_model_minus_data": gap_by_year.get(2023, blank),
            "change_start_year": evaluation_years[0],
            "change_end_year": evaluation_years[-1],
            "change_gap_model_minus_data": (
                (model_array[-1] - model_array[0])
                - (data_array[-1] - data_array[0])
            ),
            "rmse_evaluation_window": float(np.sqrt(np.mean(gaps * gaps))),
            "interval_sign_agreement_fraction": sum(matches) / len(matches),
            "interval_sign_agreement_count": sum(matches),
            "interval_count": len(matches),
        }
    )
    if tuple(evaluation_years) == CALENDAR_YEARS:
        base["rmse_five_dates"] = base["rmse_evaluation_window"]
    return base


def build_fit_statistics(
    transition: Mapping[int, Mapping[str, float]],
    moments: Mapping[int, Mapping[str, Mapping[str, Any]]],
    fertility: Mapping[int, Mapping[str, float]],
    standardized_timing: Mapping[int, Mapping[str, float]],
    acs: Mapping[int, Mapping[str, float]],
    historical: Mapping[int, Mapping[str, float]],
    nchs: Mapping[int, Mapping[str, float]],
) -> list[dict[str, Any]]:
    def values(source: Mapping[int, Mapping[str, float]], key: str) -> list[float]:
        return [float(source[year][key]) for year in CALENDAR_YEARS]

    moment_rooms = [float(moments[year][ROOM_MOMENT]["model"]) for year in CALENDAR_YEARS]
    initial_cost = transition[2007]["housing_user_cost"]
    implicit_cost = [transition[year]["housing_user_cost"] / initial_cost for year in CALENDAR_YEARS]
    model_price_to_cost = [
        transition[year]["asset_price_index"] / implicit_cost[index]
        for index, year in enumerate(CALENDAR_YEARS)
    ]
    specs: list[
        tuple[
            str,
            str,
            str,
            Sequence[float] | None,
            Sequence[float] | None,
            Sequence[int],
            bool,
            str,
        ]
    ] = [
        (
            "household_population_index",
            "Adult-household population",
            "imposed_bridge_matched_by_construction_not_fit",
            values(historical, "observed_population_index"),
            values(transition, "population_index"),
            CALENDAR_YEARS,
            False,
            "Accounting check only; never interpret as fit.",
        ),
        (
            "period_tfr_explicit_states",
            "Period fertility rate, explicit states",
            "descriptive_noncomparable_holdout",
            values(historical, "observed_tfr"),
            values(fertility, "period_tfr_explicit_states"),
            CALENDAR_YEARS,
            True,
            (
                "Scalar gap, change, RMSE, and sign claims are withheld because the "
                "model denominator is adult-household mass rather than female exposure."
            ),
        ),
        (
            "period_tfr_topcode_adjusted_sensitivity",
            "Period fertility rate, topcoded sensitivity",
            "descriptive_noncomparable_holdout_sensitivity",
            values(historical, "observed_tfr"),
            values(fertility, "period_tfr_topcode_adjusted"),
            CALENDAR_YEARS,
            True,
            (
                "Scalar claims are withheld for the denominator mismatch; this "
                "descriptive sensitivity assigns the observed mean to the 3+ parity bin."
            ),
        ),
        (
            "period_first_birth_mean_age",
            "Mean age among period first births",
            "untargeted_holdout",
            values(nchs, "period_first_birth_mean_age"),
            values(standardized_timing, "period_first_birth_mean_age"),
            CALENDAR_YEARS,
            False,
            (
                "Model age-specific first-birth rates are standardized with national "
                "ACS female exposure shares; this is not the targeted cohort object."
            ),
        ),
        (
            "period_first_birth_share_age30plus",
            "Period first births at age 30+",
            "untargeted_holdout",
            values(nchs, "period_first_birth_share_age30plus"),
            values(standardized_timing, "period_first_birth_share_age30plus"),
            CALENDAR_YEARS,
            False,
            (
                "Model age-specific first-birth rates are standardized with national "
                "ACS female exposure shares; this is not the targeted cohort object."
            ),
        ),
        (
            "aggregate_ownership_rate_18_85",
            "Homeownership, heads ages 18--85",
            "untargeted_holdout",
            values(acs, "ownership_rate"),
            values(transition, "owner_rate"),
            CALENDAR_YEARS,
            False,
            "Separate from the prime-age ownership calibration moment.",
        ),
        (
            "mean_rooms_literal_18_85",
            "Mean occupied rooms, heads ages 18--85",
            "untargeted_holdout",
            [float(acs[year]["mean_rooms_literal"]) for year in CALENDAR_YEARS[1:]],
            moment_rooms[1:],
            CALENDAR_YEARS[1:],
            False,
            (
                "Longitudinal statistics use 2011--2023 only. The isolated 2007 "
                "cross-section uses a 9+ top category and is excluded."
            ),
        ),
        (
            "real_house_price_index",
            "Real national house price",
            "untargeted_holdout",
            values(historical, "observed_real_house_price_index"),
            values(transition, "asset_price_index"),
            CALENDAR_YEARS,
            False,
            "Case-Shiller deflated by CPI; no endpoint refit.",
        ),
        (
            "real_rent_cpi_proxy_index",
            "Real rent CPI proxy",
            "untargeted_holdout_data_only",
            values(historical, "observed_real_rent_index"),
            None,
            CALENDAR_YEARS,
            False,
            "The model has no separate rent series; user cost is not substituted.",
        ),
        (
            "price_to_rent_index",
            "Price-to-rent index",
            "descriptive_noncomparable_diagnostic",
            values(historical, "observed_price_to_rent_index"),
            model_price_to_cost,
            CALENDAR_YEARS,
            True,
            (
                "Scalar claims are withheld: the observed denominator is rent CPI, "
                "whereas the model denominator is implicit rental cost."
            ),
        ),
    ]
    rows = [
        fit_statistic_row(
            series_id,
            label,
            classification,
            data,
            model,
            evaluation_years=evaluation_years,
            withhold_scalar_statistics=withhold,
            notes=notes,
        )
        for (
            series_id,
            label,
            classification,
            data,
            model,
            evaluation_years,
            withhold,
            notes,
        ) in specs
    ]
    residual_specs = (
        ("population_bridge_gap", "Population bridge gap", "population_target_gap"),
        (
            "age_bridge_max_gap",
            "Householder-age bridge maximum gap",
            "census_age_bridge_maximum_target_gap",
        ),
        ("market_clearing_residual", "Housing-market residual", "relative_market_residual"),
        ("mass_accounting_residual", "Mass-accounting residual", "mass_accounting_residual"),
    )
    zero = [0.0] * 5
    for series_id, label, key in residual_specs:
        rows.append(
            fit_statistic_row(
                series_id,
                label,
                "numerical_residual_audit_not_fit",
                zero,
                values(transition, key),
                notes="Numerical identity check; not an empirical fit statistic.",
            )
        )
    return rows


def configure_axes(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xticks(CALENDAR_YEARS)
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6, alpha=0.7)


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> list[str]:
    paths = []
    for suffix in ("png", "pdf"):
        path = output_dir / f"{stem}.{suffix}"
        fig.savefig(path, dpi=220 if suffix == "png" else None)
        paths.append(str(path.resolve()))
    plt.close(fig)
    return paths


def make_figures(
    output_dir: Path,
    transition: Mapping[int, Mapping[str, float]],
    moments: Mapping[int, Mapping[str, Mapping[str, Any]]],
    fertility: Mapping[int, Mapping[str, float]],
    standardized_timing: Mapping[int, Mapping[str, float]],
    acs: Mapping[int, Mapping[str, float]],
    historical: Mapping[int, Mapping[str, float]],
    nchs: Mapping[int, Mapping[str, float]],
) -> dict[str, list[str]]:
    years = np.asarray(CALENDAR_YEARS)
    colors = {"data": "#111111", "model": "#2f6f9f", "sensitivity": "#c56a1a"}
    figures: dict[str, list[str]] = {}

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.5))
    axes[0].plot(
        years,
        [historical[y]["observed_population_index"] for y in CALENDAR_YEARS],
        "o-",
        color=colors["data"],
        label="Imposed Census/ACS path",
    )
    axes[0].plot(
        years,
        [transition[y]["population_index"] for y in CALENDAR_YEARS],
        "s--",
        color=colors["model"],
        label="Model after bridge",
    )
    axes[0].set_title("Adult-household population index")
    axes[0].set_ylabel("2007 = 1")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(
        years,
        [
            max(abs(transition[y]["population_target_gap"]), 1e-18)
            for y in CALENDAR_YEARS
        ],
        "o-",
        label="Total household gap",
    )
    axes[1].plot(
        years,
        [
            max(abs(transition[y]["census_age_bridge_maximum_target_gap"]), 1e-18)
            for y in CALENDAR_YEARS
        ],
        "s--",
        label="Maximum age-cell gap",
    )
    axes[1].set_yscale("log")
    axes[1].set_title("Bridge accounting residuals")
    axes[1].set_ylabel("Absolute residual")
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        configure_axes(ax)
    fig.suptitle("Imposed historical inputs: matched by construction, not model fit", fontsize=12)
    figures["imposed_inputs_matched_by_construction"] = save_figure(
        fig, output_dir, "imposed_inputs_matched_by_construction"
    )

    fig, axes = plt.subplots(2, 3, figsize=(12.3, 6.7), sharex=True)
    panels = axes.ravel()
    panels[0].plot(
        years,
        [historical[y]["observed_tfr"] for y in CALENDAR_YEARS],
        "o-",
        color=colors["data"],
        label="Data",
    )
    panels[0].plot(
        years,
        [fertility[y]["period_tfr_explicit_states"] for y in CALENDAR_YEARS],
        "s--",
        color=colors["model"],
        label="Model: explicit states",
    )
    panels[0].plot(
        years,
        [fertility[y]["period_tfr_topcode_adjusted"] for y in CALENDAR_YEARS],
        ":",
        marker="^",
        color=colors["sensitivity"],
        label="Model: 3+ sensitivity",
    )
    panels[0].set_title("Period fertility rate")
    panels[0].legend(frameon=False, fontsize=7)
    panels[1].plot(
        years,
        [nchs[y]["period_first_birth_mean_age"] for y in CALENDAR_YEARS],
        "o-",
        color=colors["data"],
        label="NCHS",
    )
    panels[1].plot(
        years,
        [standardized_timing[y]["period_first_birth_mean_age"] for y in CALENDAR_YEARS],
        "s--",
        color=colors["model"],
        label="Model",
    )
    panels[1].set_title("Mean age among period first births")
    panels[1].set_ylabel("Years")
    panels[1].legend(frameon=False, fontsize=7)
    panels[2].plot(
        years,
        [nchs[y]["period_first_birth_share_age30plus"] for y in CALENDAR_YEARS],
        "o-",
        color=colors["data"],
    )
    panels[2].plot(
        years,
        [
            standardized_timing[y]["period_first_birth_share_age30plus"]
            for y in CALENDAR_YEARS
        ],
        "s--",
        color=colors["model"],
    )
    panels[2].set_title("Period first births at age 30+")
    panels[2].set_ylabel("Share")
    panels[3].plot(
        years,
        [acs[y]["ownership_rate"] for y in CALENDAR_YEARS],
        "o-",
        color=colors["data"],
        label="ACS",
    )
    panels[3].plot(
        years,
        [transition[y]["owner_rate"] for y in CALENDAR_YEARS],
        "s--",
        color=colors["model"],
        label="Model",
    )
    panels[3].set_title("Ownership, heads ages 18--85")
    panels[3].set_ylabel("Share")
    panels[3].legend(frameon=False, fontsize=7)
    rooms_years = np.asarray(CALENDAR_YEARS[1:])
    panels[4].plot(
        rooms_years,
        [acs[y]["mean_rooms_literal"] for y in CALENDAR_YEARS[1:]],
        "o-",
        color=colors["data"],
        label="ACS, 2011--23",
    )
    panels[4].plot(
        rooms_years,
        [moments[y][ROOM_MOMENT]["model"] for y in CALENDAR_YEARS[1:]],
        "s--",
        color=colors["model"],
        label="Model, 2011--23",
    )
    panels[4].scatter(
        [2007],
        [acs[2007]["mean_rooms_literal"]],
        marker="x",
        s=45,
        color=colors["data"],
        label="2007 ACS: cross-section only",
        zorder=4,
    )
    panels[4].scatter(
        [2007],
        [moments[2007][ROOM_MOMENT]["model"]],
        marker="x",
        s=45,
        color=colors["model"],
        zorder=4,
    )
    panels[4].axvspan(2006.4, 2007.6, color="#eeeeee", zorder=0)
    panels[4].text(
        2007,
        0.5,
        "2007: ROOMS 9 = 9+\ncross-section only",
        transform=panels[4].get_xaxis_transform(),
        ha="center",
        va="center",
        rotation=90,
        fontsize=6.5,
        color="#666666",
    )
    panels[4].set_title("Mean rooms: 2011--23 longitudinal")
    panels[4].set_ylabel("Literal rooms")
    panels[4].legend(frameon=False, fontsize=7)
    panels[5].axis("off")
    panels[5].text(
        0.02,
        0.92,
        "Classification\n\n"
        "Period fertility and timing: untargeted.\n"
        "Aggregate ownership: untargeted.\n"
        "National ACS rooms path: untargeted holdout.\n\n"
        "Model timing is standardized using\n"
        "ACS female age-exposure shares.\n"
        "NCHS and model timing use the same\n"
        "boundary-collapsed midpoint operator.\n"
        "Rooms: 2007 is visually isolated and\n"
        "excluded from longitudinal statistics.\n"
        "Period-TFR comparison is descriptive:\n"
        "the model denominator is adult-household mass,\n"
        "not literal female exposure; scalar fit is withheld.\n\n"
        "Scope warning: national paths and holdouts;\n"
        "the 12 maintained targets are mixed-scope.",
        va="top",
        fontsize=9,
        transform=panels[5].transAxes,
    )
    for ax in panels[:5]:
        configure_axes(ax)
    fig.suptitle("Historical fertility and housing validation", fontsize=12)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.09, top=0.88, wspace=0.28, hspace=0.35)
    figures["untargeted_fertility_housing_validation"] = save_figure(
        fig, output_dir, "untargeted_fertility_housing_validation"
    )

    initial_cost = transition[2007]["housing_user_cost"]
    cost_index = {
        year: transition[year]["housing_user_cost"] / initial_cost
        for year in CALENDAR_YEARS
    }
    fig, axes = plt.subplots(1, 3, figsize=(12.3, 3.6))
    axes[0].plot(
        years,
        [historical[y]["observed_real_house_price_index"] for y in CALENDAR_YEARS],
        "o-",
        color=colors["data"],
        label="Real Case-Shiller",
    )
    axes[0].plot(
        years,
        [transition[y]["asset_price_index"] for y in CALENDAR_YEARS],
        "s--",
        color=colors["model"],
        label="Model asset price",
    )
    axes[0].set_title("Real house price")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].plot(
        years,
        [historical[y]["observed_real_rent_index"] for y in CALENDAR_YEARS],
        "o-",
        color=colors["data"],
        label="Real rent CPI proxy",
    )
    axes[1].set_title("Real rent CPI proxy (data only)")
    axes[1].legend(frameon=False, fontsize=7)
    axes[1].text(
        0.03,
        0.08,
        "No separate model rent series;\nuser cost is not relabeled as rent.",
        transform=axes[1].transAxes,
        fontsize=8,
    )
    axes[2].plot(
        years,
        [historical[y]["observed_price_to_rent_index"] for y in CALENDAR_YEARS],
        "o-",
        color=colors["data"],
        label="Observed price / rent CPI",
    )
    axes[2].plot(
        years,
        [transition[y]["asset_price_index"] / cost_index[y] for y in CALENDAR_YEARS],
        "s--",
        color=colors["model"],
        label="Model price / implicit cost",
    )
    axes[2].set_title("Price ratios (descriptive; no scalar fit)")
    axes[2].set_ylim(0.69, 1.10)
    axes[2].legend(frameon=False, fontsize=7, loc="upper right")
    for ax in axes:
        configure_axes(ax)
        ax.set_ylabel("2007 = 1")
    fig.suptitle("Untargeted price and rent validation", fontsize=12)
    figures["untargeted_price_rent_validation"] = save_figure(
        fig, output_dir, "untargeted_price_rent_validation"
    )

    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    residuals = (
        ("Population bridge", "population_target_gap", "o-"),
        ("Age bridge maximum", "census_age_bridge_maximum_target_gap", "s--"),
        ("Housing market", "relative_market_residual", "^-"),
        ("Mass accounting", "mass_accounting_residual", "d:"),
    )
    for label, key, style in residuals:
        ax.plot(
            years,
            [max(abs(transition[y][key]), 1e-18) for y in CALENDAR_YEARS],
            style,
            label=label,
        )
    ax.set_yscale("log")
    ax.set_ylabel("Absolute residual (floor shown at $10^{-18}$)")
    ax.set_title("Numerical residual audit; these are not empirical moments")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    configure_axes(ax)
    figures["numerical_residual_audit"] = save_figure(
        fig, output_dir, "numerical_residual_audit"
    )
    return figures


def input_receipt(path: Path) -> dict[str, Any]:
    path = path.resolve()
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)}


def build_packet(args: argparse.Namespace) -> dict[str, Any]:
    case_dir = args.case_dir.resolve()
    if not case_dir.is_dir():
        raise FileNotFoundError(case_dir)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_path = find_case_summary(case_dir)
    summary = read_json(summary_path)
    model_fingerprints = validate_case_summary(summary, summary_path)
    selected_candidate = validate_selected_case(summary, case_dir)
    transition, transition_path = load_transition(case_dir)
    moments, moments_path = load_dated_moments(case_dir, selected_candidate)
    target_loss_contract = validate_target_and_loss_contract(summary, moments)
    fertility, fertility_path = load_dated_fertility(case_dir)
    acs, acs_metadata, acs_metadata_path, acs_verification_path = load_acs_housing(
        args.acs_housing_csv
    )
    standardized_timing, fertility_age_path, timing_provenance = (
        load_standardized_first_birth_timing(
            case_dir,
            args.acs_female_exposure_csv,
            fertility,
            acs_metadata,
        )
    )
    historical, historical_metadata, historical_contract = load_historical_observations(
        args.historical_comparison_csv, args.historical_metadata_json
    )
    population_contract = validate_population_bridge_contract(
        summary, transition, historical
    )
    nchs, nchs_provenance = load_nchs_first_birth_timing(
        args.nchs_counts_csv, args.nchs_manifest_csv, args.nchs_metadata_json
    )

    comparison_rows = build_comparison_rows(
        transition, moments, fertility, standardized_timing, acs, historical, nchs
    )
    fit_statistics = build_fit_statistics(
        transition, moments, fertility, standardized_timing, acs, historical, nchs
    )
    comparison_path = output_dir / "five_date_model_data_comparison.csv"
    fit_path = output_dir / "five_date_fit_statistics.csv"
    write_csv(comparison_path, comparison_rows)
    write_csv(fit_path, fit_statistics)
    figures = make_figures(
        output_dir,
        transition,
        moments,
        fertility,
        standardized_timing,
        acs,
        historical,
        nchs,
    )

    inputs = {
        "case_summary": input_receipt(summary_path),
        "transition_path": input_receipt(transition_path),
        "dated_model_moments": input_receipt(moments_path),
        "dated_period_fertility": input_receipt(fertility_path),
        "dated_period_fertility_by_age": input_receipt(fertility_age_path),
        "acs_housing_path": input_receipt(args.acs_housing_csv),
        "acs_female_exposure_path": input_receipt(args.acs_female_exposure_csv),
        "acs_metadata": input_receipt(acs_metadata_path),
        "acs_verification_report": input_receipt(acs_verification_path),
        "historical_comparison": input_receipt(args.historical_comparison_csv),
        "historical_metadata": input_receipt(args.historical_metadata_json),
        "historical_comparison_summary": input_receipt(
            Path(historical_contract["comparison_summary_path"])
        ),
        "nchs_first_birth_counts": input_receipt(args.nchs_counts_csv),
        "nchs_first_birth_manifest": input_receipt(args.nchs_manifest_csv),
        "nchs_first_birth_metadata": input_receipt(args.nchs_metadata_json),
    }
    for name, receipt in historical_contract["raw_series"].items():
        inputs[f"historical_raw__{name}"] = input_receipt(Path(receipt["path"]))
    bundle_payload = "\n".join(
        f"{name}:{receipt['sha256']}" for name, receipt in sorted(inputs.items())
    )
    input_bundle_sha256 = hashlib.sha256(bundle_payload.encode()).hexdigest()
    provenance = {
        "status": "complete_e5f_2007_2023_historical_validation_packet",
        "calendar_years": list(CALENDAR_YEARS),
        "scope_warning": SCOPE_WARNING,
        "input_bundle_sha256": input_bundle_sha256,
        "inputs": inputs,
        "model_contract": {
            **model_fingerprints,
            "selected_candidate": selected_candidate,
            "target_set": summary["target_set"],
            "task_status": summary["status"],
            "population_bridge_status": summary["population_bridge"]["status"],
            "population_validation_status": summary["population_validation_status"],
            "dated_measurement_contract": {
                "period_fertility": (
                    "untargeted current-rate diagnostic with adult-household denominator"
                ),
                "period_first_birth_timing": (
                    "untargeted age-specific-rate diagnostic standardized to ACS female exposure"
                ),
            },
            "renewal_accounting_contract": summary["renewal_accounting_contract"],
            "renewal_accounting_old_state": summary["renewal_accounting_old_state"],
            "target_and_loss_reconciliation": target_loss_contract,
            "no_policy_historical_path": True,
        },
        "data_contracts": {
            "acs": {
                "contract_sha256": acs_metadata["contract_sha256"],
                "source_sha256": acs_metadata["source"]["sha256"],
                "builder_sha256": acs_metadata["builder"]["sha256"],
                "sample": (
                    "national ACS; GQ in {1,2}; PERNUM=1; RELATE=1; HHWT>0; "
                    "householder age 18--85"
                ),
                "ownership": "OWNERSHP in {1,2}; owner=1",
                "rooms": "literal positive ROOMS; code 9 is treated as 9",
                "female_exposure": timing_provenance,
            },
            "historical_series": {
                "metadata": {
                    name: historical_metadata["series"][name]
                    for name in (
                        "consumer_price",
                        "fertility",
                        "house_price",
                        "rent",
                    )
                },
                "receipt_reconstruction": historical_contract,
            },
            "nchs_first_birth_timing": nchs_provenance,
        },
        "classification_contract": {
            "imposed_bridge": {
                "series": [
                    "household_population_index",
                    "householder_age_bridge_max_gap",
                ],
                "interpretation": "matched by construction; excluded from fit claims",
                "receipt": population_contract,
            },
            "terminal_maintained_targets": {
                "year": 2023,
                "series": [f"terminal_target__{name}" for name in sorted(moments[2023])],
                "interpretation": (
                    "the selected task's active mixed-scope SMM target system; "
                    "not a national time-series fit"
                ),
            },
            "untargeted_holdouts": {
                "series": [
                    "period_tfr_explicit_states",
                    "period_first_birth_mean_age",
                    "period_first_birth_share_age30plus",
                    "aggregate_ownership_rate_18_85",
                    "mean_rooms_literal_18_85",
                    "real_house_price_index",
                    "real_rent_cpi_proxy_index",
                    "price_to_rent_index",
                ],
                "interpretation": "never included retroactively in the transition objective",
            },
            "numerical_residuals": {
                "interpretation": "solver/accounting checks, not empirical moments"
            },
            "scope_warning": SCOPE_WARNING,
        },
        "measurement_caveats": {
            "period_tfr": (
                "Main model series uses explicit parity-state birth flows; the 3+ "
                "topcode adjustment is a sensitivity. The model denominator is "
                "pre-fertility adult-household mass, not literal national female exposure."
            ),
            "first_birth_timing": (
                "NCHS period first-birth timing is compared with model age-specific "
                "first-birth rates standardized to national ACS female exposure shares; "
                "it is not the synthetic-cohort timing object targeted at 2023."
            ),
            "ownership": (
                "The holdout uses aggregate household heads ages 18--85 and the "
                "transition's aggregate owner rate, not the maintained prime-age target."
            ),
            "rooms": (
                "This independently built national ACS literal-room path is an untargeted "
                "holdout. The 2007 code 9 means 9+, unlike later annual samples; 2007 "
                "is cross-section only and longitudinal statistics start in 2011."
            ),
            "rent": (
                "Observed rent is a rent-of-primary-residence CPI proxy. The model "
                "has no separate market-rent series, so its user-cost flow is never relabeled as rent."
            ),
            "price_to_rent": (
                "The observed denominator is rent CPI; the model denominator is its "
                "implicit rental-cost flow under the maintained user-cost identity. "
                "Scalar fit statistics are withheld."
            ),
            "scope": SCOPE_WARNING,
        },
        "normalization_contract": {
            "no_posthoc_endpoint_rescaling": True,
            "allowed_indices": (
                "source-declared 2007=1 population, Case-Shiller, real rent CPI proxy, "
                "and price/rent indices; model asset-price index is already 2007=1; "
                "model implicit rental cost is indexed to 2007 solely for the "
                "price-to-implicit-cost diagnostic"
            ),
        },
        "outputs": {
            "comparison_csv": input_receipt(comparison_path),
            "fit_statistics_csv": input_receipt(fit_path),
            "figures": {
                name: [input_receipt(Path(path)) for path in paths]
                for name, paths in figures.items()
            },
        },
        "validation_gates": {
            "five_calendar_dates_exact": True,
            "model_task_complete": True,
            "no_policy_path_confirmed": True,
            "renewal_2p1_inverse_20year_algebra_reconciled": True,
            "target_count_fingerprint_and_loss_reconciled": True,
            "model_and_data_fingerprints_present": True,
            "historical_values_reconstructed_from_receipts": True,
            "observed_population_equals_imposed_population_target": True,
            "acs_receipts_reproduce_rates": True,
            "nchs_metadata_hash_contract_verified": True,
            "nchs_boundary_collapsed_midpoint_operator_verified": True,
            "period_tfr_reconstructed_from_age_specific_rates": True,
            "first_birth_timing_standardized_to_acs_female_exposure": True,
            "terminal_target_accounting_reproduced": True,
            "imposed_paths_not_classified_as_fit": True,
            "rent_cpi_not_replaced_by_user_cost": True,
            "noncomparable_scalar_fit_claims_withheld": True,
            "rooms_2007_excluded_from_longitudinal_statistics": True,
        },
    }
    provenance_path = output_dir / "classification_and_provenance.json"
    write_json(provenance_path, provenance)
    return provenance


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    result = build_packet(args)
    print(json.dumps({"status": result["status"], "outputs": result["outputs"]}, indent=2))


if __name__ == "__main__":
    main()
