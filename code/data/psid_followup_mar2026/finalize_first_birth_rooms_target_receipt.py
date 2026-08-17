#!/usr/bin/env python3
"""Validate and fingerprint the corrected PSID first-birth rooms target."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
DEFAULT_OUTPUT = HERE / "output/sa_rooms_first_birth_household_aligned_v1"
DEFAULT_SOURCE = PROJECT.parent / "PSID/PSIDSHELF_MOBILITY.dta"
DEFAULT_DO_FILE = HERE / "sa_rooms_first_birth_household_aligned_v1.do"
CONTRACT_ID = "psid_first_birth_rooms_household_aligned_sa_v1_20260817"
EVENT_DISPLAY_TOL = 1.0e-6


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"Empty CSV: {path}")
    return rows


def number(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"{label} is not numeric") from error
    if not math.isfinite(result):
        raise RuntimeError(f"{label} is not finite")
    return result


def integer(value: Any, label: str) -> int:
    result = number(value, label)
    if abs(result - round(result)) > 1e-9:
        raise RuntimeError(f"{label} is not an integer")
    return int(round(result))


def build_contract(output: Path, source: Path, do_file: Path) -> dict[str, Any]:
    receipt_path = output / "target_receipt.csv"
    estimates_path = output / "event_study_estimates.csv"
    log_path = output / "sa_rooms_first_birth_household_aligned_v1.log"
    for path in (source, do_file, receipt_path, estimates_path, log_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    receipt_rows = read_rows(receipt_path)
    if len(receipt_rows) != 1:
        raise RuntimeError("Target receipt must contain exactly one row")
    receipt = receipt_rows[0]
    if receipt.get("moment") != "housing_increment_0to1":
        raise RuntimeError("Target receipt has the wrong moment")
    if receipt.get("status") != "corrected_primary_target":
        raise RuntimeError("Target receipt is not marked corrected primary")
    estimate = number(receipt.get("estimate"), "target estimate")
    standard_error = number(receipt.get("standard_error"), "target standard error")
    if standard_error <= 0.0:
        raise RuntimeError("Target standard error must be positive")
    if integer(receipt.get("contrast_start_time"), "contrast_start_time") != -1:
        raise RuntimeError("Target contrast must begin at event time -1")
    if integer(receipt.get("contrast_end_time"), "contrast_end_time") != 3:
        raise RuntimeError("Target contrast must end at event time +3")
    if integer(
        receipt.get("regression_omitted_time"), "regression_omitted_time"
    ) != -2:
        raise RuntimeError("Sun--Abraham regression must omit event time -2")
    observations = integer(
        receipt.get("estimation_observations"), "estimation observations"
    )
    household_years = integer(
        receipt.get("estimation_household_years"), "estimation household-years"
    )
    individuals = integer(
        receipt.get("estimation_individuals"), "estimation individuals"
    )
    if observations <= 0 or individuals <= 0 or observations != household_years:
        raise RuntimeError("Estimation sample is not one row per household-year")
    if receipt.get("weighting") != "PSID longitudinal pweight IW":
        raise RuntimeError("Target receipt changed the declared weighting")
    if receipt.get("control_group") != (
        "confirmed zero-child women from full relationship history"
    ):
        raise RuntimeError("Target receipt changed the control group")
    if integer(
        receipt.get("single_fu_only"),
        "single-family-unit restriction",
    ) != 1:
        raise RuntimeError("Target receipt does not impose the single-FID dwelling rule")
    if receipt.get("fertility_timing") != (
        "first biological child across RELCHI1-20 TYPE/BYEAR records"
    ):
        raise RuntimeError("Target receipt changed the biological-birth timing rule")
    confirmed_never_ids = integer(
        receipt.get("est_confirmed_never_ids"),
        "confirmed never-treated individuals",
    )
    excluded_multi_fu = integer(
        receipt.get("excl_multi_fu_hhyears"),
        "excluded multi-FU household-years",
    )
    excluded_known_parent = integer(
        receipt.get("excl_untimed_parent_ids"),
        "excluded known parents without biological timing",
    )
    excluded_unknown_history = integer(
        receipt.get("excl_unknown_history_ids"),
        "excluded unknown child histories",
    )
    if min(
        confirmed_never_ids,
        excluded_multi_fu,
        excluded_known_parent,
        excluded_unknown_history,
    ) <= 0:
        raise RuntimeError("Target receipt lacks the audited control/unit exclusions")

    estimates = read_rows(estimates_path)
    expected_event_times = set(range(-7, 12))
    observed_event_times = {
        integer(row["relative_time"], "relative_time") for row in estimates
    }
    if observed_event_times != expected_event_times:
        raise RuntimeError(
            "Full estimates do not contain the complete collapsed event-time grid"
        )
    row_p3 = [row for row in estimates if integer(row["relative_time"], "relative_time") == 3]
    if len(row_p3) != 1:
        raise RuntimeError("Full estimates must contain exactly one +3 row")
    row_m1 = [
        row
        for row in estimates
        if integer(row["relative_time"], "relative_time") == -1
    ]
    if len(row_m1) != 1:
        raise RuntimeError("Full estimates must contain exactly one -1 row")
    component_p3 = number(receipt.get("component_l3"), "receipt +3 component")
    component_m1 = number(receipt.get("component_f1"), "receipt -1 component")
    covariance = number(receipt.get("covariance_l3_f1"), "receipt +3/-1 covariance")
    if not math.isclose(
        number(row_p3[0]["b"], "+3 estimate"),
        component_p3,
        rel_tol=0,
        abs_tol=EVENT_DISPLAY_TOL,
    ):
        raise RuntimeError("Receipt +3 component differs from the full estimate")
    if not math.isclose(
        number(row_m1[0]["b"], "-1 estimate"),
        component_m1,
        rel_tol=0,
        abs_tol=EVENT_DISPLAY_TOL,
    ):
        raise RuntimeError("Receipt -1 component differs from the full estimate")
    implied_estimate = component_p3 - component_m1
    implied_variance = (
        number(row_p3[0]["se"], "+3 SE") ** 2
        + number(row_m1[0]["se"], "-1 SE") ** 2
        - 2.0 * covariance
    )
    if implied_variance <= 0.0:
        raise RuntimeError("The +3 minus -1 contrast has nonpositive variance")
    if not math.isclose(implied_estimate, estimate, rel_tol=0, abs_tol=5e-12):
        raise RuntimeError("Receipt estimate differs from the +3 minus -1 contrast")
    if not math.isclose(
        math.sqrt(implied_variance),
        standard_error,
        rel_tol=0,
        abs_tol=EVENT_DISPLAY_TOL,
    ):
        raise RuntimeError("Receipt SE does not reproduce from the covariance matrix")
    row_m6 = [
        row
        for row in estimates
        if integer(row["relative_time"], "relative_time") == -6
    ]
    if len(row_m6) != 1 or number(row_m6[0]["se"], "-6 SE") <= 0.0:
        raise RuntimeError("Full estimates omit the explicit F6 event coefficient")
    baseline_rows = [
        row
        for row in estimates
        if math.isclose(number(row["b"], "event coefficient"), 0.0, abs_tol=0.0)
        and math.isclose(number(row["se"], "event SE"), 0.0, abs_tol=0.0)
    ]
    if len(baseline_rows) != 1 or integer(
        baseline_rows[0]["relative_time"], "baseline relative_time"
    ) != -2:
        raise RuntimeError("K=-2 is not the unique omitted event time")

    return {
        "schema_version": 1,
        "contract_id": CONTRACT_ID,
        "status": "complete_provisional_model_mapping_target",
        "measurement_status": (
            "empirical contrast verified; structural mapping remains provisional "
            "because earlier leads reject a flat prepath"
        ),
        "target": {
            "moment": "housing_increment_0to1",
            "estimate": estimate,
            "standard_error": standard_error,
            "inverse_variance_weight": 1.0 / standard_error**2,
            "contrast_start_time": -1,
            "contrast_end_time": 3,
            "regression_omitted_time": -2,
            "component_l3": component_p3,
            "component_f1": component_m1,
            "covariance_l3_f1": covariance,
        },
        "sample": {
            "definition": receipt["sample"],
            "input_observations": integer(receipt["input_observations"], "input observations"),
            "input_individuals": integer(receipt["input_individuals"], "input individuals"),
            "treated_individuals": integer(receipt["treated_individuals"], "treated individuals"),
            "estimation_observations": observations,
            "estimation_individuals": individuals,
            "estimation_household_years": household_years,
            "one_observation_per_household_year": True,
            "single_family_unit_dwellings_only": True,
            "confirmed_never_treated_individuals": confirmed_never_ids,
            "excluded_multi_fu_household_years": excluded_multi_fu,
            "excluded_multi_fu_eligible_women": integer(
                receipt.get("excl_multi_fu_women"),
                "excluded multi-FU eligible women",
            ),
            "excluded_known_parent_without_biological_timing_individuals": excluded_known_parent,
            "excluded_unknown_child_history_individuals": excluded_unknown_history,
            "excluded_nonstandard_room_alignment_rows": integer(
                receipt.get("excl_bad_room_gap_rows"),
                "excluded nonstandard room-alignment rows",
            ),
        },
        "estimator": {
            "command": "eventstudyinteract",
            "fixed_effects_and_controls": receipt["fixed_effects"],
            "clustering": receipt["clustering"],
            "weighting": receipt["weighting"],
            "control_group": receipt["control_group"],
            "rooms_alignment": receipt["rooms_alignment"],
            "fertility_timing": receipt["fertility_timing"],
            "runtime_seconds": number(receipt["runtime_seconds"], "runtime seconds"),
        },
        "source": {
            "path": str(source.resolve()),
            "sha256": sha256(source),
            "do_file": str(do_file.resolve()),
            "do_file_sha256": sha256(do_file),
            "finalizer": str(Path(__file__).resolve()),
            "finalizer_sha256": sha256(Path(__file__).resolve()),
        },
        "outputs": {
            path.name: sha256(path)
            for path in (receipt_path, estimates_path, log_path)
        },
        "caveats": [
            "ACTUALROOMS_ is shifted forward by one observed interview before applying year-specific missing codes.",
            "Valid HHID is unavailable in 1982--1984, so those rows are excluded by the one-household-year rule.",
            "Physical dwellings containing multiple PSID family units are excluded before choosing one woman per HHID-year.",
            "The outcome is household-level; the primary sample clusters by the selected woman's stable person ID.",
            "The calibration target is the +3 minus -1 four-year contrast; the full event-study path is retained because earlier leads reject a flat prepath.",
            "The model analogue branches birth and no-birth cohorts for one four-year period; the mapping remains provisional rather than causal.",
            "The four-year contrast can include another birth and is not a mechanically pure 0-to-1 space increment.",
        ],
    }


def readme(contract: dict[str, Any]) -> str:
    target = contract["target"]
    sample = contract["sample"]
    estimator = contract["estimator"]
    return f"""# Corrected PSID first-birth rooms target

The paper-facing target is the change in occupied rooms from event time -1 to
event time +3 around a first birth. This four-year contrast is the closest data
analogue to one model period. The corrected estimate is
`{target['estimate']:.12g}` with clustered standard error
`{target['standard_error']:.12g}` (inverse-variance weight
`{target['inverse_variance_weight']:.12g}`).

The Sun--Abraham event study uses {sample['estimation_observations']:,}
woman-year observations for {sample['estimation_individuals']:,} women. Every
observation is unique at the household-year level, and the primary sample
excludes physical dwellings containing multiple PSID family units. The sample
contains current women age 18 or older who are the reference person or
spouse/partner, have observed age and positive PSID longitudinal weight, and
live in a valid single-family-unit dwelling. Women whose full relationship
history confirms zero children are the comparison group. The regression
includes person and survey-year fixed effects, age and education controls,
applies `IW` as a probability weight, and clusters by person ID.

`ACTUALROOMS_` is shifted forward by one observed interview within individual
before the vintage-specific non-room codes are removed. This repairs the
one-wave extraction error in the prior target. The earlier unweighted,
last-treated-control estimate is withdrawn.

The treatment date is the first biological child across all twenty recorded
child relationships. Known parents without a biological-birth year and people
with unknown child-count histories are excluded from the control group.
Event time -2 is the unique omitted period; event time -6 is explicitly
estimated. The reported standard error uses the full covariance between the
+3 and -1 coefficients. Earlier leads reject a flat prepath, so the complete
event-study profile remains a required diagnostic. The contrast may include
adjustment associated with a later birth and is not a mechanically pure 0-to-1
space requirement. Its use in calibration is therefore a provisional mapping
to the model's one-period birth-versus-no-birth branch, not a causal claim.

Contract ID: `{contract['contract_id']}`. Source, code, log, and output hashes
are recorded in `metadata.json`.
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--do-file", type=Path, default=DEFAULT_DO_FILE)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    contract = build_contract(output, args.source.resolve(), args.do_file.resolve())
    metadata_path = output / "metadata.json"
    readme_path = output / "README.md"
    if args.verify:
        recorded = json.loads(metadata_path.read_text(encoding="utf-8"))
        if recorded != contract:
            raise RuntimeError("Recorded metadata differs from the live receipt")
        expected_readme = readme(contract)
        if readme_path.read_text(encoding="utf-8") != expected_readme:
            raise RuntimeError("Recorded README differs from the live receipt")
        print("FIRST_BIRTH_ROOMS_RECEIPT_VERIFY_PASS")
        return
    if metadata_path.exists() or readme_path.exists():
        raise RuntimeError("Refusing to overwrite an existing finalized receipt")
    metadata_path.write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    readme_path.write_text(readme(contract), encoding="utf-8")
    print(
        "FIRST_BIRTH_ROOMS_RECEIPT_FINALIZED",
        f"estimate={contract['target']['estimate']:.12g}",
        f"se={contract['target']['standard_error']:.12g}",
    )


if __name__ == "__main__":
    main()
