"""Postprocess one completed first-birth rooms replay; never fits a regression."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parent
REPLAY = ANALYSIS.parent
OUTPUT = REPLAY / "output"
ROOT = REPLAY.parents[7]
RECEIPT = ANALYSIS / "replay_analysis.json"
WEIGHTS_OUT = ANALYSIS / "cohort_weights.csv"
EXPECTED_TARGET = 0.7202462623815278
EXPECTED_SE = 0.0852600513385958
EXPECTED_N = 49457
EXPECTED_IDS = 4112
NUMERIC_TOLERANCE = 1e-9
COVARIANCE_TOLERANCE = 1e-10
MAPPING_TOLERANCE = 1e-14


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: str | float | int) -> float:
    if value is None or str(value).strip() == "":
        return math.nan
    return float(value)


def cohort_number(value: str) -> int:
    return int(float(value.strip()))


def event_name(k: int) -> str:
    return f"F{-k}event" if k < 0 else f"L{k}event"


def stop_pending(execution: dict) -> int:
    if execution.get("status") in {"running", "starting"}:
        print(json.dumps({"status": "pending", "elapsed_seconds": execution.get("elapsed_seconds"),
                          "heartbeat_utc": execution.get("heartbeat_utc")}, indent=2))
        return 3
    return 0


def check_terminal() -> tuple[dict, dict, list[str], dict]:
    execution_path = REPLAY / "execution.json"
    execution = read_json(execution_path)
    pending = stop_pending(execution)
    if pending:
        raise PendingReplay
    manifest = read_json(REPLAY / "lead_preflight.json")
    errors: list[str] = []
    if execution.get("status") != "process_completed_pending_scientific_review":
        errors.append(f"execution status is {execution.get('status')!r}")
    if execution.get("exit_code") != 0:
        errors.append(f"Stata exit code is {execution.get('exit_code')!r}")
    if execution.get("script_sha256") != manifest.get("staged_sha256"):
        errors.append("executed script hash does not match lead-reviewed hash")
    script_hash = sha256(REPLAY / "replay_v1.do")
    if script_hash != execution.get("script_sha256"):
        errors.append("current replay script differs from executed script")
    if execution.get("elapsed_seconds", 10**9) > execution.get("cap_seconds", 0):
        errors.append("controller runtime exceeded its declared cap")
    return execution, manifest, errors, {"execution.json": sha256(execution_path),
                                       "replay_v1.do": script_hash}


class PendingReplay(Exception):
    pass


def analyze() -> tuple[dict, list[dict]]:
    execution, review, errors, hashes = check_terminal()
    names = [
        "target_receipt.csv", "event_study_estimates.csv", "reproduction_gates.csv",
        "room_code_diagnostics.csv", "cohort_interaction_coefficients.csv",
        "cohort_interaction_covariance.csv", "iw_coefficients.csv", "iw_covariance.csv",
        "estimation_cohort_event_support.csv", "eventstudyinteract_replay.ster",
        "sa_rooms_first_birth_household_aligned_v1.log", "replay_v1.log",
    ]
    paths = {name: OUTPUT / name for name in names}
    for name, path in paths.items():
        if not path.is_file():
            errors.append(f"required replay output is missing: {name}")
        else:
            hashes[name] = sha256(path)

    console = REPLAY / "stata_batch_console.log"
    if not console.is_file():
        errors.append("Stata batch console log is missing")
    else:
        hashes["stata_batch_console.log"] = sha256(console)
    if any(not p.is_file() for p in paths.values()):
        return {"status": "preanalysis_failed", "errors": errors, "input_hashes": hashes}, []

    run_log = paths["sa_rooms_first_birth_household_aligned_v1.log"].read_text(errors="replace")
    batch_log = paths["replay_v1.log"].read_text(errors="replace")
    console_text = console.read_text(errors="replace") if console.is_file() else ""
    for label, content in [("Stata fit log", run_log), ("Stata batch log", batch_log),
                           ("batch console", console_text)]:
        if re.search(r"\br\(\d+\);", content):
            errors.append(f"{label} contains a Stata return-code failure")
    if "CORRECTED_FIRST_BIRTH_ROOMS_TARGET" not in run_log:
        errors.append("fit log lacks the final target line")
    if "end of do-file" not in batch_log.lower():
        errors.append("batch log lacks normal do-file termination")

    gates = read_csv(paths["reproduction_gates.csv"])
    if len(gates) != 1:
        errors.append("reproduction gate receipt does not have exactly one row")
        gate = {}
    else:
        gate = gates[0]
        if abs(number(gate["target_observed"]) - EXPECTED_TARGET) > NUMERIC_TOLERANCE:
            errors.append("target reproduction exceeds 1e-9")
        if abs(number(gate["se_observed"]) - EXPECTED_SE) > NUMERIC_TOLERANCE:
            errors.append("standard-error reproduction exceeds 1e-9")
        if int(gate["observations_observed"]) != EXPECTED_N:
            errors.append("estimation observation count is not exactly 49,457")
        if int(gate["people_observed"]) != EXPECTED_IDS:
            errors.append("estimation person count is not exactly 4,112")
        if any(int(gate[k]) != 1 for k in
               ("target_pass", "se_pass", "observations_pass", "people_pass")):
            errors.append("one or more staged reproduction flags failed")

    target_rows = read_csv(paths["target_receipt.csv"])
    if len(target_rows) != 1:
        errors.append("target receipt must have one row")
        target = {}
    else:
        target = target_rows[0]
        if int(target["estimation_observations"]) != EXPECTED_N:
            errors.append("target receipt observation count differs from expected")
        if int(target["estimation_individuals"]) != EXPECTED_IDS:
            errors.append("target receipt individual count differs from expected")
        for field, expected in (("estimate", EXPECTED_TARGET), ("standard_error", EXPECTED_SE)):
            if abs(number(target[field]) - expected) > NUMERIC_TOLERANCE:
                errors.append(f"target receipt {field} differs beyond 1e-9")

    event_rows = read_csv(paths["event_study_estimates.csv"])
    if not any(int(number(row["relative_time"])) == -2 and number(row["b"]) == 0
               for row in event_rows):
        errors.append("original event curve lacks its explicit omitted -2 baseline")

    code_rows = read_csv(paths["room_code_diagnostics.csv"])
    code_diagnostics = {row["stage"]: {
        "observations": int(row["observations"]),
        "counts": {key.removesuffix("_count"): int(row[key]) for key in
                   ("room_code_0_count", "room_code_9_count", "room_code_98_count", "room_code_99_count")},
        "minimum_nonmissing_rooms": number(row["minimum_nonmissing_rooms"]),
    } for row in code_rows}
    if "e_sample" not in code_diagnostics:
        errors.append("room-code diagnostics lack the e(sample) row")

    # Bind coefficient and support evidence, then validate the full q x q covariance.
    coefficient_rows = read_csv(paths["cohort_interaction_coefficients.csv"])
    q = len(coefficient_rows)
    if q == 0:
        errors.append("cohort interaction coefficient export is empty")
        return {"status": "preanalysis_failed", "errors": errors, "input_hashes": hashes}, []
    by_index: dict[int, dict] = {}
    by_cohort_event: dict[tuple[int, str], dict] = {}
    for row in coefficient_rows:
        i = int(row["matrix_index"])
        g = cohort_number(row["cohort"])
        event = row["event"]
        item = {"index": i, "cohort": g, "event": event,
                "estimate": number(row["estimate"]),
                "marginal_variance": number(row["marginal_variance"]),
                "stata_name": row["stata_interaction_name"]}
        if i in by_index or (g, event) in by_cohort_event:
            errors.append("duplicate coefficient index or cohort/event mapping")
        by_index[i] = item
        by_cohort_event[g, event] = item
    if sorted(by_index) != list(range(1, q + 1)):
        errors.append("coefficient matrix indices are not contiguous from one")

    covariance = [[math.nan] * q for _ in range(q)]
    cells = 0
    for row in read_csv(paths["cohort_interaction_covariance.csv"]):
        i, j = int(row["row_index"]), int(row["column_index"])
        if not (1 <= i <= q and 1 <= j <= q) or math.isfinite(covariance[i-1][j-1]):
            errors.append("covariance index is out of range or duplicated")
            continue
        covariance[i-1][j-1] = number(row["covariance"])
        cells += 1
    if cells != q * q or any(not math.isfinite(x) for r in covariance for x in r):
        errors.append(f"full covariance has {cells} cells; expected {q*q}")
    max_symmetry_gap = 0.0
    max_diagonal_gap = 0.0
    for i in range(q):
        for j in range(i, q):
            max_symmetry_gap = max(max_symmetry_gap, abs(covariance[i][j] - covariance[j][i]))
        max_diagonal_gap = max(max_diagonal_gap,
                               abs(covariance[i][i] - by_index[i+1]["marginal_variance"]))
    if max_symmetry_gap > COVARIANCE_TOLERANCE:
        errors.append(f"full covariance symmetry gap {max_symmetry_gap:g} exceeds tolerance")
    if max_diagonal_gap > MAPPING_TOLERANCE:
        errors.append(f"full covariance diagonal does not match marginal variance map: {max_diagonal_gap:g}")

    # Aggregate e(sample) support only. A never-treated control has no relative event time.
    support_rows = read_csv(paths["estimation_cohort_event_support.csv"])
    support: dict[tuple[int, int], dict[str, float]] = {}
    for row in support_rows:
        group = row["cohort_group"].strip()
        if group == "never_treated" or row["event_time"].strip() == "":
            continue
        g, k = cohort_number(group), int(float(row["event_time"]))
        support[g, k] = {"observations": int(row["estimation_observations"]),
                         "weight": number(row["longitudinal_weight"])}
    if not support:
        errors.append("aggregate estimation cohort/event support is empty")

    if errors:
        return {"status": "preanalysis_failed", "errors": errors,
                "execution": execution, "input_hashes": hashes,
                "covariance_dimensions": [q, q],
                "covariance_max_symmetry_gap": max_symmetry_gap,
                "covariance_max_diagonal_gap": max_diagonal_gap}, []

    # Check e(b_interact) and the covariance diagonal map by matrix order.
    ncohort = len({row["cohort"] for row in by_index.values()})
    ordered_cohorts = sorted({row["cohort"] for row in by_index.values()})
    event_by_index: dict[int, str] = {}
    for row in coefficient_rows:
        idx, name = int(row["event_index"]), row["event"]
        if idx in event_by_index and event_by_index[idx] != name:
            errors.append("one event index maps to multiple event names")
        event_by_index[idx] = name
    ordered_events = [event_by_index[k] for k in sorted(event_by_index)]
    for i in range(1, q + 1):
        expected_event_index = (i - 1) // ncohort + 1
        expected_cohort_index = (i - 1) % ncohort + 1
        event = by_index[i]["event"]
        cohort = by_index[i]["cohort"]
        row = coefficient_rows[i-1]
        if int(row["event_index"]) != expected_event_index:
            errors.append("coefficient export is not event-major/cohort-minor")
        if int(row["cohort_index"]) != expected_cohort_index:
            errors.append("cohort matrix order differs from Stata ado order")
        if cohort != ordered_cohorts[expected_cohort_index-1]:
            errors.append("cohort row order differs from labeled matrix order")
        if event != ordered_events[expected_event_index-1]:
            errors.append("event column order differs from labeled matrix order")
    if errors:
        return {"status": "preanalysis_failed", "errors": errors,
                "execution": execution, "input_hashes": hashes}, []

    coefficient_index = {(x["cohort"], x["event"]): x["index"] for x in by_index.values()}
    weight_rows: list[dict] = []
    contrasts: list[dict] = []

    def make_candidate(label: str, start_k: int, end_k: int, start_event: str | None,
                       end_event: str, weight_k: int, omitted_start: bool = False) -> dict:
        start_cohorts = {g for (g, k), s in support.items() if k == start_k and s["observations"] > 0}
        end_cohorts = {g for (g, k), s in support.items() if k == end_k and s["observations"] > 0}
        common = sorted(start_cohorts & end_cohorts)
        if not common:
            raise ValueError(f"no common cohorts for {label}")
        weighted = {g: support[g, weight_k]["weight"] for g in common}
        if any(not math.isfinite(w) or w <= 0 for w in weighted.values()):
            raise ValueError(f"invalid fixed weights in {label}")
        total_weight = math.fsum(weighted.values())
        normalized = {g: weighted[g] / total_weight for g in common}
        vector = [0.0] * q
        cohort_contrasts = {}
        for g in common:
            end_key = (g, end_event)
            if end_key not in coefficient_index:
                raise ValueError(f"missing supported endpoint coefficient {end_key}")
            end_i = coefficient_index[end_key]
            end_b = by_index[end_i]["estimate"]
            if omitted_start:
                start_b = 0.0
                start_i = None
            else:
                start_key = (g, start_event)
                if start_key not in coefficient_index:
                    raise ValueError(f"missing supported start coefficient {start_key}")
                start_i = coefficient_index[start_key]
                start_b = by_index[start_i]["estimate"]
                vector[start_i-1] -= normalized[g]
            vector[end_i-1] += normalized[g]
            cohort_contrasts[g] = end_b - start_b
            weight_rows.append({"candidate": label, "cohort": g,
                                "start_event_time": start_k, "end_event_time": end_k,
                                "start_observations": support[g, start_k]["observations"],
                                "start_weight_sum": support[g, start_k]["weight"],
                                "end_observations": support[g, end_k]["observations"],
                                "end_weight_sum": support[g, end_k]["weight"],
                                "weight_endpoint": weight_k,
                                "fixed_weight_sum": weighted[g],
                                "normalized_fixed_weight": normalized[g],
                                "start_coefficient": start_b,
                                "end_coefficient": end_b,
                                "within_cohort_contrast": cohort_contrasts[g]})
        estimate = math.fsum(normalized[g] * cohort_contrasts[g] for g in common)
        variance = math.fsum(vector[i] * covariance[i][j] * vector[j]
                             for i in range(q) for j in range(q))
        # Independent cohort-contrast covariance construction checks the linear form.
        contrast_variance = math.fsum(
            normalized[g] * normalized[h] * (
                covariance[coefficient_index[g, end_event]-1][coefficient_index[h, end_event]-1]
                + (0.0 if omitted_start else covariance[coefficient_index[g, start_event]-1][coefficient_index[h, start_event]-1])
                - (0.0 if omitted_start else covariance[coefficient_index[g, end_event]-1][coefficient_index[h, start_event]-1])
                - (0.0 if omitted_start else covariance[coefficient_index[g, start_event]-1][coefficient_index[h, end_event]-1]))
            for g in common for h in common)
        linear_form_gap = abs(variance - contrast_variance)
        if variance < -1e-12 or contrast_variance < -1e-12:
            raise ValueError(f"negative contrast variance for {label}")
        if linear_form_gap > 1e-10:
            raise ValueError(f"linear-form variance mismatch for {label}: {linear_form_gap:g}")
        within_cohort_weight_sum_gap = 0.0 if omitted_start else max(
            abs((-normalized[g]) + normalized[g]) for g in common)
        candidate = {
            "name": label,
            "start_event_time": start_k,
            "end_event_time": end_k,
            "start_event_coefficient": "F2event omitted reference; normalized to zero" if omitted_start else start_event,
            "end_event_coefficient": end_event,
            "omitted_start_support_rule": "start event-time observation must be positive for every included cohort" if omitted_start else None,
            "weight_endpoint": weight_k,
            "common_cohorts": common,
            "cohort_count": len(common),
            "start_endpoint_cohort_count": len(start_cohorts),
            "end_endpoint_cohort_count": len(end_cohorts),
            "retained_start_weight_mass": math.fsum(support[g, start_k]["weight"] for g in common)
                / math.fsum(support[g, start_k]["weight"] for g in start_cohorts),
            "retained_end_weight_mass": math.fsum(support[g, end_k]["weight"] for g in common)
                / math.fsum(support[g, end_k]["weight"] for g in end_cohorts),
            "estimate": estimate,
            "variance": max(0.0, variance),
            "standard_error": math.sqrt(max(0.0, variance)),
            "variance_independent_contrast_form": max(0.0, contrast_variance),
            "linear_form_variance_abs_gap": linear_form_gap,
            "within_cohort_normalization_loading_max_abs": within_cohort_weight_sum_gap,
            "standard_error_interpretation": "conditional on these fixed cohort weights; excludes sampling variation in estimated cohort shares",
            "weight_shares": {str(g): normalized[g] for g in common},
        }
        contrasts.append(candidate)
        return candidate

    # Same 27 shared cohorts and weight conventions as the prior point-estimate review.
    pre = make_candidate("common_cohorts_fixed_prebirth_IW", -1, 3, "F1event", "L3event", -1)
    post = make_candidate("common_cohorts_fixed_postbirth_IW", -1, 3, "F1event", "L3event", 3)
    previous = read_json(REPLAY.parents[1] / "empirical_rooms" / "common_cohort_contrast.json")
    previous_candidates = {x["name"]: x for x in previous["variants"]}
    for result in (pre, post):
        old = previous_candidates.get(result["name"])
        result["previous_candidate_estimate"] = old["estimate"] if old else None
        result["gap_from_previous_candidate"] = (result["estimate"] - old["estimate"]) if old else None

    # Here -2 is the omitted reference. It is zero only for cohorts with observed -2 support.
    if any(event == "F2event" for _, event in by_cohort_event):
        raise ValueError("F2event unexpectedly appears as an estimated coefficient")
    minus2_pre = make_candidate("common_cohorts_minus2_to_plus2_fixed_minus2_IW", -2, 2,
                                None, "L2event", -2, omitted_start=True)
    minus2_post = make_candidate("common_cohorts_minus2_to_plus2_fixed_plus2_IW", -2, 2,
                                 None, "L2event", 2, omitted_start=True)
    for result in (minus2_pre, minus2_post):
        result["reference_zero_support"] = {
            str(g): {"observations_at_minus2": support[g, -2]["observations"],
                     "weight_at_minus2": support[g, -2]["weight"]} for g in result["common_cohorts"]}
        result["within_cohort_normalization_loading_max_abs"] = None

    diagnostics = read_csv(paths["room_code_diagnostics.csv"])
    codes = {row["stage"]: row for row in diagnostics}
    summary = {
        "status": "analyzed completed exact replay; candidate contrasts are not adopted",
        "execution_status": execution["status"],
        "execution_elapsed_seconds": execution.get("elapsed_seconds"),
        "processors": execution.get("processors"),
        "replay_script_sha256": execution["script_sha256"],
        "target_gate": {"expected": EXPECTED_TARGET, "observed": number(gate["target_observed"]),
                        "full_precision_abs_gap": abs(number(gate["target_observed"]) - EXPECTED_TARGET),
                        "tolerance": NUMERIC_TOLERANCE},
        "se_gate": {"expected": EXPECTED_SE, "observed": number(gate["se_observed"]),
                    "full_precision_abs_gap": abs(number(gate["se_observed"]) - EXPECTED_SE),
                    "tolerance": NUMERIC_TOLERANCE},
        "sample_gates": {"observations": int(gate["observations_observed"]),
                         "people": int(gate["people_observed"]),
                         "observations_exact_pass": True, "people_exact_pass": True},
        "cohort_interaction_matrix_dimensions": [q, q],
        "cohort_event_coefficient_count": q,
        "full_covariance_cells": cells,
        "covariance_max_symmetry_abs_gap": max_symmetry_gap,
        "covariance_max_diagonal_vs_marginal_variance_abs_gap": max_diagonal_gap,
        "room_code_diagnostics": {
            stage: {"observations": int(row["observations"]),
                    "room_code_0_count": int(row["room_code_0_count"]),
                    "room_code_9_count": int(row["room_code_9_count"]),
                    "room_code_98_count": int(row["room_code_98_count"]),
                    "room_code_99_count": int(row["room_code_99_count"]),
                    "minimum_nonmissing_rooms": number(row["minimum_nonmissing_rooms"])}
            for stage, row in codes.items()},
        "contrasts": contrasts,
        "limitations": [
            "All four contrasts use coefficients from this one regression and impose fixed common-cohort weights.",
            "Their standard errors are conditional on those chosen cohort weights and exclude sampling variation in estimated cohort shares.",
            "The -2 to +2 candidate assigns the omitted F2event coefficient zero only when the same cohort has observed support at -2 and +2.",
            "Cohort composition changes across support intersections and weight endpoints.",
            "These are descriptive alternative aggregations, not new calibration targets and not evidence of causal identification or adoption.",
        ],
    }
    hashes["lead_preflight.json"] = sha256(REPLAY / "lead_preflight.json")
    hashes["source_package_manifest.json"] = sha256(REPLAY / "source_package_manifest.json")
    baseline = REPLAY.parents[1] / "empirical_rooms" / "common_cohort_contrast.json"
    baseline_csv = baseline.with_suffix(".csv")
    hashes["previous_common_cohort_contrast.json"] = sha256(baseline)
    hashes["previous_common_cohort_contrast.csv"] = sha256(baseline_csv)
    summary["input_hashes"] = hashes
    summary["postprocessor_sha256"] = sha256(Path(__file__))
    return summary, weight_rows


def main() -> int:
    try:
        summary, weights = analyze()
    except PendingReplay:
        return 3
    except Exception as exc:
        report = {"status": "analysis_failed", "error": f"{type(exc).__name__}: {exc}"}
        RECEIPT.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        return 2
    if summary.get("status") == "preanalysis_failed":
        RECEIPT.write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2))
        return 2
    with WEIGHTS_OUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(weights[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(weights)
    RECEIPT.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
