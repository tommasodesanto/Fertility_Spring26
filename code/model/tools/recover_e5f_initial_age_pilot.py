#!/usr/bin/env python3
"""Recover report artifacts from a completed native initial-age pilot.

This driver is read-only with respect to ``age_pilot_joint``.  It accepts the
saved numerical equality certified by the frozen controller while recording
the checkpoint-provenance fields that legitimately differ across repetitions.
It never imports or calls a model solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
from typing import Any


PILOT_SOURCE_SHA256 = "51b6c6d6de51edc9387cffc76d6f5514c658aac8711655e369e73829f3b11ce1"
REBATED_HELPER_SHA256 = "d9aa97b890442d45971ec622b4b41687da10ffa26eedaf0198f352c7e6ecb790"
JOINT_ADAPTER_SHA256 = "9eee3bca39f2a98f4a58cf18196d695b6a9db3e93ef18f8eaa2bf4dbb1243bbb"
SCORED_EVALUATOR_SHA256 = "9da35b55466d74dc10a6a85ff91dabe887a807aab68417eabf63a9d7a42ca967"
AGE_HELPER_SHA256 = "5ef6be93fbf28f929d0e09f3eb4ffbe932d61aee5cd66eb707efdea6acd0d801"
CPS_AGE_DATA_SHA256 = "b415b75c916a61113f1ae054ba0baadecd35a26a00c88a21451262937271f210"
EXPECTED = {
    "baseline_original_loss": 182.6491468669,
    "baseline_extra_loss": 123.938117208,
    "selected_original_loss": 180.1870632953,
    "selected_extra_loss": 126.0058524194,
}
EXPECTED_TOLERANCE = 5e-10
EXPECTED_FAILED_CASE = "proposal_03_kappa_fert_minus"


def read(path: Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha(path: Path) -> str:
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def canonical_sha(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write(path: Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path.name}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_module(name: str, path: Path) -> Any:
    specification = importlib.util.spec_from_file_location(name, Path(path).resolve())
    if specification is None or specification.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def recursive_differences(left: Any, right: Any, path: str = "") -> list[dict[str, Any]]:
    """Return every literal JSON difference without discarding provenance fields."""
    if type(left) is not type(right):
        return [{"path": path, "left": left, "right": right, "kind": "type"}]
    if isinstance(left, dict):
        differences = []
        for key in sorted(set(left) | set(right)):
            child = f"{path}/{key}" if path else str(key)
            if key not in left:
                differences.append({"path": child, "left_missing": True, "right": right[key]})
            elif key not in right:
                differences.append({"path": child, "left": left[key], "right_missing": True})
            else:
                differences.extend(recursive_differences(left[key], right[key], child))
        return differences
    if isinstance(left, list):
        differences = []
        for index in range(max(len(left), len(right))):
            child = f"{path}/{index}" if path else str(index)
            if index >= len(left):
                differences.append({"path": child, "left_missing": True, "right": right[index]})
            elif index >= len(right):
                differences.append({"path": child, "left": left[index], "right_missing": True})
            else:
                differences.extend(recursive_differences(left[index], right[index], child))
        return differences
    return [] if left == right else [{"path": path, "left": left, "right": right}]


def allowed_checkpoint_difference(path: str) -> bool:
    parts = path.split("/")
    return (len(parts) == 3 and parts[0] == "target_fit" and parts[1].isdigit()
            and parts[2] == "model_checkpoint_sha256")


def bug_scope(score: dict[str, Any]) -> dict[str, Any]:
    """The exact overly broad scope used by the failed pilot at line 348."""
    return {"loss": score["loss"], "target_fit": score["target_fit"],
            "parameters": score["parameters"]}


def allowed_full_score_difference(path: str) -> bool:
    return (path in {"checkpoint_sha256", "evaluation_receipt_sha256",
                     "normalization/stationary_solve_seconds"}
            or allowed_checkpoint_difference(path))


def verify_extra(extra: dict[str, Any], contract: dict[str, Any]) -> None:
    rows = extra.get("rows", [])
    if len(rows) != 6 or len(extra.get("age_profiles", [])) != 5:
        raise ValueError("saved extra score is not six rows plus five complete age profiles")
    total = sum(float(row["loss_contribution"]) for row in rows)
    if not math.isclose(total, float(extra["extra_loss"]), rel_tol=1e-13, abs_tol=1e-13):
        raise ValueError("saved extra loss does not equal its six contributions")
    metadata = extra.get("metadata", {})
    if (metadata.get("empirical_source_sha256") != contract["cps_age_profile_sha256"]
            or metadata.get("weight_is_empirical_standard_error") is not False
            or metadata.get("overlap_covariance_treatment")
            != "ignored in provisional experimental objective"):
        raise ValueError("saved extra score changed its target or provisional-weight contract")
    for row in rows:
        if row.get("weight_status") != "provisional_diagonal_synthetic_5pct_of_target_not_standard_error":
            raise ValueError("saved extra row does not carry the synthetic-weight disclosure")


def verify_accounting(receipt: dict[str, Any], original: Path) -> None:
    checkpoint = Path(receipt["checkpoint"]).resolve()
    if not checkpoint.is_relative_to(original):
        raise ValueError("checkpoint escapes the immutable original pilot directory")
    if sha(checkpoint) != receipt["checkpoint_sha256"]:
        raise ValueError(f"checkpoint hash changed: {checkpoint}")
    raw_summary = read(checkpoint.parent / "summary.json")
    if raw_summary.get("checkpoint_sha256") != receipt["checkpoint_sha256"]:
        raise ValueError(f"raw summary does not pin its checkpoint: {checkpoint}")
    revenue = float(receipt["property_tax_revenue"])
    outlays = float(receipt["rebate_outlays"])
    residual = float(receipt["rebate_residual"])
    if (not all(math.isfinite(value) for value in (revenue, outlays, residual))
            or revenue <= 0 or outlays <= 0
            or not math.isclose(revenue - outlays, residual, rel_tol=1e-12, abs_tol=1e-12)
            or float(receipt["rebate_relative_gap"]) > 1e-6
            or float(receipt["pension_relative_gap"]) > 1e-6):
        raise ValueError(f"independent fiscal accounting gate failed: {checkpoint}")


def verify_case(case_dir: Path, case_id: str, controller: Any, restrictions: dict[str, Any],
                contract: dict[str, Any], original: Path) -> dict[str, Any]:
    proposal = read(case_dir / "proposal.json")
    evaluation = case_dir / "evaluation"
    candidate = read(evaluation / "candidate_result.json")
    summary = read(evaluation / "summary.json")
    launch = read(evaluation / "launch_contract.json")
    score1 = read(evaluation / "case/evaluation/scored_repetition_01/score.json")
    accounting = read(evaluation / "case/evaluation/rebate_accounting.json")
    if (proposal.get("case_id") != case_id or launch.get("proposal") != proposal
            or candidate.get("proposal") != proposal or candidate.get("status") != "verified"):
        raise ValueError(f"proposal or candidate receipt mismatch: {case_id}")
    if (summary.get("status") != "verified_rebated_initial_smoke"
            or summary.get("method") != "joint_price_fertility_rebate"
            or summary.get("loss") != score1.get("loss")
            or candidate.get("loss") != score1.get("loss")
            or candidate.get("score") != score1):
        raise ValueError(f"score/summary mismatch: {case_id}")
    if (launch.get("helper_sha256") != contract["rebated_helper_sha256"]
            or launch.get("joint_sha256") != contract["joint_adapter_sha256"]
            or launch.get("objective_sha256") != contract["original_objective_canonical_sha256"]):
        raise ValueError(f"source or target pin mismatch: {case_id}")
    if Path(candidate.get("output", "")).resolve() != (evaluation / "case/evaluation").resolve():
        raise ValueError(f"candidate output path mismatch: {case_id}")
    if len(score1.get("target_fit", [])) != 13 or sum(row.get("scored") is True
                                                       for row in score1["target_fit"]) != 12:
        raise ValueError(f"incomplete original target system: {case_id}")
    controller.validate_score(score1, proposal, restrictions)
    receipts = accounting.get("receipts", [])
    if accounting.get("status") != "verified" or receipts != candidate.get("accounting"):
        raise ValueError(f"independent accounting receipt mismatch: {case_id}")
    if len(receipts) != int(proposal.get("repetitions", 1)):
        raise ValueError(f"wrong accounting repetition count: {case_id}")
    for receipt in receipts:
        verify_accounting(receipt, original)
    scores = [score1]
    extra1 = read(case_dir / "extra_score.json")
    verify_extra(extra1, contract)
    score2 = extra2 = None
    if int(proposal.get("repetitions", 1)) == 2:
        score2 = read(evaluation / "case/evaluation/scored_repetition_02/score.json")
        extra2 = read(case_dir / "extra_score_repetition_02.json")
        controller.validate_score(score2, proposal, restrictions)
        if (candidate.get("second_signature_equal") is not True
                or controller.numeric_signature(score1) != controller.numeric_signature(score2)
                or extra1 != extra2):
            raise ValueError("two saved exact repetitions are not numerically identical")
        scores.append(score2)
    for index, (score, receipt) in enumerate(zip(scores, receipts), 1):
        scored_dir = evaluation / f"case/evaluation/scored_repetition_{index:02d}"
        evaluation_receipt = read(scored_dir / "verified_evaluation_receipt.json")
        if score.get("checkpoint_sha256") != receipt["checkpoint_sha256"]:
            raise ValueError(f"top score checkpoint pin disagrees with accounting: {case_id} repetition {index}")
        if (evaluation_receipt.get("checkpoint_sha256") != receipt["checkpoint_sha256"]
                or canonical_sha(evaluation_receipt) != score.get("evaluation_receipt_sha256")):
            raise ValueError(f"score evaluation-receipt pin is invalid: {case_id} repetition {index}")
    return {
        "case_id": case_id, "proposal": proposal, "score": score1, "score2": score2,
        "extra": extra1, "extra2": extra2, "accounting": receipts,
        "original_loss": float(score1["loss"]), "extra_loss": float(extra1["extra_loss"]),
        "augmented_loss": float(score1["loss"]) + float(extra1["extra_loss"]),
        "output": str(evaluation.resolve()),
    }


def verify_failed_case(case_dir: Path, case_id: str, contract: dict[str, Any]) -> dict[str, Any]:
    """Preserve the one known numerical failure without treating it as a score."""
    if case_id != EXPECTED_FAILED_CASE:
        raise ValueError(f"unexpected failed primary case: {case_id}")
    proposal = read(case_dir / "proposal.json")
    receipt = read(case_dir / "case_receipt.json")
    evaluation = case_dir / "evaluation"
    candidate = read(evaluation / "candidate_result.json")
    summary = read(evaluation / "summary.json")
    launch = read(evaluation / "launch_contract.json")
    if (proposal.get("case_id") != case_id or receipt.get("case_id") != case_id
            or receipt.get("status") != "failed_candidate"
            or receipt.get("proposal") != proposal or candidate.get("proposal") != proposal
            or candidate.get("status") != "failed"
            or summary.get("status") != "failed_joint_rebated_initial"):
        raise ValueError("known numerical failure lacks matching native failure receipts")
    if (launch.get("helper_sha256") != contract["rebated_helper_sha256"]
            or launch.get("joint_sha256") != contract["joint_adapter_sha256"]
            or launch.get("objective_sha256") != contract["original_objective_canonical_sha256"]):
        raise ValueError("known numerical failure used a different source or target contract")
    if Path(candidate.get("output", "")).resolve() != (evaluation / "case/evaluation").resolve():
        raise ValueError("known numerical failure output path mismatch")
    detail = str(candidate.get("failure_detail", ""))
    if not detail or receipt.get("failure_detail") != detail:
        raise ValueError("known numerical failure detail is absent or inconsistent")
    branch_path = evaluation / "case/branch_failure.json"
    if not branch_path.is_file():
        raise ValueError("known numerical failure has no native branch_failure receipt")
    branch = read(branch_path)
    if branch.get("status") != "hard_error" or branch.get("error") != detail:
        raise ValueError("known numerical failure branch receipt disagrees with candidate result")
    return {
        "case_id": case_id, "status": "failed_candidate", "proposal": proposal,
        "failure_detail": detail, "branch_failure": branch, "output": str(evaluation.resolve()),
    }


def close_to_expected(actual: float, expected: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=EXPECTED_TOLERANCE):
        raise ValueError(f"{label} differs from the recorded failed-run receipt")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--original", default="age_pilot_joint")
    parser.add_argument("--output", default="age_pilot_recovered")
    args = parser.parse_args()
    batch = args.batch.resolve()
    original = (batch / args.original).resolve()
    output = (batch / args.output).resolve()
    if not original.is_relative_to(batch) or not output.is_relative_to(batch) or original == output:
        raise ValueError("original and recovered outputs must be distinct children of --batch")
    if output.exists():
        raise FileExistsError(f"recovery output already exists: {output}")

    contract = read(original / "contract.json")
    source_paths = {
        "pilot_source_sha256": batch / "age_source_v2/run_e5f_initial_age_pilot.py",
        "rebated_helper_sha256": batch / "initial_v2/run_e5f_rebated_initial_overnight.py",
        "joint_adapter_sha256": batch / "joint_source/run_e5f_joint_rebated_initial_probe.py",
        "scored_evaluator_sha256": batch / "joint_source/run_e5f_joint_rebated_initial_scored.py",
        "extra_age_helper_sha256": batch / "age_source_v2/e5f_initial_age_profile.py",
        "cps_age_profile_sha256": batch / "age_profile_candidates.json",
    }
    fixed = {"pilot_source_sha256": PILOT_SOURCE_SHA256,
             "rebated_helper_sha256": REBATED_HELPER_SHA256,
             "joint_adapter_sha256": JOINT_ADAPTER_SHA256,
             "scored_evaluator_sha256": SCORED_EVALUATOR_SHA256,
             "extra_age_helper_sha256": AGE_HELPER_SHA256,
             "cps_age_profile_sha256": CPS_AGE_DATA_SHA256}
    verified_sources = {}
    for key, path in source_paths.items():
        actual = sha(path)
        if actual != fixed[key] or (key != "pilot_source_sha256" and contract.get(key) != actual):
            raise ValueError(f"frozen source/target pin changed: {path}")
        verified_sources[str(path)] = actual

    smoke_path = Path(contract["smoke"]).resolve()
    if smoke_path != (batch / "joint_initial_smoke/summary.json").resolve():
        raise ValueError("pilot contract points to a different initial smoke")
    if sha(smoke_path) != contract["verified_joint_rebated_smoke_summary_sha256"]:
        raise ValueError("joint initial smoke summary changed")
    smoke = read(smoke_path)
    if smoke.get("status") != "verified_rebated_initial_smoke" or smoke.get("method") != "joint_price_fertility_rebate":
        raise ValueError("joint initial smoke is not verified")

    helper = load_module("recovery_rebated_helper", source_paths["rebated_helper_sha256"])
    packet = helper.saved_packet(Path(contract["template"]))
    restrictions = helper.validate_scientific_contract(packet)
    if (digest := sha(packet["objective_path"])) != contract["original_objective_file_sha256"]:
        raise ValueError("original objective file changed")
    if packet["run"]["working_objective"]["canonical_sha256"] != contract["original_objective_canonical_sha256"]:
        raise ValueError("original objective canonical hash changed")
    controller_path = Path(packet["plan_path"]).parent / "run_capped_beta.py"
    controller = load_module("recovery_saved_controller", controller_path)
    verified_sources[str(packet["objective_path"])] = digest
    verified_sources[str(controller_path)] = sha(controller_path)

    component_keys = ("original_objective_canonical_sha256", "original_objective_file_sha256",
                      "cps_age_profile_sha256", "extra_age_helper_sha256", "rebated_helper_sha256",
                      "joint_adapter_sha256", "scored_evaluator_sha256",
                      "verified_joint_rebated_smoke_summary_sha256", "extra_weight_rule")
    components = {key: contract[key] for key in component_keys}
    fingerprint = hashlib.sha256(json.dumps(components, sort_keys=True,
                                             separators=(",", ":")).encode()).hexdigest()
    if fingerprint != contract.get("fingerprint"):
        raise ValueError("age-pilot aggregate fingerprint changed")

    primary_ids = [item["case_id"] for item in contract.get("primary_cases", [])]
    if len(primary_ids) != 6 or len(set(primary_ids)) != 6 or "baseline" not in primary_ids:
        raise ValueError("frozen pilot does not contain baseline plus five unique proposals")
    cases = []
    for case_id in primary_ids:
        case_dir = original / "cases" / case_id
        native_receipt = read(case_dir / "case_receipt.json")
        if native_receipt.get("status") == "verified":
            case = verify_case(case_dir, case_id, controller, restrictions, contract, original)
            case["status"] = "verified"
        elif native_receipt.get("status") == "failed_candidate":
            case = verify_failed_case(case_dir, case_id, contract)
        else:
            raise ValueError(f"unknown native primary-case status: {case_id}")
        cases.append(case)
    verified_cases = [case for case in cases if case["status"] == "verified"]
    failed_cases = [case for case in cases if case["status"] == "failed_candidate"]
    if len(verified_cases) != 5 or [case["case_id"] for case in failed_cases] != [EXPECTED_FAILED_CASE]:
        raise ValueError("recovery requires five verified primary cases and the one known numerical failure")
    exact = verify_case(original / "cases/selected_exact_repetitions", "selected_exact_repetitions",
                        controller, restrictions, contract, original)
    exact["status"] = "verified"
    baseline = next(case for case in verified_cases if case["case_id"] == "baseline")
    selected = min(verified_cases, key=lambda case: case["augmented_loss"])
    matching = [case for case in verified_cases
                if case["proposal"]["parameters"] == exact["proposal"]["parameters"]]
    if matching != [selected] or selected["case_id"] == "baseline":
        raise ValueError("exact repetitions do not correspond to the unique improved selected proposal")
    if (controller.numeric_signature(selected["score"]) != controller.numeric_signature(exact["score"])
            or selected["extra"] != exact["extra"]
            or not selected["augmented_loss"] < baseline["augmented_loss"]):
        raise ValueError("selected and exact saved outputs are not numerically identical or improved")

    close_to_expected(baseline["original_loss"], EXPECTED["baseline_original_loss"], "baseline original loss")
    close_to_expected(baseline["extra_loss"], EXPECTED["baseline_extra_loss"], "baseline extra loss")
    close_to_expected(selected["original_loss"], EXPECTED["selected_original_loss"], "selected original loss")
    close_to_expected(selected["extra_loss"], EXPECTED["selected_extra_loss"], "selected extra loss")

    selected_to_exact_scope = recursive_differences(bug_scope(selected["score"]), bug_scope(exact["score"]))
    exact_internal_scope = recursive_differences(bug_scope(exact["score"]), bug_scope(exact["score2"]))
    scope_differences = [{"comparison": "selected_vs_exact_repetition_01", **row}
                         for row in selected_to_exact_scope]
    scope_differences.extend({"comparison": "exact_repetition_01_vs_02", **row}
                             for row in exact_internal_scope)
    if (not selected_to_exact_scope
            or any(not allowed_checkpoint_difference(row["path"]) for row in scope_differences)):
        raise ValueError("failed-pilot comparison scope differs beyond target-row checkpoint provenance")
    selected_to_exact_full = recursive_differences(selected["score"], exact["score"])
    exact_internal_full = recursive_differences(exact["score"], exact["score2"])
    full_differences = [{"comparison": "selected_vs_exact_repetition_01", **row}
                        for row in selected_to_exact_full]
    full_differences.extend({"comparison": "exact_repetition_01_vs_02", **row}
                            for row in exact_internal_full)
    if any(not allowed_full_score_difference(row["path"]) for row in full_differences):
        raise ValueError("full saved scores contain an unexpected nonnumerical discrepancy")

    # Validation is complete. Only now create the separate recovered packet.
    output.mkdir(parents=True, exist_ok=False)
    selected_score, selected_extra = exact["score"], exact["extra"]
    original_rows = selected_score["target_fit"]
    scored_rows = [row for row in original_rows if row.get("scored") is True]
    augmented_rows = [dict(system="original_12", **row) for row in scored_rows]
    augmented_rows.extend(dict(system="experimental_extra_6", **row)
                          for row in selected_extra["rows"])
    if len(original_rows) != 13 or len(scored_rows) != 12 or len(augmented_rows) != 18:
        raise RuntimeError("recovered table dimensions are incomplete")
    write_csv(output / "selected_original_target_fit.csv", original_rows)
    write_csv(output / "selected_extra_target_fit.csv", selected_extra["rows"])
    write_csv(output / "selected_augmented_target_fit.csv", augmented_rows)
    write_csv(output / "selected_candidate_parameters.csv", selected_score["parameters"])
    write(output / "selected_extra_age_profiles.json", selected_extra["age_profiles"])
    write(output / "selected_checkpoint.json", exact["accounting"][-1])

    all_cases = [*cases, exact]
    all_original, all_extra, all_parameters, all_profiles = [], [], [], []
    for case in all_cases:
        if case["status"] != "verified":
            continue
        all_original.extend(dict(case_id=case["case_id"], **row) for row in case["score"]["target_fit"])
        all_extra.extend(dict(case_id=case["case_id"], **row) for row in case["extra"]["rows"])
        all_parameters.extend(dict(case_id=case["case_id"], **row) for row in case["score"]["parameters"])
        all_profiles.extend(dict(case_id=case["case_id"], **row) for row in case["extra"]["age_profiles"])
    write_csv(output / "all_original_target_fits.csv", all_original)
    write_csv(output / "all_extra_target_fits.csv", all_extra)
    write_csv(output / "all_candidate_parameters.csv", all_parameters)
    write(output / "all_age_profiles.json", all_profiles)
    compact = [{key: value for key, value in case.items()
                if key not in {"score", "score2", "extra", "extra2"}} for case in all_cases]
    write(output / "all_cases.json", compact)
    write(output / "cases.json", compact)
    write(output / "exact_score_discrepancies.json", {
        "status": "verified_nonnumerical_metadata_differences_only",
        "numerical_signature_equal": True,
        "extra_score_equal": True,
        "failed_pilot_scope": ["loss", "target_fit", "parameters"],
        "failed_pilot_scope_differences": scope_differences,
        "full_score_differences": full_differences,
        "top_level_checkpoint_and_evaluation_receipt_pins_verified_independently": True,
        "stationary_solve_seconds_is_recorded_but_excluded_from_numerical_signature": True,
    })
    write(output / "source_and_checkpoint_verification.json", {
        "status": "verified", "original_read_only": str(original),
        "pilot_fingerprint": fingerprint, "verified_sources": verified_sources,
        "verified_checkpoint_count": sum(len(case["accounting"]) for case in all_cases
                                         if case["status"] == "verified"),
    })
    summary = {
        "status": "recovered_verified_initial_age_pilot",
        "recovery_method": "saved_native_outputs_only_no_model_solves",
        "original_output": str(original), "selected_case": selected["case_id"],
        "baseline_original_loss": baseline["original_loss"], "baseline_extra_loss": baseline["extra_loss"],
        "baseline_augmented_loss": baseline["augmented_loss"],
        "selected_original_loss": selected["original_loss"], "selected_extra_loss": selected["extra_loss"],
        "selected_augmented_loss": selected["augmented_loss"], "selected_improves_baseline": True,
        "selected_exact_repetitions_verified": True, "selected_checkpoint": exact["accounting"][-1],
        "primary_cases": 6, "verified_primary_cases": 5, "failed_primary_cases": 1,
        "failed_case": EXPECTED_FAILED_CASE, "verified_cases_including_exact": 6,
        "original_table_rows_including_normalization": 13, "original_scored_rows": 12,
        "experimental_extra_rows": 6, "augmented_scored_rows": 18,
        "failed_pilot_scope_difference_count": len(scope_differences),
        "full_score_difference_count": len(full_differences),
        "literal_score_differences": "exact_score_discrepancies.json",
        "production_eligible": False, "pilot_fingerprint": fingerprint,
    }
    write(output / "summary.json", summary)
    print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
