"""Recover a stopped rebated-initial search by exactly verifying its saved best."""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import math
from pathlib import Path
import time


HELPER_SHA256 = "d9aa97b890442d45971ec622b4b41687da10ffa26eedaf0198f352c7e6ecb790"
JOINT_SHA256 = "9eee3bca39f2a98f4a58cf18196d695b6a9db3e93ef18f8eaa2bf4dbb1243bbb"
SCORED_SHA256 = "9da35b55466d74dc10a6a85ff91dabe887a807aab68417eabf63a9d7a42ca967"
SEARCH_SHA256 = "be2c3b2cea5723d4fa7e77b36d13a1f07e54b6561edef95c4c641d035dbaa555"
EXPECTED_SOURCE_ROOT = Path(
    "/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8"
)
OUTPUT_NAME = "initial_search_recovered"
ALLOWED_CASE_STATUSES = frozenset({
    "verified", "rejected_numerical", "rejected_equilibrium",
    "rejected_mass_gate", "timeout",
})


def load_frozen(name, path):
    spec = importlib.util.spec_from_file_location(name, Path(path).resolve())
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_original_stop(summary, records):
    """Admit only the known numerical-stop shape; source/contract failures remain fatal."""
    if summary.get("status") != "stopped_branch_hard_error":
        raise ValueError("Original search is not the stopped branch eligible for recovery")
    if not records or not any(row.get("status") == "verified" for row in records):
        raise ValueError("Original search has no verified candidate")
    bad = sorted({str(row.get("status")) for row in records
                  if row.get("status") not in ALLOWED_CASE_STATUSES})
    if bad:
        raise ValueError("Original cases contain a hard source/contract failure: " + ", ".join(bad))
    rejected = [row for row in records if str(row.get("status", "")).startswith("rejected_")]
    if not rejected:
        raise ValueError("Original hard stop has no recorded numerical/equilibrium rejection")


def validate_receipts(helper, result):
    receipts = result.get("accounting", [])
    if (result.get("receipt_status") != "verified_scored_candidate"
            or result.get("second_signature_equal") is not True
            or result.get("proposal", {}).get("repetitions") != 2
            or len(receipts) != 2):
        raise ValueError("Two exact repetition receipts are incomplete")
    if len({row.get("repetition") for row in receipts}) != 2:
        raise ValueError("Exact repetition identities are not distinct")
    for row in receipts:
        checkpoint = Path(row["checkpoint"])
        if not checkpoint.is_file() or helper.sha(checkpoint) != row["checkpoint_sha256"]:
            raise ValueError("Exact repetition checkpoint fingerprint changed")
        if (float(row["rebate_relative_gap"]) > helper.REBATE_RELATIVE_TOLERANCE
                or float(row["pension_relative_gap"]) > helper.PENSION_RELATIVE_TOLERANCE):
            raise ValueError("Exact repetition fiscal accounting gate failed")


def compact(result):
    return {key: value for key, value in result.items() if key != "score"}


def run(batch, seconds, original_search):
    if not 1 <= seconds <= 2100:
        raise ValueError("Recovery budget must be between 1 and 2100 seconds")
    batch = Path(batch).resolve()
    if Path(original_search).name != original_search:
        raise ValueError("Original search must be a direct batch child name")
    original = batch / original_search
    out = batch / OUTPUT_NAME
    if out.exists():
        raise FileExistsError(f"Recovery output already exists: {out}")

    helper_path = batch / "initial_v2/run_e5f_rebated_initial_overnight.py"
    joint_path = batch / "joint_source/run_e5f_joint_rebated_initial_probe.py"
    scored_path = batch / "joint_source/run_e5f_joint_rebated_initial_scored.py"
    search_path = batch / "joint_source/run_e5f_rebated_initial_search.py"
    search = load_frozen("frozen_rebated_initial_search", search_path)
    if search.file_sha(search_path) != SEARCH_SHA256:
        raise ValueError("Frozen rebated search source pin changed")
    search.require_source_pins(helper_path, HELPER_SHA256, joint_path, JOINT_SHA256,
                               scored_path, SCORED_SHA256)
    helper = search.load_helper(helper_path, HELPER_SHA256)

    original_contract = search.read(original / "search_contract.json")
    original_summary = search.read(original / "summary.json")
    records = search.read(original / "cases.json")
    best = search.read(original / "best_so_far.json")
    validate_original_stop(original_summary, records)
    expected_pins = {
        "helper_sha256": HELPER_SHA256,
        "joint_sha256": JOINT_SHA256,
        "scored_invoker_sha256": SCORED_SHA256,
    }
    for key, expected in expected_pins.items():
        if original_contract.get(key) != expected:
            raise ValueError(f"Original search {key} differs from frozen recovery source")
    source_root = Path(original_contract["source_root"]).resolve()
    if source_root != EXPECTED_SOURCE_ROOT:
        raise ValueError("Original search source root is not the approved saved source")
    template = source_root / "batches/capped_beta_099_20260911"
    packet = helper.saved_packet(template)
    restrictions = helper.validate_scientific_contract(packet)
    controller = search.load_controller(packet)

    if best.get("status") != "verified" or not math.isfinite(float(best.get("loss", math.nan))):
        raise ValueError("Original best-so-far is not a verified finite result")
    matching = [row for row in records if row.get("case_id") == best.get("case_id")
                and row.get("status") == "verified"]
    if len(matching) != 1 or matching[0].get("proposal") != best.get("proposal"):
        raise ValueError("Original best is not uniquely present in the case ledger")
    case_id = str(best["case_id"])
    case_folder = (original / "cases" / case_id).resolve()
    if case_folder.parent != (original / "cases").resolve():
        raise ValueError("Original best case path escapes the search folder")
    original_result = search.read(case_folder / "result.json")
    score_path = Path(best["output"]) / "scored_repetition_01/score.json"
    original_score = search.read(score_path)
    controller.validate_score(original_score, best["proposal"], restrictions)
    if (original_result.get("status") != "verified"
            or original_result.get("receipt_status") != "verified_scored_candidate"
            or float(original_score["loss"]) != float(best["loss"])
            or controller.numeric_signature(original_result["score"])
               != controller.numeric_signature(original_score)):
        raise ValueError("Original best receipt and full saved score differ")

    selected = {
        "case_id": "selected_exact_repetitions",
        "parameters": {name: controller.parameters(original_score)[name]
                       for name in controller.ALL_NAMES},
        "initial_psi": best["proposal"]["initial_psi"],
        "repetitions": 2,
    }
    controller.validate_score(original_score, selected, restrictions)
    out.mkdir(parents=True, exist_ok=False)
    recovered_contract = copy.deepcopy(original_contract)
    recovered_contract.update(
        status="recovering_exact_repetitions",
        final_exact_repetitions=2,
        recovery_seconds=seconds,
        recovery_metadata={
            "original_search": str(original),
            "original_status": original_summary["status"],
            "original_best_case_id": best["case_id"],
            "original_best_loss": best["loss"],
            "original_summary_sha256": helper.sha(original / "summary.json"),
            "original_search_contract_sha256": helper.sha(original / "search_contract.json"),
            "search_driver_sha256": SEARCH_SHA256,
            "recovery_driver": str(Path(__file__).resolve()),
            "source_note": "Nine free parameters; beta <= .99; 13 target rows (12 scored plus the 2.1 normalization); 17 raw parameter rows retained.",
        },
    )
    search.write(out / "search_contract.json", recovered_contract)
    search.write(out / "best_so_far.json", compact(best))
    search.write(out / "latest_completed.json", {
        "status": "running_exact_recovery", "case_id": selected["case_id"],
        "original_best_case_id": best["case_id"], "original_best_loss": best["loss"],
    })

    started = time.monotonic()
    exact = search.run_child(
        helper_path, HELPER_SHA256, joint_path, JOINT_SHA256,
        scored_path, SCORED_SHA256, template, selected,
        out / "cases/selected_exact_repetitions", started + seconds,
    )
    search.write(out / "latest_completed.json", compact(exact))
    search.write(out / "best_so_far.json", compact(exact))
    if exact.get("status") != "verified":
        raise RuntimeError("Recovered exact repetitions did not verify")
    validate_receipts(helper, exact)
    first_path = Path(exact["output"]) / "scored_repetition_01/score.json"
    second_path = Path(exact["output"]) / "scored_repetition_02/score.json"
    first, second = search.read(first_path), search.read(second_path)
    controller.validate_score(first, selected, restrictions)
    controller.validate_score(second, selected, restrictions)
    signature = controller.numeric_signature(original_score)
    if (controller.numeric_signature(first) != signature
            or controller.numeric_signature(second) != signature
            or controller.numeric_signature(exact["score"]) != signature):
        raise ValueError("Recovered full numeric signatures differ from the original best")
    launch = search.read(out / "cases/selected_exact_repetitions/case/launch_contract.json")
    if (launch.get("helper_sha256") != HELPER_SHA256
            or launch.get("joint_sha256") != JOINT_SHA256
            or launch.get("objective_sha256") != packet["run"]["working_objective"]["canonical_sha256"]):
        raise ValueError("Recovered exact run source or target gate differs")
    if len(first["target_fit"]) != 13 or len(first["parameters"]) != 17:
        raise ValueError("Recovered score does not preserve the full 13-row target and 17-row parameter tables")

    helper.write_selected_tables(out, exact)
    search.write(out / "selected_candidate_result.json", exact)
    search.write(out / "selected_recovered_score.json", {
        "loss": first["loss"], "target_fit": first["target_fit"],
        "parameters": first["parameters"], "contract_sha256": first["contract_sha256"],
        "free_parameter_count": first["free_parameter_count"],
        "beta_upper": 0.99, "normalization": 2.1,
        "table_note": "All 13 target rows and all 17 raw parameter rows are preserved from the exact score.",
    })
    recovered_contract["status"] = "recovered_rebated_initial_search"
    search.write(out / "search_contract.json", recovered_contract)
    summary = {
        "status": "recovered_rebated_initial_search",
        "elapsed_seconds": time.monotonic() - started,
        "selected_loss": exact["loss"],
        "selected_exact_repetitions_verified": True,
        "checkpoint": exact["accounting"][-1],
        "source_root": str(source_root),
        "smoke": original_contract["smoke"],
        "final_exact_repetitions": 2,
        "selected_exact_result_path": str(out / "cases/selected_exact_repetitions/result.json"),
        "selected_exact_output": exact["output"],
        "search_contract": str(out / "search_contract.json"),
        "original_search": str(original),
        "original_stop_status": original_summary["status"],
        "original_stop_summary": original_summary,
    }
    search.write(out / "summary.json", summary)
    print(json.dumps(summary, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--seconds", type=int, default=2100)
    parser.add_argument("--original-search", default="initial_search_joint")
    args = parser.parse_args()
    run(args.batch, args.seconds, args.original_search)


if __name__ == "__main__":
    main()
