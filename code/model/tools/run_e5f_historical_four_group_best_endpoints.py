#!/usr/bin/env python3
"""Audit closed and open stationary endpoints for the selected four-group fit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import run_e5f_historical_demographic_closure_candidate as candidate_driver
import run_e5f_open_population_transition as transition
import run_e5f_post2023_no_policy_continuations as continuation


ROOT = Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coordinate-report", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--selected-report", type=Path, default=candidate_driver.DEFAULT_REPORT)
    parser.add_argument("--source", type=Path, default=candidate_driver.DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    coordinate_report = json.loads(args.coordinate_report.resolve().read_text(encoding="utf-8"))
    if coordinate_report.get("status") != "complete_historical_four_group_coordinate_collection":
        raise RuntimeError("Coordinate report is incomplete")
    task_id = int(coordinate_report["best"]["task_id"])
    candidate_path = args.plan_dir.resolve() / f"candidate_{task_id:03d}.json"
    selected_report_path = args.selected_report.resolve()
    source_path = args.source.resolve()
    candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
    selected_report = json.loads(selected_report_path.read_text(encoding="utf-8"))
    theta, delta = candidate_driver.validate_candidate(candidate, selected_report)

    report_sha = continuation.file_sha256(selected_report_path)
    if str(candidate.get("selected_report_sha256")) != report_sha:
        raise RuntimeError("Selected-report hash differs from the coordinate candidate")
    if (
        str(coordinate_report["common_contract"]["selected_report_sha256"])
        != report_sha
    ):
        raise RuntimeError("Coordinate report and selected report are not the same run")
    source_sha = continuation.file_sha256(source_path)
    if str(coordinate_report["common_contract"]["source_sha256"]) != source_sha:
        raise RuntimeError("Coordinate report and stationary source are not the same run")

    chain, model = transition.configure_sequential_model()
    live_bundle = candidate_driver.calibration.code_fingerprint_contract(model)[
        "bundle_sha256"
    ]
    if not (
        live_bundle
        == str(candidate.get("scientific_bundle_sha256"))
        == str(coordinate_report["common_contract"]["scientific_bundle_sha256"])
    ):
        raise RuntimeError("Live scientific bundle differs from the selected coordinate")
    target_fingerprint = candidate_driver.calibration.e5_target_system().fingerprint
    if not (
        target_fingerprint
        == str(candidate.get("target_fingerprint"))
        == str(coordinate_report["common_contract"]["target_fingerprint"])
    ):
        raise RuntimeError("Live target contract differs from the selected coordinate")
    prepared, renewal_state = candidate_driver.prepare_candidate(
        chain=chain,
        model=model,
        theta=theta,
        initial_old_psi=float(selected_report["best_candidate"]["old_psi_child"]),
        outside_origin_share=float(
            selected_report["scientific_contract"]["outside_origin_entry_share"]
        ),
    )
    old_psi = float(renewal_state["old_fertility_normalization"]["psi_child"])
    new_psi = old_psi + float(delta)
    closed, schedule = continuation.solve_closed_stationary_endpoint(
        chain=prepared.chain,
        base_overrides=prepared.base_overrides,
        old_price=prepared.old_price,
        new_psi_child=new_psi,
        supply_rule=prepared.supply_rule,
        price_min_ratio=1e-4,
        price_max_ratio=3.0,
        grid_points=25,
        outdir=outdir,
    )
    continuation.write_json(outdir / "closed_stationary_endpoint.json", closed)
    continuation.write_csv(outdir / "closed_stationary_schedule.csv", schedule)

    open_dir = outdir / "open_endpoint"
    open_dir.mkdir(parents=True, exist_ok=True)
    open_raw = transition.solve_stationary_open_endpoint(
        chain=prepared.chain,
        base_overrides=prepared.base_overrides,
        old_price=prepared.old_price,
        new_psi_child=new_psi,
        outside_flow=float(renewal_state["outside_flow_M"]),
        retention=float(renewal_state["retention_rho"]),
        effective_birth_to_household_conversion=(
            transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
        ),
        supply_rule=prepared.supply_rule,
        policy_case="none",
        outdir=open_dir,
    )
    open_endpoint = continuation.verify_open_endpoint(
        open_raw,
        float(renewal_state["outside_flow_M"]),
        float(renewal_state["retention_rho"]),
    )
    continuation.write_json(outdir / "open_stationary_endpoint.json", open_endpoint)
    summary = {
        "status": "complete_historical_four_group_endpoint_audit",
        "selected_task_id": task_id,
        "selected_candidate_id": candidate["candidate_id"],
        "selected_transition_loss": coordinate_report["best"]["loss"],
        "selected_report_sha256": report_sha,
        "source_sha256": source_sha,
        "scientific_bundle_sha256": live_bundle,
        "target_fingerprint": target_fingerprint,
        "old_psi_child": old_psi,
        "new_psi_child": new_psi,
        "closed_stationary_endpoint": closed,
        "open_stationary_endpoint": open_endpoint,
        "between_steady_states_label_allowed": bool(
            closed.get("usable_closed_root", False)
        ),
    }
    continuation.write_json(outdir / "summary.json", summary)
    print(
        "HISTORICAL_FOUR_GROUP_ENDPOINT_AUDIT_COMPLETE "
        f"candidate={candidate['candidate_id']} "
        f"closed_root={summary['between_steady_states_label_allowed']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
