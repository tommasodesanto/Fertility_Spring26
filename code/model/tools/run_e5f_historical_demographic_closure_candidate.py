#!/usr/bin/env python3
"""Evaluate one calibrated-parameter candidate under the four-group closure."""

from __future__ import annotations

import argparse
import copy
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np

import audit_closed_reproductive_closure as closure
import build_e5f_transition_design_feasibility as design
import run_dynamic_population_transition as calendar
import run_e5f_historical_demographic_closure_smoke as smoke
import run_e5f_open_population_transition as transition
import run_e5f_post2023_no_policy_continuations as continuation
import run_e5f_transition_calibration as calibration


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REPORT = smoke.DEFAULT_REPORT
DEFAULT_SOURCE = smoke.DEFAULT_SOURCE
DEFAULT_AGE_PATH = smoke.DEFAULT_AGE_PATH


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-json", type=Path, required=True)
    parser.add_argument("--selected-report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--age-path", type=Path, default=DEFAULT_AGE_PATH)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--market-tol", type=float, default=2e-4)
    parser.add_argument("--market-max-iter", type=int, default=30)
    return parser.parse_args()


def prepare_candidate(
    *,
    chain: Any,
    model: Any,
    theta: dict[str, float],
    initial_old_psi: float,
    outside_origin_share: float,
) -> tuple[continuation.PreparedModel, dict[str, Any]]:
    _, profile_overrides, _ = calibration.activate_model_profile(
        calibration.REPAIRED_MODEL_PROFILE, theta
    )
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
    (
        old_solution,
        old_parameters,
        old_price_array,
        _,
        normalization,
    ) = calibration.solve_old_steady_state(
        chain,
        base,
        initial_psi=float(initial_old_psi),
        completed_fertility_target=transition.REPLACEMENT_FERTILITY,
        completed_fertility_tolerance=5e-4,
        normalize=True,
    )
    old_price = float(np.asarray(old_price_array).reshape(-1)[0])
    if int(old_parameters.J) != 17 or int(old_parameters.Nb) != 120:
        raise RuntimeError("Candidate production dimensions are invalid")
    b_grid = np.asarray(old_solution.b_grid, dtype=float)
    shared = model.precompute_shared(old_parameters, b_grid)
    old_parameters._fert2_probs = np.asarray(
        old_solution.fert2_probs, dtype=float
    ).copy()
    old_policy = calendar.policy_from_solution(
        old_solution, old_price_array, old_parameters, b_grid, shared
    )
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = (
        transition.advance_sequential_calendar_distribution
    )
    calendar.distribution_rows = transition.independent_child_distribution_rows
    stationary_g_pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(
        old_solution, old_policy, old_parameters, b_grid, shared
    )
    operator_gates = transition.operator_gates(
        old_solution,
        old_policy,
        stationary_g_pre,
        old_parameters,
        b_grid,
        shared,
    )
    operator_gates.update(reconstruction)
    gate_names = (
        "stationary_post_fertility_nesting_l1",
        "one_step_constant_path_nesting_l1",
        "mature_flow_abs_error",
        "birth_flow_abs_error",
        "topcode_adjusted_birth_flow_abs_error",
    )
    if max(abs(float(operator_gates[name])) for name in gate_names) > 5e-9:
        raise RuntimeError(f"Candidate stationary-operator gates failed: {operator_gates}")

    renewal = closure.topcode_consistent_renewal_accounting(
        old_solution, old_parameters
    )
    births_old = float(renewal["topcode_adjusted_birth_children"])
    births_old_raw = float(renewal["raw_birth_children"])
    E_old = float(old_solution.entry_rate)
    B_old = transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION * births_old
    renewal_state = transition.stationary_renewal_from_births(
        E_old,
        births_old,
        float(outside_origin_share),
        transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION,
    )
    if abs(float(renewal_state["identity_residual"])) > 2e-12:
        raise RuntimeError("Candidate old renewal identity failed")
    if abs(float(renewal_state["queue_B_over_E"]) - 1.0) > 5e-4:
        raise RuntimeError("Candidate old state is not at replacement")

    ages = (
        float(old_parameters.age_start)
        + np.arange(int(old_parameters.J), dtype=float) * float(old_parameters.da)
    )
    stationary_age_mass = smoke.age_masses(stationary_g_pre)
    age_reweight = transition.acs_2007_age_reweight_diagnostic(
        stationary_age_mass,
        ages,
        E_old,
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
    prepared = continuation.PreparedModel(
        chain=chain,
        model=model,
        base_overrides=base,
        parameters=old_parameters,
        b_grid=b_grid,
        old_price=old_price,
        old_policy=old_policy,
        initial_g_pre_2007=initial_g_pre,
        initial_mass_2007=float(np.sum(initial_g_pre)),
        births_old=births_old,
        births_old_raw=births_old_raw,
        E_old=E_old,
        B_old=B_old,
        supply_rule=supply_rule,
        supply_normalization=supply_normalization,
        operator_gates=operator_gates,
        stationary_gates={
            "old_completed_fertility": float(normalization["completed_fertility"]),
            "old_entry_flow_E": E_old,
            "old_queue_mature_flow_B": B_old,
            "old_queue_B_over_E": B_old / E_old,
            "old_price": old_price,
        },
    )
    return prepared, {
        **renewal_state,
        "old_entry_flow_E": E_old,
        "old_queue_mature_flow_B": B_old,
        "outside_flow_M": float(renewal_state["outside_flow_M"]),
        "retention_rho": float(renewal_state["retention_rho"]),
        "outside_origin_entry_share": float(outside_origin_share),
        "old_fertility_normalization": normalization,
    }


def validate_candidate(payload: dict[str, Any], report: dict[str, Any]) -> tuple[dict[str, float], float]:
    if payload.get("schema") != "e5f_historical_four_group_candidate_v1":
        raise RuntimeError("Candidate schema is invalid")
    theta = {str(key): float(value) for key, value in dict(payload["theta"]).items()}
    source_theta = dict(report["best_candidate"]["theta"])
    if set(theta) != set(source_theta):
        raise RuntimeError("Candidate theta fields differ from the selected estimate")
    if any(not math.isfinite(value) for value in theta.values()):
        raise RuntimeError("Candidate theta contains a nonfinite value")
    delta = float(payload["psi_child_change_2023"])
    if not math.isfinite(delta) or not -1.5 <= delta <= 0.2:
        raise RuntimeError("Candidate preference change is outside its declared domain")
    return theta, delta


def main() -> None:
    args = parse_args()
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "candidate": args.candidate_json.resolve(),
        "report": args.selected_report.resolve(),
        "source": args.source.resolve(),
        "age": args.age_path.resolve(),
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    started = time.perf_counter()
    candidate = json.loads(paths["candidate"].read_text(encoding="utf-8"))
    report = json.loads(paths["report"].read_text(encoding="utf-8"))
    theta, delta = validate_candidate(candidate, report)
    report_sha = continuation.file_sha256(paths["report"])
    if str(candidate.get("selected_report_sha256")) != report_sha:
        raise RuntimeError("Candidate selected-report hash is stale")

    chain, model = transition.configure_sequential_model()
    live_bundle = calibration.code_fingerprint_contract(model)["bundle_sha256"]
    expected_bundle = report["scientific_contract"]["code_fingerprints"][
        "bundle_sha256"
    ]
    if live_bundle != expected_bundle:
        raise RuntimeError(
            f"Live scientific bundle differs from selected report: {live_bundle} != {expected_bundle}"
        )
    if str(candidate.get("scientific_bundle_sha256")) != str(live_bundle):
        raise RuntimeError("Candidate scientific-bundle hash is stale")
    target_system = calibration.e5_target_system()
    if target_system.fingerprint != report["scientific_contract"]["target_fingerprint"]:
        raise RuntimeError("Live target fingerprint differs from the selected report")
    if str(candidate.get("target_fingerprint")) != str(target_system.fingerprint):
        raise RuntimeError("Candidate target fingerprint is stale")
    outside_share = float(report["scientific_contract"]["outside_origin_entry_share"])
    prepared, renewal_state = prepare_candidate(
        chain=chain,
        model=model,
        theta=theta,
        initial_old_psi=float(report["best_candidate"]["old_psi_child"]),
        outside_origin_share=outside_share,
    )
    old_psi = float(
        renewal_state["old_fertility_normalization"]["psi_child"]
    )
    contracts = {
        "report_best": {
            "old_psi_child": old_psi,
            "new_psi_child": old_psi + delta,
        },
        "renewal_old_state": renewal_state,
    }
    shares = design.load_age_shares(paths["age"])
    population_indices = transition.census_household_target_indices()
    targets = [
        float(population_indices[index]) * shares[year]
        for index, year in enumerate(design.YEARS)
    ]
    profile, profile_rows = design.fit_group_profile(
        targets, design.aging_operator(), design.FOUR_GROUPS
    )
    path_rows, age_rows, result = smoke.run_design(
        label=str(candidate["candidate_id"]),
        prepared=prepared,
        contracts=contracts,
        demographic_profile=profile,
        targets=targets,
        market_tol=float(args.market_tol),
        market_max_iter=int(args.market_max_iter),
    )
    fit_rows = list(result.pop("fit_rows"))
    smoke.write_csv(outdir / "transition_path.csv", path_rows)
    smoke.write_csv(outdir / "age_distribution_fit.csv", age_rows)
    smoke.write_csv(outdir / "target_fit.csv", fit_rows)
    smoke.write_csv(outdir / "four_group_net_flow_parameters.csv", profile_rows)
    summary = {
        "status": "complete_historical_four_group_candidate",
        "candidate": candidate,
        "candidate_sha256": continuation.file_sha256(paths["candidate"]),
        "selected_report_sha256": report_sha,
        "source_sha256": continuation.file_sha256(paths["source"]),
        "age_path_sha256": continuation.file_sha256(paths["age"]),
        "scientific_bundle_sha256": live_bundle,
        "target_fingerprint": target_system.fingerprint,
        "old_psi_child": old_psi,
        "new_psi_child": old_psi + delta,
        "renewal_old_state": renewal_state,
        **result,
        "elapsed_seconds": time.perf_counter() - started,
    }
    continuation.write_json(outdir / "summary.json", summary)
    print(
        "HISTORICAL_FOUR_GROUP_CANDIDATE_COMPLETE "
        f"candidate={candidate['candidate_id']} "
        f"loss={summary['transition_loss_at_selected_parameters']:.9g} "
        f"elapsed={summary['elapsed_seconds']:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
