#!/usr/bin/env python3
"""Collect and gate a Torch E5F dated-transition calibration panel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--expected-tasks", type=int, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-target-set", required=True)
    parser.add_argument("--expected-target-fingerprint", required=True)
    parser.add_argument("--expected-code-bundle-sha256", required=True)
    parser.add_argument("--expected-model-profile", required=True)
    parser.add_argument("--expected-panel-seed", type=int, required=True)
    parser.add_argument("--expected-center-sha256", required=True)
    parser.add_argument("--expected-panel-design", required=True)
    parser.add_argument("--expected-local-radius", type=float, required=True)
    parser.add_argument("--expected-housing-supply-elasticity", type=float, default=None)
    parser.add_argument("--expected-tenure-choice-kappa", type=float, default=None)
    parser.add_argument("--require-complete", action="store_true")
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def target_measurement_provenance(recorded: Any, target_set: str) -> dict[str, Any]:
    if not isinstance(recorded, dict) or not recorded:
        raise RuntimeError("task summary omits target-measurement definitions")
    return {
        "contract": str(target_set),
        "recorded_model_timing": recorded,
        "authoritative_source_ledgers": {
            "fertility_housing_tenure": "docs/model/e5_target_review_20260724.md",
            "model_statistic_definitions": "docs/model/intergen_target_object_audit.md",
            "wealth_level": "docs/model/intergen_wealth_target_beta_audit_20260723.md",
        },
        "first_birth_rooms_target": {
            "status": "provisional_model_mapping",
            "contract_id": "psid_first_birth_rooms_household_aligned_sa_v1_20260817",
            "receipt": (
                "code/data/psid_followup_mar2026/output/"
                "sa_rooms_first_birth_household_aligned_v1/metadata.json"
            ),
        },
    }


def population_bridge_contract(recorded: Any) -> dict[str, Any]:
    """Return only the exogenous bridge objects that must agree across candidates.

    The realized stationary age mass, survival path, and reweighted masses depend
    on a candidate's household parameters.  They are diagnostics, not provenance,
    and therefore must not be used to reject an otherwise coherent panel.
    """
    if not isinstance(recorded, dict):
        raise RuntimeError("task summary omits the population-bridge contract")
    initial = recorded.get("initial_age_reweight")
    if not isinstance(initial, dict):
        raise RuntimeError("population bridge omits the initial-age contract")
    bridge_fields = (
        "status",
        "target_indices",
        "total_households",
        "age_distribution",
    )
    initial_fields = (
        "mapping",
        "model_age_grid",
        "sample",
        "source",
        "source_path",
        "source_receipt",
        "source_sha256",
        "source_year",
        "units",
        "period_years",
        "horizon_periods",
        "horizon_calendar_year",
        "future_entry_rule",
    )
    missing_bridge = [name for name in bridge_fields if name not in recorded]
    missing_initial = [name for name in initial_fields if name not in initial]
    if missing_bridge or missing_initial:
        raise RuntimeError(
            "population bridge lacks required exogenous fields: "
            f"bridge={missing_bridge}, initial={missing_initial}"
        )
    return {
        **{name: recorded[name] for name in bridge_fields},
        "initial_age_reweight": {
            name: initial[name] for name in initial_fields
        },
    }


def validate_renewal_accounting(summary: dict[str, Any]) -> None:
    """Fail closed on the replacement and old-state renewal algebra."""
    contract = dict(summary.get("renewal_accounting_contract") or {})
    old = dict(summary.get("renewal_accounting_old_state") or {})
    replacement = float(contract.get("replacement_fertility", math.nan))
    conversion = float(
        contract.get("effective_birth_to_household_conversion", math.nan)
    )
    if not math.isclose(replacement, 2.1, rel_tol=0.0, abs_tol=1e-14):
        raise RuntimeError(f"Invalid replacement-fertility contract: {replacement}")
    if not math.isclose(conversion, 1.0 / 2.1, rel_tol=0.0, abs_tol=1e-15):
        raise RuntimeError(f"Invalid birth-to-household conversion: {conversion}")
    if int(contract.get("birth_vintage_queue_waiting_slots", -1)) != 4:
        raise RuntimeError("Renewal contract must use four queue waiting slots")
    if int(contract.get("birth_to_entry_effect_lag_dates", -1)) != 5:
        raise RuntimeError("Renewal contract must report a five-date effect lag")
    if not math.isclose(
        float(contract.get("birth_to_entry_effect_lag_years", math.nan)),
        20.0,
        rel_tol=0.0,
        abs_tol=1e-14,
    ):
        raise RuntimeError("Renewal contract must report a 20-year effect lag")
    names = (
        "old_entry_flow_E",
        "old_queue_mature_flow_B",
        "old_queue_B_over_E",
        "old_renewal_residual",
        "outside_flow_M",
        "outside_origin_entry_share",
        "retention_rho",
    )
    values = {name: float(old.get(name, math.nan)) for name in names}
    if not all(math.isfinite(value) for value in values.values()):
        raise RuntimeError("Old renewal accounting contains nonfinite values")
    E = values["old_entry_flow_E"]
    B = values["old_queue_mature_flow_B"]
    ratio = values["old_queue_B_over_E"]
    outside = values["outside_flow_M"]
    share = values["outside_origin_entry_share"]
    retention = values["retention_rho"]
    if min(E, B, outside) <= 0.0 or not 0.0 < share < 1.0 or not 0.0 <= retention <= 1.0:
        raise RuntimeError("Old renewal accounting has invalid signs or shares")
    completed = float(summary["old_model_completed_fertility"])
    identities = (
        (ratio, B / E, "reported B/E", 2e-10),
        # The completed-fertility stock and top-code-adjusted renewal flow are
        # independently measured objects.  Match the driver's established
        # stationary-measurement gate rather than imposing a tighter collector-
        # only identity that valid tasks cannot satisfy in floating point.
        (ratio, completed / replacement, "fertility-to-replacement ratio", 2e-8),
        (outside, share * E, "outside-flow formula", 2e-10),
        (retention, (1.0 - share) * E / B, "retention formula", 2e-10),
        (
            values["old_renewal_residual"],
            E - outside - retention * B,
            "renewal residual",
            2e-10,
        ),
    )
    for reported, reconstructed, label, tolerance in identities:
        if not math.isclose(
            reported, reconstructed, rel_tol=0.0, abs_tol=tolerance
        ):
            raise RuntimeError(
                f"Old renewal {label} failed: {reported} versus {reconstructed}"
            )
    if abs(ratio - 1.0) > 5e-4:
        raise RuntimeError(f"Old renewal state is not at replacement: B/E={ratio}")


def validate_calibration_scope(summary: dict[str, Any]) -> None:
    """Keep estimation distinct from post-2023 and policy experiments."""
    if int(summary.get("post_2023_periods", -1)) != 0:
        raise RuntimeError("Calibration task extends beyond the 2023 target date")
    if str(summary.get("policy_case", "missing")) != "none":
        raise RuntimeError("Calibration task contains a policy experiment")


def validate_expected_contract(
    summary: dict[str, Any], panel: dict[str, Any], args: argparse.Namespace
) -> None:
    """Bind every collected task to the predeclared production contract."""
    profile = dict(summary.get("model_profile") or {})
    code = dict(summary.get("code_fingerprints") or {})
    checks = (
        (str(summary.get("source_sha256")), args.expected_source_sha256, "source hash"),
        (str(summary.get("target_set")), args.expected_target_set, "target set"),
        (
            str(summary.get("target_fingerprint")),
            args.expected_target_fingerprint,
            "target fingerprint",
        ),
        (
            str(code.get("bundle_sha256")),
            args.expected_code_bundle_sha256,
            "scientific code bundle",
        ),
        (
            str(profile.get("name")),
            args.expected_model_profile,
            "model profile",
        ),
        (str(panel.get("panel_seed")), str(args.expected_panel_seed), "panel seed"),
        (
            str(panel.get("center_sha256")),
            args.expected_center_sha256,
            "warm-start center hash",
        ),
        (
            str(panel.get("panel_design")),
            args.expected_panel_design,
            "panel design",
        ),
    )
    for recorded, expected, label in checks:
        if recorded != str(expected):
            raise RuntimeError(
                f"Unexpected {label}: recorded={recorded!r}, expected={expected!r}"
            )
    closure = dict(summary.get("external_closure_contract") or {})
    for expected, key, label in (
        (
            args.expected_housing_supply_elasticity,
            "housing_supply_elasticity",
            "housing-supply elasticity",
        ),
        (
            args.expected_tenure_choice_kappa,
            "tenure_choice_kappa",
            "tenure-choice dispersion",
        ),
    ):
        if expected is None:
            continue
        recorded = float(closure.get(key, math.nan))
        if not math.isfinite(recorded) or not math.isclose(
            recorded, float(expected), rel_tol=0.0, abs_tol=1e-14
        ):
            raise RuntimeError(
                f"Unexpected {label}: recorded={recorded}, expected={expected}"
            )
    radius = float(panel.get("local_radius", math.nan))
    if not math.isfinite(radius) or not math.isclose(
        radius,
        float(args.expected_local_radius),
        rel_tol=0.0,
        abs_tol=1e-14,
    ):
        raise RuntimeError(
            f"Unexpected local radius: recorded={radius}, "
            f"expected={args.expected_local_radius}"
        )


def validate_target_fit(
    rows: list[dict[str, str]], summary: dict[str, Any], candidate: dict[str, Any]
) -> dict[str, float]:
    """Reconstruct the target fingerprint and every row of the SMM loss."""
    expected_count = int(summary["target_count"])
    if len(rows) != expected_count:
        raise RuntimeError("target-fit row count mismatch")
    names = [str(row.get("moment", "")) for row in rows]
    if len(set(names)) != expected_count or any(not name for name in names):
        raise RuntimeError("target-fit moment names are missing or duplicated")
    payload: list[list[Any]] = []
    model_by_name: dict[str, float] = {}
    loss = 0.0
    expected_candidate = str(candidate["candidate"])
    for row in rows:
        if str(row.get("candidate")) != expected_candidate:
            raise RuntimeError("target-fit candidate label mismatch")
        values = {
            field: float(row.get(field, math.nan))
            for field in (
                "target",
                "model",
                "gap",
                "weight",
                "loss_contribution",
                "standardized_gap",
            )
        }
        if not all(math.isfinite(value) for value in values.values()):
            raise RuntimeError(f"Nonfinite target-fit row: {row['moment']}")
        if values["weight"] <= 0.0 or values["loss_contribution"] < 0.0:
            raise RuntimeError(f"Invalid target-fit weight/loss: {row['moment']}")
        gap = values["model"] - values["target"]
        contribution = values["weight"] * gap**2
        standardized = math.copysign(math.sqrt(contribution), gap)
        identities = (
            (values["gap"], gap, "gap"),
            (values["loss_contribution"], contribution, "loss contribution"),
            (values["standardized_gap"], standardized, "standardized gap"),
        )
        for recorded, reconstructed, label in identities:
            if not math.isclose(
                recorded, reconstructed, rel_tol=2e-12, abs_tol=2e-12
            ):
                raise RuntimeError(
                    f"Target-fit {label} failed for {row['moment']}: "
                    f"{recorded} versus {reconstructed}"
                )
        payload.append([row["moment"], values["target"], values["weight"]])
        model_by_name[row["moment"]] = values["model"]
        loss += contribution
    fingerprint = hashlib.sha256(
        json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()
    if fingerprint != str(summary["target_fingerprint"]):
        raise RuntimeError(
            "Target-fit rows do not reproduce the recorded target fingerprint"
        )
    if not math.isclose(
        loss,
        float(candidate["transition_loss"]),
        rel_tol=2e-12,
        abs_tol=2e-10,
    ):
        raise RuntimeError(
            f"Target-fit loss mismatch: {loss} versus {candidate['transition_loss']}"
        )
    return model_by_name


def validate_dated_housing_ledger(
    task_dir: Path,
    summary: dict[str, Any],
    candidate: dict[str, Any],
    fit_models: dict[str, float],
) -> None:
    """Verify the dated 2019--2023 branch used by the active PSID target row."""
    contract = dict(summary.get("dated_measurement_contract") or {})
    expected_contract = {
        "first_birth_housing_horizon_periods": 1,
        "first_birth_housing_horizon_years": 4.0,
        "first_birth_housing_terminal_origin_year": 2019,
        "first_birth_housing_terminal_destination_year": 2023,
    }
    for field, expected in expected_contract.items():
        recorded = float(contract.get(field, math.nan))
        if not math.isfinite(recorded) or recorded != float(expected):
            raise RuntimeError(
                f"Invalid dated housing contract {field}: {recorded}"
            )
    case_dir = task_dir / "cases" / str(candidate["candidate"])
    ledger_path = case_dir / "dated_first_birth_housing_did.csv"
    if not ledger_path.exists():
        raise RuntimeError(f"Missing dated housing ledger: {ledger_path}")
    rows = read_csv(ledger_path)
    if len(rows) != 5 or [int(row["period"]) for row in rows] != list(range(5)):
        raise RuntimeError("Dated housing ledger must contain periods 0,...,4")
    if any(str(row.get("census_age_bridge_applied")) != "False" for row in rows):
        raise RuntimeError("Census age bridge entered the matched housing branch")
    if rows[0].get("status") != "same_policy_one_period_diagnostic_no_prior_date":
        raise RuntimeError("Unexpected period-0 housing diagnostic")
    for period, row in enumerate(rows[1:], start=1):
        if row.get("status") != "complete_one_period_dated_matched_branch_did":
            raise RuntimeError(f"Incomplete dated housing branch at period {period}")
        if int(row["origin_period"]) != period - 1 or int(row["destination_period"]) != period:
            raise RuntimeError(f"Wrong dated housing horizon at period {period}")
        values = {
            field: float(row.get(field, math.nan))
            for field in (
                "origin_mass",
                "destination_mass",
                "survival_loss",
                "origin_price",
                "destination_price",
                "treated_feasibility_projection_mass",
                "control_feasibility_projection_mass",
                "treated_continuation_births",
                "treated_mean_housing",
                "control_mean_housing",
                "housing_response",
            )
        }
        if not all(math.isfinite(value) for value in values.values()):
            raise RuntimeError(f"Nonfinite dated housing ledger at period {period}")
        if min(values["origin_mass"], values["destination_mass"]) <= 0.0:
            raise RuntimeError(f"Nonpositive branch mass at period {period}")
        if min(values["treated_mean_housing"], values["control_mean_housing"]) <= 0.0:
            raise RuntimeError(f"Nonpositive branch housing at period {period}")
        if min(
            values["treated_feasibility_projection_mass"],
            values["control_feasibility_projection_mass"],
            values["treated_continuation_births"],
        ) < 0.0:
            raise RuntimeError(f"Negative branch diagnostic at period {period}")
        if (
            values["survival_loss"] < 0.0
            or values["survival_loss"] >= values["origin_mass"]
        ):
            raise RuntimeError(f"Invalid branch survival loss at period {period}")
        if not math.isclose(
            values["origin_mass"] - values["survival_loss"],
            values["destination_mass"],
            rel_tol=0.0,
            abs_tol=2e-12,
        ):
            raise RuntimeError(f"Branch mass mismatch at period {period}")
        reconstructed_response = (
            values["treated_mean_housing"] - values["control_mean_housing"]
        )
        if not math.isclose(
            values["housing_response"],
            reconstructed_response,
            rel_tol=0.0,
            abs_tol=2e-12,
        ):
            raise RuntimeError(f"Branch housing response mismatch at period {period}")
    terminal = rows[-1]
    if int(terminal["origin_period"]) != 3 or int(terminal["destination_period"]) != 4:
        raise RuntimeError("Terminal housing row is not the 2019--2023 branch")
    response = float(terminal["housing_response"])
    expected_response = float(fit_models["housing_increment_0to1"])
    summary_response = float(candidate["terminal_first_birth_housing_response"])
    if not (
        math.isclose(response, expected_response, rel_tol=0.0, abs_tol=2e-12)
        and math.isclose(response, summary_response, rel_tol=0.0, abs_tol=2e-12)
    ):
        raise RuntimeError(
            "Terminal housing response differs across dated ledger, target fit, "
            "and candidate summary"
        )


def main() -> None:
    args = parse_args()
    results_dir = args.results_dir.resolve()
    report_dir = results_dir / "report"
    expected = int(args.expected_tasks)
    summaries: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    contracts: set[tuple[Any, ...]] = set()
    provenance_contracts: set[str] = set()

    for task in range(1, expected + 1):
        task_dir = results_dir / f"task_{task:03d}"
        summary_path = task_dir / "summary.json"
        if not summary_path.exists():
            failure_path = task_dir / "failure.json"
            if failure_path.exists():
                failure = read_json(failure_path)
                failures.append(
                    {
                        "task_id": task,
                        "reason": "classified_invalid_candidate",
                        "error_type": failure.get("error_type"),
                        "error": failure.get("error"),
                    }
                )
            else:
                failures.append({"task_id": task, "reason": "missing_summary"})
            continue
        try:
            summary = read_json(summary_path)
            candidate = dict(summary["best_candidate"])
            panel = dict(summary["panel_design"])
            profile = dict(summary.get("model_profile") or {})
            provenance = {
                "population_bridge": population_bridge_contract(
                    summary["population_bridge"]
                ),
                "population_validation_status": summary[
                    "population_validation_status"
                ],
            }
            provenance["target_measurements"] = target_measurement_provenance(
                summary["target_measurements"], str(summary["target_set"])
            )
            validate_renewal_accounting(summary)
            validate_calibration_scope(summary)
            validate_expected_contract(summary, panel, args)
            provenance_contracts.add(
                json.dumps(provenance, sort_keys=True, separators=(",", ":"))
            )
            contract = (
                str(summary["source_sha256"]),
                str(summary["target_set"]),
                str(summary["target_fingerprint"]),
                int(summary["target_count"]),
                int(panel["panel_seed"]),
                str(profile.get("name", "e5f-floor")),
                str(profile.get("profile_id", "none")),
                str(profile.get("permanent_income_log_variance", "none")),
                int(profile.get("income_state_count", 5)),
                str(profile.get("first_birth_fixed_cost_semantics", "none")),
                str(panel.get("center_sha256", "legacy_unspecified")),
                str(
                    (summary.get("code_fingerprints") or {}).get(
                        "bundle_sha256", "legacy_unspecified"
                    )
                ),
                float(summary["outside_origin_entry_share"]),
                float(summary["old_completed_fertility_reference"]),
                json.dumps(
                    summary["renewal_accounting_contract"],
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                json.dumps(
                    summary["dated_measurement_contract"],
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                int(summary["post_2023_periods"]),
                str(summary["policy_case"]),
                json.dumps(
                    summary.get("external_closure_contract") or {},
                    sort_keys=True,
                    separators=(",", ":"),
                ),
            )
            contracts.add(contract)
            if int(panel["task_id"]) != task or int(panel["panel_size"]) != expected:
                raise RuntimeError("panel task metadata mismatch")
            fit = read_csv(task_dir / "target_fit_long.csv")
            fit_models = validate_target_fit(fit, summary, candidate)
            validate_dated_housing_ledger(
                task_dir,
                summary,
                candidate,
                fit_models,
            )
            for row in fit:
                row["task_id"] = task
                target_rows.append(row)
            candidate["task_id"] = task
            candidate["design"] = panel["design"]
            candidate["model_profile"] = str(
                (summary.get("model_profile") or {}).get("name", "e5f-floor")
            )
            candidate["stationary_measurement_max_abs_gap"] = float(
                summary["stationary_measurement_max_abs_gap"]
            )
            candidate["valid"] = bool(
                math.isfinite(float(candidate["transition_loss"]))
                and float(candidate["max_market_residual"]) <= 2e-4
                and abs(float(candidate["population_target_gap"])) <= 2e-10
                and float(candidate["max_mass_residual"]) <= 2e-10
                and float(candidate["stationary_measurement_max_abs_gap"]) <= 2e-8
            )
            candidate_rows.append(candidate)
            summaries.append(summary)
        except Exception as error:
            failures.append(
                {"task_id": task, "reason": f"invalid_artifact: {type(error).__name__}: {error}"}
            )

    if len(contracts) > 1:
        raise RuntimeError(f"Mixed source/target/panel contracts: {sorted(contracts)}")
    if len(provenance_contracts) > 1:
        raise RuntimeError("Mixed population or target-measurement provenance")
    valid = [row for row in candidate_rows if bool(row["valid"])]
    best = (
        dict(min(valid, key=lambda row: float(row["transition_loss"])))
        if valid
        else None
    )
    if best is not None:
        selected_summary = next(
            summary
            for summary in summaries
            if int(summary["panel_design"]["task_id"]) == int(best["task_id"])
        )
        best.update(
            source=selected_summary["source"],
            source_sha256=selected_summary["source_sha256"],
            target_set=selected_summary["target_set"],
            target_fingerprint=selected_summary["target_fingerprint"],
            target_count=selected_summary["target_count"],
            model_profile=selected_summary.get("model_profile"),
            income_profile_gates=selected_summary.get("income_profile_gates"),
            code_fingerprints=selected_summary.get("code_fingerprints"),
            external_closure_contract=selected_summary.get(
                "external_closure_contract"
            ),
            population_bridge=selected_summary["population_bridge"],
            population_validation_status=selected_summary[
                "population_validation_status"
            ],
            target_measurements=target_measurement_provenance(
                selected_summary["target_measurements"],
                str(selected_summary["target_set"]),
            ),
            outside_origin_entry_share=float(
                selected_summary["outside_origin_entry_share"]
            ),
            old_completed_fertility_reference=float(
                selected_summary["old_completed_fertility_reference"]
            ),
            renewal_accounting_contract=selected_summary[
                "renewal_accounting_contract"
            ],
            renewal_accounting_old_state=selected_summary[
                "renewal_accounting_old_state"
            ],
            dated_measurement_contract=selected_summary[
                "dated_measurement_contract"
            ],
            old_psi_child=selected_summary["old_psi_child"],
            old_completed_fertility=selected_summary["old_model_completed_fertility"],
        )
    candidate_rows.sort(key=lambda row: float(row["transition_loss"]))
    write_csv(report_dir / "all_candidates.csv", candidate_rows)
    write_csv(report_dir / "all_target_fit.csv", target_rows)
    if best is not None:
        write_json(report_dir / "best_candidate.json", best)
        write_json(
            report_dir / "calibration_provenance.json",
            {
                "population_bridge": best["population_bridge"],
                "population_validation_status": best[
                    "population_validation_status"
                ],
                "target_measurements": best["target_measurements"],
                "renewal_accounting_contract": best[
                    "renewal_accounting_contract"
                ],
                "renewal_accounting_old_state": best[
                    "renewal_accounting_old_state"
                ],
                "dated_measurement_contract": best[
                    "dated_measurement_contract"
                ],
            },
        )
        best_fit = [
            row for row in target_rows if int(row["task_id"]) == int(best["task_id"])
        ]
        write_csv(report_dir / "best_target_fit.csv", best_fit)
        selected_case_dir = (
            results_dir
            / f"task_{int(best['task_id']):03d}"
            / "cases"
            / str(best["candidate"])
        )
        portable_artifacts = {
            "transition_path.csv": "best_transition_path.csv",
            "dated_model_moments.csv": "best_dated_model_moments.csv",
            "dated_period_fertility.csv": "best_dated_period_fertility.csv",
            "dated_period_fertility_by_age.csv": "best_dated_period_fertility_by_age.csv",
            "dated_first_birth_housing_did.csv": "best_dated_first_birth_housing_did.csv",
            "dated_cohort_timing_ledgers.json": "best_dated_cohort_timing_ledgers.json",
        }
        for source_name, destination_name in portable_artifacts.items():
            source_path = selected_case_dir / source_name
            if not source_path.exists():
                raise RuntimeError(
                    f"Selected task omits required dated artifact: {source_path}"
                )
            shutil.copy2(source_path, report_dir / destination_name)
    report = {
        "status": "complete" if not failures else "partial",
        "results_dir": str(results_dir),
        "expected_tasks": expected,
        "completed_tasks": len(candidate_rows),
        "valid_tasks": len(valid),
        "failed_or_missing_tasks": failures,
        "contract": list(next(iter(contracts))) if contracts else None,
        "outside_origin_entry_share": (
            float(best["outside_origin_entry_share"]) if best is not None else None
        ),
        "old_completed_fertility_reference": (
            float(best["old_completed_fertility_reference"])
            if best is not None
            else None
        ),
        "renewal_accounting_contract": (
            best.get("renewal_accounting_contract") if best is not None else None
        ),
        "renewal_accounting_old_state": (
            best.get("renewal_accounting_old_state") if best is not None else None
        ),
        "dated_measurement_contract": (
            best.get("dated_measurement_contract") if best is not None else None
        ),
        "code_fingerprints": (
            best.get("code_fingerprints") if best is not None else None
        ),
        "external_closure_contract": (
            best.get("external_closure_contract") if best is not None else None
        ),
        "model_profile": best.get("model_profile") if best is not None else None,
        "population_bridge": (
            best.get("population_bridge") if best is not None else None
        ),
        "population_validation_status": (
            best.get("population_validation_status") if best is not None else None
        ),
        "target_measurements": (
            best.get("target_measurements") if best is not None else None
        ),
        "best_candidate": best,
    }
    write_json(report_dir / "summary.json", report)
    print(
        "TRANSITION_PANEL_COLLECTED "
        f"completed={len(candidate_rows)}/{expected} valid={len(valid)} "
        f"best_loss={float(best['transition_loss']) if best else math.inf:.6f}",
        flush=True,
    )
    if args.require_complete and failures:
        raise SystemExit(2)
    if best is None:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
