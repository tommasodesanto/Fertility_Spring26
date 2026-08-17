#!/usr/bin/env python3
"""Focused contract tests for the canonical no-policy report builder."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import build_e5f_no_policy_transition_report as report


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def artifact_receipt(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": report.file_sha256(path),
    }


def write_sidecar(path: Path) -> None:
    path.with_suffix(".sha256").write_text(
        f"{report.file_sha256(path)}  {path.name}\n"
    )


def target_rows() -> list[dict[str, object]]:
    rows = []
    for index in range(12):
        target = 1.0 + index
        gap = 0.01 * (index + 1)
        weight = 1.0 + index
        rows.append(
            {
                "moment": f"moment_{index + 1:02d}",
                "target": target,
                "model": target + gap,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap * gap,
                "standardized_gap": math.sqrt(weight) * gap,
            }
        )
    return rows


def make_figure_pair(root: Path, stem: str) -> list[dict[str, object]]:
    records = []
    for suffix in ("png", "pdf"):
        path = root / f"{stem}.{suffix}"
        path.write_bytes(f"fixture {stem} {suffix}".encode())
        records.append(artifact_receipt(path))
    return records


def build_fixture(root: Path) -> argparse.Namespace:
    plan_dir = root / "plan"
    selected_dir = root / "selected"
    historical_dir = root / "historical"
    continuation_dir = root / "continuation"
    for directory in (plan_dir, selected_dir, historical_dir, continuation_dir):
        directory.mkdir(parents=True)

    renewal = {
        "replacement_fertility": 2.1,
        "effective_birth_to_household_conversion": 1.0 / 2.1,
        "birth_to_entry_effect_lag_years": 20.0,
        "birth_to_entry_effect_lag_dates": 5,
        "birth_vintage_queue_waiting_slots": 4,
        "birth_measure": "topcode_adjusted_birth_children",
    }
    dated = {"dated_maintained_moments": "only 2023 enters the objective"}
    population_bridge = {
        "status": "externally_normalized_not_estimated",
        "target_indices": {str(i): 1.0 + 0.02 * i for i in range(5)},
        "total_households": "Census HH-3",
        "age_distribution": "ACS household heads",
        "initial_age_reweight": {
            "mapping": "four-year bins",
            "model_age_grid": [18, 22],
            "sample": "national ACS",
            "source": "fixture",
            "source_path": "/fixture/acs.dta",
            "source_receipt": {"records": 1},
            "source_sha256": "d" * 64,
            "source_year": 2007,
            "units": "households",
            "period_years": 4.0,
            "horizon_periods": 4,
            "horizon_calendar_year": 2023.0,
            "future_entry_rule": "old queue",
        },
    }
    scientific = {
        "source": "/fixture/source.json",
        "source_sha256": "b" * 64,
        "source_metadata": {"status": "fixture"},
        "target_set": "fixture_target_set",
        "target_fingerprint": "a" * 64,
        "target_count": 12,
        "code_fingerprints": {"bundle_sha256": "c" * 64, "files": {}},
        "model_profile": {"name": "e5f-income-entry"},
        "income_profile_gates": {"income_state_count": 15},
        "target_measurements": {"contract": "fixture"},
        "renewal_accounting_contract": renewal,
        "dated_measurement_contract": dated,
        "population_bridge": population_bridge,
        "population_validation_status": "imposed bridge",
        "outside_origin_entry_share": 0.169,
        "old_completed_fertility_reference": 2.1,
        "policy_case": "none",
        "post_2023_periods": 0,
        "stationary_free_parameter_count": 10,
        "transition_free_parameter_count": 10,
    }
    domains = [
        {
            "name": f"parameter_{index:02d}",
            "lower": 0.0,
            "upper": 1.0,
            "transform": "linear",
        }
        for index in range(10)
    ]
    targets = target_rows()
    target_contract = [
        {
            "moment_index": index,
            "moment": row["moment"],
            "target": row["target"],
            "weight": row["weight"],
        }
        for index, row in enumerate(targets)
    ]
    matrix = np.vstack((np.eye(10), np.zeros((2, 10))))
    jacobian_rows = []
    for index, target in enumerate(target_contract):
        row: dict[str, object] = {
            "moment_index": index,
            "moment": target["moment"],
            "anchor_weighted_residual": 0.0,
        }
        for dimension, domain in enumerate(domains):
            row[domain["name"]] = matrix[index, dimension]
        jacobian_rows.append(row)
    write_csv(plan_dir / "weighted_jacobian.csv", jacobian_rows)
    singular = np.linalg.svd(matrix, compute_uv=False)
    write_csv(
        plan_dir / "singular_values.csv",
        [
            {
                "singular_value_index": index,
                "singular_value": value,
                "relative_to_largest": value / singular[0],
            }
            for index, value in enumerate(singular)
        ],
    )
    side = [
        {
            "dimension": index,
            "parameter": domain["name"],
            "minus_task": 2 + 2 * index,
            "plus_task": 3 + 2 * index,
            "forward_norm": 1.0,
            "backward_norm": 1.0,
            "forward_backward_dot": 1.0,
            "side_consistency_cosine": 1.0,
            "relative_side_difference": 0.0,
            "frozen_for_step": False,
            "freeze_rule": "forward_backward_dot_nonpositive",
        }
        for index, domain in enumerate(domains)
    ]
    write_csv(plan_dir / "side_consistency.csv", side)
    write_csv(
        plan_dir / "weak_directions.csv",
        [
            {
                "weak_order": order,
                "singular_value_index": 10 - order,
                "singular_value": 1.0,
                "dimension": dimension,
                "parameter": domain["name"],
                "loading": 1.0 if dimension == 10 - order else 0.0,
                "squared_loading_share": 1.0 if dimension == 10 - order else 0.0,
            }
            for order in range(1, 4)
            for dimension, domain in enumerate(domains)
        ],
    )
    plan = {
        "schema": "e5f_transition_ridge_refinement_plan_v1",
        "scientific_contract": scientific,
        "scientific_contract_sha256": report.object_sha256(scientific),
        "renewal_contract_sha256": report.object_sha256(renewal),
        "dated_contract_sha256": report.object_sha256(dated),
        "target_contract": target_contract,
        "target_contract_sha256": report.object_sha256(target_contract),
        "domain": domains,
        "jacobian_diagnostics": {
            "full": {
                "shape": [12, 10],
                "singular_values": singular.tolist(),
                "relative_ranks": {"relative_0.001": 10},
                "condition_number": 1.0,
            },
            "identification_gate": {"passed": True},
            "frozen_dimensions": [],
            "active_dimensions": list(range(10)),
            "side_consistency": side,
        },
    }
    plan_path = plan_dir / "candidate_plan.json"
    write_json(plan_path, plan)
    write_sidecar(plan_path)
    plan_sha = report.file_sha256(plan_path)

    theta = {domain["name"]: 0.5 for domain in domains}
    parameter_rows = [
        {
            "dimension": index,
            "parameter": domain["name"],
            "estimate": 0.5,
            "lower": 0.0,
            "upper": 1.0,
            "transform": "linear",
            "unit_coordinate": 0.5,
            "distance_to_nearest_unit_bound": 0.5,
            "near_bound_within_0.02": False,
            "status": "estimated_in_dated_transition_refinement",
        }
        for index, domain in enumerate(domains)
    ]
    loss = math.fsum(float(row["loss_contribution"]) for row in targets)
    selection = {
        "schema": "e5f_transition_ridge_selection_v1",
        "selected_candidate_id": 1,
        "selected_candidate_sha256": "e" * 64,
        "transition_loss": loss,
        "theta": theta,
    }
    selection_path = selected_dir / "selection.json"
    write_json(selection_path, selection)
    write_sidecar(selection_path)
    selection_sha = report.file_sha256(selection_path)
    best = {
        "candidate": "task_001",
        "theta": theta,
        "old_psi_child": 0.25,
        "new_psi_child": -0.1,
        "max_market_residual": 1e-6,
        "max_mass_residual": 1e-12,
        "population_target_gap": 1e-12,
        "terminal_synthetic_childless_consistency_gap": 1e-12,
        "renewal_accounting_contract": renewal,
        "renewal_accounting_old_state": {
            "old_entry_flow_E": 0.1,
            "old_queue_mature_flow_B": 0.1,
        },
    }
    selected_summary = {
        "schema": "e5f_transition_ridge_refinement_report_v1",
        "status": "complete_refinement_with_two_independent_identity_repeats",
        "promotion_eligible": True,
        "plan_sha256": plan_sha,
        "selection_sha256": selection_sha,
        "scientific_contract": scientific,
        "scientific_contract_sha256": report.object_sha256(scientific),
        "renewal_contract_sha256": report.object_sha256(renewal),
        "dated_contract_sha256": report.object_sha256(dated),
        "selected_candidate_id": 1,
        "selected_candidate_sha256": "e" * 64,
        "selected_transition_loss": loss,
        "selected_validation_execution_identity": "validation_identity",
        "repeat_gate": {
            "required": 2,
            "completed": 2,
            "execution_identities": ["repeat_identity_1", "repeat_identity_2"],
            "identity_gate_passed": True,
            "numerical_identity_gate_passed": True,
        },
        "eligibility_tolerances": {
            "max_market_residual": 2e-4,
            "population_target_gap": 2e-10,
            "max_mass_residual": 2e-10,
            "stationary_measurement_max_abs_gap": 2e-8,
            "terminal_synthetic_childless_consistency_gap": 2e-10,
        },
        "target_fit_table": targets,
        "estimated_parameters": parameter_rows,
        "best_candidate": best,
    }
    write_csv(selected_dir / "selected_target_fit.csv", targets)
    write_csv(selected_dir / "selected_parameters_and_bounds.csv", parameter_rows)
    write_csv(selected_dir / "selected_transition_path.csv", [{"period": 0, "value": 1.0}])
    write_csv(selected_dir / "selected_dated_model_moments.csv", [{"period": 0, "value": 1.0}])
    write_csv(selected_dir / "selected_dated_period_fertility.csv", [{"period": 0, "value": 1.0}])
    write_csv(selected_dir / "selected_dated_period_fertility_by_age.csv", [{"period": 0, "value": 1.0}])
    write_json(selected_dir / "selected_dated_cohort_timing_ledgers.json", {"period": 0})
    selected_summary_path = selected_dir / "summary.json"
    write_json(selected_summary_path, selected_summary)
    write_sidecar(selected_summary_path)
    target_provenance_path = root / "target_provenance.csv"
    write_csv(
        target_provenance_path,
        [
            {
                "moment": row["moment"],
                "target": row["target"],
                "authoritative_builder": f"code/data/builder_{index:02d}.py",
                "sample_geography_vintage": "U.S. household sample, fixture vintage",
                "estimator_measurement": "fixture measurement stated in levels",
                "fixed_effects": "not applicable",
                "clustering": "not applicable",
                "uncertainty_type": "declared calibration scale",
                "standard_error_or_scale": 1.0 / math.sqrt(float(row["weight"])),
                "weight": row["weight"],
                "status": ("measured" if index < 9 else "provisional"),
                "caveat": "Fixture provenance; replace with the authoritative empirical receipt.",
            }
            for index, row in enumerate(targets)
        ],
    )

    repeat_paths = []
    repeat_summary_hash = None
    for replicate_id in (1, 2):
        task_dir = root / f"repeat_{replicate_id:03d}"
        task_dir.mkdir()
        task_best = {
            "candidate": "task_001",
            "theta": theta,
            "new_psi_child": -0.1,
            "transition_loss": loss,
            "max_market_residual": 1e-6,
            "max_mass_residual": 1e-12,
            "population_target_gap": 1e-12,
            "terminal_synthetic_childless_consistency_gap": 1e-12,
        }
        task_summary = {
            **scientific,
            "status": "complete_transition_calibration_panel_task",
            "old_psi_child": 0.25,
            "old_model_completed_fertility": 2.1,
            "best_candidate": task_best,
            "stationary_measurement_max_abs_gap": 1e-12,
        }
        write_json(task_dir / "summary.json", task_summary)
        if replicate_id == 1:
            repeat_summary_hash = report.file_sha256(task_dir / "summary.json")
        write_csv(task_dir / "target_fit_long.csv", targets)
        wrapper = {
            "schema": "e5f_transition_ridge_task_contract_v1",
            "mode": "repeat",
            "replicate_id": replicate_id,
            "candidate_id": 1,
            "candidate_sha256": "e" * 64,
            "plan_sha256": plan_sha,
            "scientific_contract_sha256": report.object_sha256(scientific),
            "renewal_contract_sha256": report.object_sha256(renewal),
            "dated_contract_sha256": report.object_sha256(dated),
            "selection_sha256": selection_sha,
            "execution_identity": f"repeat_identity_{replicate_id}",
        }
        write_json(task_dir / "refinement_task_contract.json", wrapper)
        repeat_paths.append(task_dir)

    comparison = []
    for year in report.EXPECTED_YEARS:
        comparison.append(
            {
                "calendar_year": year,
                "series_id": "household_population_index",
                "series_label": "Household population",
                "classification": "imposed_bridge_matched_by_construction_not_fit",
                "data_value": 1.0,
                "model_value": 1.0,
            }
        )
        for series in ("market_clearing_residual", "mass_accounting_residual"):
            comparison.append(
                {
                    "calendar_year": year,
                    "series_id": series,
                    "series_label": series,
                    "classification": "numerical_residual_audit_not_fit",
                    "data_value": 0.0,
                    "model_value": 1e-12,
                }
            )
    comparison.append(
        {
            "calendar_year": 2023,
            "series_id": "period_tfr_explicit_states",
            "series_label": "Period TFR",
            "classification": "untargeted_holdout",
            "data_value": 1.6,
            "model_value": 1.7,
        }
    )
    for row in targets:
        comparison.append(
            {
                "calendar_year": 2023,
                "series_id": f"terminal_target__{row['moment']}",
                "series_label": row["moment"],
                "classification": "terminal_maintained_target",
                "data_value": row["target"],
                "model_value": row["model"],
            }
        )
    comparison_path = historical_dir / "five_date_model_data_comparison.csv"
    fit_path = historical_dir / "five_date_fit_statistics.csv"
    write_csv(comparison_path, comparison)
    fit_specs = [
        ("household_population_index", "imposed_bridge_matched_by_construction_not_fit", "computed"),
        ("period_tfr_explicit_states", "descriptive_noncomparable_holdout", "withheld_noncomparable_measurement"),
        ("period_tfr_topcode_adjusted_sensitivity", "descriptive_noncomparable_holdout_sensitivity", "withheld_noncomparable_measurement"),
        ("period_first_birth_mean_age", "untargeted_holdout", "computed"),
        ("period_first_birth_share_age30plus", "untargeted_holdout", "computed"),
        ("aggregate_ownership_rate_18_85", "untargeted_holdout", "computed"),
        ("mean_rooms_literal_18_85", "untargeted_holdout", "computed"),
        ("real_house_price_index", "untargeted_holdout", "computed"),
        ("real_rent_cpi_proxy_index", "untargeted_holdout_data_only", "not_available_data_only"),
        ("price_to_rent_index", "descriptive_noncomparable_diagnostic", "withheld_noncomparable_measurement"),
        ("population_bridge_gap", "numerical_residual_audit_not_fit", "computed"),
        ("age_bridge_max_gap", "numerical_residual_audit_not_fit", "computed"),
        ("market_clearing_residual", "numerical_residual_audit_not_fit", "computed"),
        ("mass_accounting_residual", "numerical_residual_audit_not_fit", "computed"),
    ]
    write_csv(
        fit_path,
        [
            {
                "series_id": series,
                "classification": classification,
                "evaluation_years": "2007;2011;2015;2019;2023",
                "scalar_statistics_status": status,
                "rmse_evaluation_window": "" if status != "computed" else 0.1,
                "notes": "fixture",
            }
            for series, classification, status in fit_specs
        ],
    )
    historical_figures = {
        name: make_figure_pair(historical_dir, name)
        for name in (
            "imposed_inputs_matched_by_construction",
            "untargeted_fertility_housing_validation",
            "untargeted_price_rent_validation",
            "numerical_residual_audit",
        )
    }
    historical_inputs = {
        "case_summary": {
            "path": str(repeat_paths[0] / "summary.json"),
            "sha256": repeat_summary_hash,
        }
    }
    historical_bundle = "\n".join(
        f"{name}:{record['sha256']}"
        for name, record in sorted(historical_inputs.items())
    )
    historical_provenance = {
        "status": "complete_e5f_2007_2023_historical_validation_packet",
        "calendar_years": list(report.EXPECTED_YEARS),
        "inputs": historical_inputs,
        "input_bundle_sha256": hashlib.sha256(historical_bundle.encode()).hexdigest(),
        "model_contract": {
            "source_sha256": scientific["source_sha256"],
            "target_fingerprint": scientific["target_fingerprint"],
            "code_bundle_sha256": scientific["code_fingerprints"]["bundle_sha256"],
            "target_set": scientific["target_set"],
            "renewal_accounting_contract": renewal,
            "target_and_loss_reconciliation": {
                "target_count": 12,
                "target_fingerprint": scientific["target_fingerprint"],
                "recomputed_loss": loss,
                "saved_loss": loss,
            },
            "no_policy_historical_path": True,
        },
        "outputs": {
            "comparison_csv": artifact_receipt(comparison_path),
            "fit_statistics_csv": artifact_receipt(fit_path),
            "figures": historical_figures,
        },
        "validation_gates": {"all_fixture_gates": True},
    }
    historical_provenance_path = historical_dir / "classification_and_provenance.json"
    write_json(historical_provenance_path, historical_provenance)

    continuation_artifact_names = [
        "README.md",
        "summary.json",
        "paired_continuation_path.csv",
        "shared_2007_2023_history.csv",
        "history_reproduction_audit.csv",
        "outside_share_invariance_audit.csv",
        "closed_stationary_schedule.csv",
        "closed_stationary_schedule_progress.csv",
        "closed_stationary_endpoint.json",
        "continuation_progress.csv",
        "latest_completed_period.json",
        "latest_endpoint_search.json",
        "open_stationary_endpoint.json",
        "open_endpoint/stationary_endpoint_search.csv",
        "open_endpoint/stationary_open_endpoint.csv",
        "paired_continuation_levels.png",
        "paired_continuation_levels.pdf",
        "paired_renewal_diagnostics.png",
        "paired_renewal_diagnostics.pdf",
        "closed_stationary_renewal_schedule.png",
        "closed_stationary_renewal_schedule.pdf",
    ]
    start = {
        "policy_case": "none",
        "calendar_year": 2023,
        "adult_population": 1.2,
        "population_index_2023": 1.0,
        "asset_price_index": 1.05,
        "topcode_adjusted_births_per_adult": 0.1,
        "owner_rate": 0.65,
    }
    closed_last = {
        **start,
        "calendar_year": 2183,
        "population_index_2023": 0.7,
        "asset_price_index": 0.8,
        "topcode_adjusted_births_per_adult": 0.11,
        "owner_rate": 0.62,
    }
    open_last = {
        **closed_last,
        "population_index_2023": 0.8,
        "asset_price_index": 0.85,
        "owner_rate": 0.63,
    }
    path_gate = {
        "maximum_market_residual": 1e-6,
        "maximum_mass_residual": 1e-12,
        "maximum_feasibility_projection_mass": 0.0,
        "minimum_distribution_mass": 0.0,
        "maximum_nonfinite_distribution_count": 0,
    }
    provenance = {
        "selected_report_sha256": report.file_sha256(selected_summary_path),
        "selected_case_transition_sha256": report.file_sha256(
            selected_dir / "selected_transition_path.csv"
        ),
        "source_sha256": scientific["source_sha256"],
        "renewal_contract_sha256": report.object_sha256(renewal),
        "target_fingerprint": scientific["target_fingerprint"],
        "code_bundle_sha256": scientific["code_fingerprints"]["bundle_sha256"],
        "scientific_contract_sha256": report.object_sha256(scientific),
        "selection_sha256": selection_sha,
    }
    continuation_summary = {
        "status": "complete_no_policy_post2023_continuation_pair",
        "policy_case": "none",
        "fiscal_change": "none",
        "provenance": provenance,
        "target_set": scientific["target_set"],
        "target_fingerprint": scientific["target_fingerprint"],
        "renewal_accounting_contract": renewal,
        "outside_origin_entry_share": 0.169,
        "outside_origin_entry_status": "external open sensitivity; not identified",
        "outside_share_invariance_audit": {
            "status": "passed",
            "maximum_absolute_identity_residual": 0.0,
        },
        "history_reproduction_audit": {
            "status": "passed",
            "maximum_absolute_gap": 1e-12,
            "tolerance": 5e-10,
        },
        "paired_initial_state_audit": {
            "status": "passed",
            "shared_2007_2023_history_maximum_absolute_gap": 0.0,
            "common_history_periods": 5,
            "post_2023_periods": 40,
            "matched_2023_pre_fertility_distribution_sha256": "1" * 64,
            "matched_2023_birth_vintage_queue_sha256": "2" * 64,
        },
        "paths": {
            "closed_national_benchmark": {
                "closure": "M=0, rho=1",
                "path_gates": path_gate,
                "2023": start,
                "last_date": closed_last,
            },
            "open_cbsa_sensitivity": {
                "closure": "open",
                "path_gates": path_gate,
                "2023": start,
                "last_date": open_last,
            },
        },
        "closed_stationary_endpoint": {
            "status": "no_positive_renewal_root_on_audited_grid",
            "policy_case": "none",
            "usable_closed_root": False,
            "between_steady_states_label_allowed": False,
        },
        "open_stationary_endpoint": {
            "status": "complete",
            "usable_open_endpoint": True,
            "stationary_population_scale": 0.9,
            "price_ratio": 0.84,
            "outside_flow_M": 0.09,
            "entry_households_per_adult": 0.2,
        },
        "between_steady_states_statement": (
            "The national closed finite-horizon path is a benchmark only and is not "
            "described as a transition between steady states."
        ),
    }
    write_json(continuation_dir / "summary.json", continuation_summary)
    write_csv(
        continuation_dir / "paired_continuation_path.csv",
        [
            {
                "scenario": scenario,
                "phase": phase,
                **row,
                "calendar_year": year,
            }
            for scenario, last in (
                ("closed_national_benchmark", closed_last),
                ("open_cbsa_sensitivity", open_last),
            )
            for year, phase, row in (
                [
                    (year, "shared_calibrated_history", start)
                    for year in (2007, 2011, 2015, 2019)
                ]
                + [(2023, "shared_matched_2023_state", start)]
                + [
                    (year, "post_2023_continuation", last)
                    for year in range(2027, 2184, 4)
                ]
            )
        ],
    )
    for name in continuation_artifact_names:
        path = continuation_dir / name
        if path.exists():
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".csv":
            write_csv(path, [{"status": "fixture"}])
        elif path.suffix == ".json":
            write_json(path, {"status": "fixture"})
        else:
            path.write_bytes(f"fixture {name}".encode())
    driver = Path(report.__file__).resolve().with_name(
        "run_e5f_post2023_no_policy_continuations.py"
    )
    manifest = {
        "status": "complete_no_policy_post2023_continuation_manifest",
        "driver": str(driver),
        "driver_sha256": report.file_sha256(driver),
        "input_provenance": provenance,
        "artifacts": {
            name: report.file_sha256(continuation_dir / name)
            for name in continuation_artifact_names
        },
    }
    write_json(continuation_dir / "manifest.json", manifest)

    return argparse.Namespace(
        selected_report=selected_summary_path,
        repeat_task=repeat_paths,
        refinement_plan=plan_path,
        target_provenance_csv=target_provenance_path,
        historical_packet=historical_dir,
        continuation_packet=continuation_dir,
        output_dir=root / "canonical_report",
    )


class CanonicalNoPolicyReportTest(unittest.TestCase):
    def fixture(self) -> tuple[tempfile.TemporaryDirectory[str], argparse.Namespace]:
        temporary = tempfile.TemporaryDirectory()
        return temporary, build_fixture(Path(temporary.name))

    def test_complete_fixture_builds_one_report(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        result = report.build_report(args)
        self.assertEqual(result["status"], "complete_canonical_no_policy_transition_report")
        self.assertEqual(len(read_csv_for_test(args.output_dir / "calibrated_2023_target_fit.csv")), 12)
        self.assertEqual(len(read_csv_for_test(args.output_dir / "estimated_parameters_and_bounds.csv")), 10)
        readme = (args.output_dir / "README.md").read_text().lower()
        self.assertIn("no funded policy", readme)
        self.assertIn("not described as a transition between steady states", readme)
        self.assertTrue((args.output_dir / "report_manifest.sha256").is_file())

    def test_requires_exactly_two_repeats_without_touching_output(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        args.repeat_task = args.repeat_task[:1]
        with self.assertRaisesRegex(report.ContractError, "Exactly two"):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())

    def test_rejects_changed_repeat_moment_without_touching_output(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        path = args.repeat_task[1] / "target_fit_long.csv"
        rows = read_csv_for_test(path)
        rows[0]["model"] = str(float(rows[0]["model"]) + 1e-4)
        write_csv(path, rows)
        with self.assertRaises(report.ContractError):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())

    def test_rejects_incomplete_target_provenance(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        rows = read_csv_for_test(args.target_provenance_csv)[:-1]
        write_csv(args.target_provenance_csv, rows)
        with self.assertRaisesRegex(report.ContractError, "exactly twelve"):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())

    def test_rejects_target_provenance_contract_mismatch(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        rows = read_csv_for_test(args.target_provenance_csv)
        rows[0]["weight"] = str(float(rows[0]["weight"]) + 1.0)
        write_csv(args.target_provenance_csv, rows)
        with self.assertRaisesRegex(report.ContractError, "provenance.*weight"):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())

    def test_rejects_historical_misclassification_without_touching_output(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        provenance_path = args.historical_packet / "classification_and_provenance.json"
        provenance = json.loads(provenance_path.read_text())
        comparison_path = args.historical_packet / "five_date_model_data_comparison.csv"
        rows = read_csv_for_test(comparison_path)
        rows[0]["classification"] = "fitted_population_moment"
        write_csv(comparison_path, rows)
        provenance["outputs"]["comparison_csv"] = artifact_receipt(comparison_path)
        write_json(provenance_path, provenance)
        with self.assertRaisesRegex(report.ContractError, "unknown classifications"):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())

    def test_rejects_continuation_manifest_hash_corruption(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        manifest_path = args.continuation_packet / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["artifacts"]["paired_continuation_path.csv"] = "0" * 64
        write_json(manifest_path, manifest)
        with self.assertRaisesRegex(report.ContractError, "artifact hash differs"):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())

    def test_rejects_missing_hashed_progress_artifact(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        path = args.continuation_packet / "latest_completed_period.json"
        path.unlink()
        with self.assertRaisesRegex(report.ContractError, "artifact is missing"):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())

    def test_rejects_short_continuation(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        summary_path = args.continuation_packet / "summary.json"
        summary = json.loads(summary_path.read_text())
        summary["paired_initial_state_audit"]["post_2023_periods"] = 1
        write_json(summary_path, summary)
        manifest_path = args.continuation_packet / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["artifacts"]["summary.json"] = report.file_sha256(summary_path)
        write_json(manifest_path, manifest)
        with self.assertRaisesRegex(report.ContractError, "exactly 40"):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())

    def test_rejects_fiscal_change(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        summary_path = args.continuation_packet / "summary.json"
        summary = json.loads(summary_path.read_text())
        summary["fiscal_change"] = "tax"
        write_json(summary_path, summary)
        manifest_path = args.continuation_packet / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["artifacts"]["summary.json"] = report.file_sha256(summary_path)
        write_json(manifest_path, manifest)
        with self.assertRaisesRegex(report.ContractError, "policy or a fiscal change"):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())

    def test_rejects_jacobian_corruption(self) -> None:
        temporary, args = self.fixture()
        self.addCleanup(temporary.cleanup)
        path = args.refinement_plan.parent / "weighted_jacobian.csv"
        rows = read_csv_for_test(path)
        rows[0]["parameter_00"] = "2.0"
        write_csv(path, rows)
        with self.assertRaisesRegex(report.ContractError, "singular values differ"):
            report.build_report(args)
        self.assertFalse(args.output_dir.exists())


def read_csv_for_test(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    unittest.main()
