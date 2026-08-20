#!/usr/bin/env python3
"""Pure fixture and corruption tests for the E5F ridge-refinement workflow."""

from __future__ import annotations

import copy
import csv
import json
import math
import tempfile
from pathlib import Path
from typing import Any, Callable

import numpy as np

from build_e5f_transition_ridge_refinement import (
    CENTRAL_STEP,
    ContractError,
    Expectations,
    build_identification_diagnostic,
    build_refinement,
    central_jacobian,
    file_sha256,
    load_coordinate_panel,
    object_sha256,
    transform_unit,
)
from collect_e5f_transition_ridge_refinement import (
    compare_repeat,
    load_and_validate_plan,
    validate_scientific_task,
)


SOURCE_SHA = "a" * 64
CODE_SHA = "b" * 64
CENTER_SHA = "c" * 64
TARGET_SET = "fixture_midpoint_binned_target"
TARGET_FINGERPRINT = "d" * 64
PROFILE = "e5f-income-entry"
PANEL_SEED = 2026081799


DOMAIN = [
    {"name": "beta_annual", "lower": 0.94, "upper": 0.9995, "transform": "discount"},
    {"name": "kappa_fert", "lower": 0.02, "upper": 50.0, "transform": "log"},
    {
        "name": "kappa_fert_continuation",
        "lower": 0.02,
        "upper": 50.0,
        "transform": "log",
    },
    {"name": "chi", "lower": 0.1, "upper": 5.0, "transform": "log"},
    {"name": "H0", "lower": 0.2, "upper": 80.0, "transform": "log"},
    {"name": "theta0", "lower": 0.0, "upper": 8.0, "transform": "softzero"},
    {"name": "theta1", "lower": 0.02, "upper": 16.0, "transform": "log"},
    {
        "name": "hbar_child_rooms",
        "lower": 0.1,
        "upper": 1.8,
        "transform": "log",
    },
    {
        "name": "first_birth_fixed_cost",
        "lower": 0.0,
        "upper": 8.0,
        "transform": "softzero",
    },
    {
        "name": "psi_child_change_2023",
        "lower": -1.5,
        "upper": 0.2,
        "transform": "asinh",
    },
]


def expectations(
    renewal_hash: str | None = None,
    dated_hash: str | None = None,
    *,
    dimensions: int = 10,
    central_step: float = CENTRAL_STEP,
) -> Expectations:
    return Expectations(
        source_sha256=SOURCE_SHA,
        source="/scratch/fixture/source.json",
        target_set=TARGET_SET,
        target_fingerprint=TARGET_FINGERPRINT,
        code_bundle_sha256=CODE_SHA,
        model_profile=PROFILE,
        panel_seed=PANEL_SEED,
        dimensions=dimensions,
        central_step=central_step,
        renewal_contract_sha256=renewal_hash or object_sha256(renewal_contract()),
        dated_contract_sha256=dated_hash or object_sha256(dated_contract()),
    )


def renewal_contract() -> dict[str, Any]:
    return {
        "replacement_fertility": 2.1,
        "effective_birth_to_household_conversion": 1.0 / 2.1,
        "birth_measure": "topcode_adjusted_birth_children",
        "top_bin_mean_children": 3.602,
        "birth_vintage_queue_waiting_slots": 4,
        "birth_to_entry_effect_lag_dates": 5,
        "birth_to_entry_effect_lag_years": 20.0,
        "outside_flow_formula": "M = outside_origin_entry_share * E_old",
        "retention_formula": "rho = (1 - outside_origin_entry_share) * E_old / B_old",
        "status": "fixture replacement accounting",
    }


def dated_contract() -> dict[str, Any]:
    return {
        "cohort_first_birth_timing": "fixed synthetic cohort",
        "first_birth_age_cells": "four-year cells labeled by midpoints 20,...,44",
        "period_fertility": "dated birth flows over exposure",
        "period_fertility_topcode_timing": "parity-three date",
        "period_fertility_status": "diagnostic",
        "dated_maintained_moments": "twelve moments at every date; 2023 targeted",
    }


def population_bridge() -> dict[str, Any]:
    return {
        "status": "externally_normalized_not_estimated",
        "target_indices": {"2007": 1.0, "2023": 1.13},
        "total_households": "fixture Census series",
        "age_distribution": "fixture ACS shares",
        "initial_age_reweight": {
            "mapping": "fixture mapping",
            "model_age_grid": [20, 24, 28],
            "sample": "fixture sample",
            "source": "fixture source",
            "source_path": "/scratch/fixture/pop.csv",
            "source_receipt": "fixture receipt",
            "source_sha256": "e" * 64,
            "source_year": 2007,
            "units": "households",
            "period_years": 4,
            "horizon_periods": 4,
            "horizon_calendar_year": 2023,
            "future_entry_rule": "fixture rule",
            "candidate_dependent_mass": 123.0,
        },
    }


def base_theta() -> dict[str, float]:
    return {
        "alpha_cons": 0.733,
        "beta": 0.98,
        "chi": 1.0,
        "H0": 10.0,
        "theta0": 1.0,
        "theta1": 0.5,
        "hbar_child_rooms": 0.4,
        "hbar_first_child_jump": 0.125,
        "kappa_fert": 2.0,
        "kappa_fert_continuation": 2.5,
        "first_birth_fixed_cost": 2.0,
        "psi_child": 0.0,
        "delta_alpha": 0.0,
        "delta_alpha_jump": 0.0,
        "tenure_choice_kappa": 0.0,
        "theta_n": 0.0,
    }


def fixture_domain(dimensions: int) -> list[dict[str, Any]]:
    if dimensions == 10:
        return list(DOMAIN)
    if dimensions == 11:
        return DOMAIN[:-1] + [
            {
                "name": "hbar_first_child_jump",
                "lower": 0.0,
                "upper": 0.5,
                "transform": "softzero",
            },
            DOMAIN[-1],
        ]
    raise ValueError(dimensions)


def fixture_jacobian(
    rank_deficient: bool = False, *, dimensions: int = 10
) -> np.ndarray:
    matrix = np.zeros((12, dimensions), dtype=float)
    matrix[:dimensions, :] = np.eye(dimensions)
    matrix[dimensions:, :] = np.linspace(0.1, 1.0, dimensions)
    if rank_deficient:
        matrix[:, -1] = matrix[:, 0]
    return matrix


def write_coordinate_fixture(
    root: Path,
    *,
    rank_deficient: bool = False,
    dimensions: int = 10,
    central_step: float = CENTRAL_STEP,
) -> tuple[np.ndarray, np.ndarray]:
    domain_rows = fixture_domain(dimensions)
    task_count = 1 + 2 * dimensions
    anchor = np.full(dimensions, 0.5, dtype=float)
    residual0 = np.linspace(-0.6, 0.5, 12)
    jacobian = fixture_jacobian(rank_deficient, dimensions=dimensions)
    target = np.linspace(0.2, 1.3, 12)
    weights = np.linspace(1.0, 3.2, 12)
    renewal = renewal_contract()
    dated = dated_contract()
    for task_id in range(1, task_count + 1):
        if task_id == 1:
            unit = anchor.copy()
            design = "anchor"
        else:
            offset = task_id - 2
            dimension = offset // 2
            sign = -1.0 if offset % 2 == 0 else 1.0
            unit = anchor.copy()
            unit[dimension] += sign * central_step
            design = f"coordinate_{dimension}_{'minus' if sign < 0 else 'plus'}"
        residual = residual0 + jacobian @ (unit - anchor)
        model = target + residual / np.sqrt(weights)
        loss = float(residual @ residual)
        task_dir = root / f"task_{task_id:03d}"
        task_dir.mkdir(parents=True)
        target_rows: list[dict[str, Any]] = []
        for index in range(12):
            gap = float(model[index] - target[index])
            standardized = math.sqrt(float(weights[index])) * gap
            target_rows.append(
                {
                    "candidate": f"task_{task_id:03d}",
                    "moment": f"moment_{index:02d}",
                    "target": float(target[index]),
                    "model": float(model[index]),
                    "gap": gap,
                    "weight": float(weights[index]),
                    "loss_contribution": standardized**2,
                    "standardized_gap": standardized,
                }
            )
        with (task_dir / "target_fit_long.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(target_rows[0]))
            writer.writeheader()
            writer.writerows(target_rows)
        theta = base_theta()
        terminal_change = math.nan
        for coordinate, domain in zip(unit, domain_rows, strict=True):
            value = transform_unit(
                float(coordinate),
                float(domain["lower"]),
                float(domain["upper"]),
                str(domain["transform"]),
            )
            if domain["name"] == "beta_annual":
                theta["beta"] = value**4
            elif domain["name"] == "psi_child_change_2023":
                terminal_change = value
            else:
                theta[str(domain["name"])] = value
        assert math.isfinite(terminal_change)
        best = {
            "candidate": f"task_{task_id:03d}",
            "theta": theta,
            "new_psi_child": 0.1 + terminal_change,
            "transition_loss": loss,
            "max_market_residual": 1.0e-5,
            "max_mass_residual": 1.0e-14,
            "population_target_gap": 1.0e-15,
            "terminal_synthetic_childless_consistency_gap": 1.0e-13,
        }
        summary = {
            "status": "complete_transition_calibration_panel_task",
            "source": "/scratch/fixture/source.json",
            "source_sha256": SOURCE_SHA,
            "source_metadata": {"fixture": True},
            "code_fingerprints": {
                "bundle_sha256": CODE_SHA,
                "files": {"driver": "f" * 64},
            },
            "model_profile": {
                "name": PROFILE,
                "profile_id": "fixture_income_profile",
                "income_state_count": 15,
                "permanent_income_log_variance": 0.393053,
                "first_birth_fixed_cost_semantics": "one-time cost",
            },
            "income_profile_gates": {
                "permanent_income_levels_enabled": True,
                "income_state_count": 15,
                "stationary_weight_max_abs_gap": 0.0,
            },
            "target_set": TARGET_SET,
            "target_fingerprint": TARGET_FINGERPRINT,
            "target_count": 12,
            "stationary_free_parameter_count": dimensions,
            "transition_free_parameter_count": dimensions,
            "panel_design": {
                "task_id": task_id,
                "panel_size": task_count,
                "panel_seed": PANEL_SEED,
                "design": design,
                "panel_design": "coordinate",
                "local_radius": central_step,
                "center": "/scratch/fixture/center.json",
                "center_sha256": CENTER_SHA,
                "terminal_preference_coordinate": "psi_child_change_2023",
                "unit_vector": unit.tolist(),
                "domain": domain_rows,
                "old_psi_child": 0.1,
                "new_psi_child": 0.1 + terminal_change,
                "psi_child_change_2023": terminal_change,
            },
            "policy_case": "none",
            "post_2023_periods": 0,
            "old_psi_child": 0.1,
            "old_model_completed_fertility": 2.1,
            "old_completed_fertility_reference": 2.1,
            "outside_origin_entry_share": 0.169,
            "renewal_accounting_contract": renewal,
            "renewal_accounting_old_state": {
                "old_entry_flow_E": 1.0,
                "old_queue_mature_flow_B": 1.0,
                "old_queue_B_over_E": 1.0,
                "old_renewal_residual": 0.0,
                "outside_flow_M": 0.169,
                "outside_origin_entry_share": 0.169,
                "retention_rho": 0.831,
            },
            "target_measurements": {"timing": "midpoint-binned fixture"},
            "dated_measurement_contract": dated,
            "population_bridge": population_bridge(),
            "population_validation_status": "fixture bridge",
            "stationary_measurement_max_abs_gap": 1.0e-10,
            "best_candidate": best,
        }
        (task_dir / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return jacobian, residual0


def assert_rejected(function: Callable[[], Any], contains: str | None = None) -> None:
    try:
        function()
    except ContractError as error:
        if contains is not None:
            assert contains.lower() in str(error).lower(), (contains, str(error))
        return
    raise AssertionError("Expected ContractError")


def write_scientific_validation_fixture(
    task_dir: Path,
    coordinate_dir: Path,
    plan: dict[str, Any],
    plan_path: Path,
    plan_sha: str,
    candidate_record: dict[str, Any],
) -> None:
    task_dir.mkdir(parents=True)
    summary = json.loads(
        (coordinate_dir / "task_001" / "summary.json").read_text(encoding="utf-8")
    )
    payload = json.loads(
        (plan_path.parent / candidate_record["candidate_file"]).read_text(encoding="utf-8")
    )
    terminal_delta = float(payload["best_candidate"]["new_psi_child"]) - float(
        payload["best_candidate"]["old_psi_child"]
    )
    summary["panel_design"].update(
        {
            "task_id": 1,
            "panel_size": 1,
            "panel_seed": plan["coordinate_panel_contract"]["panel_seed"],
            "design": "anchor",
            "panel_design": "mixed",
            "local_radius": 0.02,
            "center": str(plan_path.parent / candidate_record["candidate_file"]),
            "center_sha256": candidate_record["candidate_sha256"],
            "unit_vector": candidate_record["unit_vector"],
            "domain": plan["domain"],
            "old_psi_child": 0.1,
            "new_psi_child": 0.1 + terminal_delta,
            "psi_child_change_2023": terminal_delta,
        }
    )
    summary["best_candidate"]["theta"] = payload["best_candidate"]["theta"]
    summary["best_candidate"]["new_psi_child"] = 0.1 + terminal_delta
    (task_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (task_dir / "target_fit_long.csv").write_bytes(
        (coordinate_dir / "task_001" / "target_fit_long.csv").read_bytes()
    )
    wrapper = {
        "schema": "e5f_transition_ridge_task_contract_v1",
        "mode": "validation",
        "outer_task_id": int(candidate_record["candidate_id"]),
        "replicate_id": None,
        "candidate_id": int(candidate_record["candidate_id"]),
        "candidate_sha256": candidate_record["candidate_sha256"],
        "plan_sha256": plan_sha,
        "scientific_contract_sha256": plan["scientific_contract_sha256"],
        "source_sha256": plan["scientific_contract"]["source_sha256"],
        "target_set": plan["scientific_contract"]["target_set"],
        "target_fingerprint": plan["scientific_contract"]["target_fingerprint"],
        "code_bundle_sha256": plan["scientific_contract"]["code_fingerprints"][
            "bundle_sha256"
        ],
        "renewal_contract_sha256": plan["renewal_contract_sha256"],
        "dated_contract_sha256": plan["dated_contract_sha256"],
        "selection_sha256": None,
        "execution_identity": f"fixture:{candidate_record['candidate_id']}:validation",
    }
    (task_dir / "refinement_task_contract.json").write_text(
        json.dumps(wrapper, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def test_valid_fixture_builds_exact_hashed_plan() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        results = root / "coordinate"
        expected_jacobian, _ = write_coordinate_fixture(results)
        renewal_sha = object_sha256(renewal_contract())
        dated_sha = object_sha256(dated_contract())
        panel = load_coordinate_panel(results, expectations(renewal_sha, dated_sha))
        recovered, sides, frozen = central_jacobian(panel)
        assert np.allclose(recovered, expected_jacobian, rtol=0.0, atol=1.0e-11)
        assert len(sides) == 10
        assert not np.any(frozen)
        built = build_refinement(results, root / "refinement", expectations(renewal_sha, dated_sha))
        plan = built["plan"]
        assert plan["refinement_builder_sha256"] == file_sha256(
            Path(__file__).resolve().parent / "build_e5f_transition_ridge_refinement.py"
        )
        assert plan["coordinate_panel_contract"]["task_count"] == 21
        assert plan["coordinate_panel_contract"]["central_step"] == 0.02
        assert plan["jacobian_diagnostics"]["full"]["relative_ranks"]["relative_0.001"] == 10
        assert len(plan["candidates"]) == 7
        assert sum(row["kind"] == "ridge_trust" for row in plan["candidates"]) == 6
        assert sum(row["kind"] == "observed_incumbent" for row in plan["candidates"]) == 1
        plan_path = Path(built["plan_path"])
        assert file_sha256(plan_path) == built["plan_sha256"]
        loaded, loaded_sha = load_and_validate_plan(plan_path, built["plan_sha256"])
        assert loaded_sha == built["plan_sha256"]
        assert loaded["target_contract_sha256"] == plan["target_contract_sha256"]


def test_eleven_parameter_fixture_builds_and_passes_collector_plan_gate() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        results = root / "coordinate"
        expected_jacobian, _ = write_coordinate_fixture(results, dimensions=11)
        expected = expectations(dimensions=11)
        panel = load_coordinate_panel(results, expected)
        recovered, sides, frozen = central_jacobian(panel)
        assert recovered.shape == (12, 11)
        assert np.allclose(recovered, expected_jacobian, rtol=0.0, atol=1.0e-11)
        assert len(sides) == 11
        assert not np.any(frozen)
        built = build_refinement(results, root / "refinement", expected)
        plan = built["plan"]
        assert plan["coordinate_panel_contract"]["task_count"] == 23
        assert plan["coordinate_panel_contract"]["dimensions"] == 11
        assert plan["jacobian_diagnostics"]["full"]["relative_ranks"][
            "relative_0.001"
        ] == 11
        loaded, _ = load_and_validate_plan(
            Path(built["plan_path"]), built["plan_sha256"]
        )
        assert loaded["coordinate_panel_contract"]["dimensions"] == 11


def test_one_percent_coordinate_radius_builds_and_passes_collector_plan_gate() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        results = root / "coordinate"
        write_coordinate_fixture(results, dimensions=11, central_step=0.01)
        expected = expectations(dimensions=11, central_step=0.01)
        built = build_refinement(results, root / "refinement", expected)
        plan = built["plan"]
        assert plan["coordinate_panel_contract"]["central_step"] == 0.01
        loaded, _ = load_and_validate_plan(
            Path(built["plan_path"]), built["plan_sha256"]
        )
        assert loaded["coordinate_panel_contract"]["central_step"] == 0.01


def test_missing_coordinate_task_is_rejected() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        results = Path(temporary) / "coordinate"
        write_coordinate_fixture(results)
        shutil_path = results / "task_021"
        for path in shutil_path.iterdir():
            path.unlink()
        shutil_path.rmdir()
        assert_rejected(lambda: load_coordinate_panel(results, expectations()), "exactly")


def test_corrupted_central_pair_is_rejected() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        results = Path(temporary) / "coordinate"
        write_coordinate_fixture(results)
        summary_path = results / "task_002" / "summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["panel_design"]["unit_vector"][1] += 0.001
        summary_path.write_text(json.dumps(summary, sort_keys=True), encoding="utf-8")
        assert_rejected(lambda: load_coordinate_panel(results, expectations()), "unit")


def test_mixed_dated_contract_is_rejected() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        results = Path(temporary) / "coordinate"
        write_coordinate_fixture(results)
        summary_path = results / "task_013" / "summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["dated_measurement_contract"]["first_birth_age_cells"] = "wrong endpoint bins"
        summary_path.write_text(json.dumps(summary, sort_keys=True), encoding="utf-8")
        assert_rejected(
            lambda: load_coordinate_panel(results, expectations()),
            "dated contract hash",
        )


def test_target_weight_corruption_is_rejected() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        results = Path(temporary) / "coordinate"
        write_coordinate_fixture(results)
        fit_path = results / "task_008" / "target_fit_long.csv"
        rows = list(csv.DictReader(fit_path.open(newline="", encoding="utf-8")))
        rows[0]["weight"] = str(float(rows[0]["weight"]) * 2.0)
        with fit_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        assert_rejected(lambda: load_coordinate_panel(results, expectations()), "standardized")


def test_parameter_to_coordinate_corruption_is_rejected() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        results = Path(temporary) / "coordinate"
        write_coordinate_fixture(results)
        summary_path = results / "task_009" / "summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["best_candidate"]["theta"]["theta1"] *= 1.01
        summary_path.write_text(json.dumps(summary, sort_keys=True), encoding="utf-8")
        assert_rejected(
            lambda: load_coordinate_panel(results, expectations()),
            "parameter-to-unit reconstruction",
        )


def test_preference_metadata_corruption_is_rejected() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        results = Path(temporary) / "coordinate"
        write_coordinate_fixture(results)
        summary_path = results / "task_011" / "summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["best_candidate"]["new_psi_child"] += 0.001
        summary_path.write_text(json.dumps(summary, sort_keys=True), encoding="utf-8")
        assert_rejected(
            lambda: load_coordinate_panel(results, expectations()),
            "candidate new psi versus panel new psi",
        )


def test_target_candidate_label_corruption_is_rejected() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        results = Path(temporary) / "coordinate"
        write_coordinate_fixture(results)
        fit_path = results / "task_008" / "target_fit_long.csv"
        rows = list(csv.DictReader(fit_path.open(newline="", encoding="utf-8")))
        rows[3]["candidate"] = "task_999"
        with fit_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        assert_rejected(
            lambda: load_coordinate_panel(results, expectations()),
            "target row candidate",
        )


def test_rank_deficiency_blocks_candidate_generation() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        results = root / "coordinate"
        write_coordinate_fixture(results, rank_deficient=True)
        assert_rejected(
            lambda: build_refinement(results, root / "refinement", expectations()),
            "rank deficient",
        )
        assert not (root / "refinement" / "candidate_plan.json").exists()


def test_rank_deficiency_writes_diagnostic_without_candidates() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        results = root / "coordinate"
        write_coordinate_fixture(results, rank_deficient=True)
        outdir = root / "diagnostic"
        summary = build_identification_diagnostic(results, outdir, expectations())
        assert summary["status"] == "rank_gate_failed_no_candidates_created"
        assert not summary["identification_gate_passed"]
        assert not summary["candidate_centers_created"]
        assert (outdir / "diagnostic_summary.json").is_file()
        assert (outdir / "weighted_jacobian.csv").is_file()
        assert not (outdir / "candidate_plan.json").exists()


def test_nonpositive_side_consistency_freezes_column() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        results = Path(temporary) / "coordinate"
        write_coordinate_fixture(results)
        panel = load_coordinate_panel(results, expectations())
        anchor = panel.tasks[0].residual.copy()
        direction = np.zeros_like(anchor)
        direction[0] = 1.0
        panel.tasks[1].residual = anchor + CENTRAL_STEP * direction
        panel.tasks[2].residual = anchor + CENTRAL_STEP * direction
        _, sides, frozen = central_jacobian(panel)
        assert bool(frozen[0])
        assert sides[0]["forward_backward_dot"] < 0.0
        assert sides[0]["side_consistency_cosine"] < 0.0


def test_candidate_hash_corruption_is_rejected() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        results = root / "coordinate"
        write_coordinate_fixture(results)
        built = build_refinement(results, root / "refinement", expectations())
        candidate_path = root / "refinement" / "candidates" / "candidate_004.json"
        candidate_path.write_text(candidate_path.read_text(encoding="utf-8") + " ", encoding="utf-8")
        assert_rejected(
            lambda: load_and_validate_plan(Path(built["plan_path"]), built["plan_sha256"]),
            "hash mismatch",
        )


def test_scientific_validation_contract_and_wrapper_corruption() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        coordinate = root / "coordinate"
        write_coordinate_fixture(coordinate)
        built = build_refinement(coordinate, root / "refinement", expectations())
        plan = built["plan"]
        plan_path = Path(built["plan_path"])
        record = plan["candidates"][3]
        task_dir = root / "validation" / "task_004"
        write_scientific_validation_fixture(
            task_dir, coordinate, plan, plan_path, built["plan_sha256"], record
        )
        result = validate_scientific_task(
            task_dir,
            plan,
            plan_path,
            built["plan_sha256"],
            record,
            mode="validation",
            outer_task_id=4,
            replicate_id=None,
            selection_sha256=None,
        )
        assert result["candidate_id"] == 4
        wrapper_path = task_dir / "refinement_task_contract.json"
        wrapper = json.loads(wrapper_path.read_text(encoding="utf-8"))
        wrapper["code_bundle_sha256"] = "0" * 64
        wrapper_path.write_text(json.dumps(wrapper, sort_keys=True), encoding="utf-8")
        assert_rejected(
            lambda: validate_scientific_task(
                task_dir,
                plan,
                plan_path,
                built["plan_sha256"],
                record,
                mode="validation",
                outer_task_id=4,
                replicate_id=None,
                selection_sha256=None,
            ),
            "wrapper code_bundle_sha256",
        )


def test_repeat_moment_corruption_is_rejected() -> None:
    reference = {
        "candidate_id": 2,
        "candidate_sha256": "a" * 64,
        "loss": 1.0,
        "candidate": {"theta": {"beta": 0.98}},
        "target_rows": [
            {
                "moment": "tfr",
                "target": 1.9,
                "weight": 2.0,
                "model": 1.8,
                "gap": -0.1,
                "standardized_gap": -math.sqrt(2.0) * 0.1,
                "loss_contribution": 0.02,
            }
        ],
        "summary": {
            "old_psi_child": 0.1,
            "old_model_completed_fertility": 2.1,
            "panel_design": {"old_psi_child": 0.1, "new_psi_child": -0.2},
            "renewal_accounting_old_state": {
                "old_entry_flow_E": 1.0,
                "old_queue_mature_flow_B": 1.0,
            },
            "stationary_measurement_max_abs_gap": 0.0,
            "best_candidate": {
                "max_market_residual": 1.0e-5,
                "max_mass_residual": 1.0e-14,
                "population_target_gap": 0.0,
                "terminal_synthetic_childless_consistency_gap": 0.0,
            },
        },
    }
    repeat = copy.deepcopy(reference)
    repeat["replicate_id"] = 1
    repeat["target_rows"][0]["model"] += 2.0e-10
    assert_rejected(lambda: compare_repeat(reference, repeat), "model differs")


def test_launcher_resolves_the_slurm_submission_checkout() -> None:
    launcher = (
        Path(__file__).resolve().parents[3]
        / "code/cluster/submit_e5f_transition_ridge_refinement.sh"
    )
    text = launcher.read_text(encoding="utf-8")
    assert "SLURM_SUBMIT_DIR" in text
    assert 'CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-' in text
    assert "E5F_RIDGE_ESTIMATE_FIRST_CHILD_ROOM_JUMP" in text
    assert "--estimate-first-child-room-jump" in text
    assert "expected_dimensions = 11" in text


def main() -> None:
    tests = [
        test_valid_fixture_builds_exact_hashed_plan,
        test_eleven_parameter_fixture_builds_and_passes_collector_plan_gate,
        test_one_percent_coordinate_radius_builds_and_passes_collector_plan_gate,
        test_missing_coordinate_task_is_rejected,
        test_corrupted_central_pair_is_rejected,
        test_mixed_dated_contract_is_rejected,
        test_target_weight_corruption_is_rejected,
        test_parameter_to_coordinate_corruption_is_rejected,
        test_preference_metadata_corruption_is_rejected,
        test_target_candidate_label_corruption_is_rejected,
        test_rank_deficiency_blocks_candidate_generation,
        test_rank_deficiency_writes_diagnostic_without_candidates,
        test_nonpositive_side_consistency_freezes_column,
        test_candidate_hash_corruption_is_rejected,
        test_scientific_validation_contract_and_wrapper_corruption,
        test_repeat_moment_corruption_is_rejected,
        test_launcher_resolves_the_slurm_submission_checkout,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"PASS all {len(tests)} ridge-refinement tests")


if __name__ == "__main__":
    main()
