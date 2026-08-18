#!/usr/bin/env python3
"""Profile the fertility-preference change at fixed fertility-logit slopes.

For one predeclared common multiplier on the first- and continuation-birth
logit scales, this diagnostic repeatedly calls the certified dated-transition
driver.  It re-normalizes the 2007 preference intercept to replacement in every
trial and solves the 2007--2023 preference change so the model's 2023 completed
fertility matches the active target contract.  Only after that scalar profile
passes does it audit the closed stationary renewal schedule on the exact
25-price grid.

All other parameters remain at the final ridge estimate.  The exercise is a
conditional profile and root diagnostic, not a re-estimated SMM calibration.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

import run_e5f_final_closed_root_slope_frontier as fixed


SCHEMA = "e5f_target_preserving_closed_root_factor_v1"
REPORT_SCHEMA = "e5f_target_preserving_closed_root_frontier_v1"
FACTORS = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
INITIAL_DELTAS = (-0.35, -0.20, -0.50, 0.0, -0.80, 0.15, -1.20)
DELTA_LOWER = -1.50
DELTA_UPPER = 0.20
TARGET_TOLERANCE = 5.0e-4
MONOTONICITY_TOLERANCE = 5.0e-4
MAX_PROFILE_ITERATIONS = 4


class ProfileError(RuntimeError):
    """Raised when the conditional profile violates its declared contract."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ridge-report", type=Path, required=True)
    parser.add_argument("--selected-task-summary", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--factor-index", type=int)
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--results-root", type=Path)
    parser.add_argument(
        "--run-stage", choices=("smoke", "factor", "collect"), required=True
    )
    parser.add_argument("--expected-ridge-report-sha256", required=True)
    parser.add_argument("--expected-task-summary-sha256", required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-target-set", required=True)
    parser.add_argument("--expected-target-fingerprint", required=True)
    parser.add_argument("--expected-code-bundle-sha256", required=True)
    parser.add_argument("--expected-scientific-contract-sha256", required=True)
    parser.add_argument("--expected-selection-sha256", required=True)
    parser.add_argument("--expected-selected-candidate-sha256", required=True)
    parser.add_argument("--expected-selected-candidate-id", type=int, required=True)
    parser.add_argument("--expected-renewal-contract-sha256", required=True)
    parser.add_argument("--expected-dated-contract-sha256", required=True)
    parser.add_argument("--expected-frontier-driver-sha256", required=True)
    parser.add_argument("--expected-continuation-driver-sha256", required=True)
    parser.add_argument("--expected-launcher-sha256", required=True)
    parser.add_argument("--expected-diagnostic-bundle-sha256", required=True)
    parser.add_argument("--expected-profile-driver-sha256", required=True)
    parser.add_argument("--expected-profile-launcher-sha256", required=True)
    parser.add_argument("--expected-profile-bundle-sha256", required=True)
    parser.add_argument("--smoke-receipt", type=Path)
    parser.add_argument("--expected-smoke-receipt-sha256")
    return parser.parse_args()


def factor_for_index(index: int) -> float:
    if not 1 <= int(index) <= len(FACTORS):
        raise ProfileError(f"factor index must lie in [1,{len(FACTORS)}]")
    return float(FACTORS[int(index) - 1])


def profile_code_contract() -> dict[str, str]:
    driver = Path(__file__).resolve()
    launcher = (
        ROOT / "code/cluster/submit_e5f_target_preserving_closed_root_frontier.sh"
    ).resolve()
    hashes = {
        "profile_driver_sha256": fixed.file_sha256(driver),
        "profile_launcher_sha256": fixed.file_sha256(launcher),
    }
    hashes["profile_bundle_sha256"] = fixed.canonical_json_sha256(hashes)
    return hashes


def validate_profile_code(args: argparse.Namespace) -> dict[str, str]:
    hashes = profile_code_contract()
    fixed.require_equal(
        hashes["profile_driver_sha256"],
        args.expected_profile_driver_sha256,
        "profile driver hash",
    )
    fixed.require_equal(
        hashes["profile_launcher_sha256"],
        args.expected_profile_launcher_sha256,
        "profile launcher hash",
    )
    fixed.require_equal(
        hashes["profile_bundle_sha256"],
        args.expected_profile_bundle_sha256,
        "profile bundle hash",
    )
    return hashes


def base_contracts(args: argparse.Namespace) -> dict[str, Any]:
    proxy = argparse.Namespace(**vars(args))
    proxy.run_stage = "smoke"
    proxy.smoke_receipt = None
    proxy.expected_smoke_receipt_sha256 = None
    contracts = fixed.validate_inputs(proxy, require_live_code=True)
    from intergen_eqscale_seq_optimized.e5_profile import e5_target_system

    target_system = e5_target_system()
    targets = target_system.targets_dict()
    weights = target_system.weights_dict()
    if "tfr" not in targets or "tfr" not in weights:
        raise ProfileError("active target contract omits completed fertility")
    contracts["profile_target_moment"] = "tfr"
    contracts["profile_target_value"] = float(targets["tfr"])
    contracts["profile_target_weight"] = float(weights["tfr"])
    return contracts


def validate_smoke_receipt(
    args: argparse.Namespace,
    expected_provenance: dict[str, str],
    target_value: float,
) -> dict[str, str]:
    if args.run_stage == "smoke":
        if args.smoke_receipt is not None or args.expected_smoke_receipt_sha256:
            raise ProfileError("smoke stage forbids a recursive smoke receipt")
        return {}
    if args.smoke_receipt is None or not args.expected_smoke_receipt_sha256:
        raise ProfileError(f"{args.run_stage} requires an exact-loop smoke receipt")
    receipt = args.smoke_receipt.resolve()
    receipt_sha = fixed.require_file_hash(
        receipt, args.expected_smoke_receipt_sha256, "profile smoke receipt"
    )
    payload = fixed.read_json(receipt)
    fixed.require_equal(payload.get("schema"), SCHEMA, "profile smoke schema")
    fixed.require_equal(
        payload.get("status"), "complete_target_preserving_profile", "profile smoke status"
    )
    fixed.require_equal(payload.get("run_stage"), "smoke", "profile smoke stage")
    fixed.require_equal(payload.get("factor_index"), 3, "profile smoke factor index")
    fixed.require_close(payload.get("kappa_factor"), 0.75, "profile smoke factor", 1e-14)
    if not bool(payload.get("target_profile_success", False)):
        raise ProfileError("exact-loop smoke did not solve the completed-fertility profile")
    fixed.require_close(
        payload.get("profiled_terminal_completed_fertility"),
        target_value,
        "profile smoke completed fertility",
        TARGET_TOLERANCE,
    )
    provenance = dict(payload.get("input_provenance") or {})
    for name, expected in expected_provenance.items():
        fixed.require_equal(provenance.get(name), expected, f"smoke provenance {name}")
    return {"profile_smoke_receipt_sha256": receipt_sha}


def build_center(
    selected: dict[str, Any], factor: float, preference_change: float
) -> tuple[dict[str, Any], dict[str, float]]:
    _, scaled_theta, _ = fixed.build_factor_center(selected, factor)
    center = {
        "schema": "e5f_target_preserving_closed_root_center_v1",
        "kappa_factor": factor,
        "profiled_preference_change": preference_change,
        "best_candidate": {
            "theta": scaled_theta,
            "old_psi_child": 0.0,
            "new_psi_child": preference_change,
        },
    }
    return center, scaled_theta


def run_trial(
    args: argparse.Namespace,
    contracts: dict[str, Any],
    heartbeat: fixed.Heartbeat,
    outdir: Path,
    factor: float,
    preference_change: float,
    trial_id: int,
) -> dict[str, Any]:
    trial_dir = outdir / f"trial_{trial_id:02d}"
    if trial_dir.exists():
        raise ProfileError(f"refusing to overwrite trial directory: {trial_dir}")
    trial_dir.mkdir(parents=True)
    center, scaled_theta = build_center(
        contracts["selected_candidate"], factor, preference_change
    )
    center_path = trial_dir / "center.json"
    fixed.atomic_write_json(center_path, center)
    transition_dir = trial_dir / "transition"
    heartbeat.set_phase(f"profile_trial_{trial_id:02d}_delta_{preference_change:.8f}")
    command = [
        sys.executable,
        str(TOOLS_ROOT / "run_e5f_transition_calibration.py"),
        "--source",
        str(contracts["source_path"]),
        "--expected-source-sha256",
        str(args.expected_source_sha256),
        "--expected-target-set",
        str(args.expected_target_set),
        "--expected-target-fingerprint",
        str(args.expected_target_fingerprint),
        "--expected-code-bundle-sha256",
        str(args.expected_code_bundle_sha256),
        "--outdir",
        str(transition_dir),
        "--panel-task-id",
        "1",
        "--panel-size",
        "1",
        "--panel-seed",
        "2026081703",
        "--panel-local-radius",
        "0.02",
        "--panel-design",
        "mixed",
        "--panel-center-json",
        str(center_path),
        "--old-psi-child",
        str(float(contracts["selected_candidate"]["old_psi_child"])),
        "--replacement-fertility",
        "2.1",
        "--old-completed-fertility-target",
        "2.1",
        "--post-2023-periods",
        "0",
        "--policy-case",
        "none",
        "--model-profile",
        fixed.MODEL_PROFILE,
        "--estimate-first-child-room-jump",
        "--outside-origin-entry-share",
        str(float(contracts["scientific_contract"]["outside_origin_entry_share"])),
    ]
    elapsed = fixed.run_subprocess_with_heartbeat(command, trial_dir, heartbeat)
    summary, fit = fixed.validate_transition_output(
        transition_dir, contracts, scaled_theta, preference_change, factor
    )
    best = dict(summary["best_candidate"])
    completed = float(best["terminal_tfr"])
    target_value = float(contracts["profile_target_value"])
    row = {
        "trial_id": trial_id,
        "preference_change": float(preference_change),
        "old_psi_child": float(summary["old_psi_child"]),
        "new_psi_child": float(best["new_psi_child"]),
        "terminal_completed_fertility": completed,
        "target_gap": completed - target_value,
        "transition_loss": float(best["transition_loss"]),
        "elapsed_seconds": elapsed,
        "summary_path": str((transition_dir / "summary.json").resolve()),
        "summary_sha256": fixed.file_sha256(transition_dir / "summary.json"),
    }
    fixed.atomic_write_json(trial_dir / "trial_summary.json", row)
    fixed.atomic_write_csv(trial_dir / "target_fit.csv", [dict(item) for item in fit])
    return {"row": row, "summary": summary, "fit": fit, "dir": trial_dir}


def monotone_trials(trials: list[dict[str, Any]]) -> bool:
    ordered = sorted(trials, key=lambda item: float(item["row"]["preference_change"]))
    values = np.asarray(
        [float(item["row"]["terminal_completed_fertility"]) for item in ordered]
    )
    return bool(np.all(np.diff(values) >= -MONOTONICITY_TOLERANCE))


def target_bracket(
    trials: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any]] | None:
    ordered = sorted(trials, key=lambda item: float(item["row"]["preference_change"]))
    for lower, upper in zip(ordered[:-1], ordered[1:]):
        lower_gap = float(lower["row"]["target_gap"])
        upper_gap = float(upper["row"]["target_gap"])
        if lower_gap <= 0.0 <= upper_gap:
            return lower, upper
    return None


def safeguarded_secant(
    lower: dict[str, Any], upper: dict[str, Any]
) -> float:
    x0 = float(lower["row"]["preference_change"])
    x1 = float(upper["row"]["preference_change"])
    y0 = float(lower["row"]["target_gap"])
    y1 = float(upper["row"]["target_gap"])
    if not (x0 < x1 and y0 <= 0.0 <= y1):
        raise ProfileError("invalid completed-fertility bracket")
    if abs(y1 - y0) <= 1e-12:
        candidate = 0.5 * (x0 + x1)
    else:
        candidate = x0 - y0 * (x1 - x0) / (y1 - y0)
    margin = min(1e-5, 0.05 * (x1 - x0))
    return float(np.clip(candidate, x0 + margin, x1 - margin))


def choose_best(trials: list[dict[str, Any]]) -> dict[str, Any]:
    return min(trials, key=lambda item: abs(float(item["row"]["target_gap"])))


def run_factor(args: argparse.Namespace) -> None:
    if args.factor_index is None or args.collect:
        raise ProfileError("factor execution requires --factor-index and forbids --collect")
    factor_index = int(args.factor_index)
    factor = factor_for_index(factor_index)
    if args.run_stage not in {"smoke", "factor"}:
        raise ProfileError("factor execution requires smoke or factor run stage")
    if args.run_stage == "smoke" and factor_index != 3:
        raise ProfileError("exact-loop smoke is factor index 3 (0.75)")
    outdir = args.output_dir.resolve()
    if outdir.exists():
        raise ProfileError(f"refusing to overwrite output directory: {outdir}")
    outdir.mkdir(parents=True)
    heartbeat = fixed.Heartbeat(outdir, factor_index, factor)
    heartbeat.start()
    started = time.perf_counter()
    try:
        code_hashes = validate_profile_code(args)
        contracts = base_contracts(args)
        expected_provenance = {**contracts["hashes"], **code_hashes}
        target_value = float(contracts["profile_target_value"])
        smoke_hashes = validate_smoke_receipt(
            args, expected_provenance, target_value
        )
        trials: list[dict[str, Any]] = []
        for preference_change in INITIAL_DELTAS:
            if any(
                abs(float(item["row"]["preference_change"]) - preference_change)
                <= 1e-12
                for item in trials
            ):
                continue
            trial = run_trial(
                args,
                contracts,
                heartbeat,
                outdir,
                factor,
                float(preference_change),
                len(trials) + 1,
            )
            trials.append(trial)
            fixed.atomic_write_csv(
                outdir / "profile_trials.csv", [item["row"] for item in trials]
            )
            if abs(float(trial["row"]["target_gap"])) <= TARGET_TOLERANCE:
                break
            if target_bracket(trials) is not None and len(trials) >= 2:
                break

        if not monotone_trials(trials):
            raise ProfileError("completed fertility is not monotone in the preference change")
        bracket = target_bracket(trials)
        for _ in range(MAX_PROFILE_ITERATIONS):
            best = choose_best(trials)
            if abs(float(best["row"]["target_gap"])) <= TARGET_TOLERANCE:
                break
            if bracket is None:
                break
            next_delta = safeguarded_secant(*bracket)
            if any(
                abs(float(item["row"]["preference_change"]) - next_delta) <= 1e-8
                for item in trials
            ):
                next_delta = 0.5 * (
                    float(bracket[0]["row"]["preference_change"])
                    + float(bracket[1]["row"]["preference_change"])
                )
            trial = run_trial(
                args,
                contracts,
                heartbeat,
                outdir,
                factor,
                next_delta,
                len(trials) + 1,
            )
            trials.append(trial)
            fixed.atomic_write_csv(
                outdir / "profile_trials.csv", [item["row"] for item in trials]
            )
            if not monotone_trials(trials):
                raise ProfileError(
                    "completed fertility lost monotonicity during safeguarded profiling"
                )
            bracket = target_bracket(trials)

        best = choose_best(trials)
        profile_success = abs(float(best["row"]["target_gap"])) <= TARGET_TOLERANCE
        endpoint: dict[str, Any] | None = None
        schedule: list[dict[str, Any]] = []
        endpoint_seconds = 0.0
        if profile_success:
            heartbeat.set_phase("profiled_closed_stationary_schedule")
            endpoint_dir = outdir / "closed_endpoint"
            endpoint_dir.mkdir()
            endpoint, schedule, endpoint_seconds = fixed.solve_closed_endpoint(
                best["summary"], contracts, endpoint_dir
            )
            shutil.copy2(
                endpoint_dir / "closed_stationary_schedule.csv",
                outdir / "closed_stationary_schedule.csv",
            )
            shutil.copy2(
                endpoint_dir / "closed_stationary_endpoint.json",
                outdir / "closed_stationary_endpoint.json",
            )
        else:
            fixed.atomic_write_csv(
                outdir / "closed_stationary_schedule.csv",
                [{"status": "not_run_profile_target_unresolved"}],
            )
            fixed.atomic_write_json(
                outdir / "closed_stationary_endpoint.json",
                {"status": "not_run_profile_target_unresolved"},
            )

        fit_rows = [dict(row) for row in best["fit"]]
        fixed.atomic_write_csv(outdir / "profiled_target_fit.csv", fit_rows)
        parameter_rows = fixed.build_parameter_table(
            best["dir"] / "transition",
            dict(contracts["selected_candidate"]["theta"]),
            dict(best["summary"]["best_candidate"]["theta"]),
            factor,
        )
        for row in parameter_rows:
            if row["parameter"] in {"psi_child_2023", "psi_child_change_2023"}:
                row["diagnostic_status"] = (
                    "profiled_to_match_active_2023_completed_fertility_target"
                )
        fixed.atomic_write_csv(outdir / "profiled_parameters.csv", parameter_rows)
        ratios = np.asarray(
            [float(row["queue_B_over_E"]) for row in schedule]
            if profile_success
            else [],
            dtype=float,
        )
        prices = np.asarray(
            [float(row["price_ratio"]) for row in schedule]
            if profile_success
            else [],
            dtype=float,
        )
        maximum_index = int(np.argmax(ratios)) if ratios.size else None
        total_seconds = time.perf_counter() - started
        selected_summary = dict(best["summary"]["best_candidate"])
        input_provenance = {
            **contracts["hashes"],
            **code_hashes,
            **smoke_hashes,
            "ridge_report": str(contracts["report_path"]),
            "selected_task_summary": str(contracts["task_path"]),
            "source": str(contracts["source_path"]),
        }
        summary = {
            "schema": SCHEMA,
            "status": "complete_target_preserving_profile",
            "classification": "conditional_profile_diagnostic_not_calibration",
            "promotion_eligible": False,
            "run_stage": args.run_stage,
            "policy_case": "none",
            "factor_index": factor_index,
            "kappa_factor": factor,
            "factors_contract": list(FACTORS),
            "target_profile_success": profile_success,
            "target_moment": contracts["profile_target_moment"],
            "target_completed_fertility": target_value,
            "target_weight": contracts["profile_target_weight"],
            "target_tolerance": TARGET_TOLERANCE,
            "profiled_preference_change": float(best["row"]["preference_change"]),
            "profiled_terminal_completed_fertility": float(
                best["row"]["terminal_completed_fertility"]
            ),
            "profiled_completed_fertility_gap": float(best["row"]["target_gap"]),
            "transition_loss": float(selected_summary["transition_loss"]),
            "profile_trial_count": len(trials),
            "profile_monotone": monotone_trials(trials),
            "profile_bracket_found": bracket is not None,
            "profile_trials": [item["row"] for item in trials],
            "diagnostic_theta": selected_summary["theta"],
            "old_psi_child": float(best["summary"]["old_psi_child"]),
            "new_psi_child": float(selected_summary["new_psi_child"]),
            "target_fit": fit_rows,
            "parameters": parameter_rows,
            "usable_closed_root": bool(
                endpoint and endpoint.get("usable_closed_root", False)
            ),
            "between_steady_states_label_allowed": bool(
                endpoint and endpoint.get("between_steady_states_label_allowed", False)
            ),
            "closed_root_status": (
                endpoint.get("status") if endpoint else "not_run_profile_target_unresolved"
            ),
            "closed_root_price_ratio": (
                float(endpoint["price_ratio"])
                if endpoint and bool(endpoint.get("usable_closed_root", False))
                else None
            ),
            "maximum_B_over_E": (
                float(ratios[maximum_index]) if maximum_index is not None else None
            ),
            "maximum_B_over_E_price_ratio": (
                float(prices[maximum_index]) if maximum_index is not None else None
            ),
            "minimum_B_over_E": float(np.min(ratios)) if ratios.size else None,
            "closed_stationary_endpoint": endpoint,
            "runtime": {
                "profile_transition_seconds": sum(
                    float(item["row"]["elapsed_seconds"]) for item in trials
                ),
                "closed_endpoint_seconds": endpoint_seconds,
                "total_seconds": total_seconds,
            },
            "input_provenance": input_provenance,
        }
        fixed.atomic_write_json(outdir / "summary.json", summary)
        artifacts = (
            "summary.json",
            "profile_trials.csv",
            "profiled_target_fit.csv",
            "profiled_parameters.csv",
            "closed_stationary_schedule.csv",
            "closed_stationary_endpoint.json",
        )
        manifest = {
            "schema": "e5f_target_preserving_closed_root_factor_manifest_v1",
            "status": "complete_target_preserving_profile_manifest",
            "run_stage": args.run_stage,
            "factor_index": factor_index,
            "kappa_factor": factor,
            "profile_code_contract": code_hashes,
            "artifacts": fixed.artifact_hashes(outdir, artifacts),
        }
        fixed.atomic_write_json(outdir / "manifest.json", manifest)
        fixed.atomic_write_json(
            outdir / "latest_completed_case.json",
            {
                "status": "complete",
                "factor_index": factor_index,
                "kappa_factor": factor,
                "profile_success": profile_success,
                "preference_change": summary["profiled_preference_change"],
                "transition_loss": summary["transition_loss"],
                "maximum_B_over_E": summary["maximum_B_over_E"],
                "usable_closed_root": summary["usable_closed_root"],
                "elapsed_seconds": total_seconds,
            },
        )
        heartbeat.set_phase("complete")
        print(
            "TARGET_PRESERVING_ROOT_PROFILE_COMPLETE "
            f"factor={factor:g} delta={summary['profiled_preference_change']:.9f} "
            f"tfr={summary['profiled_terminal_completed_fertility']:.9f} "
            f"loss={summary['transition_loss']:.9f} "
            f"max_BE={summary['maximum_B_over_E']} root={summary['usable_closed_root']}",
            flush=True,
        )
    except Exception as error:
        fixed.atomic_write_json(
            outdir / "failure.json",
            {
                "status": "failed_target_preserving_profile",
                "factor_index": factor_index,
                "kappa_factor": factor,
                "error_type": type(error).__name__,
                "error": str(error),
                "traceback": traceback.format_exc(),
                "elapsed_seconds": time.perf_counter() - started,
            },
        )
        heartbeat.set_phase("failed")
        raise
    finally:
        heartbeat.stop()


def collect(args: argparse.Namespace) -> None:
    if not args.collect or args.factor_index is not None or args.run_stage != "collect":
        raise ProfileError("collection requires --collect, no factor index, and collect stage")
    if args.results_root is None:
        raise ProfileError("collection requires --results-root")
    outdir = args.output_dir.resolve()
    if outdir.exists():
        raise ProfileError(f"refusing to overwrite collection output: {outdir}")
    outdir.mkdir(parents=True)
    code_hashes = validate_profile_code(args)
    contracts = base_contracts(args)
    expected_provenance = {**contracts["hashes"], **code_hashes}
    target_value = float(contracts["profile_target_value"])
    smoke_hashes = validate_smoke_receipt(
        args, expected_provenance, target_value
    )
    rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    schedule_rows: list[dict[str, Any]] = []
    for index, factor in enumerate(FACTORS, start=1):
        factor_dir = args.results_root.resolve() / f"factor_{index:03d}"
        summary = fixed.read_json(factor_dir / "summary.json")
        manifest = fixed.read_json(factor_dir / "manifest.json")
        fixed.require_equal(summary.get("schema"), SCHEMA, f"factor {index} schema")
        fixed.require_equal(
            summary.get("status"),
            "complete_target_preserving_profile",
            f"factor {index} status",
        )
        fixed.require_equal(summary.get("run_stage"), "factor", f"factor {index} stage")
        fixed.require_equal(summary.get("factor_index"), index, f"factor {index} index")
        fixed.require_close(summary.get("kappa_factor"), factor, "factor value", 1e-14)
        fixed.require_close(
            summary.get("target_completed_fertility"),
            target_value,
            f"factor {index} target",
            1e-14,
        )
        fixed.require_equal(
            bool(summary.get("promotion_eligible", True)),
            False,
            f"factor {index} promotion status",
        )
        provenance = dict(summary.get("input_provenance") or {})
        for name, expected in {
            **expected_provenance,
            **smoke_hashes,
        }.items():
            fixed.require_equal(provenance.get(name), expected, f"factor {index} {name}")
        fixed.require_equal(
            manifest.get("schema"),
            "e5f_target_preserving_closed_root_factor_manifest_v1",
            f"factor {index} manifest schema",
        )
        fixed.require_equal(
            manifest.get("status"),
            "complete_target_preserving_profile_manifest",
            f"factor {index} manifest status",
        )
        fixed.require_equal(
            manifest.get("run_stage"), "factor", f"factor {index} manifest stage"
        )
        for name, expected in code_hashes.items():
            fixed.require_equal(
                dict(manifest.get("profile_code_contract") or {}).get(name),
                expected,
                f"factor {index} manifest {name}",
            )
        for name, expected in dict(manifest.get("artifacts") or {}).items():
            fixed.require_equal(
                fixed.file_sha256(factor_dir / name), expected, f"factor {index} {name}"
            )
        rows.append(
            {
                "factor_index": index,
                "kappa_factor": factor,
                "profile_success": bool(summary["target_profile_success"]),
                "preference_change": summary["profiled_preference_change"],
                "completed_fertility": summary[
                    "profiled_terminal_completed_fertility"
                ],
                "completed_fertility_gap": summary[
                    "profiled_completed_fertility_gap"
                ],
                "transition_loss": summary["transition_loss"],
                "profile_trials": summary["profile_trial_count"],
                "maximum_B_over_E": summary["maximum_B_over_E"],
                "usable_closed_root": summary["usable_closed_root"],
                "closed_root_price_ratio": summary["closed_root_price_ratio"],
                "runtime_seconds": summary["runtime"]["total_seconds"],
            }
        )
        factor_target_rows = fixed.read_csv(factor_dir / "profiled_target_fit.csv")
        if len(factor_target_rows) != fixed.TARGET_COUNT:
            raise ProfileError(
                f"factor {index} target-fit row count is {len(factor_target_rows)}"
            )
        loss_from_rows = sum(
            float(row["loss_contribution"]) for row in factor_target_rows
        )
        fixed.require_close(
            loss_from_rows,
            summary["transition_loss"],
            f"factor {index} target-loss identity",
            2e-8,
        )
        for row in factor_target_rows:
            target_rows.append({"factor_index": index, "kappa_factor": factor, **row})
        for row in fixed.read_csv(factor_dir / "closed_stationary_schedule.csv"):
            schedule_rows.append({"factor_index": index, "kappa_factor": factor, **row})
    fixed.atomic_write_csv(outdir / "frontier_summary.csv", rows)
    fixed.atomic_write_csv(outdir / "target_fit_full.csv", target_rows)
    fixed.atomic_write_csv(outdir / "closed_stationary_schedules_full.csv", schedule_rows)
    report = {
        "schema": REPORT_SCHEMA,
        "status": "complete_target_preserving_closed_root_frontier",
        "classification": "conditional_profile_diagnostic_not_calibration",
        "promotion_eligible": False,
        "target_moment": contracts["profile_target_moment"],
        "target_completed_fertility": target_value,
        "target_weight": contracts["profile_target_weight"],
        "target_tolerance": TARGET_TOLERANCE,
        "factors": list(FACTORS),
        "frontier": rows,
        "root_capable_factors": [
            float(row["kappa_factor"]) for row in rows if row["usable_closed_root"]
        ],
        "input_provenance": {
            **contracts["hashes"],
            **code_hashes,
            **smoke_hashes,
        },
    }
    fixed.atomic_write_json(outdir / "summary.json", report)
    fixed.atomic_write_json(
        outdir / "manifest.json",
        {
            "schema": "e5f_target_preserving_closed_root_report_manifest_v1",
            "status": "complete_target_preserving_closed_root_report_manifest",
            "profile_code_contract": code_hashes,
            "artifacts": fixed.artifact_hashes(
                outdir,
                (
                    "summary.json",
                    "frontier_summary.csv",
                    "target_fit_full.csv",
                    "closed_stationary_schedules_full.csv",
                ),
            ),
        },
    )


def main() -> None:
    args = parse_args()
    if args.collect:
        collect(args)
    else:
        run_factor(args)


if __name__ == "__main__":
    main()
