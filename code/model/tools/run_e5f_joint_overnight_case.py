#!/usr/bin/env python3
"""Observe and validate bounded E5F evaluations without changing scientific code.

Run this adapter from the frozen production snapshot. Every candidate invokes
the unchanged dated calibration driver. An observer copies the terminal state
after its existing measurement callback, allowing the standard graphs and a
checkpoint to accompany every completed candidate without a second solve.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
import os
import pickle
import sys
import threading
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "code/model"), str(ROOT / "code/model/tools")]

BUNDLE = "ce38de90a85de7102f4d462bd1f2618fbad17f649c5c57d840d5152e5917dff6"
SUPPORTED_BUNDLES = (BUNDLE,)
TARGET = "3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497"
SOURCE = "0afcb82d4735bd15aaa143ea04e3105a5d43df152122d02b983372102f20eef6"


SEARCH_DOMAIN = [{'lower': 0.94, 'name': 'beta_annual', 'transform': 'discount', 'upper': 0.9995}, {'lower': 0.02, 'name': 'kappa_fert', 'transform': 'log', 'upper': 50.0}, {'lower': 0.02, 'name': 'kappa_fert_continuation', 'transform': 'log', 'upper': 50.0}, {'lower': 0.1, 'name': 'chi', 'transform': 'log', 'upper': 5.0}, {'lower': 0.2, 'name': 'H0', 'transform': 'log', 'upper': 80.0}, {'lower': 0.0, 'name': 'theta0', 'transform': 'softzero', 'upper': 8.0}, {'lower': 0.02, 'name': 'theta1', 'transform': 'log', 'upper': 16.0}, {'lower': 0.1, 'name': 'hbar_child_rooms', 'transform': 'log', 'upper': 1.8}, {'lower': 0.0, 'name': 'first_birth_fixed_cost', 'transform': 'softzero', 'upper': 8.0}, {'lower': 0.0, 'name': 'hbar_first_child_jump', 'transform': 'softzero', 'upper': 2.0}, {'lower': -1.5, 'name': 'psi_child_change_2023', 'transform': 'asinh', 'upper': 0.2}]

def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def read_csv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def verify(path, expected):
    actual = digest(path)
    if actual != expected:
        raise RuntimeError(f"Hash mismatch: {path}: {actual} != {expected}")


def load_plan(path, expected):
    verify(path, expected)
    plan = read_json(path)
    if plan["schema"] != "e5f_joint_overnight_v1":
        raise RuntimeError("Unknown plan schema")
    if plan["source_sha256"] != SOURCE or plan["target_fingerprint"] != TARGET:
        raise RuntimeError("Plan changes the maintained source or target system")
    if plan["code_bundle_sha256"] not in SUPPORTED_BUNDLES:
        raise RuntimeError("Plan requests an unverified scientific bundle")
    if len(plan["cases"]) > 40 or len({c["id"] for c in plan["cases"]}) != len(plan["cases"]):
        raise RuntimeError("Unbounded or duplicated case list")
    if plan.get("search_domain") != SEARCH_DOMAIN or plan.get("first_child_jump_upper") != 2.0:
        raise RuntimeError("Plan changes the declared eleven-parameter domain")
    return plan


def validate_result(out, plan, case):
    import collect_e5f_transition_calibration as collector
    s = read_json(out / "summary.json")
    b = s["best_candidate"]
    expected = SimpleNamespace(
        expected_source_sha256=SOURCE,
        expected_target_set=plan["target_set"], expected_target_fingerprint=TARGET,
        expected_code_bundle_sha256=plan["code_bundle_sha256"], expected_model_profile="e5f-income-entry",
        expected_panel_seed=case["panel_seed"], expected_center_sha256=case["center_sha256"],
        expected_panel_design=case["panel_design"], expected_local_radius=case["radius"],
        expected_housing_supply_elasticity=.63, expected_tenure_choice_kappa=.005,
    )
    collector.validate_expected_contract(s, s["panel_design"], expected)
    collector.validate_renewal_accounting(s)
    collector.validate_calibration_scope(s)
    if s["target_count"] != 12 or s["transition_free_parameter_count"] != 11:
        raise RuntimeError("Target/parameter count changed")
    if s["panel_design"]["task_id"] != case["panel_task_id"] or s["panel_design"]["panel_size"] != case["panel_size"]:
        raise RuntimeError("Candidate generation metadata changed")
    for name, value in zip((r["name"] for r in s["panel_design"]["domain"]), s["panel_design"]["unit_vector"]):
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise RuntimeError(f"Invalid bounded coordinate {name}")
    if s["panel_design"]["domain"] != SEARCH_DOMAIN:
        raise RuntimeError("Actual search domain differs from the pinned plan")
    models = collector.validate_target_fit(read_csv(out / "target_fit_long.csv"), s, b)
    collector.validate_dated_housing_ledger(out, s, b, models)
    collector.population_bridge_contract(s["population_bridge"])
    checks = {
        "market": (b["max_market_residual"], 2e-4),
        "mass": (b["max_mass_residual"], 2e-10),
        "population": (abs(b["population_target_gap"]), 2e-10),
        "stationary_measurement": (s["stationary_measurement_max_abs_gap"], 2e-8),
        "childless_measurement": (abs(b["terminal_synthetic_childless_consistency_gap"]), 2e-10),
    }
    for name, (value, limit) in checks.items():
        if not math.isfinite(value) or value > limit:
            raise RuntimeError(f"Failed unchanged {name} gate: {value} > {limit}")
    return s, models, checks


def compare_reference(out, reference):
    """The inherited 13th row is excluded explicitly, never relabeled as 12."""
    ref = read_json(reference)
    now = read_json(out / "summary.json")
    reference_rows = {r["moment"]: r for r in read_csv(Path(reference).parent / "target_fit_long.csv")}
    for row in read_csv(out / "target_fit_long.csv"):
        for key in ("target", "model", "gap", "weight", "loss_contribution", "standardized_gap"):
            if float(row[key]) != float(reference_rows[row["moment"]][key]):
                raise RuntimeError(f"Exact reference fit mismatch: {row['moment']} / {key}")
    for key in ("theta", "new_psi_child"):
        if now["best_candidate"][key] != ref["best_candidate"][key]:
            raise RuntimeError(f"Exact reference parameter mismatch: {key}")
    if now["old_psi_child"] != ref["old_psi_child"]:
        raise RuntimeError("Exact old preference normalization mismatch")
    old_case = Path(reference).parent / "cases" / ref["best_candidate"]["candidate"]
    new_case = out / "cases" / now["best_candidate"]["candidate"]
    checked = 0
    # Reproduce all saved numeric historical path entries except elapsed times.
    old_rows, new_rows = read_csv(old_case / "transition_path.csv"), read_csv(new_case / "transition_path.csv")
    if len(old_rows) != 5 or len(new_rows) != 5:
        raise RuntimeError("Exact reference history length mismatch")
    for old, new in zip(old_rows, new_rows):
        for key in old:
            if any(word in key for word in ("elapsed", "seconds", "label", "scenario")):
                continue
            try:
                left, right = float(old[key]), float(new[key])
            except ValueError:
                continue
            if not (left == right or (math.isnan(left) and math.isnan(right))):
                raise RuntimeError(f"Exact reference historical mismatch: {key}: {left} != {right}")
            checked += 1
    return {"exact_twelve_row_fit": True, "exact_parameters": True,
            "exact_numeric_history_entries": checked, "reference_sha256": digest(reference)}


def compare_bridge(out, reference):
    """Compare fixed-jump evidence with its freely estimated-domain replay.

    Changing the square-root search transform can introduce last-bit rounding.
    Target and weight identities are exact; numerical bridge tolerances are
    1e-10 for fit/history entries and 1e-12 for physical parameters. Subsequent
    same-generator repetitions remain strictly exact.
    """
    ref, now = read_json(reference), read_json(Path(out) / "summary.json")
    reference_rows = {r["moment"]: r for r in read_csv(Path(reference).parent / "target_fit_long.csv")}
    maximum = 0.0
    for row in read_csv(Path(out) / "target_fit_long.csv"):
        old = reference_rows[row["moment"]]
        for key in ("target", "weight"):
            if float(row[key]) != float(old[key]):
                raise RuntimeError("Bridge changes empirical targets or weights")
        for key in ("model", "gap", "loss_contribution", "standardized_gap"):
            gap = abs(float(row[key]) - float(old[key])); maximum = max(maximum, gap)
            if not math.isfinite(gap) or gap > 1e-10:
                raise RuntimeError(f"Source/domain bridge fit mismatch: {row['moment']} {key} {gap}")
    left, right = ref["best_candidate"]["theta"], now["best_candidate"]["theta"]
    if set(left) != set(right):
        raise RuntimeError("Bridge physical parameter names differ")
    for key in left:
        if abs(float(left[key])-float(right[key])) > 1e-12:
            raise RuntimeError(f"Bridge parameter mismatch: {key}")
    for old, new in [(ref["old_psi_child"], now["old_psi_child"]),
                     (ref["best_candidate"]["new_psi_child"], now["best_candidate"]["new_psi_child"])]:
        if abs(float(old)-float(new)) > 1e-12:
            raise RuntimeError("Bridge preference normalization differs")
    def history(base, summary):
        return read_csv(base / "cases" / summary["best_candidate"]["candidate"] / "transition_path.csv")
    oldrows, newrows = history(Path(reference).parent, ref), history(Path(out), now)
    if len(oldrows) != len(newrows) or len(oldrows) != 5:
        raise RuntimeError("Bridge historical length mismatch")
    count = 0; history_max = 0.0
    for old, new in zip(oldrows, newrows):
        for key in old:
            if any(word in key for word in ("elapsed", "seconds", "label", "scenario")): continue
            try: a, b = float(old[key]), float(new[key])
            except ValueError: continue
            if math.isnan(a) and math.isnan(b): count += 1; continue
            gap = abs(a-b); history_max = max(history_max, gap)
            if not math.isfinite(gap) or gap > 1e-10:
                raise RuntimeError(f"Bridge historical mismatch: {key} {gap}")
            count += 1
    return {"bridge_twelve_row_fit": True, "bridge_parameters": True,
            "numeric_history_entries": count, "maximum_fit_difference": maximum,
            "maximum_history_difference": history_max,
            "reference_sha256": digest(reference), "scope": "search-domain/source bridge, not exact repeat"}


def run_case(args):
    plan_path = args.plan.resolve()
    plan = load_plan(plan_path, args.plan_sha256)
    verify(__file__, plan["adapter_sha256"])
    for name, sha in plan["helper_sha256"].items():
        verify(ROOT / "code/model/tools" / name, sha)
    for name, sha in plan.get("input_sha256", {}).items():
        verify(plan_path.parent / name, sha)
    verify(plan["source"], SOURCE)
    case = next(c for c in plan["cases"] if c["id"] == args.case_id)
    center = plan_path.parent / case["center"]
    verify(center, case["center_sha256"])
    if time.time() > plan["launch_deadline_epoch"]:
        raise RuntimeError("Predeclared launch deadline exceeded")
    out = plan_path.parent / case["output"]
    if out.exists() and any(out.iterdir()):
        raise RuntimeError(f"Refusing nonempty case output: {out}")
    out.mkdir(parents=True, exist_ok=True)
    start = time.time()
    state = {"phase": "starting", "completed_periods": 0}
    stop = threading.Event()
    def heartbeat():
        while not stop.is_set():
            write_json(out / "heartbeat.json", {**state, "elapsed_seconds": time.time()-start,
                       "epoch": time.time(), "case_id": case["id"], "slurm_job_id": os.getenv("SLURM_JOB_ID")})
            stop.wait(60)
    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        import run_e5f_transition_calibration as calibration
        import run_e5f_independent_numerical_audit as audit
        from intergen_eqscale_seq_optimized import solver
        actual = calibration.code_fingerprint_contract(solver)
        if actual["bundle_sha256"] != plan["code_bundle_sha256"]:
            raise RuntimeError("Frozen scientific bundle has drifted")
        original = calibration.transition.run_birth_vintage_scenario
        captured = {}
        def observed_scenario(**kwargs):
            original_observer = kwargs["period_observer"]
            def observer(period, evaluation, parameters, grid, shared):
                original_observer(period, evaluation, parameters, grid, shared)
                state.update(phase=f"historical_period_{period}_complete", completed_periods=period+1)
                if period == 4:
                    captured.update(copy.deepcopy(dict(parameters=parameters, b_grid=grid,
                        evaluation=evaluation, shared=shared, supply_rule=kwargs["supply_rule"])))
            kwargs["period_observer"] = observer
            return original(**kwargs)
        calibration.transition.run_birth_vintage_scenario = observed_scenario
        argv = [str(calibration.__file__), "--source", plan["source"],
            "--expected-source-sha256", SOURCE, "--expected-target-set", plan["target_set"],
            "--expected-target-fingerprint", TARGET, "--expected-code-bundle-sha256", plan["code_bundle_sha256"],
            "--outdir", str(out), "--model-profile", "e5f-income-entry",
            "--estimate-first-child-room-jump", "--first-child-room-jump-upper", "2.0",
            "--housing-supply-elasticity", "0.63",
            "--fixed-tenure-choice-kappa", "0.005", "--replacement-fertility", "2.1",
            "--old-completed-fertility-target", "2.1", "--outside-origin-entry-share", "0.169",
            "--market-tol", "0.0002", "--market-max-iter", "30", "--nb", "120",
            "--post-2023-periods", "0", "--policy-case", "none",
            "--panel-center-json", str(center), "--panel-task-id", str(case["panel_task_id"]),
            "--panel-size", str(case["panel_size"]), "--panel-design", case["panel_design"],
            "--panel-seed", str(case["panel_seed"]), "--panel-local-radius", str(case["radius"])]
        write_json(out / "execution_contract.json", {"argv": argv, "plan_sha256": args.plan_sha256,
                   "case": case, "adapter_sha256": digest(__file__), "start_epoch": start})
        state["phase"] = "unchanged_scientific_driver"
        previous_argv = sys.argv
        try:
            sys.argv = argv
            calibration.main()
        finally:
            sys.argv = previous_argv
            calibration.transition.run_birth_vintage_scenario = original
        summary, models, checks = validate_result(out, plan, case)
        reference_receipt = None
        if case.get("bridge_reference"):
            ref = plan_path.parent / case["bridge_reference"]
            verify(ref, case["bridge_reference_sha256"])
            reference_receipt = compare_bridge(out, ref)
        if case.get("reference"):
            ref = plan_path.parent / case["reference"]
            verify(ref, case["reference_sha256"])
            reference_receipt = compare_reference(out, ref)
        if not captured:
            raise RuntimeError("Terminal state observer did not run")
        state["phase"] = "checkpoint_and_standard_graphs"
        captured.update(adapter_sha256=digest(__file__), calibration_summary_sha256=digest(out / "summary.json"))
        checkpoint = out / "dated_state.pkl.gz"
        with gzip.open(checkpoint, "wb", compresslevel=1) as stream:
            pickle.dump(captured, stream, protocol=5)
        audit.standard_diagnostics(captured, out, validate_production_young=False)
        graphs = sorted((out / "standard_diagnostics").glob("*.png"))
        if len(graphs) != 17:
            raise RuntimeError(f"Expected 17 standard PNGs, found {len(graphs)}")
        budget = audit.budget_audit(captured, out)
        policy_arrays = audit.policy_array_audit(captured, out)
        for name, bounds in policy_arrays["probabilities"].items():
            if bounds["nonfinite"] or bounds["minimum"] < 0 or bounds["maximum"] > 1:
                raise RuntimeError(f"Invalid {name} probability array")
        artifacts = {str(p.relative_to(out)): digest(p) for p in [out / "summary.json", out / "target_fit_long.csv",
                     out / "parameter_table.csv", checkpoint, *graphs]}
        receipt = {"status": "complete", "case_id": case["id"], "label": case["label"],
            "plan_sha256": args.plan_sha256, "loss": summary["best_candidate"]["transition_loss"],
            "models": models, "gates": {k: {"value": v, "limit": l, "passed": True} for k,(v,l) in checks.items()},
            "reference": reference_receipt, "artifact_sha256": artifacts, "standard_graph_count": len(graphs),
            "budget_diagnostic": budget, "policy_array_diagnostic": policy_arrays,
            "elapsed_seconds": time.time()-start, "production_promoted": False}
        write_json(out / "case_receipt.json", receipt)
        state["phase"] = "complete"
    except Exception as error:
        state["phase"] = "failed"
        write_json(out / "adapter_failure.json", {"error": str(error), "type": type(error).__name__,
                   "traceback": traceback.format_exc(), "elapsed_seconds": time.time()-start})
        raise
    finally:
        stop.set()
        thread.join(timeout=2)
        write_json(out / "heartbeat.json", {**state, "elapsed_seconds": time.time()-start, "epoch": time.time()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--case-id", type=int, required=True)
    run_case(parser.parse_args())


if __name__ == "__main__":
    main()
