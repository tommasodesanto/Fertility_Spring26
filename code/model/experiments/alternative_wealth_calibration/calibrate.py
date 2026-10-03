"""Bounded alternative-clock calibration against model-matched PSID wealth.

Run one arm and one chain per fresh Python process. No model source is edited.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import importlib.util
import inspect
import json
import math
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[4]
PACKETS = ROOT / "output/model/fixed_reference_economics_20260928"
V2 = PACKETS / "normalized_calibration_v2"
TIMING = PACKETS / "purchase_timing_sandbox_v1"
SELECTION = PACKETS / "soft_timing_review_v1/soft_selected.json"
PLAN = PACKETS / "soft_timing_calibration_20261002_v1/driver_plan.json"
RESERVE = 1800.0
MAX_CALLS = 250
MAX_LIFECYCLE = 32
PENALTY = 1e12
OLD_WEALTH_TARGET = "6.92658379107299"
NEW_WEALTH_TARGET = "4.45838713455674"
WEALTH_MOMENT = "wealth_earnings"


def experimental_contract(v2, original_manifest):
    """Replace one empirical value while preserving its weight and model observer."""
    base = v2.CONFIG["base_target_contract"]
    changed = [dict(row) for row in base]
    matches = [row for row in changed if row["moment"] == WEALTH_MOMENT]
    v2.inputs.require(len(changed) == 14 and len(matches) == 1, "Wealth row missing")
    row = matches[0]
    v2.inputs.require(row["role"] == "scored" and row["target"] == OLD_WEALTH_TARGET and
                      row["weight"] == "7.595098472533724", "Original wealth contract drift")
    row["target"] = NEW_WEALTH_TARGET
    v2.inputs.require(sum(r["role"] == "scored" for r in changed) == 10,
                      "Hard target count drift")
    identity = v2.inputs.canonical(changed)
    weight = v2.inputs.canonical(dict(base_contract=changed, multipliers={}))
    v2.inputs.require(identity != original_manifest["target_fingerprint"] and
                      weight != original_manifest["weight_fingerprint"],
                      "Experimental target was not separated")
    return changed, dict(original_manifest, target_fingerprint=identity,
                         weight_fingerprint=weight,
                         original_target_fingerprint=original_manifest["target_fingerprint"],
                         original_weight_fingerprint=original_manifest["weight_fingerprint"])


def rescore_wealth(v2, rows, old_residual, out):
    """Preserve native report and rescore its one changed empirical target."""
    v2.inputs.require(v2.native.target_identity(rows) == v2.CONFIG["base_target_contract"],
                      "Native old target identity drift")
    revised = [dict(row) for row in rows]
    old = next(row for row in revised if row["moment"] == WEALTH_MOMENT)
    old_model = float(old["model"])
    old_gap = float(old["gap"])
    old_loss = float(old["loss_contribution"])
    w = float(old["weight"])
    v2.inputs.require(abs(old_gap - (old_model - float(OLD_WEALTH_TARGET))) < 1e-10,
                      "Native wealth gap sign drift")
    v2.inputs.require(abs(old_loss - w * old_gap * old_gap) < 1e-8,
                      "Native wealth loss arithmetic drift")
    old["target"] = NEW_WEALTH_TARGET
    old["gap"] = str(old_model - float(NEW_WEALTH_TARGET))
    old["loss_contribution"] = str(w * float(old["gap"]) ** 2)
    contract, experiment = experimental_contract(v2, dict(
        target_fingerprint=v2.inputs.canonical(v2.CONFIG["base_target_contract"]),
        weight_fingerprint=v2.weight_fingerprint({})))
    v2.inputs.require(v2.native.target_identity(revised) == contract,
                      "Rescored target identity drift")
    scored = [row for row in revised if row["role"] == "scored"]
    rr = np.asarray([math.sqrt(float(row["weight"])) * float(row["gap"])
                     for row in scored], float)
    v2.inputs.require(rr.shape == (10,) and np.isfinite(rr).all(),
                      "Rescored residual invalid")
    old_rr = np.asarray(old_residual, float)
    original_scored = [row for row in rows if row["role"] == "scored"]
    v2.inputs.require(np.allclose(old_rr,
                      [math.sqrt(float(row["weight"])) * float(row["gap"])
                       for row in original_scored], rtol=0, atol=1e-10),
                      "Native residual differs from native target rows")
    j = [row["moment"] for row in scored].index(WEALTH_MOMENT)
    v2.inputs.require(np.array_equal(np.delete(rr, j), np.delete(old_rr, j)),
                      "Nonwealth residual changed")
    old_total = sum(float(row["loss_contribution"] or 0) for row in rows)
    new_total = sum(float(row["loss_contribution"] or 0) for row in revised)
    expected = old_total - old_loss + w * (old_model - float(NEW_WEALTH_TARGET)) ** 2
    v2.inputs.require(abs(new_total - expected) < 1e-8 and
                      abs(new_total - float(rr @ rr)) < 1e-8,
                      "Wealth-only loss identity failed")
    v2.native.table(out / "target_fit_new_contract.csv", revised)
    v2.write(out / "wealth_rescore_receipt.json", dict(
        original_loss=old_total, experimental_loss=new_total,
        native_original_target_fingerprint=experiment["original_target_fingerprint"],
        experimental_target_fingerprint=experiment["target_fingerprint"],
        native_original_weight_fingerprint=experiment["original_weight_fingerprint"],
        experimental_weight_fingerprint=experiment["weight_fingerprint"],
        fingerprint_change_due_to_target_value_only=True,
        native_proposed_parameters_target_hash_is_original=True,
        original_wealth_target=float(OLD_WEALTH_TARGET),
        experimental_wealth_target=float(NEW_WEALTH_TARGET),
        model_wealth=old_model, unchanged_weight=w,
        nonwealth_residuals_exact=True,
        model_observer_unchanged=True,
        experimental_weight_is_control_not_new_precision=True))
    return revised, rr


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def checked_inputs(arm):
    # The original arm must never import the alternative household overlay.
    timing_driver = None
    if arm == "alternative":
        # Import order is part of the existing isolated experiment's contract.
        timing_driver = load(TIMING / "run.py", "matched_timing_driver")
        v2 = timing_driver.v2
    else:
        sys.path.insert(0, str(V2))
        v2 = load(V2 / "run_psi.py", "matched_normalized_v2")
    v2.native.verify_sources()
    for rel, digest in json.loads((V2 / "source_pins.json").read_text()).items():
        v2.inputs.require(sha(ROOT / rel) == digest, "Normalized source drift: " + rel)
    manifest = json.loads((TIMING / "manifest.json").read_text())
    selected_file = json.loads(SELECTION.read_text())
    source = ROOT / selected_file["source"]
    v2.inputs.require(sha(source) == selected_file["source_sha256"], "Soft selection source drift")
    selected = selected_file["selected"]
    v2.inputs.require(json.loads(source.read_text())[selected_file["source_key"]]["best"] == selected,
                      "Soft selected checkpoint drift")
    v2.inputs.require(v2.inputs.canonical(v2.CONFIG["base_target_contract"]) == manifest["target_fingerprint"],
                      "Target fingerprint drift")
    v2.inputs.require(v2.weight_fingerprint({}) == manifest["weight_fingerprint"] == selected["weight_fingerprint"],
                      "Weight fingerprint drift")
    v2.inputs.require(v2.native.target_identity(selected["target_fit"]) == v2.CONFIG["base_target_contract"],
                      "Selected target contract drift")
    if arm == "alternative":
        _, checked_manifest = timing_driver.verify()
        v2.inputs.require(checked_manifest == manifest, "Alternative source manifest drift")
    changed, experimental_manifest = experimental_contract(v2, manifest)
    experimental_manifest["target_contract"] = changed
    return v2, timing_driver, experimental_manifest, selected


def install_timing_observer(v2, timing_driver, out, P):
    """Use the same calendar forward map and independently adjusted audit as evaluate_selected.py."""
    original_install = v2.native.install_observer_metadata

    def timing_metadata(ge, arm, source_bounds):
        original_install(ge, arm, source_bounds)
        observe = ge.observe_price

        def timing_observe(ctx, live, label, *, final=False):
            rt = ctx["prepared"].rt
            if not getattr(rt["accounting"], "_transaction_timing_installed", False):
                cal = rt["primitive"].pf.calendar
                old_map = cal.model.build_forward_tenure_transition_maps
                new_map = timing_driver.sandbox_household.build_forward_tenure_transition_maps
                old_audit = rt["accounting"]._inherited_audit()
                before = inspect.getsource(old_audit)
                replacements = {
                    "x = grid if old == new else grid + sale[old] - costs[new]":
                        "x = grid if old == new else grid + (sale[old] - costs[new]) / float(P.R_gross)",
                    "invalid = x + y / float(P.R_gross) < floor - 1e-10":
                        "invalid = float(P.R_gross) * x + y < floor - 1e-10",
                }
                after = before
                for old, new in replacements.items():
                    v2.inputs.require(after.count(old) == 1, "Purchase audit source drift: " + old)
                    after = after.replace(old, new)
                namespace = dict(old_audit.__globals__)
                audit_path = out / "timing_purchase_audit.py"
                audit_path.write_text(after)
                (out / "timing_purchase_audit.diff").write_text("".join(difflib.unified_diff(
                    before.splitlines(True), after.splitlines(True),
                    fromfile="original", tofile="transaction_timing")))
                exec(compile(after, str(audit_path), "exec"), namespace)
                rt["accounting"]._INHERITED = namespace["audit_purchase_accounting"]
                cal.model.build_forward_tenure_transition_maps = new_map
                rt["model"].build_forward_tenure_transition_maps = new_map
                rt["accounting"]._transaction_timing_installed = True
                v2.write(out / "timing_observer_receipt.json", dict(
                    original_map_file=inspect.getsourcefile(old_map),
                    experiment_map_file=inspect.getsourcefile(new_map),
                    original_audit_sha256=hashlib.sha256(before.encode()).hexdigest(),
                    timing_audit_sha256=hashlib.sha256(after.encode()).hexdigest(),
                    acceptance_tolerances_unchanged=True,
                    audit_budget="R*b + net_sale - purchase + income - consumption - costs",
                    unrelated_audit_statements_unchanged=True))
            return observe(ctx, live, label, final=final)

        ge.observe_price = timing_observe

    v2.native.install_observer_metadata = timing_metadata


class BudgetStop(Exception):
    pass


def starts(selected, bounds, coordinates):
    """One selected point and three common deterministic, modest relative jitters."""
    base = selected["parameters"]
    seeds = [dict(base)]
    patterns = (
        {"chi": .025, "first_birth_fixed_cost": -.035, "h_P": -.015, "psi_child": .025},
        {"beta_annual": -.0015, "kappa_fert": .04, "tenure_choice_kappa": -.035,
         "child_benefit_curvature": .04},
        {"theta0": .04, "kappa_fert_continuation": -.035, "h_P": .012,
         "psi_child": -.025},
    )
    for pattern in patterns:
        point = dict(base)
        for key, shift in pattern.items():
            point[key] = min(bounds[key][1], max(bounds[key][0], base[key] * (1 + shift)))
        seeds.append(point)
    assert len(seeds) == 4 and all(set(s) == set(coordinates) for s in seeds)
    return seeds


def load_starts_table(path, expected_sha256, selected, manifest, bounds, coordinates):
    """Authenticate an external matched start table before model initialization."""
    path = Path(path).resolve()
    if sha(path) != expected_sha256:
        raise RuntimeError("Expanded start table SHA-256 drift")
    data = json.loads(path.read_text())
    if data["source_checkpoint_sha256"] != json.loads(SELECTION.read_text())["source_sha256"]:
        raise RuntimeError("Expanded starts use a different selected checkpoint")
    if (data["target_fingerprint"] != manifest["target_fingerprint"] or
            data["weight_fingerprint"] != manifest["weight_fingerprint"]):
        raise RuntimeError("Expanded starts target or weight fingerprint drift")
    if set(data["bounds"]) != set(bounds):
        raise RuntimeError("Expanded starts bound coordinates drift")
    for key in bounds:
        values = data["bounds"][key]
        if len(values) != 2 or tuple(float(v) for v in values) != tuple(bounds[key]):
            raise RuntimeError("Expanded starts bound drift: " + key)
    rows = data["starts"]
    if not isinstance(rows, list) or not rows:
        raise RuntimeError("Expanded start table is empty")
    provenance = data["start_provenance"]
    if not isinstance(provenance, list) or len(provenance) != len(rows):
        raise RuntimeError("Expanded per-start provenance missing")
    historical = json.loads((ROOT / json.loads(SELECTION.read_text())["source"]).read_text())
    result = []
    for index, row in enumerate(rows):
        if set(row) != set(coordinates):
            raise RuntimeError(f"Expanded start {index} coordinate drift")
        point = {}
        for key in coordinates:
            value = row[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise RuntimeError(f"Expanded start {index} nonfinite coordinate: {key}")
            point[key] = float(value)
        for key, value in point.items():
            low, high = bounds[key]
            if not low <= value <= high:
                raise RuntimeError(f"Expanded start {index} outside bound: {key}")
        source = provenance[index]
        if not isinstance(source, dict) or source.get("group") not in ("legacy", "historical", "near_selected", "broad"):
            raise RuntimeError(f"Expanded start {index} provenance drift")
        if source["group"] == "historical":
            key = source.get("source_key")
            if key not in historical or historical[key]["best"]["parameters"] != point:
                raise RuntimeError(f"Expanded historical start {index} source drift")
        result.append(point)
    if result[0] != selected["parameters"]:
        raise RuntimeError("Expanded first start must equal authenticated selected point")
    return result


def completion_receipt(search, status, **fields):
    """Override the provisional search status without duplicate keyword arguments."""
    result = dict(search)
    result.update(fields)
    result["status"] = status
    # Match the JSON writer's actual serialization contract before checkpointing.
    json.dumps(result, allow_nan=False)
    return result


def selected_repeat_path(report):
    report = Path(report)
    if report.name != "selected_root":
        raise RuntimeError("Native report is not the selected_root directory")
    repeat = report.parent / "selected_repeat_final"
    if not repeat.is_dir():
        raise RuntimeError("Native selected_repeat_final directory is missing")
    return repeat


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", choices=("alternative",), required=True)
    parser.add_argument("--chain", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--deadline-epoch", type=float, required=True)
    parser.add_argument("--mock-smoke", action="store_true")
    parser.add_argument("--preflight-evaluator", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--postcheck-only", action="store_true")
    parser.add_argument("--search-receipt", type=Path)
    parser.add_argument("--starts-file", type=Path)
    parser.add_argument("--starts-file-sha256")
    args = parser.parse_args()
    if args.postcheck_only != (args.search_receipt is not None):
        raise RuntimeError("Postcheck mode requires exactly one search receipt")
    if (args.starts_file is None) != (args.starts_file_sha256 is None):
        raise RuntimeError("Expanded starts require both file and SHA-256")
    if args.out.exists():
        raise RuntimeError("Refusing existing output directory")
    args.out.mkdir(parents=True)
    out = args.out.resolve()
    start = time.time()
    deadline = min(float(args.deadline_epoch), start + 21600)
    if deadline <= start + RESERVE and not (args.mock_smoke or args.postcheck_only):
        raise RuntimeError("Deadline leaves no 1800-second native reserve")
    if args.postcheck_only and deadline <= start:
        raise RuntimeError("Native postcheck child has no remaining deadline")
    v2, timing_driver, manifest, selected = checked_inputs(args.arm)
    lane = "floor_s0"
    _, bounds, _ = v2.inputs.seed_and_bounds(lane)
    bounds = {k: tuple(value) for k, value in bounds.items()}
    v2.inputs.require(bounds["h_P"] == (.1, 2.3), "Inherited h_P bound drift")
    bounds["h_P"] = (.1, 2.6)
    bounds["psi_child"] = tuple(v2.CONFIG["psi_bounds"])
    coordinates = tuple(v2.inputs.parameters(lane)) + ("psi_child",)
    v2.inputs.require(len(coordinates) == 10 and set(coordinates) == set(selected["parameters"]),
                      "Ten-coordinate contract drift")
    if args.starts_file is None:
        all_starts = starts(selected, bounds, coordinates)
    else:
        all_starts = load_starts_table(args.starts_file, args.starts_file_sha256,
                                       selected, manifest, bounds, coordinates)
    if not 0 <= args.chain < len(all_starts):
        raise RuntimeError("Chain index outside authenticated start table")
    seed = all_starts[args.chain]
    v2.inputs.check_point(seed, bounds)
    v2.inputs.LANES[lane].update(seed=seed, bounds=bounds, free_coordinates=list(coordinates))
    contract = dict(arm=args.arm, chain=args.chain, all_starts=all_starts, seed=seed,
        free_coordinates=list(coordinates), bounds=bounds,
        selected_source_sha256=json.loads(SELECTION.read_text())["source_sha256"],
        selection_receipt_sha256=sha(SELECTION),
        starts_file=str(args.starts_file.resolve()) if args.starts_file else None,
        starts_file_sha256=args.starts_file_sha256,
        starts_count=len(all_starts),
        normalized_source_pins_sha256=sha(V2 / "source_pins.json"),
        timing_manifest_sha256=sha(TIMING / "manifest.json") if args.arm == "alternative" else None,
        target_fingerprint=manifest["target_fingerprint"], weight_fingerprint=manifest["weight_fingerprint"],
        objective_calls_max=MAX_CALLS, lifecycle_solves_per_case_max=MAX_LIFECYCLE,
        reserve_seconds=RESERVE, deadline_epoch=deadline,
        economic_change=("adopted soft purchase financing; original interest clock" if args.arm == "original"
                         else "adopted soft purchase financing; experimental post-interest transaction clock"),
        other_economics="N0=1; physical parent room floor; constant alpha; no A(m); nonnegative-mean entry; phi=.8; unchanged ten coordinates; experimental matched-wealth target with original objective weight",
        no_auto_retry=True, not_adopted_calibration=True)
    v2.write(out / "start_contract.json", contract)
    lower, spans, initial = v2.simplex(seed, bounds, coordinates)
    v2.write(out / "search_contract.json", dict(method="bounded Nelder-Mead", scaled_bounds=[0., 1.],
        initial_simplex=initial.tolist(), physical_simplex_steps=v2.steps(seed, bounds, coordinates),
        max_objective_calls=2 if args.smoke else MAX_CALLS, final_reserve_seconds=RESERVE))
    def heartbeat(status, **more):
        v2.write(out / "heartbeat.json", dict(epoch=time.time(), status=status,
            arm=args.arm, chain=args.chain, **more))
    heartbeat("initialized", completed_full_ge=0, objective_calls=0)
    if args.mock_smoke:
        # Exercise the actual optimizer/checkpoint loop without imports of model runtimes.
        calls = []
        def mock(x):
            point = {k: float(lower[j] + spans[j] * x[j]) for j, k in enumerate(coordinates)}
            v2.inputs.check_point(point, bounds)
            row = dict(label=f"{len(calls):04d}_nm", parameters=point, status="passed",
                       loss=float(np.sum((np.asarray(x) - initial[0] - .01) ** 2)),
                       lifecycle_solves=0, mock=True)
            calls.append(row)
            v2.write(out / "latest_completed.json", dict(latest=row, completed_full_ge=len(calls)))
            v2.write(out / "best_so_far.json", dict(best=min(calls, key=lambda r: r["loss"])))
            heartbeat("mock_case_completed", completed_full_ge=len(calls), objective_calls=len(calls))
            if len(calls) >= 2:
                raise BudgetStop("mock_two_case_limit")
            return row["loss"]
        try:
            minimize(mock, initial[0], method="Nelder-Mead", bounds=[(0., 1.)] * 10,
                     options=dict(initial_simplex=initial, maxfev=2, maxiter=2, adaptive=True))
        except BudgetStop:
            pass
        v2.inputs.require(len(calls) == 2 and all(r["lifecycle_solves"] == 0 for r in calls),
                          "Mock did not exercise two cases")
        # Exercise both terminal receipt branches and the native sibling-path rule.
        mock_root = out / "mock_path_contract/phase_b_ge/selected_root"
        mock_repeat = mock_root.parent / "selected_repeat_final"
        mock_root.mkdir(parents=True)
        mock_repeat.mkdir()
        v2.inputs.require(selected_repeat_path(mock_root) == mock_repeat,
                          "Selected native repeat path contract drift")
        provisional = dict(status="provisional_search_finished", selected=calls[0],
                           objective_calls=len(calls), lifecycle_solves=0)
        v2.write(out / "mock_no_candidate_receipt.json",
                 completion_receipt(provisional, "no_admissible_candidate", selected=None))
        v2.write(out / "completed.json", completion_receipt(provisional,
            "mock_loop_passed_zero_solves", selected=calls[0], contract=contract))
        return
    P, grid = v2.inputs.proposal(lane)
    P, entry = v2.inputs.entry(P, grid, "nonnegative_mean")
    v2.inputs.require(P.native_purchase_income and P.native_due_stayer_credit and not P.joint_nested_choice,
                      "Purchase or choice contract drift")
    v2.inputs.require(P.N_target == 1. and P.R_gross > 1. and np.allclose(P.phi, .8, rtol=0., atol=1e-12),
                      "Population, interest, or financed-share drift")
    v2.write(out / "input_contract.json", dict(entry=entry, finance_share=np.asarray(P.phi).tolist(),
        normalized_population=float(P.N_target), original_selection_price=float(selected["price"]),
        target_contract=manifest["target_contract"], target_fingerprint=manifest["target_fingerprint"],
        weight_fingerprint=manifest["weight_fingerprint"]))
    if args.arm == "alternative":
        install_timing_observer(v2, timing_driver, out, P)
    Q = v2.native.utility_checks(P, grid, lane, out)
    if args.postcheck_only:
        search_path = args.search_receipt.resolve()
        search = json.loads(search_path.read_text())
        v2.inputs.require(search["status"] == "provisional_search_finished", "Search status drift")
        v2.inputs.require(search["arm"] == args.arm and search["chain"] == args.chain,
                          "Postcheck arm or chain drift")
        v2.inputs.require(search["target_fingerprint"] == manifest["target_fingerprint"] and
                          search["weight_fingerprint"] == manifest["weight_fingerprint"],
                          "Postcheck target or weight drift")
        v2.inputs.require(search["starts_file_sha256"] == args.starts_file_sha256,
                          "Postcheck start-table snapshot drift")
        chosen = search["selected"]
        v2.inputs.require(chosen is not None and chosen["status"] == "passed" and
                          chosen["weight_fingerprint"] == manifest["weight_fingerprint"],
                          "Missing authenticated selected candidate")
        v2.inputs.check_point(chosen["parameters"], bounds)
        v2.inputs.require(set(chosen["parameters"]) == set(coordinates),
                          "Postcheck parameter coordinate drift")
        native_evaluate = v2.normalized_objective.make_evaluator(out, lane, Q, grid,
            deadline, float(chosen.get("price", selected["price"])), native_runner=v2.native,
            exploratory=False)
        if args.preflight_evaluator:
            v2.write(out / "completed.json", dict(status="full_native_postcheck_initialized_zero_solves",
                lifecycle_solves=0, search_receipt_sha256=sha(search_path),
                starts_file_sha256=args.starts_file_sha256,
                target_fingerprint=manifest["target_fingerprint"],
                weight_fingerprint=manifest["weight_fingerprint"]))
            return
        verification = native_evaluate("selected_postcheck", chosen["parameters"], deadline)
        if verification["status"] != "passed":
            raise RuntimeError("Selected native postcheck failed: " + str(verification))
        report = Path(verification["report"])
        native_fits = v2.native.readtable(report / "target_fit.csv")
        fits, rr = rescore_wealth(v2, native_fits, verification["residual"], report)
        verification["native_original_residual"] = verification["residual"]
        verification["residual"] = rr.tolist()
        parameters = v2.native.readtable(report / "parameters.csv")
        v2.inputs.require(len(fits) == 14 and len(parameters) == 31, "Native reporting row count drift")
        v2.inputs.require(v2.native.target_identity(fits) == manifest["target_contract"],
                          "Native selected target contract drift")
        v2.inputs.require(len(list((report / "standard_diagnostics").glob("*.png"))) == 17,
                          "Native selected diagnostic plot count drift")
        native_loss = float(sum(float(r["loss_contribution"] or 0) for r in fits))
        v2.inputs.require(abs(native_loss - float(rr @ rr)) < 1e-8, "Native loss arithmetic drift")
        repeat = v2.native.compare_repeated(report, selected_repeat_path(report))
        repeat_report = selected_repeat_path(report)
        repeat_native_fits = v2.native.readtable(repeat_report / "target_fit.csv")
        repeat_fits, repeat_rr = rescore_wealth(v2, repeat_native_fits,
            v2.native.residual(repeat_native_fits), repeat_report)
        v2.inputs.require(repeat_fits == fits and np.array_equal(repeat_rr, rr),
                          "Rescored final target table differs across exact repeat")
        repeat["experimental_target_rows"] = 14
        repeat["experimental_target_fit_exact"] = True
        v2.write(out / "completed.json", dict(status="full_native_postcheck_passed",
            selected_postcheck=verification, native_loss=native_loss,
            target_fit=fits, original_target_fit=native_fits, parameters=parameters, repeat=repeat,
            search_receipt_sha256=sha(search_path), target_fingerprint=manifest["target_fingerprint"],
            starts_file_sha256=args.starts_file_sha256,
            weight_fingerprint=manifest["weight_fingerprint"],
            elapsed_seconds=time.time()-start))
        return
    evaluate = v2.normalized_objective.make_evaluator(out, lane, Q, grid, deadline,
        float(selected["price"]), native_runner=v2.native, exploratory=True)
    review_path = out / "normalization_source_review/receipt.json"
    if args.arm == "alternative":
        review = json.loads(review_path.read_text())
        review.update(household_solver_unchanged=False,
            isolated_purchase_timing_source_manifest_sha256=sha(TIMING / "manifest.json"))
        v2.write(review_path, review)
    if args.preflight_evaluator:
        v2.write(out / "completed.json", dict(status="evaluator_initialized_zero_solves",
            lifecycle_solves=0, target_fingerprint=manifest["target_fingerprint"],
            weight_fingerprint=manifest["weight_fingerprint"], source_contract=contract,
            search_evaluator="exploratory", selected_evaluator="full_native"))
        heartbeat("evaluator_initialized_zero_solves", completed_full_ge=0, objective_calls=0)
        return
    cases, cache, best, calls = [], {}, None, 0
    limit = 2 if args.smoke else MAX_CALLS
    stop = "optimizer_return"
    def objective(x):
        nonlocal best, calls
        if time.time() >= deadline - RESERVE:
            raise BudgetStop("final_native_reserve_reached")
        if calls >= limit:
            raise BudgetStop("objective_call_limit")
        if os.statvfs(out).f_bavail * os.statvfs(out).f_frsize < 20 * 1024**3:
            raise BudgetStop("free_disk_below_20GiB")
        calls += 1
        point = {k: float(lower[j] + spans[j] * x[j]) for j, k in enumerate(coordinates)}
        v2.inputs.check_point(point, bounds)
        key = tuple(point[k].hex() for k in coordinates)
        if key in cache:
            heartbeat("cache_hit", completed_full_ge=len(cases), objective_calls=calls)
            return cache[key]
        label = f"{len(cases):04d}_nm"
        heartbeat("running_full_GE", label=label, completed_full_ge=len(cases), objective_calls=calls)
        result = evaluate(label, point, deadline - RESERVE)
        if result["status"] == "passed":
            native_fits = v2.native.readtable(Path(result["report"]) / "target_fit.csv")
            fits, rr = rescore_wealth(v2, native_fits, result["residual"],
                                      Path(result["report"]))
            result["native_original_residual"] = result["residual"]
            result["residual"] = rr.tolist()
            result["target_fit"] = fits
        row = dict(label=label, parameters=point, **result)
        if result["status"] == "passed":
            rr = np.asarray(result["residual"], float)
            v2.inputs.require(rr.shape == (10,) and np.isfinite(rr).all(), "Residual contract drift")
            row["loss"] = row["objective"] = float(rr @ rr)
            row["weight_fingerprint"] = manifest["weight_fingerprint"]
            if best is None or row["loss"] < best["loss"]:
                best = row
        elif result["status"] == "inadmissible_numerical":
            row.update(objective=PENALTY, numerical_rejection=True, computed_valid_loss=False)
        elif result["status"] == "budget_exhausted":
            v2.write(out / "latest_completed.json", dict(latest=row, completed_full_ge=len(cases),
                objective_calls=calls))
            raise BudgetStop("native_evaluation_budget_exhausted")
        else:
            raise RuntimeError("Unexpected evaluator status: " + str(result["status"]))
        cases.append(row)
        cache[key] = row["objective"]
        v2.write(out / "latest_completed.json", dict(latest=row, completed_full_ge=len(cases),
            objective_calls=calls))
        v2.write(out / "best_so_far.json", dict(status="provisional_until_full_native_postcheck", best=best,
            completed_full_ge=len(cases)))
        v2.write(out / "cases.json", cases)
        heartbeat("case_completed", completed_full_ge=len(cases), objective_calls=calls,
                  best_loss=best["loss"] if best else None)
        return row["objective"]
    try:
        try:
            minimize(objective, initial[0], method="Nelder-Mead", bounds=[(0., 1.)] * 10,
                options=dict(initial_simplex=initial, maxfev=limit, maxiter=limit,
                             xatol=1e-4, fatol=1e-4, adaptive=True))
        except BudgetStop as exc:
            stop = str(exc)
        search = dict(status="provisional_search_finished", search_stop_reason=stop,
            arm=args.arm, chain=args.chain, selected=best, objective_calls=calls,
            completed_full_ge=len(cases), lifecycle_solves=sum(r.get("lifecycle_solves", 0) for r in cases),
            target_fingerprint=manifest["target_fingerprint"], weight_fingerprint=manifest["weight_fingerprint"],
            starts_file_sha256=args.starts_file_sha256,
            search_evaluator="exploratory", selected_evaluator="full_native",
            optimization_convergence_certified=False, no_auto_retry=True)
        v2.write(out / "search_completed.json", search)
        if best is None:
            v2.write(out / "completed.json", completion_receipt(search, "no_admissible_candidate"))
            return
        heartbeat("native_postcheck_running", completed_full_ge=len(cases), objective_calls=calls)
        child_out = out / "native_postcheck"
        command = [sys.executable]
        if sys.platform == "darwin":
            command.append(str(Path(__file__).with_name("local_entry.py")))
        command.extend((str(Path(__file__).resolve()), "--arm", args.arm, "--chain", str(args.chain),
                        "--out", str(child_out), "--deadline-epoch", str(deadline),
                        "--postcheck-only", "--search-receipt", str(out / "search_completed.json")))
        if args.starts_file is not None:
            command.extend(("--starts-file", str(args.starts_file.resolve()),
                            "--starts-file-sha256", args.starts_file_sha256))
        child_env = os.environ.copy()
        for key in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            child_env[key] = "1"
        remaining = min(RESERVE, deadline - time.time())
        if remaining <= 0:
            raise BudgetStop("No time remains for native postcheck child")
        child = subprocess.run(command, env=child_env, capture_output=True, text=True,
                               timeout=remaining, check=False)
        (out / "postcheck_child.stdout.log").write_text(child.stdout)
        (out / "postcheck_child.stderr.log").write_text(child.stderr)
        if child.returncode != 0:
            raise RuntimeError(f"Full native postcheck child failed with exit {child.returncode}")
        child_receipt = json.loads((child_out / "completed.json").read_text())
        v2.inputs.require(child_receipt["status"] == "full_native_postcheck_passed" and
                          child_receipt["search_receipt_sha256"] == sha(out / "search_completed.json") and
                          child_receipt["starts_file_sha256"] == args.starts_file_sha256 and
                          child_receipt["target_fingerprint"] == manifest["target_fingerprint"] and
                          child_receipt["weight_fingerprint"] == manifest["weight_fingerprint"],
                          "Native postcheck child receipt drift")
        verification = child_receipt["selected_postcheck"]
        native_loss = child_receipt["native_loss"]
        fits = child_receipt["target_fit"]
        parameters = child_receipt["parameters"]
        repeat = child_receipt["repeat"]
        report = Path(verification["report"])
        smoke_comparison = None
        if args.smoke:
            native_rows = v2.native.readtable(report / "target_fit.csv")
            selected_rows, selected_rr = rescore_wealth(v2, native_rows,
                verification["native_original_residual"], report)
            v2.inputs.require(selected_rows == fits and
                              np.allclose(selected_rr, best["residual"], rtol=0, atol=1e-10),
                              "Smoke search/full experimental-score drift")
            smoke_comparison = dict(status="search_full_new_target_exact",
                                    target_rows=14, nonwealth_rows_unchanged=True)
        status = "selected_numerically_verified"
        if abs(native_loss - best["loss"]) > 1e-8:
            status = "selected_native_passed_search_loss_differs"
        v2.write(out / "completed.json", completion_receipt(search, status,
            selected_postcheck=verification, native_loss=native_loss,
            target_fit=fits, parameters=parameters, repeat=repeat,
            smoke_fast_full_comparison=smoke_comparison,
            elapsed_seconds=time.time() - start))
        heartbeat("completed", completed_full_ge=len(cases), objective_calls=calls, native_loss=native_loss)
    except BaseException as exc:
        v2.write(out / "failure.json", dict(type=type(exc).__name__, message=str(exc),
            elapsed_seconds=time.time() - start, no_auto_retry=True))
        heartbeat("failed", completed_full_ge=len(cases), objective_calls=calls, error=str(exc))
        raise


if __name__ == "__main__":
    main()
