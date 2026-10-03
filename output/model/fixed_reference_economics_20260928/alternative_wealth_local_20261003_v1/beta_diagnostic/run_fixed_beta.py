"""One fixed-parameter, alternative-timing GE diagnostic around verified chain 2."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[5]
PACKET = Path(__file__).resolve().parent.parent
EXP = ROOT / "code/model/experiments/alternative_wealth_calibration"
sys.path.insert(0, str(EXP))
import calibrate as exp


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--beta-annual", type=float, choices=(0.95, 0.94), required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--deadline-epoch", type=float, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise RuntimeError("Refusing existing diagnostic output")
    if args.deadline_epoch - time.time() > 420 or args.deadline_epoch - time.time() < 60:
        raise RuntimeError("Seven-minute per-solve budget required")
    args.out.mkdir(parents=True)
    out = args.out.resolve()
    baseline_file = PACKET / "overnight_two_wave/chain_02/completed.json"
    baseline = json.loads(baseline_file.read_text())
    if baseline["status"] != "selected_numerically_verified" or baseline["repeat"]["experimental_target_fit_exact"] is not True:
        raise RuntimeError("Baseline chain 2 not fully verified")
    v2, timing_driver, manifest, selected = exp.checked_inputs("alternative")
    if (baseline["target_fingerprint"] != manifest["target_fingerprint"] or
            baseline["weight_fingerprint"] != manifest["weight_fingerprint"]):
        raise RuntimeError("Baseline and diagnostic target contract differ")
    point = dict(baseline["selected"]["parameters"])
    original_point = dict(point)
    point["beta_annual"] = args.beta_annual
    if set(point) != set(original_point) or sum(point[k] != original_point[k] for k in point) != 1:
        raise RuntimeError("Exactly one parameter must change")
    lane = "floor_s0"
    _, bounds, _ = v2.inputs.seed_and_bounds(lane)
    bounds = {k: tuple(value) for k, value in bounds.items()}
    bounds["h_P"] = (.1, 2.6)
    bounds["psi_child"] = tuple(v2.CONFIG["psi_bounds"])
    coordinates = tuple(v2.inputs.parameters(lane)) + ("psi_child",)
    if set(point) != set(coordinates):
        raise RuntimeError("Ten-coordinate contract drift")
    v2.inputs.check_point(point, bounds)
    v2.inputs.LANES[lane].update(seed=point, bounds=bounds,
                                 free_coordinates=list(coordinates))
    P, grid = v2.inputs.proposal(lane)
    P, entry = v2.inputs.entry(P, grid, "nonnegative_mean")
    bound_point = v2.inputs.bind(P, point, bounds, "floor")
    if (abs(float(bound_point.beta) - args.beta_annual ** float(bound_point.period_years)) > 1e-12 or
            abs(float(bound_point.psi_child) - point["psi_child"]) > 1e-12):
        raise RuntimeError("Annual-to-period discounting or free psi drift")
    if (P.N_target != 1. or P.R_gross <= 1. or
            not np.allclose(P.phi, .8, rtol=0., atol=1e-12) or
            not P.native_purchase_income or not P.native_due_stayer_credit or
            P.joint_nested_choice):
        raise RuntimeError("Non-beta economic input drift")
    source = dict(driver_sha256=sha(__file__), baseline_completion_sha256=sha(baseline_file),
                  selected_source_sha256=json.loads(exp.SELECTION.read_text())["source_sha256"],
                  timing_manifest_sha256=sha(exp.TIMING / "manifest.json"),
                  normalized_source_pins_sha256=sha(exp.V2 / "source_pins.json"))
    v2.write(out / "request.json", dict(status="fixed_beta_ge_requested",
        beta_annual=args.beta_annual, baseline_beta_annual=original_point["beta_annual"],
        baseline_native_loss=baseline["native_loss"],
        baseline_price=baseline["selected_postcheck"]["price"],
        baseline_parameters=original_point, diagnostic_parameters=point,
        changed_coordinates=["beta_annual"], bounds=bounds, entry=entry,
        arm="alternative", target_fingerprint=manifest["target_fingerprint"],
        weight_fingerprint=manifest["weight_fingerprint"], source=source,
        deadline_epoch=args.deadline_epoch, one_core=True,
        objective="illustrative_fixed_parameter_GE_sensitivity_not_calibration"))
    try:
        Q = v2.native.utility_checks(P, grid, lane, out)
        exp.install_timing_observer(v2, timing_driver, out, P)
        evaluate = v2.normalized_objective.make_evaluator(
            out, lane, Q, grid, args.deadline_epoch,
            float(baseline["selected_postcheck"]["price"]),
            native_runner=v2.native, exploratory=False)
        result = evaluate("fixed_beta", point, args.deadline_epoch)
        if result["status"] != "passed":
            raise RuntimeError("Native GE did not pass: " + str(result))
        report = Path(result["report"])
        native_rows = v2.native.readtable(report / "target_fit.csv")
        rows, residual = exp.rescore_wealth(v2, native_rows, result["residual"], report)
        params = v2.native.readtable(report / "parameters.csv")
        if len(rows) != 14 or len(params) != 31:
            raise RuntimeError("Target or parameter row count drift")
        if v2.native.target_identity(rows) != manifest["target_contract"]:
            raise RuntimeError("New target contract drift")
        repeat_report = report.parent / "selected_repeat_final"
        repeat = v2.native.compare_repeated(report, repeat_report)
        repeated_rows, repeated_residual = exp.rescore_wealth(v2,
            v2.native.readtable(repeat_report / "target_fit.csv"),
            v2.native.residual(v2.native.readtable(repeat_report / "target_fit.csv")),
            repeat_report)
        if rows != repeated_rows or not np.array_equal(residual, repeated_residual):
            raise RuntimeError("New-contract target table differs on native repeat")
        if len(list((report / "standard_diagnostics").glob("*.png"))) != 17:
            raise RuntimeError("Standard diagnostic plot count drift")
        loss = sum(float(row["loss_contribution"] or 0) for row in rows)
        if abs(loss - float(residual @ residual)) > 1e-8:
            raise RuntimeError("Experimental objective arithmetic drift")
        v2.write(out / "completed.json", dict(status="fixed_beta_native_ge_verified",
            beta_annual=args.beta_annual, baseline_beta_annual=original_point["beta_annual"],
            parameters_changed=["beta_annual"], baseline_completion_sha256=sha(baseline_file),
            source=source, target_fingerprint=manifest["target_fingerprint"],
            weight_fingerprint=manifest["weight_fingerprint"],
            native_result=result, experimental_loss=loss, target_fit=rows,
            original_diagnostic_target_fit=native_rows, parameters=params,
            repeat=repeat, new_contract_exact_repeat=True,
            selected_report=str(report), standard_plot_count=17,
            optimizer_convergence_not_applicable=True))
    except BaseException as exc:
        v2.write(out / "failure.json", dict(type=type(exc).__name__, message=str(exc),
                                         no_fallback=True))
        raise


if __name__ == "__main__":
    main()
