"""Matched fixed-price Estate-A cap-one credit diagnostic; no GE or calibration."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
OUT = ROOT / "output/model/experiments/birth_count_choice/credit_at_binary_winner_v1"
PRICE = 0.7811670615311468
PER_SOLVE_SECONDS = 600
THREADS = ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def readrows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    files = sorted((HERE / "model").glob("*.py")) + sorted((HERE / "model/engine").glob("*.py"))
    files += [HERE / "run_cap2_at_binary_winner.py", Path(__file__).resolve()]
    return {str(path.relative_to(ROOT)): sha(path) for path in files}


def verify_sources(expected):
    if source_hashes() != expected:
        raise RuntimeError("Birth-count experiment source changed after preparation")


def setup():
    from model.inputs import load_inputs
    from model.estate_contract import experiment_flags, apply_experiment_flags, contract
    from model.calibration import effective_input_fingerprint
    from run_cap2_at_binary_winner import inputs_and_receipt, SOURCE, SOURCE_SHA256

    source, point, h0, _, grid, target_pin, weight_pin = inputs_and_receipt()
    if source["selected"]["price"] != PRICE:
        raise RuntimeError("Binary-winner price drift")
    P0, g0 = load_inputs(parameters=point, external_inputs={"H0": [h0]})
    P1, g1 = load_inputs(parameters=point, external_inputs={"H0": [h0], "phi": [1.] * 4})
    for P in (P0, P1):
        apply_experiment_flags(P, experiment_flags(1))
    if not np.array_equal(g0, g1) or not np.array_equal(g0, grid):
        raise RuntimeError("Wealth grid drift")
    if set(vars(P0)) != set(vars(P1)):
        raise RuntimeError("Input field inventory drift")
    differences = sorted(k for k in vars(P0) if not np.array_equal(
        np.asarray(getattr(P0, k)), np.asarray(getattr(P1, k))))
    if differences != ["phi"] or not np.array_equal(P0.phi, [.8] * 4) or not np.array_equal(P1.phi, [1.] * 4):
        raise RuntimeError(f"Expected exactly uniform phi change: {differences}")
    if set(point) != set(contract()[3]) or len(point) != 10:
        raise RuntimeError("Ten-coordinate winner contract drift")
    from model.engine.shared import get_phi_choice_tensor, get_phi_state_matrix
    propagation = {}
    for label, P, expected in (("phi_080", P0, .8), ("phi_100", P1, 1.)):
        choice, state = get_phi_choice_tensor(P), get_phi_state_matrix(P)
        if not np.array_equal(choice, np.full_like(choice, expected)) or not np.array_equal(state, np.full_like(state, expected)):
            raise RuntimeError(f"{label}: phi does not propagate uniformly to choice and state")
        propagation[label] = dict(phi_choice_shape=list(choice.shape), phi_state_shape=list(state.shape), value=expected)
    preparation = dict(status="zero_solve_prepared", source_search_sha256=SOURCE_SHA256,
        source_selected_input_fingerprint=source["selected"]["input_fingerprint"],
        source_selected_new_contract_loss=source["selected"]["loss"],
        target_fingerprint=target_pin, weight_fingerprint=weight_pin,
        selected_root_target_fit_sha256=sha(SOURCE / "selected_root/target_fit_new_contract.csv"),
        selected_root_parameters_sha256=sha(SOURCE / "selected_root/parameters_estate_a.csv"),
        price=PRICE, birth_cap=1, fixed_H0=h0, point=point,
        input_changed_fields=differences,
        baseline_input_fingerprint=effective_input_fingerprint(P0, g0),
        relaxed_input_fingerprint=effective_input_fingerprint(P1, g1),
        source_sha256=source_hashes(),
        phi_propagation=propagation, per_solve_deadline_seconds=PER_SOLVE_SECONDS,
        outer_deadline_seconds=1500, sequential_solve_count=2,
        interpretation="fixed-price partial equilibrium; no GE root or recalibration")
    return (P0, g0), (P1, g1), preparation


def _timeout(_signum, _frame):
    raise TimeoutError("600-second fixed-price solve limit")


def solve_one(label, P, grid, preparation):
    from model.equilibrium import solve_at_price
    from model.engine.diagnostics import write_diagnostics
    from model.reporting import build_context
    from model import native_phase_b

    verify_sources(preparation["source_sha256"])
    arm = OUT / label
    arm.mkdir(exist_ok=False)
    (arm / "reporting").mkdir()
    write(arm / "start.json", dict(label=label, epoch=time.time(), price=PRICE,
                                    phi=P.phi.tolist(), preparation_sha256=sha(OUT / "prepared.json")))
    # Install authenticated observers before the solver; this is a zero-solve operation.
    context = build_context(P, grid, arm / "reporting", price_start=PRICE,
                            deadline=time.time() + PER_SOLVE_SECONDS,
                            max_lifecycle=32, closure="fixed_h0")
    started = time.monotonic()
    prior = signal.signal(signal.SIGALRM, _timeout)
    prior_timer = signal.setitimer(signal.ITIMER_REAL, PER_SOLVE_SECONDS)
    try:
        outcome = solve_at_price(P, grid, PRICE)
    finally:
        signal.setitimer(signal.ITIMER_REAL, *prior_timer)
        signal.signal(signal.SIGALRM, prior)
    elapsed = time.monotonic() - started
    sol, Q, sd = outcome["solution"], outcome["P"], outcome["shared"]
    for field, expected in (("phi_choice", float(P.phi[0])), ("phi_state", float(P.phi[0]))):
        actual = np.asarray(getattr(sd, field))
        if not np.array_equal(actual, np.full_like(actual, expected)):
            raise RuntimeError(f"{label}: executed shared.{field} differs from phi")
    arrays = {k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and not v.dtype.hasobject}
    arrays.update({"shared." + k: v for k, v in vars(sd).items()
                   if isinstance(v, np.ndarray) and not v.dtype.hasobject})
    np.savez_compressed(arm / "solution_arrays.npz", **arrays)
    write(arm / "solve_completed.json", dict(status="fixed_price_solve_complete", label=label,
        solve_elapsed_seconds=elapsed, price=PRICE, phi=Q.phi.tolist(),
        array_count=len(arrays), arrays_sha256=sha(arm / "solution_arrays.npz"),
        shared_phi_choice_values=np.unique(sd.phi_choice).tolist(),
        shared_phi_state_values=np.unique(sd.phi_state).tolist(),
        g_mass=float(np.asarray(sol.g).sum()),
        beginning_mass=float(np.asarray(sol.g_beginning_distribution).sum())))
    write_diagnostics(sol, Q, arm / "standard_diagnostics")
    if len(list((arm / "standard_diagnostics").glob("*.png"))) != 17:
        raise RuntimeError(f"{label}: standard diagnostic plot count differs from 17")
    live = dict(P=Q, b_grid=np.asarray(grid), sd=sd, sol=sol,
                price=np.asarray([PRICE]), case_deadline_epoch=time.time() + 300)
    # Baseline uses the unchanged final native observer/gates to authenticate
    # its 14 moments against the selected GE point. Relaxed arm is diagnostic:
    # final=True would impose a GE renewal root on a fixed-price counterfactual.
    context["out"] = arm / "reporting"
    try:
        observed = native_phase_b.observe_price(context, live, label,
                                                 final=(label == "phi_080"))
        write(arm / "native_observation.json", observed)
    except Exception as exc:
        write(arm / "reporting_failure.json", dict(error_type=type(exc).__name__,
            message=str(exc), traceback=traceback.format_exc(),
            solver_output_preserved=True, no_gate_change=True))
        if label == "phi_080":
            raise
    return arm


def verify_baseline(arm):
    from run_cap2_at_binary_winner import SOURCE
    from model.estate_contract import rescore_report
    native = arm / "reporting/phase_b_ge/phi_080"
    rows, _, parameters, receipt = rescore_report(native)
    saved = readrows(SOURCE / "selected_root/target_fit_new_contract.csv")
    current = readrows(native / "target_fit_new_contract.csv")
    source = json.loads((SOURCE / "provenance/search_completed.json").read_text())
    if len(rows) != 14 or len(parameters) != 31 or len(current) != 14:
        raise RuntimeError("Baseline PE full native target/parameter report missing")
    if [{k: r[k] for k in ("moment", "target", "weight", "role")} for r in current] != [
            {k: r[k] for k in ("moment", "target", "weight", "role")} for r in saved]:
        raise RuntimeError("Baseline target/weight contract differs")
    errors = {r["moment"]: abs(float(r["model"]) - float(s["model"]))
              for r, s in zip(current, saved)}
    residual = np.asarray([np.sqrt(float(r["weight"])) * float(r["gap"])
                           for r in current if r["role"] == "scored"])
    selected_residual = np.asarray(source["selected"]["residual"])
    if residual.shape != (10,) or selected_residual.shape != (10,):
        raise RuntimeError("Baseline ten scored residuals missing")
    residual_gap = float(np.max(np.abs(residual - selected_residual)))
    # The existing selected-point residual agreement gate in cluster_calibrate.py.
    if residual_gap > 1e-10:
        raise RuntimeError(f"Baseline PE scored residuals differ from verified winner: {residual_gap}; raw moments {errors}")
    result = dict(status="passed", target_rows=14, parameter_rows=31,
                  max_model_moment_abs_gap=max(errors.values()),
                  max_scored_residual_abs_gap=residual_gap,
                  existing_selected_residual_gate=1e-10,
                  source_comparison_scope="full target moments; no local winner arrays available",
                  new_contract_loss=receipt["loss"])
    write(OUT / "baseline_fit_gate.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--run", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(HERE))
    baseline, relaxed, preparation = setup()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.prepare:
        write(OUT / "prepared.json", preparation)
        print(json.dumps(preparation, indent=2, sort_keys=True))
        return
    if not all(os.environ.get(name) == "1" for name in THREADS):
        raise RuntimeError("All six Numba/BLAS/OpenMP thread limits must be one")
    if not (OUT / "prepared.json").is_file() or json.loads((OUT / "prepared.json").read_text()) != preparation:
        raise RuntimeError("Prepared exact inputs/source missing or drifted")
    if any((OUT / name).exists() for name in ("phi_080", "phi_100")):
        raise RuntimeError("Refusing duplicate fixed-price solve")
    try:
        first = solve_one("phi_080", *baseline, preparation)
        gate = verify_baseline(first)
        second = solve_one("phi_100", *relaxed, preparation)
        write(OUT / "completed.json", dict(status="two_fixed_price_solves_complete",
              phi_080=str(first), phi_100=str(second), baseline_fit_gate=gate,
              relaxed_reporting_status="observation_only" if not (second / "reporting_failure.json").exists() else "failed_preserved",
              completed_epoch=time.time()))
    except Exception as exc:
        write(OUT / "failure.json", dict(error_type=type(exc).__name__, message=str(exc),
              traceback=traceback.format_exc(), no_auto_retry=True))
        raise


if __name__ == "__main__":
    main()
