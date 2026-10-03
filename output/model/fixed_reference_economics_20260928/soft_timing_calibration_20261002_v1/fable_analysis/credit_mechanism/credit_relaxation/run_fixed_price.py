"""Two bounded canonical-production fixed-price lifecycle solves; no GE root.

Preparation is read-only for model code and performs no solve. The run phase is
intentionally gated by --run after separate lead review of prepared.json.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import signal
import sys
import threading
import time
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[6]
sys.path.insert(0, str(ROOT / "code/model"))
from production.inputs import load_inputs, DEFAULT_PRICE
from production.equilibrium import solve_at_price
from production.engine.diagnostics import write_diagnostics
from production.credit import bind_engine_credit

PIN = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/native_postcheck/selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz"
SAMEHOST = ROOT / "output/model/production_deployment_20261003/reference_unchanged/reference/phase_b_ge/selected_repeat/stage/solution_arrays.npz"
SAMEHOST_RECEIPT = ROOT / "output/model/production_deployment_20261003/reference_unchanged/receipt.json"
SAMEHOST_COMPARISON = ROOT / "output/model/production_deployment_20261003/compare_unchanged/comparison.json"
RECEIPT = ROOT / "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/revision1/saved_array_source_receipt.json"
PRICE = 0.7760569760205563
CORE = ("V", "c_pol", "hR_pol", "bp_pol", "tenure_probs", "fert_probs", "fert2_probs", "g", "g_beginning_distribution", "b_grid")
THREADS = ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
RUN_DEADLINE_SECONDS = 1500


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x.item() if isinstance(x, np.generic) else str(x)) + "\n")


def same(a, b) -> bool:
    try:
        return bool(np.array_equal(a, b))
    except Exception:
        return a == b


def setup():
    pin = json.loads(RECEIPT.read_text())["arrays"]["alternative"]
    if PIN.resolve() != Path(pin["path"]).resolve() or digest(PIN) != pin["sha256"]:
        raise RuntimeError("pinned chain-13 solution identity mismatch")
    if abs(DEFAULT_PRICE - PRICE) > 1e-14:
        raise RuntimeError("canonical price differs from chain-13 price")
    P0, grid0 = load_inputs()
    P1, grid1 = load_inputs(external_inputs={"phi": [1.0] * 4})
    if not np.array_equal(grid0, grid1):
        raise RuntimeError("wealth grids differ")
    differences = sorted(k for k in vars(P0) if not same(getattr(P0, k), getattr(P1, k)))
    if differences != ["phi"] or not np.array_equal(P0.phi, [0.8] * 4) or not np.array_equal(P1.phi, [1.0] * 4):
        raise RuntimeError(f"P change is not exclusively phi: {differences}")
    if set(vars(P0)) != set(vars(P1)):
        raise RuntimeError("P field inventory differs")
    C0, C1 = copy.deepcopy(P0), copy.deepcopy(P1)
    bind_engine_credit(C0, "corrected", float(C0.unsecured_credit_limit))
    bind_engine_credit(C1, "corrected", float(C1.unsecured_credit_limit))
    credit_differences = sorted(k for k in vars(C0) if not same(getattr(C0, k), getattr(C1, k)))
    if set(vars(C0)) != set(vars(C1)) or credit_differences != ["phi"]:
        raise RuntimeError(f"credit-bound P differs beyond phi: {credit_differences}")
    source_files = sorted((ROOT / "code/model/production").rglob("*.py")) + [
        ROOT / "code/model/production/reference_inputs/bundle.json",
        ROOT / "code/model/production/reference_inputs/arrays.npz", PIN]
    source = {str(path.relative_to(ROOT)): digest(path) for path in source_files}
    meta = {"api": "production.inputs.load_inputs; production.equilibrium.solve_at_price; production.engine.diagnostics.write_diagnostics",
            "price": PRICE, "reference": "post-interest chain13, selected_repeat saved arrays", "P_changed_fields": differences,
            "phi_baseline": P0.phi.tolist(), "phi_relaxed": P1.phi.tolist(),
            "grid_nodes": len(grid0), "income_states": int(P0.Nz), "source_sha256": source,
            "credit_bound_P_changed_fields": credit_differences,
            "solve_budget_seconds_each": 600, "solve_count": 2, "sequential": True,
            "outer_process_deadline_seconds": RUN_DEADLINE_SECONDS,
            "memory_budget_gib": 8, "array_archive_mib": round(PIN.stat().st_size / 1024**2, 2),
            "threads": {k: os.environ.get(k) for k in THREADS},
            "fixed_objects": "price, rent, payroll, H0, entry law, preferences, and all P fields other than phi",
            "credit_note": "uniform financed share 0.8 to 1; redundant purchase screen retained; no GE or calibration"}
    return (P0, grid0), (P1, grid1), meta


def verify_sources(expected: dict) -> None:
    for relative, expected_hash in expected.items():
        if digest(ROOT / relative) != expected_hash:
            raise RuntimeError(f"source changed: {relative}")


class Deadline(Exception):
    pass


def alarm_handler(_signum, _frame):
    raise Deadline("600-second per-solve deadline")


def heartbeat(case: str, stop: threading.Event):
    while not stop.is_set():
        write_json(HERE / "heartbeat.json", {"case": case, "time_epoch": time.time(), "status": "solving"})
        stop.wait(30)


def hard_timeout(label: str, stop: threading.Event, seconds: int):
    if not stop.wait(seconds):
        write_json(HERE / f"{label.replace(' ', '_')}_hard_timeout.json", {"error": f"hard deadline exceeded: {label}",
                    "time_epoch": time.time()})
        os._exit(124)


def save_case(name: str, outcome) -> dict:
    destination = HERE / name
    destination.mkdir(exist_ok=False)
    sol, P, sd = outcome["solution"], outcome["P"], outcome["shared"]
    arrays = {k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and v.dtype != object}
    arrays.update({"shared." + k: v for k, v in vars(sd).items() if isinstance(v, np.ndarray) and v.dtype != object})
    expected_phi = float(P.phi[0])
    if not np.allclose(arrays["shared.phi_choice"], expected_phi, atol=0, rtol=0):
        raise RuntimeError(f"{name}: financed share did not propagate to tenure choices")
    if not np.allclose(arrays["shared.phi_state"], expected_phi, atol=0, rtol=0):
        raise RuntimeError(f"{name}: financed share did not propagate to owner debt floors")
    np.savez_compressed(destination / "solution_arrays.npz", **arrays)
    write_json(destination / "executed_P.json", vars(P))
    write_diagnostics(sol, P, destination / "standard_diagnostics")
    figures = sorted((destination / "standard_diagnostics").glob("*.png"))
    if len(figures) != 17:
        raise RuntimeError(f"{name}: expected 17 standard diagnostic figures, got {len(figures)}")
    for key in ("tenure_probs", "fert_probs", "g", "g_beginning_distribution"):
        if not np.isfinite(arrays[key]).all():
            raise RuntimeError(f"{name}: nonfinite {key}")
    if np.min(arrays["tenure_probs"]) < -1e-8 or np.max(arrays["tenure_probs"]) > 1 + 1e-8:
        raise RuntimeError(f"{name}: invalid tenure probabilities")
    if np.min(arrays["fert_probs"]) < -1e-8 or np.max(arrays["fert_probs"]) > 1 + 1e-8:
        raise RuntimeError(f"{name}: invalid fertility probabilities")
    if abs(float(arrays["g"].sum()) - 1) > 1e-7:
        raise RuntimeError(f"{name}: realized mass does not sum to one")
    return {"arrays": str(destination / "solution_arrays.npz"), "array_count": len(arrays),
            "standard_figures": len(figures), "g_mass": float(arrays["g"].sum()),
            "phi_choice_values": np.unique(arrays["shared.phi_choice"]).tolist(),
            "source_array_sha256": digest(destination / "solution_arrays.npz")}


def run_case(name, P, grid, source_hashes):
    verify_sources(source_hashes)
    write_json(HERE / "start.json", {"case": name, "time_epoch": time.time(), "price": PRICE, "phi": P.phi})
    stop = threading.Event()
    worker = threading.Thread(target=heartbeat, args=(name, stop), daemon=True)
    worker.start()
    deadline_thread = threading.Thread(target=hard_timeout, args=(name, stop, 600), daemon=True)
    deadline_thread.start()
    started = time.monotonic()
    previous = signal.signal(signal.SIGALRM, alarm_handler)
    signal.setitimer(signal.ITIMER_REAL, 600)
    try:
        outcome = solve_at_price(P, grid, PRICE)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)
        stop.set(); worker.join(timeout=2)
        deadline_thread.join(timeout=2)
    elapsed = time.monotonic() - started
    verify_sources(source_hashes)
    result = save_case(name, outcome)
    result["solve_elapsed_seconds"] = elapsed
    write_json(HERE / f"{name}_completed.json", result)
    return result


def compare_pin():
    fresh = HERE / "phi_080/solution_arrays.npz"
    with np.load(PIN, allow_pickle=False) as old, np.load(fresh, allow_pickle=False) as new:
        table = {}
        for key in CORE:
            if key not in old or key not in new:
                raise RuntimeError(f"missing baseline core array {key}")
            if old[key].shape != new[key].shape:
                raise RuntimeError(f"baseline core shape mismatch {key}")
            error = float(np.max(np.abs(old[key].astype(float) - new[key].astype(float))))
            table[key] = error
            if not np.isfinite(error) or error > 1e-10:
                raise RuntimeError(f"baseline mismatch {key}: max abs {error}")
        write_json(HERE / "baseline_core_comparison.json", table)
    return table


def compare_samehost_control():
    if not SAMEHOST.is_file() or not SAMEHOST_RECEIPT.is_file() or not SAMEHOST_COMPARISON.is_file():
        raise RuntimeError("authenticated same-host control missing")
    receipt = json.loads(SAMEHOST_RECEIPT.read_text())
    comparison = json.loads(SAMEHOST_COMPARISON.read_text())
    if (receipt.get("status") != "passed" or receipt.get("backend") != "reference"
            or receipt.get("case") != "unchanged" or receipt.get("closure", {}).get("price") != PRICE
            or comparison.get("status") != "passed" or not comparison.get("same_host_runtime")
            or not comparison.get("all_arrays_equal") or not comparison.get("same_inputs")):
        raise RuntimeError("same-host reference authentication failed")
    with np.load(SAMEHOST, allow_pickle=False) as old, np.load(HERE / "phi_080/solution_arrays.npz", allow_pickle=False) as new:
        table = {}
        for key in CORE:
            if key not in old or key not in new or old[key].shape != new[key].shape:
                raise RuntimeError(f"same-host control lacks matched core {key}")
            equal = bool(np.array_equal(old[key], new[key]))
            table[key] = equal
            if not equal:
                raise RuntimeError(f"same-host control mismatch: {key}")
    control = {"status": "passed_exact", "samehost_core_equal": table,
               "samehost_reference_npz": str(SAMEHOST), "samehost_reference_sha256": digest(SAMEHOST),
               "samehost_receipt_sha256": digest(SAMEHOST_RECEIPT),
               "samehost_comparison_sha256": digest(SAMEHOST_COMPARISON),
               "fresh_baseline_sha256": digest(HERE / "phi_080/solution_arrays.npz"),
               "historical_cluster_pin": str(PIN), "historical_cluster_pin_sha256": digest(PIN),
               "historical_global_1e_minus10_gate": "failed; preserved in failure.json",
               "basis": "same-host authenticated execution control; no tolerance relaxation"}
    write_json(HERE / "samehost_control.json", control)
    return control


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--continue-after-samehost-control", action="store_true")
    args = parser.parse_args()
    if sum((args.prepare, args.run, args.continue_after_samehost_control)) != 1:
        parser.error("choose exactly one mode")
    if not args.prepare and not all(os.environ.get(k) == "1" for k in THREADS):
        raise RuntimeError("all Numba/BLAS/OpenMP thread limits must equal one")
    baseline, relaxed, meta = setup()
    HERE.mkdir(exist_ok=True)
    if args.prepare:
        write_json(HERE / "prepared.json", meta)
        print(json.dumps(meta, indent=2)); return
    if not (HERE / "prepared.json").is_file() or json.loads((HERE / "prepared.json").read_text())["source_sha256"] != meta["source_sha256"]:
        raise RuntimeError("reviewed preparation missing or source identity changed")
    if args.continue_after_samehost_control:
        if (HERE / "phi_100").exists():
            raise RuntimeError("relaxed case directory already exists")
        if not (HERE / "phi_080_completed.json").is_file():
            raise RuntimeError("fresh baseline receipt absent")
        verify_sources(meta["source_sha256"])
        control = compare_samehost_control()
        if (HERE / "start.json").is_file() and not (HERE / "baseline_start.json").is_file():
            (HERE / "baseline_start.json").write_bytes((HERE / "start.json").read_bytes())
        write_json(HERE / "continuation_start.json", {"time_epoch": time.time(), "control": control,
                    "action": "one phi=1 fixed-price lifecycle solve"})
        overall_stop = threading.Event()
        overall_watch = threading.Thread(target=hard_timeout, args=("continuation", overall_stop, RUN_DEADLINE_SECONDS), daemon=True)
        overall_watch.start()
        try:
            second = run_case("phi_100", *relaxed, meta["source_sha256"])
            with (HERE / "phi_080/executed_P.json").open() as stream: executed0 = json.load(stream)
            with (HERE / "phi_100/executed_P.json").open() as stream: executed1 = json.load(stream)
            input_fields = set(vars(baseline[0]))
            executed_differences = sorted(k for k in input_fields if executed0[k] != executed1[k])
            if executed_differences != ["phi"]:
                raise RuntimeError(f"executed primitive P differs beyond phi: {executed_differences}")
            write_json(HERE / "completed.json", {"status": "passed", "control": control,
                        "baseline": json.loads((HERE / "phi_080_completed.json").read_text()),
                        "relaxed": second, "executed_primitive_P_changed_fields": executed_differences,
                        "finished_epoch": time.time()})
        except Exception as exc:
            write_json(HERE / "continuation_failure.json", {"error": str(exc), "traceback": traceback.format_exc(),
                        "time_epoch": time.time()})
            raise
        finally:
            overall_stop.set(); overall_watch.join(timeout=2)
        return
    overall_stop = threading.Event()
    overall_watch = threading.Thread(target=hard_timeout, args=("whole run", overall_stop, RUN_DEADLINE_SECONDS), daemon=True)
    overall_watch.start()
    try:
        first = run_case("phi_080", *baseline, meta["source_sha256"])
        comparison = compare_pin()  # stops before relaxation on any mismatch
        second = run_case("phi_100", *relaxed, meta["source_sha256"])
        with (HERE / "phi_080/executed_P.json").open() as stream: executed0 = json.load(stream)
        with (HERE / "phi_100/executed_P.json").open() as stream: executed1 = json.load(stream)
        executed_differences = sorted(k for k in executed0 if executed0[k] != executed1[k])
        if set(executed0) != set(executed1) or executed_differences != ["phi"]:
            raise RuntimeError(f"executed P differs beyond phi: {executed_differences}")
        write_json(HERE / "completed.json", {"status": "passed", "baseline": first,
                    "baseline_core_max_abs": comparison, "relaxed": second,
                    "executed_P_changed_fields": executed_differences,
                    "finished_epoch": time.time()})
    except Exception as exc:
        write_json(HERE / "failure.json", {"error": str(exc), "traceback": traceback.format_exc(),
                    "time_epoch": time.time()})
        raise
    finally:
        overall_stop.set(); overall_watch.join(timeout=2)


if __name__ == "__main__":
    main()
