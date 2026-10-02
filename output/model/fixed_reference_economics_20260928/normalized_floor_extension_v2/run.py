"""Fixed-coordinate normalized floor sensitivity; no optimizer or model edits."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
V2 = HERE.parent / "normalized_calibration_v2"
sys.path.insert(0, str(V2))
import run_psi as v2  # noqa: E402

POINTS = (2.3, 2.4, 2.5, 2.6)
EXTENDED_UPPER = 2.6
LIMIT_SECONDS = 1800
LANE = "floor_s0"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_contract():
    incumbent = json.loads((HERE / "incumbent.json").read_text())
    manifest = json.loads((HERE / "manifest.json").read_text())
    for rel, digest in manifest["sha256"].items():
        assert sha(ROOT / rel) == digest, f"Source drift: {rel}"
    assert sha(HERE / "run.py") == manifest["runner_sha256"], "V2 runner drift"
    v2.native.verify_sources()
    for rel, digest in json.loads((V2 / "source_pins.json").read_text()).items():
        assert sha(ROOT / rel) == digest, f"V2 source drift: {rel}"
    assert incumbent["postcheck_status"] == "selected_numerically_verified"
    assert incumbent["chain"] == 2 and incumbent["case"] == "0064_nm"
    assert sha(HERE / "incumbent.json") == manifest["incumbent_sha256"]
    assert v2.inputs.canonical(v2.CONFIG["base_target_contract"]) == manifest["target_fingerprint"]
    assert v2.weight_fingerprint({}) == manifest["weight_fingerprint"]
    assert v2.CONFIG["profiles"]["base_control"] == {}, "Original weights required"
    seed, old_bounds, _ = v2.inputs.seed_and_bounds(LANE)
    old_bounds = {key: tuple(value) for key, value in old_bounds.items()}
    old_bounds["psi_child"] = tuple(v2.CONFIG["psi_bounds"])
    assert old_bounds["h_P"] == (0.1, 2.3), "Original floor bound changed"
    assert set(incumbent["parameters"]) == set(v2.inputs.parameters(LANE)) | {"psi_child"}
    assert 0.1 <= incumbent["parameters"]["h_P"] <= 2.3
    assert 0 < incumbent["selected_price"] < 80
    assert len(incumbent["postcheck_target_fit"]) == 14
    assert len(incumbent["postcheck_parameters"]) == 31
    return incumbent, old_bounds, manifest


def configure(seed, bounds):
    v2.inputs.LANES[LANE].update(seed=dict(seed), bounds=dict(bounds),
                                 free_coordinates=list(seed))


def rows(report, filename):
    with (Path(report) / filename).open(newline="") as fh:
        return list(csv.DictReader(fh))


def inspect(result, point, bounds):
    assert result["status"] == "passed", result
    report = Path(result["report"])
    target = rows(report, "target_fit.csv")
    params = rows(report, "parameters.csv")
    assert len(target) == 14 and len(params) == 31
    assert v2.native.target_identity(target) == v2.CONFIG["base_target_contract"]
    assert result["population"] == 1.0
    assert float(next(r["estimate"] for r in params if r["parameter"] == "h_P")) == point["h_P"]
    assert float(next(r["upper"] for r in params if r["parameter"] == "h_P")) == bounds["h_P"][1]
    assert len(list((report / "standard_diagnostics").glob("*.png"))) == 17
    repeat = report.parent / "selected_repeat_final"
    repeated = v2.native.compare_repeated(report, repeat)
    assert repeated["status"] == "exact_full_ge_repeat_passed"
    return {"status": "passed", "report": str(report), "price": result["price"],
            "H0_derived": result["H0_derived"], "base_loss": float(sum(float(r["loss_contribution"] or 0) for r in target)),
            "target_fit": target, "parameters": params, "repeat": repeated,
            "lifecycle_solves": result["lifecycle_solves"]}


def solve(out, point, bounds, deadline, price, label):
    out.mkdir(parents=True)
    configure(point if label == "old_bound" else BASE["parameters"], bounds)
    P, grid = v2.inputs.proposal(LANE)
    P, _ = v2.inputs.entry(P, grid, "nonnegative_mean")
    Q = v2.native.utility_checks(P, grid, LANE, out)
    evaluate = v2.normalized_objective.make_evaluator(out, LANE, Q, grid, deadline,
                                                       price, native_runner=v2.native)
    return inspect(evaluate(label, point, deadline), point, bounds)


def replay_gate(out, incumbent, old_bounds, deadline):
    point = incumbent["parameters"]
    extended = dict(old_bounds, h_P=(0.1, EXTENDED_UPPER))
    first = run_gate_arm(out / "old", "old_bound", deadline, mock_arm=MOCK_ARM)
    assert first["target_fit"] == incumbent["postcheck_target_fit"], "Old-bound target replay differs from authenticated postcheck"
    assert first["parameters"] == incumbent["postcheck_parameters"], "Old-bound parameters differ from authenticated postcheck"
    assert first["price"] == incumbent["selected_price"], "Old-bound price differs from authenticated postcheck"
    v2.write(out / "latest_completed.json", {"stage": "old_bound", "result": first})
    second = run_gate_arm(out / "extended", "extended_bound", deadline, mock_arm=MOCK_ARM)
    compare_gate_arms(first, second, incumbent)
    result = {"status": "incumbent_replay_passed", "old": first, "extended": second,
              "old_h_P_bounds": old_bounds["h_P"], "extended_h_P_bounds": extended["h_P"],
              "all_other_coordinates_fixed": True}
    v2.write(out / "completed.json", result)
    return result


def compare_gate_arms(first, second, incumbent):
    assert first["target_fit"] == second["target_fit"], "Economic target drift across bound override"
    for a, b in zip(first["parameters"], second["parameters"]):
        assert a["parameter"] == b["parameter"] and a["estimate"] == b["estimate"], "Parameter drift"
        if a["parameter"] != "h_P":
            assert a == b, f"Unrelated parameter metadata drift: {a['parameter']}"
    assert first["price"] == second["price"] and first["H0_derived"] == second["H0_derived"]
    assert abs(first["base_loss"] - incumbent["base_loss"]) < 1e-10


def run_gate_arm(out, arm, deadline, mock_arm=False):
    remaining = deadline - time.time()
    assert remaining > 0, "Shared replay-gate deadline expired"
    out.mkdir(parents=True, exist_ok=True)
    result_path = out / "arm_result.json"
    command = [sys.executable, str(HERE / "run.py"), "--mode", "arm", "--arm", arm,
               "--deadline", str(deadline), "--out", str(out)]
    if mock_arm:
        command.append("--mock-arm")
    try:
        subprocess.run(command, check=True, timeout=remaining,
                       env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
    except subprocess.TimeoutExpired as exc:
        raise TimeoutError(f"{arm} arm exceeded shared 1800-second gate deadline") from exc
    assert result_path.exists(), f"Fresh-process {arm} arm produced no result"
    result = json.loads(result_path.read_text())
    assert result["arm"] == arm and result["deadline_epoch"] == deadline
    return result["result"]


def run_arm(out, arm, deadline, incumbent, old_bounds, mock_arm):
    remaining = deadline - time.time()
    if remaining <= 0:
        raise TimeoutError("Shared replay-gate deadline expired before arm start")
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("Shared 30-minute gate deadline")))
    signal.setitimer(signal.ITIMER_REAL, remaining)
    try:
        if mock_arm:
            result = mock_arm_result(arm, incumbent, old_bounds)
        else:
            bounds = old_bounds if arm == "old_bound" else dict(old_bounds, h_P=(0.1, EXTENDED_UPPER))
            result = solve(out / "native", incumbent["parameters"], bounds, deadline,
                           incumbent["selected_price"], arm)
        v2.write(out / "arm_result.json", {"arm": arm, "deadline_epoch": deadline,
                                             "result": result, "fresh_process": True})
    except BaseException as exc:
        v2.write(out / "failure.json", {"type": type(exc).__name__, "message": str(exc),
                                         "arm": arm, "deadline_epoch": deadline})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def mock_arm_result(arm, incumbent, old_bounds):
    target = copy.deepcopy(incumbent["postcheck_target_fit"])
    params = copy.deepcopy(incumbent["postcheck_parameters"])
    if arm == "extended_bound":
        for row in params:
            if row["parameter"] == "h_P":
                row["upper"] = str(EXTENDED_UPPER)
                row["near_bound"] = str(EXTENDED_UPPER - float(row["estimate"]) <= .01 * (EXTENDED_UPPER - .1))
    return {"target_fit": target, "parameters": params, "base_loss": incumbent["base_loss"],
            "price": incumbent["selected_price"],
            "H0_derived": float(next(r["estimate"] for r in params if r["parameter"] == "H0")),
            "report": "mock", "repeat": {"status": "exact_full_ge_repeat_passed"},
            "lifecycle_solves": 0}


def mock(out, incumbent, old_bounds):
    points = [dict(incumbent["parameters"], h_P=h) for h in POINTS]
    assert all(set(x) == set(incumbent["parameters"]) for x in points)
    assert all(all(x[k] == incumbent["parameters"][k] for k in x if k != "h_P") for x in points)
    assert old_bounds["h_P"] == (0.1, 2.3)
    assert all(0.1 <= x["h_P"] <= EXTENDED_UPPER for x in points)
    global MOCK_ARM
    MOCK_ARM = True
    deadline = time.time() + LIMIT_SECONDS
    replay = replay_gate(out / "mock_gate", incumbent, old_bounds, deadline)
    assert replay["old"]["lifecycle_solves"] == replay["extended"]["lifecycle_solves"] == 0
    assert len(replay["old"]["target_fit"]) == 14 and len(replay["old"]["parameters"]) == 31
    for index, point in enumerate(points):
        assert 0.1 <= point["h_P"] <= EXTENDED_UPPER
        assert all(point[k] == incumbent["parameters"][k] for k in point if k != "h_P")
    # Keep the original target and non-floor parameter drift rejection gates.
    first = replay["old"]
    bad_target = copy.deepcopy(replay["extended"])
    bad_target["target_fit"][0]["model"] = "drift"
    try: compare_gate_arms(first, bad_target, incumbent)
    except AssertionError as exc: assert "Economic target drift" in str(exc)
    else: raise AssertionError("Target perturbation was not rejected")
    bad_parameter = copy.deepcopy(replay["extended"])
    bad_parameter["parameters"][0]["estimate"] = "drift"
    try: compare_gate_arms(first, bad_parameter, incumbent)
    except AssertionError as exc: assert "Parameter drift" in str(exc)
    else: raise AssertionError("Parameter perturbation was not rejected")
    MOCK_ARM = False
    v2.write(out / "mock.json", {"status": "passed_zero_model_solves", "points": points,
                                   "replay_calls": 2, "fresh_child_processes": 2,
                                   "production_calls": 4,
                                   "negative_target_and_parameter_gates": True,
                                   "lifecycle_solves": 0})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("mock", "gate", "point", "arm"), required=True)
    parser.add_argument("--index", type=int, choices=range(4))
    parser.add_argument("--arm", choices=("old_bound", "extended_bound"))
    parser.add_argument("--deadline", type=float)
    parser.add_argument("--mock-arm", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.mode != "arm":
        assert not args.out.exists(), "Refusing existing output"
        args.out.mkdir(parents=True)
    incumbent, old_bounds, manifest = load_contract()
    global BASE
    BASE = incumbent
    global MOCK_ARM
    MOCK_ARM = False
    if args.mode == "arm":
        assert args.arm and args.deadline and args.out.exists()
        run_arm(args.out, args.arm, args.deadline, incumbent, old_bounds, args.mock_arm)
        return
    if args.mode == "mock":
        mock(args.out, incumbent, old_bounds)
        return
    if args.mode == "point":
        assert args.index is not None
        gate = Path(os.environ["FLOOR_GATE_RECEIPT"])
        gate_data = json.loads(gate.read_text())
        assert gate_data["status"] == "incumbent_replay_passed"
        assert gate_data["source_manifest_sha256"] == sha(HERE / "manifest.json")
    start = time.time()
    deadline = start + LIMIT_SECONDS
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("30-minute worker limit")))
    signal.setitimer(signal.ITIMER_REAL, LIMIT_SECONDS)
    v2.write(args.out / "start.json", {"mode": args.mode, "index": args.index,
                                       "started_epoch": start, "deadline_epoch": deadline,
                                       "source_manifest_sha256": sha(HERE / "manifest.json"),
                                       "target_fingerprint": manifest["target_fingerprint"],
                                       "weight_fingerprint": manifest["weight_fingerprint"]})
    try:
        if args.mode == "gate":
            result = replay_gate(args.out, incumbent, old_bounds, deadline)
            result["source_manifest_sha256"] = sha(HERE / "manifest.json")
            v2.write(args.out / "completed.json", result)
        else:
            point = dict(incumbent["parameters"], h_P=POINTS[args.index])
            bounds = dict(old_bounds, h_P=(0.1, EXTENDED_UPPER))
            v2.write(args.out / "latest.json", {"status": "running", "h_P": point["h_P"],
                                                "started_epoch": start})
            result = solve(args.out / "native", point, bounds, deadline,
                           incumbent["selected_price"], f"hP_{args.index}")
            result.update(status="fixed_coordinate_verified", h_P=point["h_P"],
                          fixed_coordinates={k: v for k, v in point.items() if k != "h_P"},
                          target_fingerprint=manifest["target_fingerprint"],
                          weight_fingerprint=manifest["weight_fingerprint"],
                          source_manifest_sha256=sha(HERE / "manifest.json"),
                          experimental_not_adopted=True)
            v2.write(args.out / "latest_completed.json", result)
            v2.write(args.out / "best_so_far.json", result)
            v2.write(args.out / "completed.json", result)
    except BaseException as exc:
        v2.write(args.out / "failure.json", {"type": type(exc).__name__, "message": str(exc),
                                             "elapsed_seconds": time.time() - start})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    main()
