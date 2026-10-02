"""Fixed-coordinate normalized floor sensitivity; no optimizer or model edits."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import os
import signal
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
    first = solve(out / "old", point, old_bounds, deadline, incumbent["selected_price"], "old_bound")
    assert first["target_fit"] == incumbent["postcheck_target_fit"], "Old-bound target replay differs from authenticated postcheck"
    assert first["parameters"] == incumbent["postcheck_parameters"], "Old-bound parameters differ from authenticated postcheck"
    assert first["price"] == incumbent["selected_price"], "Old-bound price differs from authenticated postcheck"
    v2.write(out / "latest_completed.json", {"stage": "old_bound", "result": first})
    second = solve(out / "extended", point, extended, deadline, incumbent["selected_price"], "extended_bound")
    assert first["target_fit"] == second["target_fit"], "Economic target drift across bound override"
    for a, b in zip(first["parameters"], second["parameters"]):
        assert a["parameter"] == b["parameter"] and a["estimate"] == b["estimate"], "Parameter drift"
        if a["parameter"] != "h_P":
            assert a == b, f"Unrelated parameter metadata drift: {a['parameter']}"
    assert first["price"] == second["price"] and first["H0_derived"] == second["H0_derived"]
    assert abs(first["base_loss"] - incumbent["base_loss"]) < 1e-10
    result = {"status": "incumbent_replay_passed", "old": first, "extended": second,
              "old_h_P_bounds": old_bounds["h_P"], "extended_h_P_bounds": extended["h_P"],
              "all_other_coordinates_fixed": True}
    v2.write(out / "completed.json", result)
    return result


def mock(out, incumbent, old_bounds):
    points = [dict(incumbent["parameters"], h_P=h) for h in POINTS]
    assert all(set(x) == set(incumbent["parameters"]) for x in points)
    assert all(all(x[k] == incumbent["parameters"][k] for k in x if k != "h_P") for x in points)
    assert old_bounds["h_P"] == (0.1, 2.3)
    assert all(0.1 <= x["h_P"] <= EXTENDED_UPPER for x in points)
    calls = []
    original = globals()["solve"]
    def stub(directory, point, bounds, deadline, price, label):
        assert point == incumbent["parameters"] or point in points
        assert price == incumbent["selected_price"] and deadline > time.time()
        assert bounds["h_P"] in ((0.1, 2.3), (0.1, 2.6))
        calls.append((label, point["h_P"], bounds["h_P"]))
        target = copy.deepcopy(incumbent["postcheck_target_fit"])
        params = copy.deepcopy(incumbent["postcheck_parameters"])
        for row in params:
            if row["parameter"] == "h_P":
                row["estimate"] = str(point["h_P"])
                row["upper"] = str(bounds["h_P"][1])
                row["near_bound"] = str(bounds["h_P"][1] - point["h_P"] <= .01 * (bounds["h_P"][1] - bounds["h_P"][0]))
        # Preserve the original exact strings for the authenticated old-bound replay.
        if label == "old_bound": params = copy.deepcopy(incumbent["postcheck_parameters"])
        return {"target_fit": target, "parameters": params,
                "base_loss": incumbent["base_loss"], "price": incumbent["selected_price"],
                "H0_derived": float(next(r["estimate"] for r in params if r["parameter"] == "H0")),
                "report": "mock", "repeat": {"status": "exact_full_ge_repeat_passed"},
                "lifecycle_solves": 0}
    try:
        globals()["solve"] = stub
        replay_gate(out / "mock_gate", incumbent, old_bounds, time.time() + 1800)
        extended = dict(old_bounds, h_P=(0.1, 2.6))
        for index, point in enumerate(points):
            stub(out / f"mock_point{index}", point, extended, time.time() + 1800,
                 incumbent["selected_price"], f"hP_{index}")
        assert [x[0] for x in calls] == ["old_bound", "extended_bound", "hP_0", "hP_1", "hP_2", "hP_3"]
        assert [x[1] for x in calls[2:]] == list(POINTS)
        # A changed target row must block the actual-replay gate.
        def bad_target(*args):
            result = stub(*args)
            if args[-1] == "extended_bound": result["target_fit"][0]["model"] = "drift"
            return result
        globals()["solve"] = bad_target
        try: replay_gate(out / "negative_target", incumbent, old_bounds, time.time() + 1800)
        except AssertionError as exc: assert "Economic target drift" in str(exc)
        else: raise AssertionError("Target perturbation was not rejected")
        # A changed non-floor estimate must also block it.
        def bad_parameter(*args):
            result = stub(*args)
            if args[-1] == "extended_bound": result["parameters"][0]["estimate"] = "drift"
            return result
        globals()["solve"] = bad_parameter
        try: replay_gate(out / "negative_parameter", incumbent, old_bounds, time.time() + 1800)
        except AssertionError as exc: assert "Parameter drift" in str(exc)
        else: raise AssertionError("Parameter perturbation was not rejected")
    finally:
        globals()["solve"] = original
    v2.write(out / "mock.json", {"status": "passed_zero_model_solves", "points": points,
                                   "replay_calls": 2, "production_calls": 4,
                                   "negative_target_and_parameter_gates": True,
                                   "lifecycle_solves": 0})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("mock", "gate", "point"), required=True)
    parser.add_argument("--index", type=int, choices=range(4))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "Refusing existing output"
    args.out.mkdir(parents=True)
    incumbent, old_bounds, manifest = load_contract()
    global BASE
    BASE = incumbent
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
