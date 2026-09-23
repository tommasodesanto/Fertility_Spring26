#!/usr/bin/env python3
"""Persistent, bounded adaptive worker for the four utility comparison cells."""
from __future__ import annotations

import argparse
import copy
import fcntl
import gzip
import hashlib
import json
import math
import os
import pickle
import random
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

SCHEMA = "e5f_utility_overnight_manifest_v1"
SHARED = ("H0", "beta_annual", "chi", "first_birth_fixed_cost", "kappa_fert",
          "kappa_fert_continuation", "theta0", "theta1")
POSITIVE = {"H0", "chi", "kappa_fert", "kappa_fert_continuation", "theta1", "h_P"}
WIDTH = {"beta_annual": .005, "H0": .5, "chi": .15, "first_birth_fixed_cost": .08,
         "kappa_fert": .18, "kappa_fert_continuation": .18, "theta0": .04,
         "theta1": .03, "h_P": .12, "delta_alpha_jump": .012, "delta_alpha": .012}
NUMERICAL_PHRASES = ("market clearing gate", "numerical gate", "dead-mass gate",
                     "occupied value", "budget excess", "not converged", "out-of-grid mass")


class ContractError(ValueError):
    pass


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def canon(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def deadline(value):
    return datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()


def cell_task(manifest, task_id):
    if not 1 <= task_id <= 40 or manifest.get("workers") != {"B_floor": 18, "B_shares": 18, "D_floor": 2, "D_shares": 2}:
        raise ContractError("40-worker allocation changed")
    if task_id <= 18:
        return "B_floor", task_id
    if task_id <= 36:
        return "B_shares", task_id - 18
    if task_id <= 38:
        return "D_floor", task_id - 36
    return "D_shares", task_id - 38


def validate_manifest(manifest, path):
    if manifest.get("schema") != SCHEMA or manifest.get("status") != "ready":
        raise ContractError("manifest schema/status mismatch")
    if manifest.get("entry_rule") != "fixed_reference_marginal" or manifest.get("source_file_count") != 641:
        raise ContractError("entry/source contract mismatch")
    if manifest.get("wealth_grid") != {"nodes": 160, "upper": 3000}:
        raise ContractError("wealth grid contract mismatch")
    if set(manifest.get("templates", {})) != {"B_floor", "B_shares", "D_floor", "D_shares"}:
        raise ContractError("four template cells required")
    if not (deadline(manifest["production_stop_utc"]) < deadline(manifest["verification_stop_utc"]) < deadline(manifest["report_deadline_utc"])):
        raise ContractError("absolute cutoff ordering changed")
    for cell, item in manifest["templates"].items():
        p = Path(item["path"])
        if not p.is_file() or sha(p) != item["sha256"]:
            raise ContractError(f"{cell} template missing/hash mismatch")
        plan = read(p)
        if plan["objective_canonical_sha256"] != item["objective_canonical_sha256"] or plan["target_system_sha256"] != manifest["target_system_sha256"]:
            raise ContractError(f"{cell} objective/target fingerprint mismatch")
        if plan["entry_specification"]["rule"] != "fixed_reference_marginal":
            raise ContractError(f"{cell} has wrong entry rule")
        if plan["income_specification"]["constructor_arguments"]["n_persistent"] != item["income_grid"]["persistent"] or plan["income_specification"]["constructor_arguments"]["n_iid"] != item["income_grid"]["iid"]:
            raise ContractError(f"{cell} income grid mismatch")
        expected = set(SHARED) | ({"h_P"} if cell.endswith("floor") else {"delta_alpha_jump", "delta_alpha"})
        if set(plan["structural_parameters"]) != expected or set(plan["parameter_bounds"]) != expected:
            raise ContractError(f"{cell} free-coordinate count changed")
        if plan["parameter_bounds"]["beta_annual"][1] > .99 or (cell.endswith("floor") and plan["parameter_bounds"]["h_P"][1] > 2.3):
            raise ContractError("upper bound changed")
    if path and sha(path) != manifest.get("self_sha256", sha(path)):
        raise ContractError("manifest self pin changed")


def clipped(point, bounds):
    return {name: min(float(bounds[name][1]), max(float(bounds[name][0]), float(value)))
            for name, value in point.items()}


def propose(center, bounds, seed, number, *, scale=.65):
    # Common random innovations in the eight shared dimensions for paired
    # floor/share workers. Their cells can subsequently learn separate centers.
    common = random.Random(seed + 104729 * number)
    extra = random.Random(seed + 104729 * number + 2901827)
    point = {}
    for name in (*SHARED, *(n for n in center if n not in SHARED)):
        u = common.uniform(-1, 1) if name in SHARED else extra.uniform(-1, 1)
        width = WIDTH[name] * scale
        point[name] = center[name] * math.exp(width * u) if name in POSITIVE else center[name] + width * u
    if number % 5 == 0:
        # A genuine coordinate probe around the learned center.
        name = tuple(center)[(number // 5 - 1) % len(center)]
        point = dict(center)
        point[name] = center[name] + WIDTH[name] * (.5 if number % 10 else -.5)
    return clipped(point, bounds)


def make_case(template, parameters, case_id, path, seconds, proposal_metadata=None):
    plan = copy.deepcopy(template)
    bounds = plan["parameter_bounds"]
    if set(parameters) != set(bounds):
        raise ContractError("proposal parameter dimension mismatch")
    for name, value in parameters.items():
        lo, hi = map(float, bounds[name])
        if not math.isfinite(value) or not lo <= value <= hi:
            raise ContractError(f"proposal outside bound: {name}")
    plan["structural_parameters"] = dict(parameters)
    plan["starting_structural_parameters"] = dict(parameters)
    plan["comparison_draw"] = {"shared_parameters": {k: parameters[k] for k in SHARED},
        "preference_parameters": {k: v for k, v in parameters.items() if k not in SHARED}}
    plan["case_id"] = case_id
    plan["overnight_proposal"] = proposal_metadata or {}
    plan["cases"] = [{"id": case_id, "arm": "literature_income_purchase",
        "native_seconds": min(6900, int(seconds - 300)), "wrapper_seconds": min(7100, int(seconds - 100)),
        "seconds": int(seconds), "repetitions": 1}]
    write(path, plan)
    return sha(path)


def kill_group(proc):
    if proc.poll() is not None:
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()


def run_case(manifest, cell, stage, case_id, parameters, output, cutoff, timeout,
             proposal_metadata=None):
    case_dir = output / stage / cell / case_id
    case_dir.mkdir(parents=True, exist_ok=False)
    template = read(manifest["templates"][cell]["path"])
    plan_path = case_dir / "plan.json"
    planned = min(timeout, int(cutoff - time.time()))
    if planned < 600:
        raise ContractError("less than 600 seconds remain for a native case")
    plan_sha = make_case(template, parameters, case_id, plan_path, planned, proposal_metadata)
    status_path = case_dir / "status.json"
    started = time.time()
    status = {"status": "started", "cell": cell, "stage": stage, "case_id": case_id,
        "started": started, "updated": started, "parameters": parameters,
        "proposal_metadata": proposal_metadata or {},
        "plan_sha256": plan_sha, "objective_canonical_sha256": template["objective_canonical_sha256"],
        "target_system_sha256": manifest["target_system_sha256"]}
    write(status_path, status)
    command = [sys.executable, template["adapter_path"], "--plan", str(plan_path),
               "--output", str(case_dir / "result"), "--arm", "literature_income_purchase",
               "--repetitions", "1"]
    with (case_dir / "worker.log").open("w") as log:
        proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        while proc.poll() is None:
            now = time.time()
            if now >= min(cutoff, started + timeout):
                kill_group(proc)
                status.update(status="incomplete_timeout", completed=time.time(),
                              elapsed_seconds=time.time()-started)
                write(status_path, status)
                return status
            status["updated"] = now
            write(case_dir / "heartbeat.json", {"status": "running", "updated": now,
                  "case_id": case_id, "pid": proc.pid})
            time.sleep(min(30, max(.2, min(cutoff, started + timeout)-now)))
    if proc.returncode:
        tail = (case_dir / "worker.log").read_text(errors="replace")[-12000:]
        native_failure = case_dir / "result/evaluation/raw/failure.json"
        native_infeasible = native_failure.is_file() and read(native_failure).get("error_type") == "InfeasibleThetaError"
        if native_infeasible or any(phrase in tail.lower() for phrase in NUMERICAL_PHRASES):
            status.update(status="numerical_gate_rejection", completed=time.time(),
                          elapsed_seconds=time.time()-started, error_excerpt=tail[-2000:])
            write(status_path, status)
            return status
        status.update(status="fatal_contract_or_code_error", completed=time.time(),
                      error_excerpt=tail[-2000:])
        write(status_path, status)
        raise ContractError(f"{case_id}: adapter exit {proc.returncode}; inspect {case_dir / 'worker.log'}")
    result = case_dir / "result" / "evaluation"
    score_path = result / "scored_repetition_01" / "score.json"
    summary_path = result / "summary.json"
    if not score_path.is_file() or not summary_path.is_file():
        raise ContractError(f"{case_id}: missing score/summary")
    score, summary = read(score_path), read(summary_path)
    expected_count = 17 if cell.endswith("floor") else 19
    if summary.get("status") != "verified_scored_candidate" or summary.get("objective_canonical_sha256") != template["objective_canonical_sha256"]:
        raise ContractError(f"{case_id}: wrapper quality/objective mismatch")
    if len(score.get("target_fit", [])) != 13 or len(score.get("parameters", [])) != expected_count:
        raise ContractError(f"{case_id}: incomplete full fit/parameter table")
    graphs = summary.get("original_graphs", [])
    if len(graphs) != 17 or any(not (result / "raw" / row["path"]).is_file() for row in graphs):
        raise ContractError(f"{case_id}: missing standard diagnostic graphs")
    entry_path = result / "entry_wealth.json"
    if not entry_path.is_file():
        raise ContractError(f"{case_id}: missing inherited entry-wealth receipt")
    entry = read(entry_path)
    marginal = entry.get("candidate_wealth_marginal")
    reference = entry.get("reference_wealth_marginal")
    if ("rank coupled" not in entry.get("rule", "") or not isinstance(marginal, list)
            or not isinstance(reference, list) or len(marginal) != 160
            or len(reference) != 160 or max(abs(float(a)-float(b)) for a, b in zip(marginal, reference)) > 2e-14
            or abs(float(marginal[0])) > 1e-12 or abs(float(marginal[-1])) > 1e-12
            or float(entry.get("candidate_wealth_mean", 0)) <= 0):
        raise ContractError(f"{case_id}: inherited wealth marginal/rank or endpoint-censor check failed")
    loss = float(score["loss"])
    if not math.isfinite(loss):
        raise ContractError(f"{case_id}: nonfinite loss")
    status.update(status="completed", completed=time.time(), elapsed_seconds=time.time()-started,
        objective=loss, score_sha256=sha(score_path), summary_sha256=sha(summary_path),
        target_rows=13, parameter_rows=expected_count, graph_count=17,
        entry_wealth_sha256=sha(entry_path), entry_wealth_mean=entry["candidate_wealth_mean"],
        entry_wealth_marginal_max_gap=entry["wealth_marginal_max_gap"],
        score_path=str(score_path), summary_path=str(summary_path))
    write(status_path, status)
    return status


def update_cell(output, cell, result):
    root = output / "rollup" / cell
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        history = root / "cases.jsonl"
        with history.open("a") as stream:
            stream.write(json.dumps(result, sort_keys=True, allow_nan=False) + "\n")
        write(root / "latest_status.json", result)
        if result["status"] == "completed":
            write(root / "latest_completed_case.json", result)
        best_path = root / "best_so_far.json"
        best = read(best_path) if best_path.is_file() else None
        if result["status"] == "completed" and (best is None or result["objective"] < best["objective"]):
            write(best_path, result)
            best = result
        fcntl.flock(lock, fcntl.LOCK_UN)
    return best


def smoke(manifest, task_id, output):
    cell = ("B_floor", "B_shares", "D_floor", "D_shares")[task_id-1]
    cutoff = min(deadline(manifest["production_stop_utc"]), time.time()+15000)
    template = manifest["templates"][cell]
    seed = template["seed_parameters"]
    for number, point in enumerate((seed, propose(seed, template["bounds"], 20260923, 1, scale=.2)), 1):
        if time.time()+max(210, manifest["smoke_case_seconds"]*.15) >= cutoff:
            raise ContractError("smoke budget insufficient for exact two-case loop")
        result = run_case(manifest, cell, "smoke", f"smoke_{number:02d}", point,
                          output, cutoff, manifest["smoke_case_seconds"],
                          {"design": "exact_loop_anchor" if number == 1 else "exact_loop_probe",
                           "paired_seed_id": 1, "scale": 0 if number == 1 else .2})
        update_cell(output, cell, result)
        if result["status"] != "completed":
            raise ContractError(f"{cell} smoke case {number} did not pass: {result['status']}")
    write(output / "smoke" / cell / "complete.json", {"status": "verified_two_case_loop", "cell": cell})


def production(manifest, task_id, output):
    cell, pair = cell_task(manifest, task_id)
    smoke_done = output / "smoke" / cell / "complete.json"
    if not smoke_done.is_file() or read(smoke_done)["status"] != "verified_two_case_loop":
        raise ContractError(f"{cell} exact-loop smoke not verified")
    template = manifest["templates"][cell]
    bound = template["bounds"]
    seed = template["seed_parameters"]
    paired_seed = 20260923 + pair * 7907
    cutoff = deadline(manifest["production_stop_utc"])
    local_best = None
    history = []
    for number in range(1, 19):
        # 18 independent paired seeds, then each persistent worker follows its
        # own valid best. A fresh draw is made only after reading that score.
        best_path = output / "rollup" / cell / "best_so_far.json"
        global_best = read(best_path) if best_path.is_file() else None
        center = local_best["parameters"] if local_best else (global_best["parameters"] if global_best else seed)
        scale_group = ("local" if pair == 1 else "medium") if cell[0] == "D" else (
            "local" if pair <= 6 else "medium" if pair <= 12 else "broad")
        base_scale = {"local": .8, "medium": 1.6, "broad": 3.2}[scale_group]
        scale = max(.3 * base_scale, base_scale * (.95 ** (number-1)))
        if number % 6 == 0:
            scale = base_scale
        point = seed if number == 1 and pair == 1 else propose(center, bound, paired_seed, number, scale=scale)
        if any(canon(point) == canon(row["parameters"]) for row in history):
            point = propose(center, bound, paired_seed, number+37, scale=scale)
        remaining = cutoff - time.time()
        initial_estimate = manifest["observed_seconds_per_objective"]["B_15x1" if cell[0] == "B" else "D_7x3"]
        if remaining < max(600, 1.15 * (history[-1]["elapsed_seconds"] if history else initial_estimate)):
            break
        case_id = f"worker{task_id:02d}_proposal{number:02d}"
        proposal_metadata = {"paired_seed_id": pair, "scale_group": scale_group,
            "scale": scale, "proposal_number": number, "center": center,
            "design": "fixed_anchor" if number == 1 and pair == 1 else ("coordinate_probe" if number % 5 == 0 else "multivariate_adaptive"),
            "shared_random_innovations": True,
            "parameter_vectors_need_not_match_across_utility_arms": True}
        result = run_case(manifest, cell, "production", case_id, point,
                          output, cutoff, manifest["case_seconds"], proposal_metadata)
        history.append(result)
        update_cell(output, cell, result)
        write(output / "production" / "workers" / f"worker_{task_id:02d}.json",
              {"cell": cell, "pair_seed_id": pair, "task_id": task_id,
               "started": len(history), "completed": sum(r["status"] == "completed" for r in history),
               "rejected": sum(r["status"] == "numerical_gate_rejection" for r in history),
               "incomplete": sum(r["status"] == "incomplete_timeout" for r in history),
               "updated": time.time(), "last_case": case_id})
        if result["status"] == "completed" and (local_best is None or result["objective"] < local_best["objective"]):
            local_best = result
        if result["status"] == "incomplete_timeout":
            break
    return {"task_id": task_id, "cell": cell, "attempted": len(history),
            "completed": sum(r["status"] == "completed" for r in history)}


def selected_for(output, cell):
    best_path = output / "rollup" / cell / "best_so_far.json"
    if not best_path.is_file():
        raise ContractError(f"{cell} has no valid smoke/production selection")
    selected = read(best_path)
    if selected["status"] != "completed":
        raise ContractError("selected point is not completed")
    return selected


def verification(manifest, task_id, output):
    cell = ("B_floor", "B_shares", "D_floor", "D_shares")[(task_id-1)//2]
    number = 1 + (task_id-1) % 2
    cutoff = deadline(manifest["verification_stop_utc"])
    selected = selected_for(output, cell)
    if time.time()+max(600, 1.05 * selected.get("elapsed_seconds", 900)) >= cutoff:
        raise ContractError(f"{cell} selected repeat cannot fit verification deadline")
    result = run_case(manifest, cell, "verification", f"selected_repeat_{number:02d}",
                      selected["parameters"], output, cutoff, manifest["case_seconds"],
                      {"design": "exact_selected_repeat", "repeat": number,
                       "selected_case_id": selected["case_id"]})
    if result["status"] != "completed":
        raise ContractError(f"{cell} selected repeat {number} did not complete")
    return result


def verified_target_fit(score, checkpoint_sha, cell):
    rows = []
    for target in score["target_fit"]:
        if target.get("model_checkpoint_sha256") != checkpoint_sha:
            raise ContractError(f"{cell}: target row checkpoint provenance mismatch")
        rows.append({k: v for k, v in target.items()
                     if k not in ("model_checkpoint_sha256", "model_source_path")})
    return rows


def scientific_signature(manifest, cell, row):
    """Compare saved price, V, g, moments and normalized psi exactly."""
    import numpy as np
    source = Path(read(manifest["templates"][cell]["path"])["source_root"])
    sys.path[:0] = [str(source / "code/model/tools"), str(source / "code/model")]
    evaluation = Path(row["score_path"]).parents[1]
    raw = evaluation / "raw/repetition_01"
    summary = read(raw / "summary.json")
    with gzip.open(raw / "initial_state.pkl.gz", "rb") as stream:
        packet = pickle.load(stream)
    score = read(row["score_path"])
    checkpoint_sha = sha(raw / "initial_state.pkl.gz")
    if checkpoint_sha != summary["checkpoint_sha256"]:
        raise ContractError(f"{cell}: saved checkpoint hash differs from native summary")
    target_fit = verified_target_fit(score, checkpoint_sha, cell)
    def array_hash(array):
        data = np.ascontiguousarray(array)
        return {"shape": list(data.shape), "dtype": str(data.dtype),
                "sha256": hashlib.sha256(data.view(np.uint8)).hexdigest()}
    signature = {
        "price": summary["price"],
        "V": array_hash(packet["evaluation"].policy.V),
        "g": array_hash(packet["evaluation"].g_current),
        "moments": summary["legacy_stationary_moments"],
        "early_measurement": summary["early_measurement"],
        "fiscal": summary["fiscal"],
        "household_budget": summary["household_budget"],
        "psi": float(packet["parameters"].psi_child),
        "normalization": {k: v for k, v in score.get("normalization", {}).items()
                          if k != "stationary_solve_seconds"},
        "target_fit": target_fit,
        "parameters": score["parameters"],
        "loss": row["objective"],
    }
    return signature


def finalize(manifest, task_id, output):
    cell = ("B_floor", "B_shares", "D_floor", "D_shares")[task_id-1]
    selected = selected_for(output, cell)
    repeats = []
    for number in (1, 2):
        path = output / "verification" / cell / f"selected_repeat_{number:02d}" / "status.json"
        if not path.is_file() or read(path).get("status") != "completed":
            raise ContractError(f"{cell} repeat {number} missing/incomplete")
        repeats.append(read(path))
    signature = scientific_signature(manifest, cell, selected)
    comparison = [scientific_signature(manifest, cell, row) for row in repeats]
    exact = all(item == signature for item in comparison)
    receipt = {"cell": cell, "status": "verified_exact_twice" if exact else "requires_review",
        "selected": selected, "repeats": repeats,
        "signature_sha256": canon(signature),
        "repeat_signature_sha256": [canon(item) for item in comparison],
        "fields_compared": ["price", "V", "g", "moments", "early_measurement", "fiscal", "household_budget", "psi", "normalization", "target_fit", "parameters", "loss"],
        "completed": time.time()}
    write(output / "verification" / cell / "receipt.json", receipt)
    if not exact:
        raise ContractError(f"{cell} independent scientific-array comparison differs")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--stage", choices=("smoke", "production", "verification", "finalize"), required=True)
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    if sha(args.manifest) != args.manifest_sha256:
        raise ContractError("manifest byte hash mismatch")
    manifest = read(args.manifest)
    validate_manifest(manifest, args.manifest)
    if sha(Path(__file__)) != manifest["runner"]["sha256"]:
        raise ContractError("running controller does not match manifest code pin")
    adapter_path = Path(manifest["templates"]["B_floor"]["path"])
    adapter_pin = manifest["adapter_source_change"]["bundled_sha256"]
    if sha(read(adapter_path)["adapter_path"]) != adapter_pin:
        raise ContractError("bundled preference adapter hash changed")
    args.output_root.mkdir(parents=True, exist_ok=True)
    try:
        if args.stage == "smoke":
            if not 1 <= args.task_id <= 4: raise ContractError("smoke task index out of range")
            result = smoke(manifest, args.task_id, args.output_root)
        elif args.stage == "production":
            result = production(manifest, args.task_id, args.output_root)
        elif args.stage == "verification":
            if not 1 <= args.task_id <= 8: raise ContractError("verification task index out of range")
            result = verification(manifest, args.task_id, args.output_root)
        else:
            if not 1 <= args.task_id <= 4: raise ContractError("finalize task index out of range")
            result = finalize(manifest, args.task_id, args.output_root)
    except BaseException as exc:
        write(args.output_root / "failures" / f"{args.stage}_{args.task_id:02d}.json",
              {"status": "fatal_or_incomplete_stage", "stage": args.stage,
               "task_id": args.task_id, "error_type": type(exc).__name__,
               "error": str(exc), "updated": time.time()})
        raise
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
