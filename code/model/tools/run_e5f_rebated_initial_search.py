"""Bounded search after a verified equal-rebate initial smoke."""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import copy
import csv
import importlib.util
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import threading
import time

HELPER_SHA256 = "d9aa97b890442d45971ec622b4b41687da10ffa26eedaf0198f352c7e6ecb790"
JOINT_SHA256 = "9eee3bca39f2a98f4a58cf18196d695b6a9db3e93ef18f8eaa2bf4dbb1243bbb"
SCORED_SHA256 = "9da35b55466d74dc10a6a85ff91dabe887a807aab68417eabf63a9d7a42ca967"
DEFAULT_WORKERS = 6
MAXIMUM_WORKERS = 18
MAXIMUM_COORDINATES = 18
MAXIMUM_JOINT = 18
MAXIMUM_WALL_SECONDS = 10800
FINAL_RESERVE_SECONDS = 2100


def read(path): return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def file_sha(path):
    import hashlib
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_source_pins(helper_path, helper_sha256, joint_path, joint_sha256,
                        scored_path, scored_sha256):
    expected = ((helper_path, helper_sha256, HELPER_SHA256, "helper"),
                (joint_path, joint_sha256, JOINT_SHA256, "joint adapter"),
                (scored_path, scored_sha256, SCORED_SHA256, "scored invoker"))
    for path, supplied, approved, label in expected:
        if supplied != approved or file_sha(path) != supplied:
            raise ValueError(f"Frozen {label} source pin changed")


def load_helper(path, expected_sha256=HELPER_SHA256):
    path = Path(path).resolve()
    spec = importlib.util.spec_from_file_location("rebated_initial_helper", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    if expected_sha256 != HELPER_SHA256 or module.sha(path) != expected_sha256:
        raise ValueError("Frozen rebated initial helper fingerprint changed")
    return module


def load_controller(packet):
    path = Path(packet["plan_path"]).parent / "run_capped_beta.py"
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location("saved_rebated_search_controller", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    if module.sha(path) != packet["plan"]["controller_sha256"]:
        raise ValueError("Saved proposal controller fingerprint changed")
    return module


def smoke_seed(helper, smoke, packet, restrictions, joint_sha256=JOINT_SHA256):
    smoke = Path(smoke).resolve()
    summary = read(smoke / "summary.json")
    if (summary.get("status") != "verified_rebated_initial_smoke"
            or summary.get("method") != "joint_price_fertility_rebate"):
        raise ValueError("Search requires a verified equal-rebate smoke")
    launch = read(smoke / "launch_contract.json")
    if launch.get("helper_sha256") != HELPER_SHA256 or launch.get("joint_sha256") != joint_sha256:
        raise ValueError("Smoke numerical adapter pins differ from this search")
    accounting = read(smoke / "case/evaluation/rebate_accounting.json")
    if accounting.get("status") != "verified" or len(accounting.get("receipts", [])) != 1:
        raise ValueError("Smoke rebate accounting receipt is incomplete")
    score = read(smoke / "case/evaluation/scored_repetition_01/score.json")
    proposal = copy.deepcopy(launch["proposal"])
    controller = load_controller(packet)
    controller.validate_score(score, proposal, restrictions)
    if score["loss"] != summary["loss"]:
        raise ValueError("Smoke score and summary loss differ")
    return dict(case_id=proposal["case_id"], status="verified", loss=score["loss"],
                output=str(smoke / "case/evaluation"), proposal=proposal,
                score=score, accounting=accounting["receipts"])


def failure_status(detail):
    text = str(detail).lower()
    hard_markers = ("fingerprint", "source pin", "source manifest", "target changed",
                    "objective source", "objective canonical", "mixed score contract",
                    "structural contract", "scientific contract", "saved full twelve-row",
                    "saved first-birth rooms target", "complete nine-coordinate")
    if any(marker in text for marker in hard_markers):
        return "failed"
    if "mass" in text or "feasibility" in text:
        return "rejected_mass_gate"
    if "equilibrium" in text or "market" in text or "converg" in text:
        return "rejected_equilibrium"
    return "rejected_numerical"


def bounded_failure_evidence(folder, maximum_characters=12000):
    folder = Path(folder)
    paths = [folder / "candidate.log", folder / "case/summary.json",
             folder / "case/candidate_result.json", folder / "case/case/branch_failure.json",
             folder / "case/case/evaluation/failure.json",
             folder / "case/case/evaluation/raw/failure.json",
             folder / "case/case/evaluation/initial_solve.log"]
    pieces = []
    remaining = maximum_characters
    for path in paths:
        if remaining <= 0 or not path.is_file(): continue
        text = path.read_text(errors="replace")
        excerpt = text[-min(len(text), remaining, 3000):]
        pieces.append(f"{path.name}: {excerpt}"); remaining -= len(excerpt)
    return "\n".join(pieces) or "candidate exited without a readable receipt"


def error_signature(result):
    text = str(result.get("failure_detail", "")).lower()
    text = re.sub(r"/[^\s:]+", "<path>", text)
    text = re.sub(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?", "<n>", text)
    return re.sub(r"\s+", " ", text).strip()[:500]


def scored_command(scored_path, helper_path, helper_sha256, joint_path, joint_sha256,
                   template, item_path, output):
    return [sys.executable, "-B", str(scored_path),
            "--helper", str(helper_path), "--helper-sha256", helper_sha256,
            "--joint", str(joint_path), "--joint-sha256", joint_sha256,
            "--template", str(template), "--output", str(output),
            "--proposal", str(item_path)]


def candidate_mode(helper_path, helper_sha256, joint_path, joint_sha256,
                   scored_path, scored_sha256, template, item_path, output):
    require_source_pins(helper_path, helper_sha256, joint_path, joint_sha256,
                        scored_path, scored_sha256)
    command = scored_command(scored_path, helper_path, helper_sha256, joint_path,
                             joint_sha256, template, item_path, output)
    completed = subprocess.run(command)
    receipt = Path(output) / "candidate_result.json"
    if not receipt.exists():
        raise RuntimeError(f"Joint scored candidate exited {completed.returncode} without a receipt")
    result = read(receipt)
    if result["status"] != "verified":
        detail = str(result.get("failure_detail", "")) + "\n" + bounded_failure_evidence(output)
        result["failure_detail"] = detail[-12000:]
        result["status"] = failure_status(result["failure_detail"])
    write(Path(output).parent / "result.json", result)
    if result["status"] == "failed": raise SystemExit(2)


def run_child(helper_path, helper_sha256, joint_path, joint_sha256,
              scored_path, scored_sha256, template, item, folder, deadline):
    folder.mkdir(parents=True, exist_ok=False)
    write(folder / "proposal.json", item)
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        return dict(case_id=item["case_id"], status="timeout", proposal=item,
                    output=str(folder), failure_detail="branch wall clock exhausted")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMBA_NUM_THREADS="1", NUMBA_DISABLE_JIT="0",
               PYTHONOPTIMIZE="0", MPLCONFIGDIR=str(folder / "mpl"))
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "candidate",
               "--helper", str(helper_path), "--helper-sha256", helper_sha256,
               "--joint", str(joint_path), "--joint-sha256", joint_sha256,
               "--scored", str(scored_path), "--scored-sha256", scored_sha256,
               "--template", str(template),
               "--item", str(folder / "proposal.json"), "--output", str(folder / "case")]
    try:
        with (folder / "candidate.log").open("w") as log:
            child = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                     start_new_session=True)
            try:
                code = child.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL); child.wait(); raise
    except subprocess.TimeoutExpired:
        return dict(case_id=item["case_id"], status="timeout", proposal=item,
                    output=str(folder), failure_detail="candidate wall clock exhausted")
    if not (folder / "result.json").exists():
        detail = bounded_failure_evidence(folder)
        return dict(case_id=item["case_id"], status=failure_status(detail), proposal=item,
                    output=str(folder), failure_detail=detail)
    return read(folder / "result.json")


def bounded_batch(items, worker, workers, completed):
    """Keep at most workers live and stop submitting after a hard branch error."""
    pending = iter(items); active = {}; results = []; hard = False; failures = {}
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in range(min(workers, len(items))):
            item = next(pending); active[pool.submit(worker, item)] = item
        while active:
            done, _ = cf.wait(active, return_when=cf.FIRST_COMPLETED)
            for future in done:
                active.pop(future)
                result = future.result(); results.append(result); completed(result)
                hard = hard or result["status"] == "failed"
                if result["status"].startswith("rejected_"):
                    signature = error_signature(result)
                    failures[signature] = failures.get(signature, 0) + 1
                    if signature and failures[signature] >= 3:
                        hard = True
                if not hard:
                    try:
                        item = next(pending); active[pool.submit(worker, item)] = item
                    except StopIteration:
                        pass
    return results, hard


def write_csv(path, rows):
    if not rows: return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(rows)


def aggregate_tables(out, records):
    fits, parameters = [], []
    for record in records:
        if record["status"] != "verified": continue
        score = record["score"]
        fits.extend(dict(case_id=record["case_id"], case_loss=record["loss"], **row)
                    for row in score["target_fit"])
        parameters.extend(dict(case_id=record["case_id"], case_loss=record["loss"], **row)
                          for row in score["parameters"])
    write_csv(Path(out) / "all_target_fits.csv", fits)
    write_csv(Path(out) / "all_parameters.csv", parameters)


def run_search(helper_path, helper_sha256, joint_path, joint_sha256,
               scored_path, scored_sha256, template, smoke, output, workers, wall_seconds):
    if not 1 <= workers <= MAXIMUM_WORKERS or not 1 <= wall_seconds <= MAXIMUM_WALL_SECONDS:
        raise ValueError("Worker or wall-clock budget exceeds the approved bound")
    require_source_pins(helper_path, helper_sha256, joint_path, joint_sha256,
                        scored_path, scored_sha256)
    started = time.monotonic(); deadline = started + wall_seconds
    search_deadline = deadline - FINAL_RESERVE_SECONDS
    helper = load_helper(helper_path, helper_sha256); packet = helper.saved_packet(template)
    restrictions = helper.validate_scientific_contract(packet); controller = load_controller(packet)
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    seed = smoke_seed(helper, smoke, packet, restrictions, joint_sha256)
    records = [seed]; best = seed; lock = threading.Lock(); hard_error = False
    write(out / "search_contract.json", dict(status="running", helper_sha256=helper_sha256,
        joint_sha256=joint_sha256, scored_invoker_sha256=scored_sha256,
        source_root=str(packet["source_root"]), smoke=str(Path(smoke).resolve()), workers=workers,
        maximum_pool=MAXIMUM_WORKERS, wall_seconds=wall_seconds,
        maximum_coordinate_cases=MAXIMUM_COORDINATES, maximum_joint_cases=MAXIMUM_JOINT,
        search_seconds=wall_seconds-FINAL_RESERVE_SECONDS,
        final_verification_reserve_seconds=FINAL_RESERVE_SECONDS,
        final_exact_repetitions=2, scored_moments=12, free_parameters=9,
        beta_upper=0.99, normalization=2.1, rooms_target=helper.ROOMS_TARGET))

    def persist(latest):
        nonlocal best
        if latest["status"] == "verified" and latest["loss"] < best["loss"]: best = latest
        compact = lambda x: {k: v for k, v in x.items() if k != "score"}
        write(out / "latest_completed.json", compact(latest))
        write(out / "best_so_far.json", compact(best))
        write(out / "cases.json", [compact(row) for row in records])

    def completed(result):
        with lock:
            records.append(result); persist(result)

    persist(seed)
    def search_worker(item): return run_child(helper_path, helper_sha256, joint_path,
        joint_sha256, scored_path, scored_sha256, template, item,
        out / "cases" / item["case_id"], search_deadline)

    center = dict(seed)
    coordinates = controller.derivative_proposals(center, restrictions, 0)
    if len(coordinates) > MAXIMUM_COORDINATES:
        raise RuntimeError("Coordinate proposal cap exceeded")
    coordinate_results, hard_error = bounded_batch(coordinates, search_worker, workers, completed)
    joint_results = []
    if not hard_error and time.monotonic() < search_deadline:
        try:
            jacobian = controller.jacobian(coordinate_results, center)
            write(out / "coordinate_jacobian.json", dict(names=list(controller.ALL_NAMES),
                                                           weighted_jacobian=jacobian.tolist()))
            joints = controller.joint_proposals(center, jacobian, restrictions, 0)[:MAXIMUM_JOINT]
            joint_results, hard_error = bounded_batch(joints, search_worker, workers, completed)
        except Exception as exc:
            write(out / "joint_stage_failure.json", dict(error_type=type(exc).__name__, error=str(exc)))
    if hard_error:
        aggregate_tables(out, records)
        write(out / "summary.json", dict(status="stopped_branch_hard_error",
            elapsed_seconds=time.monotonic()-started, completed_cases=len(records), best_loss=best["loss"]))
        raise SystemExit(2)
    if time.monotonic() >= deadline:
        aggregate_tables(out, records)
        write(out / "summary.json", dict(status="stopped_before_exact_repetitions",
            elapsed_seconds=time.monotonic()-started, completed_cases=len(records), best_loss=best["loss"]))
        raise SystemExit(2)
    selected = dict(case_id="selected_exact_repetitions",
        parameters={name: controller.parameters(best["score"])[name] for name in controller.ALL_NAMES},
        initial_psi=best["proposal"]["initial_psi"], repetitions=2)
    exact = run_child(helper_path, helper_sha256, joint_path, joint_sha256,
        scored_path, scored_sha256, template, selected,
        out / "cases" / selected["case_id"], deadline)
    completed(exact)
    if (exact["status"] != "verified" or not exact.get("second_signature_equal")
            or controller.numeric_signature(exact["score"]) != controller.numeric_signature(best["score"])):
        aggregate_tables(out, records)
        write(out / "summary.json", dict(status="failed_exact_repetitions",
            elapsed_seconds=time.monotonic()-started, completed_cases=len(records), best_loss=best["loss"]))
        raise SystemExit(2)
    best = exact; persist(exact); helper.write_selected_tables(out, exact); aggregate_tables(out, records)
    summary = dict(status="completed_rebated_initial_search", elapsed_seconds=time.monotonic()-started,
        workers=workers, coordinate_cases=len(coordinate_results), joint_cases=len(joint_results),
        completed_cases=len(records), selected_loss=exact["loss"], selected_exact_repetitions_verified=True,
        checkpoint=exact["accounting"][-1])
    write(out / "summary.json", summary); print(json.dumps(summary), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest="mode", required=True)
    candidate = sub.add_parser("candidate")
    candidate.add_argument("--helper", type=Path, required=True); candidate.add_argument("--helper-sha256", required=True)
    candidate.add_argument("--joint", type=Path, required=True); candidate.add_argument("--joint-sha256", required=True)
    candidate.add_argument("--scored", type=Path, required=True); candidate.add_argument("--scored-sha256", required=True)
    candidate.add_argument("--template", type=Path, required=True)
    candidate.add_argument("--item", type=Path, required=True); candidate.add_argument("--output", type=Path, required=True)
    run = sub.add_parser("run")
    run.add_argument("--helper", type=Path, required=True); run.add_argument("--helper-sha256", required=True)
    run.add_argument("--joint", type=Path, required=True); run.add_argument("--joint-sha256", required=True)
    run.add_argument("--scored", type=Path, required=True); run.add_argument("--scored-sha256", required=True)
    run.add_argument("--template", type=Path, required=True); run.add_argument("--smoke", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True); run.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    run.add_argument("--wall-seconds", type=int, default=MAXIMUM_WALL_SECONDS)
    args = parser.parse_args()
    if args.mode == "candidate": candidate_mode(args.helper, args.helper_sha256,
        args.joint, args.joint_sha256, args.scored, args.scored_sha256,
        args.template, args.item, args.output)
    else: run_search(args.helper, args.helper_sha256, args.joint, args.joint_sha256,
        args.scored, args.scored_sha256, args.template, args.smoke, args.output,
        args.workers, args.wall_seconds)


if __name__ == "__main__": main()
