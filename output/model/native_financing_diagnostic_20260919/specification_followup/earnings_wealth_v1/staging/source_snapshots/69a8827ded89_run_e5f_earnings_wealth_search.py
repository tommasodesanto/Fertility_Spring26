"""Bounded adaptive search for the earnings and wealth diagnostic.

The final plan owns the arm, objective/source fingerprints, nine-coordinate
bounds, and all budgets.  This controller only proposes points and invokes the
existing earnings-wealth adapter; it makes no calibration or adoption claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import signal
import subprocess
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable

PARAMETERS = (
    "beta_annual", "kappa_fert", "kappa_fert_continuation", "chi", "H0",
    "theta0", "theta1", "first_birth_fixed_cost", "h_P",
)
OBJECTIVE = "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1"
FATAL_WORDS = ("contract", "fingerprint", "hash mismatch", "source inventory",
               "code error", "unexpected schema", "objective mismatch", "traceback")
GATE_WORDS = ("market clearing gate", "numerical gate", "infeasible choice",
              "nonfinite loss", "not converged")
STOP = threading.Event()


class ContractError(ValueError):
    pass


def read(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def write(path: Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def score(receipt: dict[str, Any]) -> float:
    for key in ("loss", "objective", "score", "weighted_loss"):
        value = receipt.get(key)
        if isinstance(value, (int, float)) and math.isfinite(float(value)):
            return float(value)
    raise ContractError("receipt has no finite objective")


def validate_plan(plan: dict[str, Any]) -> None:
    if plan.get("objective_canonical_sha256") != OBJECTIVE:
        raise ContractError("objective canonical fingerprint mismatch")
    bounds = plan.get("parameter_bounds")
    start = plan.get("starting_structural_parameters", plan.get("structural_parameters"))
    if not isinstance(bounds, dict) or set(bounds) != set(PARAMETERS):
        raise ContractError("plan must pin bounds for all nine coordinates")
    if not isinstance(start, dict) or set(start) != set(PARAMETERS):
        raise ContractError("plan must pin starting structural parameters")
    for name in PARAMETERS:
        lo, hi = map(float, bounds[name])
        value = float(start[name])
        if not (math.isfinite(lo) and math.isfinite(hi) and lo <= value <= hi):
            raise ContractError(f"invalid bound/start for {name}")
    if float(bounds["beta_annual"][1]) > .99 + 1e-12:
        raise ContractError("beta upper bound exceeds .99")
    if float(bounds["h_P"][1]) > 2.3 + 1e-12:
        raise ContractError("h_P upper bound exceeds 2.3")
    cfg = plan.get("search_config", {})
    if int(cfg.get("max_proposals", 64)) > 64:
        raise ContractError("proposal cap exceeds 64")
    if int(cfg.get("workers", 4)) < 1:
        raise ContractError("invalid worker count")
    if not plan.get("files") or not plan.get("source_root") or not plan.get("source_fingerprints"):
        raise ContractError("runtime source/file pins are required")


def _bounds(plan: dict[str, Any]) -> dict[str, tuple[float, float]]:
    return {n: tuple(map(float, plan["parameter_bounds"][n])) for n in PARAMETERS}


def _clip(point: dict[str, float], bounds: dict[str, tuple[float, float]]) -> dict[str, float]:
    return {n: min(bounds[n][1], max(bounds[n][0], float(point[n]))) for n in PARAMETERS}


def proposal_batch(center: dict[str, float], plan: dict[str, Any], rng: random.Random,
                   count: int, seen: set[str], scale: float) -> list[dict[str, float]]:
    """Create deterministic multivariate proposals; every coordinate moves."""
    bounds = _bounds(plan)
    positive = {"kappa_fert", "kappa_fert_continuation", "chi", "H0", "theta1", "h_P"}
    widths = plan.get("search_config", {}).get("coordinate_widths", {})
    out = []
    for _ in range(count * 8):
        point = {}
        moved = True
        for name in PARAMETERS:
            width = float(widths.get(name, .01 if name == "beta_annual" else .05))
            u = rng.uniform(-1.0, 1.0)
            if abs(u) < 1e-9:
                moved = False
            point[name] = center[name] * math.exp(scale * width * u) if name in positive else center[name] + scale * width * u
        point = _clip(point, bounds)
        key = fingerprint(point)
        if moved and key not in seen and any(abs(point[n] - center[n]) > 1e-14 for n in PARAMETERS):
            seen.add(key)
            out.append(point)
            if len(out) >= count:
                break
    return out


def _kill_process_group(proc: subprocess.Popen[Any]) -> None:
    try:
        import psutil
        parent = psutil.Process(proc.pid)
        children = parent.children(recursive=True)
        for child in reversed(children):
            child.terminate()
        _, alive = psutil.wait_procs(children, timeout=3)
        for child in alive:
            child.kill()
    except Exception:
        pass
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(timeout=5)
    except (ProcessLookupError, OSError):
        pass


def classify_failure(error: BaseException) -> str:
    text = str(error).lower()
    if any(word in text for word in FATAL_WORDS):
        return "fatal_contract_error"
    if any(word in text for word in GATE_WORDS):
        return "numerical_gate_rejection"
    return "fatal_unexpected_error"


def _numeric_fit(receipt: dict[str, Any]) -> Any:
    return {
        "loss": receipt.get("loss"),
        "target_fit": [(r.get("restriction_id", r.get("target")), r.get("target"), r.get("model"), r.get("gap"), r.get("loss_contribution")) for r in receipt.get("target_fit", [])],
        "parameters": [(r.get("parameter"), r.get("estimate")) for r in receipt.get("parameters", [])],
        "price": receipt.get("_native_price"),
        "normalization": receipt.get("normalization"),
    }


def _validate_receipt(receipt: dict[str, Any], point: dict[str, float], *, repetitions: int) -> None:
    if receipt.get("schema") != "e5f_initial_minimum_distance_result_v1":
        raise ContractError("native objective/schema mismatch")
    if len(receipt.get("target_fit", [])) != 13 or len(receipt.get("parameters", [])) != 17:
        raise ContractError("receipt must contain 13 targets and 17 parameters")
    actual = {r.get("parameter"): float(r.get("estimate")) for r in receipt["parameters"] if r.get("structural_coordinate")}
    if set(actual) != set(PARAMETERS) or any(not math.isclose(actual[n], point[n], rel_tol=0, abs_tol=1e-10) for n in PARAMETERS):
        raise ContractError("receipt structural parameters differ from requested point")
    summary = receipt.get("_summary", receipt.get("summary", {}))
    if summary.get("status") != "verified_scored_candidate":
        raise ContractError("native summary is not verified")
    if summary.get("objective_canonical_sha256") != OBJECTIVE:
        raise ContractError("native objective fingerprint mismatch")
    if summary.get("repetitions") != repetitions:
        raise ContractError("native repetition count mismatch")
    if repetitions == 2 and (summary.get("repetitions") != 2 or summary.get("exact_loss_equality") is not True):
        raise ContractError("two repetition exactness receipt missing")
    score(receipt)


def subprocess_evaluator(plan_path: Path, output: Path, point: dict[str, float], *, psi: float,
                         repetitions: int, timeout: float, case_id: str) -> dict[str, Any]:
    """Invoke the existing adapter in a fresh process and kill its full group on timeout."""
    plan = read(plan_path)
    derived = dict(plan)
    derived["parent_plan_sha256"] = hashlib.sha256(plan_path.read_bytes()).hexdigest()
    derived["structural_parameters"] = dict(point)
    derived["starting_structural_parameters"] = dict(point)
    derived["initial_psi"] = float(psi)
    arm = plan.get("search_config", {}).get("arm") or plan.get("search_arm")
    if not arm:
        raise ContractError("search arm is missing")
    derived["cases"] = [{"id": case_id, "arm": arm, "repetitions": repetitions,
                          "native_seconds": min(int(plan.get("search_config", {}).get("native_seconds", 1800)) * repetitions, int(timeout) - 100),
                          "seconds": int(timeout), "wrapper_seconds": int(timeout) - 50}]
    derived_path = output.parent / (output.name + ".plan.json")
    write(derived_path, derived)
    # The adapter creates the case directory itself with exist_ok=False.
    log = output.parent / (output.name + ".log")
    command = [str(plan.get("python", sys.executable)), str(plan["adapter_path"]), "--plan", str(derived_path),
               "--output", str(output), "--arm", arm, "--repetitions", str(repetitions)]
    with log.open("w") as stream:
        proc = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        deadline = time.monotonic() + timeout
        next_heartbeat = 0.
        try:
            while proc.poll() is None:
                if STOP.is_set():
                    raise ContractError("batch stopped after an unexpected failure")
                if time.monotonic() >= deadline:
                    raise TimeoutError(f"case exceeded {timeout:.0f}s")
                if time.monotonic() >= next_heartbeat:
                    write(output.parent / (output.name + ".heartbeat.json"),
                          {"status": "running", "pid": proc.pid, "updated": time.time(), "case": case_id})
                    next_heartbeat = time.monotonic() + 30.
                time.sleep(.5)
        except BaseException:
            _kill_process_group(proc)
            raise
    if proc.returncode:
        detail = log.read_text()[-3000:]
        raise RuntimeError(detail or f"adapter exited {proc.returncode}")
    rep = output / "evaluation" / "scored_repetition_01" / "score.json"
    receipt = read(rep)
    receipt["_evaluation"] = str(output / "evaluation")
    receipt["_summary"] = read(output / "evaluation" / "summary.json")
    receipt["_native_price"] = read(output / "evaluation/raw/repetition_01/summary.json")["price"]
    receipt["_runtime_contract"] = read(output / "runtime_contract.json")
    if receipt["_runtime_contract"].get("plan_sha256") != hashlib.sha256(derived_path.read_bytes()).hexdigest():
        raise ContractError("runtime economic plan fingerprint mismatch")
    receipt["_parent_plan_sha256"] = derived["parent_plan_sha256"]
    return receipt


def run_search(plan: dict[str, Any], output: Path, evaluator: Callable[..., dict[str, Any]], *,
               anchor: dict[str, Any], search_seconds: float | None = None) -> dict[str, Any]:
    validate_plan(plan)
    STOP.clear()
    _validate_receipt(anchor, plan.get("starting_structural_parameters", plan["structural_parameters"]), repetitions=2)
    cfg = plan.get("search_config", {})
    budget = float(search_seconds if search_seconds is not None else cfg.get("search_seconds", 14400))
    max_proposals = min(64, int(cfg.get("max_proposals", 64)))
    workers = max(1, int(cfg.get("workers", 4)))
    timeout = float(cfg.get("case_timeout", 1800))
    verification_budget = float(cfg.get("verification_seconds", 3600))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    write(output / "anchor.json", anchor)
    best = {"status": "anchor", "objective": score(anchor), "parameters": plan.get("starting_structural_parameters", plan["structural_parameters"]), "receipt": anchor}
    write(output / "best_so_far.json", best)
    rows: list[dict[str, Any]] = []
    seen = {fingerprint(best["parameters"])}
    rng = random.Random(int(cfg.get("seed", 20260921)))
    center = dict(best["parameters"])
    scale = 1.0
    index = 0
    while index < max_proposals and time.monotonic() - start < budget - verification_budget - timeout:
        batch = proposal_batch(center, plan, rng, min(workers, max_proposals - index), seen, scale)
        if not batch:
            break
        def evaluate(item: tuple[int, dict[str, float]]) -> dict[str, Any]:
            i, point = item
            case = output / "cases" / f"case_{i:03d}"
            try:
                receipt = evaluator(point, case, psi=float(plan.get("initial_psi", cfg.get("initial_psi", 0.0))), repetitions=1, timeout=timeout, case_id=f"earnings_wealth_search_{i:03d}")
                _validate_receipt(receipt, point, repetitions=1)
                return {"case": i, "status": "completed", "objective": score(receipt), "parameters": point, "receipt": receipt}
            except Exception as exc:
                kind = classify_failure(exc)
                if kind != "numerical_gate_rejection":
                    STOP.set()
                    write(output / "failure.json", {"case": i, "status": kind, "error": str(exc)})
                    raise ContractError(str(exc)) from exc
                return {"case": i, "status": kind, "parameters": point, "error": str(exc)}
        with ThreadPoolExecutor(max_workers=min(workers, len(batch))) as pool:
            futures = [pool.submit(evaluate, (index + j + 1, p)) for j, p in enumerate(batch)]
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                write(output / "latest.json", row)
                write(output / "heartbeat.json", {"status": "running", "completed": len(rows), "updated": time.time()})
                with (output / "cases.jsonl").open("a") as stream:
                    stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
                if row["status"] == "completed" and row["objective"] < best["objective"]:
                    best = {"status": "improved", **row}
                    center = dict(row["parameters"])
                    scale = .5
                    write(output / "best_so_far.json", best)
        if batch and all(row["status"] == "numerical_gate_rejection" for row in rows[-len(batch):]):
            write(output / "stop.json", {"status": "all_cases_gate_rejected", "completed": len(rows)})
            raise ContractError("all dispatched cases rejected by the native numerical gate")
        index += len(batch)
    verification: dict[str, Any]
    try:
        check = evaluator(best["parameters"], output / "selected_verification", psi=float(plan.get("initial_psi", cfg.get("initial_psi", 0.0))), repetitions=2, timeout=verification_budget, case_id="earnings_wealth_selected_exact_repetition")
        _validate_receipt(check, best["parameters"], repetitions=2)
        summary = check["_summary"]
        if len(summary.get("original_graphs", [])) != 17:
            raise ContractError("selected verification must retain 17 native diagnostic plots")
        if check.get("_evaluation"):
            evaluation = Path(check["_evaluation"])
            for rel in plan.get("native_artifact_paths", []):
                if not (evaluation / "raw" / "repetition_02" / rel).is_file():
                    raise ContractError(f"missing pinned native artifact: {rel}")
        verification = {"status": "verified" if score(check) == best["objective"] and _numeric_fit(check) == _numeric_fit(best["receipt"]) else "mismatch", "receipt": check, "exact_objective": score(check) == best["objective"], "numeric_fit_equal": _numeric_fit(check) == _numeric_fit(best["receipt"])}
    except Exception as exc:
        verification = {"status": "failed", "error": str(exc)}
    result = {"schema": "e5f_earnings_wealth_search_v1", "status": "verified_selection" if verification["status"] == "verified" else "requires_review", "selected": best, "verification": verification, "cases": rows, "proposal_count": len(rows), "max_proposals": max_proposals, "objective_canonical_sha256": OBJECTIVE, "seed": cfg.get("seed", 20260921), "elapsed_seconds": time.monotonic()-start, "total_stage_budget_seconds": budget}
    write(output / "summary.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("smoke", "search"), required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verified-smoke", type=Path)
    args = parser.parse_args()
    plan = read(args.plan)
    validate_plan(plan)
    import run_e5f_earnings_wealth_candidate as adapter
    adapter.verify_plan(plan)
    cfg = plan.get("search_config", {})
    arm = cfg.get("arm") or plan.get("search_arm")
    if not arm:
        raise ContractError("search arm is missing")
    start = plan.get("starting_structural_parameters", plan["structural_parameters"])
    evaluator = lambda point, case, **kw: subprocess_evaluator(args.plan, case, point, **kw)
    if args.mode == "smoke":
        anchor = evaluator(start, args.output / "smoke_anchor", psi=float(plan["initial_psi"]), repetitions=2, timeout=float(cfg.get("smoke_anchor_seconds", 2 * cfg.get("case_timeout", 1800))), case_id="earnings_wealth_smoke_anchor")
        _validate_receipt(anchor, start, repetitions=2)
        probe = dict(start); bounds = _bounds(plan)
        name = "beta_annual" if start["beta_annual"] > bounds["beta_annual"][0] else "h_P"
        probe[name] = max(bounds[name][0], float(start[name]) - (.005 if name == "beta_annual" else .05))
        probe_receipt = evaluator(probe, args.output / "smoke_probe", psi=float(plan["initial_psi"]), repetitions=1, timeout=float(cfg.get("case_timeout", 1800)), case_id="earnings_wealth_smoke_probe")
        _validate_receipt(probe_receipt, probe, repetitions=1)
        write(args.output / "smoke_receipt.json", {"status": "verified_smoke", "anchor": anchor, "probe": probe_receipt, "plan_sha256": hashlib.sha256(args.plan.read_bytes()).hexdigest()})
        return
    if not args.verified_smoke:
        raise ContractError("production requires --verified-smoke")
    smoke = read(args.verified_smoke)
    if smoke.get("status") != "verified_smoke":
        raise ContractError("smoke receipt is not verified")
    if smoke.get("plan_sha256") != hashlib.sha256(args.plan.read_bytes()).hexdigest():
        raise ContractError("smoke plan hash differs")
    anchor = smoke["anchor"]
    result = run_search(plan, args.output, evaluator, anchor=anchor)
    if result["status"] != "verified_selection":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
