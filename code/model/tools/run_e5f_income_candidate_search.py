"""Bounded, deterministic native-income coordinate poll.

This is a diagnostic local controller.  It evaluates a finite set of fresh
plans through the existing candidate adapter and makes no SMM or optimum claim.
The evaluator hook is intentionally small so the loop can be tested without a
native solve.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import signal
import subprocess
import sys
import time
import threading
import csv
import shutil
import gzip
import pickle
import psutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable

PARAMETERS = ("beta_annual", "kappa_fert", "kappa_fert_continuation", "chi", "H0",
              "theta0", "theta1", "first_birth_fixed_cost", "h_P")
OBJECTIVE = "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1"
DEFAULT_PSI = 0.2429719621740803
ORIGINAL_PSI = 0.1489153145785918
STOP = threading.Event()

class ContractError(ValueError):
    pass


def read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                  allow_nan=False).encode()).hexdigest()


def _numbers(value: Any):
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if math.isfinite(float(value)):
            yield float(value)
    elif isinstance(value, dict):
        for item in value.values():
            yield from _numbers(item)
    elif isinstance(value, list):
        for item in value:
            yield from _numbers(item)


def extract_score(payload: dict[str, Any]) -> float:
    """Read the scalar objective without depending on a report-only layout."""
    for key in ("objective", "loss", "score", "weighted_loss", "objective_value"):
        if key in payload and isinstance(payload[key], (int, float)):
            value = float(payload[key])
            if math.isfinite(value):
                return value
    for key in ("summary", "result", "fit", "scored"):
        if isinstance(payload.get(key), dict):
            try:
                return extract_score(payload[key])
            except ValueError:
                pass
    raise ValueError("score receipt has no finite objective scalar")


def validate_score_contract(payload, plan):
    if payload.get("schema") != "e5f_initial_minimum_distance_result_v1" or payload.get("contract_sha256") != OBJECTIVE:
        raise ContractError("native objective/schema mismatch")
    if not math.isfinite(float(payload.get("loss", math.nan))):
        raise ContractError("nonfinite native loss")
    if {p["parameter"] for p in payload.get("parameters", []) if p.get("structural_coordinate")} != set(PARAMETERS):
        raise ContractError("missing structural coordinates")
    if len(payload["parameters"]) != 17:
        raise ContractError("incomplete parameter table")
    n = payload.get("normalization", {})
    if n.get("target") != 2.1 or not math.isfinite(float(n.get("absolute_gap", math.nan))) or float(n["absolute_gap"]) > 5e-4:
        raise ContractError("normalization gate mismatch")
    summary = payload.get("_summary", {})
    if summary.get("status") != "verified_scored_candidate" or summary.get("objective_canonical_sha256") != OBJECTIVE or summary.get("loss") != payload["loss"]:
        raise ContractError("unverified native score")
    if payload.get("candidate_payload_fingerprint") != plan["candidate_payload_fingerprint"]:
        raise ContractError("income identity mismatch")
    rows = payload.get("target_fit", [])
    if len(rows) != 13 or sum(bool(x.get("scored")) for x in rows) != 12:
        raise ContractError("incomplete target fit")
    signature = [(x["restriction_id"], x["target"], x.get("actual_weight"), x.get("scored")) for x in rows]
    if plan.get("target_signature") is not None and fingerprint(signature) != plan["target_signature"]:
        raise ContractError("target/weight fingerprint mismatch")


def parameter_bounds(plan):
    raw = plan["parameter_bounds"]
    if set(raw) != set(PARAMETERS): raise ContractError("all nine bounds required")
    return {n: tuple(map(float, raw[n])) for n in PARAMETERS}


def validate_parameters(receipt, expected):
    actual = {r["parameter"]: r["estimate"] for r in receipt["parameters"]}
    if any(not math.isclose(float(actual[n]), v, rel_tol=0., abs_tol=1e-12)
           for n, v in expected.items()):
        raise ContractError("scored parameters differ from requested parameters")


def load_verified_score(score_path, plan):
    score_path = Path(score_path)
    evaluation = score_path.parent.parent
    score = read(score_path)
    score["_summary"] = read(evaluation / "summary.json")
    initial = read(evaluation.parent / "initial_contract.json")
    score["candidate_payload_fingerprint"] = initial["income_candidate_payload_fingerprint"]
    if initial["income_candidate_json_sha256"] != plan["candidate_json_sha256"] or initial["income_constructor_sha256"] != plan["constructor_sha256"]:
        raise ContractError("native income constructor/candidate mismatch")
    score["_evaluation"] = str(evaluation)
    score["_initial_psi"] = initial["initial_psi"]
    score["_native_price"] = read(evaluation / "raw/repetition_01/summary.json")["price"]
    validate_score_contract(score, plan)
    return score


def proposals(seed: dict[str, float], plan: dict[str, Any], max_proposals: int = 16) -> list[dict[str, float]]:
    """Return one feasible direction per coordinate, then its opposites.

    The base point is the incumbent and is therefore not returned as a proposal.
    Clipping is only performed against bounds pinned in the plan; duplicate
    points are removed by exact serialized coordinates.
    """
    steps = {"beta_annual": .005, "kappa_fert": .20, "kappa_fert_continuation": .20,
             "chi": .10, "H0": .10, "theta0": .04, "theta1": .25,
             "first_birth_fixed_cost": .075, "h_P": .15}
    bounds = parameter_bounds(plan)
    seen = {fingerprint(seed)}
    first: list[dict[str, float]] = []
    second: list[dict[str, float]] = []
    # beta starts toward the interior because the pilot is at its upper bound.
    signs = {name: (-1 if name in ("beta_annual", "h_P") else 1) for name in PARAMETERS}
    for phase, target in ((signs, first), ({n: -s for n, s in signs.items()}, second)):
        for name in PARAMETERS:
            value = float(seed[name]); step = steps[name]
            delta = step if name in ("beta_annual", "theta0", "h_P", "first_birth_fixed_cost") else value * step
            candidate = dict(seed); candidate[name] = value + phase[name] * delta
            if name in bounds:
                lo, hi = bounds[name]
                candidate[name] = min(hi, max(lo, candidate[name]))
            if not math.isfinite(candidate[name]) or fingerprint(candidate) in seen:
                continue
            seen.add(fingerprint(candidate)); target.append(candidate)
    return (first + second)[:max_proposals]


def _kill_process_group(proc: subprocess.Popen) -> None:
    # The frozen wrapper creates a fresh session for its native child.  Kill
    # descendants first so a timed-out case cannot leave a solver behind.
    try:
        import psutil
        parent = psutil.Process(proc.pid)
        descendants = parent.children(recursive=True)
        for child in reversed(descendants):
            try:
                child.terminate()
            except psutil.Error:
                pass
        _, alive = psutil.wait_procs(descendants, timeout=3)
        for child in alive:
            try:
                child.kill()
            except psutil.Error:
                pass
    except (psutil.NoSuchProcess, ProcessLookupError):
        pass
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(timeout=5)
    except ProcessLookupError:
        pass


def adapter_evaluator(plan_path: Path, output: Path, parameters: dict[str, float],
                      *, psi: float, repetitions: int, timeout: float,
                      case_id: str) -> dict[str, Any]:
    param_path = output.parent / (output.name + ".parameters.json")
    write(param_path, parameters)
    command = [sys.executable, str(Path(__file__).resolve()), "--mode", "adapter-child",
               "--plan", str(plan_path), "--output", str(output),
               "--parameters-json", str(param_path), "--initial-psi", str(psi),
               "--repetitions", str(repetitions), "--case-id", case_id]
    log = output.parent / (output.name + ".log")
    with log.open("w") as stream:
        proc = subprocess.Popen(command, start_new_session=True, stdout=stream, stderr=subprocess.STDOUT)
        try:
            deadline = time.monotonic() + timeout
            last_heartbeat = 0.
            while proc.poll() is None:
                if STOP.is_set(): raise RuntimeError("search interrupted")
                if time.monotonic() >= deadline: raise TimeoutError(f"candidate exceeded {timeout:.0f}s")
                if time.monotonic() - last_heartbeat >= 30:
                    write(output.parent / (output.name + ".heartbeat.json"), {"status":"running", "pid":proc.pid,"updated":time.time()})
                    last_heartbeat = time.monotonic()
                time.sleep(.5)
        except BaseException:
            _kill_process_group(proc)
            raise
    if proc.returncode:
        detail = log.read_text()[-2500:]
        if any(word in detail.lower() for word in ("fingerprint", "hash mismatch", "contracterror", "source inventory", "preflight")):
            raise ContractError(detail)
        raise RuntimeError(f"adapter failed ({proc.returncode}): {detail}")
    return load_verified_score(output / "evaluation/scored_repetition_01/score.json", read(plan_path))


def _array_payload(payload: dict[str, Any]) -> Any:
    for key in ("native_arrays", "arrays", "repeated_arrays", "price_arrays"):
        if key in payload:
            return payload[key]
    return None


def run_search(plan: dict[str, Any], output: Path, evaluator: Callable[..., dict[str, Any]],
               *, incumbent: dict[str, Any], search_seconds: float = 2100,
               max_proposals: int = 16, case_timeout: float = 900,
               verification_seconds: float = 1200, workers: int = 4) -> dict[str, Any]:
    validate_score_contract(incumbent, plan)
    incumbent_score = extract_score(incumbent)
    seed = {name: float(plan["pilot_parameters"][name]) for name in PARAMETERS}
    validate_parameters(incumbent, seed)
    output.mkdir(parents=True, exist_ok=False)
    write(output / "incumbent_score.json", incumbent)
    write(output / "latest_completed.json", {"status":"awaiting_first_case"})
    best = {"status": "incumbent", "objective": incumbent_score, "parameters": seed,
            "psi": ORIGINAL_PSI, "source": "immutable_pilot", "receipt": incumbent}
    write(output / "best_so_far.json", best)
    completed: list[dict[str, Any]] = []
    started = time.monotonic()
    all_proposals = proposals(seed, plan, max_proposals)
    for batch_start in range(0, len(all_proposals), max(1, workers)):
        batch = all_proposals[batch_start:batch_start + max(1, workers)]
        remaining = search_seconds - (time.monotonic() - started)
        if remaining < 500:
            write(output / "stop.json", {"status": "budget_stop", "remaining_seconds": remaining,
                                          "completed": len(completed)})
            break
        rows = {batch_start + j + 1: {"case": batch_start + j + 1, "parameters": params,
                                      "status": "started", "started": time.time()}
                for j, params in enumerate(batch)}
        for row in rows.values():
            write(output / "progress.json", row)
        def evaluate(index_params):
            index, params = index_params
            case = output / "cases" / f"case_{index:02d}"
            try:
                receipt = evaluator(params, case, psi=DEFAULT_PSI, repetitions=1,
                                     timeout=min(case_timeout, remaining, 900),
                                     case_id=f"income_search_{index:02d}")
                validate_score_contract(receipt, plan)
                validate_parameters(receipt, params)
                score = extract_score(receipt)
                return index, dict(rows[index], status="completed", objective=score, receipt=receipt)
            except TimeoutError as exc:
                return index, dict(rows[index], status="timed_out", error=str(exc))
            except ContractError:
                STOP.set()
                raise
            except (ValueError, RuntimeError, OSError) as exc:
                return index, dict(rows[index], status="failed", error=str(exc))
        with ThreadPoolExecutor(max_workers=min(workers, len(batch))) as pool:
            futures = [pool.submit(evaluate, (batch_start + j + 1, params))
                       for j, params in enumerate(batch)]
            for future in as_completed(futures):
                index, row = future.result()
                completed.append(row)
                if row.get("status") == "completed" and (row["objective"],index) < (best["objective"],best.get("case",0)):
                    best = {"status":"improved","case":index,"objective":row["objective"],"parameters":row["parameters"],"psi":DEFAULT_PSI,"receipt":row["receipt"]}
                    write(output / "best_so_far.json", best)
                write(output / "latest_completed.json", row)
                write(output / "cases.json", completed)
    selected = best
    verification = {"status": "not_run"}
    if selected.get("status") in ("improved", "incumbent"):
        remaining = verification_seconds
        try:
            verify_psi = DEFAULT_PSI if selected.get("status") == "improved" else ORIGINAL_PSI
            check = evaluator(selected["parameters"], output / "selected_verification",
                              psi=verify_psi, repetitions=2, timeout=remaining,
                              case_id="income_selected_verification")
            validate_score_contract(check, plan)
            validate_parameters(check, selected["parameters"])
            score = extract_score(check)
            exact = score == float(selected["objective"])
            summary_ok = check.get("_summary", {}).get("repetitions") == 2 and \
                check.get("_summary", {}).get("exact_loss_equality") is True
            def numeric_fit(receipt):
                return {"loss": receipt.get("loss"), "price": receipt.get("_native_price"),
                        "parameters": [(x.get("parameter"), x.get("estimate"))
                                       for x in receipt.get("parameters", [])
                                       ],
                        "target_fit": [(x.get("restriction_id", x.get("target")), x.get("target"),
                                        x.get("model"), x.get("gap"), x.get("actual_weight"), x.get("loss_contribution"))
                                       for x in receipt.get("target_fit", [])],
                        "normalization": (receipt.get("normalization", {}).get("completed_fertility"),
                                          receipt.get("normalization", {}).get("psi_child"))}
            numeric_equal = numeric_fit(selected.get("receipt", {})) == numeric_fit(check)
            verification = {"status": "verified" if exact and summary_ok and numeric_equal else "mismatch",
                             "objective": score, "exact_objective": exact,
                             "native_exact_loss_equality": summary_ok,
                             "numeric_fit_equal": numeric_equal, "receipt": check}
        except Exception as exc:  # preserve a reviewable failure receipt
            verification = {"status": "failed", "error": str(exc)}
    summary = {"status": "verified_selection" if verification["status"] == "verified" else "requires_review", "poll_complete": len(completed)==len(all_proposals), "interpretation": "bounded_coordinate_poll",
               "objective_canonical_sha256": OBJECTIVE, "incumbent_objective": incumbent_score,
               "selected": selected, "verification": verification, "cases": completed,
               "proposal_count": len(completed), "max_proposals": max_proposals,
               "search_seconds": search_seconds}
    chosen = verification.get("receipt") if verification.get("status") == "verified" else selected.get("receipt")
    if chosen:
        for name, rows in (("selected_target_fit", chosen["target_fit"]), ("selected_parameters", chosen["parameters"])):
            table = [dict(r) for r in rows]
            if name == "selected_parameters":
                for row in table:
                    if row["parameter"] == "beta_annual": row.update(active_upper=.99,active_near_bound=abs(row["estimate"]-.99)<=.0005)
            fields = list(dict.fromkeys(k for r in table for k in r))
            with (output / (name + ".csv")).open("w",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");writer.writeheader();writer.writerows(table)
        if chosen.get("_evaluation"):
            source = Path(chosen["_evaluation"])/"raw"/("repetition_02" if verification["status"]=="verified" else "repetition_01")/"standard_diagnostics"
            shutil.copytree(source,output/"selected_standard_diagnostics")
            if len(list((output/"selected_standard_diagnostics").glob("*.png"))) != 17: raise ContractError("missing selected graphs")
    write(output / "summary.json", summary)
    return summary


def preflight(plan, output):
    """Zero-solve wrapper and actual incumbent income-payload checks."""
    from run_e5f_income_candidate_calibration import run_scored_pilot, candidate_overrides, sha
    import numpy as np
    source = Path(plan["source_root"])
    sys.path[:0] = [str(source / "code/model/tools"), str(source / "code/model")]
    checkpoint = Path(plan["incumbent_checkpoint"])
    if sha(checkpoint) != plan["incumbent_checkpoint_sha256"]:
        raise ContractError("incumbent checkpoint hash mismatch")
    with gzip.open(checkpoint, "rb") as stream:
        packet = pickle.load(stream)
    P = packet["parameters"]
    candidate = read(Path(plan["candidate_json"]))
    expected = candidate_overrides(candidate)
    errors = {k: float(np.max(np.abs(np.asarray(getattr(P,k)) - np.asarray(expected[k]))))
              for k in ("z_grid", "z_weights", "Pi_z")}
    z, w, pi = map(np.asarray, (P.z_grid,P.z_weights,P.Pi_z))
    logz = np.log(z); centered = logz - w @ logz
    annual = candidate["annual_coefficients_recovered_from_nested_fitted_covariances"]
    vp = float(annual["persistent_variance"])
    ve = float(candidate["four_year_diagnostic_mapping"]["transitory_log_sd_period"])**2
    rho = float(annual["rho_annual"])**4
    moments = {"mean":float(w@z),"log_variance":float(w@centered**2),
               "log_covariance_lag1":float((w*centered)@pi@centered),
               "log_covariance_lag2":float((w*centered)@pi@pi@centered)}
    desired = dict(mean=1.,log_variance=vp+ve,log_covariance_lag1=vp*rho,log_covariance_lag2=vp*rho**2)
    if max(errors.values()) > 1e-12 or P.permanent_income_levels_enabled or len(z)!=15:
        raise ContractError("checkpoint income payload mismatch")
    if P.income_candidate_fingerprint != plan["candidate_payload_fingerprint"]:
        raise ContractError("checkpoint income identity mismatch")
    if any(abs(moments[k]-desired[k])>1e-12 for k in desired):
        raise ContractError("checkpoint income log moments mismatch")
    actual = {n: float(P.beta)**.25 if n=="beta_annual" else
              float(P.hbar_first_child_jump)+float(P.hbar_child_rooms) if n=="h_P" else
              float(np.asarray(getattr(P,n)).reshape(-1)[0]) for n in PARAMETERS}
    if any(abs(actual[n]-plan["pilot_parameters"][n])>1e-12 for n in PARAMETERS):
        raise ContractError("checkpoint structural parameters mismatch")
    del packet, P
    point = proposals(plan["pilot_parameters"],plan)[0]
    receipt = run_scored_pilot(plan,output,preflight_only=True,parameters=point,
                              initial_psi=DEFAULT_PSI,repetitions=2,case_id="income_search_preflight")
    write(output/"income_payload_audit.json",dict(status="passed",max_errors=errors,
          moments=moments,expected_moments=desired,parameters=actual,solves=0))
    return receipt


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=("search", "adapter-child", "preflight"), required=True)
    p.add_argument("--plan", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--incumbent-score", type=Path)
    p.add_argument("--parameters-json", type=Path)
    p.add_argument("--initial-psi", type=float)
    p.add_argument("--repetitions", type=int, default=1)
    p.add_argument("--case-id", default="income_candidate_pilot")
    args = p.parse_args()
    plan = read(args.plan)
    plan["plan_path"] = str(args.plan.resolve())
    from run_e5f_income_candidate_calibration import validate_plan, sha
    validate_plan(plan,require_source=True)
    if sha(Path(__file__)) != plan["controller_sha256"]: raise ContractError("controller hash mismatch")
    def handle_stop(signum, frame):
        STOP.set()
        raise KeyboardInterrupt("controller terminated")
    signal.signal(signal.SIGTERM,handle_stop)
    signal.signal(signal.SIGINT,handle_stop)
    if args.mode == "preflight":
        print(json.dumps(preflight(plan,args.output)))
        return
    if args.mode == "adapter-child":
        from run_e5f_income_candidate_calibration import run_scored_pilot
        plan["plan_path"] = str(args.plan.resolve())
        params = read(args.parameters_json)
        result = run_scored_pilot(plan, args.output, parameters=params, initial_psi=args.initial_psi,
                                  repetitions=args.repetitions, case_id=args.case_id)
        evaluation = Path(result["output"])
        score = evaluation / "scored_repetition_01" / "score.json"
        if not score.exists():
            score = evaluation / "score.json"
        if not score.exists():
            raise FileNotFoundError("adapter completed without score.json")
        receipt = read(score)
        summary_path = evaluation / "summary.json"
        if summary_path.exists():
            receipt["_summary"] = read(summary_path)
        write(args.output / "score.json", receipt)
        return
    if args.incumbent_score is None:
        raise ValueError("search requires --incumbent-score")
    if sha(args.incumbent_score) != plan["incumbent_score_sha256"]: raise ContractError("incumbent hash mismatch")
    incumbent = load_verified_score(args.incumbent_score,plan)
    bounds = parameter_bounds(plan)
    objective = read(Path(plan["source_manifest_path"]))
    for row in objective["parameter_restrictions"]:
        name=row["parameter"]; expected=(row["lower"], .99 if name=="beta_annual" else row["upper"])
        if bounds[name] != expected: raise ContractError("bounds changed: "+name)
    result = run_search(plan, args.output, lambda params, case, **kw: adapter_evaluator(
        args.plan, case, params, **kw), incumbent=incumbent,
        max_proposals=min(16, int(plan.get("search_max_proposals", 16) or 16)))
    if result["status"] != "verified_selection":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
