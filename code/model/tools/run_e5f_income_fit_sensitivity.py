"""Controlled finite-difference panel around the selected new-income fit.

The panel deliberately reuses the frozen income candidate adapter and scorer.
It does not alter the model, objective, target contract, or numerical gates.
The old controller is loaded from the pinned ``income_overnight_v1`` runtime so
its adapter-child dispatch continues to use the already validated code path.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import importlib.util
import json
import math
import os
import pickle
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable

PARAMETERS = (
    "beta_annual", "kappa_fert", "kappa_fert_continuation", "chi", "H0",
    "theta0", "theta1", "first_birth_fixed_cost", "h_P",
)
OBJECTIVE = "4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1"
ANCHOR_LOSS = 353.6588729140903
ANCHOR_PSI = 0.23949950404168222
MAX_STATIONARY_SOLVES = 8
FULL_PROBE_COUNT = 16
HALF_PROBE_COUNT = 8
MAX_EVALUATIONS = 28
EXPECTED_REJECTION_WORDS = (
    "nonfinite", "infeasible", "no feasible", "did not converge",
    "failed to converge", "line search", "bracket failure", "root solve",
)
FORBIDDEN_FAILURE_WORDS = (
    "source", "fingerprint", "hash mismatch", "target", "objective",
    "accounting", "budget", "validation", "contract", "preflight",
)
ARRAY_FILES = ("raw/repetition_01/policy_array_summary.json", "raw/repetition_01/initial_state.pkl.gz")


class ContractError(ValueError):
    pass


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def write_json(path: Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def load_old_controller(plan: dict[str, Any]):
    """Load the old controller by its pinned path, preserving its ``__file__``."""
    path = Path(plan["controller_path"])
    if not path.exists():
        raise ContractError(f"pinned old controller is unavailable: {path}")
    if plan.get("controller_sha256") and sha256(path) != plan["controller_sha256"]:
        raise ContractError("pinned old controller hash mismatch")
    spec = importlib.util.spec_from_file_location("e5f_income_candidate_search_pinned", path)
    if spec is None or spec.loader is None:
        raise ContractError("could not load pinned old controller")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_old_adapter(plan: dict[str, Any]):
    path = Path(plan["adapter_path"])
    if not path.exists():
        raise ContractError(f"pinned old adapter is unavailable: {path}")
    if plan.get("adapter_sha256") and sha256(path) != plan["adapter_sha256"]:
        raise ContractError("pinned old adapter hash mismatch")
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location("e5f_income_candidate_calibration_pinned", path)
    if spec is None or spec.loader is None:
        raise ContractError("could not load pinned old adapter")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def structural_parameters(receipt: dict[str, Any]) -> dict[str, float]:
    rows = [row for row in receipt.get("parameters", []) if row.get("structural_coordinate")]
    if {row.get("parameter") for row in rows} != set(PARAMETERS):
        raise ContractError("receipt does not contain all nine structural coordinates")
    return {str(row["parameter"]): float(row["estimate"]) for row in rows}


def numeric_signature(receipt: dict[str, Any]) -> dict[str, Any]:
    return {
        "loss": receipt.get("loss"),
        "parameters": [(row.get("parameter"), row.get("estimate")) for row in receipt.get("parameters", [])],
        "target_fit": [
            (row.get("restriction_id", row.get("target")), row.get("target"), row.get("model"),
             row.get("gap"), row.get("actual_weight"), row.get("loss_contribution"))
            for row in receipt.get("target_fit", [])
        ],
        "normalization": (
            receipt.get("normalization", {}).get("completed_fertility"),
            receipt.get("normalization", {}).get("psi_child"),
        ),
    }


def _hash_arrays(value: Any, digest: Any, path: str = "root", seen: set[int] | None = None) -> int:
    """Hash ndarray paths, dtypes, shapes and C-order bytes from a native checkpoint."""
    import numpy as np
    if seen is None: seen = set()
    if isinstance(value, np.ndarray):
        arr = np.ascontiguousarray(value)
        digest.update(path.encode()); digest.update(str(arr.dtype).encode()); digest.update(repr(arr.shape).encode()); digest.update(arr.tobytes())
        return 1
    if value is None or isinstance(value, (str, bytes, int, float, bool)):
        return 0
    marker = id(value)
    if marker in seen: return 0
    seen.add(marker)
    total = 0
    if isinstance(value, dict):
        for key in sorted(value, key=lambda x: repr(x)):
            total += _hash_arrays(value[key], digest, f"{path}.{key!r}", seen)
    elif isinstance(value, (list, tuple)):
        for i, item in enumerate(value): total += _hash_arrays(item, digest, f"{path}[{i}]", seen)
    elif hasattr(value, "__dict__"):
        total += _hash_arrays(vars(value), digest, path, seen)
    return total


def artifact_signature(evaluation: Path, source_root: str | Path | None = None, *, require_files: bool = True) -> dict[str, str]:
    """Hash the native checkpoint arrays and policy gate used for anchor checks."""
    out: dict[str, str] = {}
    for relative in ARRAY_FILES:
        path = Path(evaluation) / relative
        if not path.exists():
            if require_files:
                raise ContractError(f"missing repeated-array artifact: {path}")
            continue
        if relative.endswith("initial_state.pkl.gz"):
            if source_root:
                sys.path[:0] = [str(Path(source_root) / "code/model/tools"), str(Path(source_root) / "code/model")]
            digest = hashlib.sha256()
            with gzip.open(path, "rb") as stream:
                packet = pickle.load(stream)
            count = _hash_arrays(packet, digest)
            if count < 3:
                raise ContractError("native checkpoint lacks the required policy/stationary arrays")
            out[relative] = digest.hexdigest()
        else:
            out[relative] = sha256(path)
    if require_files and set(out) != set(ARRAY_FILES):
        raise ContractError("incomplete repeated-array artifact set")
    return out


def resolve_anchor(summary_path: Path, plan: dict[str, Any], old: Any, *, require_files: bool = True) -> dict[str, Any]:
    """Resolve the selected score through the saved summary, never a guessed path."""
    summary = read_json(summary_path)
    if summary.get("status") != "verified_selection":
        raise ContractError("selected search summary is not verified")
    verification = summary.get("verification", {})
    if verification.get("status") != "verified":
        raise ContractError("selected search verification is not exact")
    reference = verification.get("receipt", {})
    evaluation_text = reference.get("_evaluation")
    if not evaluation_text:
        raise ContractError("selected summary lacks verification evaluation path")
    evaluation = Path(evaluation_text)
    score_path = evaluation / "scored_repetition_01" / "score.json"
    if not score_path.exists():
        raise ContractError(f"summary-resolved selected score is missing: {score_path}")
    score = old.load_verified_score(score_path, plan)
    if not math.isclose(float(score["loss"]), ANCHOR_LOSS, rel_tol=0.0, abs_tol=1e-12):
        raise ContractError("selected loss differs from the frozen case-60 anchor")
    psi = float(summary.get("selected", {}).get("psi", score.get("_initial_psi", math.nan)))
    if not math.isclose(psi, ANCHOR_PSI, rel_tol=0.0, abs_tol=1e-15):
        raise ContractError("selected psi differs from the frozen case-60 anchor")
    selected = structural_parameters(score)
    summary_selected = summary.get("selected", {}).get("parameters", {})
    if summary_selected and any(abs(selected[k] - float(summary_selected[k])) > 1e-12 for k in PARAMETERS):
        raise ContractError("selected summary and score parameter tables disagree")
    arrays = artifact_signature(evaluation, plan.get("source_root"), require_files=require_files)
    return {
        "summary": summary,
        "score": score,
        "parameters": selected,
        "psi": psi,
        "evaluation": str(evaluation),
        "score_path": str(score_path),
        "numeric_signature": numeric_signature(score),
        "array_signature": arrays,
    }


def parameter_bounds(plan: dict[str, Any]) -> dict[str, tuple[float, float]]:
    raw = plan.get("parameter_bounds")
    if not isinstance(raw, dict) or set(raw) != set(PARAMETERS):
        raise ContractError("plan does not contain the complete nine-coordinate bound map")
    return {name: (float(raw[name][0]), float(raw[name][1])) for name in PARAMETERS}


def build_probe_specs(anchor: dict[str, float], bounds: dict[str, tuple[float, float]]) -> list[dict[str, Any]]:
    """Build 16 legal full steps and eight legal half steps in fixed order."""
    increments = {
        "beta_annual": 0.001,
        "h_P": 0.05,
        "first_birth_fixed_cost": 0.01,
    }
    for name in PARAMETERS:
        increments.setdefault(name, 0.05 * abs(anchor[name]))
    # Zero-valued coordinates are absent here; none of the six 5% coordinates
    # is zero at the selected point. Keep the failure explicit if that changes.
    if any(increments[name] <= 0.0 for name in PARAMETERS):
        raise ContractError("selected coordinate gives a zero finite-difference step")
    specs: list[dict[str, Any]] = []

    def add(name: str, sign: int, scale: float, group: str) -> None:
        step = increments[name] * scale
        proposed = anchor[name] + sign * step
        lo, hi = bounds[name]
        if proposed < lo - 1e-14 or proposed > hi + 1e-14:
            return
        point = dict(anchor)
        point[name] = proposed
        specs.append({"probe_id": f"{group}_{name}_{'plus' if sign > 0 else 'minus'}", "coordinate": name,
                      "sign": sign, "scale": scale, "step": step, "group": group, "parameters": point})

    for name in PARAMETERS:
        add(name, -1, 1.0, "full")
        add(name, 1, 1.0, "full")
    for name in ("beta_annual", "h_P", "chi", "H0", "kappa_fert"):
        add(name, -1, 0.5, "half")
        add(name, 1, 0.5, "half")
    if len([s for s in specs if s["group"] == "full"]) != FULL_PROBE_COUNT:
        raise ContractError("full-step panel does not contain 16 legal probes")
    if len([s for s in specs if s["group"] == "half"]) != HALF_PROBE_COUNT:
        raise ContractError("half-step panel does not contain 8 legal probes")
    return specs


def expected_rejection(error: BaseException) -> bool:
    text = str(error).lower()
    if any(word in text for word in FORBIDDEN_FAILURE_WORDS):
        return False
    return any(word in text for word in EXPECTED_REJECTION_WORDS)


def solve_count(receipt: dict[str, Any]) -> int | None:
    value = receipt.get("normalization", {}).get("stationary_solves")
    return int(value) if isinstance(value, (int, float)) else None


def compare_to_reference(receipt: dict[str, Any], evaluation: Path, reference: dict[str, Any], source_root: str | Path, *, arrays: bool = True) -> None:
    if numeric_signature(receipt) != reference["numeric_signature"]:
        raise ContractError("anchor objective/target/parameter table differs from saved case-60 reference")
    if arrays and artifact_signature(evaluation, source_root) != reference["array_signature"]:
        raise ContractError("anchor repeated-array artifacts differ from saved case-60 reference")


def _heartbeat(case_dir: Path, stop: threading.Event, payload: dict[str, Any]) -> None:
    while not stop.wait(30.0):
        write_json(case_dir / "heartbeat.json", {**payload, "status": "running", "updated": time.time()})


def evaluate_one(old: Any, plan_path: Path, plan: dict[str, Any], output: Path, *, parameters: dict[str, float],
                psi: float, case_id: str, anchor: dict[str, Any] | None = None,
                timeout: float = 900.0) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=False)
    bounds = parameter_bounds(plan)
    for name, value in parameters.items():
        if name not in bounds or not bounds[name][0] - 1e-14 <= float(value) <= bounds[name][1] + 1e-14:
            raise ContractError(f"requested coordinate {name} lies outside the active bound")
    write_json(output / "parameters.json", parameters)
    stop = threading.Event()
    thread = threading.Thread(target=_heartbeat, args=(output, stop, {"case_id": case_id}), daemon=True)
    thread.start()
    started = time.monotonic()
    try:
        receipt = old.adapter_evaluator(plan_path, output / "native", parameters, psi=psi, repetitions=1,
                                        timeout=timeout, case_id=case_id)
        old.validate_score_contract(receipt, plan)
        actual = structural_parameters(receipt)
        if any(abs(actual[name] - float(parameters[name])) > 1e-12 for name in PARAMETERS):
            raise ContractError("adapter returned a parameter table different from the requested point")
        evaluation = Path(receipt["_evaluation"])
        result = {"case_id": case_id, "status": "valid", "parameters": parameters,
                  "psi": psi, "objective": float(receipt["loss"]),
                  "stationary_solves": solve_count(receipt),
                  "elapsed_seconds": time.monotonic() - started,
                  "evaluation": str(evaluation), "numeric_signature": numeric_signature(receipt),
                  "array_signature": artifact_signature(evaluation, plan.get("source_root")), "receipt": receipt}
        write_json(output / "result.json", result)
        return result
    except Exception as exc:
        if expected_rejection(exc):
            result = {"case_id": case_id, "status": "expected_numerical_rejection", "parameters": parameters,
                      "psi": psi, "error": str(exc), "elapsed_seconds": time.monotonic() - started}
            write_json(output / "result.json", result)
            return result
        raise
    finally:
        stop.set()
        thread.join(timeout=1.0)


def write_tables(output: Path, rows: list[dict[str, Any]]) -> None:
    valid = [row for row in rows if row.get("status") == "valid"]
    fields = ["case_id", "status", "coordinate", "group", "sign", "scale", "step", "objective",
              "stationary_solves", "elapsed_seconds", "evaluation"]
    with (output / "cases.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)
    for row in valid:
        receipt = row.get("receipt", {})
        case_dir = Path(row["case_dir"])
        for filename, values in (("target_fit.csv", receipt.get("target_fit", [])),
                                 ("parameters.csv", receipt.get("parameters", []))):
            if not values: continue
            keys = list(dict.fromkeys(key for value in values for key in value))
            with (case_dir / filename).open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore")
                writer.writeheader(); writer.writerows(values)


def finite_difference_panel(rows: list[dict[str, Any]], anchor: dict[str, Any], specs: list[dict[str, Any]], output: Path) -> None:
    anchor_fit = anchor["score"].get("target_fit", [])
    base = {row.get("restriction_id", row.get("target")): float(row["model"]) for row in anchor_fit if row.get("scored")}
    weights = {row.get("restriction_id", row.get("target")): float(row.get("actual_weight") or 0.0) for row in anchor_fit if row.get("scored")}
    by_id = {row.get("probe_id"): row for row in rows}
    columns = []
    for spec in specs:
        row = by_id.get(spec["probe_id"])
        if not row or row.get("status") != "valid":
            continue
        receipt = row["receipt"]
        moments = {x.get("restriction_id", x.get("target")): float(x["model"]) for x in receipt.get("target_fit", []) if x.get("scored")}
        if set(moments) != set(base):
            continue
        scaled = {name: (moments[name] - base[name]) / (spec["sign"] * spec["step"]) for name in base}
        columns.append({"probe_id": spec["probe_id"], "coordinate": spec["coordinate"], "group": spec["group"],
                        "sign": spec["sign"], "step": spec["step"], "scaled_response": scaled,
                        "weighted_scaled_response": {name: value * math.sqrt(max(weights[name], 0.0)) for name, value in scaled.items()}})
    write_json(output / "finite_difference_jacobian.json", {"status": "descriptive_local_map", "columns": columns,
              "anchor_objective": ANCHOR_LOSS, "note": "Rank and conditioning are local diagnostics, not identification proofs."})


def run_panel(args: argparse.Namespace) -> dict[str, Any]:
    plan = read_json(args.plan)
    plan["plan_path"] = str(Path(args.plan).resolve())
    if args.expected_plan_sha256 and sha256(args.plan) != args.expected_plan_sha256:
        raise ContractError("sensitivity plan hash mismatch")
    if plan.get("objective_canonical_sha256") != OBJECTIVE:
        raise ContractError("plan objective fingerprint differs from the frozen objective")
    old = load_old_controller(plan)
    adapter = load_old_adapter(plan) if args.mode != "plan" else None
    if args.mode != "plan":
        adapter.validate_plan(plan, require_source=True)
    anchor = resolve_anchor(Path(args.selected_summary), plan, old, require_files=args.mode != "plan")
    bounds = parameter_bounds(plan)
    specs = build_probe_specs(anchor["parameters"], bounds)
    output = Path(args.output)
    if args.mode == "plan":
        output.mkdir(parents=True, exist_ok=False)
        write_json(output / "panel_plan.json", {"status": "prepared_not_submitted", "anchor": anchor,
                  "probe_count": len(specs), "full_probe_count": FULL_PROBE_COUNT, "half_probe_count": HALF_PROBE_COUNT,
                  "max_objective_evaluations": MAX_EVALUATIONS, "max_nested_stationary_solves": MAX_EVALUATIONS * MAX_STATIONARY_SOLVES,
                  "probes": specs, "target_fingerprint": plan.get("target_signature"), "objective": OBJECTIVE})
        return {"status": "prepared_not_submitted", "probe_count": len(specs)}
    output.mkdir(parents=True, exist_ok=False)
    finish = dt.datetime.fromisoformat(args.finish_utc.replace("Z", "+00:00"))
    seconds_to_finish = (finish - dt.datetime.now(dt.timezone.utc)).total_seconds()
    if seconds_to_finish <= 0:
        raise ContractError("absolute UTC finish cutoff has passed")
    global_deadline = time.monotonic() + min(float(args.global_timeout), seconds_to_finish)
    write_json(output / "panel_plan.json", {"status": "running", "anchor": {k: v for k, v in anchor.items() if k != "score"},
              "probes": specs, "max_objective_evaluations": MAX_EVALUATIONS,
              "max_nested_stationary_solves": MAX_EVALUATIONS * MAX_STATIONARY_SOLVES})
    rows: list[dict[str, Any]] = []
    beta_inward = next(spec for spec in specs if spec["probe_id"] == "full_beta_annual_minus")
    def run_case(spec: dict[str, Any], label: str) -> dict[str, Any]:
        if time.monotonic() >= global_deadline:
            raise ContractError("global panel deadline exceeded before dispatch")
        case_dir = output / "cases" / label
        result = evaluate_one(old, Path(args.plan), plan, case_dir, parameters=spec["parameters"], psi=anchor["psi"],
                              case_id=f"income_fit_sensitivity_{label}", timeout=args.case_timeout)
        result.update({k: spec.get(k) for k in ("probe_id", "coordinate", "group", "sign", "scale", "step")}, case_dir=str(case_dir))
        return result
    # Smoke is intentionally sequential: two fresh anchors and the first probe.
    if args.mode == "smoke":
        smoke_rows: list[dict[str, Any]] = []
        for label in ("smoke_anchor_01", "smoke_anchor_02"):
            result = evaluate_one(old, Path(args.plan), plan, output / "smoke" / label, parameters=anchor["parameters"],
                                  psi=anchor["psi"], case_id=f"income_fit_sensitivity_{label}", timeout=args.case_timeout)
            compare_to_reference(result["receipt"], Path(result["evaluation"]), anchor, plan.get("source_root"))
            result.update(case_id=label, case_dir=str(output / "smoke" / label), coordinate="anchor", group="anchor", step=0.0)
            smoke_rows.append(result)
            write_json(output / "latest_completed.json", result)
            write_json(output / "best_so_far.json", result)
        smoke_probe = run_case(beta_inward, "smoke_beta_annual_minus")
        if smoke_probe.get("status") != "valid":
            raise ContractError("smoke beta inward probe was rejected")
        smoke_rows.append(smoke_probe)
        write_json(output / "latest_completed.json", smoke_probe)
        write_json(output / "best_so_far.json", min(smoke_rows, key=lambda row: float(row["objective"])))
        write_json(output / "smoke_summary.json", {"status": "passed", "rows": smoke_rows, "anchor_loss": ANCHOR_LOSS,
                  "reference_score_path": anchor["score_path"], "reference_array_signature": anchor["array_signature"]})
        write_tables(output, smoke_rows)
        return {"status": "passed", "evaluations": len(smoke_rows), "stationary_solves": sum(r.get("stationary_solves") or 0 for r in smoke_rows)}

    if not args.smoke_summary:
        raise ContractError("production requires --smoke-summary")
    smoke = read_json(args.smoke_summary)
    if smoke.get("status") != "passed":
        raise ContractError("successful smoke summary is required")
    smoke_rows = list(smoke.get("rows", []))
    if len(smoke_rows) != 3 or any(row.get("status") != "valid" for row in smoke_rows):
        raise ContractError("smoke summary is incomplete")
    anchor_rows = [row for row in smoke_rows if row.get("group") == "anchor"]
    beta_rows = [row for row in smoke_rows if row.get("probe_id") == beta_inward["probe_id"]]
    if {row.get("case_id") for row in anchor_rows} != {"smoke_anchor_01", "smoke_anchor_02"} or len(beta_rows) != 1:
        raise ContractError("smoke summary does not contain two anchors and the named beta probe")
    if beta_rows[0].get("parameters") != beta_inward["parameters"]:
        raise ContractError("smoke beta probe differs from the declared panel")
    for row in smoke_rows:
        receipt = row.get("receipt")
        if not isinstance(receipt, dict):
            raise ContractError("smoke row lacks its complete score receipt")
        old.validate_score_contract(receipt, plan)
        actual = structural_parameters(receipt)
        if any(abs(actual[name] - float(row["parameters"][name])) > 1e-12 for name in PARAMETERS):
            raise ContractError("smoke row parameter table does not match its requested point")
        evaluation = Path(row["evaluation"])
        observed_arrays = artifact_signature(evaluation, plan.get("source_root"))
        if observed_arrays != row.get("array_signature"):
            raise ContractError("smoke row array receipt changed before production")
    if any(row.get("group") == "anchor" for row in smoke_rows):
        for row in smoke_rows:
            if row.get("group") == "anchor":
                compare_to_reference(row["receipt"], Path(row["evaluation"]), anchor, plan.get("source_root"))
    rows.extend(smoke_rows)
    anchor_row = {"case_id": "anchor_reference", "status": "valid", "parameters": anchor["parameters"], "objective": ANCHOR_LOSS,
                  "evaluation": anchor["evaluation"], "numeric_signature": anchor["numeric_signature"], "array_signature": anchor["array_signature"],
                  "group": "anchor", "coordinate": "anchor", "step": 0.0, "case_dir": str(output / "anchor_reference")}
    write_json(output / "best_so_far.json", anchor_row)
    remaining_specs = [spec for spec in specs if spec["probe_id"] != beta_inward["probe_id"]]
    def run_indexed(item: tuple[int, dict[str, Any]]) -> dict[str, Any]:
        index, spec = item
        return run_case(spec, f"probe_{index:02d}_{spec['probe_id']}")
    with ThreadPoolExecutor(max_workers=max(1, int(args.workers))) as pool:
        futures = [pool.submit(run_indexed, (index, spec)) for index, spec in enumerate(remaining_specs)]
        try:
            for future in as_completed(futures):
                result = future.result()
                rows.append(result)
                write_json(output / "latest_completed.json", result)
                valid_now = [row for row in rows if row.get("status") == "valid" and row.get("group") != "anchor"]
                if valid_now:
                    write_json(output / "best_so_far.json", min(valid_now, key=lambda row: (float(row["objective"]), row["case_id"])))
        except BaseException as exc:
            if hasattr(old, "STOP"):
                old.STOP.set()
            for future in futures:
                future.cancel()
            write_json(output / "failure.json", {"status": "fatal", "error": str(exc), "stop_reason": "unexpected_or_contract_failure"})
            raise
    valid = [row for row in rows if row.get("status") == "valid" and row.get("group") != "anchor"]
    valid.append(anchor_row)
    if not valid: raise ContractError("no valid probe or anchor remains for final repetition")
    selected = min(valid, key=lambda row: (float(row["objective"]), str(row["case_id"])))
    final_rows = []
    final_started = time.monotonic()
    for rep in ("final_rep_01", "final_rep_02"):
        if time.monotonic() >= global_deadline:
            raise ContractError("global panel deadline exceeded before final repetition")
        spec = {"parameters": selected["parameters"]}
        remaining = max(1.0, min(900.0, 1800.0 - (time.monotonic() - final_started)))
        result = evaluate_one(old, Path(args.plan), plan, output / "final" / rep, parameters=spec["parameters"], psi=anchor["psi"],
                              case_id=f"income_fit_sensitivity_{rep}", timeout=remaining)
        if result["status"] != "valid": raise ContractError("final selected repetition was rejected")
        if result["numeric_signature"] != selected["numeric_signature"] or result["array_signature"] != selected["array_signature"]:
            raise ContractError("final selected repetition differs from selected candidate")
        result.update(case_id=rep, case_dir=str(output / "final" / rep), group="final", coordinate="selected", step=0.0)
        final_rows.append(result)
    all_rows = rows + final_rows
    if len(all_rows) != MAX_EVALUATIONS:
        raise ContractError(f"objective evaluation count {len(all_rows)} differs from declared maximum {MAX_EVALUATIONS}")
    finite_difference_panel(all_rows, anchor, specs, output)
    write_tables(output, all_rows)
    summary = {"status": "verified_selection", "anchor_loss": ANCHOR_LOSS, "selected": selected,
               "rows": all_rows, "objective_evaluations": len(all_rows),
               "stationary_solves_actual": sum(int(r.get("stationary_solves") or 0) for r in all_rows if r.get("status") == "valid"),
               "stationary_solves_maximum": MAX_EVALUATIONS * MAX_STATIONARY_SOLVES,
               "target_fingerprint": plan.get("target_signature"), "objective_canonical_sha256": OBJECTIVE}
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("plan", "smoke", "production"), required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--selected-summary", type=Path, required=True)
    parser.add_argument("--smoke-summary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case-timeout", type=float, default=900.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--global-timeout", type=float, default=10800.0)
    parser.add_argument("--finish-utc", default="2026-09-21T11:30:00Z")
    parser.add_argument("--expected-plan-sha256")
    args = parser.parse_args()
    result = run_panel(args)
    print(json.dumps({k: v for k, v in result.items() if k not in {"rows", "probes", "anchor"}}, sort_keys=True))


if __name__ == "__main__":
    main()
