"""Bounded local income refinement around one verified native anchor.

This controller owns scheduling and receipts only.  The income adapter and all
scientific score gates remain in ``run_e5f_income_candidate_search``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import signal
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable

import run_e5f_income_candidate_search as old

PARAMETERS = old.PARAMETERS
OBJECTIVE = old.OBJECTIVE
SEED = 20260921
DEFAULT_ANCHOR_LOSS = 326.9831988727637
STOP = old.STOP


def read(path: Path) -> Any:
    return json.loads(Path(path).read_text())


def write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                  allow_nan=False).encode()).hexdigest()


def numeric_fit(receipt: dict[str, Any]) -> dict[str, Any]:
    return {"loss": receipt.get("loss"), "price": receipt.get("_native_price"),
            "parameters": [(x.get("parameter"), x.get("estimate")) for x in receipt.get("parameters", [])],
            "target_fit": [(x.get("restriction_id", x.get("target")), x.get("target"),
                            x.get("model"), x.get("gap"), x.get("actual_weight"),
                            x.get("loss_contribution")) for x in receipt.get("target_fit", [])],
            "normalization": (receipt.get("normalization", {}).get("completed_fertility"),
                              receipt.get("normalization", {}).get("psi_child"))}


def anchor_receipt(payload: dict[str, Any]) -> dict[str, Any]:
    """Accept a native score JSON or an enriched selected/verification wrapper."""
    for key in ("receipt", "selected", "verification"):
        if isinstance(payload.get(key), dict) and ("loss" not in payload or key != "receipt"):
            try:
                return anchor_receipt(payload[key])
            except ValueError:
                pass
    if "loss" not in payload or "parameters" not in payload:
        raise ValueError("anchor score must contain native loss and parameter rows")
    return payload


def structural_parameters(receipt: dict[str, Any]) -> dict[str, float]:
    out = {r["parameter"]: float(r["estimate"]) for r in receipt.get("parameters", [])
           if r.get("structural_coordinate")}
    if set(out) != set(PARAMETERS):
        raise old.ContractError("anchor lacks all nine structural coordinates")
    return out


def bounds(plan: dict[str, Any]) -> dict[str, tuple[float, float]]:
    raw = plan.get("parameter_bounds") or plan.get("pilot_box")
    if not isinstance(raw, dict) or set(raw) != set(PARAMETERS):
        raise old.ContractError("complete parameter bounds required")
    return {n: (float(raw[n][0]), float(raw[n][1])) for n in PARAMETERS}


def proposals(anchor: dict[str, float], plan: dict[str, Any], count: int = 20) -> list[dict[str, float]]:
    """Five seeded antithetic directions at 5% and 2% normalized steps."""
    import random
    box = bounds(plan); rng = random.Random(SEED)
    scales = (0.05, 0.02)
    seen = {digest(anchor)}; out = []
    # Positive coordinates use relative steps; beta uses small absolute steps.
    width = {n: (max(abs(anchor[n]), 1e-12) if n not in ("beta_annual",)
                 else .02) for n in PARAMETERS}
    directions = []
    for _ in range(5):
        v = [rng.uniform(-1., 1.) for _ in PARAMETERS]
        directions.extend((v, [-x for x in v]))
    for scale in scales:
        for direction in directions:
            point = {n: float(anchor[n]) + scale * width[n] * u
                     for n, u in zip(PARAMETERS, direction)}
            point = {n: min(box[n][1], max(box[n][0], point[n])) for n in PARAMETERS}
            key = digest(point)
            if key not in seen:
                seen.add(key); out.append(point)
            if len(out) == count:
                return out
    return out


def adapter_evaluator(plan_path: Path, output: Path, parameters: dict[str, float], **kwargs: Any) -> dict[str, Any]:
    return old.adapter_evaluator(plan_path, output, parameters, **kwargs)


def run_refinement(plan: dict[str, Any], plan_path: Path, output: Path,
                   anchor: dict[str, Any], evaluator: Callable[..., dict[str, Any]],
                   *, anchor_psi: float, seconds: float = 1200., workers: int = 10) -> dict[str, Any]:
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    started = time.monotonic(); deadline = started + seconds
    STOP.clear()
    params = structural_parameters(anchor); old.validate_score_contract(anchor, plan)
    anchor_loss = old.extract_score(anchor)
    design = proposals(params, plan, 20)
    write(output / "design.json", {"schema": "e5f_income_local_refinement_design_v1", "seed": SEED,
                                    "anchor_parameters": params, "anchor_loss": anchor_loss,
                                    "proposals": design, "bounds": bounds(plan), "max_proposals": 20,
                                    "steps": [0.05, 0.02], "step_convention": "relative to anchor except beta absolute .001/.0004; clipped to frozen bounds"})
    best = {"status": "anchor", "objective": anchor_loss, "parameters": params,
            "psi": anchor_psi, "receipt": anchor}
    write(output / "best_so_far.json", best)
    write(output / "latest_completed.json", {"status": "anchor_smoke_pending"})
    # The watchdog ensures a native child sees STOP even while a worker is in I/O.
    def watchdog() -> None:
        while not STOP.is_set():
            wait = min(30., max(0., deadline - time.monotonic()))
            if STOP.wait(wait):
                return
            write(output / "heartbeat.json", {"status": "running", "updated": time.time(),
                                                "elapsed_seconds": time.monotonic() - started})
            if time.monotonic() >= deadline:
                STOP.set(); return
    threading.Thread(target=watchdog, daemon=True).start()
    try:
        remaining = max(0., deadline - time.monotonic())
        if remaining <= 0:
            raise TimeoutError("no budget for anchor smoke")
        smoke = evaluator(params, output / "anchor_smoke", psi=anchor_psi, repetitions=2,
                          timeout=min(600., remaining), case_id="income_local_refinement_anchor_smoke")
        old.validate_score_contract(smoke, plan); old.validate_parameters(smoke, params)
        smoke_loss = old.extract_score(smoke)
        exact = (smoke.get("_summary", {}).get("exact_loss_equality") is True and smoke.get("_summary", {}).get("repetitions") == 2)
        same = numeric_fit(anchor) == numeric_fit(smoke)
        rel = abs(smoke_loss - anchor_loss) / max(1., abs(anchor_loss))
        smoke_gate = exact and same and rel <= 1e-8
        write(output / "anchor_smoke.json", {"loss": smoke_loss, "relative_loss_gap": rel,
                                               "exact_loss_equality": exact, "numeric_fit_equal": same,
                                               "gate": smoke_gate, "receipt": smoke})
        if not smoke_gate:
            raise old.ContractError("anchor smoke gate failed")
        completed = []
        for start in range(0, len(design), max(1, workers)):
            remaining = deadline - time.monotonic()
            if remaining <= 1.0 or STOP.is_set():
                break
            batch = design[start:start + workers]
            def one(item: tuple[int, dict[str, float]]) -> dict[str, Any]:
                index, point = item; case = output / "cases" / f"case_{index:02d}"; case.parent.mkdir(parents=True, exist_ok=True)
                try:
                    receipt = evaluator(point, case, psi=anchor_psi, repetitions=1,
                                        timeout=max(1., min(remaining, 900.)),
                                        case_id=f"income_local_refinement_{index:02d}")
                    old.validate_score_contract(receipt, plan); old.validate_parameters(receipt, point)
                    return {"case": index, "status": "completed", "objective": old.extract_score(receipt),
                            "parameters": point, "receipt": receipt}
                except old.ContractError:
                    STOP.set(); raise
                except (TimeoutError, RuntimeError) as exc:
                    detail = str(exc).lower()
                    forbidden = ('fingerprint', 'hash mismatch', 'source', 'accounting', 'budget', 'validation', 'contract', 'preflight')
                    expected = ('infeasible', 'did not converge', 'failed to converge', 'bracket failure', 'root solve', 'nonfinite', 'candidate exceeded')
                    if any(word in detail for word in forbidden) or not (isinstance(exc, TimeoutError) or STOP.is_set() or any(word in detail for word in expected)):
                        STOP.set()
                        raise
                    return {"case": index, "status": "numerical_rejected", "parameters": point,
                            "error": str(exc)}
                except Exception:
                    STOP.set(); raise
            with ThreadPoolExecutor(max_workers=min(workers, len(batch))) as pool:
                for future in as_completed([pool.submit(one, (start + i + 1, p)) for i, p in enumerate(batch)]):
                    row = future.result(); completed.append(row)
                    write(output / "latest_completed.json", row); write(output / "cases.json", completed)
                    if row["status"] == "completed" and row["objective"] < best["objective"]:
                        best = {"status": "improved", **{k: row[k] for k in ("case", "objective", "parameters", "receipt")}, "psi": anchor_psi}
                        write(output / "best_so_far.json", best)
        verification = {"status": "not_run"}
        if time.monotonic() + 600. <= deadline and not STOP.is_set():
            try:
                check = evaluator(best["parameters"], output / "selected_verification", psi=anchor_psi,
                                  repetitions=2, timeout=min(600., deadline - time.monotonic()),
                                  case_id="income_local_refinement_selected_verification")
                old.validate_score_contract(check, plan); old.validate_parameters(check, best["parameters"])
                verification = {"status": "verified" if old.extract_score(check) == best["objective"] and
                                check.get("_summary", {}).get("exact_loss_equality") is True and
                                numeric_fit(check) == numeric_fit(best["receipt"]) else "mismatch",
                                "receipt": check}
            except Exception as exc:
                verification = {"status": "failed", "error": str(exc)}
        result = {"schema": "e5f_income_local_refinement_v1",
                  "status": "verified_selection" if verification.get("status") == "verified" else "provisional",
                  "objective_canonical_sha256": OBJECTIVE, "anchor_loss": anchor_loss,
                  "anchor_psi": anchor_psi, "selected": best, "verification": verification,
                  "cases": completed, "proposal_count": len(completed), "max_proposals": len(design),
                  "workers": workers, "elapsed_seconds": time.monotonic() - started,
                  "time_budget_seconds": seconds, "bounds": bounds(plan),
                  "interpretation": "bounded local numerical refinement; no economic claim",
                  "anchor_comparison": "exact native loss and numeric-fit equality; relative tolerance is cross-platform documentation only"}
        write(output / "summary.json", result)
        return result
    except BaseException as exc:
        write(output / "failure.json", {"status": "failed", "error": repr(exc),
                                         "elapsed_seconds": time.monotonic() - started})
        raise
    finally:
        STOP.set()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--anchor-score", type=Path, required=True)
    parser.add_argument("--anchor-psi", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=1200.)
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("NUMBA_NUM_THREADS", "1")
    plan = read(args.plan); anchor = anchor_receipt(read(args.anchor_score))
    old.validate_score_contract(anchor, plan)
    result = run_refinement(plan, args.plan, args.output, anchor,
                            lambda p, c, **kw: adapter_evaluator(args.plan, c, p, **kw),
                            anchor_psi=args.anchor_psi, seconds=args.seconds, workers=args.workers)
    if result["status"] not in ("verified_selection", "provisional"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
