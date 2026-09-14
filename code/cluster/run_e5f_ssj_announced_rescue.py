"""Experiment-only warm restart of the 104-date announced root with a supplied Jacobian.

Reuses the announced batch's own driver module (``run_e5f_announced_original_queue``)
for context loading, the fresh endpoint check, announced preference routing,
the exact ``run_path`` root call, plots and receipts.  The only injected
changes are (i) the start guess, taken from the finished announced root's
exactly reproduced best iterate, (ii) the Broyden ``initial_jacobian``
(``broyden_final`` from that root's receipt, or ``toeplitz`` extrapolated from
the ten-date measured lag profiles), and (iii) optionally the isolated
direction-preserving step rule.  Nothing in production changes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import threading
import time
from unittest.mock import patch

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from e5f_ssj_toeplitz_jacobian import assemble_from_receipt
from e5f_ssj_scaled_step_root import solve_price_path_scaled


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=lambda v: v.tolist()) + "\n")
    tmp.replace(path)


def select_jacobian(mode, announced_receipt, toeplitz_receipt, horizon):
    """Return ``(matrix, description)`` for the requested warm-start Jacobian."""
    n = 3 * horizon
    if mode == "broyden_final":
        J = np.asarray(announced_receipt["final_jacobian"], dtype=float)
        info = dict(mode=mode, source="announced root receipt final Broyden matrix")
    elif mode == "toeplitz":
        J, info = assemble_from_receipt(toeplitz_receipt, horizon)
        info["mode"] = mode
    else:
        raise ValueError("jacobian_mode must be 'broyden_final' or 'toeplitz'")
    if J.shape != (n, n) or not np.isfinite(J).all():
        raise ValueError("Warm-start Jacobian has the wrong shape or is not finite")
    info["condition_number"] = float(np.linalg.cond(J))
    return J, info


def warm_coordinates(announced_receipt, horizon):
    best = announced_receipt.get("best")
    if best is None or announced_receipt.get("final_reproduction_max_abs") != 0.0:
        raise ValueError("Warm start requires an exactly reproduced announced best iterate")
    x = np.asarray(best["prices"], dtype=float)
    if x.shape != (3 * horizon,) or not np.isfinite(x).all() or np.any(x <= 0):
        raise ValueError("Announced best iterate has the wrong shape")
    return x.reshape(3, horizon), float(best["score"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()
    m = read(args.manifest)
    for path, digest in m["file_sha256"].items():
        if sha(path) != digest:
            raise ValueError("Changed pinned input: " + path)
    announced = read(m["announced_manifest"])
    announced_receipt = read(m["announced_root_receipt"])
    toeplitz_receipt = read(m["toeplitz_receipt"]) if m.get("toeplitz_receipt") else None
    out = Path(m["output"])
    folder = out / "run"
    if folder.exists():
        raise ValueError("Refusing to overwrite a started rescue")
    sys.path.insert(0, str(Path(m["announced_source_dir"])))
    import run_e5f_announced_original_queue as ann
    spec = read(announced["spec"])
    sys.path.insert(0, str(Path(spec["batch"]) / "source"))
    import run_e5f_original_queue_experiments as runner
    c = runner.load_context(announced["spec"])
    c.spec_path = Path(announced["spec"])
    endpoint, receipt = ann.endpoint_from_manifest(announced)
    horizon = 104
    guess, previous_best = warm_coordinates(announced_receipt, horizon)
    jacobian, jacobian_info = select_jacobian(m["jacobian_mode"], announced_receipt, toeplitz_receipt, horizon)
    folder.mkdir(parents=True)
    end = time.time() + float(m["seconds"])
    deadline = time.monotonic() + (end - time.time())
    save(folder / "experiment_contract.json", dict(
        manifest_sha256=sha(args.manifest), horizon=horizon, warm_start_from=m["announced_root_receipt"],
        previous_best_score=previous_best, previous_evaluations=announced_receipt.get("evaluations"),
        jacobian=jacobian_info, step_rule=m.get("step_rule", "clipped"), seconds=m["seconds"], deadline_unix=end,
        root_mapping_cap=8, announced_deadline_not_applied=True, production_eligible=False,
        fake_news_derivatives_constructed=False))
    save(folder / "warm_start_jacobian.json", dict(jacobian=jacobian.tolist(), **jacobian_info))
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(60):
            save(folder / "controller_heartbeat.json", dict(remaining_seconds=end - time.time()))
            if time.time() >= end:
                save(folder / "controller_failure.json", dict(error="Rescue deadline"))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    try:
        ann.fresh_endpoint_check(c, endpoint, receipt, folder, min(deadline, time.monotonic() + 600))
        psi_path = np.r_[announced["psi_levels"], np.full(100, ann.FINAL_PSI)]
        evaluate, first = ann.announced_queue_path(c, psi_path)
        native_solve = c.rebated.solve_rebated_forecast

        def solve_with_jacobian(**kwargs):
            controls = dict(kwargs["root_controls"])
            if controls.get("max_evaluations") != 8:
                raise ValueError("Rescue keeps the retained eight-mapping budget")
            controls["initial_jacobian"] = jacobian
            return native_solve(**dict(kwargs, root_controls=controls))
        step_solver = solve_price_path_scaled if m.get("step_rule", "clipped") == "scaled" else c.rebated._path_root_solver()
        with c.queue.original_queue_adapter(), patch.object(c.rebated, "evaluate_forecast", evaluate), \
                patch.object(c.rebated, "first_period_state", first), \
                patch.object(c.rebated, "_path_root_solver", lambda: step_solver), \
                patch.object(c.rebated, "solve_rebated_forecast", solve_with_jacobian), \
                c.cache.policy_cache(c.joined.pf, max_bytes=12 * 1024**3):
            result = ann.run_path(c, endpoint, announced, folder, deadline, guess=guess)
        history = [dict(evaluation=e["evaluation"], phase=e["phase"], score=e["score"], evaluation_seconds=e.get("evaluation_seconds"))
                   for e in result.root_receipt.get("history", []) if "evaluation" in e]
        previous = [dict(evaluation=e["evaluation"], phase=e["phase"], score=e["score"], evaluation_seconds=e.get("evaluation_seconds"))
                    for e in announced_receipt.get("history", []) if "evaluation" in e]
        save(folder / "comparison.json", dict(previous_history=previous, rescue_history=history,
             previous_best_score=previous_best, rescue_best_score=(result.root_receipt.get("best") or {}).get("score"),
             rescue_converged=bool(result.root_receipt.get("finite_horizon_market_fiscal_converged")),
             jacobian=jacobian_info, step_rule=m.get("step_rule", "clipped"), production_eligible=False))
        save(folder / "controller_complete.json", dict(completed=True,
             finite_horizon_market_fiscal_converged=bool(result.root_receipt.get("finite_horizon_market_fiscal_converged")),
             stationary_endpoint_verified=True, horizon_verified=False, production_eligible=False))
    except BaseException as exc:
        save(folder / "controller_failure.json", dict(error_type=type(exc).__name__, error=str(exc)))
        raise
    finally:
        stop.set()


if __name__ == "__main__":
    main()
