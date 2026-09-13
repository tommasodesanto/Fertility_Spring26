"""Run the authorized six-mapping, horizon-two original-queue SSJ smoke.

This is an isolated measurement bridge.  It calls the pinned original queue
operator and reconstructs its *existing* three scaled residual equations; it
does not solve a transition or alter a scientific kernel.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace as NS

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--contract", type=Path, required=True)
    ap.add_argument("--seconds", type=float, default=540.0)
    args = ap.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("Refusing to overwrite native-smoke output")
    contract = json.loads(args.contract.read_text())
    for name, expected in contract["sha256"].items():
        path = Path(name)
        if sha(path) != expected:
            raise ValueError("Pinned input changed: " + str(path))
    started = time.monotonic()
    deadline = started + min(float(args.seconds), 540.0)
    args.output.mkdir(parents=True, exist_ok=True)
    # Source is placed first only to make this connection driver importable;
    # load_context subsequently pins the scientific helpers and verifies hashes.
    sys.path.insert(0, str(args.spec.parent / "source"))
    import run_e5f_original_queue_experiments as runner

    c = runner.load_context(args.spec)
    c.spec_path = args.spec
    rebated = c.rebated
    save(args.output / "started.json", dict(
        status="load_context_passed", spec_sha256=sha(args.spec),
        runner_module=str(Path(runner.__file__).resolve()),
        rebated_module=str(Path(rebated.__file__).resolve()),
        queue_module=str(Path(c.queue.__file__).resolve()),
        horizon=2, mapping_cap=6, production_eligible=False))
    if time.time() >= float(c.spec["absolute_deadline_unix"]):
        raise TimeoutError("Frozen experiment deadline has passed")
    q = float(c.packet["evaluation"].policy.price[0])
    pension = float(c.old.parameters.pension)
    rebate = float(c.old.parameters.property_tax_lump_sum_transfer)
    terminal = NS(parameters=c.old.parameters, policy=c.packet["evaluation"].policy, asset_price=q)
    inherited = rebated.InheritedState(2007, c.old.initial_state)
    base = np.vstack((np.full(2, np.log(q)), np.full(2, pension), np.full(2, rebate)))
    direction = np.zeros((3, 2)); direction[0, 0] = 1.0
    h = 1e-5
    evaluations: list[dict] = []
    baseline_state = None
    stationary_drift = None
    replay_state_gap = None

    def state_gap(actual, reference) -> dict:
        actual_g = np.asarray(actual.g_pre, dtype=float)
        reference_g = np.asarray(reference.g_pre, dtype=float)
        scale = max(float(reference_g.sum()), 1e-15)
        actual_q = np.asarray(actual.scheduled_entries, dtype=float)
        reference_q = np.asarray(reference.scheduled_entries, dtype=float)
        actual_raw = np.asarray(actual.scheduled_raw_entries, dtype=float)
        reference_raw = np.asarray(reference.scheduled_raw_entries, dtype=float)
        return dict(distribution_relative_l1=float(np.abs(actual_g-reference_g).sum()/scale),
            population_relative_gap=float(abs(actual_g.sum()/scale-1.0)),
            queue_relative_max=float(np.max(np.abs(actual_q/reference_q-1.0))),
            raw_queue_relative_max=float(np.max(np.abs(actual_raw/reference_raw-1.0))))

    def mapping(label: str, unknowns: np.ndarray) -> np.ndarray:
        nonlocal baseline_state, stationary_drift, replay_state_gap
        if len(evaluations) >= 6 or time.monotonic() >= deadline:
            raise TimeoutError("Native mapping cap or time budget reached")
        log_q, b, r = unknowns
        result = c.queue.queue_path(inherited=inherited, old_state=c.old,
            prices=np.exp(log_q), pensions=b, transfers=r, psi=float(c.old.parameters.psi_child),
            terminal=terminal, demographics=None)
        if len(result.values) != 3 or len(result.rows) != 2:
            raise ValueError("Native horizon-two path did not retain T+1 values and T rows")
        blocks = []
        for row in result.rows:
            blocks.append(rebated.dated_residual(
                demand=row["housing_demand"], supply=row["housing_supply"],
                payroll_accounts=row, tax_accounts=row))
        residual = rebated.stack_dated_residuals(blocks)
        gates = dict(maximum_mass_accounting_error=float(result.maximum_mass_accounting_error),
            maximum_policy_reproduction_error=float(result.maximum_policy_reproduction_error),
            maximum_feasibility_projection_mass=float(result.maximum_feasibility_projection_mass))
        state = result.person_tail.terminal_state
        if label == "baseline":
            stationary_drift = state_gap(state, c.old.initial_state)
            baseline_state = state
        elif label == "exact_baseline_replay":
            if baseline_state is None:
                raise RuntimeError("Exact replay evaluated before baseline state")
            replay_state_gap = state_gap(state, baseline_state)
            baseline_state = None
        evaluations.append(dict(label=label, residual=residual.tolist(), max_abs=float(np.max(np.abs(residual))),
            rows=result.rows, native_g_pre_mass=float(result.person_tail.terminal_state.g_pre.sum()),
            scheduled_entries=np.asarray(result.person_tail.terminal_state.scheduled_entries).tolist(),
            scheduled_raw_entries=np.asarray(result.person_tail.terminal_state.scheduled_raw_entries).tolist(), gates=gates))
        save(args.output / "latest_completed.json", dict(status="mapping_complete", mappings=len(evaluations),
            latest_label=label, latest_max_abs=float(np.max(np.abs(residual))),
            elapsed_seconds=time.monotonic()-started))
        return residual

    try:
        with c.queue.original_queue_adapter(), c.cache.policy_cache(c.joined.pf, max_bytes=6 * 1024**3):
            f0 = mapping("baseline", base)
            f0_replay = mapping("exact_baseline_replay", base)
            fp_h = mapping("plus_h", base + h * direction)
            fm_h = mapping("minus_h", base - h * direction)
            fp_half = mapping("plus_h_over_2", base + 0.5 * h * direction)
            fm_half = mapping("minus_h_over_2", base - 0.5 * h * direction)
        baseline_gap = float(np.max(np.abs(f0)))
        replay_gap = float(np.max(np.abs(f0 - f0_replay)))
        deriv_h = (fp_h - fm_h) / (2 * h)
        deriv_half = (fp_half - fm_half) / h
        if stationary_drift is None or replay_state_gap is None:
            raise RuntimeError("Baseline/replay state receipts are missing")
        save(args.output / "native_smoke.json", dict(status="complete", horizon=2, mappings=len(evaluations),
            unknowns="(log_house_price, pension, rebate)", direction="log_house_price[0]",
            h=h, baseline_max_abs=baseline_gap, exact_replay_max_abs=replay_gap,
            derivative_h=deriv_h.tolist(), derivative_h_over_2=deriv_half.tolist(),
            derivative_max_abs_difference=float(np.max(np.abs(deriv_h-deriv_half))),
            baseline_gate_passed=baseline_gap <= 2e-4, exact_replay_passed=replay_gap <= 2e-10,
            baseline_stationary_state_drift=stationary_drift,
            baseline_stationary_state_drift_limit=1e-5,
            baseline_stationary_state_drift_passed=all(v <= 1e-5 for v in stationary_drift.values()),
            exact_replay_state_gap=replay_state_gap,
            exact_replay_state_gap_passed=all(v <= 2e-10 for v in replay_state_gap.values()),
            elapsed_seconds=time.monotonic()-started, evaluations=evaluations,
            production_eligible=False, fast_news_claimed=False))
    except BaseException as exc:
        save(args.output / "native_smoke_failure.json", dict(error_type=type(exc).__name__, error=str(exc),
            mappings=len(evaluations), elapsed_seconds=time.monotonic()-started, evaluations=evaluations))
        raise


if __name__ == "__main__":
    main()
