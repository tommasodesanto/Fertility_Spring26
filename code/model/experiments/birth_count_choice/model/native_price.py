"""Single at-price native Bellman+KFE with one counted lifecycle call."""
from __future__ import annotations

import copy
import signal
import time
from pathlib import Path

import numpy as np


class CaseDeadline(TimeoutError):
    pass


def _alarm(_signum, _frame):
    raise CaseDeadline("300-second per-case deadline")


def solve_fixed_price(context, d_bar, q, budget, label, stage_dir):
    from . import credit
    from .engine import solver
    from .credit import DEAD_MASS_TOL

    stage_dir = Path(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=False)
    budget.claim_lifecycle(label)
    started = time.monotonic()
    case_deadline_epoch = min(time.time() + 300.0, budget.deadline_epoch)
    P = copy.deepcopy(context["P"])
    credit.bind_engine_credit(P, "corrected", d_bar)
    P.native_inherited_distribution_evidence_dir = str(stage_dir / "inherited_state_failures")
    grid = context["b_grid"]
    sd = solver.precompute_shared(P, grid)
    prior_handler = signal.getsignal(signal.SIGALRM)
    remaining = case_deadline_epoch - time.time()
    if remaining <= 0:
        raise CaseDeadline("Global deadline reached before solve")
    signal.signal(signal.SIGALRM, _alarm)
    prior_timer = signal.setitimer(signal.ITIMER_REAL, remaining)
    try:
        sol = solver.solve_markov_income_at_prices(np.array([q], dtype=float), P, grid,
                                                  SD=sd, verbose=False, fast_stats=False)
    finally:
        signal.setitimer(signal.ITIMER_REAL, *prior_timer)
        signal.signal(signal.SIGALRM, prior_handler)
    elapsed = time.monotonic() - started
    if elapsed > 300 or time.time() > case_deadline_epoch:
        raise CaseDeadline("Case or global deadline overrun")
    if float(getattr(P, "_entry_censored_mass", 0.0)) > DEAD_MASS_TOL:
        raise RuntimeError("Inherited entry censoring would remove occupied mass")
    if not np.isfinite(float(sol.mean_parity)):
        raise RuntimeError("Nonfinite native solution")
    save_arrays = any(token in label for token in ("phase_a", "selected", "repeat", "final"))
    arrays = {}
    if save_arrays:
        arrays = {k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and v.dtype != object}
        arrays.update({"shared." + k: v for k, v in vars(sd).items()
                       if isinstance(v, np.ndarray) and v.dtype != object})
        np.savez_compressed(stage_dir / "solution_arrays.npz", **arrays)
    if time.time() > case_deadline_epoch:
        raise CaseDeadline("Case deadline exceeded during stage serialization")
    if "prepared" in context:
        tfr = float(context["prepared"].rt["chain"].extract_moments(sol, P)["tfr"])
    else:
        tfr = None
    summary = dict(label=label, d_bar=d_bar, q=q, elapsed_seconds=elapsed,
                   native_array_count=len(arrays), tfr=tfr,
                   own_rate=float(sol.own_rate),
                   entry_censored_mass=float(getattr(P, "_entry_censored_mass", 0.0)),
                   adult_entry_relative_gap=float(getattr(sol, "adult_entry_stationary_relative_gap", np.nan)),
                   stage_arrays=str(stage_dir / "solution_arrays.npz") if save_arrays else None)
    context["write_json"](stage_dir / "summary.json", summary)
    budget.progress("completed_case", summary=summary)
    return dict(sol=sol, sd=sd, P=P, b_grid=grid, price=np.array([q]),
                stage_dir=stage_dir, summary=summary,
                case_deadline_epoch=case_deadline_epoch)
