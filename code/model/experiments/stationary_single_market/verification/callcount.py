"""Verification-only call counter (sys.setprofile; no function replacement).

Counts calls and inclusive wall time of named Python functions, matched by
code object. Used identically for the lab and old engine in the GE benchmark,
so both headline times carry the same overhead. Numba kernels are not Python
frames and are not counted. Never enabled by default on the normal path.
"""
from __future__ import annotations

import signal
import sys
import time


class BudgetExceeded(RuntimeError):
    """Raised BEFORE an at-price solve beyond the call budget starts, or at the
    solve-stage deadline. A budget stop is a failed case, never a certificate."""

STAGES = ("solve_markov_income_at_prices", "solve_bellman_full_markov_income",
          "forward_distribution_markov_income", "upgrade_fast_markov_solution",
          "pack_solution_markov_income", "pack_fast_solution_markov_income",
          "refine_one_market_markov_income", "solve_markov_income_equilibrium",
          "attach_markov_market_accounting")


class CallCounter:
    def __init__(self, module, *, max_lifecycle=None, deadline_seconds=None):
        self.max_lifecycle, self.deadline_seconds = max_lifecycle, deadline_seconds
        self.started = self.ended = None
        self.stop_reason = None
        self.codes = {getattr(module, n).__code__: n for n in STAGES if hasattr(module, n)}
        self.calls = {n: 0 for n in self.codes.values()}
        self.seconds = {n: 0.0 for n in self.codes.values()}
        self.events = []          # (name, 'call'|'return', monotonic time), in order
        self._open = {}

    def _profile(self, frame, event, arg):
        name = self.codes.get(frame.f_code)
        if name is None:
            return
        now = time.perf_counter()
        if event == "call" and self.deadline_seconds is not None and now - self.started > self.deadline_seconds:
            self.stop_reason = f"solve-stage deadline {self.deadline_seconds} s reached before {name}"
            raise BudgetExceeded(self.stop_reason)
        if (event == "call" and name == "solve_markov_income_at_prices" and self.max_lifecycle is not None
                and self.calls[name] >= self.max_lifecycle):
            self.stop_reason = f"lifecycle budget {self.max_lifecycle} reached; call {self.calls[name] + 1} not started"
            raise BudgetExceeded(self.stop_reason)
        if event == "call":
            self.calls[name] += 1
            self._open.setdefault(name, []).append(now)
            self.events.append((name, "call", now))
        elif event == "return":
            start = self._open.get(name, [now]).pop() if self._open.get(name) else now
            self.seconds[name] += now - start
            self.events.append((name, "return", now))

    def _alarm(self, signum, frame):
        self.stop_reason = f"solve-stage deadline {self.deadline_seconds} s (hard alarm)"
        raise BudgetExceeded(self.stop_reason)

    def __enter__(self):
        self.started = time.perf_counter()
        if self.deadline_seconds is not None:   # backstop between counted calls (fires after a running kernel returns)
            signal.signal(signal.SIGALRM, self._alarm)
            signal.setitimer(signal.ITIMER_REAL, float(self.deadline_seconds))
        sys.setprofile(self._profile)
        return self

    def __exit__(self, *exc):
        sys.setprofile(None)
        self.ended = time.perf_counter()
        if self.deadline_seconds is not None:
            signal.setitimer(signal.ITIMER_REAL, 0.0)
            signal.signal(signal.SIGALRM, signal.SIG_DFL)

    def summary(self) -> dict:
        full = self.calls.get("solve_markov_income_at_prices", 0)
        return dict(instrumentation="sys.setprofile, verification-only, identical for both engines",
                    budget=dict(max_lifecycle=self.max_lifecycle, deadline_seconds=self.deadline_seconds,
                                stop_reason=self.stop_reason),
                    profiled_solve_stage_seconds=(self.ended - self.started) if self.started and self.ended else None,
                    calls=self.calls, inclusive_seconds=self.seconds,
                    lifecycle_evaluations=full,
                    bellman_calls=self.calls.get("solve_bellman_full_markov_income", 0),
                    kfe_calls=self.calls.get("forward_distribution_markov_income", 0),
                    final_payload_upgrades=self.calls.get("upgrade_fast_markov_solution", 0),
                    final_solution=self.final_solution(),
                    note="lifecycle_evaluations counts every at-price solve (fast and full); "
                         "unique_fast_price_evaluations from the engine already includes refinement entries.")

    def final_solution(self) -> str:
        """Resolve the GE's final full solution from the ordered event log."""
        tail = [(n, e) for n, e, _ in self.events
                if n in ("upgrade_fast_markov_solution", "solve_markov_income_at_prices", "refine_one_market_markov_income")]
        after = []
        for n, e in reversed(tail):
            if n == "refine_one_market_markov_income" and e == "return":
                break
            after.append((n, e))
        names = {n for n, e in after if e == "call"}
        if "upgrade_fast_markov_solution" in names:
            return "payload_upgrade (no additional household solve)"
        if "solve_markov_income_at_prices" in names:
            return "full_household_replay (one additional at-price solve)"
        return "unresolved"
