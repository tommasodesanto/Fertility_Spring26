"""Isolated financing-path binding for the existing native dated observer.

The caller supplies a source-authenticated, freshly repeated native runtime.
This module never constructs a calibration reference or changes a model default.
"""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import threading
import numpy as np

import dated_phi

_LOCK = threading.RLock()


def financing_path(kind: str, periods: int) -> np.ndarray:
    """One surprise date, or an announced permanent financing expansion."""
    if periods < 1:
        raise ValueError("At least one date is required")
    if kind == "control":
        return np.full(periods, 0.8)
    if kind == "temporary":
        path = np.full(periods, 0.8)
        path[0] = 1.0
        return path
    if kind == "permanent":
        return np.ones(periods)
    raise ValueError("Unknown financing experiment")


@contextmanager
def bind_phi_path(pf, phi_path):
    """Install the dated engine only while the reviewed mapping is evaluated.

    The existing mapping keeps its household, fiscal, estate, population, and
    policy-array audits.  Its exact policy cache hashes the date-specific P.phi.
    """
    phi = np.asarray(phi_path, dtype=float).reshape(-1).copy()
    if pf is not dated_phi.base:
        raise RuntimeError("Dated financing and native policy cache use different PF modules")
    if not len(phi) or not np.isfinite(phi).all() or np.any((phi < 0) | (phi > 1)):
        raise ValueError("Invalid financing path")
    with _LOCK:
        original = pf.evaluate_path_at_prices
        def evaluate_with_phi(**kwargs):
            if len(np.asarray(kwargs["prices"]).reshape(-1)) != len(phi):
                raise RuntimeError("Financing and price paths differ")
            return dated_phi.evaluate_path_at_prices(phi_path=phi, **kwargs)
        pf.evaluate_path_at_prices = evaluate_with_phi
        try:
            yield
        finally:
            pf.evaluate_path_at_prices = original


def map_case(runtime, *, kind, terminal, endpoint, prices, pensions, folder,
             start_year=2007, initial_state=None):
    """Run one unchanged audited mapping from the authenticated baseline state."""
    prices = np.asarray(prices, dtype=float).reshape(-1)
    pensions = np.asarray(pensions, dtype=float).reshape(-1)
    phi = financing_path(kind, len(prices))
    if len(pensions) != len(prices):
        raise ValueError("Pension and price paths differ")
    terminal_phi = float(np.asarray(terminal["parameters"].phi).reshape(-1)[0])
    required_terminal = 1.0 if kind == "permanent" else 0.8
    if terminal_phi != required_terminal:
        raise RuntimeError("Terminal financing rule differs from dated policy")
    if float(np.asarray(runtime.P.phi).reshape(-1)[0]) != 0.8:
        raise RuntimeError("Initial fitted baseline must use phi=0.8")
    if not runtime.reference_verified:
        raise RuntimeError("Native initial-state repeat has not passed")
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    psi = np.full(len(prices), float(runtime.P.psi_child))
    with bind_phi_path(runtime.pf, phi):
        native, record = runtime.mapping(
            terminal, endpoint, prices, pensions, psi, folder,
            initial_state=runtime.initial_state if initial_state is None else initial_state,
            start_year=start_year,
        )
    if not all(record["gates"].values()):
        raise RuntimeError("Existing dated mapping gate failed")
    if [float(row["phi"]) for row in native.rows] != list(phi):
        raise RuntimeError("Dated financing audit differs from requested path")
    record["phi_path"] = phi.tolist()
    record["case_kind"] = kind
    return native, record
