"""Borrowing limits and the entrant feasibility census.

Two explicitly named rules. Neither is a default: callers must pick one.

ReferenceCredit  the block0506 rule, restated from solver.py/parameters.py:
  renter   b' >= min(s_{j+1} * min(b, 0), -D_{j+1})       (taper rollover)
  buyer    b' >= m_{j+1} * (-phi * p * H')                 (purchase income on)
  stayer   DUE: b' >= max(min(b, m_{j+1}(-phi p H)), death floor) when
           native_due_stayer_credit; death floor -(1-psi) p H if death is
           possible at j (terminal age or survival < 1), else -inf.

CorrectedCredit is a historical September 29 original-timing fixture:
  renter   b' >= -d_bar; b' >= 0 at terminal age or positive death risk
  sale     fixture uses b + (1 - psi) p H >= 0 (raw, pre-clip)
  owners   unchanged from the reference rule.
The fixture classes and credit_rule are not production transaction-feasibility
callbacks. The adopted post-interest engine instead requires R*b + S >= 0
for owner-to-renter sale solvency, where S = (1 - psi) p H.
`bind_engine_credit` selects the active engine mode without invoking fixtures.
The byte-preserved engine provenance is recorded in engine_inventory.json;
native solve/root provenance is recorded in native_core_inventory.json.

The census reports, never repairs: no truncation, transfer, deletion,
renormalization or positive-credit fallback.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

DEAD_VALUE_CUTOFF = -1e9   # solver.DEAD_VALUE_CUTOFF
DEAD_MASS_TOL = 1e-12      # solver.DEAD_MASS_TOL


def death_possible(P, j: int) -> bool:
    return j == int(P.J) - 1 or (bool(P.use_age_survival) and float(P.survival_probs[j]) < 1.0)


def collateral_floor(P, j: int, price: float, house: float) -> float:
    """Next-period financed-share floor with the (inert here) LTV multiplier."""
    m = np.asarray(P.owner_ltv_multipliers, dtype=float)
    return float(m[min(max(j + 1, 0), m.size - 1)]) * (-float(np.asarray(P.phi)[0]) * price * house)


def stayer_death_floor(P, j: int, price: float, house: float) -> float:
    return -(1.0 - float(P.psi)) * price * house if death_possible(P, j) else -np.inf


@dataclass(frozen=True)
class ReferenceCredit:
    """Historical block0506 original-timing fixture; not an engine callback."""
    name: str = "reference_block0506"

    def renter_floor(self, P, b, j: int) -> np.ndarray:
        s = float(np.asarray(P.debt_taper_weights)[j + 1])
        D = float(np.asarray(P.debt_caps)[j + 1])
        return np.minimum(s * np.minimum(np.asarray(b, dtype=float), 0.0), -D)

    def buyer_floor(self, P, j: int, price: float, house: float) -> float:
        if not bool(P.native_purchase_income):
            raise ValueError("Only the purchase-income reference contract is supported")
        return collateral_floor(P, j, price, house)

    def stayer_floor(self, P, b, j: int, price: float, house: float) -> np.ndarray:
        if not bool(P.native_due_stayer_credit):
            raise ValueError("Reference stayer rule requires native_due_stayer_credit")
        return np.maximum(np.minimum(np.asarray(b, dtype=float), collateral_floor(P, j, price, house)),
                          stayer_death_floor(P, j, price, house))

    def sale_allowed(self, P, b, price: float, house: float) -> np.ndarray:
        return np.ones(np.shape(b), dtype=bool)


@dataclass(frozen=True)
class CorrectedCredit(ReferenceCredit):
    """Historical original-timing fixture, including its uncapitalized sale test."""
    d_bar: float = 0.0
    name: str = "corrected_explicit_d_bar"

    def __post_init__(self):
        if not (np.isfinite(self.d_bar) and self.d_bar >= 0.0):
            raise ValueError("d_bar must be finite and nonnegative")

    def renter_floor(self, P, b, j: int) -> np.ndarray:
        floor = 0.0 if death_possible(P, j) else -float(self.d_bar)
        return np.full(np.shape(b), floor, dtype=float)

    def sale_allowed(self, P, b, price: float, house: float) -> np.ndarray:
        """Original-timing fixture only; active production tests R*b + S."""
        return np.asarray(b, dtype=float) + (1.0 - float(P.psi)) * price * house >= 0.0


def credit_rule(mode: str, d_bar: float | None = None):
    """Explicit selection; no fallback between modes."""
    if mode == "reference":
        if d_bar is not None:
            raise ValueError("Reference mode takes no d_bar")
        return ReferenceCredit()
    if mode == "corrected":
        if d_bar is None:
            raise ValueError("Corrected mode requires an explicit d_bar")
        return CorrectedCredit(d_bar=float(d_bar))
    raise ValueError("Unknown credit mode: " + mode)


def entrant_feasibility_census(entry_mass: np.ndarray, entry_value: np.ndarray, b_grid: np.ndarray) -> dict:
    """Positive-mass entrant nodes whose value is Bellman-dead.

    entry_mass, entry_value: arrays over identical entrant states with the
    wealth node on axis 0. Returns every dead occupied node; `feasible` is
    False when their mass exceeds DEAD_MASS_TOL. Nothing is altered.
    """
    mass = np.asarray(entry_mass, dtype=float)
    value = np.asarray(entry_value, dtype=float)
    if mass.shape != value.shape or mass.shape[0] != len(b_grid):
        raise ValueError("Census arrays must share shape with wealth on axis 0")
    dead = (mass > 0.0) & ~(value > DEAD_VALUE_CUTOFF)
    rows = [dict(index=[int(k) for k in idx], wealth=float(b_grid[idx[0]]), mass=float(mass[idx]))
            for idx in zip(*np.nonzero(dead))]
    total = float(mass[dead].sum())
    return dict(feasible=total <= DEAD_MASS_TOL, dead_mass=total,
                dead_share=total / float(mass.sum()) if mass.sum() > 0 else float("nan"),
                cells=rows)


def bind_engine_credit(P, mode: str, d_bar: float | None = None):
    """Set the engine's credit field for an explicit mode; returns P.

    reference: `unsecured_credit_limit` must be absent or None (legacy
    rollover/taper rule, bit-for-bit). corrected: validation copied from the
    reviewed upstream overlay `parameters.py` (lines 508-520): scalar,
    finite, nonnegative, float-cast, and incompatible with
    native_solvency_credit. The engine then applies b' >= -D for renters, the
    zero floor under death risk, and the post-interest raw sale gate
    R_gross*b + (1-psi) p H >= 0 before grid clipping. Its source hashes are
    recorded in engine_inventory.json, not a materialize_receipt file.
    """
    if mode == "reference":
        if d_bar is not None or getattr(P, "unsecured_credit_limit", None) is not None:
            raise ValueError("Reference mode requires unsecured_credit_limit None")
        return P
    if mode != "corrected" or d_bar is None:
        raise ValueError("Corrected mode requires an explicit d_bar")
    raw_credit = d_bar
    if not np.isscalar(raw_credit):
        raise ValueError("unsecured_credit_limit must be None or a finite non-negative scalar")
    try:
        value = float(raw_credit)
    except (TypeError, ValueError) as exc:
        raise ValueError("unsecured_credit_limit must be None or a finite non-negative scalar") from exc
    if not np.isfinite(value) or value < 0.0:
        raise ValueError("unsecured_credit_limit must be None or a finite non-negative scalar")
    if bool(getattr(P, "native_solvency_credit", False)):
        raise ValueError("unsecured_credit_limit cannot be combined with native_solvency_credit")
    P.unsecured_credit_limit = value
    return P
