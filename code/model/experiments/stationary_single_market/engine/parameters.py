"""Supported-path parameters; provenance and exact changes are in source_provenance.json/review.diff."""
from __future__ import annotations
import copy
import math
from types import SimpleNamespace
from typing import Any, Mapping
import numpy as np


def unsecured_debt_floor(current_unsecured: Any, s_next: float, D_next: float) -> np.ndarray:
    """Return the next-period unsecured floor for a current unsecured position."""

    u = np.asarray(current_unsecured, dtype=float)
    return np.minimum(float(s_next) * np.minimum(u, 0.0), -float(D_next))


def get_fecundity_by_age(P: SimpleNamespace) -> np.ndarray:
    """Per-period conception probability by age index j (length J).

    omega1 == 0 -> all ones (production behavior, including beyond the
    terminal age): this exact rule is the bitwise-nesting guarantee.
    """
    J = int(P.J)
    w1 = float(getattr(P, "fecundity_omega1", 0.0))
    if w1 == 0.0:
        return np.ones(J, dtype=float)
    w2 = float(getattr(P, "fecundity_omega2", 0.0))
    terminal = float(getattr(P, "fecundity_terminal_age", 45.0))
    ages = float(P.age_start) + np.arange(J, dtype=float) * float(P.da)
    pi = 1.0 - w1 * np.exp(w2 * (ages - float(P.age_start)))
    pi = np.clip(pi, 0.0, 1.0)
    terminal_decay = float(getattr(P, "fecundity_terminal_decay", 0.0))
    if terminal_decay < 0.0 or not np.isfinite(terminal_decay):
        raise ValueError("fecundity_terminal_decay must be finite and nonnegative.")
    if terminal_decay > 0.0:
        tail_start = float(getattr(P, "fecundity_tail_start_age", 40.0))
        if not np.isfinite(tail_start):
            raise ValueError("fecundity_tail_start_age must be finite.")
        pi *= np.exp(-terminal_decay * np.maximum(ages - tail_start, 0.0))
    pi[ages >= terminal] = 0.0
    return pi


def readiness_gate_active(P: SimpleNamespace) -> bool:
    """Whether the default-off E6c childless readiness state is active."""
    return bool(getattr(P, "readiness_gate_enabled", False))


def readiness_cumulative_probability(P: SimpleNamespace, age: float) -> float:
    """Unconditional probability that readiness has arrived by ``age``."""
    location = float(getattr(P, "readiness_location_age", 14.0))
    spread = float(getattr(P, "readiness_spread_years", 2.0))
    if not np.isfinite(location):
        raise ValueError("readiness_location_age must be finite.")
    if not np.isfinite(spread) or spread <= 0.0:
        raise ValueError("readiness_spread_years must be finite and positive.")
    x = float(np.clip((float(age) - location) / spread, -40.0, 40.0))
    return float(1.0 / (1.0 + np.exp(-x)))


def readiness_transition_hazard(
    P: SimpleNamespace,
    current_age: float,
    next_age: float,
) -> float:
    """Conditional unsettled-to-settled probability over one age interval."""
    if float(next_age) < float(current_age):
        raise ValueError("readiness transition ages must be weakly increasing.")
    current = readiness_cumulative_probability(P, current_age)
    nxt = readiness_cumulative_probability(P, next_age)
    return float(np.clip((nxt - current) / max(1.0 - current, 1e-14), 0.0, 1.0))


def readiness_childless_states(P: SimpleNamespace) -> tuple[int, ...]:
    """Child-state indices representing childless households."""
    return (0, 1) if readiness_gate_active(P) else (0,)


def readiness_settled_state(P: SimpleNamespace) -> int:
    """State from which first-child entry is available."""
    return 1 if readiness_gate_active(P) else 0


def independent_child_maturation_active(P: SimpleNamespace) -> bool:
    """This engine supports independent child counts only."""
    return True


def child_earnings_penalty_active(P: SimpleNamespace) -> bool:
    """Whether any children-at-home earnings penalty entry is nonzero."""
    return bool(np.any(np.asarray(getattr(P, "child_earnings_penalty", np.zeros(4)), dtype=float) != 0.0))


def mortgage_stay_floor_active(P: SimpleNamespace) -> bool:
    """Whether stayer-specific mortgage floors must be computed."""
    if bool(getattr(P, "mortgage_origination_only", False)):
        return True
    return float(getattr(P, "mortgage_amortization", 0.0)) > 0.0


def rental_wedge_active(P: SimpleNamespace) -> bool:
    """Whether the size-dependent rental wedge changes renter budgets."""
    return (
        float(getattr(P, "rental_wedge_intercept", 0.0)) != 0.0
        or float(getattr(P, "rental_wedge_slope", 0.0)) != 0.0
    )


def estate_receiver_active(P: SimpleNamespace) -> bool:
    """Whether aggregate estates are paid as lump sums to ages 45-65."""
    return str(getattr(P, "estate_receiver", "none")).strip().lower() == "ages_45_65"


def bequest_utility_net_active(P: SimpleNamespace) -> bool:
    """Whether the bequest-utility estate is valued net of the selling cost.

    True when the utility-side boolean is set or when the estate transfer is
    on (the transfer values estates net in both utility and accounting).
    """
    if bool(getattr(P, "bequest_net_of_selling_cost", False)):
        return True
    return estate_receiver_active(P)


def estate_housing_value(P: SimpleNamespace, price: float, rooms: float, *, for_accounting: bool) -> float:
    """Housing leg of an estate: gross, or net of the selling cost.

    The utility side nets out the selling cost when ``bequest_utility_net``
    holds; the aggregate-flow accounting nets it out when the estate transfer
    is on. Otherwise the gross value nests the current model bit for bit.
    """
    gross = float(price) * float(rooms)
    net = bool(for_accounting and estate_receiver_active(P)) or bool(
        (not for_accounting) and bequest_utility_net_active(P)
    )
    if net:
        return float(1.0 - float(getattr(P, "psi", 0.0))) * gross
    return gross


def estate_transfer_at_age(P: SimpleNamespace, j: int) -> float:
    """Per-household estate transfer at age index j (0 outside 45-65)."""
    if not estate_receiver_active(P):
        return 0.0
    transfer = float(getattr(P, "estate_lump_sum_transfer", 0.0))
    if transfer == 0.0:
        return 0.0
    age = float(getattr(P, "age_start", 18.0)) + float(int(j)) * float(getattr(P, "da", 4.0))
    if 45.0 <= age <= 65.0:
        return transfer
    return 0.0


def children_at_home_count(nn: int, cs: int, P: SimpleNamespace) -> int:
    """Children currently at home, bounded by children ever born and three."""
    return int(min(max(int(cs), 0), int(nn), 3))


def child_earnings_multiplier(P: SimpleNamespace, j: int, m: int) -> float:
    """After-tax earnings multiplier for age index j and m children at home.

    Working ages scale by (1 - penalty); retirement ages are untouched, as
    are pensions, the payroll-tax base, and the PAYGO balance (those use
    unpenalized earnings by construction).
    """
    if int(j) >= int(getattr(P, "J_R", 0)):
        return 1.0
    penalty = np.asarray(getattr(P, "child_earnings_penalty", np.zeros(4)), dtype=float).reshape(-1)
    if penalty.size == 0:
        return 1.0
    return float(1.0 - penalty[min(max(int(m), 0), penalty.size - 1)])


def parent_age_maturation_active(P: SimpleNamespace) -> bool:
    """Whether the default-off parent-age maturation law is active."""
    return str(getattr(P, "child_maturation_mode", "constant")).strip().lower() == "parent_age"
