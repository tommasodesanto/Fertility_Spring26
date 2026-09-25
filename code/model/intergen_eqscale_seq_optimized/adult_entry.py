"""Aggregate birth-cohort entry, separate from parental child dependency.

One model period is four years. ``births_t`` counts children born during date t,
including the extra children represented by the 3+ top bin. The queue returns
potential entrant *households* for date t+1: half of cohort t-3 (age 16 at
entry) and half of cohort t-4 (age 20). Retention, if any, acts after this
queue; parental survival and the child-departure transition never enter it.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


REPLACEMENT_FERTILITY = 2.1


def adjusted_births(raw_births: float, top_bin_entries: float,
                    top_bin_weight: float, top_state: int = 3) -> float:
    """Return child units after adding the excess represented by the top bin."""
    raw = float(raw_births)
    top = float(top_bin_entries)
    weight = float(top_bin_weight)
    if not all(map(math.isfinite, (raw, top, weight))):
        raise ValueError("Birth accounting inputs must be finite")
    if raw < 0 or top < 0 or top > raw + 1e-12 or weight < top_state:
        raise ValueError("Invalid birth or top-bin accounting inputs")
    return raw + (weight - top_state) * top


def potential_entry_households(adjusted_birth_children: float) -> float:
    """Apply the retained birth-to-household conversion exactly once."""
    births = float(adjusted_birth_children)
    if not math.isfinite(births) or births < 0:
        raise ValueError("Adjusted births must be finite and nonnegative")
    return births / REPLACEMENT_FERTILITY


@dataclass(frozen=True)
class SplitBirthEntryQueue:
    """Pending potential household flows, before any geographic retention.

    The first tuple has three waiting slots and the second four. At date t,
    popping each first slot supplies date t+1 entry. Appending date-t births
    therefore makes their effects arrive at dates t+4 and t+5, respectively.
    """

    due_in_16: tuple[float, float, float]
    due_in_20: tuple[float, float, float, float]

    def __post_init__(self) -> None:
        if len(self.due_in_16) != 3 or len(self.due_in_20) != 4:
            raise ValueError("Split entry requires three and four waiting slots")
        if any(not math.isfinite(x) or x < 0 for x in self.due_in_16 + self.due_in_20):
            raise ValueError("Pending cohort flows must be finite and nonnegative")

    @classmethod
    def constant_prehistory(cls, adjusted_birth_children: float) -> "SplitBirthEntryQueue":
        half = 0.5 * potential_entry_households(adjusted_birth_children)
        return cls((half,) * 3, (half,) * 4)

    @property
    def stock(self) -> float:
        return sum(self.due_in_16) + sum(self.due_in_20)

    def step(self, adjusted_birth_children: float) -> tuple[float, "SplitBirthEntryQueue"]:
        """Return potential entry at t+1 and the queue after date-t births."""
        half = 0.5 * potential_entry_households(adjusted_birth_children)
        due = self.due_in_16[0] + self.due_in_20[0]
        next_queue = SplitBirthEntryQueue(
            self.due_in_16[1:] + (half,),
            self.due_in_20[1:] + (half,),
        )
        return due, next_queue


def require_closed_stationary_renewal(entry_households: float,
                                      adjusted_birth_children: float,
                                      fertility_tolerance: float) -> dict[str, float]:
    """Gate B=E using the existing completed-fertility tolerance in child units.

    Call this only after the fertility-intercept normalization has completed.
    Since B = births/2.1, the equivalent gate is
    ``abs(births/E - 2.1) <= fertility_tolerance``. No outside entrants or
    retention coefficient are introduced by the closed stationary model.
    """
    E = float(entry_households)
    tolerance = float(fertility_tolerance)
    B = potential_entry_households(adjusted_birth_children)
    if not math.isfinite(E) or E <= 0 or not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("Invalid closed-renewal gate inputs")
    completed = float(adjusted_birth_children) / E
    gap = completed - REPLACEMENT_FERTILITY
    if abs(gap) > tolerance:
        raise ValueError(
            f"Closed birth-entry renewal failed: births/E={completed:.12g}, "
            f"target={REPLACEMENT_FERTILITY:.12g}, tolerance={tolerance:.12g}"
        )
    return {"entry_E": E, "potential_B": B, "births_per_entry": completed,
            "fertility_gap": gap, "entry_residual": E - B}
