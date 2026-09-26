"""Small, auditable observation operators for the first-birth bridge.

This module maps an event-year interview onto the housing date in the model's
four-year grid. It does not estimate an event study, simulate households, or
load model code/checkpoints. ``birth_offset_years`` is an explicit observation
convention measured from the model decision date; it is not estimated.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import floor, isfinite
from typing import Mapping


@dataclass(frozen=True)
class ObserverBridgeContract:
    """Pinned definitions for a descriptive, education-free auxiliary bridge."""

    empirical_design: str = "PSID H-v2"
    model_period_years: int = 4
    model_age_start: int = 18
    model_age_max: int = 82
    age_cell_width_years: int = 4
    baseline_event_years: tuple[int, int] = (-3, -2)
    outcome_event_years: tuple[int, int] = (3, 4)
    main_birth_offset_years: float = 0.0
    midpoint_sensitivity_years: float = 2.0
    uniform_sensitivity: str = "seeded household-level U[0,4); separate from main rule"
    within_period_housing: str = "held at the model-period choice"
    education_rule: str = "omit in both data bridge and model; auxiliary only"
    data_weights: str = "PSID longitudinal IW"
    model_weights: str = "model probability mass"
    adoption_status: str = "not adopted; does not replace H-v2 or experimental 1.465"

    def validate(self) -> None:
        if self.model_period_years <= 0 or self.age_cell_width_years <= 0:
            raise ValueError("period length and age-cell width must be positive")
        if self.model_age_start < 0 or self.model_age_max < self.model_age_start:
            raise ValueError("model age bounds must be ordered and nonnegative")
        if (self.model_age_max - self.model_age_start) % self.age_cell_width_years:
            raise ValueError("model age bounds must align with the age-cell grid")
        if not 0 <= self.main_birth_offset_years < self.model_period_years:
            raise ValueError("main birth offset must lie within the model period")
        if not 0 <= self.midpoint_sensitivity_years < self.model_period_years:
            raise ValueError("midpoint offset must lie within the model period")
        for window in (self.baseline_event_years, self.outcome_event_years):
            if len(window) != 2 or window[0] > window[1]:
                raise ValueError("event windows must be ordered two-endpoint ranges")


def model_period_for_event_year(
    birth_period: int,
    event_year: int,
    birth_offset_years: float,
    period_years: int = 4,
) -> int:
    """Return the model housing period observed at an integer event year.

    Birth occurs ``birth_offset_years`` after the beginning of ``birth_period``.
    Housing is held constant within each model period. Integer event year 0 is
    the calendar birth year, so negative event years can map to the prior model
    period. ``floor`` (not truncation toward zero) is essential for that case.
    """

    if isinstance(birth_period, bool) or not isinstance(birth_period, int):
        raise TypeError("birth_period must be an integer model-period index")
    if isinstance(event_year, bool) or not isinstance(event_year, int):
        raise TypeError("event_year must be an integer number of calendar years")
    if isinstance(period_years, bool) or not isinstance(period_years, int):
        raise TypeError("period_years must be an integer")
    if period_years <= 0:
        raise ValueError("period_years must be positive")
    if not isfinite(birth_offset_years) or not 0 <= birth_offset_years < period_years:
        raise ValueError("birth_offset_years must be finite and in [0, period_years)")
    return birth_period + floor((event_year + birth_offset_years) / period_years)


def model_age_cell_lower_bound(
    age_years: int,
    age_start: int = 18,
    cell_width_years: int = 4,
    age_max: int = 82,
) -> int:
    """Map an empirical age to the lower endpoint of its model-aligned cell."""

    if any(isinstance(x, bool) or not isinstance(x, int) for x in (age_years, age_start, cell_width_years, age_max)):
        raise TypeError("age and cell parameters must be integers")
    if cell_width_years <= 0:
        raise ValueError("cell_width_years must be positive")
    if age_years < age_start or age_years > age_max:
        raise ValueError("age_years must lie within the model's adult age support")
    return age_start + cell_width_years * ((age_years - age_start) // cell_width_years)


def observe_model_rooms(
    rooms_by_period: Mapping[int, float],
    birth_period: int,
    event_year: int,
    birth_offset_years: float,
    period_years: int = 4,
) -> float:
    """Read housing at the model date corresponding to one event-time row."""

    model_period = model_period_for_event_year(
        birth_period, event_year, birth_offset_years, period_years
    )
    if model_period not in rooms_by_period:
        raise KeyError(f"model history does not contain housing at period {model_period}")
    return float(rooms_by_period[model_period])
