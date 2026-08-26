"""Couple annual person cohorts to the four-year household distribution."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .four_year_bridge import PersonBlockLedger, advance_person_state_block
from .household_head_bridge import (
    HeadMassRebalanceLedger,
    aggregate_heads_to_model_age_cells,
    rebalance_household_distribution_by_age,
)
from .person_cohort_law import CohortState


Array = np.ndarray


@dataclass(frozen=True)
class HouseholdPersonBlockLedger:
    """Exact accounting records for one coupled model-period transition."""

    person: PersonBlockLedger
    household_heads: HeadMassRebalanceLedger
    target_head_mass: float
    coupled_household_mass: float
    household_person_head_gap: float


def advance_household_person_block(
    raw_next_household_distribution: Array,
    person_state: CohortState,
    *,
    model_period_births: float,
    period_years: int,
    birth_sex_shares: Mapping[int, Array],
    survival: Mapping[int, Array],
    net_migration: Mapping[int, Array],
    headship_rates: Mapping[int, Array] | Array,
    model_age_start: int,
    model_age_cell_width: int,
    number_of_model_age_cells: int,
    age_axis: int = 3,
    empty_age_templates: Mapping[int, Array] | None = None,
    topcoded_last_age: bool = True,
    tolerance: float = 1e-8,
) -> tuple[Array, CohortState, HouseholdPersonBlockLedger]:
    """Advance persons, then impose their head stocks on household ages.

    The household transition supplied here determines the conditional economic
    state distribution within each age cell. The person law determines the
    cell's total household-head mass. Rescaling is explicit because the current
    model does not yet contain separate economic states for head formation,
    dissolution, or migration.
    """

    next_person_state, person_ledger = advance_person_state_block(
        person_state,
        model_period_births=model_period_births,
        period_years=period_years,
        birth_sex_shares=birth_sex_shares,
        survival=survival,
        net_migration=net_migration,
        headship_rates=headship_rates,
        topcoded_last_age=topcoded_last_age,
        tolerance=tolerance,
    )
    target_age_mass = aggregate_heads_to_model_age_cells(
        next_person_state,
        age_start=model_age_start,
        cell_width=model_age_cell_width,
        number_of_cells=number_of_model_age_cells,
    )
    next_households, head_ledger = rebalance_household_distribution_by_age(
        raw_next_household_distribution,
        target_age_mass,
        age_axis=age_axis,
        empty_age_templates=empty_age_templates,
        tolerance=tolerance,
    )
    target_total = float(np.sum(target_age_mass))
    household_total = float(np.sum(next_households))
    gap = household_total - target_total
    scale = max(1.0, target_total)
    if abs(gap) > tolerance * scale:
        raise RuntimeError(
            "Coupled household mass does not equal person-law head mass: "
            f"gap={gap:.6g}"
        )
    return next_households, next_person_state, HouseholdPersonBlockLedger(
        person=person_ledger,
        household_heads=head_ledger,
        target_head_mass=target_total,
        coupled_household_mass=household_total,
        household_person_head_gap=gap,
    )
