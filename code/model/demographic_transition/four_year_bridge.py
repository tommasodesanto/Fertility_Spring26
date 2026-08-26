"""Four-year interface between model birth flows and annual person cohorts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .person_cohort_law import CohortState, FlowLedger, advance_one_year


Array = np.ndarray


@dataclass(frozen=True)
class PersonBlockLedger:
    """Annual ledgers and block totals for one model period."""

    start_year: int
    end_year: int
    model_period_births: float
    annual_ledgers: tuple[FlowLedger, ...]
    person_identity_max_abs: float
    head_identity_max_abs: float
    total_net_migration: float
    total_new_heads_from_nonheads: float
    total_head_dissolutions: float
    total_net_migrant_heads: float


def advance_person_state_block(
    state: CohortState,
    *,
    model_period_births: float,
    period_years: int,
    birth_sex_shares: Mapping[int, Array],
    survival: Mapping[int, Array],
    net_migration: Mapping[int, Array],
    headship_rates: Mapping[int, Array] | Array,
    topcoded_last_age: bool = True,
    tolerance: float = 1e-8,
) -> tuple[CohortState, PersonBlockLedger]:
    """Spread a model-period birth flow uniformly and advance annual cohorts.

    The active model period is four years. Its expected birth count is a flow
    over that interval, so the initial bridge allocates one equal quarter to
    each calendar year while using that year's observed/projected sex ratio.
    This timing convention is explicit and can later be replaced by a measured
    within-period birth profile without changing the accounting boundary.
    """

    years = int(period_years)
    births_total = float(model_period_births)
    if years < 1:
        raise ValueError("period_years must be positive")
    if not np.isfinite(births_total) or births_total < 0.0:
        raise ValueError("model_period_births must be finite and nonnegative")
    current = state.validated(tolerance=tolerance)
    ledgers: list[FlowLedger] = []
    annual_total = births_total / years
    for year in range(current.year + 1, current.year + years + 1):
        if year not in birth_sex_shares:
            raise KeyError(f"missing birth sex shares for {year}")
        shares = np.asarray(birth_sex_shares[year], dtype=float).reshape(-1)
        if shares.shape != (current.persons.shape[0],):
            raise ValueError(
                f"birth sex shares for {year} have shape {shares.shape}; "
                f"expected {(current.persons.shape[0],)}"
            )
        if np.any(~np.isfinite(shares)) or np.any(shares < 0.0):
            raise ValueError(f"birth sex shares for {year} are invalid")
        if not np.isclose(float(np.sum(shares)), 1.0, atol=1e-12, rtol=0.0):
            raise ValueError(f"birth sex shares for {year} do not sum to one")
        if year not in survival or year not in net_migration:
            raise KeyError(f"missing survival or migration for {year}")
        rates = (
            headship_rates[year]
            if isinstance(headship_rates, Mapping)
            else headship_rates
        )
        current, ledger = advance_one_year(
            current,
            births=annual_total * shares,
            survival=survival[year],
            net_migration=net_migration[year],
            headship_rates=rates,
            topcoded_last_age=topcoded_last_age,
            tolerance=tolerance,
        )
        ledgers.append(ledger)

    return current, PersonBlockLedger(
        start_year=state.year,
        end_year=current.year,
        model_period_births=births_total,
        annual_ledgers=tuple(ledgers),
        person_identity_max_abs=max(
            abs(ledger.person_identity_residual) for ledger in ledgers
        ),
        head_identity_max_abs=max(
            abs(ledger.head_identity_residual) for ledger in ledgers
        ),
        total_net_migration=sum(
            float(np.sum(ledger.net_migration)) for ledger in ledgers
        ),
        total_new_heads_from_nonheads=sum(
            float(np.sum(ledger.new_heads_from_nonheads)) for ledger in ledgers
        ),
        total_head_dissolutions=sum(
            float(np.sum(ledger.head_dissolutions)) for ledger in ledgers
        ),
        total_net_migrant_heads=sum(
            float(np.sum(ledger.net_migrant_heads)) for ledger in ledgers
        ),
    )
