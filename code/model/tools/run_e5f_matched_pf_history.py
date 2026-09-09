"""Conditional historical/person PF composition; no equilibrium or calibration.

The caller supplies the full dated price/preference/transfer paths, a terminal
value boundary, and frozen 2023 person/head primitives. No roots or empirical
targets are inferred here. The 2023 household decision is evaluated only once.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Callable, Sequence

import numpy as np

import run_e5f_perfect_foresight_transition as pf
import run_e5f_perfect_foresight_person_demography as person_pf


@dataclass
class ConditionalHistoryEvaluation:
    history: pf.PathEvaluation
    person_tail: person_pf.PersonPathEvaluation
    rows: list[dict[str, Any]]
    values: list[np.ndarray]
    bellman_solves: int
    initial_2023_age_head_gap: float
    scope: str = "Conditional PF path with supplied terminal boundary; not a converged equilibrium"


def evaluate_history_and_person_tail(
    *, years: Sequence[int], prices: Sequence[float], psi_path: Sequence[float],
    transfer_path: Sequence[float], terminal_price: float, terminal_V: np.ndarray,
    base_parameters: Any, b_grid: np.ndarray, initial_state: pf.PFInitialState,
    historical_conditioning: pf.HistoricalConditioning,
    initial_2023_persons: person_pf.CohortState,
    demographic_primitives: person_pf.AnnualDemographicPrimitives,
    supply_rule: Any, birth_to_entry_conversion: float,
    observer: Callable | None = None,
) -> ConditionalHistoryEvaluation:
    """Evaluate 2007--2019 history followed by the person tail from 2023.

    The optional observer receives the global date index, current evaluation,
    dated parameters, wealth grid, and shared inputs before each advancement.
    """
    dates = np.asarray(years)
    p = np.asarray(prices, dtype=float)
    psi = np.asarray(psi_path, dtype=float)
    transfers = np.asarray(transfer_path, dtype=float)
    if (dates.ndim != 1 or len(dates) < 5
            or not np.array_equal(dates, 2007 + 4 * np.arange(len(dates)))
            or float(base_parameters.period_years) != 4.):
        raise ValueError("Joined PF evaluation requires complete four-year dates from 2007 through at least 2023")
    if (any(a.shape != dates.shape or not np.isfinite(a).all() for a in (p, psi, transfers))
            or np.any(p <= 0) or np.any(transfers < 0)):
        raise ValueError("Explicit price, preference and nonnegative transfer paths must match all dates")
    if historical_conditioning.start_year != 2007:
        raise ValueError("Joined history must begin in 2007")
    historical_conditioning.validate(base_parameters, 4, initial_state, birth_to_entry_conversion)
    if observer is not None and not callable(observer):
        raise ValueError("Joined observer must be callable")
    if (observer is not None and historical_conditioning.observer is not None
            and historical_conditioning.observer is not observer):
        raise ValueError("Supply one common observer, not different historical and joined observers")
    callback = observer if observer is not None else historical_conditioning.observer
    people = initial_2023_persons.validated()
    frozen_people = demographic_primitives.initial_person_state
    if (people.year != 2023 or frozen_people.year != 2023
            or not np.array_equal(people.persons, frozen_people.persons)
            or not np.array_equal(people.heads, frozen_people.heads)):
        raise ValueError("Initial 2023 persons/heads must match the supplied frozen demographic primitives")

    tail_prices, tail_psi, tail_transfers = p[4:], psi[4:], transfers[4:]
    tail_rents = pf.rents_from_asset_prices(tail_prices, terminal_price, base_parameters)
    tail_values, tail_backward_solves = pf.backward_value_path(
        prices=tail_prices, rents=tail_rents, psi_path=tail_psi,
        terminal_V=terminal_V, base_parameters=base_parameters,
        b_grid=b_grid, transfer_path=tail_transfers,
    )
    history = pf.evaluate_path_at_prices(
        prices=p[:4], psi_path=psi[:4], transfer_path=transfers[:4],
        terminal_price=float(p[4]), terminal_V=tail_values[0],
        base_parameters=base_parameters, b_grid=b_grid, initial_state=initial_state,
        supply_rule=supply_rule, birth_to_entry_conversion=birth_to_entry_conversion,
        historical_conditioning=replace(historical_conditioning, observer=callback),
    )
    g_2023 = history.terminal_state.g_pre
    heads_2023 = person_pf.aggregate_heads_to_model_age_cells(
        people, age_start=int(base_parameters.age_start),
        cell_width=int(base_parameters.da), number_of_cells=int(base_parameters.J),
    )
    age_gap = float(np.max(np.abs(g_2023.sum(axis=(0, 1, 2, 4, 5, 6)) - heads_2023)))
    if not np.isfinite(age_gap) or age_gap > 2e-9:
        raise RuntimeError(f"Historical 2023 household/person head-age identity fails: {age_gap}")

    def tail_observer(period, evaluation, parameters, grid, shared):
        callback(period + 4, evaluation, parameters, grid, shared)

    tail = person_pf.evaluate_path_at_prices_person_demography(
        prices=tail_prices, psi_path=tail_psi, transfer_path=tail_transfers,
        terminal_price=terminal_price, terminal_V=terminal_V,
        base_parameters=base_parameters, b_grid=b_grid,
        initial_state=person_pf.PersonPFState(g_pre=g_2023.copy(), persons=people),
        demographic_primitives=demographic_primitives, supply_rule=supply_rule,
        precomputed_value_path=tail_values,
        observer=tail_observer if callback is not None else None,
    )
    rows = [dict(row) for row in history.rows]
    rows.extend(dict(row, period=int(row['period']) + 4) for row in tail.rows)
    if [row['calendar_year'] for row in rows] != dates.tolist():
        raise RuntimeError("Joined path duplicated or omitted a calendar date")
    count = history.bellman_solves + tail_backward_solves + tail.bellman_solves
    if count != 2 * len(dates):
        raise RuntimeError("Joined path repeated or omitted a backward/forward household solve")
    return ConditionalHistoryEvaluation(
        history=history, person_tail=tail, rows=rows,
        values=history.values[:-1] + tail.values, bellman_solves=count,
        initial_2023_age_head_gap=age_gap,
    )
