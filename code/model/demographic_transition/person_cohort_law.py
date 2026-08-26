"""Exact annual person-to-household-head cohort accounting.

The active household model indexes adult household heads.  This module keeps
the demographic objects that produce those heads separate and explicit:

* births create age-zero persons;
* age-specific survival advances every existing cohort;
* an age-by-sex net-migration residual enters after survival;
* nonheads form households and existing heads dissolve households; and
* household-head mass is an accounting stock, not ``births / 2.1``.

The headship adjustment is deliberately an accounting closure rather than a
behavioral household-formation model.  Within every age-sex cell it uses the
minimum one-way gross flow required to attain an externally supplied headship
rate: either nonheads become heads or heads cease to be heads, never both in
the same cell and year.  This makes the missing margin visible and produces an
exact flow ledger that can later be replaced by an estimated hazard model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np


Array = np.ndarray


def _finite_array(value: Array, name: str) -> Array:
    array = np.asarray(value, dtype=float)
    if array.ndim != 2:
        raise ValueError(f"{name} must have shape (sex, age); got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _rate_array(value: Array, name: str, shape: tuple[int, int]) -> Array:
    array = _finite_array(value, name)
    if array.shape != shape:
        raise ValueError(f"{name} has shape {array.shape}; expected {shape}")
    if np.any((array < 0.0) | (array > 1.0)):
        raise ValueError(f"{name} must lie in [0, 1]")
    return array


@dataclass(frozen=True)
class CohortState:
    """Resident persons and household heads by sex and single year of age."""

    year: int
    persons: Array
    heads: Array

    def validated(self, *, tolerance: float = 1e-10) -> "CohortState":
        persons = _finite_array(self.persons, "persons")
        heads = _finite_array(self.heads, "heads")
        if heads.shape != persons.shape:
            raise ValueError(
                f"heads has shape {heads.shape}; expected {persons.shape}"
            )
        if np.min(persons) < -tolerance:
            raise ValueError("persons contains negative mass")
        if np.min(heads) < -tolerance:
            raise ValueError("heads contains negative mass")
        if np.max(heads - persons) > tolerance:
            raise ValueError("household heads cannot exceed resident persons")
        return CohortState(
            year=int(self.year),
            persons=np.maximum(persons, 0.0).copy(),
            heads=np.maximum(heads, 0.0).copy(),
        )


@dataclass(frozen=True)
class FlowLedger:
    """One-year demographic and household-formation flows."""

    year: int
    births: Array
    newborn_deaths: Array
    existing_person_deaths: Array
    net_migration: Array
    existing_head_deaths: Array
    newborn_heads: Array
    new_heads_from_nonheads: Array
    head_dissolutions: Array
    net_migrant_heads: Array
    person_identity_residual: float
    head_identity_residual: float
    headship_mapping_residual: float


@dataclass(frozen=True)
class DemographicSteadyState:
    """Stationary person and head stocks under a births-per-head rule."""

    persons: Array
    heads: Array
    annual_births: float
    annual_births_per_head: float
    birth_head_multiplier: float
    fixed_migration_head_base: float
    renewal_ratio: float
    population_fixed_point_residual: float
    headship_mapping_residual: float
    response_iterations: int


@dataclass(frozen=True)
class DemographicSteadyStateComponents:
    """Linear stationary responses to one birth and fixed net migration."""

    persons_per_annual_birth: Array
    persons_from_fixed_migration: Array
    birth_head_multiplier: float
    fixed_migration_head_base: float
    response_iterations: int


def cohort_survivors(
    persons: Array,
    births: Array,
    survival: Array,
    *,
    topcoded_last_age: bool = True,
) -> tuple[Array, Array, Array]:
    """Advance persons one year before migration.

    ``survival[:, a]`` is the probability of reaching destination age ``a``.
    The last cell is interpreted as an open-ended top code when
    ``topcoded_last_age`` is true, so survivors from both the penultimate and
    the last cells enter the last cell.
    """

    stock = _finite_array(persons, "persons")
    if np.any(stock < 0.0):
        raise ValueError("persons contains negative mass")
    rates = _rate_array(survival, "survival", stock.shape)
    births_array = np.asarray(births, dtype=float).reshape(-1)
    if births_array.shape != (stock.shape[0],):
        raise ValueError(
            f"births has shape {births_array.shape}; expected {(stock.shape[0],)}"
        )
    if np.any(~np.isfinite(births_array)) or np.any(births_array < 0.0):
        raise ValueError("births must be finite and nonnegative")

    survivors = np.zeros_like(stock)
    survivors[:, 0] = births_array * rates[:, 0]
    if stock.shape[1] > 2:
        survivors[:, 1:-1] = stock[:, :-2] * rates[:, 1:-1]
    if stock.shape[1] == 1:
        if topcoded_last_age:
            survivors[:, 0] += stock[:, 0] * rates[:, 0]
    elif topcoded_last_age:
        survivors[:, -1] = (
            stock[:, -2] + stock[:, -1]
        ) * rates[:, -1]
    else:
        survivors[:, -1] = stock[:, -2] * rates[:, -1]

    newborn_deaths = births_array * (1.0 - rates[:, 0])
    existing_survivors = np.sum(survivors, axis=1) - survivors[:, 0]
    existing_deaths = np.sum(stock, axis=1) - existing_survivors
    scale = np.maximum(1.0, np.sum(stock, axis=1))
    if np.any(existing_deaths < -1e-12 * scale):
        raise RuntimeError("survival accounting created negative deaths")
    return survivors, newborn_deaths, np.maximum(existing_deaths, 0.0)


def _advance_existing_linear(
    stock: Array,
    survival: Array,
    *,
    topcoded_last_age: bool,
) -> Array:
    """Apply the existing-cohort survival operator to a possibly signed stock."""

    values = _finite_array(stock, "stock")
    rates = _rate_array(survival, "survival", values.shape)
    advanced = np.zeros_like(values)
    if values.shape[1] > 2:
        advanced[:, 1:-1] = values[:, :-2] * rates[:, 1:-1]
    if values.shape[1] == 1:
        if topcoded_last_age:
            advanced[:, 0] = values[:, 0] * rates[:, 0]
    elif topcoded_last_age:
        advanced[:, -1] = (values[:, -2] + values[:, -1]) * rates[:, -1]
    else:
        advanced[:, -1] = values[:, -2] * rates[:, -1]
    return advanced


def _stationary_linear_response(
    forcing: Array,
    survival: Array,
    *,
    topcoded_last_age: bool,
    tolerance: float,
    max_iterations: int,
) -> tuple[Array, int]:
    """Solve ``x = A x + forcing`` by contraction iteration."""

    forcing_array = _finite_array(forcing, "forcing")
    rates = _rate_array(survival, "survival", forcing_array.shape)
    current = np.zeros_like(forcing_array)
    for iteration in range(1, int(max_iterations) + 1):
        updated = _advance_existing_linear(
            current,
            rates,
            topcoded_last_age=topcoded_last_age,
        ) + forcing_array
        gap = float(np.max(np.abs(updated - current)))
        scale = max(1.0, float(np.max(np.abs(updated))))
        current = updated
        if gap <= float(tolerance) * scale:
            return current, iteration
    raise RuntimeError(
        "Stationary cohort response did not converge within "
        f"{int(max_iterations)} iterations"
    )


def stationary_demographic_components(
    *,
    birth_sex_shares: Array,
    survival: Array,
    net_migration: Array,
    headship_rates: Array,
    topcoded_last_age: bool = True,
    tolerance: float = 1e-12,
    max_iterations: int = 20_000,
) -> DemographicSteadyStateComponents:
    """Return the two linear components of the stationary cohort law."""

    rates = _finite_array(survival, "survival")
    rates = _rate_array(rates, "survival", rates.shape)
    migration = _finite_array(net_migration, "net_migration")
    headship = _rate_array(headship_rates, "headship_rates", rates.shape)
    if migration.shape != rates.shape:
        raise ValueError(
            f"net_migration has shape {migration.shape}; expected {rates.shape}"
        )
    shares = np.asarray(birth_sex_shares, dtype=float).reshape(-1)
    if shares.shape != (rates.shape[0],):
        raise ValueError(
            f"birth_sex_shares has shape {shares.shape}; expected {(rates.shape[0],)}"
        )
    if np.any(~np.isfinite(shares)) or np.any(shares < 0.0):
        raise ValueError("birth_sex_shares must be finite and nonnegative")
    if not np.isclose(float(np.sum(shares)), 1.0, rtol=0.0, atol=1e-12):
        raise ValueError("birth_sex_shares must sum to one")

    one_birth_forcing = np.zeros_like(rates)
    one_birth_forcing[:, 0] = shares * rates[:, 0]
    birth_response, birth_iterations = _stationary_linear_response(
        one_birth_forcing,
        rates,
        topcoded_last_age=topcoded_last_age,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )
    migration_response, migration_iterations = _stationary_linear_response(
        migration,
        rates,
        topcoded_last_age=topcoded_last_age,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )
    return DemographicSteadyStateComponents(
        persons_per_annual_birth=birth_response,
        persons_from_fixed_migration=migration_response,
        birth_head_multiplier=float(np.sum(birth_response * headship)),
        fixed_migration_head_base=float(np.sum(migration_response * headship)),
        response_iterations=max(birth_iterations, migration_iterations),
    )


def solve_demographic_steady_state(
    *,
    annual_births_per_head: float,
    birth_sex_shares: Array,
    survival: Array,
    net_migration: Array,
    headship_rates: Array,
    topcoded_last_age: bool = True,
    tolerance: float = 1e-12,
    max_iterations: int = 20_000,
) -> DemographicSteadyState:
    """Solve a stationary cohort law with fixed absolute net migration.

    Let ``A`` be the survival/aging operator, ``f`` the newborn-survivor
    vector generated by one annual birth, ``m`` fixed annual net migration,
    and ``h`` the headship-rate vector.  The person stock satisfies

    ``P = A P + f * b * (h'P) + m``,

    where ``b`` is annual births per household head.  The solution is obtained
    from the two linear responses ``(I-A)^-1 f`` and ``(I-A)^-1 m``.  A finite
    positive-migration steady state requires the renewal ratio
    ``b * h'(I-A)^-1 f`` to be strictly below one.
    """

    birth_rate = float(annual_births_per_head)
    if not np.isfinite(birth_rate) or birth_rate < 0.0:
        raise ValueError("annual_births_per_head must be finite and nonnegative")
    rates = _finite_array(survival, "survival")
    rates = _rate_array(rates, "survival", rates.shape)
    migration = _finite_array(net_migration, "net_migration")
    headship = _rate_array(headship_rates, "headship_rates", rates.shape)
    shares = np.asarray(birth_sex_shares, dtype=float).reshape(-1)
    components = stationary_demographic_components(
        birth_sex_shares=shares,
        survival=rates,
        net_migration=migration,
        headship_rates=headship,
        topcoded_last_age=topcoded_last_age,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )
    birth_response = components.persons_per_annual_birth
    migration_response = components.persons_from_fixed_migration
    birth_head_multiplier = components.birth_head_multiplier
    migration_head_base = components.fixed_migration_head_base
    renewal_ratio = birth_rate * birth_head_multiplier
    if renewal_ratio >= 1.0 - 10.0 * float(tolerance):
        raise RuntimeError(
            "No finite fixed-migration demographic steady state: "
            f"renewal_ratio={renewal_ratio:.12g}"
        )
    annual_births = birth_rate * migration_head_base / (1.0 - renewal_ratio)
    persons = migration_response + annual_births * birth_response
    heads = persons * headship
    if np.min(persons) < -1e-8 * max(1.0, float(np.max(np.abs(persons)))):
        raise RuntimeError(
            "The fixed migration vector implies negative stationary person mass"
        )
    persons = np.maximum(persons, 0.0)
    heads = np.maximum(heads, 0.0)

    survivor_stock = _advance_existing_linear(
        persons,
        rates,
        topcoded_last_age=topcoded_last_age,
    )
    newborns = np.zeros_like(persons)
    newborns[:, 0] = annual_births * shares * rates[:, 0]
    fixed_point_residual = float(
        np.max(np.abs(persons - survivor_stock - newborns - migration))
    )
    mapping_residual = float(np.max(np.abs(heads - persons * headship)))
    scale = max(1.0, float(np.max(persons)))
    if fixed_point_residual > 50.0 * float(tolerance) * scale:
        raise RuntimeError(
            "Demographic steady-state identity failed by "
            f"{fixed_point_residual:.6g}"
        )
    return DemographicSteadyState(
        persons=persons,
        heads=heads,
        annual_births=float(annual_births),
        annual_births_per_head=birth_rate,
        birth_head_multiplier=birth_head_multiplier,
        fixed_migration_head_base=migration_head_base,
        renewal_ratio=renewal_ratio,
        population_fixed_point_residual=fixed_point_residual,
        headship_mapping_residual=mapping_residual,
        response_iterations=components.response_iterations,
    )


def _surviving_heads(
    heads: Array,
    births: Array,
    survival: Array,
    headship_rates: Array,
    *,
    topcoded_last_age: bool,
) -> tuple[Array, Array, Array]:
    """Advance existing heads and report deaths before headship adjustment."""

    head_stock = _finite_array(heads, "heads")
    rates = _rate_array(survival, "survival", head_stock.shape)
    headship = _rate_array(headship_rates, "headship_rates", head_stock.shape)
    births_array = np.asarray(births, dtype=float).reshape(-1)

    survivors = np.zeros_like(head_stock)
    newborn_heads = births_array * rates[:, 0] * headship[:, 0]
    survivors[:, 0] = newborn_heads
    if head_stock.shape[1] > 2:
        survivors[:, 1:-1] = head_stock[:, :-2] * rates[:, 1:-1]
    if head_stock.shape[1] == 1:
        if topcoded_last_age:
            survivors[:, 0] += head_stock[:, 0] * rates[:, 0]
    elif topcoded_last_age:
        survivors[:, -1] = (
            head_stock[:, -2] + head_stock[:, -1]
        ) * rates[:, -1]
    else:
        survivors[:, -1] = head_stock[:, -2] * rates[:, -1]

    existing_survivors = np.sum(survivors, axis=1) - newborn_heads
    deaths = np.sum(head_stock, axis=1) - existing_survivors
    scale = np.maximum(1.0, np.sum(head_stock, axis=1))
    if np.any(deaths < -1e-12 * scale):
        raise RuntimeError("head-survival accounting created negative deaths")
    return survivors, newborn_heads, np.maximum(deaths, 0.0)


def advance_one_year(
    state: CohortState,
    *,
    births: Array,
    survival: Array,
    net_migration: Array,
    headship_rates: Array,
    topcoded_last_age: bool = True,
    tolerance: float = 1e-8,
) -> tuple[CohortState, FlowLedger]:
    """Advance a cohort state and return an exact flow ledger."""

    current = state.validated(tolerance=tolerance)
    shape = current.persons.shape
    migration = _finite_array(net_migration, "net_migration")
    if migration.shape != shape:
        raise ValueError(
            f"net_migration has shape {migration.shape}; expected {shape}"
        )
    headship = _rate_array(headship_rates, "headship_rates", shape)

    person_survivors, newborn_deaths, existing_person_deaths = cohort_survivors(
        current.persons,
        births,
        survival,
        topcoded_last_age=topcoded_last_age,
    )
    head_survivors, newborn_heads, existing_head_deaths = _surviving_heads(
        current.heads,
        births,
        survival,
        headship,
        topcoded_last_age=topcoded_last_age,
    )

    target_heads_before_migration = person_survivors * headship
    new_heads = np.maximum(target_heads_before_migration - head_survivors, 0.0)
    dissolutions = np.maximum(head_survivors - target_heads_before_migration, 0.0)
    migrant_heads = migration * headship

    next_persons = person_survivors + migration
    next_heads = (
        head_survivors + new_heads - dissolutions + migrant_heads
    )
    next_state = CohortState(
        year=current.year + 1,
        persons=next_persons,
        heads=next_heads,
    ).validated(tolerance=tolerance)

    births_array = np.asarray(births, dtype=float).reshape(-1)
    person_rhs = (
        float(np.sum(current.persons))
        + float(np.sum(births_array))
        - float(np.sum(newborn_deaths))
        - float(np.sum(existing_person_deaths))
        + float(np.sum(migration))
    )
    head_rhs = (
        float(np.sum(current.heads))
        - float(np.sum(existing_head_deaths))
        + float(np.sum(newborn_heads))
        + float(np.sum(new_heads))
        - float(np.sum(dissolutions))
        + float(np.sum(migrant_heads))
    )
    person_residual = float(np.sum(next_state.persons)) - person_rhs
    head_residual = float(np.sum(next_state.heads)) - head_rhs
    mapping_residual = float(
        np.max(np.abs(next_state.heads - next_state.persons * headship))
    )
    scale = max(1.0, float(np.sum(next_state.persons)))
    if abs(person_residual) > tolerance * scale:
        raise RuntimeError(
            f"person accounting identity failed by {person_residual:.6g}"
        )
    if abs(head_residual) > tolerance * scale:
        raise RuntimeError(f"head accounting identity failed by {head_residual:.6g}")
    if mapping_residual > tolerance * max(1.0, float(np.max(next_state.persons))):
        raise RuntimeError(
            f"headship mapping identity failed by {mapping_residual:.6g}"
        )

    return next_state, FlowLedger(
        year=next_state.year,
        births=births_array.copy(),
        newborn_deaths=newborn_deaths,
        existing_person_deaths=existing_person_deaths,
        net_migration=migration.copy(),
        existing_head_deaths=existing_head_deaths,
        newborn_heads=newborn_heads,
        new_heads_from_nonheads=new_heads,
        head_dissolutions=dissolutions,
        net_migrant_heads=migrant_heads,
        person_identity_residual=person_residual,
        head_identity_residual=head_residual,
        headship_mapping_residual=mapping_residual,
    )


def infer_net_migration_path(
    target_persons: Mapping[int, Array],
    births: Mapping[int, Array],
    survival: Mapping[int, Array],
    *,
    topcoded_last_age: bool = True,
) -> dict[int, Array]:
    """Back out the net-migration residual that reproduces a target path."""

    years = sorted(int(year) for year in target_persons)
    if len(years) < 2 or any(b != a + 1 for a, b in zip(years, years[1:])):
        raise ValueError("target_persons must contain consecutive annual years")
    first = _finite_array(target_persons[years[0]], "target_persons")
    if np.any(first < 0.0):
        raise ValueError("target_persons contains negative mass")
    result: dict[int, Array] = {}
    for previous_year, year in zip(years, years[1:]):
        if year not in births or year not in survival:
            raise KeyError(f"missing births or survival for {year}")
        previous = _finite_array(
            target_persons[previous_year], f"target_persons[{previous_year}]"
        )
        target = _finite_array(target_persons[year], f"target_persons[{year}]")
        if previous.shape != first.shape or target.shape != first.shape:
            raise ValueError("target population arrays have inconsistent shapes")
        survivors, _, _ = cohort_survivors(
            previous,
            births[year],
            survival[year],
            topcoded_last_age=topcoded_last_age,
        )
        result[year] = target - survivors
    return result


def simulate_path(
    initial_state: CohortState,
    *,
    births: Mapping[int, Array],
    survival: Mapping[int, Array],
    net_migration: Mapping[int, Array],
    headship_rates: Mapping[int, Array] | Array,
    final_year: int,
    topcoded_last_age: bool = True,
    tolerance: float = 1e-8,
) -> tuple[dict[int, CohortState], dict[int, FlowLedger]]:
    """Simulate a consecutive annual demographic path."""

    initial = initial_state.validated(tolerance=tolerance)
    if int(final_year) < initial.year:
        raise ValueError("final_year precedes the initial state")
    states = {initial.year: initial}
    ledgers: dict[int, FlowLedger] = {}
    current = initial
    for year in range(initial.year + 1, int(final_year) + 1):
        if year not in births or year not in survival or year not in net_migration:
            raise KeyError(f"missing demographic input for {year}")
        rates = (
            headship_rates[year]
            if isinstance(headship_rates, Mapping)
            else headship_rates
        )
        current, ledger = advance_one_year(
            current,
            births=births[year],
            survival=survival[year],
            net_migration=net_migration[year],
            headship_rates=rates,
            topcoded_last_age=topcoded_last_age,
            tolerance=tolerance,
        )
        if current.year != year:
            raise RuntimeError("calendar-time state update is inconsistent")
        states[year] = current
        ledgers[year] = ledger
    return states, ledgers
