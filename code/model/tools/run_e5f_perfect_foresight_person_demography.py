#!/usr/bin/env python3
"""Isolated PF evaluator with a coherent person-to-household-head law.

This driver leaves the existing birth-vintage-queue solver untouched. It uses
the same household Bellman and within-age economic-state transition, but births
first create annual person cohorts. Survival, fixed net migration, and explicit
headship rates then determine the household-head mass supplied to every
four-year model age cell.

The command-line entry point is an accounting smoke, not an equilibrium solve.
It verifies the coupled forward operator at the reconstructed 2023 state before
that operator is used in a jointly solved transition and terminal fixed point.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from demographic_transition import (  # noqa: E402
    CohortState,
    HouseholdPersonBlockLedger,
    advance_household_person_block,
    aggregate_heads_to_model_age_cells,
    infer_net_migration_path,
    rebalance_household_distribution_by_age,
    solve_demographic_steady_state,
)
import build_e5f_coherent_person_cohort_path as demographic_source  # noqa: E402
import run_dynamic_population_transition as calendar  # noqa: E402
import run_e5f_open_population_transition as transition  # noqa: E402
import run_e5f_perfect_foresight_transition as pf  # noqa: E402
import run_e5f_post2023_no_policy_continuations as continuation  # noqa: E402


DEFAULT_SOURCE_DIR = demographic_source.DEFAULT_SOURCE_DIR
DEFAULT_HEADSHIP_DIR = demographic_source.DEFAULT_HEADSHIP_DIR
DEFAULT_OUTPUT = (
    ROOT / "output/model/e5f_pf_person_demography_accounting_smoke_20260826a"
)


@dataclass(frozen=True)
class AnnualDemographicPrimitives:
    """Annual person-law inputs in model household-mass units."""

    start_year: int
    last_empirical_year: int
    scale_model_units_per_person: float
    initial_person_state: CohortState
    birth_sex_shares: dict[int, np.ndarray]
    survival: dict[int, np.ndarray]
    net_migration: dict[int, np.ndarray]
    headship_rates: np.ndarray
    raw_acs_headship_rates: np.ndarray
    model_age_headship_alignment_factors: np.ndarray
    initial_household_age_mass: np.ndarray
    initial_person_head_age_mass: np.ndarray
    source_paths: dict[str, Path]

    def _value_for_year(
        self, values: Mapping[int, np.ndarray], year: int
    ) -> np.ndarray:
        requested = int(year)
        if requested in values:
            return np.asarray(values[requested], dtype=float)
        if requested > self.last_empirical_year:
            return np.asarray(values[self.last_empirical_year], dtype=float)
        raise KeyError(f"No demographic primitive is available for {requested}")

    def block_inputs(
        self, start_year: int, period_years: int
    ) -> tuple[
        dict[int, np.ndarray],
        dict[int, np.ndarray],
        dict[int, np.ndarray],
    ]:
        years = range(int(start_year) + 1, int(start_year) + int(period_years) + 1)
        return (
            {year: self._value_for_year(self.birth_sex_shares, year) for year in years},
            {year: self._value_for_year(self.survival, year) for year in years},
            {year: self._value_for_year(self.net_migration, year) for year in years},
        )


@dataclass
class PersonPFState:
    """Joint economic-household and person-cohort state at one model date."""

    g_pre: np.ndarray
    persons: CohortState


@dataclass
class PersonPathEvaluation:
    prices: np.ndarray
    rents: np.ndarray
    values: list[np.ndarray]
    rows: list[dict[str, Any]]
    terminal_state: PersonPFState
    maximum_market_residual: float
    maximum_policy_reproduction_error: float
    maximum_person_identity_error: float
    maximum_head_identity_error: float
    maximum_household_person_head_gap: float
    maximum_age_head_gap: float
    maximum_feasibility_projection_mass: float
    bellman_solves: int
    elapsed_seconds: float


@dataclass
class TerminalHouseholdPersonFixedPoint:
    """Stationary household composition and demographic scale at fixed policy."""

    converged: bool
    iterations: int
    g_pre: np.ndarray
    persons: CohortState
    annual_births_per_head: float
    model_period_births: float
    housing_demand: float
    housing_supply: float
    housing_excess_demand_relative: float
    housing_market_relative_residual: float
    government_budget_residual: float
    scaled_government_budget_residual: float
    implied_equal_transfer: float
    equal_transfer_gap: float
    distribution_mapping_relative_l1: float
    annual_birth_rate_relative_gap: float
    person_one_step_relative_l1: float
    head_one_step_relative_l1: float
    age_head_one_step_max_abs: float
    household_person_head_gap: float
    renewal_ratio: float
    response_iterations: int
    history: list[dict[str, float]]


def _source_paths(source_dir: Path, headship_dir: Path) -> dict[str, Path]:
    return {
        "population_mid": source_dir / "population_mid.csv",
        "births_mid": source_dir / "births_mid.csv",
        "survival": source_dir / "survival.csv",
        "vintage_2025": source_dir / "vintage_2025_age_sex.csv",
        "acs_headship": headship_dir / "acs_headship_profiles.csv",
    }


def build_annual_demographic_primitives(
    household_g_pre: np.ndarray,
    parameters: SimpleNamespace,
    *,
    source_dir: Path = DEFAULT_SOURCE_DIR,
    headship_dir: Path = DEFAULT_HEADSHIP_DIR,
    start_year: int = pf.CALENDAR_START_YEAR,
    last_empirical_year: int = 2100,
) -> AnnualDemographicPrimitives:
    """Align the 2023 person state to the model's household-head age masses.

    A single global scale converts Census persons into the model's household
    mass units. Within each four-year age cell, one multiplicative alignment
    factor preserves the ACS sex-age pattern while making the implied 2023 head
    stock exactly equal the reconstructed model head stock. The factors and raw
    rates are retained for inspection.
    """

    source_root = source_dir.resolve()
    headship_root = headship_dir.resolve()
    paths = _source_paths(source_root, headship_root)
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    population = demographic_source.population_path(
        demographic_source.read_csv(paths["population_mid"])
    )
    vintage = demographic_source.vintage_2025(
        demographic_source.read_csv(paths["vintage_2025"])
    )
    target = demographic_source.cohort_anchor_projection(
        population,
        vintage,
        start_year=int(start_year),
        maximum_year=int(last_empirical_year),
    )
    census_births = demographic_source.births_path(
        demographic_source.read_csv(paths["births_mid"])
    )
    survival = demographic_source.survival_path(
        demographic_source.read_csv(paths["survival"]), nativity=0
    )
    raw_headship = demographic_source.fixed_2023_headship(
        demographic_source.read_csv(paths["acs_headship"])
    )
    migration_unscaled = infer_net_migration_path(
        target,
        census_births,
        survival,
        topcoded_last_age=True,
    )

    distribution = np.asarray(household_g_pre, dtype=float)
    if distribution.ndim != 7:
        raise ValueError(
            f"Expected the seven-dimensional household distribution, got {distribution.shape}"
        )
    if distribution.shape[3] != int(parameters.J):
        raise ValueError("Household age dimension does not match parameters.J")
    household_age_mass = np.sum(distribution, axis=(0, 1, 2, 4, 5, 6))
    model_age_start = int(round(float(parameters.age_start)))
    model_age_width = int(round(float(parameters.da)))
    model_age_stop = model_age_start + model_age_width * int(parameters.J)
    empirical_heads = target[int(start_year)] * raw_headship
    empirical_model_heads = float(
        np.sum(empirical_heads[:, model_age_start:model_age_stop])
    )
    household_total = float(np.sum(household_age_mass))
    if empirical_model_heads <= 0.0 or household_total <= 0.0:
        raise RuntimeError("Initial household or empirical head mass is nonpositive")
    scale = household_total / empirical_model_heads

    aligned_headship = raw_headship.copy()
    alignment = np.ones(int(parameters.J), dtype=float)
    for age_index in range(int(parameters.J)):
        lower = model_age_start + age_index * model_age_width
        upper = lower + model_age_width
        base_mass = scale * float(
            np.sum(target[int(start_year)][:, lower:upper] * raw_headship[:, lower:upper])
        )
        desired_mass = float(household_age_mass[age_index])
        if base_mass <= 0.0:
            raise RuntimeError(f"Empirical head mass is zero in model age cell {age_index}")
        alignment[age_index] = desired_mass / base_mass
        aligned_headship[:, lower:upper] *= alignment[age_index]
    if np.min(aligned_headship) < 0.0 or np.max(aligned_headship) > 1.0:
        raise RuntimeError(
            "Model-age alignment pushed a headship rate outside [0,1]: "
            f"min={np.min(aligned_headship):.6g}, max={np.max(aligned_headship):.6g}"
        )

    initial_persons = scale * target[int(start_year)]
    initial_state = CohortState(
        year=int(start_year),
        persons=initial_persons,
        heads=initial_persons * aligned_headship,
    ).validated()
    person_head_age_mass = aggregate_heads_to_model_age_cells(
        initial_state,
        age_start=model_age_start,
        cell_width=model_age_width,
        number_of_cells=int(parameters.J),
    )
    maximum_gap = float(np.max(np.abs(person_head_age_mass - household_age_mass)))
    if maximum_gap > 2e-10:
        raise RuntimeError(
            f"Initial person/head alignment failed by {maximum_gap:.6g} model units"
        )

    birth_shares: dict[int, np.ndarray] = {}
    for year, births in census_births.items():
        total = float(np.sum(births))
        if total <= 0.0:
            raise RuntimeError(f"Census births are nonpositive in {year}")
        birth_shares[int(year)] = np.asarray(births, dtype=float) / total
    scaled_migration = {
        int(year): scale * np.asarray(value, dtype=float)
        for year, value in migration_unscaled.items()
    }
    return AnnualDemographicPrimitives(
        start_year=int(start_year),
        last_empirical_year=int(last_empirical_year),
        scale_model_units_per_person=scale,
        initial_person_state=initial_state,
        birth_sex_shares=birth_shares,
        survival={int(year): np.asarray(value, dtype=float) for year, value in survival.items()},
        net_migration=scaled_migration,
        headship_rates=aligned_headship,
        raw_acs_headship_rates=raw_headship,
        model_age_headship_alignment_factors=alignment,
        initial_household_age_mass=np.asarray(household_age_mass, dtype=float),
        initial_person_head_age_mass=person_head_age_mass,
        source_paths=paths,
    )


def _owner_rate(distribution: np.ndarray) -> float:
    total = float(np.sum(distribution))
    return float(np.sum(distribution[:, 1:, ...])) / max(total, 1e-15)


def _age_masses(distribution: np.ndarray) -> np.ndarray:
    return np.sum(np.asarray(distribution, dtype=float), axis=(0, 1, 2, 4, 5, 6))


def solve_terminal_household_person_fixed_point(
    *,
    policy: calendar.PolicyBundle,
    parameters: SimpleNamespace,
    b_grid: np.ndarray,
    initial_g_pre: np.ndarray,
    demographic_primitives: AnnualDemographicPrimitives,
    supply_rule: calendar.HousingSupplyRule,
    distribution_tolerance: float = 1e-9,
    birth_rate_tolerance: float = 1e-10,
    one_step_tolerance: float = 1e-8,
    maximum_iterations: int = 250,
    damping: float = 0.50,
) -> TerminalHouseholdPersonFixedPoint:
    """Solve the demographic scale and household composition at fixed policy.

    For each candidate household distribution, the policy implies annual births
    per head. The linear annual cohort law then solves the associated finite
    demographic steady state analytically under fixed terminal survival,
    migration, and headship. The household economic-state distribution is
    advanced once and reweighted to those head-age stocks. Iteration continues
    until this mapping and the birth rate both stop moving; an independent
    four-year one-step check then verifies the resulting fixed point.
    """

    if int(maximum_iterations) < 1:
        raise ValueError("maximum_iterations must be positive")
    if not 0.0 < float(damping) <= 1.0:
        raise ValueError("damping must lie in (0,1]")
    for label, value in (
        ("distribution_tolerance", distribution_tolerance),
        ("birth_rate_tolerance", birth_rate_tolerance),
        ("one_step_tolerance", one_step_tolerance),
    ):
        if not math.isfinite(float(value)) or float(value) <= 0.0:
            raise ValueError(f"{label} must be finite and positive")

    P = copy.deepcopy(parameters)
    shared = calendar.model.precompute_shared(P, b_grid)
    price = float(np.asarray(policy.price, dtype=float).reshape(-1)[0])
    period_years = int(round(float(P.period_years)))
    age_start = int(round(float(P.age_start)))
    age_width = int(round(float(P.da)))
    age_stop = age_start + age_width * int(P.J)
    headship = np.asarray(demographic_primitives.headship_rates, dtype=float)
    outside_model_ages = headship.copy()
    outside_model_ages[:, age_start:age_stop] = 0.0
    if float(np.max(np.abs(outside_model_ages))) > 1e-14:
        raise ValueError(
            "Terminal headship must be zero outside the household model's age "
            "support so births per model head and births per demographic head "
            "are the same object"
        )
    terminal_year = int(demographic_primitives.last_empirical_year)
    birth_shares = demographic_primitives._value_for_year(
        demographic_primitives.birth_sex_shares, terminal_year
    )
    survival = demographic_primitives._value_for_year(
        demographic_primitives.survival, terminal_year
    )
    migration = demographic_primitives._value_for_year(
        demographic_primitives.net_migration, terminal_year
    )
    entry_shares = np.asarray(P.entry_shares, dtype=float).reshape(-1)
    entry_shares /= float(np.sum(entry_shares))
    entrant_template = calendar.entrant_cohort(entry_shares, P, b_grid)

    g_pre = np.asarray(initial_g_pre, dtype=float).copy()
    if g_pre.ndim != 7 or g_pre.shape[3] != int(P.J):
        raise ValueError("initial_g_pre has the wrong household-state shape")
    if np.any(~np.isfinite(g_pre)) or np.min(g_pre) < -1e-13:
        raise ValueError("initial_g_pre must be finite and nonnegative")
    g_pre = np.maximum(g_pre, 0.0)
    previous_birth_rate: float | None = None
    history: list[dict[str, float]] = []
    final_ss = None
    final_evaluation = None
    final_births = math.nan
    iteration_birth_rate = math.nan
    mapping_gap = math.inf
    birth_rate_gap = math.inf

    for iteration in range(1, int(maximum_iterations) + 1):
        evaluation = calendar.evaluate_period(
            np.array([price]),
            g_pre,
            P,
            b_grid,
            shared,
            calendar.SolveCounter(),
            supply_rule=supply_rule,
            supplied_policy=policy,
        )
        accounting = transition.calendar_topcode_birth_accounting(
            evaluation.g_pre,
            evaluation.g_post_fertility,
            float(evaluation.births),
            P,
        )
        adjusted_births = float(accounting["topcode_adjusted_birth_children"])
        head_mass = float(np.sum(evaluation.g_current))
        annual_birth_rate = adjusted_births / max(
            head_mass * period_years, 1e-15
        )
        iteration_birth_rate = annual_birth_rate
        demographic_ss = solve_demographic_steady_state(
            annual_births_per_head=annual_birth_rate,
            birth_sex_shares=birth_shares,
            survival=survival,
            net_migration=migration,
            headship_rates=headship,
            topcoded_last_age=True,
            tolerance=1e-11,
            max_iterations=20000,
        )
        demographic_state = CohortState(
            year=terminal_year,
            persons=demographic_ss.persons,
            heads=demographic_ss.heads,
        )
        target_age_mass = aggregate_heads_to_model_age_cells(
            demographic_state,
            age_start=age_start,
            cell_width=age_width,
            number_of_cells=int(P.J),
        )
        raw_next, _, _, raw_mass_residual = (
            transition.advance_sequential_calendar_distribution(
                evaluation,
                np.zeros(int(P.I)),
                P,
                b_grid,
                shared,
            )
        )
        if abs(float(raw_mass_residual)) > 2e-8:
            raise RuntimeError(
                f"Raw household transition mass residual is {raw_mass_residual:.6g}"
            )
        proposed, _ = rebalance_household_distribution_by_age(
            raw_next,
            target_age_mass,
            age_axis=3,
            empty_age_templates={0: entrant_template},
            tolerance=1e-10,
        )
        scale = max(float(np.sum(proposed)), 1e-15)
        mapping_gap = float(np.sum(np.abs(proposed - g_pre))) / scale
        birth_rate_gap = (
            math.inf
            if previous_birth_rate is None
            else abs(annual_birth_rate - previous_birth_rate)
            / max(abs(annual_birth_rate), 1e-15)
        )
        history.append(
            {
                "iteration": float(iteration),
                "annual_births_per_head": annual_birth_rate,
                "household_heads": head_mass,
                "resident_persons": float(np.sum(demographic_ss.persons)),
                "renewal_ratio": float(demographic_ss.renewal_ratio),
                "distribution_mapping_relative_l1": mapping_gap,
                "annual_birth_rate_relative_gap": birth_rate_gap,
                "housing_market_relative_residual": float(
                    evaluation.relative_market_residual
                ),
            }
        )
        final_ss = demographic_ss
        final_evaluation = evaluation
        final_births = adjusted_births
        if (
            mapping_gap <= float(distribution_tolerance)
            and birth_rate_gap <= float(birth_rate_tolerance)
        ):
            g_pre = proposed
            break
        if iteration == 1 or float(damping) == 1.0:
            g_pre = proposed
        else:
            mixed = (1.0 - float(damping)) * g_pre + float(damping) * proposed
            g_pre, _ = rebalance_household_distribution_by_age(
                mixed,
                target_age_mass,
                age_axis=3,
                empty_age_templates={0: entrant_template},
                tolerance=1e-10,
            )
        previous_birth_rate = annual_birth_rate
    else:
        iteration = int(maximum_iterations)

    if final_ss is None or final_evaluation is None:
        raise RuntimeError("Terminal fixed-point iteration produced no evaluation")

    final_evaluation = calendar.evaluate_period(
        np.array([price]),
        g_pre,
        P,
        b_grid,
        shared,
        calendar.SolveCounter(),
        supply_rule=supply_rule,
        supplied_policy=policy,
    )
    final_accounting = transition.calendar_topcode_birth_accounting(
        final_evaluation.g_pre,
        final_evaluation.g_post_fertility,
        float(final_evaluation.births),
        P,
    )
    final_births = float(final_accounting["topcode_adjusted_birth_children"])
    final_head_mass = float(np.sum(final_evaluation.g_current))
    final_rate = final_births / max(final_head_mass * period_years, 1e-15)
    final_ss = solve_demographic_steady_state(
        annual_births_per_head=final_rate,
        birth_sex_shares=birth_shares,
        survival=survival,
        net_migration=migration,
        headship_rates=headship,
        topcoded_last_age=True,
        tolerance=1e-11,
        max_iterations=20000,
    )
    raw_next, _, _, raw_mass_residual = transition.advance_sequential_calendar_distribution(
        final_evaluation,
        np.zeros(int(P.I)),
        P,
        b_grid,
        shared,
    )
    if abs(float(raw_mass_residual)) > 2e-8:
        raise RuntimeError(
            f"Final raw household transition mass residual is {raw_mass_residual:.6g}"
        )
    shares = {
        year: birth_shares
        for year in range(terminal_year + 1, terminal_year + period_years + 1)
    }
    survival_block = {
        year: survival
        for year in range(terminal_year + 1, terminal_year + period_years + 1)
    }
    migration_block = {
        year: migration
        for year in range(terminal_year + 1, terminal_year + period_years + 1)
    }
    stationary_state = CohortState(
        year=terminal_year,
        persons=final_ss.persons,
        heads=final_ss.heads,
    )
    one_step_g, one_step_persons, one_step_ledger = advance_household_person_block(
        raw_next,
        stationary_state,
        model_period_births=final_births,
        period_years=period_years,
        birth_sex_shares=shares,
        survival=survival_block,
        net_migration=migration_block,
        headship_rates=headship,
        model_age_start=age_start,
        model_age_cell_width=age_width,
        number_of_model_age_cells=int(P.J),
        age_axis=3,
        empty_age_templates={0: entrant_template},
        tolerance=1e-10,
    )
    household_scale = max(float(np.sum(g_pre)), 1e-15)
    mapping_gap = float(np.sum(np.abs(one_step_g - g_pre))) / household_scale
    person_scale = max(float(np.sum(final_ss.persons)), 1e-15)
    head_scale = max(float(np.sum(final_ss.heads)), 1e-15)
    person_one_step_gap = float(
        np.sum(np.abs(one_step_persons.persons - final_ss.persons))
    ) / person_scale
    head_one_step_gap = float(
        np.sum(np.abs(one_step_persons.heads - final_ss.heads))
    ) / head_scale
    final_target_age = aggregate_heads_to_model_age_cells(
        stationary_state,
        age_start=age_start,
        cell_width=age_width,
        number_of_cells=int(P.J),
    )
    age_gap = float(np.max(np.abs(_age_masses(one_step_g) - final_target_age)))
    birth_rate_gap = (
        abs(final_rate - iteration_birth_rate) / max(abs(final_rate), 1e-15)
    )
    property_tax_revenue = float(
        calendar.model.property_tax_revenue_from_distribution(
            final_evaluation.g_current,
            final_evaluation.policy.hR_pol,
            final_evaluation.policy.price,
            P,
        )
    )
    transfer = float(getattr(P, "property_tax_lump_sum_transfer", 0.0))
    fiscal_residual = property_tax_revenue - transfer * final_head_mass
    housing_demand = float(final_evaluation.demand_by_loc[0])
    housing_supply = float(final_evaluation.supply_by_loc[0])
    housing_excess = (housing_demand - housing_supply) / max(
        abs(housing_supply), 1e-15
    )
    implied_transfer = property_tax_revenue / max(final_head_mass, 1e-15)
    scaled_fiscal_residual = fiscal_residual / max(
        abs(property_tax_revenue),
        abs(transfer * final_head_mass),
        1e-15,
    )
    converged = bool(
        mapping_gap <= float(distribution_tolerance)
        and birth_rate_gap <= max(float(birth_rate_tolerance), 1e-8)
        and person_one_step_gap <= float(one_step_tolerance)
        and head_one_step_gap <= float(one_step_tolerance)
        and age_gap <= float(one_step_tolerance) * max(1.0, household_scale)
        and abs(one_step_ledger.household_person_head_gap)
        <= float(one_step_tolerance) * max(1.0, household_scale)
    )
    return TerminalHouseholdPersonFixedPoint(
        converged=converged,
        iterations=int(iteration),
        g_pre=g_pre,
        persons=stationary_state,
        annual_births_per_head=final_rate,
        model_period_births=final_births,
        housing_demand=housing_demand,
        housing_supply=housing_supply,
        housing_excess_demand_relative=housing_excess,
        housing_market_relative_residual=float(
            final_evaluation.relative_market_residual
        ),
        government_budget_residual=fiscal_residual,
        scaled_government_budget_residual=scaled_fiscal_residual,
        implied_equal_transfer=implied_transfer,
        equal_transfer_gap=implied_transfer - transfer,
        distribution_mapping_relative_l1=mapping_gap,
        annual_birth_rate_relative_gap=birth_rate_gap,
        person_one_step_relative_l1=person_one_step_gap,
        head_one_step_relative_l1=head_one_step_gap,
        age_head_one_step_max_abs=age_gap,
        household_person_head_gap=one_step_ledger.household_person_head_gap,
        renewal_ratio=float(final_ss.renewal_ratio),
        response_iterations=int(final_ss.response_iterations),
        history=history,
    )


def evaluate_path_at_prices_person_demography(
    *,
    prices: Sequence[float],
    psi_path: Sequence[float],
    terminal_price: float,
    terminal_V: np.ndarray,
    base_parameters: SimpleNamespace,
    b_grid: np.ndarray,
    initial_state: PersonPFState,
    demographic_primitives: AnnualDemographicPrimitives,
    supply_rule: calendar.HousingSupplyRule,
    transfer_path: Sequence[float] | None = None,
) -> PersonPathEvaluation:
    """Evaluate a price path with endogenous births and coherent head stocks."""

    started = time.perf_counter()
    price_path = np.asarray(prices, dtype=float).reshape(-1)
    psi_values = np.asarray(psi_path, dtype=float).reshape(-1)
    if price_path.shape != psi_values.shape or len(price_path) < 1:
        raise ValueError("Price and preference paths must have the same positive length")
    transfers = (
        np.zeros_like(price_path)
        if transfer_path is None
        else np.asarray(transfer_path, dtype=float).reshape(-1)
    )
    if transfers.shape != price_path.shape:
        raise ValueError("The equal-transfer path must match the price path")
    if np.any(~np.isfinite(transfers)) or np.any(transfers < 0.0):
        raise ValueError("Equal transfers must be finite and nonnegative")

    rents = pf.rents_from_asset_prices(price_path, terminal_price, base_parameters)
    values, backward_solves = pf.backward_value_path(
        prices=price_path,
        rents=rents,
        psi_path=psi_values,
        terminal_V=terminal_V,
        base_parameters=base_parameters,
        b_grid=b_grid,
        transfer_path=transfers,
    )
    state = PersonPFState(
        g_pre=np.asarray(initial_state.g_pre, dtype=float).copy(),
        persons=CohortState(
            year=int(initial_state.persons.year),
            persons=np.asarray(initial_state.persons.persons, dtype=float).copy(),
            heads=np.asarray(initial_state.persons.heads, dtype=float).copy(),
        ).validated(),
    )
    rows: list[dict[str, Any]] = []
    maximum_reproduction = 0.0
    maximum_person_identity = 0.0
    maximum_head_identity = 0.0
    maximum_household_head_gap = 0.0
    maximum_age_gap = 0.0
    maximum_projection = 0.0
    forward_solves = 0

    for period, (price, rent, psi, transfer_value) in enumerate(
        zip(price_path, rents, psi_values, transfers, strict=True)
    ):
        parameters = copy.deepcopy(base_parameters)
        parameters.psi_child = float(psi)
        parameters.property_tax_lump_sum_transfer = float(transfer_value)
        period_years = int(round(float(parameters.period_years)))
        expected_year = pf.CALENDAR_START_YEAR + period * period_years
        if state.persons.year != expected_year:
            raise RuntimeError(
                f"Person state year {state.persons.year} does not match model date {expected_year}"
            )

        current_target_heads = aggregate_heads_to_model_age_cells(
            state.persons,
            age_start=int(round(float(parameters.age_start))),
            cell_width=int(round(float(parameters.da))),
            number_of_cells=int(parameters.J),
        )
        current_age_gap = float(
            np.max(np.abs(_age_masses(state.g_pre) - current_target_heads))
        )
        maximum_age_gap = max(maximum_age_gap, current_age_gap)
        if current_age_gap > 2e-9:
            raise RuntimeError(
                f"Household/person head-age identity fails at {expected_year}: {current_age_gap:.6g}"
            )

        shared = calendar.model.precompute_shared(parameters, b_grid)
        policy = pf.solve_date_policy(
            price=float(price),
            rent=float(rent),
            P=parameters,
            b_grid=b_grid,
            shared=shared,
            continuation_V=values[period + 1],
        )
        forward_solves += 1
        reproduction = float(np.max(np.abs(policy.V - values[period])))
        maximum_reproduction = max(maximum_reproduction, reproduction)
        evaluation = calendar.evaluate_period(
            np.array([float(price)]),
            state.g_pre,
            parameters,
            b_grid,
            shared,
            calendar.SolveCounter(),
            supply_rule=supply_rule,
            supplied_policy=policy,
        )
        accounting = transition.calendar_topcode_birth_accounting(
            evaluation.g_pre,
            evaluation.g_post_fertility,
            float(evaluation.births),
            parameters,
        )
        adjusted_births = float(accounting["topcode_adjusted_birth_children"])

        raw_next, model_mature_by_loc, raw_exit_mass, raw_mass_residual = (
            transition.advance_sequential_calendar_distribution(
                evaluation,
                np.zeros(int(parameters.I)),
                parameters,
                b_grid,
                shared,
            )
        )
        entry_shares = np.asarray(parameters.entry_shares, dtype=float).reshape(-1)
        entry_shares /= float(np.sum(entry_shares))
        entrant_template = calendar.entrant_cohort(
            entry_shares, parameters, b_grid
        )
        shares, survival, migration = demographic_primitives.block_inputs(
            state.persons.year, period_years
        )
        next_g, next_persons, coupled = advance_household_person_block(
            raw_next,
            state.persons,
            model_period_births=adjusted_births,
            period_years=period_years,
            birth_sex_shares=shares,
            survival=survival,
            net_migration=migration,
            headship_rates=demographic_primitives.headship_rates,
            model_age_start=int(round(float(parameters.age_start))),
            model_age_cell_width=int(round(float(parameters.da))),
            number_of_model_age_cells=int(parameters.J),
            age_axis=3,
            empty_age_templates={0: entrant_template},
        )
        maximum_person_identity = max(
            maximum_person_identity, coupled.person.person_identity_max_abs
        )
        maximum_head_identity = max(
            maximum_head_identity, coupled.person.head_identity_max_abs
        )
        maximum_household_head_gap = max(
            maximum_household_head_gap, abs(coupled.household_person_head_gap)
        )
        maximum_projection = max(
            maximum_projection, float(evaluation.feasibility_projection_mass)
        )
        health = calendar.distribution_health(
            {
                "pre": evaluation.g_pre,
                "post_fertility": evaluation.g_post_fertility,
                "current": evaluation.g_current,
                "raw_next_pre": raw_next,
                "coupled_next_pre": next_g,
            }
        )
        if (
            int(health["nonfinite_distribution_count"]) != 0
            or health["min_distribution_mass"] is None
            or float(health["min_distribution_mass"]) < -1e-13
        ):
            raise RuntimeError(f"Distribution-health gate failed: {health}")

        current_heads = float(np.sum(evaluation.g_current))
        current_persons = float(np.sum(state.persons.persons))
        property_tax_revenue = float(
            calendar.model.property_tax_revenue_from_distribution(
                evaluation.g_current,
                evaluation.policy.hR_pol,
                evaluation.policy.price,
                parameters,
            )
        )
        transfer_outlays = float(transfer_value) * current_heads
        rows.append(
            {
                "period": period,
                "calendar_year": expected_year,
                "psi_child": float(psi),
                "asset_price": float(price),
                "renter_price": float(rent),
                "housing_demand": float(evaluation.demand_by_loc[0]),
                "housing_supply": float(evaluation.supply_by_loc[0]),
                "relative_market_residual": float(evaluation.relative_market_residual),
                "resident_persons": current_persons,
                "household_heads": current_heads,
                "persons_per_household_head": current_persons / max(current_heads, 1e-15),
                "birth_children": float(evaluation.births),
                "birth_children_topcode_adjusted": adjusted_births,
                "births_per_household_head": adjusted_births / max(current_heads, 1e-15),
                "next_resident_persons": float(np.sum(next_persons.persons)),
                "next_household_heads": float(np.sum(next_g)),
                "annual_net_migration_over_period": coupled.person.total_net_migration,
                "new_heads_from_nonheads_over_period": coupled.person.total_new_heads_from_nonheads,
                "head_dissolutions_over_period": coupled.person.total_head_dissolutions,
                "net_migrant_heads_over_period": coupled.person.total_net_migrant_heads,
                "household_bridge_added_mass": float(
                    np.sum(coupled.household_heads.added_mass)
                ),
                "household_bridge_removed_mass": float(
                    np.sum(coupled.household_heads.removed_mass)
                ),
                "household_person_head_gap": coupled.household_person_head_gap,
                "person_identity_max_abs": coupled.person.person_identity_max_abs,
                "head_identity_max_abs": coupled.person.head_identity_max_abs,
                "raw_household_exit_mass": float(raw_exit_mass),
                "raw_household_mass_residual": float(raw_mass_residual),
                "model_state_same_period_mature_flow_B": float(
                    parameters.entrant_conversion_factor
                )
                * float(np.sum(model_mature_by_loc)),
                "property_tax_revenue": property_tax_revenue,
                "equal_transfer_period_units": float(transfer_value),
                "equal_transfer_outlays": transfer_outlays,
                "government_budget_residual": property_tax_revenue - transfer_outlays,
                "owner_rate": _owner_rate(evaluation.g_current),
                "policy_reproduction_max_abs": reproduction,
                "feasibility_frontier_projection_mass": float(
                    evaluation.feasibility_projection_mass
                ),
                "minimum_distribution_mass": float(health["min_distribution_mass"]),
                "nonfinite_distribution_count": int(
                    health["nonfinite_distribution_count"]
                ),
            }
        )
        state = PersonPFState(g_pre=next_g, persons=next_persons)

    return PersonPathEvaluation(
        prices=price_path,
        rents=rents,
        values=values,
        rows=rows,
        terminal_state=state,
        maximum_market_residual=max(
            abs(float(row["relative_market_residual"])) for row in rows
        ),
        maximum_policy_reproduction_error=maximum_reproduction,
        maximum_person_identity_error=maximum_person_identity,
        maximum_head_identity_error=maximum_head_identity,
        maximum_household_person_head_gap=maximum_household_head_gap,
        maximum_age_head_gap=maximum_age_gap,
        maximum_feasibility_projection_mass=maximum_projection,
        bellman_solves=backward_solves + forward_solves,
        elapsed_seconds=time.perf_counter() - started,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--headship-dir", type=Path, default=DEFAULT_HEADSHIP_DIR)
    parser.add_argument("--selected-report", type=Path, default=pf.DEFAULT_REPORT)
    parser.add_argument(
        "--selected-transition", type=Path, default=pf.DEFAULT_SELECTED_TRANSITION
    )
    parser.add_argument("--source", type=Path, default=pf.DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--periods", type=int, default=1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if int(args.periods) < 1:
        raise ValueError("--periods must be positive")
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    chain, model = transition.configure_sequential_model()
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    calendar.distribution_rows = transition.independent_child_distribution_rows
    contracts, _ = pf.load_diagnostic_contracts(
        args.selected_report.resolve(),
        args.selected_transition.resolve(),
        args.source.resolve(),
    )
    prepared = continuation.prepare_model(
        contracts, args.source.resolve(), chain, model
    )
    initial_old, _, history_gate, impact_price = pf.reconstruct_2023_state(
        prepared,
        contracts,
        args.selected_transition.resolve(),
        market_tolerance=2e-4,
    )
    primitives = build_annual_demographic_primitives(
        initial_old.g_pre,
        prepared.parameters,
        source_dir=args.source_dir,
        headship_dir=args.headship_dir,
    )
    best = contracts["report_best"]
    current_psi = float(best["new_psi_child"])
    parameters = copy.deepcopy(prepared.parameters)
    parameters.psi_child = current_psi
    periods = int(args.periods)
    evaluation = evaluate_path_at_prices_person_demography(
        prices=np.full(periods, float(impact_price)),
        psi_path=np.full(periods, current_psi),
        terminal_price=float(impact_price),
        terminal_V=prepared.old_policy.V,
        base_parameters=parameters,
        b_grid=prepared.b_grid,
        initial_state=PersonPFState(
            g_pre=initial_old.g_pre,
            persons=primitives.initial_person_state,
        ),
        demographic_primitives=primitives,
        supply_rule=prepared.supply_rule,
    )
    pf.write_csv(output_dir / "accounting_smoke_path.csv", evaluation.rows)
    initial_total_people = float(np.sum(primitives.initial_person_state.persons))
    initial_total_heads = float(np.sum(primitives.initial_person_state.heads[:, 18:86]))
    payload = {
        "status": "passed_accounting_smoke_not_equilibrium",
        "interpretation": (
            "The household continuation value and constant price path are temporary "
            "boundary inputs. This packet tests only the coupled forward demographic "
            "operator and must not be read as a solved equilibrium transition."
        ),
        "periods": periods,
        "history_reproduction_status": history_gate.get("status"),
        "initial_resident_persons_model_units": initial_total_people,
        "initial_household_heads_model_units": initial_total_heads,
        "initial_actual_resident_persons": initial_total_people
        / primitives.scale_model_units_per_person,
        "model_units_per_actual_person": primitives.scale_model_units_per_person,
        "model_age_headship_alignment_factors": (
            primitives.model_age_headship_alignment_factors.tolist()
        ),
        "maximum_market_residual_not_gated": evaluation.maximum_market_residual,
        "maximum_policy_reproduction_error": (
            evaluation.maximum_policy_reproduction_error
        ),
        "maximum_person_identity_error": evaluation.maximum_person_identity_error,
        "maximum_head_identity_error": evaluation.maximum_head_identity_error,
        "maximum_household_person_head_gap": (
            evaluation.maximum_household_person_head_gap
        ),
        "maximum_age_head_gap": evaluation.maximum_age_head_gap,
        "maximum_feasibility_projection_mass": (
            evaluation.maximum_feasibility_projection_mass
        ),
        "source_hashes": {
            name: pf.file_sha256(path)
            for name, path in primitives.source_paths.items()
        },
        "source_code_hashes": {
            "person_demography_driver": pf.file_sha256(Path(__file__).resolve()),
            "demographic_source_builder": pf.file_sha256(
                TOOLS_ROOT / "build_e5f_coherent_person_cohort_path.py"
            ),
            "demographic_satellite_helper": pf.file_sha256(
                TOOLS_ROOT / "build_e5f_persons_demographic_satellite.py"
            ),
            "perfect_foresight_driver": pf.file_sha256(
                TOOLS_ROOT / "run_e5f_perfect_foresight_transition.py"
            ),
            "household_solver": pf.file_sha256(
                MODEL_ROOT / "intergen_eqscale_seq_optimized/solver.py"
            ),
            "person_cohort_law": pf.file_sha256(
                MODEL_ROOT / "demographic_transition/person_cohort_law.py"
            ),
            "four_year_bridge": pf.file_sha256(
                MODEL_ROOT / "demographic_transition/four_year_bridge.py"
            ),
            "household_head_bridge": pf.file_sha256(
                MODEL_ROOT / "demographic_transition/household_head_bridge.py"
            ),
            "household_person_coupling": pf.file_sha256(
                MODEL_ROOT / "demographic_transition/household_person_coupling.py"
            ),
        },
    }
    gates = {
        "policy_reproduction": evaluation.maximum_policy_reproduction_error <= 2e-8,
        "person_identity": evaluation.maximum_person_identity_error <= 2e-8,
        "head_identity": evaluation.maximum_head_identity_error <= 2e-8,
        "household_person_head_identity": (
            evaluation.maximum_household_person_head_gap <= 2e-9
        ),
        "age_head_identity": evaluation.maximum_age_head_gap <= 2e-9,
        "distribution_feasibility": (
            evaluation.maximum_feasibility_projection_mass <= 2e-8
        ),
    }
    payload["accounting_gates"] = gates
    if not all(gates.values()):
        payload["status"] = "failed_accounting_smoke"
    (output_dir / "summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if not all(gates.values()):
        raise RuntimeError(f"Person-demography accounting smoke failed: {gates}")


if __name__ == "__main__":
    main()
