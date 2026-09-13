#!/usr/bin/env python3
"""Isolated demographic adapters for the E5F overnight A0/B0 exercises.

The surviving-maturation adapter deliberately does not choose a conversion
from model child units to persons or from persons to household heads.  Those
are measurement/formation primitives, not numerical conventions.  Production
use therefore requires a labelled external contract; diagnostic use still
requires the caller to provide the two factors explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

import numpy as np


Array = np.ndarray
AdvanceDistribution = Callable[[Any, Array, Any, Array, Any], tuple[Array, Array, float, float]]
EntrantCohort = Callable[[Array, Any, Array], Array]

_PRODUCTION_CLASSIFICATIONS = {
    "estimated",
    "empirically_normalized",
    "externally_fixed",
}
_ALL_CLASSIFICATIONS = _PRODUCTION_CLASSIFICATIONS | {"diagnostic"}


@dataclass(frozen=True)
class FormationUnitContract:
    """Explicit map from capped model child units to new household heads.

    ``mature_persons_per_model_child_unit`` owns any mapping required because
    the representative 3+ fertility state is coded as three model child units.
    ``household_heads_per_mature_person`` is the household-formation map.  The
    adapter never infers either number from completed fertility or replacement.
    """

    mature_persons_per_model_child_unit: float
    household_heads_per_mature_person: float
    classification: str
    source: str
    model_child_unit_definition: str
    mature_person_definition: str
    household_head_definition: str

    def validated(self, *, production: bool) -> "FormationUnitContract":
        values = (
            self.mature_persons_per_model_child_unit,
            self.household_heads_per_mature_person,
        )
        if any(not np.isfinite(value) or value < 0.0 for value in values):
            raise ValueError("formation conversion factors must be finite and nonnegative")
        classification = str(self.classification).strip().lower()
        allowed = _PRODUCTION_CLASSIFICATIONS if production else _ALL_CLASSIFICATIONS
        if classification not in allowed:
            if production and classification == "diagnostic":
                raise ValueError(
                    "production surviving-maturation entry requires an estimated, "
                    "empirically_normalized, or externally_fixed formation contract; "
                    "a diagnostic conversion cannot be promoted"
                )
            raise ValueError(f"unsupported formation-contract classification: {classification!r}")
        labels = {
            "source": self.source,
            "model_child_unit_definition": self.model_child_unit_definition,
            "mature_person_definition": self.mature_person_definition,
            "household_head_definition": self.household_head_definition,
        }
        missing = [name for name, value in labels.items() if not str(value).strip()]
        if missing:
            raise ValueError("formation contract is missing: " + ", ".join(missing))
        return self

    @property
    def household_entries_per_model_child_unit(self) -> float:
        return float(self.mature_persons_per_model_child_unit) * float(
            self.household_heads_per_mature_person
        )


@dataclass(frozen=True)
class OutsideEntryContract:
    """Evidence and state-allocation definition for a migration-on B+ flow."""

    classification: str
    source: str
    input_unit_definition: str
    household_state_allocation_definition: str

    def validated(self, *, production: bool) -> "OutsideEntryContract":
        classification = str(self.classification).strip().lower()
        allowed = _PRODUCTION_CLASSIFICATIONS if production else _ALL_CLASSIFICATIONS
        if classification not in allowed:
            raise ValueError(
                "migration-on production requires an estimated, empirically_normalized, "
                "or externally_fixed person-to-household allocation contract"
            )
        labels = {
            "source": self.source,
            "input_unit_definition": self.input_unit_definition,
            "household_state_allocation_definition": self.household_state_allocation_definition,
        }
        missing = [name for name, value in labels.items() if not str(value).strip()]
        if missing:
            raise ValueError("outside-entry contract is missing: " + ", ".join(missing))
        return self


@dataclass(frozen=True)
class SurvivingMaturationLedger:
    """Auditable child and household identities for one four-year transition."""

    post_birth_dependent_units: float
    dependent_deaths: float
    surviving_maturations_model_child_units: float
    next_dependent_units: float
    child_identity_residual: float
    starting_households: float
    household_exits: float
    domestic_formation_entries: float
    outside_migration_entries: float
    next_households: float
    household_identity_residual: float
    base_transition_mass_residual: float
    mature_model_child_units_by_location: Array
    domestic_entries_by_location: Array
    outside_entries_by_location: Array
    total_entries_by_location: Array
    formation_contract_classification: str
    formation_contract_source: str


def dependent_units_by_location(distribution: Array) -> Array:
    """Count children currently at home in each location.

    The expected E5F axes are
    ``wealth, tenure, location, age, income, children-ever-born, children-at-home``.
    In the 3+ state this is the model's capped child-unit count, not a claim that
    the family literally has exactly three children.
    """

    values = np.asarray(distribution, dtype=float)
    if values.ndim != 7:
        raise ValueError(f"household distribution must be seven-dimensional; got {values.shape}")
    if not np.all(np.isfinite(values)) or np.min(values) < -1.0e-13:
        raise ValueError("household distribution must be finite and nonnegative")
    result = np.zeros(values.shape[2], dtype=float)
    for children_ever_born in range(1, values.shape[5]):
        largest_at_home = min(children_ever_born, values.shape[6] - 1)
        for children_at_home in range(1, largest_at_home + 1):
            result += children_at_home * np.sum(
                values[:, :, :, :, :, children_ever_born, children_at_home],
                axis=(0, 1, 3, 4),
            )
    return result


def _dependent_deaths_by_location(distribution: Array, parameters: Any) -> Array:
    values = np.asarray(distribution, dtype=float)
    locations = values.shape[2]
    ages = values.shape[3]
    deaths = np.zeros(locations, dtype=float)
    survival = np.ones(max(ages - 1, 0), dtype=float)
    if bool(getattr(parameters, "use_age_survival", False)):
        survival = np.asarray(parameters.survival_probs, dtype=float).reshape(-1)
        if survival.size != ages - 1:
            raise ValueError(
                f"survival_probs has length {survival.size}; expected {ages - 1}"
            )
        if np.any(~np.isfinite(survival)) or np.any((survival < 0.0) | (survival > 1.0)):
            raise ValueError("survival probabilities must lie in [0, 1]")
    for age in range(ages - 1):
        deaths += (1.0 - survival[age]) * dependent_units_by_location(
            values[:, :, :, age : age + 1, :, :, :]
        )
    # The last model-age household exits with all dependents, as in the frozen
    # sequential calendar operator.
    deaths += dependent_units_by_location(values[:, :, :, -1:, :, :, :])
    return deaths


def _nonnegative_vector(value: Array | None, *, length: int, name: str) -> Array:
    if value is None:
        return np.zeros(length, dtype=float)
    vector = np.asarray(value, dtype=float).reshape(-1)
    if vector.size != length:
        raise ValueError(f"{name} has length {vector.size}; expected {length}")
    if np.any(~np.isfinite(vector)) or np.any(vector < 0.0):
        raise ValueError(f"{name} must be finite and nonnegative")
    return vector.copy()


def _validate_stochastic_maturation(parameters: Any, post_birth: Array) -> None:
    if not bool(getattr(parameters, "use_stochastic_aging", False)):
        raise ValueError("surviving-maturation closure B requires stochastic child maturation")
    if str(getattr(parameters, "child_state_mode", "")).strip().lower() != "independent_count":
        raise ValueError("closure B requires the independent_count children-at-home process")
    transition = np.asarray(getattr(parameters, "Pi_child", None), dtype=float)
    expected_shape = (post_birth.shape[6], post_birth.shape[6], post_birth.shape[5])
    if transition.shape != expected_shape:
        raise ValueError(f"Pi_child has shape {transition.shape}; expected {expected_shape}")
    if np.any(~np.isfinite(transition)) or np.any((transition < 0.0) | (transition > 1.0)):
        raise ValueError("Pi_child probabilities must lie in [0, 1]")
    for children_ever_born in range(post_birth.shape[5]):
        for children_at_home in range(children_ever_born + 1):
            if not np.isclose(
                float(np.sum(transition[children_at_home, :, children_ever_born])),
                1.0,
                rtol=0.0,
                atol=1.0e-12,
            ):
                raise ValueError("Pi_child has a non-stochastic reachable row")


def advance_surviving_maturation_distribution(
    evaluation: Any,
    parameters: Any,
    b_grid: Array,
    shared: Any,
    *,
    advance_distribution: AdvanceDistribution,
    entrant_cohort: EntrantCohort,
    formation_contract: FormationUnitContract,
    outside_entries_by_location: Array | None = None,
    outside_entry_contract: OutsideEntryContract | None = None,
    zero_migration: bool = True,
    production: bool = True,
    tolerance: float = 1.0e-9,
) -> tuple[Array, Array, SurvivingMaturationLedger]:
    """Advance B using only survivor-conditioned maturation plus declared entry.

    ``advance_distribution`` should be the frozen
    ``advance_sequential_calendar_distribution`` and ``entrant_cohort`` should
    be the matching calendar entrant builder.  The first call uses exactly zero
    entry, so no legacy queue or outside-origin residual can enter implicitly.
    """

    contract = formation_contract.validated(production=production)
    post_birth = np.asarray(evaluation.g_post_fertility, dtype=float)
    locations = int(getattr(parameters, "I"))
    if post_birth.ndim != 7 or post_birth.shape[2] != locations:
        raise ValueError("evaluation.g_post_fertility is incompatible with parameters.I")
    _validate_stochastic_maturation(parameters, post_birth)
    outside = _nonnegative_vector(
        outside_entries_by_location,
        length=locations,
        name="outside_entries_by_location",
    )
    if zero_migration and np.any(outside != 0.0):
        raise ValueError("zero-migration B0 requires every outside-entry cell to be exactly zero")
    if zero_migration and outside_entry_contract is not None:
        raise ValueError("zero-migration B0 cannot carry an outside-entry contract")
    if not zero_migration and outside_entries_by_location is None:
        raise ValueError("migration-on B+ requires an explicitly mapped outside-entry vector")
    if not zero_migration:
        if outside_entry_contract is None:
            raise ValueError(
                "migration-on B+ requires an explicit person-to-household state-allocation contract"
            )
        outside_entry_contract.validated(production=production)

    zero_entry = np.zeros(locations, dtype=float)
    incumbent_next, mature_by_loc, household_exits, base_residual = advance_distribution(
        evaluation, zero_entry, parameters, np.asarray(b_grid, dtype=float), shared
    )
    incumbent_next = np.asarray(incumbent_next, dtype=float).copy()
    mature = _nonnegative_vector(
        mature_by_loc,
        length=locations,
        name="surviving maturation flow",
    )
    domestic = mature * contract.household_entries_per_model_child_unit
    total_entry = domestic + outside
    if float(np.sum(np.abs(incumbent_next[:, :, :, 0, :, :, :]))) > tolerance:
        raise RuntimeError("zero-entry incumbent transition unexpectedly populated entrant age")
    entrant_mass = np.asarray(
        entrant_cohort(total_entry, parameters, np.asarray(b_grid, dtype=float)),
        dtype=float,
    )
    if entrant_mass.shape != incumbent_next[:, :, :, 0, :, :, :].shape:
        raise ValueError("entrant cohort shape does not match the youngest-age household slice")
    incumbent_next[:, :, :, 0, :, :, :] = entrant_mass

    post_dependents = float(np.sum(dependent_units_by_location(post_birth)))
    dependent_deaths = float(np.sum(_dependent_deaths_by_location(post_birth, parameters)))
    next_dependents = float(np.sum(dependent_units_by_location(incumbent_next)))
    mature_total = float(np.sum(mature))
    child_residual = next_dependents - (post_dependents - dependent_deaths - mature_total)

    starting_households = float(np.sum(post_birth))
    exits = float(household_exits)
    entries = float(np.sum(total_entry))
    next_households = float(np.sum(incumbent_next))
    household_residual = next_households - (starting_households - exits + entries)
    scale = max(1.0, post_dependents, starting_households)
    if abs(float(base_residual)) > tolerance * scale:
        raise RuntimeError(f"base zero-entry transition mass residual is {base_residual:.6g}")
    if abs(child_residual) > tolerance * scale:
        raise RuntimeError(f"dependent-stock identity failed: residual={child_residual:.6g}")
    if abs(household_residual) > tolerance * scale:
        raise RuntimeError(f"household-stock identity failed: residual={household_residual:.6g}")

    ledger = SurvivingMaturationLedger(
        post_birth_dependent_units=post_dependents,
        dependent_deaths=dependent_deaths,
        surviving_maturations_model_child_units=mature_total,
        next_dependent_units=next_dependents,
        child_identity_residual=child_residual,
        starting_households=starting_households,
        household_exits=exits,
        domestic_formation_entries=float(np.sum(domestic)),
        outside_migration_entries=float(np.sum(outside)),
        next_households=next_households,
        household_identity_residual=household_residual,
        base_transition_mass_residual=float(base_residual),
        mature_model_child_units_by_location=mature,
        domestic_entries_by_location=domestic,
        outside_entries_by_location=outside,
        total_entries_by_location=total_entry,
        formation_contract_classification=str(contract.classification).strip().lower(),
        formation_contract_source=str(contract.source),
    )
    return incumbent_next, total_entry, ledger


def zero_post_origin_migration(
    net_migration_by_year: Mapping[int, Array],
    outside_entry_by_year: Mapping[int, Array],
    *,
    forecast_origin_year: int = 2023,
) -> tuple[dict[int, Array], dict[int, Array]]:
    """Return copied paths with every post-origin migration/entry cell zero."""

    origin = int(forecast_origin_year)
    migration: dict[int, Array] = {}
    outside: dict[int, Array] = {}
    for year, value in net_migration_by_year.items():
        array = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError(f"net migration for {year} is not finite")
        migration[int(year)] = np.zeros_like(array) if int(year) >= origin else array.copy()
    for year, value in outside_entry_by_year.items():
        array = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError(f"outside entry for {year} is not finite")
        outside[int(year)] = np.zeros_like(array) if int(year) >= origin else array.copy()
    return migration, outside


def assert_zero_post_origin_migration(
    net_migration_by_year: Mapping[int, Array],
    outside_entry_by_year: Mapping[int, Array],
    *,
    forecast_origin_year: int = 2023,
) -> None:
    """Reject aggregate-zero paths that conceal offsetting cell flows."""

    origin = int(forecast_origin_year)
    for label, path in (
        ("net migration", net_migration_by_year),
        ("outside entry", outside_entry_by_year),
    ):
        for year, value in path.items():
            if int(year) >= origin and np.any(np.asarray(value, dtype=float) != 0.0):
                raise ValueError(f"{label} is nonzero in at least one cell for {year}")
