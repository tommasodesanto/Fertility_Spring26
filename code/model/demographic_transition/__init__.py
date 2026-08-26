"""Person-cohort accounting for demographic–housing transitions."""

from .person_cohort_law import (
    CohortState,
    DemographicSteadyState,
    DemographicSteadyStateComponents,
    FlowLedger,
    advance_one_year,
    cohort_survivors,
    infer_net_migration_path,
    simulate_path,
    solve_demographic_steady_state,
    stationary_demographic_components,
)
from .household_head_bridge import (
    HeadMassRebalanceLedger,
    aggregate_heads_to_model_age_cells,
    rebalance_household_distribution_by_age,
)
from .four_year_bridge import (
    PersonBlockLedger,
    advance_person_state_block,
)
from .household_person_coupling import (
    HouseholdPersonBlockLedger,
    advance_household_person_block,
)

__all__ = [
    "CohortState",
    "DemographicSteadyState",
    "DemographicSteadyStateComponents",
    "FlowLedger",
    "advance_one_year",
    "cohort_survivors",
    "infer_net_migration_path",
    "simulate_path",
    "solve_demographic_steady_state",
    "stationary_demographic_components",
    "HeadMassRebalanceLedger",
    "aggregate_heads_to_model_age_cells",
    "rebalance_household_distribution_by_age",
    "PersonBlockLedger",
    "advance_person_state_block",
    "HouseholdPersonBlockLedger",
    "advance_household_person_block",
]
