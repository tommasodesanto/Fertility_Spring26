"""Timing- and wealth-definition-repaired fourteen-moment July 2026 contract."""

from __future__ import annotations

from typing import Any

from .m5_profile import M5_THETA, m5_overrides
from .target_system import TargetSystem


NEW_MOMENT_PROFILE_NAME = "intergen_new_moment_ledger_gross_wealth_income_20260723"
NEW_MOMENT_PROFILE_RUNNABLE = True

# Cross-sectional wealth and living-old wealth use the beginning-of-period
# state. The bequest flow uses post-saving resources at the death node. The
# aggregate saving target is the matched 2005-2019 PSID gross/gross ratio from
# code/data/psid_followup_mar2026/output/aggregate_wealth_earnings_audit/.
NEW_MOMENT_TARGETS = {
    "aggregate_wealth_to_annual_gross_labor_earnings": 6.87307721116358,
    "annual_bequest_flow_to_aggregate_wealth": 0.0088,
    "old_total_wealth_to_annual_income_p90_p50_7684": 3.44811075444552,
    "childless_renter_rent_expenditure_slope": 0.264688845713763,
    "childless_renter_intercept_at_mean_price_model_units": 0.106415445844303,
    "bottom_quintile_childless_renter_mean_rooms": 3.89742817329093,
    "one_shot_parity_consumption_coefficient": 0.0752004729804767,
    "housing_increment_0to1": 0.664435,
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": 0.367700,
    "tfr": 1.918,
    "childless_rate": 0.188,
    "own_rate": 0.575472,
    "model_feasible_four_year_tenure_brier": 0.117612603915696,
    "aggregate_mean_occupied_rooms_18_85": 5.779970,
}

# Equal squared percentage deviations.  This is deliberately transparent for
# the preliminary three-hour search and avoids mixing units across the new
# auxiliary estimates before a joint covariance matrix is available.
NEW_MOMENT_WEIGHTS = {
    name: 1.0 / (target * target) for name, target in NEW_MOMENT_TARGETS.items()
}

NEW_MOMENT_SEED = {
    **M5_THETA,
    # Timing-repaired v2 winner, used only as the continuation seed for the
    # gross/gross target revision.
    "beta": 0.997459215934036,
    "alpha_cons": 0.6274604033546914,
    "c_bar_0": 0.5168123052802193,
    "c_bar_n": 0.1581965720667883,
    "h_bar_0": 0.5772612042418263,
    "h_bar_jump": 0.04817446833955657,
    "h_bar_n": 0.9259334338411668,
    "psi_child": -1.038325531838713,
    "kappa_fert": 5.492875182493313,
    "chi": 1.0802049637994615,
    "H0": 5.137297717207506,
    "theta0": 0.020734580571044135,
    "theta1": 0.12370835450345574,
    "tenure_choice_kappa": 0.018730093416335268,
}


def new_moment_target_system() -> TargetSystem:
    system = TargetSystem.from_mappings(
        "new_moment_ledger_14_v3_recent_gross_wealth_income_relative_loss",
        NEW_MOMENT_TARGETS,
        NEW_MOMENT_WEIGHTS,
    )
    system.require_identified(free_parameter_count=14)
    return system


def new_moment_overrides(*, tight: bool = True, optimized: bool = True) -> dict[str, Any]:
    """Use the current M5 model and numerical contract with the new seed."""

    return {
        **m5_overrides(tight=tight, optimized=optimized),
        **NEW_MOMENT_SEED,
    }
