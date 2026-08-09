"""E5 sequential-fertility target contract, revised 2026-08-09.

``tfr`` is retained as a legacy key, but denotes completed fertility (mean
children ever born for women ages 40--44), with the literal 3+ top-bin weight
used by the E4 L4 convention.  Rows with an available standard error use the
diagonal inverse-variance weight; other rows use provisional target-scale
weights.  The first-birth rooms target and standard error come from the
household-and-year-fixed-effect PSID specification restored on 2026-08-09.
The borrowed bequest-flow row uses the explicitly synthetic provisional
target-scale weight ``1 / 0.0088**2`` pending the weights pass.
"""
from __future__ import annotations

from typing import Any

from .calibration import OLD_NONLOCATION_TARGETS
from .externals import flhsv_income_overrides
from .target_system import TargetSystem

E5_PROFILE_NAME = "eqscale_seq_e5_idfe_20260809"
E5_TARGET_SET = "e5_idfe_review_20260809"
E5_TARGET_PROVENANCE = {
    "housing_increment_0to1": {
        "source": "PSID first-birth Sun--Abraham event study",
        "builder": "code/data/psid_followup_mar2026/sa_rooms_first_birth_one_variant_v1.do",
        "event_time": 3,
        "fixed_effects": ["ID", "year"],
        "covariates": ["i.AGEREP", "i.EDUYEAR"],
        "cluster": "ID",
        "estimate": 0.80494368,
        "standard_error": 0.16728361,
        "regenerated": "2026-08-09",
    }
}
E5_TARGETS: dict[str, float | None] = {
    "tfr": 1.918,
    "childless_rate": 0.188,
    # NCHS natality archive, cohorts 1979-84 (code/data/nchs_natality_timing)
    "mean_age_first_birth": 25.310560799362,
    "share_first_births_age30plus": 0.270062376851342,
    # PSID Sun--Abraham first-birth rooms response at k=+3, with household
    # and year fixed effects; regenerated 2026-08-09 from the active builder.
    "housing_increment_0to1": 0.80494368,
    "prime30_55_parent_3plus_minus_1to2_mean_rooms": 0.36769955881,
    "own_family_gap": OLD_NONLOCATION_TARGETS["own_family_gap"],
    "own_rate": 0.575472,
    "aggregate_mean_occupied_rooms_18_85": 5.779970481941968,
    "aggregate_wealth_to_annual_gross_labor_earnings": 6.8731,
    "annual_bequest_flow_to_aggregate_wealth": 0.0088,
    "old_total_wealth_to_annual_income_p90_p50_7684": 3.44811075444552,
}


# One weight rule (provisional, revisited at the freeze): every row enters the
# loss as (gap / SE)^2.  Measured bootstrap SEs where they exist; the two
# timing SEs are declared from the primary-vs-sensitivity window spread; all
# remaining rows carry a declared synthetic SE of 5 percent of the target.
E5_MEASURED_SES: dict[str, float] = {
    "tfr": 0.026483780917,                     # CPS Jun 2024 bootstrap
    "childless_rate": 0.007629200330,          # CPS Jun 2024 bootstrap
    "mean_age_first_birth": 0.25,              # declared, window spread 0.23
    "share_first_births_age30plus": 0.008,     # declared, window spread 0.0076
    "housing_increment_0to1": 0.16728361,      # PSID ID-FE Sun--Abraham, clustered by household
    "aggregate_wealth_to_annual_gross_labor_earnings": 0.3988,   # PSID bootstrap
    "old_total_wealth_to_annual_income_p90_p50_7684": 0.1325,    # PSID bootstrap
}


def _weight(name: str) -> float:
    value = E5_TARGETS[name]
    if value is None:
        return 1.0
    se = E5_MEASURED_SES.get(name, 0.05 * abs(float(value)))
    return 1.0 / se**2


E5_WEIGHTS: dict[str, float] = {name: _weight(name) for name in E5_TARGETS}

E5_DOMAIN: tuple[tuple[str, float, float, str], ...] = (
    ("beta_annual", 0.94, 0.9995, "discount"),
    ("delta_alpha", 0.0, 0.25, "softzero"),
    ("delta_alpha_jump", 0.0, 0.25, "softzero"),
    ("psi_child", -3.0, 3.0, "asinh"),
    ("kappa_fert", 0.02, 50.0, "log"),
    ("kappa_fert_continuation", 0.02, 50.0, "log"),
    ("chi", 0.10, 5.0, "log"),
    ("H0", 0.20, 80.0, "log"),
    ("theta0", 0.0, 8.0, "softzero"),
    ("theta1", 0.02, 16.0, "log"),
)
E5_FIXED = {"theta_n": 0.0, "alpha_cons": 0.733, "tenure_choice_kappa": 0.0}
E5_EXTERNAL_METADATA = {
    "alpha_cons": {"value": 0.733, "source": "CEX childless-cash-renter expenditure slope, 2019-2023"},
    "tenure_choice_kappa": {"value": 0.0, "source": "author-pinned deterministic tenure choice"},
    "theta_n": {"value": 0.0, "source": "author-pinned child-invariant bequest utility"},
}


def e5_target_system() -> TargetSystem:
    """Build E5, refusing any target row that is still pending."""
    pending = [name for name, value in E5_TARGETS.items() if value is None]
    if pending:
        raise ValueError(f"E5 target values are pending: {pending}")
    targets = {name: float(value) for name, value in E5_TARGETS.items()}
    system = TargetSystem.from_mappings(E5_TARGET_SET, targets, E5_WEIGHTS)
    system.require_identified(10)
    return system


def e5_overrides() -> dict[str, Any]:
    """Return the signed E5 external conventions (E4 L4, power, FL x HSV)."""
    return {
        **E5_FIXED,
        "n_parity": 4,
        "fertility_units": "literal_topcode",
        "tfr_top_bin_weight": 3.602359422009,  # CPS Jun 2024 capped top-bin mean (SE 0.0279)
        "entrant_conversion_factor": 0.5,
        "child_bin_high_cutoff": 3,
        "eqscale_form": "power",
        "gamma_e": 0.0,
        **flhsv_income_overrides(),
    }
