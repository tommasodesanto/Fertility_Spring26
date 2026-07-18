"""Named, reproducible production configuration for the July 9 repair cycle."""

from __future__ import annotations

from typing import Any

import numpy as np

from .calibration import PERIOD_YEARS
from .parameters import make_persistent_transition_matrix


PRODUCTION_PROFILE_NAME = "intergen_july9_repair_v1"
PRODUCTION_J = 17
PRODUCTION_SEARCH_NB = 120
PRODUCTION_VERIFY_NB = 240
PRODUCTION_H_OWN = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
PRODUCTION_RENTER_CAP = 6.0
PRODUCTION_INCOME_STATES = 5
PRODUCTION_TARGET_SET = "candidate_replacement_post_audit_v1"
PRODUCTION_MAX_ITER_EQ = 10
PRODUCTION_HOUSING_EVENT_HORIZON = 0
PRODUCTION_INCOME_AGE_BREAKS = np.array([22.0, 26.0, 34.0, 46.0, 58.0])
PRODUCTION_INCOME_AGE_VALUES = np.array([0.650, 0.850, 1.000, 0.985, 0.935])
PRODUCTION_Z_GRID = np.array([0.60, 0.80, 1.00, 1.20, 1.40])
PRODUCTION_Z_WEIGHTS = np.array([0.10, 0.20, 0.40, 0.20, 0.10])
PRODUCTION_INCOME_PERSISTENCE = 0.85
PRODUCTION_B_MIN = -12.0
PRODUCTION_B_MAX = 30.0
PRODUCTION_B_CORE_LO = -5.0
PRODUCTION_B_CORE_HI = 7.0
PRODUCTION_B_MID_HI = 15.0
PRODUCTION_B_FRAC_LOW = 0.08
PRODUCTION_B_FRAC_CORE = 0.72
PRODUCTION_B_FRAC_MID = 0.12
PRODUCTION_B_GRID_POWER = 1.5

# These are the source-controlled bounds that generated the reported July runs.
PRODUCTION_SEARCH_BOUNDS: tuple[tuple[str, float, float], ...] = (
    ("beta_annual", 0.940, 0.995),
    ("alpha_cons", 0.400, 0.950),
    ("c_bar_0", 0.020 * PERIOD_YEARS, 0.320 * PERIOD_YEARS),
    ("c_bar_n", 0.050, 1.500),
    ("h_bar_0", 1.000, 6.000),
    ("h_bar_jump", 0.050, 2.500),
    ("h_bar_n", 0.020, 2.000),
    ("psi_child", 0.000, 0.350),
    ("kappa_fert", 1.000, 12.000),
    ("tenure_choice_kappa", 0.000, 0.120),
    ("chi", 0.400, 1.150),
    ("theta0", 0.000, 2.000),
    ("theta_n", 0.000, 1.500),
)

FROZEN_SOURCE_THETA: dict[str, float] = {
    "alpha_cons": 0.6195578033469402,
    "beta": 0.7808047585225533,
    "c_bar_0": 1.2793341787106423,
    "c_bar_n": 0.11617574349164181,
    "chi": 1.1496942798885048,
    "h_bar_0": 1.0,
    "h_bar_jump": 2.1198232149690575,
    "h_bar_n": 1.3819945583480897,
    "kappa_fert": 1.0188742355464353,
    "psi_child": 0.3490340757304799,
    "tenure_choice_kappa": 0.0,
    "theta0": 0.0014981114599317271,
    "theta_n": 0.9722009696445806,
}
FROZEN_SOURCE_RECORD_NAME = "nb120_cap60_onehour_best.json"
FROZEN_SOURCE_RANK_LOSS = 6.3052679727526595
FROZEN_SOURCE_MARKET_RESIDUAL = 0.00019303196737969033
FROZEN_SOURCE_STRICT_CONVERGED = False

FROZEN_SOURCE_REPRO_SWITCHES: dict[str, Any] = {
    "use_postdecision_current_distribution": False,
    "legacy_entry_income_peak": False,
    "propagate_birth_entry_grant": False,
    "housing_event_horizon": PRODUCTION_HOUSING_EVENT_HORIZON,
}
REPAIRED_TIMING_SWITCHES: dict[str, Any] = {
    "use_postdecision_current_distribution": True,
    "legacy_entry_income_peak": False,
    "propagate_birth_entry_grant": True,
    "housing_event_horizon": PRODUCTION_HOUSING_EVENT_HORIZON,
}

PRELIMINARY_FIXED_BETA_ANNUAL = 0.94
PRELIMINARY_FIXED_THETA_N = FROZEN_SOURCE_THETA["theta_n"]
PRELIMINARY_CHI_RUNGS: tuple[float, ...] = (0.40, 0.55, 0.70, 0.85, 1.00, 1.15)


def production_profile_metadata() -> dict[str, Any]:
    runtime = {
        key: value.tolist() if isinstance(value, np.ndarray) else value
        for key, value in production_profile_overrides().items()
    }
    return {
        "name": PRODUCTION_PROFILE_NAME,
        "J": PRODUCTION_J,
        "search_Nb": PRODUCTION_SEARCH_NB,
        "verification_Nb": PRODUCTION_VERIFY_NB,
        "H_own": PRODUCTION_H_OWN.tolist(),
        "renter_cap": PRODUCTION_RENTER_CAP,
        "income_states": PRODUCTION_INCOME_STATES,
        "target_set": PRODUCTION_TARGET_SET,
        "max_iter_eq": PRODUCTION_MAX_ITER_EQ,
        "runtime_overrides": runtime,
        "frozen_source": {
            "record_name": FROZEN_SOURCE_RECORD_NAME,
            "theta": dict(FROZEN_SOURCE_THETA),
            "rank_loss": FROZEN_SOURCE_RANK_LOSS,
            "market_residual": FROZEN_SOURCE_MARKET_RESIDUAL,
            "strict_converged": FROZEN_SOURCE_STRICT_CONVERGED,
        },
        "bounds": [
            {"name": name, "lower": lower, "upper": upper}
            for name, lower, upper in PRODUCTION_SEARCH_BOUNDS
        ],
    }


def production_profile_overrides() -> dict[str, Any]:
    """Return model objects that must not drift with generic run arguments."""

    return {
        "H_own": PRODUCTION_H_OWN.copy(),
        "hR_max": PRODUCTION_RENTER_CAP,
        "child_housing_spec": "jump_plus_linear",
        "income_age_breaks": PRODUCTION_INCOME_AGE_BREAKS.copy(),
        "income_age_values": PRODUCTION_INCOME_AGE_VALUES.copy(),
        "normalize_income_profile": True,
        "legacy_entry_income_peak": False,
        "use_income_types": True,
        "income_type_transition": "markov",
        "z_grid": PRODUCTION_Z_GRID.copy(),
        "z_weights": PRODUCTION_Z_WEIGHTS.copy(),
        "income_shock_persistence": PRODUCTION_INCOME_PERSISTENCE,
        "Pi_z": make_persistent_transition_matrix(
            PRODUCTION_Z_WEIGHTS, PRODUCTION_INCOME_PERSISTENCE
        ),
        "b_min": PRODUCTION_B_MIN,
        "b_max": PRODUCTION_B_MAX,
        "b_core_lo": PRODUCTION_B_CORE_LO,
        "b_core_hi": PRODUCTION_B_CORE_HI,
        "b_mid_hi": PRODUCTION_B_MID_HI,
        "b_frac_low": PRODUCTION_B_FRAC_LOW,
        "b_frac_core": PRODUCTION_B_FRAC_CORE,
        "b_frac_mid": PRODUCTION_B_FRAC_MID,
        "b_grid_power": PRODUCTION_B_GRID_POWER,
        "entry_wealth_spread_nodes": 1,
        "interp_method": "linear",
        "max_iter_eq": PRODUCTION_MAX_ITER_EQ,
        "tol_eq": 1e-4,
        "scalar_market_refine": True,
        "scalar_market_refine_method": "brent",
        "scalar_market_refine_iter": 16,
        "scalar_market_refine_max_expand": 8,
        "housing_event_horizon": PRODUCTION_HOUSING_EVENT_HORIZON,
        "use_pti_constraint": False,
        "parent_dp_waiver": False,
        "birth_dp_grant": False,
        "birth_entry_grant": False,
        "birth_entry_grant_amount": 0.0,
        "propagate_birth_entry_grant": True,
        "use_postdecision_current_distribution": True,
    }


def comparison_arm_switches(label: str) -> dict[str, Any]:
    arms = {
        "frozen_source_repro": FROZEN_SOURCE_REPRO_SWITCHES,
        "repaired_timing": REPAIRED_TIMING_SWITCHES,
    }
    try:
        return dict(arms[str(label)])
    except KeyError as exc:
        raise ValueError(f"unknown comparison arm: {label}") from exc


def validate_frozen_source_theta(theta: dict[str, float]) -> None:
    received = {str(key): float(value) for key, value in theta.items()}
    if received != FROZEN_SOURCE_THETA:
        mismatches = {
            key: {"expected": FROZEN_SOURCE_THETA.get(key), "received": received.get(key)}
            for key in sorted(set(FROZEN_SOURCE_THETA) | set(received))
            if FROZEN_SOURCE_THETA.get(key) != received.get(key)
        }
        raise ValueError(
            f"comparison requires the exact frozen source theta; mismatches: {mismatches}"
        )


def validate_production_profile(
    profile_name: str,
    *,
    J: int,
    Nb: int,
    n_house: int,
    income_states: int,
    target_set: str,
    max_iter_eq: int,
    stage: str = "search",
) -> None:
    if str(profile_name) != PRODUCTION_PROFILE_NAME:
        raise ValueError(f"unknown production profile: {profile_name}")
    stage_clean = str(stage).strip().lower()
    expected_nb = PRODUCTION_SEARCH_NB if stage_clean == "search" else PRODUCTION_VERIFY_NB
    if stage_clean not in {"search", "verify"}:
        raise ValueError(f"unknown production profile stage: {stage}")
    received = {
        "J": int(J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "income_states": int(income_states),
        "target_set": str(target_set),
        "max_iter_eq": int(max_iter_eq),
    }
    expected = {
        "J": PRODUCTION_J,
        "Nb": expected_nb,
        "n_house": len(PRODUCTION_H_OWN),
        "income_states": PRODUCTION_INCOME_STATES,
        "target_set": PRODUCTION_TARGET_SET,
        "max_iter_eq": PRODUCTION_MAX_ITER_EQ,
    }
    if received != expected:
        raise ValueError(
            f"{PRODUCTION_PROFILE_NAME}/{stage_clean} configuration drift: "
            f"expected {expected}, received {received}"
        )
