"""Complete, package-owned M5 benchmark contract."""

from __future__ import annotations

from typing import Any

import numpy as np

from .calibration import base_overrides, external_entry_wealth_overrides_1824, get_target_system
from .local_panel import income_process_overrides
from .production_profile import production_profile_overrides


M5_PROFILE_NAME = "intergen_m5_income_disciplined_20260717"
M5_TARGET_SET = "candidate_replacement_income_disciplined_v1"
M5_J = 17
M5_NB = 120
M5_PRICE = 0.7589039340587314
M5_LOSS = 9.044422069071352
# Same theta under the optimized strict-root evaluator in the completed Torch
# promotion panel. The gap from M5_LOSS comes from nearby accepted prices.
M5_OPTIMIZED_STRICT_ROOT_LOSS = 9.14507565057346
M5_OPTIMIZED_STRICT_ROOT_PRICE = 0.7588693705046422
M5_ANNUAL_RHO = 0.9601845894041878
M5_ANNUAL_INNOVATION_SD = 0.06453733259357768
M5_POSTRETIREMENT_SURVIVAL = {
    66.0: 0.9391263063710125,
    70.0: 0.9184976343249724,
    74.0: 0.8849521927812863,
    78.0: 0.8300468061015381,
}
M5_THETA = {
    "H0": 8.645577738994191,
    "alpha_cons": 0.5912263839788428,
    "beta": 0.9653857228201739,
    "c_bar_0": 1.2595937731478881,
    "c_bar_n": 0.39294466990660537,
    "chi": 1.113265753289088,
    "h_bar_0": 0.3915418909982497,
    "h_bar_jump": 1.6021182503934808,
    "h_bar_n": 0.16436541806228172,
    "kappa_fert": 2.0433862404293337,
    "psi_child": 0.19791383267507004,
    "tenure_choice_kappa": 0.010017488787185433,
    "theta0": 0.31176482840616865,
    "theta1": 0.3972817386975299,
    "theta_n": 0.0,
}


def _survival_schedule() -> np.ndarray:
    schedule = np.ones(M5_J - 1, dtype=float)
    for transition in range(M5_J - 1):
        age = 18.0 + 4.0 * transition
        if age in M5_POSTRETIREMENT_SURVIVAL:
            schedule[transition] = M5_POSTRETIREMENT_SURVIVAL[age]
    return schedule


def m5_overrides(*, tight: bool = True, optimized: bool = True) -> dict[str, Any]:
    """Return the complete M5 solver contract without importing a runner."""

    max_iter_eq = 40 if tight else 10
    tol_eq = 2.5e-5 if tight else 1e-4
    overrides = {
        **base_overrides(J=M5_J, Nb=M5_NB, n_house=5, max_iter_eq=max_iter_eq),
        **production_profile_overrides(),
        **income_process_overrides(
            5,
            "rouwenhorst",
            M5_ANNUAL_INNOVATION_SD,
            M5_ANNUAL_RHO,
        ),
        **external_entry_wealth_overrides_1824(),
        **M5_THETA,
        "bequest_spec": "linear_child_scale",
        "normalize_bequest_utility": True,
        "owner_ltv_taper": False,
        "owner_ltv_taper_start_age": 66.0,
        "owner_ltv_taper_end_age": 82.0,
        "owner_ltv_terminal_share": 0.0,
        "use_age_survival": True,
        "entry_wealth_censor_to_frontier": True,
        "max_iter_eq": max_iter_eq,
        "tol_eq": tol_eq,
        "q": (1.0 + 0.02) ** 4.0 - 1.0,
        "delta": 1.0 - (1.0 - 0.011) ** 4.0,
        "eta_supply": np.array([1.75]),
        "lambda_d": 0.0,
        "debt_taper_start_age": 42.0,
        "debt_taper_end_age": 62.0,
        "survival_probs": _survival_schedule(),
        "markov_equilibrium_method": "direct_brent" if optimized else "legacy_damped",
        "use_numba_scatter": bool(optimized),
    }
    return overrides


def m5_target_system():
    system = get_target_system(M5_TARGET_SET)
    system.require_identified(free_parameter_count=14)
    return system
