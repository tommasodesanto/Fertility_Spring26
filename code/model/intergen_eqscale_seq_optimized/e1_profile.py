"""Package-owned E1 eqscale/sequential benchmark contract."""
from __future__ import annotations

from types import SimpleNamespace
from typing import Any
import numpy as np

from .calibration import base_overrides, external_entry_wealth_overrides_1824, get_target_system
from .local_panel import income_process_overrides
from .production_profile import production_profile_overrides

E1_PROFILE_NAME = "eqscale_seq_e1_20260718"
E1_TARGET_SET = "candidate_replacement_income_disciplined_v1"
E1_J, E1_NB = 17, 120
E1_PRICE, E1_LOSS = 0.8312301242822234, 12.608491237749933
E1_THETA = {"H0": 7.4067615181119395, "alpha_cons": 0.6936410095385107,
    "beta": 0.9522043265603773, "chi": 1.063516193644096,
    "delta_alpha": 0.04229109060592314, "delta_alpha_jump": 0.062385084056824774,
    "gamma_e": 0.3830886809204491, "kappa_fert": 10.769451597500348,
    "psi_child": -0.023903107073990792, "tenure_choice_kappa": 0.01994450789165068,
    "theta0": 0.21496360295871036, "theta1": 0.19165791253321096, "theta_n": 0.0}
_SURVIVAL = {66.0: .9391263063710125, 70.0: .9184976343249724, 74.0: .8849521927812863, 78.0: .8300468061015381}

def _survival_schedule() -> np.ndarray:
    return np.array([_SURVIVAL.get(18.0 + 4.0*j, 1.0) for j in range(E1_J - 1)])

def e1_overrides(*, tight: bool = True, optimized: bool = True) -> dict[str, Any]:
    max_iter, tol = (40, 2.5e-5) if tight else (10, 1e-4)
    return {**base_overrides(J=E1_J, Nb=E1_NB, n_house=5, max_iter_eq=max_iter),
        **production_profile_overrides(),
        **income_process_overrides(5, "rouwenhorst", .20, .9601845894041878),
        **external_entry_wealth_overrides_1824(), **E1_THETA,
        "preference_spec": "eqscale", "sequential_births": True,
        "fecundity_omega1": .02, "fecundity_omega2": .134, "fecundity_terminal_age": 45.0,
        "transfer_floor_G0": 0.0, "transfer_floor_Gn": 0.0, "use_age_survival": True,
        "entry_wealth_censor_to_frontier": True, "bequest_spec": "linear_child_scale",
        "normalize_bequest_utility": True, "owner_ltv_taper": False,
        "max_iter_eq": max_iter, "tol_eq": tol, "q": (1.02 ** 4) - 1.0,
        "delta": 1.0 - (.989 ** 4), "eta_supply": np.array([1.75]), "lambda_d": 0.0,
        "debt_taper_start_age": 42.0, "debt_taper_end_age": 62.0,
        "survival_probs": _survival_schedule(),
        "markov_equilibrium_method": "direct_brent" if optimized else "legacy_damped",
        "use_numba_scatter": bool(optimized)}

def e1_target_system():
    system = get_target_system(E1_TARGET_SET)
    system.require_identified(13)
    return system
