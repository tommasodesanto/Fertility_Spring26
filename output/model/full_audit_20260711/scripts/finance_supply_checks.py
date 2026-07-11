#!/usr/bin/env python3
"""Finance/supply spec-audit checks (audit task C).

Part 1 (no solve): defaults vs combined-spec finance parameters; exact 4-year
conversions; supply-curve identity on the saved repro records; required H0 to
hit the rooms target at the candidate equilibrium price.

Part 2 (tiny solve, Nb=32): pushes the current_bound_best theta through the
exact overnight evaluation chain (production profile + combined fixed spec +
Rouwenhorst income) at a small grid and verifies that q, delta, eta_supply,
H_own, hR_max, phi landed on P; that the renter cap and owner debt floor bind
as specified; and that the rooms moment equals aggregate demand / total mass.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import intergen_housing_fertility.local_panel as lp
from intergen_housing_fertility.calibration import base_overrides, extract_moments
from intergen_housing_fertility.parameters import apply_overrides, setup_parameters
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_TARGET_SET,
    production_profile_overrides,
)
from intergen_housing_fertility.solver import run_model_cp_dt

ROOT = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
AUDIT = ROOT / "output/model/full_audit_20260711"
OUT = AUDIT / "spec_audit"
OUT.mkdir(parents=True, exist_ok=True)

MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
PERIOD_YEARS = 4.0

report: dict = {}

# ---------------------------------------------------------------- Part 1
q_new = (1.0 + 0.02) ** PERIOD_YEARS - 1.0
delta_new = 1.0 - (1.0 - 0.011) ** PERIOD_YEARS

P0 = setup_parameters()
report["defaults"] = {
    "q_default": float(P0.q),
    "delta_default": float(P0.delta),
    "tau_H": float(P0.tau_H),
    "user_cost_rate_default": float(P0.user_cost_rate),
    "R_gross_default": float(P0.R_gross),
    "r_bar_default": float(P0.r_bar[0]),
    "H0_default": float(P0.H0[0]),
    "eta_supply_default": float(P0.eta_supply[0]),
    "H_own_default": P0.H_own.tolist(),
    "hR_max_default": float(P0.hR_max),
    "phi_default": np.atleast_1d(P0.phi).tolist(),
    "implied_annual_r_default": float((1.0 + P0.q) ** 0.25 - 1.0),
    "implied_annual_delta_default": float(1.0 - (1.0 - P0.delta) ** 0.25),
}
report["combined_spec"] = {
    "q_new": q_new,
    "delta_new": delta_new,
    "user_cost_rate_new": q_new + delta_new + float(P0.tau_H),
    "R_gross_new": 1.0 + q_new,
    "q_ratio_default_over_new": float(P0.q) / q_new,
    "delta_ratio_default_over_new": float(P0.delta) / delta_new,
    "delta4_exact_repr": repr(delta_new),
}

# apply overrides the way the overnight chain does and confirm recomputation
Pc = setup_parameters()
Pc = apply_overrides(
    Pc,
    {"q": q_new, "delta": delta_new, "eta_supply": np.array([1.75]), "H0": 9.999672059746388},
)
report["apply_overrides_check"] = {
    "user_cost_rate": float(Pc.user_cost_rate),
    "user_cost_matches_sum": bool(
        abs(Pc.user_cost_rate - (q_new + delta_new + Pc.tau_H)) < 1e-15
    ),
    "R_gross": float(Pc.R_gross),
    "xi_supply_follows_eta": Pc.xi_supply.tolist(),
    "H0_scalar_override_shape": list(np.shape(Pc.H0)),
    "H0_value": np.atleast_1d(Pc.H0).tolist(),
}

# supply identity on the saved repro records
uc_rate = q_new + delta_new + 0.01 * PERIOD_YEARS
records = {}
for name in ["current_repro_exact", "housing_repro_exact", "current_eq30", "current_pinit_lo"]:
    path = AUDIT / "repro" / f"{name}.json"
    if not path.exists():
        continue
    r = json.loads(path.read_text())
    p = float(np.asarray(r["p_eq"]).reshape(-1)[0])
    H0 = float(r["theta"]["H0"])
    supply_pred = H0 * (uc_rate * p / 0.16) ** 1.75
    rooms = float(r["moments"]["aggregate_mean_occupied_rooms_18_85"])
    req_H0 = ROOMS_TARGET / ((uc_rate * p / 0.16) ** 1.75)
    records[name] = {
        "p_eq": p,
        "H0": H0,
        "user_cost": uc_rate * p,
        "supply_predicted_from_formula": supply_pred,
        "rooms_moment": rooms,
        "abs_rel_gap_supply_vs_rooms": abs(supply_pred - rooms) / rooms,
        "H0_needed_for_target_at_this_price": req_H0,
        "rooms_gap_vs_target": rooms - ROOMS_TARGET,
    }
report["supply_identity_on_records"] = records

# ---------------------------------------------------------------- Part 2
targets, weights = lp.get_target_set(PRODUCTION_TARGET_SET)
targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
weights["aggregate_mean_occupied_rooms_18_85"] = 6.0
income = lp.income_process_overrides(
    5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO
)
extra = production_profile_overrides()
extra.update(
    {
        "q": q_new,
        "delta": delta_new,
        "eta_supply": np.array([1.75]),
        "normalize_bequest_utility": True,
    }
)
theta = json.loads((AUDIT / "repro" / "current_repro_exact.json").read_text())["theta"]

NB = 32
overrides = {
    **base_overrides(J=PRODUCTION_J, Nb=NB, n_house=5, max_iter_eq=PRODUCTION_MAX_ITER_EQ),
    **extra,
    **income,
    **theta,
}
sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
moments = extract_moments(sol, P)

p = float(np.asarray(p_eq).reshape(-1)[0])
hR = np.asarray(sol.hR_pol)
renter_slice = hR[:, 0, ...]
owner_slices = hR[:, 1:, ...]
bp = np.asarray(sol.bp_pol)
phi = float(np.atleast_1d(P.phi)[0])
floor_check = {}
for ten in range(1, 1 + P.n_house):
    floor = -phi * p * float(P.H_own[ten - 1])
    floor_check[f"rung_{ten}_H{P.H_own[ten - 1]:.0f}"] = {
        "debt_floor_-phi_pH": floor,
        "min_bp_pol": float(np.min(bp[:, ten, ...])),
        "violation": bool(np.min(bp[:, ten, ...]) < floor - 1e-9),
    }

supply_formula = float(P.H0[0]) * (P.user_cost_rate * p / float(P.r_bar[0])) ** float(
    P.xi_supply[0]
)
report["tiny_solve_nb32"] = {
    "P_q": float(P.q),
    "P_delta": float(P.delta),
    "P_user_cost_rate": float(P.user_cost_rate),
    "P_R_gross": float(P.R_gross),
    "P_eta_supply": np.atleast_1d(P.eta_supply).tolist(),
    "P_xi_supply": np.atleast_1d(P.xi_supply).tolist(),
    "P_H_own": np.asarray(P.H_own).tolist(),
    "P_hR_max": float(P.hR_max),
    "P_phi": np.atleast_1d(P.phi).tolist(),
    "P_r_bar": np.atleast_1d(P.r_bar).tolist(),
    "P_H0": np.atleast_1d(P.H0).tolist(),
    "p_eq": p,
    "max_renter_rooms_policy": float(np.max(renter_slice)),
    "renter_cap_respected": bool(np.max(renter_slice) <= P.hR_max + 1e-9),
    "owner_slice_of_hR_pol_all_zero": bool(np.max(np.abs(owner_slices)) == 0.0),
    "owner_debt_floor_checks": floor_check,
    "rooms_moment": float(moments["aggregate_mean_occupied_rooms_18_85"]),
    "aggregate_housing_demand": float(sol.aggregate_housing_demand),
    "total_mass": float(sol.total_mass),
    "rooms_equals_demand_over_mass": bool(
        abs(
            moments["aggregate_mean_occupied_rooms_18_85"]
            - sol.aggregate_housing_demand / sol.total_mass
        )
        < 1e-12
    ),
    "aggregate_housing_supply": float(sol.aggregate_housing_supply),
    "supply_formula_recomputed": supply_formula,
    "supply_matches_formula": bool(
        abs(supply_formula - float(sol.aggregate_housing_supply)) < 1e-10
    ),
    "market_residual": float(sol.best_max_abs_rel_excess),
}

# default-parameter drift demonstration: same theta WITHOUT the fixed spec
overrides_nofix = {
    **base_overrides(J=PRODUCTION_J, Nb=NB, n_house=5, max_iter_eq=PRODUCTION_MAX_ITER_EQ),
    **production_profile_overrides(),
    **income,
    **theta,
}
P_probe = apply_overrides(setup_parameters(), overrides_nofix)
report["silent_default_drift_probe"] = {
    "q_without_fixed_spec": float(P_probe.q),
    "delta_without_fixed_spec": float(P_probe.delta),
    "user_cost_rate_without_fixed_spec": float(P_probe.user_cost_rate),
    "note": "production_profile_overrides() alone leaves the old 4%/2% annual finance defaults in place",
}

(OUT / "finance_supply_checks.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
