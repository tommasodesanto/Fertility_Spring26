#!/usr/bin/env python3
"""INC-2 verifier: quantify contamination if the July-10 combined candidate is
evaluated through the checked-in tool path (production profile + default
'current' income process, old q/delta/eta/bequest) instead of the combined spec.
Tiny grid Nb=40 for both arms; same theta, same target system.
Audit-only; writes JSON next to this script's output dir.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

import intergen_housing_fertility.local_panel as lp
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_TARGET_SET,
    production_profile_overrides,
)

MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
PERIOD_YEARS = 4.0
NB = 40
MAX_ITER_EQ = 10

THETA = {
    "H0": 9.999672059746388,
    "alpha_cons": 0.6731064483929385,
    "beta": 0.8411423842686492,
    "c_bar_0": 1.28,
    "c_bar_n": 0.4559787014242588,
    "chi": 1.15,
    "h_bar_0": 1.0,
    "h_bar_jump": 1.4758939442666474,
    "h_bar_n": 0.984106767291167,
    "kappa_fert": 1.8906304000000005,
    "psi_child": 0.26907999999999993,
    "tenure_choice_kappa": 0.0,
    "theta0": 0.1318301350511569,
    "theta_n": 0.7679202774006243,
}


def targets_weights():
    targets, weights = lp.get_target_set(PRODUCTION_TARGET_SET)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = 6.0
    return targets, weights


def run(label, income, extra, targets, weights):
    t0 = time.perf_counter()
    rec = lp.run_local_panel_case(
        0,
        {"label": label, "theta": dict(THETA)},
        PRODUCTION_J,
        NB,
        5,
        MAX_ITER_EQ,
        income,
        targets,
        weights,
        extra,
    )
    rec["wall_seconds"] = time.perf_counter() - t0
    return rec


def main():
    targets, weights = targets_weights()

    # Arm A: correct combined spec (matched Rouwenhorst + fixed spec).
    income_a = lp.income_process_overrides(
        5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO
    )
    extra_a = production_profile_overrides()
    extra_a.update(
        {
            "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
            "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
            "eta_supply": np.array([1.75]),
            "normalize_bequest_utility": True,
        }
    )

    # Arm B: what checked-in tools/cli defaults would use — production profile
    # overrides with the old 5-point z grid and NO combined fixed spec.
    income_b = lp.income_process_overrides(5)  # default process="current"
    extra_b = production_profile_overrides()

    # Arm C: cli.py bare "--income-process rouwenhorst" defaults (rho=0.90, sd=0.20),
    # still no combined fixed spec.
    income_c = lp.income_process_overrides(5, "rouwenhorst")
    extra_c = production_profile_overrides()

    out = {}
    for label, income, extra in [
        ("A_combined_correct", income_a, extra_a),
        ("B_tool_default_current", income_b, extra_b),
        ("C_cli_bare_rouwenhorst", income_c, extra_c),
    ]:
        rec = run(label, income, extra, targets, weights)
        keep = {
            "rank_loss": rec.get("rank_loss"),
            "market_residual": rec.get("market_residual"),
            "strict_converged": rec.get("strict_converged"),
            "status": rec.get("status"),
            "p_eq": rec.get("p_eq"),
            "moments": {k: rec.get("moments", {}).get(k) for k in sorted(targets)},
            "wall_seconds": rec.get("wall_seconds"),
        }
        out[label] = keep
        print(label, "loss=", keep["rank_loss"], "resid=", keep["market_residual"],
              "strict=", keep["strict_converged"], f"wall={keep['wall_seconds']:.1f}s")

    # per-moment side-by-side
    rows = []
    for name in sorted(targets):
        row = {"moment": name, "target": targets[name], "weight": weights[name]}
        for label in out:
            row[label] = out[label]["moments"].get(name)
        rows.append(row)
    out["_fit_table"] = rows
    out["_theta"] = THETA
    out["_config"] = {"Nb": NB, "max_iter_eq": MAX_ITER_EQ, "J": PRODUCTION_J}

    dest = Path(__file__).resolve().parent.parent / "spec_audit" / "inc2_severity_spec_gap_nb40.json"
    dest.write_text(json.dumps(lp.jsonable(out), indent=2, sort_keys=True))
    print("wrote", dest)


if __name__ == "__main__":
    main()
