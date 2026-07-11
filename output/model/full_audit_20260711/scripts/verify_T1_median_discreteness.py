#!/usr/bin/env python3
"""Verify finding T1: old nonhousing wealth/income median is b-grid-discrete.

No model solve. Builds P exactly as the overnight evaluation path does
(base_overrides + production profile + combined fixed spec + rouwenhorst
income + candidate theta), then:
  1. checks the retirement pension is constant across income states i and
     ages 65-75 (so ratio values = bg / const);
  2. checks the reproduced candidate medians (1.028539..., 0.378935...)
     coincide EXACTLY with bg-node / annual-pension values;
  3. reports local grid spacing at those nodes = the minimum possible jump
     size of the moment, and the implied objective jump
     |dL| = w * |2*gap*step + step^2|.
Writes JSON next to the target_audit outputs.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import intergen_housing_fertility.local_panel as lp
from intergen_housing_fertility.parameters import setup_parameters, apply_overrides
from intergen_housing_fertility.utils import make_grid
from intergen_housing_fertility.solver import annual_gross_income_at_state
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    production_profile_overrides,
)

MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
PERIOD_YEARS = 4.0
NB = 120

REPRO_DIR = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/full_audit_20260711/repro"
)
OUT = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/full_audit_20260711/target_audit/verify_T1_median_discreteness.json"
)

TARGET = 2.23046078
WEIGHT = 0.8


def build_P(theta: dict) -> object:
    income = lp.income_process_overrides(
        5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO
    )
    extra = production_profile_overrides()
    extra.update(
        {
            "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
            "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
            "eta_supply": np.array([1.75]),
            "normalize_bequest_utility": True,
            "max_iter_eq": PRODUCTION_MAX_ITER_EQ,
        }
    )
    overrides = {
        **lp.base_overrides(J=PRODUCTION_J, Nb=NB, n_house=5, max_iter_eq=PRODUCTION_MAX_ITER_EQ),
        **extra,
        **income,
        **theta,
    }
    P = setup_parameters()
    P = apply_overrides(P, overrides)
    return P


def main() -> None:
    results: dict = {}
    for label, fname in [
        ("current_bound_best", "current_repro_exact.json"),
        ("housing_first_best", "housing_repro_exact.json"),
    ]:
        rec = json.loads((REPRO_DIR / fname).read_text())
        theta = rec.get("theta") or rec.get("candidate", {}).get("theta")
        if theta is None:
            # reproduce_candidate stores the full record; theta key may vary
            theta = rec["theta"] if "theta" in rec else {}
        moments = rec.get("moments", {})
        med = float(moments["old_nonhousing_wealth_to_income_median_6575"])

        P = build_P(dict(theta))
        bg = make_grid(P)
        # ages 65-75 window as in solver stats: indices for da=4 lifecycle
        da = float(getattr(P, "da", getattr(P, "period_years", 4.0)))
        age0 = float(getattr(P, "age_start", 18.0))
        ages = [age0 + j * da for j in range(P.J)]
        jwin = [j for j in range(P.J) if ages[j] >= 65.0 - 1e-9 and ages[j] <= 75.0 + 1e-9]
        pens = {}
        for j in jwin:
            for i in range(P.I):
                pens[(i, j)] = annual_gross_income_at_state(P, i, j, 1.0)
        pen_vals = sorted(set(round(v, 12) for v in pens.values()))
        yj = float(pen_vals[0])

        ratios = bg / yj
        # nearest node to reproduced median
        k = int(np.argmin(np.abs(ratios - med)))
        node_ratio = float(ratios[k])
        exact = abs(node_ratio - med)
        step_up = float(ratios[k + 1] - ratios[k]) if k + 1 < len(ratios) else float("nan")
        step_dn = float(ratios[k] - ratios[k - 1]) if k > 0 else float("nan")
        gap = med - TARGET
        dloss_up = WEIGHT * abs(2 * gap * step_up + step_up**2)
        dloss_dn = WEIGHT * abs(-2 * gap * step_dn + step_dn**2)
        results[label] = {
            "reported_median": med,
            "pension_values_in_window_unique": pen_vals,
            "pension_constant_across_i_and_j": len(pen_vals) == 1,
            "j_window": jwin,
            "ages_window": [ages[j] for j in jwin],
            "nearest_grid_node_index": k,
            "nearest_grid_node_ratio": node_ratio,
            "abs_diff_median_vs_node": exact,
            "median_is_exact_grid_node": exact < 1e-10,
            "ratio_step_up": step_up,
            "ratio_step_down": step_dn,
            "gap_vs_target": gap,
            "loss_contribution": WEIGHT * gap * gap,
            "objective_jump_if_median_moves_one_node_up": dloss_up,
            "objective_jump_if_median_moves_one_node_down": dloss_dn,
            "Nb": NB,
            "b_max": float(getattr(P, "b_max", np.nan)),
        }
        # distinct ratio values overall
        results[label]["n_distinct_ratio_values"] = int(len(np.unique(np.round(ratios, 12))))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
