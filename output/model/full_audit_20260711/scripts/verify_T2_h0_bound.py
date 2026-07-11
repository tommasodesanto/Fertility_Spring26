#!/usr/bin/env python3
"""Verify finding T2: is the H0=10 upper bound locally binding, and does more
supply close the aggregate rooms gap?

Evaluates the current_bound_best theta with H0 overridden on a small grid
(Nb=40, audit-only) holding all other parameters fixed. Reports loss, the
rooms moment, and the big ownership moments for each H0.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

import intergen_housing_fertility.local_panel as lp
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_TARGET_SET,
    production_profile_overrides,
)

MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
PERIOD_YEARS = 4.0
NB = 40

REPORT = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/"
    "combined_recalibration/overnight_20260710_report/current_bound_best.json"
)
OUTDIR = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/"
    "full_audit_20260711/verify_T2"
)

KEY_MOMENTS = [
    "aggregate_mean_occupied_rooms_18_85",
    "own_rate",
    "own_rate_2534",
    "old_age_own_rate",
    "prime30_55_childless_renter_mean_rooms",
    "prime30_55_childless_owner_minus_renter_mean_rooms",
]


def build_env():
    targets, weights = lp.get_target_set(PRODUCTION_TARGET_SET)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = 6.0
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
            "max_iter_eq": int(PRODUCTION_MAX_ITER_EQ),
        }
    )
    return targets, weights, income, extra


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    payload = json.loads(REPORT.read_text())
    base_theta = dict(payload["theta"])
    targets, weights, income, extra = build_env()

    results = []
    for h0 in [8.56, 10.0, 12.0, 16.0]:
        theta = dict(base_theta)
        theta["H0"] = h0
        t0 = time.perf_counter()
        rec = lp.run_local_panel_case(
            0,
            {"label": f"T2_H0_{h0}", "theta": theta},
            PRODUCTION_J,
            NB,
            5,
            int(PRODUCTION_MAX_ITER_EQ),
            income,
            targets,
            weights,
            extra,
        )
        wall = time.perf_counter() - t0
        moments = rec.get("moments") or {}
        row = {
            "H0": h0,
            "loss": float(rec["rank_loss"]),
            "strict": bool(rec.get("strict_converged")),
            "residual": float(rec.get("market_residual", float("nan"))),
            "wall_s": wall,
        }
        for m in KEY_MOMENTS:
            row[m] = float(moments.get(m, float("nan")))
        rooms_gap = row["aggregate_mean_occupied_rooms_18_85"] - ROOMS_TARGET
        row["rooms_gap"] = rooms_gap
        row["rooms_loss_contrib"] = 6.0 * rooms_gap**2
        results.append(row)
        print(json.dumps(row))

    (OUTDIR / "h0_sweep_nb40.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
