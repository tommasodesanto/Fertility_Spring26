#!/usr/bin/env python3
"""BQ-1 verifier: quantify the effect of normalize_bequest_utility alone.

Solves the July-10 combined-spec candidate theta twice at a tiny grid
(Nb=32), identical in every respect except normalize_bequest_utility
(True = estimated model, False = what the un-patched tools/CLI would use),
and reports per-moment divergence. Audit-only; writes JSON to the audit
output folder.
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
NB = 32

REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")
RECORD = REPO / "output/model/full_audit_20260711/repro/current_repro_exact.json"
OUT = REPO / "output/model/full_audit_20260711/spec_audit/verify_bq1_bequest_flag_nb32.json"


def build_env(normalize: bool):
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
            "normalize_bequest_utility": bool(normalize),
        }
    )
    extra["max_iter_eq"] = int(PRODUCTION_MAX_ITER_EQ)
    return targets, weights, income, extra


def main() -> None:
    payload = json.loads(RECORD.read_text())
    theta = dict(payload["theta"])

    results = {}
    for label, normalize in (("normalized_true_estimated", True), ("normalized_false_tools", False)):
        targets, weights, income, extra = build_env(normalize)
        t0 = time.perf_counter()
        rec = lp.run_local_panel_case(
            0,
            {"label": label, "theta": theta},
            PRODUCTION_J,
            NB,
            5,
            int(PRODUCTION_MAX_ITER_EQ),
            income,
            targets,
            weights,
            extra,
        )
        rec["wall_seconds"] = time.perf_counter() - t0
        results[label] = rec

    targets, weights, _, _ = build_env(True)
    rows = []
    m_true = results["normalized_true_estimated"].get("moments") or {}
    m_false = results["normalized_false_tools"].get("moments") or {}
    for name in sorted(targets):
        vt = m_true.get(name)
        vf = m_false.get(name)
        rows.append(
            {
                "moment": name,
                "target": float(targets[name]),
                "weight": float(weights[name]),
                "model_normalized": None if vt is None else float(vt),
                "model_unnormalized": None if vf is None else float(vf),
                "diff_unnorm_minus_norm": (
                    None if vt is None or vf is None else float(vf) - float(vt)
                ),
            }
        )

    summary = {
        "Nb": NB,
        "theta": theta,
        "loss_normalized": float(results["normalized_true_estimated"]["rank_loss"]),
        "loss_unnormalized": float(results["normalized_false_tools"]["rank_loss"]),
        "p_eq_normalized": results["normalized_true_estimated"].get("p_eq"),
        "p_eq_unnormalized": results["normalized_false_tools"].get("p_eq"),
        "strict_normalized": results["normalized_true_estimated"].get("strict_converged"),
        "strict_unnormalized": results["normalized_false_tools"].get("strict_converged"),
        "moment_table": rows,
        "wall_normalized": results["normalized_true_estimated"]["wall_seconds"],
        "wall_unnormalized": results["normalized_false_tools"]["wall_seconds"],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(lp.jsonable(summary), indent=2, sort_keys=True))
    print(json.dumps(lp.jsonable(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
