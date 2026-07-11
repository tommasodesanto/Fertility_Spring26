#!/usr/bin/env python3
"""Audit-only reproduction of the July 10 combined-spec overnight candidates.

Rebuilds the evaluation environment of tmp/overnight_combined_20260710/
run_wave2_lateral.py exactly (production profile overrides + combined fixed
spec + Rouwenhorst income process + rooms target) and re-solves a saved
candidate theta through the same lp.run_local_panel_case path the overnight
used. Writes one JSON record per run. Does not touch production files.
"""

from __future__ import annotations

import argparse
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


def build_env(nb: int, max_iter_eq: int, tol_eq: float | None, p_init: float | None):
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
        }
    )
    extra["max_iter_eq"] = int(max_iter_eq)
    if tol_eq is not None:
        extra["tol_eq"] = float(tol_eq)
    if p_init is not None:
        extra["p_init_override"] = np.array([float(p_init)])
    return targets, weights, income, extra


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, required=True,
                        help="extracted report JSON holding theta (current_bound_best.json etc.)")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--label", type=str, required=True)
    parser.add_argument("--nb", type=int, default=120)
    parser.add_argument("--max-iter-eq", type=int, default=PRODUCTION_MAX_ITER_EQ)
    parser.add_argument("--tol-eq", type=float, default=None)
    parser.add_argument("--p-init", type=float, default=None)
    args = parser.parse_args()

    payload = json.loads(args.record.read_text())
    theta = dict(payload["theta"])
    stored_loss = float(payload.get("loss", float("nan")))
    stored_residual = float(payload.get("market_residual", float("nan")))

    targets, weights, income, extra = build_env(
        args.nb, args.max_iter_eq, args.tol_eq, args.p_init
    )
    t0 = time.perf_counter()
    record = lp.run_local_panel_case(
        0,
        {"label": args.label, "theta": theta},
        PRODUCTION_J,
        args.nb,
        5,
        int(args.max_iter_eq),
        income,
        targets,
        weights,
        extra,
    )
    record["audit_wall_seconds"] = time.perf_counter() - t0
    record["audit_config"] = {
        "nb": args.nb,
        "max_iter_eq": args.max_iter_eq,
        "tol_eq": args.tol_eq,
        "p_init": args.p_init,
    }
    record["stored_loss"] = stored_loss
    record["stored_residual"] = stored_residual
    record["loss_gap_vs_stored"] = float(record["rank_loss"]) - stored_loss

    # per-moment fit table from the reproduced moments
    fit = []
    moments = record.get("moments") or {}
    for name in sorted(targets):
        model_value = moments.get(name)
        target_value = targets[name]
        weight = weights[name]
        if model_value is None:
            fit.append({"moment": name, "error": "missing from extract_moments"})
            continue
        gap = float(model_value) - float(target_value)
        fit.append(
            {
                "moment": name,
                "target": float(target_value),
                "model": float(model_value),
                "gap": gap,
                "weight": float(weight),
                "loss_contribution": float(weight) * gap * gap,
            }
        )
    record["audit_fit"] = fit

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(lp.jsonable(record), indent=2, sort_keys=True))
    print(
        f"{args.label}: loss={record['rank_loss']:.9f} stored={stored_loss:.9f} "
        f"gap={record['loss_gap_vs_stored']:.3e} residual={record['market_residual']:.3e} "
        f"strict={record['strict_converged']} wall={record['audit_wall_seconds']:.1f}s"
    )


if __name__ == "__main__":
    main()
