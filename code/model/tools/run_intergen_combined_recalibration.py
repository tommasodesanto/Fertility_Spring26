#!/usr/bin/env python3
"""Bounded combined-specification calibration after the feasibility repair."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from intergen_housing_fertility.local_panel import (
    GLOBAL_DE_BOUNDS,
    global_unit_from_theta,
    run_local_polish,
    theta_from_global_unit,
)
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
)


ROOT = Path(__file__).resolve().parents[3]
MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
H0_SEARCH_UPPER = 20.0
MAX_MINUTES = 235.0
MAX_EVALS = 2000


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-theta", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--max-evals", type=int, default=int(os.getenv("COMBINED_MAX_EVALS", "60")))
    parser.add_argument("--minutes", type=float, default=float(os.getenv("COMBINED_MINUTES", "220")))
    parser.add_argument("--seed", type=int, default=int(os.getenv("COMBINED_SEED", "20260710")))
    parser.add_argument("--initial-step", type=float, default=float(os.getenv("COMBINED_INITIAL_STEP", "0.04")))
    parser.add_argument(
        "--method",
        choices=("pattern", "nelder-mead"),
        default=os.getenv("COMBINED_METHOD", "pattern"),
    )
    parser.add_argument("--start-mix", type=float, default=float(os.getenv("COMBINED_START_MIX", "0")))
    args = parser.parse_args()

    payload = json.loads(args.seed_theta.read_text())
    theta = dict(payload.get("theta", payload))
    theta["H0"] = float(payload.get("H0", theta.get("H0", 4.0)))
    start_mix = float(np.clip(args.start_mix, 0.0, 1.0))
    search_bounds = [*GLOBAL_DE_BOUNDS, ("H0", 1.0, H0_SEARCH_UPPER)]
    original_unit = global_unit_from_theta(theta, bounds=search_bounds)
    if original_unit is None:
        raise ValueError("seed record does not contain the complete 14-parameter vector")
    if start_mix > 0.0:
        rng = np.random.default_rng(args.seed)
        start_unit = (1.0 - start_mix) * original_unit + start_mix * rng.random(len(search_bounds))
        theta = theta_from_global_unit(np.clip(start_unit, 0.0, 1.0), bounds=search_bounds)
    else:
        start_unit = original_unit
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "combined_driver_config.json").write_text(
        json.dumps(
            {
                "method": args.method,
                "seed": args.seed,
                "start_mix": start_mix,
                "original_unit_vector": original_unit.tolist(),
                "start_unit_vector": np.asarray(start_unit, dtype=float).tolist(),
                "seed_theta_path": str(args.seed_theta),
            },
            indent=2,
            sort_keys=True,
        )
    )
    period_years = 4.0
    run_j = int(os.getenv("COMBINED_J", str(PRODUCTION_J)))
    run_nb = int(os.getenv("COMBINED_NB", str(PRODUCTION_SEARCH_NB)))
    profile_name = PRODUCTION_PROFILE_NAME if os.getenv("COMBINED_PROFILE", "1") == "1" else None
    fix_tenure_kappa = os.getenv("COMBINED_FIX_TENURE_KAPPA", "0").strip().lower() in {
        "1", "true", "yes", "on",
    }
    fixed_spec = {
        "q": (1.0 + 0.02) ** period_years - 1.0,
        "delta": 1.0 - (1.0 - 0.011) ** period_years,
        "eta_supply": np.array([1.75]),
        "normalize_bequest_utility": True,
        "lambda_d": 0.0,
        "debt_taper_start_age": 42.0,
        "debt_taper_end_age": 62.0,
    }
    summary = run_local_polish(
        args.outdir,
        max_evals=min(args.max_evals, MAX_EVALS),
        seed=args.seed,
        J=run_j,
        Nb=run_nb,
        n_house=5,
        max_iter_eq=PRODUCTION_MAX_ITER_EQ,
        minutes=min(args.minutes, MAX_MINUTES),
        income_states=5,
        income_process="rouwenhorst",
        income_annual_rho=MATCHED_ANNUAL_RHO,
        income_innovation_sd=MATCHED_ANNUAL_INNOVATION_SD,
        target_set=PRODUCTION_TARGET_SET,
        seed_theta=theta,
        method=args.method,
        initial_step=args.initial_step,
        min_step=0.003,
        shrink=0.5,
        profile_name=profile_name,
        fixed_theta={"tenure_choice_kappa": 0.0} if fix_tenure_kappa else None,
        normalize_bequest_utility=True,
        additional_search_bounds=[("H0", 1.0, H0_SEARCH_UPPER)],
        additional_targets={"aggregate_mean_occupied_rooms_18_85": ROOMS_TARGET},
        additional_weights={"aggregate_mean_occupied_rooms_18_85": 6.0},
        fixed_model_overrides=fixed_spec,
        progress=True,
    )
    best = summary.get("best")
    if not isinstance(best, dict) or not best.get("strict_converged", False):
        raise RuntimeError("combined calibration produced no strict-converged candidate")


if __name__ == "__main__":
    main()
