"""Small diagnostic calibration driver for the one-market model scaffold."""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np

from .solver import run_model_cp_dt


PERIOD_YEARS = 4.0
AGE_START = 22.0
FERTILITY_START_AGE = 26.0
FERTILITY_END_AGE = 42.0
RETIREMENT_AGE = 66.0
CHILD_MATURITY_AGE = 18.0


CORE_TARGETS = {
    "own_rate": 0.575,
    "young_owner_rate": 0.25,
    "old_owner_rate": 0.764,
    "mean_completed_fertility": 0.95,
    "childless_rate": 0.20,
}


OLD_NONLOCATION_TARGETS = {
    "tfr": 1.70,
    "childless_rate": 0.15,
    "mean_age_first_birth": 26.0,
    "own_rate": 0.57547241,
    "own_family_gap": 0.16766167,
    "housing_increment_0to1": 0.66443467,
    "housing_increment_1to2": 0.56581378,
    "young_liquid_wealth_to_income": 0.60,
    "old_age_own_rate": 0.76426097,
    "old_age_parent_childless_gap": 0.070,
}


OLD_NONLOCATION_NO_TIMING_TARGETS = {
    k: v for k, v in OLD_NONLOCATION_TARGETS.items() if k != "mean_age_first_birth"
}


CANDIDATE_NO_TIMING_V0_TARGETS = {
    **OLD_NONLOCATION_NO_TIMING_TARGETS,
    "liquid_wealth_to_income": 1.20,
    "housing_user_cost_share": 0.24,
    "prime_childless_renter_median_rooms": 4.0,
    "prime_childless_owner_median_rooms": 6.0,
}


CANDIDATE_REPLACEMENT_V1_TARGETS = {
    "tfr": 1.70,
    "childless_rate": 0.15,
    "own_rate": 0.57547241,
    "own_family_gap": 0.16766167,
    "housing_increment_0to1": 0.66443467,
    "housing_increment_1to2": 0.48803143,
    "young_liquid_wealth_to_income": 0.17922556,
    "old_nonhousing_wealth_to_income_6575": 6.41854289,
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575": 1.00744952,
    "prime30_55_childless_renter_mean_rooms": 3.80528810,
    "prime30_55_childless_owner_mean_rooms": 6.22404983,
    "prime30_55_childless_renter_share_rooms_ge6": 0.13765221,
    "prime30_55_childless_owner_share_rooms_ge6": 0.59613112,
}


CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS = {
    **{k: v for k, v in CANDIDATE_REPLACEMENT_V1_TARGETS.items() if k != "old_nonhousing_wealth_to_income_6575"},
    "old_nonhousing_wealth_to_income_median_6575": 2.23046078,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS,
    "old_minus_young_owner_rate_6575_2534": 0.42309487,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_OWNGAP_V1_TARGETS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS.items()
        if k != "old_parent_childless_nonhousing_wealth_to_income_gap_6575"
    },
    "old_age_parent_childless_gap": 0.08254406,
    "old_minus_young_owner_rate_6575_2534": 0.42309487,
}


CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS,
    "old_age_own_rate": 0.76426097,
    "old_age_parent_childless_gap": 0.08254406,
}


CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_TARGETS,
    "own_rate_2534": 0.34116609,
}


CANDIDATE_REPLACEMENT_YOUNG_OLD_ROOMGAP_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_TARGETS,
    "prime30_55_childless_owner_minus_renter_mean_rooms": 2.41876173,
}


CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_TARGETS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_V1_TARGETS.items()
        if k
        not in {
            "old_nonhousing_wealth_to_income_6575",
            "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
        }
    },
    "old_total_wealth_to_income_median_6575": 5.26375360,
    "old_parent_childless_total_wealth_to_income_gap_6575": 1.28413467,
}


CANDIDATE_REPLACEMENT_TOTAL_DUE_LIFECYCLE_V1_TARGETS = {
    **CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_TARGETS,
    "old_minus_young_owner_rate_6575_2534": 0.42309487,
}


CORE_FEASIBILITY_V1_TARGETS = {
    k: CANDIDATE_NO_TIMING_V0_TARGETS[k]
    for k in [
        "tfr",
        "childless_rate",
        "own_rate",
        "own_family_gap",
        "housing_increment_0to1",
        "housing_increment_1to2",
        "prime_childless_renter_median_rooms",
        "prime_childless_owner_median_rooms",
    ]
}


COST_TEST_V1_TARGETS = {
    **CORE_FEASIBILITY_V1_TARGETS,
    "housing_user_cost_share": CANDIDATE_NO_TIMING_V0_TARGETS["housing_user_cost_share"],
}


OLDAGE_TEST_V1_TARGETS = {
    **CORE_FEASIBILITY_V1_TARGETS,
    "old_age_own_rate": CANDIDATE_NO_TIMING_V0_TARGETS["old_age_own_rate"],
    "old_age_parent_childless_gap": CANDIDATE_NO_TIMING_V0_TARGETS["old_age_parent_childless_gap"],
}


ROOMCOST_TEST_V1_TARGETS = {
    k: CANDIDATE_NO_TIMING_V0_TARGETS[k]
    for k in [
        "tfr",
        "childless_rate",
        "own_rate",
        "own_family_gap",
        "old_age_own_rate",
        "old_age_parent_childless_gap",
        "housing_user_cost_share",
        "housing_increment_0to1",
        "housing_increment_1to2",
    ]
}


CORE_WEIGHTS = {
    "own_rate": 8.0,
    "young_owner_rate": 8.0,
    "old_owner_rate": 6.0,
    "mean_completed_fertility": 10.0,
    "childless_rate": 3.0,
}


OLD_NONLOCATION_WEIGHTS = {
    "tfr": 12.0,
    "childless_rate": 12.0,
    "mean_age_first_birth": 12.0,
    "own_rate": 12.0,
    "own_family_gap": 10.0,
    "housing_increment_0to1": 8.0,
    "housing_increment_1to2": 2.0,
    "young_liquid_wealth_to_income": 12.0,
    "old_age_own_rate": 10.0,
    "old_age_parent_childless_gap": 10.0,
}


OLD_NONLOCATION_NO_TIMING_WEIGHTS = {
    k: v for k, v in OLD_NONLOCATION_WEIGHTS.items() if k != "mean_age_first_birth"
}


CANDIDATE_NO_TIMING_V0_WEIGHTS = {
    **OLD_NONLOCATION_NO_TIMING_WEIGHTS,
    "liquid_wealth_to_income": 12.0,
    "housing_user_cost_share": 250.0,
    "prime_childless_renter_median_rooms": 10.0,
    "prime_childless_owner_median_rooms": 10.0,
}


CANDIDATE_REPLACEMENT_V1_WEIGHTS = {
    "tfr": 20.0,
    "childless_rate": 20.0,
    "own_rate": 100.0,
    "own_family_gap": 45.0,
    "housing_increment_0to1": 14.0,
    "housing_increment_1to2": 8.0,
    "young_liquid_wealth_to_income": 12.0,
    "old_nonhousing_wealth_to_income_6575": 0.8,
    "old_parent_childless_nonhousing_wealth_to_income_gap_6575": 2.0,
    "prime30_55_childless_renter_mean_rooms": 6.0,
    "prime30_55_childless_owner_mean_rooms": 6.0,
    "prime30_55_childless_renter_share_rooms_ge6": 25.0,
    "prime30_55_childless_owner_share_rooms_ge6": 25.0,
}


CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS = {
    **{k: v for k, v in CANDIDATE_REPLACEMENT_V1_WEIGHTS.items() if k != "old_nonhousing_wealth_to_income_6575"},
    "old_nonhousing_wealth_to_income_median_6575": 0.8,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS,
    "old_minus_young_owner_rate_6575_2534": 80.0,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_SOFT_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS,
    "old_minus_young_owner_rate_6575_2534": 30.0,
}


CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_OWNGAP_V1_WEIGHTS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS.items()
        if k != "old_parent_childless_nonhousing_wealth_to_income_gap_6575"
    },
    "old_age_parent_childless_gap": 40.0,
    "old_minus_young_owner_rate_6575_2534": 80.0,
}


CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS,
    "old_age_own_rate": 160.0,
    "old_age_parent_childless_gap": 40.0,
}


CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_WEIGHTS,
    "own_rate_2534": 80.0,
}


CANDIDATE_REPLACEMENT_YOUNG_OLD_ROOMGAP_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_WEIGHTS,
    "prime30_55_childless_owner_minus_renter_mean_rooms": 12.0,
}


CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_WEIGHTS = {
    **{
        k: v
        for k, v in CANDIDATE_REPLACEMENT_V1_WEIGHTS.items()
        if k
        not in {
            "old_nonhousing_wealth_to_income_6575",
            "old_parent_childless_nonhousing_wealth_to_income_gap_6575",
        }
    },
    "old_total_wealth_to_income_median_6575": 0.8,
    "old_parent_childless_total_wealth_to_income_gap_6575": 2.0,
}


CANDIDATE_REPLACEMENT_TOTAL_DUE_LIFECYCLE_V1_WEIGHTS = {
    **CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_WEIGHTS,
    "old_minus_young_owner_rate_6575_2534": 80.0,
}


CORE_FEASIBILITY_V1_WEIGHTS = {
    "tfr": 20.0,
    "childless_rate": 20.0,
    "own_rate": 90.0,
    "own_family_gap": 45.0,
    "housing_increment_0to1": 14.0,
    "housing_increment_1to2": 8.0,
    "prime_childless_renter_median_rooms": 8.0,
    "prime_childless_owner_median_rooms": 8.0,
}


COST_TEST_V1_WEIGHTS = {
    **CORE_FEASIBILITY_V1_WEIGHTS,
    "housing_user_cost_share": 320.0,
}


OLDAGE_TEST_V1_WEIGHTS = {
    **CORE_FEASIBILITY_V1_WEIGHTS,
    "old_age_own_rate": 160.0,
    "old_age_parent_childless_gap": 40.0,
}


ROOMCOST_TEST_V1_WEIGHTS = {
    "tfr": 20.0,
    "childless_rate": 20.0,
    "own_rate": 100.0,
    "own_family_gap": 45.0,
    "old_age_own_rate": 160.0,
    "old_age_parent_childless_gap": 40.0,
    "housing_user_cost_share": 320.0,
    "housing_increment_0to1": 14.0,
    "housing_increment_1to2": 8.0,
}


CANDIDATE_NO_TIMING_OWNHEAVY_V1_WEIGHTS = {
    **CANDIDATE_NO_TIMING_V0_WEIGHTS,
    "own_rate": 120.0,
    "own_family_gap": 50.0,
    "old_age_own_rate": 80.0,
    "old_age_parent_childless_gap": 20.0,
    "housing_user_cost_share": 120.0,
    "liquid_wealth_to_income": 8.0,
    "young_liquid_wealth_to_income": 8.0,
    "housing_increment_0to1": 6.0,
    "prime_childless_renter_median_rooms": 5.0,
    "prime_childless_owner_median_rooms": 5.0,
}


CANDIDATE_NO_TIMING_REFINEMENT_V1_WEIGHTS = {
    **CANDIDATE_NO_TIMING_V0_WEIGHTS,
    "tfr": 20.0,
    "childless_rate": 20.0,
    "own_rate": 110.0,
    "own_family_gap": 45.0,
    "old_age_own_rate": 160.0,
    "old_age_parent_childless_gap": 40.0,
    "housing_user_cost_share": 320.0,
    "liquid_wealth_to_income": 20.0,
    "young_liquid_wealth_to_income": 14.0,
    "housing_increment_0to1": 14.0,
    "housing_increment_1to2": 8.0,
    "prime_childless_renter_median_rooms": 6.0,
    "prime_childless_owner_median_rooms": 6.0,
}


TARGET_SETS = {
    "core": (CORE_TARGETS, CORE_WEIGHTS),
    "old_nonlocation": (OLD_NONLOCATION_TARGETS, OLD_NONLOCATION_WEIGHTS),
    "old_nonlocation_no_timing": (OLD_NONLOCATION_NO_TIMING_TARGETS, OLD_NONLOCATION_NO_TIMING_WEIGHTS),
    "candidate_no_timing_v0": (CANDIDATE_NO_TIMING_V0_TARGETS, CANDIDATE_NO_TIMING_V0_WEIGHTS),
    "candidate_replacement_v1": (CANDIDATE_REPLACEMENT_V1_TARGETS, CANDIDATE_REPLACEMENT_V1_WEIGHTS),
    "candidate_replacement_nh_median_v1": (
        CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_TARGETS,
        CANDIDATE_REPLACEMENT_NH_MEDIAN_V1_WEIGHTS,
    ),
    "candidate_replacement_due_lifecycle_v1": (
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_TARGETS,
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_WEIGHTS,
    ),
    "candidate_replacement_due_lifecycle_soft_v1": (
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_V1_TARGETS,
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_SOFT_V1_WEIGHTS,
    ),
    "candidate_replacement_due_lifecycle_owngap_v1": (
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_OWNGAP_V1_TARGETS,
        CANDIDATE_REPLACEMENT_DUE_LIFECYCLE_OWNGAP_V1_WEIGHTS,
    ),
    "candidate_replacement_old_retention_v1": (
        CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_TARGETS,
        CANDIDATE_REPLACEMENT_OLD_RETENTION_V1_WEIGHTS,
    ),
    "candidate_replacement_young_old_own_v1": (
        CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_TARGETS,
        CANDIDATE_REPLACEMENT_YOUNG_OLD_OWN_V1_WEIGHTS,
    ),
    "candidate_replacement_young_old_roomgap_v1": (
        CANDIDATE_REPLACEMENT_YOUNG_OLD_ROOMGAP_V1_TARGETS,
        CANDIDATE_REPLACEMENT_YOUNG_OLD_ROOMGAP_V1_WEIGHTS,
    ),
    "candidate_replacement_total_median_v1": (
        CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_TARGETS,
        CANDIDATE_REPLACEMENT_TOTAL_MEDIAN_V1_WEIGHTS,
    ),
    "candidate_replacement_total_due_lifecycle_v1": (
        CANDIDATE_REPLACEMENT_TOTAL_DUE_LIFECYCLE_V1_TARGETS,
        CANDIDATE_REPLACEMENT_TOTAL_DUE_LIFECYCLE_V1_WEIGHTS,
    ),
    "candidate_no_timing_core_feasibility_v1": (
        CORE_FEASIBILITY_V1_TARGETS,
        CORE_FEASIBILITY_V1_WEIGHTS,
    ),
    "candidate_no_timing_cost_test_v1": (
        COST_TEST_V1_TARGETS,
        COST_TEST_V1_WEIGHTS,
    ),
    "candidate_no_timing_oldage_test_v1": (
        OLDAGE_TEST_V1_TARGETS,
        OLDAGE_TEST_V1_WEIGHTS,
    ),
    "candidate_no_timing_roomcost_test_v1": (
        ROOMCOST_TEST_V1_TARGETS,
        ROOMCOST_TEST_V1_WEIGHTS,
    ),
    "candidate_no_timing_ownheavy_v1": (
        CANDIDATE_NO_TIMING_V0_TARGETS,
        CANDIDATE_NO_TIMING_OWNHEAVY_V1_WEIGHTS,
    ),
    "candidate_no_timing_refinement_v1": (
        CANDIDATE_NO_TIMING_V0_TARGETS,
        CANDIDATE_NO_TIMING_REFINEMENT_V1_WEIGHTS,
    ),
}


def run_small_calibration(
    outdir: Path,
    *,
    n_cases: int = 24,
    seed: int = 1234,
    J: int = 12,
    Nb: int = 40,
    n_house: int = 4,
    max_iter_eq: int = 35,
    progress: bool = True,
    target_set: str = "old_nonlocation",
) -> dict[str, Any]:
    """Run a small random-search diagnostic calibration.

    This is intentionally not a formal SMM routine. Its purpose is to find a
    locally more interpretable parameter point for inspecting the scaffold.
    """

    rng = np.random.default_rng(seed)
    targets, weights = get_target_set(target_set)
    outdir.mkdir(parents=True, exist_ok=True)
    cases_path = outdir / "cases.jsonl"
    best_path = outdir / "best.json"
    meta = {
        "target_set": str(target_set),
        "targets": targets,
        "weights": weights,
        "n_cases": int(n_cases),
        "seed": int(seed),
        "J": int(J),
        "period_years": PERIOD_YEARS,
        "age_start": AGE_START,
        "fertility_start_age": FERTILITY_START_AGE,
        "fertility_end_age": FERTILITY_END_AGE,
        "retirement_age": RETIREMENT_AGE,
        "child_maturity_age": CHILD_MATURITY_AGE,
        "stage_durations": [CHILD_MATURITY_AGE / PERIOD_YEARS],
        "expected_child_years": [CHILD_MATURITY_AGE],
        "Nb": int(Nb),
        "n_house": int(n_house),
        "max_iter_eq": int(max_iter_eq),
        "progress": bool(progress),
        "status": "diagnostic_only_not_formal_smm",
    }
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True))

    base = base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq)
    best: dict[str, Any] | None = None
    start = time.perf_counter()
    cases_path.write_text("")
    for idx in range(n_cases):
        theta = draw_candidate(rng, idx)
        overrides = {**base, **theta}
        t0 = time.perf_counter()
        try:
            sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
            moments = extract_moments(sol, P)
            loss = diagnostic_loss(moments, targets=targets, weights=weights)
            status = "ok"
            err = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
            timings = getattr(sol, "timings", {})
        except Exception as exc:  # noqa: BLE001 - calibration should checkpoint failures.
            moments = {}
            loss = math.inf
            status = f"failed: {type(exc).__name__}: {exc}"
            p_eq = np.array([np.nan])
            err = math.inf
            timings = {}
        elapsed = time.perf_counter() - t0
        record = {
            "case": int(idx),
            "status": status,
            "loss": float(loss),
            "theta": jsonable(theta),
            "moments": jsonable(moments),
            "p_eq": jsonable(p_eq),
            "market_residual": float(err),
            "elapsed_sec": float(elapsed),
            "timings": jsonable(timings),
        }
        with cases_path.open("a") as fh:
            fh.write(json.dumps(record, sort_keys=True) + "\n")
        if best is None or record["loss"] < best["loss"]:
            best = record
            best_path.write_text(json.dumps(best, indent=2, sort_keys=True))
        if progress:
            best_loss = float(best["loss"]) if best is not None else math.inf
            print(
                f"case {idx + 1}/{n_cases}: loss={record['loss']:.4g}, "
                f"resid={record['market_residual']:.2e}, best={best_loss:.4g}, "
                f"elapsed={elapsed:.1f}s",
                flush=True,
            )

    summary = {
        "best": best,
        "elapsed_sec": float(time.perf_counter() - start),
        "metadata": meta,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def run_informed_smoke(
    outdir: Path,
    *,
    J: int = 16,
    Nb: int = 40,
    n_house: int = 6,
    max_iter_eq: int = 25,
    labels: list[str] | None = None,
    case_limit: int | None = None,
    fixed_price: float | None = None,
    progress: bool = True,
) -> dict[str, Any]:
    """Run a deterministic production-parameter smoke panel.

    This is a parameter-led smoke check, not a formal calibration. It varies
    economically meaningful internal blocks around the current scaffold and
    keeps finance, supply, and menu objects fixed at their benchmark inputs.
    """

    targets, weights = get_target_set("old_nonlocation")
    rank_targets = {k: v for k, v in targets.items() if k != "mean_age_first_birth"}
    rank_weights = {k: v for k, v in weights.items() if k != "mean_age_first_birth"}
    outdir.mkdir(parents=True, exist_ok=True)
    cases_path = outdir / "cases.jsonl"
    best_path = outdir / "best.json"
    base = base_overrides(J=J, Nb=Nb, n_house=n_house, max_iter_eq=max_iter_eq)
    if fixed_price is not None:
        base.update(
            {
                "solve_mode": "pe",
                "p_fixed": np.array([float(fixed_price)]),
                "w_fixed": np.array([1.0]),
                "entry_shares_fixed": np.array([1.0]),
            }
        )
    candidates = informed_smoke_candidates()
    if labels:
        wanted = {str(label) for label in labels}
        candidates = [candidate for candidate in candidates if str(candidate["label"]) in wanted]
        missing = sorted(wanted - {str(candidate["label"]) for candidate in candidates})
        if missing:
            raise ValueError(f"Unknown informed-smoke labels: {missing}")
    if case_limit is not None:
        candidates = candidates[: max(0, int(case_limit))]
    meta = {
        "status": "deterministic_informed_smoke_not_formal_calibration",
        "case_count": len(candidates),
        "labels": [str(candidate["label"]) for candidate in candidates],
        "fixed_price": None if fixed_price is None else float(fixed_price),
        "J": int(J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "max_iter_eq": int(max_iter_eq),
        "fixed_external_or_first_stage": [
            "q",
            "delta",
            "tau_H",
            "phi",
            "pti_limit",
            "psi",
            "income_age_profile",
            "z_process",
            "H_own",
            "hR_max",
            "H0",
            "r_bar",
            "eta_supply",
            "theta1",
        ],
        "varied_blocks": [
            "beta",
            "alpha_cons",
            "b_entry_fixed",
            "c_bar_0",
            "c_bar_n",
            "h_bar_0",
            "h_bar_jump",
            "h_bar_n",
            "psi_child",
            "kappa_fert",
            "chi",
            "theta0",
            "theta_n",
        ],
        "missing_internal_parameter": "lambda_a_first_birth_timing_shifter",
        "ranking_loss": "old_nonlocation_excluding_mean_age_first_birth",
        "targets": targets,
        "weights": weights,
        "ranking_targets": rank_targets,
        "ranking_weights": rank_weights,
    }
    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True))
    cases_path.write_text("")
    best: dict[str, Any] | None = None
    start = time.perf_counter()
    for idx, candidate in enumerate(candidates):
        label = str(candidate["label"])
        theta = dict(candidate["theta"])
        overrides = {**base, **theta}
        t0 = time.perf_counter()
        try:
            sol, P, p_eq = run_model_cp_dt(overrides, verbose=False)
            moments = extract_moments(sol, P)
            full_loss = diagnostic_loss(moments, targets=targets, weights=weights)
            rank_loss = diagnostic_loss(moments, targets=rank_targets, weights=rank_weights)
            status = "ok"
            err = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
            timings = getattr(sol, "timings", {})
        except Exception as exc:  # noqa: BLE001 - smoke panel should checkpoint failures.
            moments = {}
            full_loss = math.inf
            rank_loss = math.inf
            status = f"failed: {type(exc).__name__}: {exc}"
            p_eq = np.array([np.nan])
            err = math.inf
            timings = {}
        elapsed = time.perf_counter() - t0
        record = {
            "case": int(idx),
            "label": label,
            "status": status,
            "rank_loss": float(rank_loss),
            "full_old_nonlocation_loss": float(full_loss),
            "theta": jsonable(theta),
            "moments": jsonable(moments),
            "p_eq": jsonable(p_eq),
            "market_residual": float(err),
            "elapsed_sec": float(elapsed),
            "timings": jsonable(timings),
        }
        with cases_path.open("a") as fh:
            fh.write(json.dumps(record, sort_keys=True) + "\n")
        if best is None or record["rank_loss"] < best["rank_loss"]:
            best = record
            best_path.write_text(json.dumps(best, indent=2, sort_keys=True))
        if progress:
            print(
                f"case {idx + 1}/{len(candidates)} {label}: "
                f"rank={record['rank_loss']:.4g}, full={record['full_old_nonlocation_loss']:.4g}, "
                f"resid={record['market_residual']:.2e}, elapsed={elapsed:.1f}s",
                flush=True,
            )

    summary = {
        "best": best,
        "elapsed_sec": float(time.perf_counter() - start),
        "metadata": meta,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def informed_smoke_candidates() -> list[dict[str, Any]]:
    beta_hi = 0.98**PERIOD_YEARS
    beta_lo = 0.94**PERIOD_YEARS
    return [
        {"label": "baseline", "theta": {}},
        {"label": "patient_beta", "theta": {"beta": beta_hi}},
        {"label": "impatient_beta", "theta": {"beta": beta_lo}},
        {
            "label": "more_housing_share",
            "theta": {"alpha_cons": 0.60, "h_bar_0": 3.5, "c_bar_0": 0.08 * PERIOD_YEARS},
        },
        {
            "label": "less_housing_share",
            "theta": {"alpha_cons": 0.78, "h_bar_0": 4.5, "c_bar_0": 0.12 * PERIOD_YEARS},
        },
        {
            "label": "higher_family_space_need",
            "theta": {"h_bar_jump": 1.10, "h_bar_n": 0.90, "h_bar_0": 3.6},
        },
        {
            "label": "lower_child_cost_higher_taste",
            "theta": {"c_bar_n": 0.32, "psi_child": 0.10, "kappa_fert": 5.5},
        },
        {
            "label": "higher_child_cost_lower_taste",
            "theta": {"c_bar_n": 0.72, "psi_child": 0.045, "kappa_fert": 4.0},
        },
        {"label": "higher_owner_premium", "theta": {"chi": 1.30}},
        {
            "label": "stronger_bequest_parent_tilt",
            "theta": {"theta0": 0.80, "theta_n": 0.45},
        },
        {
            "label": "combined_internal_push",
            "theta": {
                "beta": beta_hi,
                "alpha_cons": 0.62,
                "b_entry_fixed": 2.0,
                "c_bar_0": 0.08 * PERIOD_YEARS,
                "c_bar_n": 0.40,
                "h_bar_0": 3.5,
                "h_bar_jump": 1.00,
                "h_bar_n": 0.80,
                "psi_child": 0.09,
                "kappa_fert": 5.5,
                "chi": 1.25,
                "theta0": 0.75,
                "theta_n": 0.40,
            },
        },
    ]


def base_overrides(*, J: int, Nb: int, n_house: int, max_iter_eq: int) -> dict[str, Any]:
    lifecycle = lifecycle_overrides(J)
    return {
        "I": 1,
        "w_hat": np.array([1.0]),
        "E_loc": np.array([0.0]),
        "N_0": np.array([1.0]),
        "entry_shares": np.array([1.0]),
        "entry_by_loc": np.array([1.0 / J]),
        "r_bar": np.array([0.16]),
        "H0": np.array([4.0]),
        "eta_supply": np.array([1.0]),
        **lifecycle,
        "J": int(J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "H_own": np.linspace(2.0, 10.0, int(n_house)),
        "max_iter_eq": int(max_iter_eq),
        "tol_eq": 1e-4,
        "use_pti_constraint": False,
        "scalar_market_refine": True,
        "scalar_market_refine_method": "brent",
        "scalar_market_refine_iter": 16,
        "scalar_market_refine_max_expand": 8,
    }


def lifecycle_overrides(J: int) -> dict[str, Any]:
    j = int(J)
    period_years = PERIOD_YEARS
    fertile_start = age_to_period_number(FERTILITY_START_AGE, period_years)
    fertile_end = age_to_period_number(FERTILITY_END_AGE, period_years)
    retirement_idx = int(round((RETIREMENT_AGE - AGE_START) / period_years))
    return {
        "period_years": period_years,
        "da": period_years,
        "age_start": AGE_START,
        "A_m": CHILD_MATURITY_AGE,
        "stage_durations": np.array([CHILD_MATURITY_AGE / period_years]),
        "J_R": max(2, min(j - 1, retirement_idx)),
        "A_f_start": max(1, min(j, fertile_start)),
        "A_f_end": max(1, min(j, fertile_end)),
    }


def age_to_period_number(age: float, period_years: float) -> int:
    return int(round((float(age) - AGE_START) / period_years)) + 1


def draw_candidate(rng: np.random.Generator, idx: int) -> dict[str, Any]:
    if idx == 0:
        return {}
    if idx == 1:
        return {"hR_max": 4.0}
    if idx == 2:
        return {"hR_max": 5.0}
    if idx == 3:
        return {"b_entry_fixed": 5.0, "phi": np.array([0.95, 0.95, 0.95])}
    if idx == 4:
        return {"hR_max": 4.0, "H0": np.array([3.0])}
    phi = rng.uniform(0.82, 0.97)
    return {
        "alpha_cons": rng.uniform(0.62, 0.82),
        "phi": np.array([phi, phi, phi]),
        "b_entry_fixed": rng.uniform(0.0, 7.0),
        "c_bar_n": rng.uniform(0.24, 0.72),
        "chi": rng.uniform(0.75, 1.50),
        "kappa_fert": rng.uniform(3.0, 7.0),
        "psi": rng.uniform(0.00, 0.08),
        "psi_child": rng.uniform(0.03, 0.12),
        "theta0": rng.uniform(0.12, 0.85),
        "theta_n": rng.uniform(0.00, 0.50),
        "theta1": rng.uniform(-0.05, 0.08),
        "h_bar_jump": rng.uniform(0.35, 1.10),
        "h_bar_n": rng.uniform(0.25, 0.90),
        "hR_max": rng.uniform(3.6, 6.5),
        "owner_h_bar_scale": rng.uniform(0.75, 1.15),
        "owner_size_cost": rng.uniform(0.0, 0.05),
        "owner_size_cost_ref": rng.uniform(2.0, 8.0),
        "H0": np.array([rng.uniform(3.0, 5.8)]),
    }


def extract_moments(sol: Any, P: Any | None = None) -> dict[str, float]:
    household_parity = float(getattr(sol, "mean_completed_fertility", np.nan))
    parity = np.asarray(getattr(sol, "parity_dist", np.array([np.nan])), dtype=float).reshape(-1)
    renter_rooms = float(getattr(sol, "prime_childless_renter_median_rooms", np.nan))
    owner_rooms = float(getattr(sol, "prime_childless_owner_median_rooms", np.nan))
    own_rate_2534 = float(getattr(sol, "own_rate_2534", np.nan))
    old_age_own_rate = float(getattr(sol, "old_age_own_rate_6575", np.nan))
    return {
        "tfr": 2.0 * household_parity,
        "own_rate": float(getattr(sol, "own_rate_3055", np.nan)),
        "aggregate_own_rate": float(getattr(sol, "own_rate", np.nan)),
        "young_owner_rate": float(getattr(sol, "young_owner_rate", np.nan)),
        "old_owner_rate": float(getattr(sol, "old_owner_rate", np.nan)),
        "old_age_own_rate": old_age_own_rate,
        "mean_completed_fertility": household_parity,
        "childless_rate": float(getattr(sol, "childless_rate", np.nan)),
        "own_rate_3055": float(getattr(sol, "own_rate_3055", np.nan)),
        "own_rate_2534": own_rate_2534,
        "old_minus_young_owner_rate_6575_2534": old_age_own_rate - own_rate_2534,
        "own_rate_3544": float(getattr(sol, "own_rate_3544", np.nan)),
        "own_rate_nonparents_3055": float(getattr(sol, "own_rate_nonparents_3055", np.nan)),
        "own_rate_newparents_3055": float(getattr(sol, "own_rate_newparents_3055", np.nan)),
        "own_gap_newparent_nonparent_3055": float(getattr(sol, "own_gap_newparent_nonparent_3055", np.nan)),
        "own_family_gap": float(
            getattr(sol, "own_gap_newparent_nonparent_3055", getattr(sol, "own_family_gap", np.nan))
        ),
        "old_age_parent_childless_gap_6575": float(getattr(sol, "old_age_parent_childless_gap_6575", np.nan)),
        "old_age_parent_childless_gap": float(getattr(sol, "old_age_parent_childless_gap_6575", np.nan)),
        "old_nonhousing_wealth_to_income_6575": float(getattr(sol, "old_nonhousing_wealth_to_income_6575", np.nan)),
        "old_total_wealth_to_income_6575": float(getattr(sol, "old_total_wealth_to_income_6575", np.nan)),
        "old_nonhousing_wealth_to_income_median_6575": float(
            getattr(sol, "old_nonhousing_wealth_to_income_median_6575", np.nan)
        ),
        "old_total_wealth_to_income_median_6575": float(
            getattr(sol, "old_total_wealth_to_income_median_6575", np.nan)
        ),
        "old_parent_childless_nonhousing_wealth_to_income_gap_6575": float(
            getattr(sol, "old_parent_childless_nonhousing_wealth_to_income_gap_6575", np.nan)
        ),
        "old_parent_childless_total_wealth_to_income_gap_6575": float(
            getattr(sol, "old_parent_childless_total_wealth_to_income_gap_6575", np.nan)
        ),
        "old_parent_childless_nonhousing_wealth_to_income_median_gap_6575": float(
            getattr(sol, "old_parent_childless_nonhousing_wealth_to_income_median_gap_6575", np.nan)
        ),
        "old_parent_childless_total_wealth_to_income_median_gap_6575": float(
            getattr(sol, "old_parent_childless_total_wealth_to_income_median_gap_6575", np.nan)
        ),
        "mean_age_first_birth": float(getattr(sol, "mean_age_first_birth", np.nan)),
        "housing_increment_0to1": float(getattr(sol, "housing_increment_0to1_eventstudy_t3", np.nan)),
        "housing_increment_1to2": float(
            getattr(sol, "housing_increment_1to2_proxy_t3", getattr(sol, "housing_increment_1to2", np.nan))
        ),
        "young_liquid_wealth_to_income": float(getattr(sol, "young_liquid_wealth_to_income", np.nan)),
        "liquid_wealth_to_income": float(getattr(sol, "liquid_wealth_to_income", np.nan)),
        "wealth_to_income": float(getattr(sol, "wealth_to_income", np.nan)),
        "parity_share_0": float(parity[0]) if parity.size > 0 else np.nan,
        "parity_share_1": float(parity[1]) if parity.size > 1 else np.nan,
        "parity_share_2plus": float(np.sum(parity[2:])) if parity.size > 2 else np.nan,
        "prime_childless_renter_median_rooms": renter_rooms,
        "prime_childless_owner_median_rooms": owner_rooms,
        "prime_childless_owner_minus_renter_rooms": owner_rooms - renter_rooms,
        "prime30_55_childless_renter_mean_rooms": float(
            getattr(sol, "prime30_55_childless_renter_mean_rooms", np.nan)
        ),
        "prime30_55_childless_owner_mean_rooms": float(
            getattr(sol, "prime30_55_childless_owner_mean_rooms", np.nan)
        ),
        "prime30_55_childless_owner_minus_renter_mean_rooms": float(
            getattr(sol, "prime30_55_childless_owner_minus_renter_mean_rooms", np.nan)
        ),
        "prime30_55_childless_renter_share_rooms_ge6": float(
            getattr(sol, "prime30_55_childless_renter_share_rooms_ge6", np.nan)
        ),
        "prime30_55_childless_owner_share_rooms_ge6": float(
            getattr(sol, "prime30_55_childless_owner_share_rooms_ge6", np.nan)
        ),
        "prime30_55_parent_owner_minus_renter_mean_rooms": float(
            getattr(sol, "prime30_55_parent_owner_minus_renter_mean_rooms", np.nan)
        ),
        "housing_user_cost_share": housing_user_cost_share(sol, P),
        "renter25_45_all_cap_share": float(getattr(sol, "renter25_45_all_cap_share", np.nan)),
        "owner_neg_liquid_share_2534": float(getattr(sol, "owner_neg_liquid_share_2534", np.nan)),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
    }


def housing_user_cost_share(sol: Any, P: Any | None) -> float:
    if P is None or not hasattr(sol, "g") or not hasattr(sol, "hR_pol"):
        return np.nan
    g = np.asarray(sol.g, dtype=float)
    hR = np.asarray(sol.hR_pol, dtype=float)
    user_cost = np.asarray(getattr(sol, "owner_user_cost", np.nan), dtype=float).reshape(-1)
    if g.ndim != 7 or hR.shape != g.shape or user_cost.size != int(P.I):
        return np.nan
    z_grid = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    total_income = 0.0
    total_housing = 0.0
    for i in range(P.I):
        for j in range(P.J):
            for zz, z_value in enumerate(z_grid):
                yj = float(P.income[i, j]) * (float(z_value) if j < int(P.J_R) else 1.0)
                mass_ijz = float(np.sum(g[:, :, i, j, zz, :, :]))
                total_income += yj * mass_ijz
                total_housing += float(np.sum(g[:, 0, i, j, zz, :, :] * hR[:, 0, i, j, zz, :, :])) * user_cost[i]
                for ten in range(1, 1 + int(P.n_house)):
                    mass = float(np.sum(g[:, ten, i, j, zz, :, :]))
                    total_housing += mass * user_cost[i] * float(P.H_own[ten - 1])
    return float(total_housing / max(total_income, 1e-12))


def get_target_set(target_set: str) -> tuple[dict[str, float], dict[str, float]]:
    key = str(target_set).lower()
    if key not in TARGET_SETS:
        raise ValueError(f"Unknown target_set {target_set!r}; choose from {sorted(TARGET_SETS)}")
    targets, weights = TARGET_SETS[key]
    return dict(targets), dict(weights)


def diagnostic_loss(
    moments: dict[str, float],
    *,
    targets: dict[str, float] | None = None,
    weights: dict[str, float] | None = None,
) -> float:
    if targets is None or weights is None:
        targets, weights = get_target_set("core")
    loss = 0.0
    for name, target in targets.items():
        value = float(moments.get(name, np.nan))
        if not np.isfinite(value):
            return math.inf
        loss += float(weights.get(name, 1.0)) * (value - target) ** 2
    residual = float(moments.get("market_residual", np.nan))
    if not np.isfinite(residual) or residual > 5e-3:
        loss += 100.0
    return float(loss)


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value
