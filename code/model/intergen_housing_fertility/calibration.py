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


TARGET_SETS = {
    "core": (CORE_TARGETS, CORE_WEIGHTS),
    "old_nonlocation": (OLD_NONLOCATION_TARGETS, OLD_NONLOCATION_WEIGHTS),
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
            moments = extract_moments(sol)
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
        "use_pti_constraint": True,
        "scalar_market_refine": True,
        "scalar_market_refine_iter": 16,
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
        return {"b_entry_fixed": 5.0, "phi": np.array([0.95, 0.95, 0.95]), "pti_limit": 0.45}
    if idx == 4:
        return {"hR_max": 4.0, "H0": np.array([3.0])}
    phi = rng.uniform(0.82, 0.97)
    return {
        "alpha_cons": rng.uniform(0.62, 0.82),
        "phi": np.array([phi, phi, phi]),
        "b_entry_fixed": rng.uniform(0.0, 7.0),
        "pti_limit": rng.uniform(0.28, 0.65),
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


def extract_moments(sol: Any) -> dict[str, float]:
    household_parity = float(getattr(sol, "mean_completed_fertility", np.nan))
    return {
        "tfr": 2.0 * household_parity,
        "own_rate": float(getattr(sol, "own_rate_3055", np.nan)),
        "aggregate_own_rate": float(getattr(sol, "own_rate", np.nan)),
        "young_owner_rate": float(getattr(sol, "young_owner_rate", np.nan)),
        "old_owner_rate": float(getattr(sol, "old_owner_rate", np.nan)),
        "old_age_own_rate": float(getattr(sol, "old_age_own_rate_6575", np.nan)),
        "mean_completed_fertility": household_parity,
        "childless_rate": float(getattr(sol, "childless_rate", np.nan)),
        "own_rate_3055": float(getattr(sol, "own_rate_3055", np.nan)),
        "own_rate_2534": float(getattr(sol, "own_rate_2534", np.nan)),
        "own_rate_3544": float(getattr(sol, "own_rate_3544", np.nan)),
        "own_rate_nonparents_3055": float(getattr(sol, "own_rate_nonparents_3055", np.nan)),
        "own_rate_newparents_3055": float(getattr(sol, "own_rate_newparents_3055", np.nan)),
        "own_gap_newparent_nonparent_3055": float(getattr(sol, "own_gap_newparent_nonparent_3055", np.nan)),
        "own_family_gap": float(
            getattr(sol, "own_gap_newparent_nonparent_3055", getattr(sol, "own_family_gap", np.nan))
        ),
        "old_age_parent_childless_gap_6575": float(getattr(sol, "old_age_parent_childless_gap_6575", np.nan)),
        "old_age_parent_childless_gap": float(getattr(sol, "old_age_parent_childless_gap_6575", np.nan)),
        "mean_age_first_birth": float(getattr(sol, "mean_age_first_birth", np.nan)),
        "housing_increment_0to1": float(getattr(sol, "housing_increment_0to1_eventstudy_t3", np.nan)),
        "housing_increment_1to2": float(
            getattr(sol, "housing_increment_1to2_proxy_t3", getattr(sol, "housing_increment_1to2", np.nan))
        ),
        "young_liquid_wealth_to_income": float(getattr(sol, "young_liquid_wealth_to_income", np.nan)),
        "renter25_45_all_cap_share": float(getattr(sol, "renter25_45_all_cap_share", np.nan)),
        "owner_neg_liquid_share_2534": float(getattr(sol, "owner_neg_liquid_share_2534", np.nan)),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
    }


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
