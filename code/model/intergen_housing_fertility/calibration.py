"""Small diagnostic calibration driver for the one-market model scaffold."""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np

from .solver import run_model_cp_dt


DIAGNOSTIC_TARGETS = {
    "own_rate": 0.575,
    "young_owner_rate": 0.25,
    "old_owner_rate": 0.764,
    "mean_completed_fertility": 0.95,
    "childless_rate": 0.20,
}


DIAGNOSTIC_WEIGHTS = {
    "own_rate": 8.0,
    "young_owner_rate": 8.0,
    "old_owner_rate": 6.0,
    "mean_completed_fertility": 10.0,
    "childless_rate": 3.0,
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
) -> dict[str, Any]:
    """Run a small random-search diagnostic calibration.

    This is intentionally not a formal SMM routine. Its purpose is to find a
    locally more interpretable parameter point for inspecting the scaffold.
    """

    rng = np.random.default_rng(seed)
    outdir.mkdir(parents=True, exist_ok=True)
    cases_path = outdir / "cases.jsonl"
    best_path = outdir / "best.json"
    meta = {
        "targets": DIAGNOSTIC_TARGETS,
        "weights": DIAGNOSTIC_WEIGHTS,
        "n_cases": int(n_cases),
        "seed": int(seed),
        "J": int(J),
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
            loss = diagnostic_loss(moments)
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
        "stage_durations": np.array([1.0]),
        "J": int(J),
        "J_R": max(2, min(J - 2, 8)),
        "A_f_start": 2,
        "A_f_end": min(5, J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "H_own": np.linspace(2.0, 10.0, int(n_house)),
        "max_iter_eq": int(max_iter_eq),
        "tol_eq": 1e-4,
        "use_pti_constraint": True,
        "scalar_market_refine": True,
        "scalar_market_refine_iter": 16,
    }


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
        "phi": np.array([phi, phi, phi]),
        "b_entry_fixed": rng.uniform(0.0, 7.0),
        "pti_limit": rng.uniform(0.28, 0.65),
        "chi": rng.uniform(0.75, 1.50),
        "kappa_fert": rng.uniform(3.0, 7.0),
        "theta0": rng.uniform(0.12, 0.85),
        "h_bar_jump": rng.uniform(0.35, 1.10),
        "h_bar_n": rng.uniform(0.25, 0.90),
        "hR_max": rng.uniform(3.6, 6.5),
        "H0": np.array([rng.uniform(3.0, 5.8)]),
    }


def extract_moments(sol: Any) -> dict[str, float]:
    return {
        "own_rate": float(getattr(sol, "own_rate", np.nan)),
        "young_owner_rate": float(getattr(sol, "young_owner_rate", np.nan)),
        "old_owner_rate": float(getattr(sol, "old_owner_rate", np.nan)),
        "mean_completed_fertility": float(getattr(sol, "mean_completed_fertility", np.nan)),
        "childless_rate": float(getattr(sol, "childless_rate", np.nan)),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
    }


def diagnostic_loss(moments: dict[str, float]) -> float:
    loss = 0.0
    for name, target in DIAGNOSTIC_TARGETS.items():
        value = float(moments.get(name, np.nan))
        if not np.isfinite(value):
            return math.inf
        loss += DIAGNOSTIC_WEIGHTS[name] * (value - target) ** 2
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
