#!/usr/bin/env python3
"""One-sided finite-difference moment Jacobian at the current-bound candidate.

Tight-equilibrium solves (max_iter_eq=40, tol_eq=2.5e-5) to keep FD signal above
the equilibrium-acceptance noise measured in the repeatability battery.
Steps are 2% of each parameter's box width, flipped inward at bounds.
Audit-only; writes JSON/CSV under the audit folder.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

import intergen_housing_fertility.local_panel as lp
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_SEARCH_BOUNDS,
    PRODUCTION_TARGET_SET,
    production_profile_overrides,
)

MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
PERIOD_YEARS = 4.0
NB = 120
MAX_ITER_EQ = 40
TOL_EQ = 2.5e-5
STEP_FRAC = 0.02

AUDIT = Path(__file__).resolve().parents[1]
REPORT = AUDIT.parent / "combined_recalibration/overnight_20260710_report"


def build_env():
    targets, weights = lp.get_target_set(PRODUCTION_TARGET_SET)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = 6.0
    income = lp.income_process_overrides(5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO)
    extra = production_profile_overrides()
    extra.update({
        "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
        "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
        "eta_supply": np.array([1.75]),
        "normalize_bequest_utility": True,
        "max_iter_eq": MAX_ITER_EQ,
        "tol_eq": TOL_EQ,
    })
    return targets, weights, income, extra


def solve(theta, label, income, targets, weights, extra):
    rec = lp.run_local_panel_case(0, {"label": label, "theta": theta}, PRODUCTION_J, NB, 5,
                                  MAX_ITER_EQ, income, targets, weights, extra)
    if rec["status"] != "ok":
        raise RuntimeError(f"{label}: {rec['status']}")
    return rec


def main() -> None:
    payload = json.loads((REPORT / "current_bound_best.json").read_text())
    theta0 = dict(payload["theta"])

    bounds = [*[tuple(b) for b in PRODUCTION_SEARCH_BOUNDS], ("H0", 1.0, 10.0)]
    targets, weights, income, extra = build_env()
    moment_names = sorted(targets)

    t0 = time.perf_counter()
    base = solve(theta0, "jac_base", income, targets, weights, extra)
    m0 = np.array([float(base["moments"][k]) for k in moment_names])
    print(f"base loss={base['rank_loss']:.6f} residual={base['market_residual']:.2e} "
          f"({time.perf_counter()-t0:.0f}s)", flush=True)

    rows = []
    Jmat = np.zeros((len(moment_names), len(bounds)))
    for pi, (name, lo, hi) in enumerate(bounds):
        width = hi - lo
        key = "beta" if name == "beta_annual" else name
        cur_search = float(theta0[key]) ** (1.0 / PERIOD_YEARS) if name == "beta_annual" else float(theta0[key])
        step = STEP_FRAC * width
        if cur_search + step > hi:
            step = -step
        new_search = cur_search + step
        theta = dict(theta0)
        theta[key] = new_search ** PERIOD_YEARS if name == "beta_annual" else new_search
        t1 = time.perf_counter()
        rec = solve(theta, f"jac_{name}", income, targets, weights, extra)
        m1 = np.array([float(rec["moments"][k]) for k in moment_names])
        Jmat[:, pi] = (m1 - m0) / (step / width)  # derivative w.r.t. unit-cube coordinate
        rows.append({
            "parameter": name, "search_value": cur_search, "step": step,
            "loss": rec["rank_loss"], "residual": rec["market_residual"],
            "strict": rec["strict_converged"],
        })
        print(f"{name}: dloss={rec['rank_loss']-base['rank_loss']:+.4f} step={step:+.4f} "
              f"res={rec['market_residual']:.1e} ({time.perf_counter()-t1:.0f}s)", flush=True)

    w = np.array([float(weights[k]) for k in moment_names])
    WJ = np.sqrt(w)[:, None] * Jmat
    sv = np.linalg.svd(WJ, compute_uv=False)
    out = {
        "moment_names": moment_names,
        "parameters": [b[0] for b in bounds],
        "jacobian_unit_cube": Jmat.tolist(),
        "weighted_singular_values": sv.tolist(),
        "condition_number": float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf"),
        "rank_tol_relative_1e6": int(np.sum(sv > sv[0] * 1e-6)),
        "rank_tol_relative_1e3": int(np.sum(sv > sv[0] * 1e-3)),
        "rank_tol_relative_1e2": int(np.sum(sv > sv[0] * 1e-2)),
        "base": {"loss": base["rank_loss"], "residual": base["market_residual"],
                 "moments": {k: float(base["moments"][k]) for k in moment_names}},
        "fd_rows": rows,
        "config": {"Nb": NB, "max_iter_eq": MAX_ITER_EQ, "tol_eq": TOL_EQ, "step_frac": STEP_FRAC},
    }
    outpath = AUDIT / "diagnostics/jacobian_current_bound.json"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text(json.dumps(out, indent=2))
    print("singular values:", np.array2string(sv, precision=4))
    print(f"condition number: {out['condition_number']:.3e}; "
          f"rank@1e-3: {out['rank_tol_relative_1e3']}/{len(bounds)}")


if __name__ == "__main__":
    main()
