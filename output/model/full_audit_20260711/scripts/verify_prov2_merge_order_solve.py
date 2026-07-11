"""PROV-2 decisive test (verifier): solve the current_bound_best theta at Nb=40
under the dirty-tree merge order (income wins -> Rouwenhorst) and the HEAD merge
order (profile wins -> stale 5-point grid); compare losses and moments.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import intergen_housing_fertility.local_panel as lp
from intergen_housing_fertility.calibration import run_model_cp_dt
from intergen_housing_fertility.local_panel import base_overrides, diagnostic_loss, extract_moments
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

record_path = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/"
    "full_audit_20260711/repro/current_repro_exact.json"
)
theta = dict(json.loads(record_path.read_text())["theta"])

targets, weights = lp.get_target_set(PRODUCTION_TARGET_SET)
targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
weights["aggregate_mean_occupied_rooms_18_85"] = 6.0

income = lp.income_process_overrides(5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO)
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
base = base_overrides(J=PRODUCTION_J, Nb=NB, n_house=5, max_iter_eq=PRODUCTION_MAX_ITER_EQ)

ov_dirty = {**base, **extra, **income, **theta}   # dirty tree: income wins
ov_head = {**base, **income, **extra, **theta}    # HEAD: profile wins

diff_keys = sorted(
    k for k in set(ov_dirty) | set(ov_head)
    if not np.array_equal(np.asarray(ov_dirty.get(k)), np.asarray(ov_head.get(k)))
)
print("keys differing between merged dicts:", diff_keys)

out = {"theta_source": str(record_path), "Nb": NB, "diff_keys": diff_keys}
for label, ov in [("dirty_income_wins", ov_dirty), ("head_profile_wins", ov_head)]:
    sol, P, p_eq = run_model_cp_dt(ov, verbose=False)
    moments = extract_moments(sol, P)
    loss = diagnostic_loss(moments, targets=targets, weights=weights)
    resid = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
    out[label] = {
        "loss": float(loss),
        "market_residual": resid,
        "p_eq": [float(x) for x in np.atleast_1d(p_eq)],
        "z_grid_effective": [float(x) for x in np.asarray(P.z_grid)],
        "moments": {k: float(v) for k, v in moments.items() if k in targets},
    }
    print(f"{label}: loss={loss:.6f} residual={resid:.3e} "
          f"z_grid_max={max(out[label]['z_grid_effective']):.4f}")

print("\nmoment                                            target    dirty     head   dcontrib")
for name in sorted(targets):
    md = out["dirty_income_wins"]["moments"].get(name)
    mh = out["head_profile_wins"]["moments"].get(name)
    w, t = weights[name], targets[name]
    dc = w * (mh - t) ** 2 - w * (md - t) ** 2
    print(f"  {name:45s} {t:9.4f} {md:9.4f} {mh:9.4f} {dc:+9.4f}")
print(f"\nloss(head) - loss(dirty) = "
      f"{out['head_profile_wins']['loss'] - out['dirty_income_wins']['loss']:+.6f}")

outfile = Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/"
    "full_audit_20260711/spec_audit/prov2_merge_order_solve_nb40.json"
)
outfile.parent.mkdir(parents=True, exist_ok=True)
outfile.write_text(json.dumps(out, indent=2, sort_keys=True))
print("wrote", outfile)
