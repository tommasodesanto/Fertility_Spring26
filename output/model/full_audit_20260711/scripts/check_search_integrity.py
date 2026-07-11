#!/usr/bin/env python3
"""Search-integrity audit checks (no model solves).

Verifies, using only importable code and the local report JSON:
 1. bounds_for_arm('baseline_pattern') == PRODUCTION_SEARCH_BOUNDS exactly.
 2. Enumerates the relaxed-arm boxes.
 3. global_unit_from_theta silently CLIPS out-of-box seed thetas (no error).
 4. beta bound-scale identity: annual^4 == reported 4-year beta estimate.
 5. Shared-coordinate / pattern-step arithmetic between the two wave-2 winners.
 6. Wave-2 global-stage budget arithmetic.
 7. jitter_theta determinism and distinctness across task seeds.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/overnight_combined_20260710")
import intergen_housing_fertility.local_panel as lp
from intergen_housing_fertility.production_profile import PRODUCTION_SEARCH_BOUNDS
from run_overnight_combined import bounds_for_arm, jitter_theta

out = {}

# 1. baseline bounds identity
base = bounds_for_arm("baseline_pattern")
out["baseline_equals_production_bounds"] = (
    [tuple(b) for b in base] == [tuple(b) for b in PRODUCTION_SEARCH_BOUNDS]
)

# 2. relaxed boxes
out["preference_relaxed"] = bounds_for_arm("preference_relaxed")
out["housing_relaxed"] = bounds_for_arm("housing_relaxed")

# 3. clipping behavior: theta outside the production box
theta_out = {name: hi * 1.5 for name, lo, hi in base}
theta_out["beta"] = 0.999**4
unit = lp.global_unit_from_theta(theta_out, bounds=base)
out["out_of_box_seed_silently_clipped"] = bool(unit is not None and float(unit.max()) == 1.0)

# 4. beta identity from report
rep = json.loads(Path(
    "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/"
    "combined_recalibration/overnight_20260710_report/summary.json").read_text())
cb = rep["current_bound_best"]
hr = rep["housing_relaxed_best"]
out["beta_annual_pow4_matches"] = abs(0.9576732996982973**4 - cb["theta"]["beta"]) < 1e-15

# 5. shared coordinates and pattern-step arithmetic
shared = [k for k in cb["theta"] if cb["theta"][k] == hr["theta"][k]]
out["identical_coordinates"] = shared
d = {k: hr["theta"][k] - cb["theta"][k] for k in cb["theta"] if k not in shared}
out["coordinate_diffs"] = d
# pattern steps on the unit cube: 0.08 * 0.6^k, scaled by (hi-lo)
steps = [0.08 * 0.6**k for k in range(8)]
def step_multiple(diff, rng):
    u = abs(diff) / rng
    for s in steps:
        m = u / s
        if abs(m - round(m)) < 1e-9 and round(m) != 0:
            return f"{round(m)} x step({s:.6g})"
    return None
ranges = {"H0": 9.0, "kappa_fert": 11.0, "h_bar_jump": None, "alpha_cons": 0.55,
          "psi_child": 0.35, "theta_n": 1.5, "chi": None, "h_bar_0": None, "h_bar_n": None}
out["step_decomposition"] = {
    "H0": step_multiple(d.get("H0", 0), 9.0),
    "kappa_fert": step_multiple(d.get("kappa_fert", 0), 11.0),
    "alpha_cons": step_multiple(d.get("alpha_cons", 0), 0.55),
    "psi_child": step_multiple(d.get("psi_child", 0), 0.35),
    "theta_n": step_multiple(d.get("theta_n", 0), 1.5),
}

# 6. wave-2 budget arithmetic
minutes = 330.0
global_deadline_min = min(100.0, minutes * 0.28)
out["wave2_global_deadline_min"] = global_deadline_min
out["draws_at_37s"] = global_deadline_min * 60 / 37.0
out["draws_at_74s"] = global_deadline_min * 60 / 74.0
out["draws_at_150s"] = global_deadline_min * 60 / 150.0

# 7. jitter determinism/distinctness
seed_theta = dict(cb["theta"])
t1 = jitter_theta(seed_theta, base, 2026071105, 0.05)
t1b = jitter_theta(seed_theta, base, 2026071105, 0.05)
t2 = jitter_theta(seed_theta, base, 2026071111, 0.05)
t0 = jitter_theta(seed_theta, base, 2026071101, 0.0)
out["jitter_deterministic"] = t1 == t1b
out["jitter_distinct_across_seeds"] = t1 != t2
out["jitter_zero_returns_seed"] = all(
    abs(t0[k] - v) < 1e-12 for k, v in lp.theta_from_global_unit(
        lp.global_unit_from_theta(seed_theta, bounds=[*base, ("H0", 1.0, 10.0)]),
        bounds=[*base, ("H0", 1.0, 10.0)]).items())

print(json.dumps(out, indent=2, default=str))
